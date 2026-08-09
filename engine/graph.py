"""Building the task graph in two stages, because its shape is not known in advance.

Step 0 decides, per sample, whether a supplied matrix is accepted or must be rebuilt from FASTQ.
A graph assembled before that answer would either place an alignment that is not needed or omit
one that is, and both are silent: the first wastes hours, the second produces a deliverable built
on a matrix nobody verified.

So the run is staged. Stage one is ingest alone. Its results decide stage two, which is
everything else. The alternative - putting both branches in the graph and skipping one - was
rejected because a skipped task and a task that was never required record identically in the
manifest, and the manifest is what a later reader trusts.
"""

from __future__ import annotations

from pathlib import Path

from . import steps
from .task import Refusal, Task, TaskFailure


def ingest_stage(pipeline, python_exe: str) -> list[Task]:
    """Stage one: validate and verify every sample. Nothing is produced, only decided."""
    out = []
    for row in pipeline.samples:
        s = row["sample"]
        out.append(Task(
            key=f"00_ingest/{s}", step="00_ingest", sample=s, fn=steps._ingest,
            inputs=tuple(str(x) for x in (row.get("matrix"), row.get("fastq_r1")) if x),
            params={"row": row, "python_exe": python_exe,
                    "expected_genes": row.get("expected_genes")},
        ))
    return out


def main_stage(pipeline, python_exe: str, tools: dict, ingest: dict) -> list[Task]:
    """Stage two, built from what stage one decided.

    `ingest` maps sample -> the TaskResult metrics of its ingest task, which carry `mode`
    ("accept" or "run") and the processor that would rebuild it.
    """
    tasks: list[Task] = []
    design = steps._design(pipeline.samples)
    by_sample = {r["sample"]: r for r in pipeline.samples}
    assay = {s: str(r.get("assay") or "").lower() for s, r in by_sample.items()}

    raw_of: dict = {}
    for s, row in by_sample.items():
        mode = (ingest.get(s) or {}).get("mode")
        if mode == "run":
            raw = pipeline.work / f"{s}_align" / "raw"
            tasks.append(Task(
                key=f"00_align/{s}", step="00_align", sample=s, fn=steps._align,
                inputs=(str(row.get("fastq_r1") or ""),),
                params={"row": row, "tools": tools, "processor":
                        (ingest.get(s) or {}).get("processor")},
                outputs=(str(raw),), cpus=20, memory_gb=64, walltime_h=12,
            ))
            raw_of[s] = raw
        else:
            raw_of[s] = Path(str(row.get("matrix")))

    # --- step 1: ambient. Mandatory for nuclei, optional for cells; the module decides, not us.
    ambient_keys = []
    for s in by_sample:
        needs = (f"00_align/{s}",) if f"00_align/{s}" in {t.key for t in tasks} \
            else (f"00_ingest/{s}",)
        h5 = pipeline.results / "objects" / f"{s}_ambient.h5"
        k = f"01_ambient/{s}"
        tasks.append(Task(
            key=k, step="01_ambient", sample=s, fn=steps._ambient, needs=needs,
            inputs=(str(raw_of[s]),),
            params={"raw": str(raw_of[s]), "assay": assay.get(s),
                    "exe": tools.get("cellbender", "cellbender"),
                    "env_bin": tools.get("cellbender_bin"),
                    "device": tools.get("device", "cuda"),
                    "fpr": 0.0, "learning_rate": tools.get("learning_rate")},
            outputs=(str(h5),), cpus=4, memory_gb=64, walltime_h=8, gpu=True,
        ))
        ambient_keys.append(k)

    tasks.append(Task(
        key="01_ambient_audit", step="01_ambient", fn=steps._ambient_audit,
        needs=tuple(ambient_keys),
        params={"per_sample": {s: {"h5": str(pipeline.results / "objects" / f"{s}_ambient.h5"),
                                   "raw": str(raw_of[s])} for s in by_sample},
                "design": design},
    ))
    tasks.append(Task(
        key="02_cells", step="02_cells", fn=steps._cellcall, needs=("01_ambient_audit",),
        params={"design": design, "samples": list(by_sample)},
    ))

    # --- steps 3-4: the light floor, then doublet scoring on what clears it.
    dbl_keys = []
    for s in by_sample:
        k = f"04_doublets/{s}"
        csv = pipeline.results / "tables" / f"{s}_doublets.csv"
        tasks.append(Task(
            key=k, step="04_doublets", sample=s, fn=steps._doublets, needs=("02_cells",),
            params={"sample": s, "h5": str(pipeline.results / "objects" / f"{s}_ambient.h5"),
                    "rscript": tools.get("rscript", "Rscript"),
                    "light_floor": tools.get("light_floor", 200),
                    "dbr": tools.get("dbr"), "dbr_sd": tools.get("dbr_sd"),
                    "seed": tools.get("seed", 0)},
            outputs=(str(csv),), cpus=4, memory_gb=32, walltime_h=4,
        ))
        dbl_keys.append(k)

    tasks.append(Task(
        key="04_doublet_health", step="04_doublets", fn=steps._doublet_health,
        needs=tuple(dbl_keys), params={"design": design, "samples": list(by_sample)},
    ))

    # --- step 5: thresholds, derived per library and applied as one cohort constant.
    tasks.append(Task(
        key="05_quality", step="05_quality", fn=steps._quality_stage,
        needs=("04_doublet_health",),
        params={"samples": list(by_sample), "python_exe": python_exe,
                "light_floor": tools.get("light_floor", 200),
                "decisions": pipeline.decisions},
    ))

    # --- step 6: cluster and profile, then flag.
    clus_keys = []
    for s in by_sample:
        k = f"06_cluster/{s}"
        tasks.append(Task(
            key=k, step="06_cluster_check", sample=s, fn=steps._cluster, needs=("05_quality",),
            params={"sample": s, "python_exe": python_exe,
                    "resolution": tools.get("resolution", 1.0), "seed": tools.get("seed", 0)},
            cpus=8, memory_gb=64, walltime_h=6,
        ))
        clus_keys.append(k)
    tasks.append(Task(
        key="06_cluster_flags", step="06_cluster_check", fn=steps._cluster_flags,
        needs=tuple(clus_keys), params={"design": design, "decisions": pipeline.decisions},
    ))

    last = "06_cluster_flags"
    if pipeline.mode == "apply":
        # Placed ONLY in apply mode. In evidence mode this task does not exist, so there is no
        # code path from `--mode evidence` to a deletion - not a flag that defaults to safe.
        tasks.append(Task(
            key="07_apply", step="07_apply", fn=steps._apply, needs=("06_cluster_flags",),
            params={"decisions": pipeline.decisions, "samples": list(by_sample)},
        ))
        last = "07_apply"

    tasks.append(Task(key="report", step="report", fn=steps._report, needs=(last,),
                      params={"extra": {}}))
    return tasks
