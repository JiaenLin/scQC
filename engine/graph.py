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



def _num(*candidates):
    """The first candidate that is actually a number, else None.

    A samplesheet cell arrives as a string, and an empty one arrives as "" - which is not a
    value and must not shadow the cohort-wide flag behind it. `float("")` raises and `float(None)`
    raises, so both are treated as absent rather than as zero; a doublet rate of 0.0 declared by
    accident would silently disable the detector.
    """
    for c in candidates:
        if c is None:
            continue
        s = str(c).strip()
        if not s or s.lower() in ("na", "nan", "none", "null"):
            continue
        try:
            return float(s)
        except (TypeError, ValueError):
            continue
    return None


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
    #
    # A sample may also arrive ALREADY CORRECTED - a samplesheet `ambient_h5` column, with the
    # provenance columns beside it. That is a third state, not a skip (modules/01_ambient), and
    # it is how the pipeline is re-entered at step 2 with denoised matrices. No CellBender task
    # is created for such a sample and none is pretended: the object is used where step 1's
    # output would have gone, and the audit is told it cannot measure what it did not see.
    ambient_keys, supplied_of = [], {}
    for s, row in by_sample.items():
        given = str(row.get("ambient_h5") or "").strip()
        if not given:
            continue
        supplied_of[s] = {
            "path": given,
            "tool": str(row.get("ambient_tool") or "").strip(),
            "version": str(row.get("ambient_version") or "").strip(),
            "params": str(row.get("ambient_params") or "").strip(),
            "produced_by": str(row.get("ambient_produced_by") or "").strip()}

    for s in by_sample:
        h5 = pipeline.results / "objects" / f"{s}_ambient.h5"
        if s in supplied_of:
            continue
        # `needs` may only name tasks in THIS stage. Stage one - ingest - has already run to
        # completion before this function is called, so a dependency on it is both unresolvable
        # and unnecessary: the ordering it would express is already guaranteed.
        #
        # This read `else (f"00_ingest/{s}",)`, which is a cross-stage edge and stops the run with
        # "depends on unknown task(s)". It never fired because it only applies to samples whose
        # matrix is ACCEPTED rather than aligned - and until this cohort, no run had ever supplied
        # a matrix. A latent break on the path a reusable pipeline is most likely to be used on.
        needs = (f"00_align/{s}",) if f"00_align/{s}" in {t.key for t in tasks} else ()
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

    if supplied_of:
        # Runs BEFORE the audit and before anything consumes an ambient object: it validates the
        # provenance of every supplied matrix through the same module that decides whether
        # CellBender runs, so an unattributed object stops the run here rather than surfacing as
        # a blank field in the report.
        tasks.append(Task(
            key="01_ambient_supplied", step="01_ambient", fn=steps._ambient_supplied,
            # No `needs`: stage one has already completed, and a supplied matrix depends on
            # nothing built in this stage.
            params={"supplied": supplied_of,
                    "assay": {s: assay.get(s) for s in supplied_of}},
            outputs=tuple(str(pipeline.results / "objects" / f"{s}_ambient.h5")
                          for s in supplied_of),
        ))
        ambient_keys.append("01_ambient_supplied")

    tasks.append(Task(
        key="01_ambient_audit", step="01_ambient", fn=steps._ambient_audit,
        needs=tuple(ambient_keys),
        # EVERY LIBRARY WHOSE RAW MATRIX IS IN HAND IS AUDITED, supplied or not.
        #
        # This excluded every supplied library, on the stated grounds that a supplied object
        # arrives without the raw counts the removal fraction is measured against. That is true
        # of a bare object and false of this samplesheet: `matrix` names the raw counts for every
        # row, supplied or not, so `raw_of` is populated for all of them. The audit was declining
        # to measure something it had. On a cohort where all ten were supplied it covered NOTHING
        # and said "0 libraries" - honest about a fact it did not have to accept.
        #
        # A library with genuinely no raw is still excluded and still named, below: unmeasurable
        # and unmeasured must not read the same.
        params={"per_sample": {s: {"h5": str(pipeline.results / "objects" / f"{s}_ambient.h5"),
                                   "raw": str(raw_of[s]),
                                   "supplied": s in supplied_of,
                                   # Optional and DECLARED, never inferred from a sibling path.
                                   # CellBender's learning curve lives in its own log or metrics
                                   # file, and a converted object does not carry it.
                                   "log": str(by_sample[s].get("ambient_log") or "").strip(),
                                   "metrics_csv": str(
                                       by_sample[s].get("ambient_metrics") or "").strip(),
                                   # WHICH BARCODES THE FRACTION IS SUMMED OVER. A converted
                                   # object does not carry CellBender's cell-probability latent,
                                   # so on a supplied library the cell call has to be declared -
                                   # and it already is, in the column step 2 reads. Without it
                                   # parse_metrics refuses, correctly: summing over every droplet
                                   # instead would answer a different question in the same units.
                                   "cell_barcodes": str(
                                       by_sample[s].get("cellbender_barcodes") or "").strip()}
                               for s in by_sample if raw_of.get(s)},
                "no_raw": [s for s in by_sample if not raw_of.get(s)],
                "supplied": supplied_of,
                "design": design},
    ))
    # Step 2 needs each caller's BARCODES, not two counts: the number the gate turns on is how
    # many the aligner called that the denoiser did not.
    #
    # Where they come from depends on who ran what. When CellBender runs here it writes
    # `<stem>_cell_barcodes.csv` beside its output; when the object was SUPPLIED, this run never
    # saw that file and the samplesheet must say where it is. The aligner's call is its filtered
    # matrix directory - `outs/filtered` for CeleScope, `filtered_feature_bc_matrix` for
    # CellRanger - which the samplesheet names for an accepted matrix.
    #
    # Neither is defaulted or inferred from the raw matrix. The whole point of this step is to
    # catch a population lost at the boundary between the two callers, and a step that quietly
    # compares something else is worse than one that refuses.
    call_paths = {}
    for s, row in by_sample.items():
        cb_csv = str(row.get("cellbender_barcodes") or "").strip()
        if not cb_csv and s not in supplied_of:
            stem = pipeline.results / "objects" / f"{s}_cellbender"
            cb_csv = str(stem.parent / f"{stem.name}_cell_barcodes.csv")
        call_paths[s] = {"aligner": str(row.get("aligner_cells") or "").strip(),
                         "cellbender": cb_csv}

    tasks.append(Task(
        key="02_cells", step="02_cells", fn=steps._cellcall,
        # The denoiser's cell call comes from the measurement of <sample>_ambient.h5, which is the
        # object every step after step 1 reads. It used to come from a second file named in the
        # samplesheet, from a run nothing verified was this one.
        needs=tuple(f"05_quality/{s}" for s in by_sample),
        params={"design": design, "samples": list(by_sample), "call_paths": call_paths},
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
                    # Per-sample first: loading concentration can differ between libraries, and
                    # a cohort-wide rate imposed over that is a number describing no library.
                    # Falls back to the cohort flags; neither is defaulted, and the adapter
                    # refuses None rather than choosing a rate nobody declared.
                    "dbr": _num(by_sample[s].get("dbr"), tools.get("dbr")),
                    "dbr_sd": _num(by_sample[s].get("dbr_sd"), tools.get("dbr_sd")),
                    "seed": tools.get("seed", 0)},
            outputs=(str(csv),), cpus=4, memory_gb=32, walltime_h=4,
        ))
        dbl_keys.append(k)

    tasks.append(Task(
        key="04_doublet_health", step="04_doublets", fn=steps._doublet_health,
        needs=tuple(dbl_keys), params={"design": design, "samples": list(by_sample)},
    ))

    # --- step 5: thresholds, derived per library and applied as one cohort constant.
    # Step 5 measures each library independently and then proposes ONE cohort constant from the
    # ten results. The measurement was a serial loop inside the single 05_quality task, which no
    # amount of --jobs could parallelise - concurrency releases tasks, not iterations. It is now
    # ten tasks and a barrier, the same shape as step 4.
    for s in by_sample:
        tasks.append(Task(
            key=f"05_quality/{s}", step="05_quality", sample=s, fn=steps._quality_measure,
            # Runs right after step 1, not after the doublets: it MEASURES and removes nothing,
            # and both step 2 and step 5 need what it produces. Ordering constraints in this
            # pipeline are about what may FILTER before what, never about what may be looked at.
            needs=("01_ambient_audit",),
            params={"sample": s, "python_exe": python_exe,
                    "mt_prefix": str(by_sample[s].get("mt_prefix") or "").strip(),
                    "ribo_pattern": str(by_sample[s].get("ribo_pattern") or "").strip(),
                    # The floor the MITOCHONDRIAL quartiles are taken above. The count valleys in
                    # the same pass do not use it and must not: they need the debris mode this
                    # would delete.
                    "light_floor": tools.get("light_floor", 200)},
            cpus=4, memory_gb=32, walltime_h=2,
        ))

    tasks.append(Task(
        key="05_quality", step="05_quality", fn=steps._quality_stage,
        needs=tuple(f"05_quality/{s}" for s in by_sample) + ("04_doublet_health",),
        params={"samples": list(by_sample), "python_exe": python_exe,
                "light_floor": tools.get("light_floor", 200),
                "decisions": pipeline.decisions},
    ))

    # --- step 6: cluster and profile, then flag.
    clus_keys = []
    for s in by_sample:
        k = f"06_cluster/{s}"
        tasks.append(Task(
            key=k, step="06_cluster_check", sample=s, fn=steps._cluster,
            # 04_doublets/<s> is already upstream through 05_quality, and is named anyway: this
            # task now READS that library's doublet CSV, and a dependency that exists only by
            # transitivity is one a later reordering can remove without anything noticing.
            needs=("05_quality", f"04_doublets/{s}"),
            params={"sample": s, "python_exe": python_exe,
                    "resolution": tools.get("resolution", 1.0), "seed": tools.get("seed", 0),
                    "mt_prefix": str(by_sample[s].get("mt_prefix") or "").strip(),
                    "ribo_pattern": str(by_sample[s].get("ribo_pattern") or "").strip(),
                    "doublet_csv": str(pipeline.results / "tables" / f"{s}_doublets.csv")},
            # 48 rather than 64: a common workq ceiling is 50 gb, and one library's
            # clustering does not need more. check_resources() would refuse the
            # graph otherwise - correctly, but before doing any work.
            cpus=8, memory_gb=48, walltime_h=6,
        ))
        clus_keys.append(k)
    tasks.append(Task(
        key="06_cluster_flags", step="06_cluster_check", fn=steps._cluster_flags,
        needs=tuple(clus_keys), params={"design": design, "decisions": pipeline.decisions,
                                        "samples": list(by_sample)},
    ))

    last = "06_cluster_flags"
    if pipeline.mode == "apply":
        # Placed ONLY in apply mode. In evidence mode these tasks do not exist, so there is no
        # code path from `--mode evidence` to a deletion - not a flag that defaults to safe.
        #
        # TWO STAGES, because the gate has to run between them. `07_measure/<s>` computes every
        # criterion per barcode and writes a table; `07_apply` builds the ledger, puts it through
        # the gate, audits the result and only then materialises a filtered object. An op that
        # measured and wrote in one pass would have removed before the approval was checked.
        # The thresholds are NOT resolved here. With `per_library` they come from step 5's table,
        # which on a clean run does not exist while the graph is being built - resolving now
        # would refuse a correct run for a file that is three steps from being written. Each task
        # resolves its own when it runs, through the same function, so they cannot disagree.
        measure_keys = []
        for s in by_sample:
            k = f"07_measure/{s}"
            tasks.append(Task(
                key=k, step="07_apply", sample=s, fn=steps._apply_measure,
                needs=("06_cluster_flags", f"04_doublets/{s}"),
                params={"sample": s, "python_exe": python_exe,
                        "decisions": pipeline.decisions,
                        "mt_prefix": str(by_sample[s].get("mt_prefix") or "").strip(),
                        "ribo_pattern": str(by_sample[s].get("ribo_pattern") or "").strip(),
                        "light_floor": tools.get("light_floor", 200),
                        "doublet_csv": str(pipeline.results / "tables"
                                           / f"{s}_doublets.csv")},
                cpus=4, memory_gb=32, walltime_h=2,
            ))
            measure_keys.append(k)
        tasks.append(Task(
            key="07_apply", step="07_apply", fn=steps._apply,
            needs=tuple(measure_keys) + ("06_cluster_flags",),
            params={"decisions": pipeline.decisions, "samples": list(by_sample),
                    "python_exe": python_exe, "light_floor": tools.get("light_floor", 200)},
            cpus=8, memory_gb=48, walltime_h=4,
        ))
        last = "07_apply"

    tasks.append(Task(key="report", step="report", fn=steps._report, needs=(last,),
                      params={"extra": {}}))
    return tasks
