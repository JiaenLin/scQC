# Running scQC on PBS Pro

How to run a multi-library cohort as PBS jobs, end to end.

scQC's `--executor pbs` turns each per-library step into its own PBS job and waits for it. One
orchestrator job submits and supervises; the heavy per-library work fans out. This document is the
procedure: what the environment must provide, how to size the allocation, what the job script has
to do, and how to check the run afterwards.

For what the steps compute, see [`WORKFLOW.md`](WORKFLOW.md); for what they write,
[`OUTPUTS.md`](OUTPUTS.md).

---

## 1. Before the first run — check the interpreter

**This is the single most important check, and it is one line.**

```bash
python -c "import pandas as pd; print(pd.get_option('mode.string_storage'))"
```

pandas 3.0 defaults this to `auto`, which means *use pyarrow if it can be imported*. Under
`--executor pbs`, scQC's worker threads read `.h5ad` objects while sibling threads fork for
`qsub`/`qstat`, and building Arrow-backed string arrays in that situation crashes the process with
a SIGSEGV and no Python traceback.

scQC sets `mode.string_storage = "python"` for you at startup (`engine/native_io.harden()`, called
from `main()` before argument parsing), so a current version is safe. The check matters because it
tells you what you are looking at if you ever meet the crash on an older build: **the same code
crashes or does not depending on whether pyarrow is installed in the environment**, so an
environment change can break a pipeline that has not been touched. Bisecting the tool will not find
it. Check the environment first.

The rest of the environment:

| | |
|---|---|
| Python | the project's own interpreter, called by absolute path — no activation needed |
| `Rscript` | required for the doublet step; pass it with `--rscript` |
| PBS commands | `qsub`/`qstat` must be on `PATH` **inside the job**. scQC searches `PATH` and the usual install prefixes and discovers `PBS_CONF_FILE`/`PBS_EXEC`/`PBS_SERVER` itself, then probes with `qstat -B` and refuses up front if PBS is present but unconfigured — rather than failing once per task |

## 2. Size the orchestrator for the whole pipeline

> **`--executor pbs` sends only the SHELL half of a task to PBS. The Python half runs in the
> orchestrator's allocation.**

A task's declared `cpus`/`memory_gb`/`walltime_h` govern **the child job it submits**. Everything
scQC does in Python runs where you submitted it: reading each `.h5ad`, exporting the matrix the R
step consumes, every cohort-level step, and the whole report build.

So do not submit the orchestrator as a one-core helper.

| | orchestrator | per-library children |
|---|---|---|
| cpus | enough for the widest concurrent step, not for "just submitting" | scQC declares them |
| memory | enough for the largest single object plus the merge and the report | scQC declares them |
| walltime | **the whole run, including every child's queue wait** | scQC declares them |

The orchestrator is asleep in `qstat` polling for most of its life, but it holds the allocation the
whole time. Give it a walltime that covers the slowest possible path, not the compute time.

A workable starting point is *a few cores and enough memory to hold your largest library's object
several times over*, then adjust from what the run actually used — `qstat -xf <jobid>` reports
`resources_used.mem` and `resources_used.walltime` once it finishes. Measure your own cohort rather
than copying a number from someone else's.

**Discover your site's queues; do not assume any name.** `workq`, `short`, `long` and `batch` are
conventions, not standards, and their ceilings differ everywhere:

```bash
qstat -Q                 # one line per queue
qstat -Q -f <queue>      # its resources_max: ncpus, mem, walltime
```

Pick the queue from the requirement, never the requirement from the queue: over-requesting costs
some queue time, under-requesting costs the whole run and you find out at the end. `select` counts
**chunks** and `ncpus`/`mem` are **per chunk**, so `select=2:ncpus=8` is sixteen cores — the most
common reason `qsub` refuses a job, with a message that does not say so.

`--jobs N` sizes the orchestrator's thread pool, so it bounds both how many children exist at once
and how many threads are doing Python work concurrently. If your site caps queued jobs per user,
`--jobs` must stay under that cap.

## 3. The job script

Start from your site's PBS template if it has one. Every job script needs these regardless:

```bash
#!/bin/bash
#PBS -N scqc_run
#PBS -q <your-queue>                       # from `qstat -Q`, not from a guess
#PBS -l select=1:ncpus=<N>:mem=<M>gb      # per CHUNK; check against the queue ceiling
#PBS -l walltime=<HH:MM:SS>              # the server default is 5 minutes
#PBS -j oe

set -euo pipefail            # without it a failed step runs on and PBS records exit 0
cd "$PBS_O_WORKDIR"          # jobs start in $HOME, not where you submitted

RUNDIR="${RUNDIR:?pass RUNDIR with qsub -v}"
mkdir -p "$RUNDIR"/{logs,cache}

# Keep every tool's scratch inside the run directory. R, pip, matplotlib and numba all write
# under $HOME otherwise. Leave TMPDIR alone and never reassign $HOME.
export PYTHONNOUSERSITE=1 \
       XDG_CACHE_HOME="$RUNDIR/cache" \
       MPLCONFIGDIR="$RUNDIR/cache/mpl" \
       NUMBA_CACHE_DIR="$RUNDIR/cache/numba" \
       R_LIBS_USER="$RUNDIR/cache/R"

export OMP_NUM_THREADS="${NCPUS:-1}" MKL_NUM_THREADS="${NCPUS:-1}"

# qsub/qstat must be on PATH INSIDE the job, and a batch job does not inherit the login shell's
# environment. The module name is site-specific - `module avail` to find yours.
module load <pbs-module> 2>/dev/null || true

"$PY" "$TOOL/scqc_cli.py" run \
    --project "$RUNDIR" --samplesheet "$RUNDIR/samplesheet.csv" \
    --mode apply \
    --executor pbs --queue <child-queue> \
    --rscript "$RSCRIPT" --cpu \
    --dbr <rate> --dbr-sd <sd> \
    --jobs <N> --seed 0 \
  2>&1 | tee "$RUNDIR/logs/run.log"
```

`--project` is the run directory; scQC writes `results/`, `work/` and `logs/` beneath it. Never let
it fall back to a default output location.

### Logs

`#PBS -o` and `-e` belong inside the run directory, and **the submitting command must create that
directory first** — PBS resolves those paths at job exit, and if the directory is missing the copy
fails with `Exit_status = 0` and no file at all. Never point them at `/tmp`, which is node-local: a
log written on one node is invisible from every other.

### Sealing the run directory

Seal from an exit trap, so a crash seals too. **Turn `set -e` off inside the trap first:**

```bash
finish() { s=$?
  set +e                                  # REQUIRED — see below
  miss=""
  D=$(ls -d "$RUNDIR"/results/*/ 2>/dev/null | head -1)
  if [ -z "$D" ]; then
    miss="$miss results/<digest>"
  else
    # THE PATHS MUST BE THE ONES THE TOOL ACTUALLY WRITES. `report.json` is under `reports/`,
    # not at the top of the digest directory; `INPUTS.json` is at the top. Getting either wrong
    # fails a run that produced everything — see the note below.
    [ -s "$D/reports/report.json" ] || miss="$miss reports/report.json"
    [ -s "$D/INPUTS.json" ]         || miss="$miss INPUTS.json"
    n=$(ls "$D"/objects/cohort_per_sample/*.filtered.h5ad 2>/dev/null | wc -l)
    [ "$n" -gt 0 ] || miss="$miss filtered_objects"
  fi
  if [ -n "$miss" ] && [ "$s" -eq 0 ]; then s=97; fi        # `if`, not an && chain
  { echo "exit=$s"; echo "jobid=${PBS_JOBID:-none}"; hostname
    date -u +%Y-%m-%dT%H:%M:%SZ
    [ -n "$miss" ] && echo "missing_products:$miss"; } \
    > "$RUNDIR/$([ "$s" -eq 0 ] && echo SEALED.txt || echo FAILED.txt)"
  # Files, not -R on the tree, and leave logs/ writable.
  find "$RUNDIR" -type f -not -path "$RUNDIR/logs/*" -exec chmod a-w {} +
  return 0; }
trap finish EXIT
```

Three things this encodes:

1. **`[ -n "$miss" ] && [ "$s" -eq 0 ] && s=97` under `set -e` kills the trap.** When a run both
   fails and loses products — the case the seal exists for — that `&&` list returns non-zero and
   the function dies before writing anything, leaving a directory that reads as unfinished rather
   than failed.
2. **`chmod -R a-w` on the whole tree races PBS's own copy of `-o`/`-e` into `logs/`** and destroys
   the log. Chmod files and leave `logs/` writable.
3. **A seal must check the products, not just the exit status.** A job without `set -e` ends on a
   successful `echo`, and a trap reading `$?` then writes a green seal on a run that died.
4. **…but check for the paths the tool really writes, and verify them against a finished run
   before trusting the seal.** A product check listing `<digest>/report.json` — one directory
   level off from the `<digest>/reports/report.json` scQC actually writes — turned a completely
   successful run into `FAILED.txt`, with every object and figure present on disk. **A seal that
   fails a correct run is worse than no seal**, because the next person learns to ignore it. Run
   the check against a known-good run directory once, and confirm it says nothing is missing.

## 4. Submitting

```bash
RUNDIR=<project>/runs/scqc/<stage>/<UTCSTAMP>__scqc-<commit>__<stage>
mkdir "$RUNDIR"; mkdir "$RUNDIR/logs"      # bare mkdir: a key collision must abort, not merge
cp <samplesheet> "$RUNDIR/samplesheet.csv"

qsub -v RUNDIR="$RUNDIR",TOOLHEAD="$(git -C <tool> rev-parse HEAD)" \
     -o "$RUNDIR/logs/driver.log" \
     <job>.pbs
```

**Record the job id the same day.** PBS keeps job history for 24 hours; after that the log file in
the run directory is the only evidence the job existed.

Note also that a **compute node has no `git`** — pass the commit in with `-v` and write it into the
run directory rather than resolving it inside the job.

## 5. Watching a run

```bash
qstat -u "$USER"                       # the orchestrator plus its children
qstat -s <jobid>                       # the scheduler's own comment on why something is not starting
grep -E '^  (RUN|REFUSE|FAIL)' "$RUNDIR/logs/run.log"
```

Children appear and disappear as steps complete; the orchestrator stays in `R` throughout. A run
that has stopped making progress but still holds its allocation is usually waiting on a child that
cannot be scheduled — `qstat -s` on the child says why.

## 6. Checking the result

A finished run leaves `SEALED.txt` (or `FAILED.txt`) and `results/<digest>/`. Check the seal first,
then the products it names.

**Comparing two runs.** `tables/<sample>.percell.csv` carries `barcode` and `keep`, so barcode
identity is a text comparison — no `anndata`, no script, seconds on a login node.

> **scQC writes CSV through Python's `csv` module, whose default line terminator is `\r\n`.** Every
> table is CRLF, including on Linux. A comparison that does not strip it fails **silently**:
>
> ```bash
> grep -c ',True$' tables/<sample>.percell.csv                        # 0 — WRONG
> cut -d, -f16 tables/<sample>.percell.csv | tr -d '\r' | grep -cx True   # right
> ```
>
> The first form reports zero kept barcodes in every library and raises nothing.

**What should be identical between two runs of the same commit** — the kept-barcode set. The apply
criteria are `fail_not_cellbender_cell`, `fail_umi_floor`, `fail_gene_floor`, `fail_mito_ceiling`
and `fail_doublet`; they read counts, genes, mitochondrial percentage and a seeded doublet call, and
the count valleys and Tukey-calibrated MAD k are deterministic. None of them depends on the
executor. **A one-barcode difference between two runs of the same commit is a finding, not noise.**

**What may legitimately differ** — `cluster_FLAG`. Leiden is not bit-reproducible. The flag removes
nothing (`uns["scqc"]["removed_on_flag"]` is `0`), so it does not change the deliverable's
membership; do not use it as a reproducibility check.

A run key records **provenance, not results**. Where two runs at the same commit disagree, say in
writing which one the deliverable is.

## 7. Known limitations

- **A finished child's log can be reported missing on NFS.** The executor waits for the log with a
  `Path.exists()` retry loop; on NFS a `stat()` for a name the client has already looked up and not
  found is answered from the negative dentry cache, revalidated against the parent directory rather
  than by asking about the file. Repeated `exists()` calls therefore re-read one cached answer. It
  needs the orchestrator to be on a different node from the task — the normal condition under
  `--executor pbs`, and never reproduced by a local-executor test. If you see *"PBS job N finished
  but <log> does not exist"* while the log is plainly on disk, this is it. Listing the parent
  directory on each attempt is what makes the retry a retry.
- **`pgrep -f '<pattern>'` issued over SSH matches its own command line**, so a check for "is scQC
  running on a login node" reports a violation every time. Put such checks in a script file and
  exclude `$$`.
