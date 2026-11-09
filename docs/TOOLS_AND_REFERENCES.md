# Tools and references

## References — resolved by name, through the project's registry

A samplesheet declares a reference by a `<species>/<build>` key, and the key is resolved through the
**project's** registry: a tab-separated file kept with the project and passed as `scqc run
--registry FILE` (and `scqc validate --registry FILE`). Paths are not accepted in the samplesheet,
because a path records where an index happened to sit on one machine and says nothing about what is
in it. `setup/init_project.sh` writes the file's header; the indexes themselves are large and live
outside any repository.

**scQC ships no registry and falls back to none.** Ingest refuses a run without `--registry`. Until
2026-11-09 the repository shipped `references/_registry/registry.tsv` with one row - the first
cohort's mouse reference - and ingest read it for any project that brought none, so such a project
was resolved against another cohort's reference without a word.

```
species        build                 path              built_by                   star_version  notes
mus_musculus   ensembl_112_filtered  <the index's path> celescope 2.7.3 rna mkref  2.7.11a       ...
```

That is a worked example of one complete row (tab-separated in the file). It was built by
`celescope 2.7.3 rna mkref` with STAR 2.7.11a from the Ensembl 112 primary assembly and a GTF
filtered by `celescope utils mkgtf` with its **default attributes, which retain introns**.
`conf/env/build_reference.sh` is the build, start to finish. It takes the organism, assembly and
release from `SPECIES`, `ASSEMBLY` and `ENSEMBL_RELEASE`, and refuses to start without all three: none
has a default.

That row records **34,290 genes**, and the count is the point of recording it. A gene count is the
cheapest available check that an index is the one the registry claims: rebuild from the same
release with the same `mkgtf` attributes and it must come out identical. A different count means
the release, the filter attributes or the annotation differs from what was declared, and every
downstream matrix inherits that difference without saying so.

**Intron retention is the reference's most consequential property.** For single-nuclei data most
reads are intronic; an exon-only reference silently discards the majority of the signal and still
produces a matrix that looks like data. Any reference added to the registry must record whether
introns are retained, because it cannot be recovered from the index by inspection.

## Tools — the recipe is owned, the environment is not

**A conda environment cannot be moved.** Every console script in it carries an absolute
interpreter path in its shebang, written at install time; in a CeleScope environment that is over a
hundred files under `bin` alone. Relocating the directory — copying it to another machine, moving
it to a shared filesystem, renaming a parent directory — breaks all of them at once.

The failure is worse than it sounds, because of how CeleScope invokes STAR: **by bare name on
`PATH`**. A relocated environment therefore fails at the shebang with `STAR: command not found`,
which on a batch scheduler means a job that disappears having produced nothing and having consumed
its allocation. Nothing in the output says the environment was the cause.

So the pipeline owns the **recipe**, never the directory:

```
conf/env/
  requirements.lock.txt       core: the analysis stack, 55 exact pins, python 3.12.12
  cellbender.lock.txt         the denoiser: 165 exact pins, python 3.10.22, linux x86_64 only
  lock.sh                     reads, checks and applies a lock (sourced by the installers)
  install_cellbender.sh       the denoiser, from its lock: cellbender 0.3.2, torch 2.5.1, numpy 1.26.4
  install_cs_cb.sh            celescope 2.7.3; its cellbender component calls install_cellbender.sh
  install_rdoublet.sh         R 4.3.3, scDblFinder 1.16.0, xgboost, Matrix — the supported route
  install_scdblfinder.sh      superseded: the same packages built from source
  build_reference.sh          the reference build, start to finish
  conda_pkgs_v2.7.3.txt       the celescope solve
```

An environment is a build artifact. It is rebuilt from a lock file, not copied — which is the only
form in which it is portable to another machine at all. A lock holds every package of an
environment's closure at one exact version, and names its interpreter to the patch;
`setup/install_env.sh` applies it with `pip install --no-deps`, so nothing is resolved on the day of
the build, and then compares every installed version with it. `--check` reads the locks and builds
nothing. The installer builds under the prefix you give it; a component that already exists there is
checked against its pins and never modified, so re-running it cannot disturb an environment something
is currently using — and an environment that does not match is refused, not reused.

The denoiser's lock must pin torch 2.5.1 and numpy 1.26.4: the installer refuses a lock naming any other
version, or naming torch unpinned, before anything is built. Torch 2.5.1 is on record for the environment
that made the calibration's denoised objects; numpy 1.26.4 is the pin drafted with it, and `--check` says
so. The R environment's conda half is pinned only where it decides the scores (R, scDblFinder and its
Bioconductor release, xgboost, Matrix); the other CRAN packages and the BLAS are solved at build time.
A pinned R package is checked as R reports it, which is not always conda's number: conda-forge's
`r-xgboost` 1.7.6 is R's xgboost 1.7.6.1, and `r-matrix` 1.6_5 is Matrix 1.6.5.

## One version fact worth carrying

STAR **2.7.11a**, supplied by the celescope conda environment. It is not pinned independently
anywhere, which means a rebuilt environment could silently arrive with a different STAR. The
registry records the version each reference was built with, so an index and an aligner that no
longer match can at least be detected.
