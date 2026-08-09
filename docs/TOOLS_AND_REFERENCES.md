# Tools and references

## References — resolved by name, owned here

```
references/
  _registry/registry.tsv              resolve a reference by NAME, not by path
  <species>/<build>/                  the index itself — large, and not tracked by git
```

A samplesheet declares a reference by its registry name. Paths are not accepted, because a path
records where an index happened to sit on one machine and says nothing about what is in it.

The registry ships one worked example, `mus_musculus / ensembl_112_filtered`, built by
`celescope 2.7.3 rna mkref` with STAR 2.7.11a from the Ensembl 112 primary assembly and a GTF
filtered by `celescope utils mkgtf` with its **default attributes, which retain introns**. It is an
example of a complete registry row, not a reference you are expected to use; the index itself is
not distributed. `conf/env/build_reference.sh` is the build, start to finish.

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
  install_cs_cb.sh            celescope 2.7.3 + cellbender 0.3.2
  install_rdoublet.sh         R with scDblFinder and scds — the supported route
  install_scdblfinder.sh      superseded: the same packages built from source
  build_reference.sh          the reference build, start to finish
  conda_pkgs_v2.7.3.txt       the celescope solve
  requirements.lock.txt       workstation
  requirements.hpc.lock.txt   HPC
```

An environment is a build artifact. It is rebuilt from a lock file, not copied — which is the only
form in which it is portable to another machine at all. `setup/install_env.sh` builds them under
the prefix you give it, and skips any component that already exists there, so re-running it cannot
disturb an environment something is currently using.

## One version fact worth carrying

STAR **2.7.11a**, supplied by the celescope conda environment. It is not pinned independently
anywhere, which means a rebuilt environment could silently arrive with a different STAR. The
registry records the version each reference was built with, so an index and an aligner that no
longer match can at least be detected.
