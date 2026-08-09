#!/usr/bin/env Rscript
# Execution adapter: score one sample with scDblFinder and write the calls. It removes nothing,
# it selects nothing, and it refuses rather than repairing anything it did not expect.
#
# WHAT THIS SCRIPT IS FOR
#
# scDblFinder is an R package and the pipeline that needs it is Python. The obvious bridge -
# rpy2, or an R session that reads .h5ad - couples the two runtimes and makes the doublet step
# fail for reasons that have nothing to do with doublets. This script is the other bridge: a
# MatrixMarket triple in, a three-column CSV out, one process, no shared address space. The
# Python side writes the triple and reads the CSV; nothing else crosses.
#
# WHY IT FAILS INSTEAD OF COPING
#
# Every check below stops the run with a non-zero status rather than continuing on a repaired
# input, because the failure this pipeline exists to catch is the one that completes. Three of
# them are worth naming:
#
#   ORIENTATION. The export is genes x cells; AnnData is cells x genes. A transposed matrix is
#   still a valid matrix, scDblFinder still returns scores, and the scores are of nothing. The
#   only cheap detector is the barcode and feature counts, so both are checked against the
#   dimensions before anything is computed.
#
#   ZERO-COUNT BARCODES. scDblFinder's vignette records the failure directly: "Size factors
#   should be positive" if there are cells "that have zero reads (or a very low read count,
#   leading to zero after feature selection)". That error arrives from deep inside normalisation
#   and reads as a bug in the tool. It is an input that should have been floored, so it is
#   refused here, by name and with a count.
#
#   THE xgboost SHIM. scDblFinder 1.16.0 passes max_depth, eta, subsample, nthread and
#   eval_metric to xgboost as top-level arguments. xgboost 2.x and later accept them only
#   through a deprecation shim that folds them into `params`. The run completes and returns
#   scores that look entirely normal and were not computed the way the method was
#   characterised - the mismatch is visible only because a sibling package, scds, crashes on
#   the same change. The version is therefore printed always, and refused when the caller
#   passes a ceiling. conf/env/install_rdoublet.sh pins r-xgboost=1.7.6 for the same reason.
#
# WHY IT PRINTS ITS VERSIONS
#
# A version that was not observed is a fabricated provenance record, and it cannot be told apart
# from a real one. So every package version written into the report is one this process asked
# the installed library for and printed on stdout, where the calling adapter parses it back.
# Nothing is inferred from an environment name or a lock file.
#
# WHAT IT DOES NOT CLAIM TO HAVE OBSERVED
#
# The same rule cuts the other way and it caught this script out. It used to print
# `dbr_sd_package_default`, deparsed from `formals(scDblFinder)$dbr.sd`, under a name that says
# it is the value the installed package would have used had dbr.sd been omitted. It is not: the
# installed signature declares `dbr.sd = NULL` and resolves the effective value inside the
# function body, so what was being printed - and carried into the run's metrics, and available to
# be quoted in a report - was the string "NULL" wearing the name of a measurement. The formals
# ARE observable and are still printed, under `dbr_sd_formal_default`, beside
# `dbr_sd_formal_default_route`, which states in words whether that string is the value the
# package uses or only the placeholder in its signature. Where the effective default cannot be
# seen from here it is reported as not observed, rather than approximated by the placeholder.
#
# WHY IT DELETES THE OUTPUT BEFORE IT SCORES
#
# A calls file left by an earlier run satisfies every check a caller can make after this script
# returns. It exists, it parses, its barcodes match the export - because the export is the same -
# and the previous parameters are then recorded under the new ones. So out_csv, and the .partial
# it is renamed from, are removed at startup: after this script exits 0 the file at out_csv was
# written by THIS run, and after any other exit there is no file there at all.
#
# USAGE
#
#   Rscript scdblfinder.R <mtx_dir> <out_csv> <dbr> <dbr.sd|default> <seed> [key=value ...]
#
#     mtx_dir   directory holding matrix.mtx[.gz], barcodes.tsv[.gz], features.tsv[.gz],
#               written genes x cells by adapters/doublets.py
#     out_csv   destination; written atomically via a .partial file and renamed
#     dbr       expected doublet rate. DECLARED by the platform; there is no default here and
#               there is no fallback to scDblFinder's 10x loading formula, which on a Singleron
#               cohort tracked library size at r = 0.872
#     dbr.sd    the uncertainty on that rate, or the literal token `default` to omit the
#               argument entirely and let the installed package apply its own
#     seed       integer; set.seed() is called with it immediately before scDblFinder
#
#   optional key=value tokens:
#     threads=<n>        BiocParallel workers; >1 uses MulticoreParam(RNGseed = seed)
#     xgboost_max=<ver>  refuse (status 4) if the installed xgboost is at or above this version
#     features_col=<n>   column of features.tsv to use as row names (default 1)
#
# EXIT STATUS
#   0  scored, written, verified
#   1  a package is missing, or an input is unreadable, malformed or unsafe to score
#   2  the command line is wrong
#   4  the installed xgboost is at or above the refused version

options(warn = 1, stringsAsFactors = FALSE)

die <- function(status, ...) {
    cat("scdblfinder.R: FATAL: ", paste0(..., collapse = ""), "\n", sep = "", file = stderr())
    quit(save = "no", status = status, runLast = FALSE)
}

vline <- function(name, value) cat(sprintf("##scqc-version\t%s\t%s\n", name, value))
mline <- function(name, value) cat(sprintf("##scqc-metric\t%s\t%s\n", name, value))

# ---------------------------------------------------------------- command line

USAGE <- paste("usage: scdblfinder.R <mtx_dir> <out_csv> <dbr> <dbr.sd|default> <seed>",
               "[threads=n] [xgboost_max=ver] [features_col=n]")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) die(2, USAGE)

mtx_dir <- args[[1]]
out_csv <- args[[2]]
dbr_arg <- args[[3]]
sd_arg  <- args[[4]]
seed_arg <- args[[5]]

opts <- list(threads = "1", xgboost_max = "", features_col = "1")
for (tok in args[-(1:5)]) {
    if (!grepl("=", tok, fixed = TRUE)) die(2, "expected key=value, got '", tok, "'. ", USAGE)
    k <- sub("=.*$", "", tok)
    v <- sub("^[^=]*=", "", tok)
    if (!(k %in% names(opts))) {
        die(2, "unknown option '", k, "'. Known: ", paste(names(opts), collapse = ", "),
            ". An option the script does not recognise is silently ignored by most tools, ",
            "which is how a run acquires a setting nobody applied.")
    }
    opts[[k]] <- v
}

# A blank cell read out of a table arrives as the string "NA", and as.numeric("NA") is NA, not
# an error. Every numeric argument goes through this, so a missing value stops the run here
# rather than becoming a threshold further in.
as_num <- function(x, nm) {
    if (!nzchar(x)) die(2, nm, " is empty. Unknown is not a value; supply it or do not run.")
    v <- suppressWarnings(as.numeric(x))
    if (length(v) != 1L || is.na(v) || !is.finite(v)) {
        die(2, nm, " is not a finite number: '", x, "'. Values such as NA, NaN, Inf or a blank ",
            "table cell are refused rather than defaulted.")
    }
    v
}

dbr <- as_num(dbr_arg, "dbr")
if (dbr <= 0 || dbr >= 1) {
    die(2, "dbr must be a fraction strictly between 0 and 1, got ", dbr,
        ". A percentage such as 6 rather than 0.06 lands here.")
}

sd_is_default <- identical(tolower(sd_arg), "default")
dbr_sd <- if (sd_is_default) NULL else as_num(sd_arg, "dbr.sd")
if (!is.null(dbr_sd) && dbr_sd < 0) die(2, "dbr.sd must be >= 0, got ", dbr_sd)

seed_num <- as_num(seed_arg, "seed")
if (seed_num != trunc(seed_num)) die(2, "seed must be an integer, got ", seed_arg)
seed <- as.integer(seed_num)

threads <- as_num(opts$threads, "threads")
if (threads < 1 || threads != trunc(threads)) die(2, "threads must be an integer >= 1")
threads <- as.integer(threads)

features_col <- as_num(opts$features_col, "features_col")
if (features_col < 1 || features_col != trunc(features_col)) {
    die(2, "features_col must be an integer >= 1")
}
features_col <- as.integer(features_col)

# ---------------------------------------------------------------- no leftovers
#
# Removed here, at the top, and not just overwritten at the end. Overwriting makes the file
# correct only when the run reaches the write; every other exit - a missing package, an
# unreadable matrix, a killed job - leaves the PREVIOUS run's calls sitting at out_csv, where
# they satisfy every check a caller can make afterwards and are read as this run's result. After
# this point the file at out_csv either was written by this run or does not exist.
partial <- paste0(out_csv, ".partial")
for (leftover in c(out_csv, partial)) {
    if (file.exists(leftover)) {
        if (!file.remove(leftover)) {
            die(1, "a file from an earlier run is in the way and could not be removed: ",
                leftover, ". It is removed before scoring rather than overwritten after, ",
                "because a file that survives a failed run is indistinguishable from one this ",
                "run wrote.")
        }
    }
}

# ---------------------------------------------------------------- packages and versions

NEEDED <- c("Matrix", "SingleCellExperiment", "scDblFinder", "xgboost")
for (p in NEEDED) {
    if (!requireNamespace(p, quietly = TRUE)) {
        die(1, "R package '", p, "' is not installed in this library tree (",
            paste(.libPaths(), collapse = ", "), "). scQC does not substitute a different ",
            "detector or a different version - the result would not be the one the report ",
            "describes. Build the environment with conf/env/install_rdoublet.sh.")
    }
}

vline("R", R.version.string)
for (p in NEEDED) vline(p, as.character(utils::packageVersion(p)))

xgb <- utils::packageVersion("xgboost")
if (nzchar(opts$xgboost_max)) {
    if (xgb >= opts$xgboost_max) {
        die(4, "xgboost ", as.character(xgb), " is at or above the refused version ",
            opts$xgboost_max, ". scDblFinder passes max_depth, eta, subsample, nthread and ",
            "eval_metric as top-level arguments; from 2.0.0 xgboost accepts those only through ",
            "a deprecation shim, and the run then completes and returns scores that were not ",
            "computed the way the method was characterised. Pin r-xgboost=1.7.6 ",
            "(conf/env/install_rdoublet.sh), or pass xgboost_max= empty to score anyway and ",
            "record that the scores went through the shim.")
    }
}

# What the installed scDblFinder() DECLARES for dbr.sd, and what that declaration is worth.
#
# This used to be printed as `dbr_sd_package_default` and described as the value the package
# would have used had dbr.sd been omitted. It is not that, and it never was: the installed
# signature declares `dbr.sd = NULL` and the effective value is resolved inside the function
# body, so the string being reported as an observed default was "NULL". An observation of the
# formals is a real observation - of the formals - and it is now labelled as one, with a second
# metric saying in words whether it determines the value used. Where it does not, this script
# says the effective default was NOT observed rather than offering the placeholder in its place;
# reading it out of the function body would be parsing source to guess at a number.
sd_formals <- formals(scDblFinder::scDblFinder)
# The formal is never bound to a variable. An argument declared WITHOUT a default is the empty
# symbol, and binding that to a name makes every later use raise 'argument is missing, with no
# default' - from inside this reporting block, about a function it is only describing.
if (!("dbr.sd" %in% names(sd_formals))) {
    sd_formal_txt <- "absent"
    sd_formal_route <- paste0("the installed scDblFinder() has no dbr.sd argument at all, so ",
                              "omitting it selects nothing; check the package version above")
} else if (identical(sd_formals[["dbr.sd"]], quote(expr = ))) {
    sd_formal_txt <- "none"
    sd_formal_route <- paste0("the argument is declared with no default, so omitting it is ",
                              "an error rather than a choice of default")
} else {
    sd_formal_txt <- paste(deparse(sd_formals[["dbr.sd"]]), collapse = " ")
    if (identical(trimws(sd_formal_txt), "NULL")) {
        sd_formal_route <- paste0("NOT the value used: the signature default is NULL and the ",
                                  "effective dbr.sd is computed inside scDblFinder(); this run ",
                                  "did not observe it and does not guess it")
    } else {
        sd_formal_route <- paste0("the value used: the signature default is a literal, so it ",
                                  "is what scDblFinder() applies when dbr.sd is omitted")
    }
}
# The calling adapter parses '##scqc-metric<TAB>name<TAB>value' and requires exactly three
# fields, so a deparsed default carrying a tab or a newline would look like a malformed line
# rather than an unusual default.
sd_formal_txt <- gsub("[\t\r\n]+", " ", sd_formal_txt)
sd_formal_route <- gsub("[\t\r\n]+", " ", sd_formal_route)

# ---------------------------------------------------------------- input

if (!dir.exists(mtx_dir)) die(1, "input directory does not exist: ", mtx_dir)

pick <- function(stem) {
    plain <- file.path(mtx_dir, stem)
    gzd <- paste0(plain, ".gz")
    if (file.exists(plain)) return(list(path = plain, gz = FALSE))
    if (file.exists(gzd)) return(list(path = gzd, gz = TRUE))
    die(1, "neither ", plain, " nor ", gzd, " exists. The export written by ",
        "adapters/doublets.py holds matrix.mtx, barcodes.tsv and features.tsv, each optionally ",
        "gzipped.")
}

mtx_f <- pick("matrix.mtx")
bc_f <- pick("barcodes.tsv")
ft_f <- pick("features.tsv")

read_col <- function(f, col, label) {
    con <- if (f$gz) gzfile(f$path, "rt") else file(f$path, "rt")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    tab <- tryCatch(
        utils::read.delim(con, header = FALSE, sep = "\t", quote = "", comment.char = "",
                          colClasses = "character", check.names = FALSE),
        error = function(e) die(1, "could not read ", label, " from ", f$path, ": ",
                                conditionMessage(e)))
    if (ncol(tab) < col) {
        die(1, label, " (", f$path, ") has ", ncol(tab), " column(s); column ", col,
            " was requested.")
    }
    trimws(tab[[col]])
}

barcodes <- read_col(bc_f, 1L, "barcodes.tsv")
features <- read_col(ft_f, features_col, "features.tsv")

read_mtx <- function(f) {
    # A connection opened at script top level is never closed by on.exit - on.exit only fires
    # when a function returns - so the open and the close are wrapped in one here.
    if (!f$gz) {
        return(tryCatch(Matrix::readMM(f$path),
                        error = function(e) die(1, "could not read ", f$path,
                                                " as MatrixMarket: ", conditionMessage(e))))
    }
    con <- gzfile(f$path, "rb")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    tryCatch(Matrix::readMM(con),
             error = function(e) die(1, "could not read ", f$path, " as MatrixMarket: ",
                                     conditionMessage(e)))
}

m <- read_mtx(mtx_f)

# Orientation, checked before anything is computed. A transposed export is a valid matrix and
# produces scores of nothing; the dimensions are the only cheap way to see it.
if (nrow(m) != length(features) || ncol(m) != length(barcodes)) {
    die(1, "shape mismatch: matrix is ", nrow(m), " x ", ncol(m), " but features.tsv has ",
        length(features), " row(s) and barcodes.tsv has ", length(barcodes), " row(s). The ",
        "export must be genes x cells. A cells x genes matrix - AnnData's own orientation - ",
        "produces exactly this mismatch, and produces no error at all when the two counts ",
        "happen to be equal.")
}
if (length(barcodes) == 0L) die(1, "barcodes.tsv is empty; there is nothing to score")
if (anyDuplicated(barcodes) != 0L) {
    d <- unique(barcodes[duplicated(barcodes)])
    die(1, length(d), " duplicated barcode(s) in barcodes.tsv, e.g. ",
        paste(utils::head(d, 5), collapse = ", "), ". A duplicated key silently merges two ",
        "nuclei's calls on the way back into Python.")
}
if (any(!nzchar(barcodes))) die(1, sum(!nzchar(barcodes)), " empty barcode(s) in barcodes.tsv")
bad_chr <- grepl("[,\"\r\n]", barcodes)
if (any(bad_chr)) {
    die(1, sum(bad_chr), " barcode(s) contain a comma, quote or newline, e.g. '",
        barcodes[which(bad_chr)[1]], "'. The output is unquoted CSV; such a barcode would ",
        "split into two fields and shift every column after it.")
}

cs <- Matrix::colSums(m)
if (any(!is.finite(cs))) die(1, "non-finite column sums in the matrix; the export is not counts")
zero <- which(cs <= 0)
if (length(zero) > 0L) {
    die(1, length(zero), " barcode(s) have zero total counts, e.g. ",
        paste(utils::head(barcodes[zero], 5), collapse = ", "), ". scDblFinder's vignette ",
        "records the consequence: 'Size factors should be positive' if there are cells with ",
        "zero reads or a read count low enough to reach zero after feature selection. Apply ",
        "the light floor (step 3) before scoring; it selects the scoring set and removes ",
        "nothing from the analysis.")
}

# ---------------------------------------------------------------- score

m <- methods::as(m, "CsparseMatrix")
sce <- SingleCellExperiment::SingleCellExperiment(list(counts = m))
colnames(sce) <- barcodes
rownames(sce) <- features

set.seed(seed)

call_args <- list(sce, dbr = dbr, verbose = FALSE)
if (!is.null(dbr_sd)) call_args$dbr.sd <- dbr_sd
if (threads > 1L) {
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
        die(1, "threads=", threads, " was requested but BiocParallel is not installed. A ",
            "parallel run without a seeded backend is not reproducible, so this is refused ",
            "rather than silently run on one core.")
    }
    vline("BiocParallel", as.character(utils::packageVersion("BiocParallel")))
    call_args$BPPARAM <- BiocParallel::MulticoreParam(workers = threads, RNGseed = seed)
}

sce <- tryCatch(do.call(scDblFinder::scDblFinder, call_args),
                error = function(e) die(1, "scDblFinder failed: ", conditionMessage(e)))

out_bc <- colnames(sce)
if (length(out_bc) != length(barcodes) || !setequal(out_bc, barcodes)) {
    die(1, "scDblFinder returned ", length(out_bc), " column(s) for an input of ",
        length(barcodes), ". The detector contract forbids subsetting the input: a nucleus ",
        "that vanishes here would be recorded downstream as not a doublet, which is not the ",
        "same as never having been examined.")
}

score <- sce$scDblFinder.score
cls <- sce$scDblFinder.class
if (is.null(score) || is.null(cls)) {
    die(1, "scDblFinder returned no scDblFinder.score / scDblFinder.class column. Available ",
        "colData: ", paste(colnames(SummarizedExperiment::colData(sce)), collapse = ", "))
}
score <- as.numeric(score)
cls <- as.character(cls)
if (length(score) != length(out_bc) || length(cls) != length(out_bc)) {
    die(1, "scDblFinder returned ", length(score), " score(s) and ", length(cls),
        " class(es) for ", length(out_bc), " barcode(s)")
}
if (any(is.na(score)) || any(!is.finite(score))) {
    die(1, sum(is.na(score) | !is.finite(score)), " non-finite doublet score(s). A missing ",
        "score written out as a number reads downstream as a confident low one.")
}
if (any(is.na(cls))) die(1, sum(is.na(cls)), " missing doublet class(es)")
unexpected <- setdiff(unique(cls), c("singlet", "doublet"))
if (length(unexpected) > 0L) {
    die(1, "unexpected doublet class value(s): ", paste(unexpected, collapse = ", "),
        ". The reader accepts 'singlet' and 'doublet' only, rather than mapping an unknown ",
        "label onto one of them.")
}

# ---------------------------------------------------------------- write

out_dir <- dirname(out_csv)
if (nzchar(out_dir) && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(out_dir)) die(1, "could not create output directory ", out_dir)
}

df <- data.frame(barcode = out_bc,
                 doublet_score = sprintf("%.10g", score),
                 doublet_class = cls)

# Written to a sibling and renamed, so a job killed mid-write leaves no file rather than a
# short one that parses cleanly and is missing the last few thousand nuclei. Both this path and
# out_csv were removed at startup, so neither can be a survivor of an earlier run.
ok <- tryCatch({
    utils::write.csv(df, file = partial, row.names = FALSE, quote = FALSE)
    TRUE
}, error = function(e) die(1, "could not write ", partial, ": ", conditionMessage(e)))
if (!file.rename(partial, out_csv)) {
    die(1, "could not rename ", partial, " to ", out_csv)
}

n_doublet <- sum(cls == "doublet")
mline("n_cells", length(out_bc))
mline("n_genes", nrow(sce))
mline("n_doublets", n_doublet)
mline("n_singlets", length(out_bc) - n_doublet)
mline("dbr_used", sprintf("%.10g", dbr))
mline("dbr_sd_used", if (is.null(dbr_sd)) "package-default" else sprintf("%.10g", dbr_sd))
mline("dbr_sd_formal_default", sd_formal_txt)
mline("dbr_sd_formal_default_route", sd_formal_route)
mline("seed_used", format(seed))
mline("threads_used", format(threads))
mline("out_csv", out_csv)

quit(save = "no", status = 0, runLast = FALSE)
