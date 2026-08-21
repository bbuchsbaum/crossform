#!/usr/bin/env Rscript
# 00-download.R -- fetch and checksum the Haxby 2001 subject tarballs.
#
# Fetches subj1 .. subj6 from the PyMVPA distribution. Idempotent at two
# levels: a tarball that is already present and hash-verified is never
# refetched, and a subject directory that already holds the expected files is
# never re-unpacked. A tarball whose digest does not match the published one
# is a hard stop, not a warning: the script refuses to hand an unknown archive
# to the rest of the exemplar.
#
# Partial downloads resume. Interrupt it and run it again.
#
# Environment:
#   SUBJECTS       space/comma separated subject ids to fetch, e.g. "subj2 subj5".
#                  Default: all six.
#   EXEMPLAR_DIR   override the exemplar directory (normally inferred).
#   KEEP_TARBALL   "0" to delete each tarball once it is verified and unpacked
#                  (halves the disk footprint; provenance keeps the digests).
#
# Provenance is recorded per subject in data/download-provenance-<subj>.rds and
# collected into data/download-provenance-all.rds. The pre-existing single-file
# data/download-provenance.rds (subject 1, written 2026-08-13) is left alone.

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  Sys.getenv("EXEMPLAR_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
data_dir <- file.path(exemplar_dir, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

BASE_URL <- "http://data.pymvpa.org/datasets/haxby2001/"
STAMP <- "2010.01.14"

## ---- Published MD5s ------------------------------------------------------
# Verbatim from the distribution's own checksum file,
#   http://data.pymvpa.org/datasets/haxby2001/MD5SUMS
# fetched 2026-08-20. These are the dataset maintainers' hashes, not ours: the
# file is the authority, and subject 1's entry below is the same value the
# tarball downloaded here on 2026-08-13 already hashes to. Do not edit a value
# to make a download pass; refetch the tarball.
EXPECTED_MD5 <- c(
  subj1 = "03e6865ec33bd5ebccdad3f13ef1d77f",
  subj2 = "56902b0583c0329b8364cadc1abb3ed5",
  subj3 = "fa6e21e52d161f108e5b45f9b6571d76",
  subj4 = "5c0a0e45b562112ed4f15da9f27c7176",
  subj5 = "2b983a6f4d16bad9053654e57c00483a",
  subj6 = "2fecf6b76e781b349e616c594a74064e"
)
MD5SUMS_URL <- paste0(BASE_URL, "MD5SUMS")
MD5SUMS_FETCHED <- "2026-08-20"

# SHA-256 of subj1-2010.01.14.tar.gz, recorded from the first verified
# download on 2026-08-13 from data.pymvpa.org (300.2 MB / 314803244 bytes,
# Last-Modified 2010-01-15). The published MD5SUMS carries no SHA-256 column,
# so this is the one digest of our own; the other five are recorded into
# provenance on first download rather than pinned here.
EXPECTED_SHA256 <- c(
  subj1 = "3c14bd7fad6c869e5d9b81739e24e21bb9feb6a7abef27682bbf131b6b4bec5c"
)

SUBJECTS_ALL <- names(EXPECTED_MD5)

subjects <- Sys.getenv("SUBJECTS")
subjects <- if (nzchar(subjects)) {
  trimws(strsplit(subjects, "[,[:space:]]+")[[1L]])
} else {
  SUBJECTS_ALL
}
subjects <- subjects[nzchar(subjects)]
unknown <- setdiff(subjects, SUBJECTS_ALL)
if (length(unknown)) {
  stop("Unknown subject id(s): ", paste(unknown, collapse = ", "),
       ". Known: ", paste(SUBJECTS_ALL, collapse = ", "))
}
keep_tarball <- !identical(Sys.getenv("KEEP_TARBALL"), "0")

## ---- Digests -------------------------------------------------------------
digest_of <- function(path, algo) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(path, algo = algo, file = TRUE))
  }
  tool <- switch(algo, md5 = c("md5", "-q"), sha256 = c("shasum", "-a", "256"),
                 stop("unsupported algo ", algo))
  out <- system2(tool[1L], c(tool[-1L], shQuote(path)), stdout = TRUE)
  sub("\\s.*$", "", trimws(out[1L]))
}
md5_of <- function(path) digest_of(path, "md5")
sha256_of <- function(path) digest_of(path, "sha256")

## ---- Download ------------------------------------------------------------
# curl resumes a partial file with -C -; download.file() cannot, so it is the
# fallback and starts over. Either way the digest is what decides.
have_curl <- nzchar(Sys.which("curl"))

fetch <- function(url, dest) {
  part <- paste0(dest, ".part")
  ok <- FALSE
  if (have_curl) {
    status <- suppressWarnings(system2(
      "curl", c("-fsSL", "-C", "-", "--retry", "5", "--retry-delay", "5",
                "--max-time", "3600", "-o", shQuote(part), shQuote(url))))
    ok <- identical(as.integer(status), 0L)
    if (!ok) message("  curl exit status ", status)
  } else {
    ok <- tryCatch({
      utils::download.file(url, destfile = part, mode = "wb", quiet = FALSE,
                           method = "libcurl")
      TRUE
    }, error = function(e) { message("  failed: ", conditionMessage(e)); FALSE })
  }
  if (ok && file.exists(part) && file.size(part) > 1e8) {
    file.rename(part, dest)
    return(TRUE)
  }
  FALSE
}

## ---- What a usable subject directory must contain ------------------------
# bold.nii.gz, mask4_vt.nii.gz and labels.txt are what 01-prepare-data.R reads;
# the mask8_* functional ROIs are what 07's ledger groups by. anat.nii.gz is
# not read by any script here. Everything beyond the three required files is
# reported rather than required, because the six subjects are not identical:
# see the per-subject notes the provenance records.
REQUIRED <- c("bold.nii.gz", "mask4_vt.nii.gz", "labels.txt")
OPTIONAL <- c("anat.nii.gz", "mask8_face_vt.nii.gz", "mask8_house_vt.nii.gz",
              "mask8b_face_vt.nii.gz", "mask8b_house_vt.nii.gz")

inventory <- function(subj_dir) {
  files <- list.files(subj_dir)
  list(
    required_present = REQUIRED[file.exists(file.path(subj_dir, REQUIRED))],
    required_missing = REQUIRED[!file.exists(file.path(subj_dir, REQUIRED))],
    optional_present = OPTIONAL[file.exists(file.path(subj_dir, OPTIONAL))],
    optional_missing = OPTIONAL[!file.exists(file.path(subj_dir, OPTIONAL))],
    extra = setdiff(files, c(REQUIRED, OPTIONAL)),
    n_files = length(files)
  )
}

## ---- Per subject ---------------------------------------------------------
provenance <- list()
failures <- character(0)

for (subj in subjects) {
  tarball <- sprintf("%s-%s.tar.gz", subj, STAMP)
  url <- paste0(BASE_URL, tarball)
  dest <- file.path(data_dir, tarball)
  subj_dir <- file.path(data_dir, subj)
  prov_path <- file.path(data_dir, sprintf("download-provenance-%s.rds", subj))

  message("\n== ", subj, " ==")

  unpacked_ok <- dir.exists(subj_dir) &&
    all(file.exists(file.path(subj_dir, REQUIRED)))

  # A verified tarball that has already been unpacked and then deleted is a
  # complete state, not a missing one: trust the recorded provenance.
  if (!file.exists(dest) && unpacked_ok && file.exists(prov_path)) {
    prov <- readRDS(prov_path)
    message("  already unpacked and verified (tarball removed); md5 ", prov$md5)
    provenance[[subj]] <- prov
    next
  }

  if (file.exists(dest)) {
    got <- md5_of(dest)
    message("  tarball present, md5 ", got)
    if (!identical(EXPECTED_MD5[[subj]], got)) {
      message("  MD5 MISMATCH\n    expected: ", EXPECTED_MD5[[subj]],
              "\n    found:    ", got,
              "\n  Refusing an unverified archive. Delete it to refetch.")
      failures[[subj]] <- "md5 mismatch against the published MD5SUMS"
      next
    }
    downloaded_now <- FALSE
  } else {
    message("  downloading ", url)
    t0 <- Sys.time()
    if (!fetch(url, dest)) {
      message("  DOWNLOAD FAILED for ", tarball)
      failures[[subj]] <- "download failed (network or truncated response)"
      next
    }
    message("  downloaded in ",
            round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
    got <- md5_of(dest)
    if (!identical(EXPECTED_MD5[[subj]], got)) {
      message("  MD5 MISMATCH after download\n    expected: ",
              EXPECTED_MD5[[subj]], "\n    found:    ", got)
      failures[[subj]] <- "md5 mismatch after download"
      next
    }
    downloaded_now <- TRUE
  }

  sha <- sha256_of(dest)
  pinned_sha <- unname(EXPECTED_SHA256[subj])
  if (!is.na(pinned_sha) && !identical(pinned_sha, sha)) {
    message("  SHA-256 MISMATCH against the recorded pin\n    expected: ",
            pinned_sha, "\n    found:    ", sha)
    failures[[subj]] <- "sha256 mismatch against the pinned digest"
    next
  }
  size <- file.size(dest)
  message("  md5 ok (", EXPECTED_MD5[[subj]], "); sha256 ", sha,
          "; ", format(size, big.mark = ","), " bytes")

  if (!unpacked_ok) {
    message("  unpacking")
    utils::untar(dest, exdir = data_dir)
    unpacked_ok <- dir.exists(subj_dir) &&
      all(file.exists(file.path(subj_dir, REQUIRED)))
  }
  inv <- inventory(subj_dir)
  if (length(inv$required_missing)) {
    message("  MISSING REQUIRED FILES under ", subj_dir, ": ",
            paste(inv$required_missing, collapse = ", "))
    failures[[subj]] <- paste0("missing required file(s): ",
                               paste(inv$required_missing, collapse = ", "))
    next
  }
  quirk <- if (length(inv$optional_missing)) {
    paste0("missing optional file(s): ",
           paste(inv$optional_missing, collapse = ", "))
  } else {
    ""
  }
  message("  contents ok: ", inv$n_files, " files",
          if (nzchar(quirk)) paste0("; ", quirk) else "")

  prov <- list(
    subject = subj, source = "pymvpa", url = url, tarball = tarball,
    bytes = size, md5 = EXPECTED_MD5[[subj]], md5_source = MD5SUMS_URL,
    md5_source_fetched = MD5SUMS_FETCHED, sha256 = sha,
    downloaded_now = downloaded_now,
    required_present = inv$required_present,
    optional_present = inv$optional_present,
    optional_missing = inv$optional_missing,
    extra_files = inv$extra, n_files = inv$n_files,
    quirk = quirk, when = Sys.time()
  )
  saveRDS(prov, prov_path)
  provenance[[subj]] <- prov

  if (!keep_tarball) {
    unlink(dest)
    message("  removed tarball (KEEP_TARBALL=0); digests kept in ",
            basename(prov_path))
  }
}

## ---- Summary -------------------------------------------------------------
if (length(provenance)) {
  saveRDS(provenance, file.path(data_dir, "download-provenance-all.rds"))
}

message("\n== summary ==")
for (subj in subjects) {
  if (!is.null(provenance[[subj]])) {
    p <- provenance[[subj]]
    message(sprintf("  %s  OK    %s bytes  md5 %s%s", subj,
                    format(p$bytes, big.mark = ","), p$md5,
                    if (nzchar(p$quirk)) paste0("  [", p$quirk, "]") else ""))
  } else {
    message(sprintf("  %s  FAIL  %s", subj,
                    if (!is.null(failures[[subj]])) failures[[subj]] else "unknown"))
  }
}

if (length(failures)) {
  message("\n", length(failures), " subject(s) unavailable: ",
          paste(names(failures), collapse = ", "),
          "\nThe downstream scripts skip subjects that are absent; they do not ",
          "invent them. Rerun this script to retry (partial files resume).")
} else {
  message("\nAll ", length(subjects), " requested subject(s) verified.")
}
