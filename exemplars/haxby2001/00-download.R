#!/usr/bin/env Rscript
# 00-download.R -- fetch and checksum the Haxby 2001 subject 1 tarball.
#
# Idempotent: if data/subj1-2010.01.14.tar.gz already exists with the recorded
# SHA-256, nothing is downloaded. If it exists with a different digest the
# script stops rather than silently analyzing an unknown file.
#
# Primary source: PyMVPA dataset mirror. Fallbacks are listed in SOURCES and
# tried in order; the source actually used is recorded in data/download-provenance.rds.

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  Sys.getenv("EXEMPLAR_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
data_dir <- file.path(exemplar_dir, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

TARBALL <- "subj1-2010.01.14.tar.gz"

# SHA-256 of subj1-2010.01.14.tar.gz, recorded from the first verified
# download on 2026-08-13 from data.pymvpa.org (300.2 MB / 314803244 bytes,
# Last-Modified 2010-01-15). See data/download-provenance.rds.
EXPECTED_SHA256 <- "3c14bd7fad6c869e5d9b81739e24e21bb9feb6a7abef27682bbf131b6b4bec5c"

SOURCES <- c(
  pymvpa = "http://data.pymvpa.org/datasets/haxby2001/subj1-2010.01.14.tar.gz"
)

sha256_of <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(path, algo = "sha256", file = TRUE))
  }
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("\\s.*$", "", out[1])
}

dest <- file.path(data_dir, TARBALL)

if (file.exists(dest)) {
  got <- sha256_of(dest)
  message("Existing tarball SHA-256: ", got)
  if (!identical(EXPECTED_SHA256, got) &&
      !identical(EXPECTED_SHA256, "PENDING-FIRST-DOWNLOAD")) {
    stop("Tarball digest mismatch.\n  expected: ", EXPECTED_SHA256,
         "\n  found:    ", got,
         "\nRefusing to proceed with an unverified archive. Delete it to refetch.")
  }
} else {
  ok <- FALSE
  used <- NA_character_
  for (nm in names(SOURCES)) {
    url <- SOURCES[[nm]]
    message("Downloading from ", nm, ": ", url)
    part <- paste0(dest, ".part")
    status <- tryCatch({
      utils::download.file(url, destfile = part, mode = "wb", quiet = FALSE,
                           method = "libcurl")
      0L
    }, error = function(e) {
      message("  failed: ", conditionMessage(e)); 1L
    })
    if (identical(status, 0L) && file.exists(part) && file.size(part) > 1e8) {
      file.rename(part, dest)
      ok <- TRUE
      used <- nm
      break
    }
    unlink(part)
  }
  if (!ok) {
    stop("All sources failed. Fetch ", TARBALL, " manually into ", data_dir)
  }
  got <- sha256_of(dest)
  message("Downloaded from '", used, "'; SHA-256: ", got)
  if (identical(EXPECTED_SHA256, "PENDING-FIRST-DOWNLOAD")) {
    message("Record this digest as EXPECTED_SHA256 in 00-download.R.")
  } else if (!identical(EXPECTED_SHA256, got)) {
    stop("Digest mismatch against recorded value: ", got)
  }
  saveRDS(list(source = used, url = SOURCES[[used]], sha256 = got,
               when = Sys.time()),
          file.path(data_dir, "download-provenance.rds"))
}

# Unpack (idempotent: skip if the expected subject directory is present).
subj_dir <- file.path(data_dir, "subj1")
if (!dir.exists(subj_dir)) {
  message("Unpacking ", TARBALL, " into ", data_dir)
  utils::untar(dest, exdir = data_dir)
}

expected <- c("bold.nii.gz", "mask4_vt.nii.gz", "labels.txt")
present <- file.exists(file.path(subj_dir, expected))
if (!all(present)) {
  stop("Missing expected files under ", subj_dir, ": ",
       paste(expected[!present], collapse = ", "))
}
message("OK: ", subj_dir, " contains ", paste(expected, collapse = ", "))
