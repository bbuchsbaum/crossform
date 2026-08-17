#!/usr/bin/env Rscript
## Reproduce the file-level call graph of R/ and its strongly connected
## components, the numbers design/architecture.md quotes in "What remains".
##
##   Rscript benchmarks/call-graph-scc.R [R-dir]
##   SCC_FOCUS=receipt.R Rscript benchmarks/call-graph-scc.R
##
## The extraction rule is the one tests/testthat/test-architecture.R uses:
## top-level `name <- function(...)` definitions, then every symbol a file
## references that resolves to a definition in another file. `SCC_FOCUS`
## names the file whose edges into and out of the largest component are
## printed (default `receipt.R`, the file ticket C1 removed from it).
## Components are found with Tarjan.

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[[1L]] else "R"

internal_call_graph <- function(dir) {
  files <- list.files(dir, pattern = "[.]R$", full.names = TRUE)
  defs <- new.env(parent = emptyenv())
  parsed <- lapply(files, parse, keep.source = FALSE)
  names(parsed) <- basename(files)
  for (name in names(parsed)) {
    for (e in parsed[[name]]) {
      if (is.call(e) && (identical(e[[1L]], quote(`<-`)) ||
          identical(e[[1L]], quote(`=`))) && is.name(e[[2L]]) &&
          is.call(e[[3L]]) && identical(e[[3L]][[1L]], quote(`function`))) {
        assign(as.character(e[[2L]]), name, envir = defs)
      }
    }
  }
  collect <- function(x, out) {
    if (is.name(x)) {
      n <- as.character(x)
      if (nzchar(n) && exists(n, envir = defs, inherits = FALSE)) {
        assign(n, TRUE, envir = out)
      }
      return(invisible())
    }
    if (is.call(x) || is.pairlist(x)) {
      for (i in seq_along(x)) {
        el <- x[[i]]
        if (missing(el) || is.null(el)) next
        collect(el, out)
      }
    }
    invisible()
  }
  rows <- list()
  for (name in names(parsed)) {
    out <- new.env(parent = emptyenv())
    for (e in parsed[[name]]) collect(e, out)
    syms <- ls(out, all.names = TRUE)
    if (!length(syms)) next
    to <- vapply(syms, get, character(1), envir = defs)
    keep <- to != name
    if (!any(keep)) next
    rows[[length(rows) + 1L]] <- data.frame(from = name, to = unname(to[keep]),
      symbol = syms[keep], stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

## Tarjan, iterative-free (recursion depth is bounded by the file count).
tarjan_scc <- function(nodes, edges) {
  adj <- split(edges$to, factor(edges$from, levels = nodes))
  index <- setNames(rep(NA_integer_, length(nodes)), nodes)
  low <- index
  on_stack <- setNames(rep(FALSE, length(nodes)), nodes)
  stack <- character()
  counter <- 0L
  comps <- list()
  strongconnect <- function(v) {
    counter <<- counter + 1L
    index[[v]] <<- counter
    low[[v]] <<- counter
    stack <<- c(stack, v)
    on_stack[[v]] <<- TRUE
    for (w in unique(adj[[v]])) {
      if (is.na(index[[w]])) {
        strongconnect(w)
        low[[v]] <<- min(low[[v]], low[[w]])
      } else if (on_stack[[w]]) {
        low[[v]] <<- min(low[[v]], index[[w]])
      }
    }
    if (low[[v]] == index[[v]]) {
      pos <- match(v, stack)
      comp <- stack[pos:length(stack)]
      stack <<- stack[seq_len(pos - 1L)]
      on_stack[comp] <<- FALSE
      comps[[length(comps) + 1L]] <<- sort(comp)
    }
    invisible()
  }
  for (v in nodes) if (is.na(index[[v]])) strongconnect(v)
  comps
}

edges <- internal_call_graph(dir)
nodes <- sort(unique(c(edges$from, edges$to,
  basename(list.files(dir, pattern = "[.]R$")))))
comps <- tarjan_scc(nodes, edges)
sizes <- lengths(comps)
comps <- comps[order(-sizes)]
cat("files:", length(nodes), " edges:", nrow(unique(edges[c("from", "to")])),
  "\n")
cat("largest SCC size:", length(comps[[1L]]), "\n")
for (cm in comps[lengths(comps) > 1L]) {
  cat("\nSCC (", length(cm), "):\n", sep = "")
  cat(strwrap(paste(cm, collapse = "  "), width = 76, prefix = "  "), sep = "\n")
}

scc_files <- comps[[1L]]
focus <- if (nchar(Sys.getenv("SCC_FOCUS"))) Sys.getenv("SCC_FOCUS") else "receipt.R"
out <- edges[edges$from == focus & edges$to %in% scc_files, ]
cat("\n", focus, " out-edges into the largest SCC (", nrow(out), "):\n", sep = "")
if (nrow(out)) print(out[order(out$to, out$symbol), ], row.names = FALSE)
inn <- edges[edges$to == focus & edges$from %in% scc_files, ]
cat("\nin-edges to ", focus, " from the largest SCC (", nrow(inn), "):\n", sep = "")
if (nrow(inn)) print(unique(inn[order(inn$from, inn$symbol), ]), row.names = FALSE)
