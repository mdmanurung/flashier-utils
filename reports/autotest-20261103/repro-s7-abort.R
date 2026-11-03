#!/usr/bin/env Rscript
# Purpose: reproduce the autotest 0.2.0.3 abort on ggplot2 4.x (S7) return values.
# Run:     Rscript reports/autotest-20261103/repro-s7-abort.R
# Inputs:  none; builds a one-function package in a temporary directory
# Outputs: the autotest error on stderr

# Build a one-function package --------------------------------------------------

pkg <- file.path(tempdir(), "s7plot")
dir.create(file.path(pkg, "R"), recursive = TRUE)
dir.create(file.path(pkg, "man"))
writeLines(c(
  "Package: s7plot", "Title: Minimal Autotest Abort on S7 ggplot Objects", "Version: 0.0.1",
  "Authors@R: person(\"Repro\", \"Author\", role = c(\"aut\", \"cre\"), email = \"repro@example.org\")",
  "Description: One function that takes an integer and returns a ggplot.",
  "License: MIT", "Imports: ggplot2", "Encoding: UTF-8"
), file.path(pkg, "DESCRIPTION"))
writeLines("export(bar_plot)", file.path(pkg, "NAMESPACE"))
writeLines(c(
  "bar_plot <- function(n = 3L) {",
  "  ggplot2::ggplot(data.frame(x = seq_len(n)), ggplot2::aes(x, x)) + ggplot2::geom_col()",
  "}"
), file.path(pkg, "R", "bar_plot.R"))
writeLines(c(
  "\\name{bar_plot}", "\\alias{bar_plot}", "\\title{Bar plot}", "\\usage{bar_plot(n = 3L)}",
  "\\arguments{\\item{n}{Number of bars, an integer.}}", "\\value{A ggplot object.}",
  "\\description{Draws n bars.}", "\\examples{", "p <- bar_plot(3L)", "}"
), file.path(pkg, "man", "bar_plot.Rd"))

# Install it and run autotest ---------------------------------------------------

lib <- file.path(tempdir(), "lib")
dir.create(lib)
install.packages(pkg, lib = lib, repos = NULL, type = "source", quiet = TRUE)
# autotest runs subprocesses, which read R_LIBS rather than this session's .libPaths().
Sys.setenv(R_LIBS = paste(c(lib, .libPaths()), collapse = .Platform$path.sep))
.libPaths(c(lib, .libPaths()))
cat(R.version.string, "| autotest", format(packageVersion("autotest")), "| ggplot2", format(packageVersion("ggplot2")), "\n")
autotest::autotest_package(pkg, test = TRUE, progress = "none")
