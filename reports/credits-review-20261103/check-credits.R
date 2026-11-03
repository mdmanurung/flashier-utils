#!/usr/bin/env Rscript
# Purpose: check that every package the code calls is declared and credited, and
#          that NOTICE, provenance.csv and DESCRIPTION pins agree.
# Run:     scripts/run-r.sh reports/credits-review-20261103/check-credits.R reports/credits-review-20261103
# Inputs:
#   DESCRIPTION            Imports, Suggests, Remotes pins [package]
#   R/*.R                  pkg:: calls and .check_engine() declarations [package]
#   CREDITS.md             credited packages and prose credit sections [package]
#   inst/NOTICE            upstream commits named in the notice [package]
#   inst/provenance.csv    one row per audited function [scripts/generate_provenance.py]
#   inst/licenses/         bundled third-party licence texts [package]
# Outputs (in <run_dir>):
#   packages.tsv           per package: declared field, called in R/, credited
#   provenance.tsv         per provenance row: function defined in R/, source commit named
#   pins.tsv               per pinned remote: commit named in NOTICE or CREDITS
#   session_info.txt       R and package versions

run_dir <- commandArgs(trailingOnly = TRUE)[1]

# Read the package sources ------------------------------------------------------

description <- read.dcf("DESCRIPTION")[1, ]
code <- unlist(lapply(list.files("R", pattern = "[.]R$", full.names = TRUE), readLines))
credits <- readLines("CREDITS.md")
notice <- readLines("inst/NOTICE")
provenance <- utils::read.csv("inst/provenance.csv", stringsAsFactors = FALSE)

# Declared, called and credited packages ----------------------------------------

# package names from a DESCRIPTION field, without version constraints
field_packages <- function(field) trimws(sub("\\(.*", "", strsplit(description[[field]], ",")[[1]]))
imports <- field_packages("Imports")
suggests <- field_packages("Suggests")
# namespaces called as pkg::fn in executable lines (comments and strings excluded)
executable <- sub("#.*$", "", gsub('"[^"]*"', '""', code))
called <- unique(sub("::$", "", unlist(regmatches(executable, gregexpr("\\b[A-Za-z][A-Za-z0-9.]+::", executable)))))
# engines checked by name before use
engines <- unique(unlist(regmatches(code, gregexpr('(?<=\\.check_engine\\(")[A-Za-z0-9.]+', code, perl = TRUE))))
base_packages <- c("base", "stats", "utils", "methods", "parallel")
packages <- sort(unique(setdiff(c(imports, suggests, called, engines), base_packages)))
credited_heading <- sub("^### ([A-Za-z0-9.]+).*", "\\1", grep("^### ", credits, value = TRUE))
package_table <- data.frame(
  package = packages,
  declared = ifelse(packages %in% imports, "Imports", ifelse(packages %in% suggests, "Suggests", NA_character_)),
  called_in_R = packages %in% c(called, engines),
  credited_heading = packages %in% credited_heading,
  credited_text = vapply(packages, function(p) any(grepl(paste0("\\b", p, "\\b"), credits)), logical(1))
)
utils::write.table(package_table, file.path(run_dir, "packages.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Provenance rows against R/ definitions and pins -------------------------------

# commit SHAs pinned in DESCRIPTION Remotes, by repository
remotes <- trimws(strsplit(description[["Remotes"]], ",")[[1]])
pins <- stats::setNames(sub(".*@", "", remotes), sub("@.*", "", remotes))
defined <- unique(sub("^([A-Za-z0-9_.]+) <- function.*", "\\1", grep("^[A-Za-z0-9_.]+ <- function", code, value = TRUE)))
provenance_table <- data.frame(
  provenance_id = provenance$provenance_id,
  function_name = provenance$`function`,
  defined_in_R = provenance$`function` %in% defined,
  source_commit = provenance$source_commit,
  # an audited source commit must be citable from CREDITS.md or NOTICE
  commit_named = vapply(provenance$source_commit, function(sha) {
    !grepl("^[0-9a-f]{40}$", sha) || any(grepl(sha, c(credits, notice), fixed = TRUE))
  }, logical(1))
)
utils::write.table(provenance_table, file.path(run_dir, "provenance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Pinned remotes against NOTICE and CREDITS -------------------------------------

pin_table <- data.frame(
  remote = names(pins), commit = unname(pins),
  in_notice = vapply(pins, function(sha) any(grepl(sha, notice, fixed = TRUE)), logical(1)),
  in_credits = vapply(pins, function(sha) any(grepl(sha, credits, fixed = TRUE)), logical(1))
)
utils::write.table(pin_table, file.path(run_dir, "pins.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Session info ------------------------------------------------------------------

writeLines(utils::capture.output(utils::sessionInfo()), file.path(run_dir, "session_info.txt"))
