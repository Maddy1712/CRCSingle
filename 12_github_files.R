# 12_github_files.R -- README / citation files for the CRCsingle package repo, and an upload-ready copy.
source("R/utils.R")
pkgdir <- PKGDIR
gh     <- "https://github.com/Maddy1712/CRCsingle"
stopifnot(dir.exists(pkgdir))

writeLines(c("# CRCsingle", "",
  "Author: Madhuri Tiwari (17madhuritiwari@gmail.com)", "",
  "Single-sample molecular subtyping of colorectal cancer. Each sample is assigned to one of three subtypes (S1-S3) using within-sample gene ranks, with a margin-based confidence call. A cohort-adjusted mode and a 100-gene lite model are also available.", "",
  "## Installation", "",
  "    install.packages(\"remotes\")",
  "    remotes::install_github(\"Maddy1712/CRCsingle\")", "",
  "## Usage", "",
  "    library(CRCsingle)",
  "    x <- crc_example()                    # example log-scale expression matrix",
  "    res <- predict_subtype(x)             # genes (HGNC symbols) x samples",
  "    table(res$call)", "",
  "Input: log-scale expression, HGNC gene symbols as row names, no duplicated symbols.", "",
  "## Licence", "", "MIT"), file.path(pkgdir, "README.md"))

writeLines(c("cff-version: 1.2.0",
  "message: \"If you use this software, please cite it as below.\"",
  "title: \"CRCsingle: single-sample molecular subtyping of colorectal cancer\"",
  "version: 0.1.0", "date-released: 2026-10-21", "authors:",
  "  - family-names: Tiwari", "    given-names: Madhuri",
  "    email: 17madhuritiwari@gmail.com", "license: MIT",
  paste0("repository-code: \"", gh, "\"")), file.path(pkgdir, "CITATION.cff"))
writeLines("^CITATION\\.cff$", file.path(pkgdir, ".Rbuildignore"))
dir.create(file.path(pkgdir, "inst"), showWarnings = FALSE)
writeLines(c("bibentry(", "  bibtype = \"Manual\",",
  "  title   = \"CRCsingle: Single-Sample Molecular Subtyping of Colorectal Cancer\",",
  "  author  = person(\"Madhuri\", \"Tiwari\"),", "  year    = \"2026\",",
  "  note    = \"R package version 0.1.0\",", paste0("  url     = \"", gh, "\""), ")"),
  file.path(pkgdir, "inst", "CITATION"))

desc <- readLines(file.path(pkgdir, "DESCRIPTION"))
if (!any(grepl("^URL:", desc)))
  writeLines(c(desc, paste0("URL: ", gh), paste0("BugReports: ", gh, "/issues")),
             file.path(pkgdir, "DESCRIPTION"))

# upload-ready copy (hidden files such as .Rbuildignore included)
dest <- paste0(pkgdir, "_upload")
dir.create(dest, showWarnings = TRUE)
ok <- file.copy(list.files(pkgdir, full.names = TRUE, all.files = TRUE, no.. = TRUE),
                dest, recursive = TRUE, overwrite = TRUE)
stopifnot(all(ok))
print(list.files(dest, recursive = TRUE, all.files = TRUE))   # expect 17 files

# After upload, test the install from a fresh R session (not one where CRCsingle is already loaded):
#   remotes::install_github("Maddy1712/CRCsingle"); library(CRCsingle)
#   table(predict_subtype(crc_example())$call)
