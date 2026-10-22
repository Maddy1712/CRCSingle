# 10_build_package.R -- add cohort-adjusted mode, compare the four modes, build + test the CRCsingle package.
source("R/utils.R"); library(mclust)
set.seed(2026)

d          <- readRDS("data/processed/ranked.rds")
cohorts    <- readRDS("data/processed/geo_cohorts.rds")
cms_calls  <- readRDS("data/processed/cms_calls.rds")
tcga       <- readRDS("data/processed/tcga_expr.rds")
model_full <- readRDS("results/models/crc_single_model_v2.rds")
model_lite <- readRDS("results/models/crc_single_model_lite100.rds")
y <- d$y; lv <- levels(y); ids <- names(d$R_geo); tune <- "GSE17536"; test <- setdiff(ids, tune)

# ---------- 1. centroids and threshold for the cohort-adjusted mode ----------
add_z <- function(m) {
  tz <- t(scale(t(tcga$expr[m$panel, ])))                 # genes x patients, z-scored per gene
  cz <- sapply(lv, function(k) rowMeans(tz[, y == k, drop = FALSE]))
  cr <- cor(tz, cz)
  o  <- t(apply(cr, 1, sort, decreasing = TRUE))
  m$centroids_z <- cz
  m$margin_threshold_z <- unname(quantile(o[, 1] - o[, 2], 0.15))
  m }
crc_models <- list(full = add_z(model_full), lite = add_z(model_lite))
saveRDS(crc_models, "results/models/crc_models_pkg.rds")

# ---------- 2. evaluate the four modes against CMS ----------
source("R/predict_subtype_pkg.R")      # uses crc_models from the global environment
ex_c <- function(id) cohorts[[id]]$expr[, colnames(d$R_geo[[id]]), drop = FALSE]
run_mode <- function(model, adj) lapply(setNames(ids, ids), function(id)
  predict_subtype(ex_c(id), model, cohort_adjust = adj))
modes <- list(`full, single-sample` = run_mode("full", FALSE), `full, cohort-adjusted` = run_mode("full", TRUE),
              `lite, single-sample` = run_mode("lite", FALSE), `lite, cohort-adjusted` = run_mode("lite", TRUE))
eval_mode <- function(pl) {
  s  <- unlist(lapply(test, function(id) pl[[id]]$subtype))
  m  <- unlist(lapply(test, function(id) unname(cms_calls[[id]][as.character(pl[[id]]$sample)])))
  kp <- unlist(lapply(test, function(id) pl[[id]]$call != "Unclassified"))
  c(ARI_all = ari_f(s, m), ARI_called = ari_f(s[kp], m[kp]),
    CMS4_found = mean(s[which(m == "CMS4")] == "S3"), unclassified = mean(!kp)) }
tab <- do.call(rbind, lapply(modes, eval_mode))
cat("\n==== Pooled test cohorts vs CMS (single-sample ARI should be 0.421 / 0.363) ====\n")
print(round(tab, 3))
print(round(sapply(ids, function(id) sapply(modes, function(pl)
  ari_f(pl[[id]]$subtype, cms_calls[[id]][as.character(pl[[id]]$sample)]))), 3))
write.csv(round(tab, 4), "results/tables/modes_comparison.csv")

# ---------- 3. write the R package ----------
pkgdir <- PKGDIR
for (s in c("R", "man", "tests/testthat", "vignettes", "inst/extdata"))
  dir.create(file.path(pkgdir, s), recursive = TRUE, showWarnings = FALSE)
save(crc_models, file = file.path(pkgdir, "R", "sysdata.rda"), compress = "xz")
file.copy("R/predict_subtype_pkg.R", file.path(pkgdir, "R", "predict_subtype.R"), overwrite = TRUE)

set.seed(1)
ex_small <- cohorts$GSE39582$expr[d$universe, sample(ncol(cohorts$GSE39582$expr), 40)]
saveRDS(ex_small, file.path(pkgdir, "inst", "extdata", "example_expr.rds"), compress = "xz")

writeLines(c(
  "crc_panel <- function(model = c(\"full\", \"lite\")) {",
  "  model <- match.arg(model)",
  "  m <- crc_models[[model]]",
  "  data.frame(rank = seq_along(m$panel), gene = m$panel,",
  "             round(m$centroids[m$panel, , drop = FALSE], 4), row.names = NULL)",
  "}",
  "crc_example <- function() {",
  "  f <- system.file(\"extdata\", \"example_expr.rds\", package = \"CRCsingle\")",
  "  if (!nzchar(f)) stop(\"example data not found\")",
  "  readRDS(f)",
  "}"), file.path(pkgdir, "R", "utils.R"))

writeLines(c(
  "Package: CRCsingle", "Type: Package",
  "Title: Single-Sample Molecular Subtyping of Colorectal Cancer", "Version: 0.1.0",
  "Authors@R: person(\"Madhuri\", \"Tiwari\", email = \"17madhuritiwari@gmail.com\", role = c(\"aut\", \"cre\"))",
  "Description: Assigns one of three molecular subtypes (S1-S3) to individual colorectal cancer expression profiles using within-sample gene ranks, a platform-stable gene set and a margin-based confidence call. An optional cohort-adjusted mode is available for datasets with many samples.",
  "License: MIT + file LICENSE", "Encoding: UTF-8", "Depends: R (>= 4.1.0)", "Imports: stats, utils",
  "Suggests: testthat (>= 3.0.0), knitr, rmarkdown", "VignetteBuilder: knitr"),
  file.path(pkgdir, "DESCRIPTION"))
writeLines(c("YEAR: 2026", "COPYRIGHT HOLDER: Madhuri Tiwari"), file.path(pkgdir, "LICENSE"))
writeLines(c("export(predict_subtype)", "export(crc_panel)", "export(crc_example)",
             "importFrom(stats, cor)"), file.path(pkgdir, "NAMESPACE"))

writeLines(c(
  "\\name{predict_subtype}", "\\alias{predict_subtype}",
  "\\title{Predict the molecular subtype of colorectal cancer samples}",
  "\\description{Assigns subtype S1, S2 or S3 to each sample with a margin-based confidence call.}",
  "\\usage{predict_subtype(expr, model = \"full\", cohort_adjust = FALSE, min_gene_frac = 0.8)}",
  "\\arguments{",
  "\\item{expr}{numeric matrix, genes (HGNC symbols as row names) by samples, on log scale.}",
  "\\item{model}{\"full\" (500 genes), \"lite\" (100 genes) or a model object.}",
  "\\item{cohort_adjust}{if TRUE, each gene is standardised across the supplied samples (needs a cohort of several samples).}",
  "\\item{min_gene_frac}{minimum fraction of the reference gene set that must be present.}}",
  "\\value{A data frame with one row per sample: subtype, margin, call (Unclassified when the margin is below the threshold) and correlations with the subtype centroids.}",
  "\\examples{", "x <- crc_example()", "head(predict_subtype(x))", "}"),
  file.path(pkgdir, "man", "predict_subtype.Rd"))
writeLines(c("\\name{crc_panel}", "\\alias{crc_panel}", "\\title{Genes and centroids of a model}",
             "\\description{Returns the gene panel with subtype centroid values.}",
             "\\usage{crc_panel(model = c(\"full\", \"lite\"))}",
             "\\arguments{\\item{model}{\"full\" or \"lite\".}}",
             "\\value{A data frame of genes and centroid values.}",
             "\\examples{", "head(crc_panel())", "}"), file.path(pkgdir, "man", "crc_panel.Rd"))
writeLines(c("\\name{crc_example}", "\\alias{crc_example}", "\\title{Example expression matrix}",
             "\\description{A small log-scale expression matrix (40 public colorectal tumours) for examples.}",
             "\\usage{crc_example()}", "\\value{A numeric matrix, genes by samples.}",
             "\\examples{", "dim(crc_example())", "}"), file.path(pkgdir, "man", "crc_example.Rd"))

writeLines(c("library(testthat)", "library(CRCsingle)", "test_check(\"CRCsingle\")"),
           file.path(pkgdir, "tests", "testthat.R"))
writeLines(c(
  "test_that(\"single sample equals batch\", {",
  "  x <- crc_example(); b <- predict_subtype(x); s <- predict_subtype(x[, 3, drop = FALSE])",
  "  expect_equal(s$subtype, b$subtype[3])",
  "})",
  "test_that(\"output structure\", {",
  "  x <- crc_example(); b <- predict_subtype(x)",
  "  expect_true(all(c(\"sample\", \"subtype\", \"margin\", \"call\") %in% names(b)))",
  "  expect_equal(nrow(b), ncol(x))",
  "  expect_true(all(b$subtype %in% c(\"S1\", \"S2\", \"S3\")))",
  "})",
  "test_that(\"lite and cohort modes run\", {",
  "  x <- crc_example()",
  "  expect_equal(nrow(predict_subtype(x, model = \"lite\")), ncol(x))",
  "  expect_equal(nrow(predict_subtype(x, cohort_adjust = TRUE)), ncol(x))",
  "})",
  "test_that(\"bad input fails clearly\", {",
  "  x <- crc_example()",
  "  expect_error(predict_subtype(x[1:100, ]))",
  "  expect_error(predict_subtype(x[, 1:2], cohort_adjust = TRUE))",
  "})"), file.path(pkgdir, "tests", "testthat", "test-predict.R"))

fence <- strrep("`", 3)
writeLines(c(
  "---", "title: \"CRCsingle: single-sample colorectal cancer subtyping\"", "author: \"Madhuri Tiwari\"",
  "output: rmarkdown::html_vignette", "vignette: >", "  %\\VignetteIndexEntry{CRCsingle}",
  "  %\\VignetteEngine{knitr::rmarkdown}", "  %\\VignetteEncoding{UTF-8}", "---", "",
  "## Overview", "",
  "CRCsingle assigns one of three subtypes (S1, S2, S3) to each colorectal cancer sample. By default each sample is scored on its own, using within-sample gene ranks, so no cohort is needed. Samples with a low margin between the best and second-best subtype are called `Unclassified`.", "",
  "## Single-sample prediction", "",
  paste0(fence, "{r}"), "library(CRCsingle)", "x <- crc_example()", "dim(x)",
  "res <- predict_subtype(x)", "head(res)", "table(res$call)", fence, "",
  "## One sample", "",
  paste0(fence, "{r}"), "predict_subtype(x[, 1, drop = FALSE])", fence, "",
  "## Lite model and cohort-adjusted mode", "",
  "The lite model uses 100 genes. With a cohort of at least 20 samples, `cohort_adjust = TRUE` standardises genes across the cohort, which improves agreement with CMS but is no longer single-sample.", "",
  paste0(fence, "{r}"), "table(predict_subtype(x, model = \"lite\")$call)",
  "table(predict_subtype(x, cohort_adjust = TRUE)$call)", fence, "",
  "## Input requirements", "",
  "Genes as HGNC symbols in the row names, log-scale expression (log2 microarray or log-CPM), no duplicated symbols, no missing values in the reference genes."),
  file.path(pkgdir, "vignettes", "CRCsingle.Rmd"))

# ---------- 4. install and test (fresh session objects removed so the package version is used) ----------
rm(predict_subtype)
if ("CRCsingle" %in% loadedNamespaces()) unloadNamespace("CRCsingle")
install.packages(pkgdir, repos = NULL, type = "source")
library(CRCsingle)
print(environment(predict_subtype))                  # must be <environment: namespace:CRCsingle>
x <- crc_example(); print(dim(x)); print(table(predict_subtype(x)$call))
testthat::test_dir(file.path(pkgdir, "tests", "testthat"), package = "CRCsingle", load_package = "installed")
cat("Package written to", pkgdir, "\n")
