# R/utils.R -- shared helpers. Every script runs from the repository root.
suppressPackageStartupMessages({
  library(tidyverse); library(glmnet); library(ranger); library(caret); library(survival)
})

for (p in c("data/raw/geo", "data/processed", "results/models", "results/tables",
            "results/figures", "results/tables/paper", "results/figures/paper"))
  dir.create(p, recursive = TRUE, showWarnings = FALSE)

# Package output folder (sibling of the repository)
PKGDIR <- file.path(dirname(normalizePath(getwd(), winslash = "/")), "CRCsingle")

# Collapse duplicate gene symbols, keeping the row with the highest mean signal
collapse_by_symbol <- function(mat, symbols) {
  ok <- !is.na(symbols) & symbols != ""
  mat <- mat[ok, , drop = FALSE]; symbols <- symbols[ok]
  o <- order(rowMeans(mat), decreasing = TRUE)
  mat <- mat[o, , drop = FALSE]; symbols <- symbols[o]
  keep <- !duplicated(symbols)
  out <- mat[keep, , drop = FALSE]
  rownames(out) <- symbols[keep]
  out
}

# Single-sample rank transform: each sample is ranked on its own
rank_transform <- function(m, genes) {
  m <- m[genes, , drop = FALSE]
  apply(m, 2, function(x) rank(x, ties.method = "average") / length(x))
}

cramers_v <- function(tab) {
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  as.numeric(sqrt(chi / (sum(tab) * (min(dim(tab)) - 1))))
}

# Normalised mutual information
nmi <- function(a, b) {
  p <- table(a, b) / length(a); pa <- rowSums(p); pb <- colSums(p)
  mi <- sum(p[p > 0] * log((p / outer(pa, pb))[p > 0]))
  mi / sqrt(-sum(pa[pa > 0] * log(pa[pa > 0])) * -sum(pb[pb > 0] * log(pb[pb > 0])))
}

# ARI that ignores missing labels
ari_f <- function(a, b) { ok <- !is.na(a) & !is.na(b); mclust::adjustedRandIndex(a[ok], b[ok]) }

# Survival helpers
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
to_event <- function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(is.na(x) | x %in% c("", "na", "n/a", "nd", "--"), NA_integer_,
         ifelse(grepl("^(0|no|alive|censor|disease.free|none)", x), 0L, 1L))
}
cap <- function(df, tau = 60) { df$e <- ifelse(df$t > tau, 0L, df$e); df$t <- pmin(df$t, tau); df }

# Glmnet helpers
train_glm <- function(R, y) cv.glmnet(t(R), y, family = "multinomial", alpha = 0.5,
                                      nfolds = 5, type.measure = "class")
pred_glm  <- function(m, Rn) {
  pr <- predict(m, t(Rn), s = "lambda.1se", type = "response")
  matrix(pr, nrow = ncol(Rn), dimnames = list(colnames(Rn), dimnames(pr)[[2]]))
}
