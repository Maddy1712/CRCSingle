# 07_platform_robustness.R -- diagnose S3 under-calling of v0.1, compare variants with stricter platform filters.
# The tuning cohort (GSE17536) alone is used to choose the variant; other cohorts are test cohorts.
source("R/utils.R"); source("R/predict_subtype_v1.R")
library(limma); library(mclust)
set.seed(2026)

d         <- readRDS("data/processed/ranked.rds")
cohorts   <- readRDS("data/processed/geo_cohorts.rds")
model     <- readRDS("results/models/crc_single_model.rds")
cms_calls <- readRDS("data/processed/cms_calls.rds")      # from 04_cms_calls.R
y <- d$y; lv <- levels(y); ids <- names(d$R_geo)
tune <- "GSE17536"; test <- setdiff(ids, tune)

# ---------- A. S3 calibration of the current model (V1) ----------
diag_tab <- do.call(rbind, lapply(ids, function(id) {
  Rg <- d$R_geo[[id]]; p <- pred_glm(model$fit, Rg[model$panel, ])
  sub <- colnames(p)[max.col(p, ties.method = "first")]; cm <- cms_calls[[id]][colnames(Rg)]
  i4 <- which(cm == "CMS4")
  data.frame(cohort = id, S3_called = mean(sub == "S3"), CMS4_share = mean(cm == "CMS4", na.rm = TRUE),
             S3_is_CMS4 = mean(cm[sub == "S3"] == "CMS4", na.rm = TRUE),
             CMS4_found_by_S3 = mean(sub[i4] == "S3"), mean_prob_S3 = mean(p[, "S3"])) }))
cat("\n==== A. S3 calibration of the current model ====\n"); print(diag_tab, digits = 2, row.names = FALSE)

pm <- sapply(c(list(TCGA = d$R_tcga), d$R_geo), function(m) apply(m[model$panel, ], 1, median))
shift <- apply(pm[, -1], 1, median) - pm[, "TCGA"]
cat("\nPanel genes: median rank shift, arrays minus TCGA (largest first)\n")
print(round(shift[order(-abs(shift))], 2)[1:12])

# ---------- B. variants with stricter platform filters ----------
train_cent <- function(R, y) sapply(levels(y), function(k) rowMeans(R[, y == k, drop = FALSE]))
pred_cent  <- function(m, Rn) { s <- cor(Rn, m, method = "spearman") * 20
  e <- exp(s - apply(s, 1, max)); e / rowSums(e) }
train_tsp  <- function(R, y) lapply(levels(y), function(k) {
  cls <- factor(ifelse(y == k, k, "rest"), levels = c("rest", k))
  switchBox::SWAP.KTSP.Train(R, cls, krange = c(3, 5, 7, 9), featureNo = 500) })
pred_tsp   <- function(m, Rn) {
  v <- sapply(seq_along(m), function(i)
    as.integer(as.character(switchBox::SWAP.KTSP.Classify(Rn, m[[i]])) == lv[i]))
  v <- matrix(v, nrow = ncol(Rn), dimnames = list(colnames(Rn), lv))
  pr <- v / pmax(1, rowSums(v)); pr[rowSums(v) == 0, ] <- 1 / length(lv); pr }
stab_panel <- function(R, y, n_panel = 30, B = 100) {
  X <- t(R)
  lam <- cv.glmnet(X, y, family = "multinomial", alpha = 1, nfolds = 10, type.measure = "class")$lambda.1se
  cnt <- setNames(numeric(ncol(X)), colnames(X))
  for (b in seq_len(B)) {
    idx <- unlist(lapply(split(seq_along(y), y), function(i) sample(i, floor(length(i) * 0.5))))
    fit <- glmnet(X[idx, ], y[idx], family = "multinomial", alpha = 1, lambda = lam)
    nz  <- setdiff(unique(unlist(lapply(coef(fit), function(m) rownames(m)[as.vector(m != 0)]))), "(Intercept)")
    cnt[nz] <- cnt[nz] + 1 }
  names(sort(cnt, decreasing = TRUE))[1:n_panel] }

# platform filter: median-rank difference between TCGA and the TUNING cohort only
dmed <- abs(apply(d$R_tcga[d$universe, ], 1, median) - apply(d$R_geo[[tune]][d$universe, ], 1, median))
st <- function(thr) names(which(dmed <= thr))
cat("\nGenes passing platform filter:  <=0.20:", length(st(.20)), " <=0.10:", length(st(.10)),
    " <=0.05:", length(st(.05)), "\n")

vars <- list(); keep_models <- list()
vars[["V1 glmnet-30, filter 0.20 (current)"]] <- function(Rg) {
  p <- pred_glm(model$fit, Rg[model$panel, ]); colnames(p)[max.col(p, ties.method = "first")] }
for (thr in c(0.10, 0.05)) if (length(st(thr)) >= 300) local({
  thr <- thr; g <- st(thr)
  pan <- stab_panel(d$R_tcga[g, ], y, 30, 100)
  fit <- train_glm(d$R_tcga[pan, ], y)
  keep_models[[paste0("filter", thr)]] <<- list(fit = fit, panel = pan, universe = d$universe,
                                                classes = lv, conf_threshold = 0.5, version = "0.2.0")
  vars[[sprintf("V2 glmnet-30, filter %.2f", thr)]] <<- function(Rg) {
    p <- pred_glm(fit, Rg[pan, ]); colnames(p)[max.col(p, ties.method = "first")] }
})
g10 <- st(0.10)
if (length(g10) >= 300) {
  cent <- train_cent(d$R_tcga[g10, ], y)
  vars[["V3 centroid, all genes, filter 0.10"]] <- function(Rg) {
    p <- pred_cent(cent, Rg[g10, ]); colnames(p)[max.col(p, ties.method = "first")] }
  if (requireNamespace("switchBox", quietly = TRUE)) {
    tsp <- train_tsp(d$R_tcga[g10, ], y)
    vars[["V4 k-TSP, filter 0.10"]] <- function(Rg) {
      p <- pred_tsp(tsp, Rg[g10, ]); colnames(p)[max.col(p, ties.method = "first")] }
  }
}

# ---------- C. evaluate every variant against CMS ----------
metrics <- function(sub, cms) {
  ok <- !is.na(cms); sub <- sub[ok]; cms <- cms[ok]
  c(ARI = adjustedRandIndex(sub, cms), NMI = nmi(sub, cms), S3_share = mean(sub == "S3"),
    CMS4_share = mean(cms == "CMS4"),
    S3_is_CMS4 = if (any(sub == "S3")) mean(cms[sub == "S3"] == "CMS4") else NA,
    CMS4_found = if (any(cms == "CMS4")) mean(sub[cms == "CMS4"] == "S3") else NA) }
calls <- lapply(vars, function(f) lapply(setNames(ids, ids), function(id) {
  s <- f(d$R_geo[[id]]); names(s) <- colnames(d$R_geo[[id]]); s }))
get <- function(v, set) {
  s <- unlist(lapply(set, function(id) unname(calls[[v]][[id]])))
  m <- unlist(lapply(set, function(id) unname(cms_calls[[id]][names(calls[[v]][[id]])])))
  list(s = s, m = m) }
tab <- function(set) do.call(rbind, lapply(names(vars), function(v) {
  z <- get(v, set); data.frame(variant = v, t(round(metrics(z$s, z$m), 3))) }))
cat("\n==== C1. TUNING cohort (", tune, ") vs CMS ====\n", sep = "")
t_tune <- tab(tune); print(t_tune, row.names = FALSE)
cat("\n==== C2. TEST cohorts pooled (", length(test), " cohorts) vs CMS ====\n", sep = "")
t_test <- tab(test); print(t_test, row.names = FALSE)
ari <- sapply(ids, function(id) sapply(names(vars), function(v)
  round(metrics(calls[[v]][[id]], cms_calls[[id]][names(calls[[v]][[id]])])["ARI"], 3)))
cat("\n==== C3. ARI vs CMS per cohort ====\n"); print(ari)
cat("\nSelected on the tuning cohort (highest ARI):", t_tune$variant[which.max(t_tune$ARI)], "\n")
write.csv(t_test, "results/tables/variants_test_pooled.csv", row.names = FALSE)
write.csv(t_tune, "results/tables/variants_tuning.csv", row.names = FALSE)
saveRDS(keep_models, "results/models/variant_models.rds")
