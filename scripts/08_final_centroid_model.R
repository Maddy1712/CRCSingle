# 08_final_centroid_model.R -- final v0.2 single-sample centroid model.
# Gene ranking: TCGA only (moderated F-test). Platform filter and panel size: tuning cohort GSE17536 only.
source("R/utils.R")
library(limma); library(mclust)
set.seed(2026)

d         <- readRDS("data/processed/ranked.rds")
cohorts   <- readRDS("data/processed/geo_cohorts.rds")
cms_calls <- readRDS("data/processed/cms_calls.rds")
preds_v1  <- readRDS("results/tables/external_predictions.rds")
y <- d$y; lv <- levels(y); ids <- names(d$R_geo); tune <- "GSE17536"; test <- setdiff(ids, tune)

metrics <- function(sub, cms) { ok <- !is.na(cms) & !is.na(sub); sub <- sub[ok]; cms <- cms[ok]
  c(ARI = adjustedRandIndex(sub, cms), S3_share = mean(sub == "S3"),
    CMS4_found = if (any(cms == "CMS4")) mean(sub[cms == "CMS4"] == "S3") else NA) }

# ---------- 1. rank genes (TCGA only) and define the platform filter ----------
fit0 <- eBayes(lmFit(d$R_tcga[d$universe, ], model.matrix(~ y)))
rank_genes <- rownames(topTable(fit0, coef = 2:3, number = Inf, sort.by = "F"))
dmed <- abs(apply(d$R_tcga[d$universe, ], 1, median) - apply(d$R_geo[[tune]][d$universe, ], 1, median))
st   <- function(thr) names(which(dmed <= thr))
cent_fn <- function(genes) sapply(lv, function(k) rowMeans(d$R_tcga[genes, y == k, drop = FALSE]))
sub_of <- function(id, genes, cent) {
  Rg <- d$R_geo[[id]]
  s  <- lv[max.col(cor(Rg[genes, , drop = FALSE], cent, method = "spearman"), ties.method = "first")]
  names(s) <- colnames(Rg); s }

# ---------- 2. grid: filter strength x panel size ----------
grid <- expand.grid(thr = c(0.05, 0.10, 0.15, 0.20), n = c(30, 50, 100, 200, 500, 1000, Inf))
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g_ok  <- rank_genes[rank_genes %in% st(grid$thr[i])]
  n_use <- if (is.infinite(grid$n[i])) length(g_ok) else min(grid$n[i], length(g_ok))
  if (n_use < 30) return(NULL)
  genes <- g_ok[1:n_use]; cent <- cent_fn(genes)
  tt <- metrics(sub_of(tune, genes, cent), cms_calls[[tune]][colnames(d$R_geo[[tune]])])
  s  <- unlist(lapply(test, function(id) unname(sub_of(id, genes, cent))))
  m  <- unlist(lapply(test, function(id) unname(cms_calls[[id]][colnames(d$R_geo[[id]])])))
  te <- metrics(s, m)
  data.frame(thr = grid$thr[i], n_genes = n_use, ARI_tune = round(unname(tt["ARI"]), 3),
             ARI_test = round(unname(te["ARI"]), 3), S3_test = round(unname(te["S3_share"]), 3),
             CMS4_found_test = round(unname(te["CMS4_found"]), 3)) }))
res <- res[!duplicated(res[, c("thr", "n_genes")]), ]
cat("\n==== Grid: filter strength x panel size ====\n"); print(res, row.names = FALSE)
write.csv(res, "results/tables/centroid_grid.csv", row.names = FALSE)

# pre-specified rule: smallest panel within 0.03 ARI of the best TUNING result
ok <- res$ARI_tune >= max(res$ARI_tune) - 0.03
ch <- res[ok, ]; ch <- ch[order(ch$n_genes, -ch$ARI_tune), ][1, ]
cat("\nChosen (tuning cohort only): filter <=", ch$thr, ", genes =", ch$n_genes, "\n")

# ---------- 3. final single-sample centroid model ----------
genes <- rank_genes[rank_genes %in% st(ch$thr)][1:ch$n_genes]
cent  <- cent_fn(genes)
cr_tr <- cor(d$R_tcga[genes, ], cent, method = "spearman")
marg  <- function(cr) { o <- t(apply(cr, 1, sort, decreasing = TRUE)); o[, 1] - o[, 2] }
thr_m <- unname(quantile(marg(cr_tr), 0.15))      # lowest 15% TCGA margins = Unclassified
model2 <- list(type = "centroid", panel = genes, centroids = cent, universe = d$universe,
               classes = lv, margin_threshold = thr_m, filter_thr = ch$thr, version = "0.2.0")
saveRDS(model2, "results/models/crc_single_model_v2.rds")
cat("Margin threshold (15th percentile in TCGA):", round(thr_m, 4), "\n")
source("R/predict_subtype_v2.R")

# ---------- 4. external evaluation ----------
preds2 <- lapply(setNames(ids, ids), function(id)
  predict_subtype(cohorts[[id]]$expr[, colnames(d$R_geo[[id]]), drop = FALSE], model2))
summ <- do.call(rbind, lapply(ids, function(id) {
  p <- preds2[[id]]; cm <- cms_calls[[id]][as.character(p$sample)]; kp <- p$call != "Unclassified"
  data.frame(cohort = id, n = nrow(p), S1 = round(mean(p$subtype == "S1"), 2),
             S2 = round(mean(p$subtype == "S2"), 2), S3 = round(mean(p$subtype == "S3"), 2),
             unclassified = round(mean(!kp), 3), ARI_v2_all = round(ari_f(p$subtype, cm), 3),
             ARI_v2_called = round(ari_f(p$subtype[kp], cm[kp]), 3),
             ARI_v1 = round(ari_f(preds_v1[[id]]$subtype, cm), 3)) }))
cat("\n==== Final model: external cohorts (", tune, "= tuning) ====\n", sep = ""); print(summ, row.names = FALSE)
write.csv(summ, "results/tables/v2_external_summary.csv", row.names = FALSE)

pool <- function(f) { s <- unlist(lapply(test, function(id) f(id)$s))
  m <- unlist(lapply(test, function(id) f(id)$m)); round(metrics(s, m), 3) }
cat("\nPooled TEST cohorts vs CMS (ARI, S3 share, share of CMS4 found as S3):\n")
print(rbind(
  V1_glmnet30 = pool(function(id) list(s = preds_v1[[id]]$subtype,
                                       m = cms_calls[[id]][as.character(preds_v1[[id]]$sample)])),
  V2_centroid = pool(function(id) list(s = preds2[[id]]$subtype,
                                       m = cms_calls[[id]][as.character(preds2[[id]]$sample)]))))

# CIT in GSE39582
ph <- cohorts$GSE39582$pheno; cc <- grep("cit.*molecular", colnames(ph), ignore.case = TRUE, value = TRUE)[1]
cit <- as.character(ph[[cc]])[match(preds2$GSE39582$sample, rownames(ph))]
v <- !is.na(cit) & !(toupper(cit) %in% c("NC", "NA", "N/A", ""))
cat("\nGSE39582 vs CIT (V2): ARI =", round(ari_f(preds2$GSE39582$subtype[v], cit[v]), 3),
    " NMI =", round(nmi(preds2$GSE39582$subtype[v], cit[v]), 3),
    " CramersV =", round(cramers_v(table(preds2$GSE39582$subtype[v], cit[v])), 3), "\n")
print(table(Predicted = preds2$GSE39582$subtype[v], CIT = cit[v]))

# held-out-gene reproducibility
heldout <- setdiff(st(0.20), genes)
tstat <- function(Rm, g, gg) sapply(lv, function(k) {
  grp <- factor(ifelse(g == k, "in", "out"), levels = c("out", "in"))
  if (min(table(grp)) < 3) return(rep(NA_real_, length(gg)))
  eBayes(lmFit(Rm[gg, ], model.matrix(~ grp)))$t[, 2] })
if (length(heldout) >= 100) {
  t_tc <- tstat(d$R_tcga, y, heldout)
  rp <- do.call(rbind, lapply(ids, function(id) {
    tc <- tstat(d$R_geo[[id]], factor(preds2[[id]]$subtype, lv), heldout)
    data.frame(cohort = id, t(round(sapply(seq_along(lv), function(j)
      cor(tc[, j], t_tc[, j], method = "spearman", use = "complete.obs")), 2))) }))
  colnames(rp)[-1] <- lv
  cat("\nReproducibility on", length(heldout), "genes not used by the model:\n"); print(rp, row.names = FALSE)
}

# ---------- 5. single-sample consistency and noise ----------
ex   <- cohorts$GSE39582$expr[d$universe, colnames(d$R_geo$GSE39582)]
base <- predict_subtype(ex, model2); idx <- sample(ncol(ex), 100)
single <- sapply(idx, function(i) predict_subtype(ex[, i, drop = FALSE], model2)$subtype)
cat("\nSingle-sample agreement (expect 1.00):", mean(single == base$subtype[idx]), "\n")
gsd <- apply(ex, 1, sd)
noise <- data.frame(sd_mult = c(0.1, 0.25, 0.5, 1), agreement = sapply(c(0.1, 0.25, 0.5, 1), function(s) {
  ex2 <- ex + matrix(rnorm(length(ex), 0, rep(gsd * s, ncol(ex))), nrow(ex))
  mean(predict_subtype(ex2, model2)$subtype == base$subtype) }))
print(noise, digits = 3)
