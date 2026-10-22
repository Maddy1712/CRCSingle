# 09_ablation_survival_lite.R -- rank-representation ablation, survival with final model, 100-gene lite model.
source("R/utils.R"); source("R/predict_subtype_v2.R")
library(survival); library(mclust)
set.seed(2026)

d         <- readRDS("data/processed/ranked.rds")
cohorts   <- readRDS("data/processed/geo_cohorts.rds")
cms_calls <- readRDS("data/processed/cms_calls.rds")
tcga      <- readRDS("data/processed/tcga_expr.rds")
model2    <- readRDS("results/models/crc_single_model_v2.rds")
y <- d$y; lv <- levels(y); ids <- names(d$R_geo); tune <- "GSE17536"; test <- setdiff(ids, tune)
genes <- model2$panel
preds2 <- lapply(setNames(ids, ids), function(id)
  predict_subtype(cohorts[[id]]$expr[, colnames(d$R_geo[[id]]), drop = FALSE], model2))

metr <- function(sub, cms) { ok <- !is.na(cms) & !is.na(sub); sub <- sub[ok]; cms <- cms[ok]
  c(ARI = round(adjustedRandIndex(sub, cms), 3), S3_share = round(mean(sub == "S3"), 3),
    CMS4_found = round(mean(sub[cms == "CMS4"] == "S3"), 3)) }
pooled <- function(get_sub) {
  s <- unlist(lapply(test, function(id) unname(get_sub(id))))
  m <- unlist(lapply(test, function(id) unname(cms_calls[[id]][colnames(d$R_geo[[id]])])))
  metr(s, m) }
per_cohort <- function(get_sub) sapply(ids, function(id)
  metr(get_sub(id), cms_calls[[id]][colnames(d$R_geo[[id]])])["ARI"])

# ---------- A. ablation: does the single-sample rank representation work? ----------
tz       <- t(scale(t(tcga$expr[genes, ])))
cent_z   <- sapply(lv, function(k) rowMeans(tz[, y == k, drop = FALSE]))
cent_raw <- sapply(lv, function(k) rowMeans(tcga$expr[genes, y == k, drop = FALSE]))
ex_g <- function(id) cohorts[[id]]$expr[genes, colnames(d$R_geo[[id]]), drop = FALSE]
get_rank <- function(id) setNames(preds2[[id]]$subtype, preds2[[id]]$sample)
get_cohz <- function(id) { z <- t(scale(t(ex_g(id))))        # needs the whole cohort
  setNames(lv[max.col(cor(z, cent_z), ties.method = "first")], colnames(z)) }
get_raw  <- function(id) { e <- ex_g(id)                     # no normalisation
  setNames(lv[max.col(cor(e, cent_raw), ties.method = "first")], colnames(e)) }
abl <- rbind(`Rank, single-sample (model)`    = pooled(get_rank),
             `Cohort z-score (needs cohort)`  = pooled(get_cohz),
             `Raw expression (single-sample)` = pooled(get_raw))
cat("\n==== A. Ablation, pooled test cohorts vs CMS (same 500 genes) ====\n"); print(abl)
print(round(rbind(rank = per_cohort(get_rank), cohort_z = per_cohort(get_cohz), raw = per_cohort(get_raw)), 3))
write.csv(abl, "results/tables/ablation_pooled.csv")

# ---------- B. survival with the final model ----------
get_surv <- function(id, tcol, ecol) {
  ph <- cohorts[[id]]$pheno; p <- preds2[[id]]; ix <- match(p$sample, rownames(ph))
  df <- data.frame(cohort = id, sub = factor(p$subtype, lv), called = p$call != "Unclassified",
                   t = num(ph[[tcol]])[ix], e = to_event(ph[[ecol]])[ix])
  if (median(df$t, na.rm = TRUE) > 200) df$t <- df$t / 30.4375
  cap(df[!is.na(df$t) & !is.na(df$e) & df$t > 0, ]) }
os17 <- c("overall survival follow-up time:ch1", "overall_event (death from any cause):ch1")
df17 <- c("dfs_time:ch1", "dfs_event (disease free survival; cancer recurrence):ch1")
spec <- list(list("GSE39582", "OS",  "os.delay (months):ch1", "os.event:ch1"),
             list("GSE39582", "DFS", "rfs.delay:ch1", "rfs.event:ch1"),
             list("GSE17536", "OS",  os17[1], os17[2]), list("GSE17536", "DFS", df17[1], df17[2]),
             list("GSE17537", "OS",  os17[1], os17[2]), list("GSE17537", "DFS", df17[1], df17[2]),
             list("GSE33113", "DFS", "time to meta or recurrence:ch1", "meta or recurrence within 3 years:ch1"))
sv <- do.call(rbind, lapply(spec, function(s) { x <- get_surv(s[[1]], s[[3]], s[[4]]); x$ep <- s[[2]]; x }))

cox_row <- function(df, label, ep, strat = FALSE) {
  df <- droplevels(df)
  if (nrow(df) < 30 || sum(df$e) < 10 || nlevels(df$sub) < 2) return(NULL)
  f   <- if (strat) Surv(t, e) ~ sub + strata(cohort) else Surv(t, e) ~ sub
  fit <- coxph(f, data = df); sm <- summary(fit); cx <- sm$conf.int
  g <- function(nm) if (nm %in% rownames(cx)) sprintf("%.2f (%.2f-%.2f)", cx[nm, 1], cx[nm, 3], cx[nm, 4]) else NA
  data.frame(set = label, endpoint = ep, n = nrow(df), events = sum(df$e),
             p = signif(unname(sm$sctest["pvalue"]), 3), HR_S2 = g("subS2"), HR_S3 = g("subS3")) }
rows <- list()
for (id in unique(sv$cohort)) for (ep in unique(sv$ep[sv$cohort == id])) {
  x <- sv[sv$cohort == id & sv$ep == ep, ]
  rows[[length(rows) + 1]] <- cox_row(x, paste(id, "all calls", if (id == tune) "(tuning)" else ""), ep)
  rows[[length(rows) + 1]] <- cox_row(x[x$called, ], paste(id, "called only", if (id == tune) "(tuning)" else ""), ep) }
for (ep in c("OS", "DFS")) {
  rows[[length(rows) + 1]] <- cox_row(sv[sv$ep == ep & sv$cohort != tune, ], "POOLED test cohorts, strata=cohort", ep, TRUE)
  rows[[length(rows) + 1]] <- cox_row(sv[sv$ep == ep, ], "POOLED incl. tuning cohort, strata=cohort", ep, TRUE) }
surv_tab <- do.call(rbind, rows)
cat("\n==== B. Survival by predicted subtype (5-year cap; HR vs S1) ====\n"); print(surv_tab, row.names = FALSE)
write.csv(surv_tab, "results/tables/survival_v2.csv", row.names = FALSE)

# stage- and age-adjusted Cox in GSE39582 (stage coded I-II / III / IV)
ph <- cohorts$GSE39582$pheno; p <- preds2$GSE39582; ix <- match(p$sample, rownames(ph))
st <- num(ph[["tnm.stage:ch1"]])[ix]
sg <- factor(ifelse(st %in% 1:2, "I-II", ifelse(st == 3, "III", ifelse(st == 4, "IV", NA))),
             levels = c("I-II", "III", "IV"))
agecol <- grep("age.at.diagnosis", colnames(ph), value = TRUE)[1]
for (e in list(c("OS", "os.delay (months):ch1", "os.event:ch1"), c("DFS", "rfs.delay:ch1", "rfs.event:ch1"))) {
  df <- data.frame(sub = factor(p$subtype, lv), stage = sg, age = num(ph[[agecol]])[ix],
                   t = num(ph[[e[2]]])[ix], e = to_event(ph[[e[3]]])[ix])
  df <- cap(df[complete.cases(df) & df$t > 0, ])
  cat("\nGSE39582", e[1], "adjusted for stage + age (n =", nrow(df), ", events =", sum(df$e), ")\n")
  print(round(summary(coxph(Surv(t, e) ~ sub + stage + age, data = df))$conf.int[, c(1, 3, 4)], 2)) }

# ---------- C. lightweight 100-gene model ----------
lite_genes <- model2$panel[1:100]                    # top of the F-ranked list
cent_l <- sapply(lv, function(k) rowMeans(d$R_tcga[lite_genes, y == k, drop = FALSE]))
cr_l   <- cor(d$R_tcga[lite_genes, ], cent_l, method = "spearman")
thr_l  <- unname(quantile(apply(cr_l, 1, function(v) { s <- sort(v, decreasing = TRUE); s[1] - s[2] }), 0.15))
model_lite <- list(type = "centroid", panel = lite_genes, centroids = cent_l, universe = d$universe,
                   classes = lv, margin_threshold = thr_l, filter_thr = model2$filter_thr, version = "0.2.0-lite")
saveRDS(model_lite, "results/models/crc_single_model_lite100.rds")
pl <- lapply(setNames(ids, ids), function(id)
  predict_subtype(cohorts[[id]]$expr[, colnames(d$R_geo[[id]]), drop = FALSE], model_lite))
get_lite <- function(id) setNames(pl[[id]]$subtype, pl[[id]]$sample)
cat("\n==== C. Full (500 genes) vs lite (100 genes), pooled test cohorts vs CMS ====\n")
print(rbind(full_500 = pooled(get_rank), lite_100 = pooled(get_lite)))
cat("Unclassified (mean over test cohorts):  full =",
    round(mean(sapply(test, function(id) mean(preds2[[id]]$call == "Unclassified"))), 3),
    "  lite =", round(mean(sapply(test, function(id) mean(pl[[id]]$call == "Unclassified"))), 3), "\n")

# ---------- D. exports for the supplement ----------
write.csv(data.frame(rank = seq_along(genes), gene = genes, round(model2$centroids, 4)),
          "results/tables/final_model_500_genes_centroids.csv", row.names = FALSE)
write.csv(data.frame(rank = 1:100, gene = lite_genes), "results/tables/lite_model_100_genes.csv", row.names = FALSE)
