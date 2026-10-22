# 05_baseline_classifier.R -- v0.1: model comparison, stability selection, 30-gene glmnet panel.
source("R/utils.R")
library(limma); library(mclust)
set.seed(2026)

d       <- readRDS("data/processed/ranked.rds")
cohorts <- readRDS("data/processed/geo_cohorts.rds")
R <- d$R_tcga[d$stable, ]; y <- d$y; lv <- levels(y)
cat("Training matrix (genes x patients):", dim(R), "\n")

# ---------- A. model functions + leak-free cross-validation ----------
train_rf <- function(R, y) ranger(x = as.data.frame(t(R)), y = y, probability = TRUE,
                                  num.trees = 300, mtry = 60)
pred_rf  <- function(m, Rn) predict(m, as.data.frame(t(Rn)))$predictions
train_cent <- function(R, y) sapply(levels(y), function(k) rowMeans(R[, y == k, drop = FALSE]))
pred_cent  <- function(m, Rn) {
  s <- cor(Rn, m, method = "spearman") * 20
  e <- exp(s - apply(s, 1, max)); e / rowSums(e)
}
train_tsp <- function(R, y) lapply(levels(y), function(k) {
  cls <- factor(ifelse(y == k, k, "rest"), levels = c("rest", k))
  switchBox::SWAP.KTSP.Train(R, cls, krange = c(3, 5, 7, 9), featureNo = 500)
})
pred_tsp <- function(m, Rn) {
  votes <- sapply(seq_along(m), function(i)
    as.integer(as.character(switchBox::SWAP.KTSP.Classify(Rn, m[[i]])) == lv[i]))
  votes <- matrix(votes, nrow = ncol(Rn), dimnames = list(colnames(Rn), lv))
  pr <- votes / pmax(1, rowSums(votes)); pr[rowSums(votes) == 0, ] <- 1 / length(lv); pr
}
cv_generic <- function(R, y, train_fun, pred_fun, k = 5, reps = 2) {
  out <- list()
  for (r in seq_len(reps)) for (te in createFolds(y, k = k)) {
    m  <- train_fun(R[, -te, drop = FALSE], y[-te])
    pr <- pred_fun(m, R[, te, drop = FALSE])
    out[[length(out) + 1]] <- data.frame(truth = y[te], pred = colnames(pr)[max.col(pr)], pr)
  }
  do.call(rbind, out)
}
score <- function(cv) {
  cm <- confusionMatrix(factor(cv$pred, lv), factor(cv$truth, lv))
  c(acc = unname(cm$overall["Accuracy"]), kappa = unname(cm$overall["Kappa"]),
    bal_acc = mean(cm$byClass[, "Balanced Accuracy"]))
}

models <- list(glmnet = list(train_glm, pred_glm), ranger = list(train_rf, pred_rf),
               centroid = list(train_cent, pred_cent))
if (requireNamespace("switchBox", quietly = TRUE)) models$kTSP <- list(train_tsp, pred_tsp)
perf <- NULL
for (nm in names(models)) {
  t0 <- Sys.time(); cat("CV:", nm, "...\n")
  r <- tryCatch(cv_generic(R, y, models[[nm]][[1]], models[[nm]][[2]]),
                error = function(e) { message(nm, " failed: ", conditionMessage(e)); NULL })
  if (!is.null(r)) perf <- rbind(perf, data.frame(model = nm, t(score(r)),
        minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)))
}
cat("\n==== Model comparison (5-fold x 2, inside TCGA) ====\n"); print(perf, digits = 3)
write.csv(perf, "results/tables/model_comparison_cv.csv", row.names = FALSE)

# ---------- B. stability selection and minimal panel ----------
cvl    <- cv.glmnet(t(R), y, family = "multinomial", alpha = 1, nfolds = 10, type.measure = "class")
lambda <- cvl$lambda.1se
X <- t(R); cnt <- setNames(numeric(ncol(X)), colnames(X)); B <- 100
for (b in seq_len(B)) {
  idx <- unlist(lapply(split(seq_along(y), y), function(i) sample(i, floor(length(i) * 0.5))))
  fit <- glmnet(X[idx, ], y[idx], family = "multinomial", alpha = 1, lambda = lambda)
  nz  <- unique(unlist(lapply(coef(fit), function(m) rownames(m)[as.vector(m != 0)])))
  nz  <- setdiff(nz, "(Intercept)"); cnt[nz] <- cnt[nz] + 1
}
freq   <- sort(cnt / B, decreasing = TRUE); ranked <- names(freq)
write.csv(data.frame(gene = ranked, freq = as.numeric(freq)),
          "results/tables/stability_frequency.csv", row.names = FALSE)
cat("\nGenes selected in >=80% of resamples:", sum(freq >= 0.8), "\n")

sizes <- c(5, 10, 15, 20, 30, 50, 100)
sweep <- do.call(rbind, lapply(sizes, function(n)
  data.frame(n_genes = n, t(score(cv_generic(R[ranked[1:n], ], y, train_glm, pred_glm))))))
cat("\n==== Panel size vs accuracy (optimistic: selection used all TCGA) ====\n"); print(sweep, digits = 3)
write.csv(sweep, "results/tables/panel_size_sweep.csv", row.names = FALSE)
png("results/figures/panel_size.png", 900, 650, res = 130)
plot(sweep$n_genes, sweep$bal_acc, type = "b", log = "x", pch = 19, ylim = c(0.5, 1),
     xlab = "Genes in panel", ylab = "Balanced accuracy (CV)"); dev.off()

n_star <- min(sweep$n_genes[sweep$bal_acc >= max(sweep$bal_acc) - 0.01])
panel  <- ranked[1:n_star]
cat("Chosen panel size:", n_star, "\n"); print(panel)

# ---------- C. confidence threshold and final v0.1 model ----------
cvp  <- cv_generic(R[panel, ], y, train_glm, pred_glm)
conf <- apply(cvp[, lv], 1, max); correct <- cvp$pred == cvp$truth
thr_tab <- do.call(rbind, lapply(seq(0.40, 0.95, 0.05), function(t) {
  k <- conf >= t; data.frame(thr = t, coverage = mean(k), accuracy = mean(correct[k])) }))
cat("\n==== Confidence: coverage vs accuracy (CV) ====\n"); print(thr_tab, digits = 3)
ok  <- thr_tab$accuracy >= 0.90 & thr_tab$coverage >= 0.70
thr <- if (any(ok)) min(thr_tab$thr[ok]) else 0.5
write.csv(thr_tab, "results/tables/confidence_threshold_table.csv", row.names = FALSE)

fit_final <- cv.glmnet(t(d$R_tcga[panel, ]), y, family = "multinomial", alpha = 0.5,
                       nfolds = 10, type.measure = "class")
model <- list(fit = fit_final, panel = panel, universe = d$universe, classes = lv,
              conf_threshold = thr, version = "0.1.0")
saveRDS(model, "results/models/crc_single_model.rds")
cat("\nConfidence threshold:", thr, "\n")

# ---------- D. external cohorts: predictions ----------
source("R/predict_subtype_v1.R")
preds <- lapply(names(d$R_geo), function(id)
  predict_subtype(cohorts[[id]]$expr[, colnames(d$R_geo[[id]]), drop = FALSE], model))
names(preds) <- names(d$R_geo)
ext <- do.call(rbind, lapply(names(preds), function(id) {
  p <- preds[[id]]
  data.frame(cohort = id, n = nrow(p),
             S1 = mean(p$subtype == "S1"), S2 = mean(p$subtype == "S2"), S3 = mean(p$subtype == "S3"),
             median_conf = median(p$confidence), unclassified = mean(p$call == "Unclassified")) }))
cat("\n==== External cohorts: predicted subtype mix and confidence ====\n"); print(ext, digits = 2)
write.csv(ext, "results/tables/external_prediction_summary.csv", row.names = FALSE)
saveRDS(preds, "results/tables/external_predictions.rds")

# ---------- E. concordance with CIT in GSE39582 ----------
ph  <- cohorts$GSE39582$pheno
col <- grep("cit.*molecular", colnames(ph), ignore.case = TRUE, value = TRUE)
if (length(col)) {
  cit <- as.character(ph[[col[1]]])[match(preds$GSE39582$sample, rownames(ph))]
  v   <- !is.na(cit) & !(toupper(cit) %in% c("NC", "NA", "N/A", ""))
  print(table(Predicted = preds$GSE39582$subtype[v], CIT = cit[v]))
  cat("ARI:", round(ari_f(preds$GSE39582$subtype[v], cit[v]), 3),
      "  Cramer's V:", round(cramers_v(table(preds$GSE39582$subtype[v], cit[v])), 3), "\n")
}

# ---------- F. held-out-gene reproducibility ----------
heldout <- setdiff(d$stable, panel)
tstat <- function(Rm, g, genes) sapply(lv, function(k) {
  grp <- factor(ifelse(g == k, "in", "out"), levels = c("out", "in"))
  if (min(table(grp)) < 3) return(rep(NA_real_, length(genes)))
  eBayes(lmFit(Rm[genes, ], model.matrix(~ grp)))$t[, 2] })
t_tcga <- tstat(d$R_tcga, y, heldout)
rep_tab <- do.call(rbind, lapply(names(preds), function(id) {
  tc <- tstat(d$R_geo[[id]], factor(preds[[id]]$subtype, lv), heldout)
  data.frame(cohort = id, subtype = lv, spearman = sapply(seq_along(lv), function(j)
    cor(tc[, j], t_tcga[, j], method = "spearman", use = "complete.obs"))) }))
cat("\n==== Reproducibility on genes NOT in the panel (Spearman vs TCGA) ====\n")
print(reshape(rep_tab, idvar = "cohort", timevar = "subtype", direction = "wide"), digits = 2)
write.csv(rep_tab, "results/tables/heldout_gene_reproducibility.csv", row.names = FALSE)

# ---------- G. possible duplicated patients across cohorts ----------
allR <- do.call(cbind, lapply(names(d$R_geo), function(id) {
  m <- d$R_geo[[id]]; colnames(m) <- paste(id, colnames(m), sep = "|"); m }))
cm  <- cor(allR[d$universe, ]); coh <- sub("\\|.*", "", colnames(cm))
best <- t(sapply(seq_len(ncol(cm)), function(i) {
  o <- cm[i, coh != coh[i]]; c(max(o), names(o)[which.max(o)]) }))
dup <- data.frame(sample = colnames(cm), cohort = coh,
                  max_other = as.numeric(best[, 1]), partner = best[, 2])
print(round(quantile(dup$max_other, c(.5, .9, .99, 1)), 3))
print(table(dup$cohort[dup$max_other > 0.98]))
write.csv(dup[dup$max_other > 0.98, ], "results/tables/possible_duplicates.csv", row.names = FALSE)

# ---------- H. single-sample consistency and noise robustness ----------
ex   <- cohorts$GSE39582$expr[d$universe, colnames(d$R_geo$GSE39582)]
base <- predict_subtype(ex, model)
idx  <- sample(ncol(ex), 100)
single <- sapply(idx, function(i) predict_subtype(ex[, i, drop = FALSE], model)$subtype)
cat("\nSingle-sample agreement (expect 1.00):", mean(single == base$subtype[idx]), "\n")
gsd <- apply(ex, 1, sd)
noise <- data.frame(sd_mult = c(0.1, 0.25, 0.5, 1), agreement = sapply(c(0.1, 0.25, 0.5, 1), function(s) {
  ex2 <- ex + matrix(rnorm(length(ex), 0, rep(gsd * s, ncol(ex))), nrow(ex))
  mean(predict_subtype(ex2, model)$subtype == base$subtype) }))
print(noise, digits = 3)
write.csv(noise, "results/tables/robust_noise.csv", row.names = FALSE)
