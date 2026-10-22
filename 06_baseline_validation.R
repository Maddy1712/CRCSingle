# 06_baseline_validation.R -- unbiased internal accuracy, CMS/CIT concordance, external survival (v0.1 model).
source("R/utils.R")
library(survival); library(mclust)
set.seed(2026)

d         <- readRDS("data/processed/ranked.rds")
cohorts   <- readRDS("data/processed/geo_cohorts.rds")
model     <- readRDS("results/models/crc_single_model.rds")
preds     <- readRDS("results/tables/external_predictions.rds")
cms_calls <- readRDS("data/processed/cms_calls.rds")
y <- d$y; lv <- levels(y); lab <- d$lab; panel <- model$panel
R <- d$R_tcga[d$stable, ]

# ---------- A. nested CV: panel selection INSIDE each fold ----------
RUN_NESTED <- TRUE
if (RUN_NESTED) {
  nested_fold <- function(Rtr, ytr, Rte, B = 30, n_panel = 30) {
    X <- t(Rtr)
    lam <- cv.glmnet(X, ytr, family = "multinomial", alpha = 1, nfolds = 5,
                     type.measure = "class")$lambda.1se
    cnt <- setNames(numeric(ncol(X)), colnames(X))
    for (b in seq_len(B)) {
      idx <- unlist(lapply(split(seq_along(ytr), ytr), function(i) sample(i, floor(length(i) * 0.5))))
      fit <- glmnet(X[idx, ], ytr[idx], family = "multinomial", alpha = 1, lambda = lam)
      nz  <- setdiff(unique(unlist(lapply(coef(fit), function(m) rownames(m)[as.vector(m != 0)]))), "(Intercept)")
      cnt[nz] <- cnt[nz] + 1
    }
    g <- names(sort(cnt, decreasing = TRUE))[1:n_panel]
    pred_glm(train_glm(Rtr[g, ], ytr), Rte[g, , drop = FALSE])
  }
  out <- list()
  for (r in 1:2) for (te in caret::createFolds(y, k = 5)) {
    pr <- nested_fold(R[, -te], y[-te], R[, te, drop = FALSE])
    out[[length(out) + 1]] <- data.frame(truth = y[te], pred = colnames(pr)[max.col(pr)])
  }
  nc <- do.call(rbind, out)
  cm <- caret::confusionMatrix(factor(nc$pred, lv), factor(nc$truth, lv))
  cat("\n==== A. Nested CV (panel selected inside folds) ====\n")
  print(round(c(acc = unname(cm$overall["Accuracy"]), kappa = unname(cm$overall["Kappa"]),
                bal_acc = mean(cm$byClass[, "Balanced Accuracy"])), 3))
}

# ---------- B. does the 30-gene panel reproduce a full-gene classifier externally? ----------
cent <- sapply(lv, function(k) rowMeans(R[, y == k, drop = FALSE]))
agree <- do.call(rbind, lapply(names(preds), function(id) {
  cr <- cor(d$R_geo[[id]][d$stable, ], cent, method = "spearman")
  fp <- factor(lv[max.col(cr)], lv); pp <- factor(preds[[id]]$subtype, lv)
  data.frame(cohort = id, agreement = round(mean(fp == pp), 3),
             kappa = round(unname(caret::confusionMatrix(pp, fp)$overall["Kappa"]), 3)) }))
cat("\n==== B. 30-gene panel vs full-gene centroid, external cohorts ====\n"); print(agree)

# ---------- C. CMS and CIT concordance ----------
conc <- function(a, b) { ok <- !is.na(a) & !is.na(b); tab <- table(a[ok], b[ok])
  c(n = sum(ok), ARI = round(adjustedRandIndex(a[ok], b[ok]), 3),
    NMI = round(nmi(a[ok], b[ok]), 3), CramersV = round(cramers_v(tab), 3)) }
res <- list(TCGA = conc(as.character(lab$subtype), cms_calls$TCGA[lab$sample]))
for (id in names(preds))
  res[[id]] <- conc(preds[[id]]$subtype, cms_calls[[id]][as.character(preds[[id]]$sample)])
cat("\n==== C. Concordance with CMScaller calls ====\n"); print(do.call(rbind, res))

ph <- cohorts$GSE39582$pheno
cc <- grep("cit.*molecular", colnames(ph), ignore.case = TRUE, value = TRUE)[1]
cit <- as.character(ph[[cc]])[match(preds$GSE39582$sample, rownames(ph))]
v <- !is.na(cit) & !(toupper(cit) %in% c("NC", "NA", "N/A", ""))
cat("\nGSE39582 vs CIT:", paste(names(conc(preds$GSE39582$subtype[v], cit[v])),
    conc(preds$GSE39582$subtype[v], cit[v]), sep = "=", collapse = "  "), "\n")

# ---------- D. survival in external cohorts ----------
pick <- function(ph, pats) { for (p in pats) { h <- grep(p, colnames(ph), ignore.case = TRUE, value = TRUE)
  if (length(h)) return(h[1]) }; NA_character_ }
ends <- list(
  OS  = list(t = c("^os\\.delay", "overall.*(survival|follow).*time", "overall.*time", "os.*time"),
             e = c("^os\\.event", "overall.*event", "death|dead|vital")),
  DFS = list(t = c("^rfs\\.delay", "(dfs|disease.free|relapse.free|recurrence.free).*(time|delay)",
                   "(dfs|disease.free).*(survival|follow)", "time.*(relapse|recurr)"),
             e = c("^rfs\\.event", "(dfs|disease.free|relapse|recurrence).*event", "relapse|recurr")))
cat("\n==== D. Survival columns detected (CHECK THESE) ====\n")
surv_rows <- list()
for (id in names(preds)) for (ep in names(ends)) {
  ph <- cohorts[[id]]$pheno
  tcol <- pick(ph, ends[[ep]]$t); ecol <- pick(ph, ends[[ep]]$e)
  cat(sprintf("%-9s %-3s time: %-45s event: %s\n", id, ep, tcol, ecol))
  if (is.na(tcol) || is.na(ecol)) next
  p  <- preds[[id]]; ix <- match(p$sample, rownames(ph))
  df <- data.frame(sub = factor(p$subtype, lv), t = num(ph[[tcol]])[ix], e = to_event(ph[[ecol]])[ix])
  if (median(df$t, na.rm = TRUE) > 200) df$t <- df$t / 30.4375          # days -> months
  df <- df[!is.na(df$t) & !is.na(df$e) & df$t > 0, ]
  df <- cap(df); df$sub <- droplevels(df$sub)
  if (nrow(df) < 30 || sum(df$e) < 10 || nlevels(df$sub) < 2) next
  lr <- survdiff(Surv(t, e) ~ sub, data = df)
  cx <- summary(coxph(Surv(t, e) ~ sub, data = df))$conf.int
  g  <- function(nm) if (nm %in% rownames(cx)) sprintf("%.2f (%.2f-%.2f)", cx[nm, 1], cx[nm, 3], cx[nm, 4]) else NA
  surv_rows[[length(surv_rows) + 1]] <- data.frame(cohort = id, endpoint = ep, n = nrow(df),
    events = sum(df$e), p_logrank = signif(1 - pchisq(lr$chisq, length(lr$n) - 1), 3),
    HR_S2_vs_S1 = g("subS2"), HR_S3_vs_S1 = g("subS3"))
}
sv <- do.call(rbind, surv_rows)
cat("\n==== D. External survival by subtype (5-year cap) ====\n"); print(sv, row.names = FALSE)
write.csv(sv, "results/tables/external_survival.csv", row.names = FALSE)

# ---------- E. TCGA: are subtypes tied to tumour site? ----------
tcga <- readRDS("data/processed/tcga_expr.rds")
print(round(prop.table(table(subtype = lab$subtype, project = tcga$clin[lab$sample, "project"]), 1), 2))
