# 02_tcga_subtypes.R -- consensus clustering (k = 3), clinical associations, markers.
source("R/utils.R")
library(ConsensusClusterPlus); library(limma)
set.seed(2026)

tcga <- readRDS("data/processed/tcga_expr.rds")
expr <- tcga$expr; clin <- tcga$clin

# ---- 1. consensus clustering (k = 3 pre-specified; k = 2..6 reported) ----
mads <- apply(expr, 1, mad)
top  <- names(sort(mads, decreasing = TRUE))[1:2000]
mat  <- sweep(expr[top, ], 1, apply(expr[top, ], 1, median))     # median-centre genes
cc <- ConsensusClusterPlus(mat, maxK = 6, reps = 1000, pItem = 0.8, pFeature = 0.8,
                           clusterAlg = "hc", distance = "pearson", seed = 2026,
                           plot = "png", title = "results/figures/cc")

# stability: PAC = share of ambiguous consensus values (lower = cleaner)
kdf <- data.frame(
  k   = 2:6,
  PAC = sapply(2:6, function(k) { cm <- cc[[k]]$consensusMatrix
    v <- cm[lower.tri(cm)]; round(mean(v > 0.1 & v < 0.9), 3) }),
  smallest_cluster = sapply(2:6, function(k) min(table(cc[[k]]$consensusClass))))
print(kdf); write.csv(kdf, "results/tables/k_comparison.csv", row.names = FALSE)

# name clusters S1, S2, ... from largest to smallest
relabel <- function(cls) {
  cls <- unname(cls); ord <- names(sort(table(cls), decreasing = TRUE))
  factor(paste0("S", match(cls, ord)), levels = paste0("S", seq_along(ord)))
}
member_conf <- function(k, cls) {            # mean consensus with own cluster members
  cm <- cc[[k]]$consensusMatrix
  sapply(seq_along(cls), function(i) { s <- which(cls == cls[i]); s <- s[s != i]; mean(cm[i, s]) })
}
sub3 <- relabel(cc[[3]]$consensusClass)
sub4 <- relabel(cc[[4]]$consensusClass)
lab <- data.frame(sample = colnames(mat), subtype = sub3,
                  member_conf = round(member_conf(3, sub3), 3),
                  subtype_k4 = sub4, stringsAsFactors = FALSE)
cat("\nSubtype sizes (k = 3):\n"); print(table(lab$subtype))
cat("\nMembership confidence by subtype:\n"); print(tapply(lab$member_conf, lab$subtype, summary))

# ---- 2. clinical associations (descriptive) ----
d <- cbind(lab, clin[lab$sample, c("project", "sex", "age", "stage", "os_months", "os_event")])
d$stage_grp <- sub("^Stage ([IV]+).*", "\\1", d$stage)
set.seed(1)
assoc <- c(
  project = chisq.test(table(d$subtype, d$project),   simulate.p.value = TRUE, B = 10000)$p.value,
  sex     = chisq.test(table(d$subtype, d$sex),       simulate.p.value = TRUE, B = 10000)$p.value,
  stage   = chisq.test(table(d$subtype, d$stage_grp), simulate.p.value = TRUE, B = 10000)$p.value,
  age     = kruskal.test(age ~ subtype, data = d)$p.value)
cat("\nAssociation p-values:\n"); print(signif(assoc, 3))
print(table(d$subtype, d$stage_grp, useNA = "ifany"))

# ---- 3. overall survival by subtype (5-year cap) ----
sv <- d[!is.na(d$os_months) & !is.na(d$os_event) & d$os_months > 0, ]
sv$t <- pmin(sv$os_months, 60); sv$e <- ifelse(sv$os_months > 60, 0, sv$os_event)
lr <- survdiff(Surv(t, e) ~ subtype, data = sv)
cat("\nLog-rank p =", signif(1 - pchisq(lr$chisq, length(lr$n) - 1), 3),
    " (n =", nrow(sv), ", events =", sum(sv$e), ")\n")
print(summary(coxph(Surv(t, e) ~ subtype, data = sv))$conf.int)
png("results/figures/km_subtypes.png", 1000, 750, res = 130)
plot(survfit(Surv(t, e) ~ subtype, data = sv), col = 1:3, lwd = 2,
     xlab = "Months", ylab = "Overall survival", main = "TCGA COAD+READ by subtype")
legend("bottomleft", levels(sv$subtype), col = 1:3, lwd = 2, bty = "n")
dev.off()

# ---- 4. marker genes (one vs rest) and heatmap ----
mk <- do.call(rbind, lapply(levels(d$subtype), function(k) {
  g   <- factor(ifelse(d$subtype == k, "in", "out"), levels = c("out", "in"))
  fit <- eBayes(lmFit(expr, model.matrix(~ g)))
  tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "t")
  tt$gene <- rownames(tt); tt$subtype <- k
  head(tt[tt$logFC > 0, ], 50)
}))
write.csv(mk, "results/tables/subtype_markers_top50.csv", row.names = FALSE)

hm_genes <- unlist(lapply(levels(d$subtype), function(k) head(mk$gene[mk$subtype == k], 10)))
o <- order(d$subtype)
z <- t(scale(t(expr[hm_genes, o]))); z[z > 3] <- 3; z[z < -3] <- -3
png("results/figures/subtype_heatmap.png", 1400, 900, res = 130)
heatmap(z, Rowv = NA, Colv = NA, scale = "none", labCol = NA, margins = c(2, 7),
        ColSideColors = c("#E41A1C", "#377EB8", "#4DAF4A")[as.integer(d$subtype[o])],
        col = colorRampPalette(c("navy", "white", "firebrick"))(60))
dev.off()

# ---- 5. save labels ----
saveRDS(lab, "data/processed/tcga_labels.rds")
write.csv(lab, "results/tables/tcga_subtype_labels.csv", row.names = FALSE)
