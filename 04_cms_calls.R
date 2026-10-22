# 04_cms_calls.R -- CMScaller calls for TCGA and every GEO cohort (cached for later scripts).
source("R/utils.R")
library(edgeR); library(AnnotationDbi); library(org.Hs.eg.db)
stopifnot(requireNamespace("CMScaller", quietly = TRUE))

tcga    <- readRDS("data/processed/tcga_expr.rds")
lab     <- readRDS("data/processed/tcga_labels.rds")
d       <- readRDS("data/processed/ranked.rds")
cohorts <- readRDS("data/processed/geo_cohorts.rds")

run_cms <- function(em, rnaseq) {
  em <- em[rowSums(is.na(em)) == 0, , drop = FALSE]
  r <- NULL
  invisible(capture.output(suppressMessages(
    r <- tryCatch(CMScaller::CMScaller(em, rowNames = "symbol", RNAseq = rnaseq,
                                       FDR = 0.05, doPlot = FALSE),
                  error = function(e) NULL))))
  if (is.null(r)) {                                   # fallback: Entrez IDs
    ez <- mapIds(org.Hs.eg.db, keys = rownames(em), column = "ENTREZID",
                 keytype = "SYMBOL", multiVals = "first")
    ok <- !is.na(ez) & !duplicated(ez)
    em2 <- em[ok, , drop = FALSE]; rownames(em2) <- ez[ok]
    invisible(capture.output(suppressMessages(
      r <- CMScaller::CMScaller(em2, RNAseq = rnaseq, FDR = 0.05, doPlot = FALSE))))
  }
  setNames(as.character(r$prediction), rownames(r))
}

cms_calls <- list()
# TCGA: non-log CPM (CMScaller logs it itself)
cms_calls$TCGA <- run_cms(cpm(calcNormFactors(DGEList(tcga$counts[rownames(tcga$expr), ]))), TRUE)
for (id in names(d$R_geo))
  cms_calls[[id]] <- run_cms(cohorts[[id]]$expr[, colnames(d$R_geo[[id]]), drop = FALSE], FALSE)
saveRDS(cms_calls, "data/processed/cms_calls.rds")

# TCGA subtype vs CMS
lab$CMS <- cms_calls$TCGA[match(lab$sample, names(cms_calls$TCGA))]
print(table(lab$CMS, useNA = "ifany"))
print(table(lab$subtype, lab$CMS, useNA = "ifany"))
cat("Adjusted Rand index:", round(ari_f(as.character(lab$subtype), lab$CMS), 3), "\n")
saveRDS(lab, "data/processed/tcga_labels.rds")
write.csv(lab, "results/tables/tcga_subtype_labels.csv", row.names = FALSE)
