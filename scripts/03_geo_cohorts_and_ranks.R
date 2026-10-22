# 03_geo_cohorts_and_ranks.R -- load GPL570 cohorts, build single-sample ranks.
# Inputs: GEO series matrix files in data/raw/geo/ (GSE39582, GSE17536, GSE17537,
#         GSE14333, GSE33113, GSE37892, GSE35896), e.g. GSE39582_series_matrix.txt
source("R/utils.R")
library(GEOquery); library(Biobase); library(AnnotationDbi); library(hgu133plus2.db)

ids   <- c("GSE39582", "GSE17536", "GSE17537", "GSE14333", "GSE33113", "GSE37892", "GSE35896")
found <- sapply(ids, function(id) {
  h <- list.files("data/raw/geo", pattern = paste0("^", id, ".*series_matrix"),
                  full.names = TRUE, ignore.case = TRUE)
  if (length(h)) h[1] else NA_character_
})
print(found)
if (all(is.na(found))) stop("No series matrix files found in data/raw/geo")

# ---- 1. load each cohort, map probes to gene symbols ----
load_cohort <- function(f) {
  g <- getGEO(filename = f, getGPL = FALSE)
  plat <- annotation(g)
  if (!identical(plat, "GPL570")) return(list(skip = paste("platform", plat)))
  ex <- exprs(g)
  if (max(ex, na.rm = TRUE) > 50) ex <- log2(pmax(ex, 1))
  sym <- mapIds(hgu133plus2.db, keys = rownames(ex), column = "SYMBOL",
                keytype = "PROBEID", multiVals = "first")
  list(expr = collapse_by_symbol(ex, as.character(sym[rownames(ex)])),
       pheno = pData(g), platform = plat)
}
cohorts <- list()
for (id in ids[!is.na(found[ids])]) {
  r <- tryCatch(load_cohort(found[[id]]), error = function(e) list(skip = conditionMessage(e)))
  if (!is.null(r$skip)) { message(id, " skipped: ", r$skip); next }
  cohorts[[id]] <- r
  cat(sprintf("%-9s %6d genes x %4d samples   NA fraction %.3f\n", id,
              nrow(r$expr), ncol(r$expr), mean(is.na(r$expr))))
}
saveRDS(cohorts, "data/processed/geo_cohorts.rds")

# ---- 2. gene universe and single-sample ranks ----
tcga <- readRDS("data/processed/tcga_expr.rds")
lab  <- readRDS("data/processed/tcga_labels.rds")
common <- Reduce(intersect, c(list(rownames(tcga$expr)), lapply(cohorts, function(x) rownames(x$expr))))
cat("Genes shared by TCGA and all cohorts:", length(common), "\n")
mad_t    <- apply(tcga$expr[common, ], 1, mad)
universe <- names(sort(mad_t, decreasing = TRUE))[1:min(5000, length(common))]

R_tcga <- rank_transform(tcga$expr, universe)
R_geo  <- lapply(cohorts, function(x) {
  m <- x$expr[universe, , drop = FALSE]
  m <- m[, colSums(is.na(m)) == 0, drop = FALSE]      # drop samples with missing genes
  rank_transform(m, universe)
})
y <- factor(lab$subtype[match(colnames(R_tcga), lab$sample)])

# platform-stable genes: median rank, TCGA vs the TUNING cohort (GSE17536) only
stable <- universe
if ("GSE17536" %in% names(R_geo)) {
  d_med  <- abs(apply(R_tcga, 1, median) - apply(R_geo[["GSE17536"]], 1, median))
  stable <- names(which(d_med <= 0.20))
}
cat("Platform-stable genes:", length(stable), "of", length(universe), "\n")

# ---- 3. quick sanity check: nearest centroid in rank space ----
cent <- sapply(levels(y), function(k) rowMeans(R_tcga[stable, y == k, drop = FALSE]))
cat("\nTCGA subtype proportions:\n"); print(round(prop.table(table(y)), 2))
for (id in names(R_geo)) {
  cr  <- cor(R_geo[[id]][stable, ], cent, method = "spearman")
  top <- levels(y)[max.col(cr)]
  cat(sprintf("\n%s (n = %d): predicted proportions\n", id, ncol(R_geo[[id]])))
  print(round(prop.table(table(factor(top, levels(y)))), 2))
}

saveRDS(list(universe = universe, stable = stable, R_tcga = R_tcga,
             R_geo = R_geo, y = y, lab = lab), "data/processed/ranked.rds")
cat("\nSaved data/processed/ranked.rds\n")
