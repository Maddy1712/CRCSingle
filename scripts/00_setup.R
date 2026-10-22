# 00_setup.R -- install dependencies. Run once, from the repository root.
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
# Windows only, if downloads time out:  options(download.file.method = "wininet")

cran <- c("tidyverse", "data.table", "glmnet", "ranger", "caret", "pROC", "survival",
          "survminer", "mclust", "remotes", "testthat", "BiocManager")
need <- setdiff(cran, rownames(installed.packages()))
if (length(need)) install.packages(need)

bioc <- c("TCGAbiolinks", "SummarizedExperiment", "GEOquery", "limma", "edgeR", "sva",
          "ConsensusClusterPlus", "switchBox", "GSVA", "Biobase",
          "org.Hs.eg.db", "hgu133plus2.db", "AnnotationDbi")
need <- setdiff(bioc, rownames(installed.packages()))
if (length(need)) BiocManager::install(need, update = FALSE, ask = FALSE)

# CMS classifier (GitHub)
if (!requireNamespace("CMScaller", quietly = TRUE))
  remotes::install_github("peterawe/CMScaller")

for (p in c("data/raw", "data/raw/geo", "data/processed", "results/models",
            "results/tables", "results/figures"))
  dir.create(p, recursive = TRUE, showWarnings = FALSE)

set.seed(2026)
sessionInfo()
