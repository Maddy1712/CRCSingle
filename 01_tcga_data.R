# 01_tcga_data.R -- TCGA COAD+READ: expression (TMM log-CPM) and overall survival.
# Inputs (put in data/raw/), downloaded from the GDC portal:
#   * *Merged_Count_Matrix*.tsv   STAR gene-count matrix (columns = GDC File IDs)
#   * gdc_sample_sheet*.tsv       GDC sample sheet for the same files
#   * clinical.cart*.tar.gz       GDC clinical cart (clinical.tsv, follow_up.tsv)
source("R/utils.R")
library(data.table); library(edgeR)

find_file <- function(pattern) {
  hits <- list.files("data/raw", pattern = pattern, full.names = TRUE,
                     recursive = TRUE, ignore.case = TRUE)
  if (!length(hits)) stop("No file matching '", pattern, "' in data/raw")
  hits[order(file.mtime(hits), decreasing = TRUE)][1]
}
f_mat   <- find_file("Merged_Count_Matrix.*\\.tsv$")
f_sheet <- find_file("sample.sheet.*\\.tsv$")
f_clin  <- find_file("clinical.cart.*(tar\\.gz|tar_gz|gz)$")

# ---- 1. one primary tumour per patient ----
ss <- read.delim(f_sheet, check.names = FALSE, stringsAsFactors = FALSE)
ss <- ss[grepl("^TCGA-", ss$`Project ID`), ]
ss$code <- sub(".*-(\\d{2})[A-Z]$", "\\1", ss$`Sample ID`)
ss$vial <- sub(".*-\\d{2}([A-Z])$", "\\1", ss$`Sample ID`)
ss <- ss[ss$code == "01", ]                           # 01 = primary solid tumour
hdr <- names(fread(f_mat, nrows = 0))
ss  <- ss[ss$`File ID` %in% hdr, ]
ss  <- ss[order(ss$`Case ID`, ss$vial), ]
ss  <- ss[!duplicated(ss$`Case ID`), ]                # 01A first
cat("Patients kept:", nrow(ss), "(621 expected)\n")

# ---- 2. read needed columns ----
mat <- fread(f_mat, select = c("gene_id", "gene_name", ss$`File ID`),
             data.table = FALSE, check.names = FALSE)
gene_id <- mat$gene_id; gene_name <- mat$gene_name
counts <- as.matrix(mat[, ss$`File ID`]); storage.mode(counts) <- "numeric"
rm(mat); invisible(gc())
keep <- !grepl("^N_", gene_id)                        # drop STAR summary rows
counts <- counts[keep, ]; gene_id <- gene_id[keep]; gene_name <- gene_name[keep]
colnames(counts) <- ss$`Case ID`

# ---- 3. protein-coding genes only ----
ens <- sub("\\..*", "", gene_id)
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  gt <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = unique(ens),
                              keytype = "ENSEMBL", column = "GENETYPE", multiVals = "first")
  pc <- !is.na(gt[ens]) & gt[ens] == "protein-coding"
} else {
  message("org.Hs.eg.db not installed: using an approximate symbol-based filter.")
  pc <- !grepl("^(MIR|LINC|SNOR|RNU|RNY|SCARNA|MT-)|-AS[0-9]$|-IT[0-9]$|^RP[0-9]+-|^[A-Z]{1,2}[0-9]{5,6}\\.[0-9]",
               gene_name)
}
cat("Protein-coding genes:", sum(pc), "\n")
counts <- collapse_by_symbol(counts[pc, ], gene_name[pc])

# ---- 4. filter, TMM, log-CPM ----
y <- DGEList(counts = counts)
y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
logcpm <- cpm(y, log = TRUE, prior.count = 2)

# ---- 5. clinical data and overall survival ----
tmp <- file.path("data/raw", "clinical_extract"); dir.create(tmp, showWarnings = FALSE)
untar(f_clin, exdir = tmp)
rd <- function(f) read.delim(file.path(tmp, f), check.names = FALSE,
                             na.strings = c("'--", "", "NA"), stringsAsFactors = FALSE)
cl <- rd("clinical.tsv"); fu <- rd("follow_up.tsv")
cl <- cl[cl$`cases.submitter_id` %in% colnames(logcpm), ]
fu <- fu[fu$`cases.submitter_id` %in% colnames(logcpm), ]

mx <- function(x) { x <- num(x); if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE) }
fs <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[1] else NA }

clin <- cl %>%
  group_by(case = `cases.submitter_id`) %>%
  summarise(project = fs(`project.project_id`),
            vital   = fs(`demographic.vital_status`),
            d_death = mx(`demographic.days_to_death`),
            d_last  = mx(`diagnoses.days_to_last_follow_up`),
            age     = mx(`demographic.age_at_index`),
            sex     = fs(`demographic.sex_at_birth`),
            stage   = fs(`diagnoses.ajcc_pathologic_stage`),
            .groups = "drop") %>% as.data.frame()

fu_max <- tapply(num(fu$`follow_ups.days_to_follow_up`), fu$`cases.submitter_id`,
                 function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
clin$d_fu      <- as.numeric(fu_max[clin$case])
clin$last_seen <- pmax(clin$d_last, clin$d_fu, na.rm = TRUE)
clin$os_event  <- ifelse(clin$vital == "Dead", 1L, ifelse(clin$vital == "Alive", 0L, NA))
clin$os_days   <- ifelse(clin$os_event == 1, clin$d_death, clin$last_seen)
clin$os_months <- clin$os_days / 30.4375
clin <- clin[match(colnames(logcpm), clin$case), ]
rownames(clin) <- clin$case

# ---- 6. save and report ----
saveRDS(list(expr = logcpm, counts = counts, clin = clin, samples = ss),
        "data/processed/tcga_expr.rds")
cat("\n==== SUMMARY ====\n")
print(dim(logcpm))
print(table(clin$project, useNA = "ifany"))
print(table(clin$os_event, useNA = "ifany"))
print(summary(clin$os_months))
print(table(clin$stage, useNA = "ifany"))
