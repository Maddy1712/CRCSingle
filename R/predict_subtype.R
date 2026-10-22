
predict_subtype <- function(expr, model = "full", cohort_adjust = FALSE, min_gene_frac = 0.8) {
  if (is.character(model)) {
    model <- match.arg(model, c("full", "lite"))
    model <- crc_models[[model]]
  }
  expr <- as.matrix(expr)
  if (!is.numeric(expr)) stop("'expr' must be a numeric matrix (genes x samples, log scale)")
  if (is.null(rownames(expr))) stop("'expr' needs HGNC gene symbols as row names")
  if (anyDuplicated(rownames(expr))) stop("Duplicated gene symbols in row names; collapse them first")
  if (max(expr, na.rm = TRUE) > 50) warning("Values look unlogged; supply log-scale expression")
  have <- intersect(model$universe, rownames(expr))
  if (length(have) < min_gene_frac * length(model$universe))
    stop("Too few reference genes found in 'expr'.")
  miss <- setdiff(model$panel, have)
  if (length(miss)) stop("Panel genes missing from 'expr': ", paste(utils::head(miss, 5), collapse = ", "))
  if (anyNA(expr[have, , drop = FALSE])) stop("Missing values in reference genes; remove or impute first.")
  sid <- if (is.null(colnames(expr))) paste0("sample", seq_len(ncol(expr))) else colnames(expr)
  if (cohort_adjust) {
    if (ncol(expr) < 3) stop("cohort_adjust = TRUE needs several samples from one cohort")
    if (ncol(expr) < 20) warning("cohort_adjust = TRUE works best with at least 20 samples")
    z <- t(scale(t(expr[model$panel, , drop = FALSE])))
    z[is.na(z)] <- 0
    cr  <- cor(z, model$centroids_z[model$panel, , drop = FALSE])
    thr <- model$margin_threshold_z
  } else {
    r <- apply(expr[have, , drop = FALSE], 2,
               function(x) rank(x, ties.method = "average") / length(x))
    r <- matrix(r, nrow = length(have), dimnames = list(have, NULL))
    cr  <- cor(r[model$panel, , drop = FALSE], model$centroids[model$panel, , drop = FALSE],
               method = "spearman")
    thr <- model$margin_threshold
  }
  o      <- t(apply(cr, 1, sort, decreasing = TRUE))
  margin <- o[, 1] - o[, 2]
  sub    <- colnames(cr)[max.col(cr, ties.method = "first")]
  colnames(cr) <- paste0("r_", colnames(cr))
  data.frame(sample = sid, subtype = sub, margin = round(margin, 4),
             call = ifelse(margin >= thr, sub, "Unclassified"),
             round(cr, 4), row.names = NULL, check.names = FALSE)
}

