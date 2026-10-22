# v0.2 predictor: centroid model (type = "centroid") with glmnet fallback (v0.1 models).
predict_subtype <- function(expr, model, min_gene_frac = 0.8) {
  expr <- as.matrix(expr)
  if (is.null(rownames(expr))) stop("'expr' needs gene symbols as rownames")
  have <- intersect(model$universe, rownames(expr))
  if (length(have) < min_gene_frac * length(model$universe))
    stop("Too few reference genes found in 'expr'.")
  if (!all(model$panel %in% have)) stop("Panel genes missing from 'expr'.")
  r <- apply(expr[have, , drop = FALSE], 2,
             function(x) rank(x, ties.method = "average") / length(x))
  r <- matrix(r, nrow = length(have), dimnames = list(have, colnames(expr)))
  sid <- if (is.null(colnames(expr))) seq_len(ncol(expr)) else colnames(expr)
  if (identical(model$type, "centroid")) {
    cr <- cor(r[model$panel, , drop = FALSE], model$centroids[model$panel, , drop = FALSE],
              method = "spearman")
    o  <- t(apply(cr, 1, sort, decreasing = TRUE))
    margin <- o[, 1] - o[, 2]
    sub <- colnames(cr)[max.col(cr, ties.method = "first")]
    colnames(cr) <- paste0("r_", colnames(cr))
    return(data.frame(sample = sid, subtype = sub, margin = margin,
                      call = ifelse(margin >= model$margin_threshold, sub, "Unclassified"),
                      cr, row.names = NULL))
  }
  X  <- t(r[model$panel, , drop = FALSE])
  pr <- predict(model$fit, X, s = "lambda.1se", type = "response")
  p  <- matrix(pr, nrow = nrow(X), dimnames = list(rownames(X), dimnames(pr)[[2]]))
  conf <- apply(p, 1, max); sub <- colnames(p)[max.col(p)]
  data.frame(sample = sid, subtype = sub, confidence = conf,
             call = ifelse(conf >= model$conf_threshold, sub, "Unclassified"), p, row.names = NULL)
}
