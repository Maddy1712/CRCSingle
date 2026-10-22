# v0.1 glmnet predictor (30-gene panel). Kept for comparison with the final centroid model.
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
  X  <- t(r[model$panel, , drop = FALSE])
  pr <- predict(model$fit, X, s = "lambda.1se", type = "response")
  p  <- matrix(pr, nrow = nrow(X), dimnames = list(rownames(X), dimnames(pr)[[2]]))
  conf <- apply(p, 1, max); sub <- colnames(p)[max.col(p)]
  sid  <- if (is.null(rownames(X))) seq_len(nrow(X)) else rownames(X)
  data.frame(sample = sid, subtype = sub, confidence = conf,
             call = ifelse(conf >= model$conf_threshold, sub, "Unclassified"),
             p, row.names = NULL)
}
