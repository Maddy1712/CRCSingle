crc_panel <- function(model = c("full", "lite")) {
  model <- match.arg(model)
  m <- crc_models[[model]]
  data.frame(rank = seq_along(m$panel), gene = m$panel,
             round(m$centroids[m$panel, , drop = FALSE], 4), row.names = NULL)
}
crc_example <- function() {
  f <- system.file("extdata", "example_expr.rds", package = "CRCsingle")
  if (!nzchar(f)) stop("example data not found")
  readRDS(f)
}
