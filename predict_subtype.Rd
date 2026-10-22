\name{predict_subtype}
\alias{predict_subtype}
\title{Predict the molecular subtype of colorectal cancer samples}
\description{Assigns subtype S1, S2 or S3 to each sample with a margin-based confidence call.}
\usage{predict_subtype(expr, model = "full", cohort_adjust = FALSE, min_gene_frac = 0.8)}
\arguments{
\item{expr}{numeric matrix, genes (HGNC symbols as row names) by samples, on log scale.}
\item{model}{"full" (500 genes), "lite" (100 genes) or a model object.}
\item{cohort_adjust}{if TRUE, each gene is standardised across the supplied samples (needs a cohort of several samples).}
\item{min_gene_frac}{minimum fraction of the reference gene set that must be present.}}
\value{A data frame with one row per sample: subtype, margin, call (Unclassified when the margin is below the threshold) and correlations with the subtype centroids.}
\examples{
x <- crc_example()
head(predict_subtype(x))
}
