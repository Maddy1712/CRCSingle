# CRCsingle

Author: Madhuri Tiwari (17madhuritiwari@gmail.com)

Single-sample molecular subtyping of colorectal cancer. Each sample is assigned to one of three subtypes (S1-S3) using within-sample gene ranks, with a margin-based confidence call. A cohort-adjusted mode and a 100-gene lite model are also available.

## Installation

    install.packages("remotes")
    remotes::install_github("Maddy1712/CRCsingle")

## Usage

    library(CRCsingle)
    x <- crc_example()                    # example log-scale expression matrix
    res <- predict_subtype(x)             # genes (HGNC symbols) x samples
    table(res$call)

Input: log-scale expression, HGNC gene symbols as row names, no duplicated symbols.

## Licence

MIT
