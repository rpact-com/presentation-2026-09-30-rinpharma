# crmPack Workshop at R/Pharma 2026

Slides for the crmPack workshop at R/Pharma on 30 September 2026.

## Generated R code

`scripts/generate-code.R` extracts the R chunks from `parts/workshop.qmd` into
`code/workshop.R` without running them. Quarto runs it before each render, so the
Code page and the Getting Started slide offer the current script for download.
Run `Rscript scripts/generate-code.R` from the project root to regenerate it
without rendering. The script requires `knitr`.
