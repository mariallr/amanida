## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(123)

## ----install,eval=FALSE-------------------------------------------------------
#  install.packages("amanida")

## -----------------------------------------------------------------------------
library(amanida)

## -----------------------------------------------------------------------------
coln = c("Compound Name", "P-value", "Fold-change", "N total", "References")
input_file <- getsampleDB()
datafile <- amanida_read(input_file, mode = "quan", coln, separator=";")

## -----------------------------------------------------------------------------
amanida_result <- compute_amanida(datafile, comp.inf = F)

## -----------------------------------------------------------------------------
head(amanida_result@stat)

## ----fig.width = 7, fig.height = 9--------------------------------------------
volcano_plot(amanida_result, cutoff = c(0.05,3.5))

## ----fig.width = 5, fig.height = 6--------------------------------------------
vote_plot(amanida_result, counts = 1)

## ----fig.width = 6.5, fig.height = 8.5----------------------------------------
explore_plot(sample_data, type = "mix", counts = 1)

## ----eval = F-----------------------------------------------------------------
#  column_id = c("Compound Name", "P-value", "Fold-change", "N total", "References")
#  input_file <- getsampleDB()
#  
#  amanida_report(input_file,
#                 separator = ";",
#                 analysis_type = "quan-qual",
#                 column_id = column_id,
#                 pvalue_cutoff = 0.05,
#                 fc_cutoff = 4,
#                 votecount_lim = 2,
#                 comp_inf = F)

## -----------------------------------------------------------------------------
sessionInfo()

