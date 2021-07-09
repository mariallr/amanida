
input_file <- system.file("extdata", "dataset2.csv", package = "amanida")

rmarkdown::render(
  input = "report.Rmd",
  output_file = "prova1.html",
  output_dir = "/Users/maria/OneDrive - URV/metaanalysis_pack",
  params = list(
    file_name = input_file,
    separator = ";",
    analysis_type = "quan",
    column_id = c("Compound Name", "P-value", "Fold-change", "N total", "References"),
    pvalue_cutoff = 0.05,
    fc_cutoff = 4,
    votecount_lim = 1,
    show_code = FALSE
  )
)
