pkgname <- "amanida"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "amanida-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('amanida')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("amanida_read")
### * amanida_read

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: amanida_read
### Title: Import data
### Aliases: amanida_read

### ** Examples

coln <-  c("Compound Name", "P-value", "Fold-change", "N total", "References")
input_file <- getsampleDB()
datafile <- amanida_read(input_file, mode = "quan", coln, separator=";")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("amanida_read", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("amanida_report")
### * amanida_report

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: amanida_report
### Title: Report
### Aliases: amanida_report

### ** Examples

## Not run: 
##D column_id = c("Compound Name", "P-value", "Fold-change", "N total", "References")
##D input_file <- getsampleDB()
##D 
##D amanida_report(input_file, separator = ";", column_id, analysis_type = "quan", 
##D                 pvalue_cutoff = 0.05, fc_cutoff = 4, votecount_lim = 2, 
##D                 comp_inf = F)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("amanida_report", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("amanida_vote")
### * amanida_vote

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: amanida_vote
### Title: Qualitative meta-analysis
### Aliases: amanida_vote

### ** Examples

## Not run: 
##D coln = c("Compound Name", "Behaviour", "References")
##D input_file <- system.file("extdata", "dataset2.csv", package = "amanida")
##D data_votes <- amanida_read(input_file, mode = "qual", coln, separator = ";")
##D 
##D vote_result <- amanida_vote(data_votes)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("amanida_vote", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_amanida")
### * compute_amanida

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_amanida
### Title: Combine statistical results and compute vote-counting
### Aliases: compute_amanida

### ** Examples

## Not run: 
##D data("sample_data")
##D 
##D compute_amanida(sample_data)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_amanida", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("explore_plot")
### * explore_plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: explore_plot
### Title: Plot for compounds divergence in reports
### Aliases: explore_plot

### ** Examples

data("sample_data")

explore_plot(sample_data, type = "mix", counts = 1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("explore_plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("volcano_plot")
### * volcano_plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: volcano_plot
### Title: Volcano plot of combined results
### Aliases: volcano_plot

### ** Examples

## Not run: 
##D    data("sample_data")
##D    
##D    amanida_result <- compute_amanida(sample_data)
##D    volcano_plot(amanida_result)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("volcano_plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("vote_plot")
### * vote_plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: vote_plot
### Title: Bar-plot for compounds vote-counting
### Aliases: vote_plot

### ** Examples

## Not run: 
##D     data("sample_data")
##D     result <- compute_amanida(sample_data)
##D     vote_plot(result)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("vote_plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
