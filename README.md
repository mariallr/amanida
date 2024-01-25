
# amanida: a R package for meta-analysis with non-integral data

[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/) 

## Description

`Amanida` package contains a collection of functions for computing a meta-analysis in R only using significance and effect size. It covers the lack of data provided on metabolomic studies, where is rare to have error or variance disclosed. With this adaptation, only using p-value and fold-change, global significance and effect size for compounds or metabolites are obtained. 

Furthermore, `Amanida` also computes qualitative meta-analysis performing a vote-counting for compounds, including the option of only using identifier and trend labels.  


## Documentation

The following computations are included:

* P-value combination: Fisher's method weighted by number of participants on the study. 
* Fold-change combination: logarithmic transformation for average with weighting by number of participants. 
* Compound vote-counting: votes are +1 for up-regulation, -1 for down-regulation and 0 if no trend. The total votes are divided by the number of reports. 

The following plots are included to visualize the results: 

* Volcano plot of meta-analysis results: showing compounds labels for over the selected cut-off. 
* Bar plot of qualitative results.
* Bar plot of reports divided by trend including the total vote-counting.

## Installation

#### Beta/Github release:

Installation using R package devtools:

```r
install.packages("devtools")
devtools::install_github("mariallr/amanida")
```

#### CRAN:

```r
install.packages("amanida")
```

## Usage

You can use `Amanida` package in RStudio or R. After installation (explained before) follow this steps: 

**1. Load package in your script:**

```r
library(amanida)
```

**2. Read your data: `amanida_read`**

Supported files are csv, xls/xlsx and txt. 

For quantitative meta-analysis include the following parameters:

* Indicate mode = "quan"
* coln: vector containing the column names, which need to be in this order:
  * Id: compound name or unique identification
  * P-value
  * Fold-change
  * N: number of individuals in the study
  * Reference: bibliographic reference of the results

```r
coln = c("Compound Name", "P-value", "Fold-change", "N total", "References")
input_file <- system.file("extdata", "dataset2.csv", package = "amanida")
datafile <- amanida_read(input_file, mode = "quan", coln, separator=";")
```

For qualitative meta-analysis include the following parameters:

* Indicate mode = "qual"
* coln: vector containing the column names, which need to be in this order:
  * Id: compound name or unique identification
  * Trend: can be up-regulated or down-regulated
  * Reference: bibliographic reference of the results

```r
coln = c("Compound Name", "Behaviour", "References")
input_file <- system.file("extdata", "dataset2.csv", package = "amanida")
datafile <- amanida_read(input_file, mode = "qual", coln, separator=";")
```

Before the meta-analysis the IDs can be checked using public databases information. The IDs in format chemical name, InChI, InChIKey, and SMILES are searched in PubChem to transform all into a common nomenclature using `webchem` package. Harmonization names process is based in *Villalba H, Llambrich M, Gumà J, Brezmes J, Cumeras R. A Metabolites Merging Strategy (MMS): Harmonization to Enable Studies’ Intercomparison. Metabolites. 2023; 13(12):1167. https://doi.org/10.3390/metabo13121167*

```{r}
datafile <- check_names(datafile)
```


**3. Perform adapted meta-analysis: `compute_amanida`**

```r
amanida_result <- compute_amanida(datafile, comp.inf = F)
```

In this step you will obtain an S4 object with two tables:

* adapted meta-analysis acces by `amanida_result@stat`
* vote-counting acces by `amanida_results@vote`

Selecting the option `comp.inf = T` the package need the previous use of `check_names`. Then using PubChem ID duplicates are checked. Results are returned including the following information: PubChem ID, Molecular Formula, Molecular Weight, SMILES, InChIKey, KEGG, ChEBI, HMDB, Drugbank. 

**4. Perform qualitative meta-analysis: `amanida_vote`**


```r
coln = c("Compound Name", "Behaviour", "References")
input_file <- system.file("extdata", "dataset2.csv", package = "amanida")
data_votes <- amanida_read(input_file, mode = "qual", coln, separator = ";")

vote_result <- amanida_vote(data_votes)
```

For qualitative analysis the `check_names` can be also used, following the same procedure explained in Section 2. 

In this step you will obtain an S4 object with one table:

* vote-counting access by `vote_results@vote`

#### Plots

**Graphical visualization for adapted meta-analysis results: `volcano_plot`**

```r
volcano_plot(amanida_result, cutoff = c(0.05,4))
```

**Graphical visualization of compounds vote-counting: `vote_plot`**

Data can be subset for better visualization using counts parameter to indicate the vote-counting cut-off. 

```r
vote_plot(amanida_result)
```

**Graphical visualization of compounds vote-counting and reports divided trend: `explore_plot`**

Data can be shown in three types:
* type = "all": show all data
* type = "sub": subset the data by a cut-off value indicated by the counts parameter 
* type = "mix": subset the data by a cut-off value indicated by the counts parameter and show compounds with discrepancies (reports up-regulated and down-regulated)

```r
explore_plot(sample_data, type = "mix", counts = 1)
```

## Report 

All results using Amanida can be obtained in a single step using `amanida_report` function. It only requires the following parameters for qualitative analysis report:
* file: path to the dataset
* separator: separator used in the dataset
* analysis_type: specify "quan"
* column_id: nomes of columns to be used, see `amanida_read` documentation for more information
* pvalue_cutoff: numeric value where the p-value will be considered as significant, usually 0.05
* fc_cutoff: numeric value where the fold-change will be considered as significant, usually 2
* votecount_lim: numeric value set as minimum to show vote-counting results
* comp_inf: to include name checking and IDs retrieval.

And for quantitative analysis report:
* file: path to the dataset
* separator: separator used in the dataset
* analysis_type: specify "qual"
* column_id: nomes of columns to be used, see `amanida_read` documentation for more information
* votecount_lim: numeric value set as minimum to show vote-counting results
* comp_inf: to include name checking and IDs retrieval.

```r
column_id = c("Compound Name", "P-value", "Fold-change", "N total", "References")
input_file <- system.file("extdata", "dataset2.csv", package = "amanida")
amanida_report(input_file, 
                separator = ";", 
                column_id, 
                analysis_type = "quan", 
                pvalue_cutoff = 0.05, 
                fc_cutoff = 4, 
                votecount_lim = 2, 
                comp_inf = F)
  
```


## Examples

There is an example dataset installed, to run examples please load:

```r
data("sample_data")
```

The dataset consist in a short list of compounds extracted from *Comprehensive Volatilome and Metabolome Signatures of Colorectal Cancer in Urine: A Systematic Review and Meta-Analysis* Mallafré et al. Cancers 2021, 13(11), 2534; https://doi.org/10.3390/cancers13112534


Please fill an issue if you have any question or problem :)

