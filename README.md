
# amanida: a R package for meta-analysis with non-integral data

[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/) 

## Description

`Amanida` package contains a collection of functions for computing a meta-analysis in R only using significance and effect size. It covers the lack of data provided on metabolomic studies, where is rare to have error or variance disclosed. With this adaptation, only using p-value and fold-change, global significance and effect size for compounds or metabolites are obtained. 

Furthermore, `Amanida` also computes quanlitative meta-analysis performing a vote-counting for compounds, including the option of only using identifier and trend labels.  


## Documentation

The following computations are included:

* Trend division: values are divided in tables to mantain the trend, up-regulated (fold-change > 1) and down-regulated (fold-change < 1).
* P-value combination: Fisher's method weigthed by number of participants on the study. 
* Fold-change combination: logarithmic transformation for average with weigthing by number of participants. 
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


**3. Perform adapted meta-analysis: `compute_amanida`**

```r
amanida_result <- compute_amanida(datafile)
```

In this step you will obtain an S4 object with two tables:

* adapted meta-analysis acces by `amanida_result@stat`
* vote-counting acces by `amanida_results@vote`

**4. Perform quanlitative meta-analysis: `amanida_vote`**


```r
coln = c("Compound Name", "Behaviour", "References")
input_file <- system.file("extdata", "dataset2.csv", package = "amanida")
data_votes <- amanida_read(input_file, mode = "qual", coln, separator = ";")

vote_result <- amanida_vote(data_votes)
```

In this step you will obtain an S4 object with one table:

* vote-counting acces by `vote_results@vote`

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


## Examples

There is an example dataset installed, to run examples please load:

```r
data("sample_data")
```

The dataset consits in a short list of compounds extracted from *Comprehensive volatilome and metabolome signatures of colorectal cancer in urine: A systematic review and meta-analysis.* Mallafré et al. 2021 Article in revision.


Please fill an issue if you have any question or problem :)

