---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# FixFRM

<!-- badges: start -->
<!-- badges: end -->

The goal of FixFRM is to fit functional regression-based model (FRM) to test gene-based association between set of SNPs and traits. This version can analyze quantitative, dichotomous, and survival time trait for population data. The idea of FRM is to treat multiple genetic variants as a realization of a stochastic process, which means that the region is modeled as a continuum of sequence data rather than discrete observations. The genome of an individual is viewed as a stochastic function that contains both linkage and LD or correlation information of the genetic markers. Details can be found in the following references:

* Fan RZ, Wang YF, Mills JL, Wilson AF, Bailey-Wilson JE, and Xiong MM (2013) Functional linear models for association analysis of quantitative traits. Genetic Epidemiology, 37:726-742.
* Fan RZ, Wang YF, Mills JL, Carter TC, Lobach I, Wilson AF, Bailey-Wilson JE, Weeks DE, and Xiong MM (2014) Generalized functional linear models for case-control association studies.
Genetic Epidemiology, 38:622-637.
* Fan RZ, Wang YF, Qi Y, Ding Y, Weeks DE, Lu ZH, Ren HB, Cook RJ, Xiong MM and Chen W (2016) Gene-based association analysis for censored traits via functional regressions.
Genetic Epidemiology, 40(2):133-143.

## Installation

This package is uploaded to GitHub, you can install the latest version of FixFRM by running the following R codes:

``` r
devtools::install_github("Vk1n9/flm")
library(FixFRM)
```

## Example

All FRM models can be divided into two types: beta-smooth-only model and fixed model, regardless of the type of traits. In FRM framework, genetic effect of a region is treated as a smooth function, which genetic variants can be treated as either smooth function or just discrete numbers. 

There are four datasets that needed for all functions: 

1. phenotype
2. genotype
3. covariates
4. variant position information

The phenotype data contains at least three columns, they are pedigree ID, person ID, and traits (for right censored data, traits are recorded in the last two columns). These columns should be named as ped, person, and trait, respectively. Below is an example of phenotype data for right censored cases.
``` {r echo = FALSE}
pheno <- read.csv("/Users/Beal/Box\ Sync/广东医科大学项目分类/教学相关/本科毕业论文指导/2024年/阮文谦/测试数据/ped.csv")
head(pheno)
```

The genotype data only contains SNP data, the number of columns is the total number of SNPs to be analyzed.
``` {r echo = FALSE}
geno <- read.csv("/Users/Beal/Box\ Sync/广东医科大学项目分类/教学相关/本科毕业论文指导/2024年/阮文谦/测试数据/geno.csv")
geno[1:6, 3:10]
```

The covariates data contrains at least three columns, they are pedigree ID, person ID, and any other covariates. The first two columns follows the same naming regulation which other covariates can be named any allowed forms.
``` {r echo = FALSE}
covariate <- read.csv("/Users/Beal/Box\ Sync/广东医科大学项目分类/教学相关/本科毕业论文指导/2024年/阮文谦/测试数据/covariate.csv")
head(covariate)
```

The position information data is the physical location information of the SNPs to be analyzed, the length of this data should be equal to the number of columns of genetic variants data.
``` {r echo = FALSE}
pos <- read.csv("/Users/Beal/Box\ Sync/广东医科大学项目分类/教学相关/本科毕业论文指导/2024年/阮文谦/测试数据/pos.csv")
head(pos)
```

For bete-smooth-only function, the codes to run quantitative model is

``` r
flm_beta_smooth_only(pheno, mode = "Additive", geno, pos, order, bbasis, covariate, base = "bspline", interaction = FALSE)
```

This will give p values of testing association by likelihood ratio test (LRT), chi-squared test, and F-distributed test, respectively, as shown below.

``` r
$LRT
[1] 0.6598117
$Chisq
[1] 0.6598117
$F
[1] 0.6601542
```
For fixed model, the codes to run quantitative model is

``` r
flm_fixed_model(pheno, mode = "Additive", geno, pos, order, bbasis, fbasis, gfasis, covariate, base = "fspline", interaction = FALSE)
```

These two functions share similar structure except that fixed model allows different set of basis function to fit genetic variant as a smooth function, so does gfasis arguement include. The testing results are similar, as shown below.

``` r
$LRT
[1] 0.6598117
$Chisq
[1] 0.6598117
$F
[1] 0.660154
```

