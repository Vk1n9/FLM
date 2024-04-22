library(fda)
library(MASS)
library(Matrix)
#library(globaltest)
gt
### This Function is modified by Ruzong Fan, March 20, 2013 ###
#' Title
#'
#' @param pheno a data frame of phenotype data. This data frame should contain three columns: pedigree id, person id, and phenotype value. These three columns should be named as ped, person, and trait, respectively.
#' @param mode a character value specifying the discrete realization of genetic variant function. The default of “Additive” assumes additive effect of minor alleles, “Dom” assume dominant effect, and “Rec” assumes recessive effect.
#' @param geno  a matrix of genotype data only. The number of rows should be equal to the number of individuals in pheno. The number of columns is the total number of SNPs.
#' @param pos a vector of physical positions of the SNP used in geno. The length of pos should be equal to the column number of geno.
#' @param order an integer specifying the order of b-splines, which is one higher than their degree. The default of 4 gives cubic splines. If base=”fspline”, order is ignored.
#' @param basis an integer variable specifying the number of basis functions used to expand genetic effect function. If base=”fspline”, basis should be an odd integer.
#' @param covariate a data frame of covariate data. This data frame should contrain at least three columns: pedigree id, person id, and at least one covariate. The first two colums should be named as ped and person, respectively.
#' @param base a character value to specify basis function system to expand the genetic effect function. Options are “bspline” for b-splines and “fspline” for Fourier splines.
#' @param interaction a logic value indicating whether gene-environment interaction effect should be considered in the model. The default FALSE assumes no interaction.
#' @importFrom fda create.bspline.basis
#' @importFrom fda create.fourier.basis
#' @importFrom fda eval.basis
#' @importFrom stats glm
#' @importFrom stats anova
#'
#' @return a list of p values of testing the genetic effect, given by likelihood ratio test, chi-squared test, and F-distributed test, respectively.
#' @export
#'
flm_beta_smooth_only=function(pheno, mode = "Additive", geno, pos, order, basis, covariate, base = "bspline", interaction = FALSE)
{
  geno[is.na(geno)]=0
  covariate[is.na(covariate)] = 0
  ### define genotyping matrix for Dom and Rec modes ###
  ### For Dom mode, redefine geno[i,j] = 1 if geno[i,j] = 1 or 2
  ### For Rec mode, redefine geno[i,j] = 1 if geno[i,j] = 2
  geno_X = geno * 0
  for (i in 1:nrow(geno))
    for (j in 1:ncol(geno))
    {
      if (mode == "Dom")
      {
        if (geno[i, j] == 1 || geno[i, j] == 2)
          geno_X[i,j] = 1
      }
      else if ( mode == "Rec")
        if (geno[i, j] == 2)
          geno_X[i, j] = 1
    }

  if  ( mode == "Rec" || mode == "Dom")
    geno = geno_X

  idx     = is.na(pheno)
  pheno   = pheno[!idx]
  geno    = geno[!idx,]
  covariate = covariate[!idx,]
  dqr     = qr(geno)
  index   = dqr$pivot[1:dqr$rank]
  geno    = geno[, index]
  pos     = pos[index]
  nsample = nrow(geno)
  nsnp    = ncol(geno)

  if(max(pos) > 1) {
    pos = (pos - min(pos)) / (max(pos) - min(pos))
  }

  if (base ==  "bspline"){
    betabasis  = create.bspline.basis(norder = order, nbasis = basis)
  } else if (base == "fspline"){
    betabasis  = create.fourier.basis(c(0,1), nbasis = basis)
  }else { }

  B = eval.basis(pos, betabasis)

  UJ = geno %*% B
  #########
  pval = list()

  if (interaction == FALSE)
  {
    fit        = glm (pheno ~ covariate + UJ, family = "gaussian")
    pval$LRT   = anova(fit, test = "LRT")[3,5]
    pval$Chisq = anova(fit, test = "Chisq")[3,5]
    pval$F     = anova(fit, test = "F")[3,6]

    #gtfit      = gt(pheno ~ covariate, pheno ~ covariate + UJ, model = "linear")
    #pval$gt    =  p.value(gtfit)
  }else
  {
    Interaction = matrix (0, ncol =  ncol(UJ), nrow =nrow(UJ) )
    for ( i in 1:nrow(Interaction) ) { Interaction[i,] = covariate[i,1] * UJ[i,] }

    fitInt     = glm (pheno ~ covariate + UJ + Interaction, family = "gaussian")
    pval$LRT   = anova(fitInt, test = "LRT")[4,5]
    pval$Chisq = anova(fitInt, test = "Chisq")[4,5]
    pval$F     = anova(fitInt, test = "F")[4,6]

    #gtfitInt   = gt(pheno ~ covariate + UJ, pheno ~ covariate + UJ + Interaction, model = "linear")
    #pval$gt    =  p.value(gtfitInt)
  }

  pval
}

