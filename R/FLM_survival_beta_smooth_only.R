library(fda)
library(MASS)
library(Matrix)
library(survival)

### This Function is modified by Ruzong Fan, March 20, 2013 ###
#' Title
#'
#' @param time describe
#' @param delta describe
#' @param mode describe
#' @param geno describe
#' @param pos describe
#' @param order describe
#' @param basis describe
#' @param covariate describe
#' @param base describe
#'
#' @importFrom survival coxph
#' @importFrom stats pchisq
#' @return a value
#' @export
#'
#'
flm_surv_beta_smooth_only = function(time, delta, mode = "Additive", geno, pos, order, basis, covariate, base = "bspline")
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

  idx     = is.na(time)
  time    = time[!idx]
  delta   = delta[!idx]
  geno    = geno[!idx,]
  if (is.vector(covariate)==FALSE )
  {
    covariate = covariate[!idx,]
  } else if (is.vector(covariate))
  {
    covariate[!idx]
  }
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

  ### Make sure UJ has full rank of bbasis or fbasis ###
  UJdqr   = qr(UJ)
  UJindex = UJdqr$pivot[1:UJdqr$rank]
  UJ      = UJ[, UJindex]
  ###

  #fitNull   =  coxph(Surv(time, delta) ~ as.matrix( covariate ) )
  pval = list()

  fit        = coxph(Surv(time, delta) ~ as.matrix( covariate ) + as.matrix( UJ ) )
  #pval$LRT   = anova(fit, test = "LRT")[3,5]
  pval$Chisq = anova(fit, test = "Chisq")[3, 4]
  #pval$Rao   = anova(fit, test = "Rao")[3,6]

  if (pval$Chisq == 0)
  {
    Stat = anova(fit, test = "Chisq")[3, 2]
    df = anova(fit, test = "Chisq")[3, 3]
    pval$Chisq = pchisq(Stat, df, lower.tail=FALSE)
  }

  pval
}

