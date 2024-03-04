library(fda)
library(MASS)
library(Matrix)
#library(globaltest)

### This Function is modified by Ruzong Fan, March 20, 2013 ###



#' Title
#'
#' @param pheno describe
#' @param mode describe
#' @param geno describe
#' @param pos describe
#' @param order describe
#' @param bbasis describe
#' @param fbasis describe
#' @param gbasis describe
#' @param covariate describe
#' @param base describe
#' @param interaction describe
#'
#' @return a value
#' @export
#'

gflm_fixed_model=function(pheno, mode = "Additive", geno, pos, order, bbasis, fbasis, gbasis, covariate, base = "bspline", interaction = FALSE)
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
    betabasis  = create.bspline.basis(norder = order, nbasis = bbasis)
    genobasis  = create.bspline.basis(norder = order, nbasis = gbasis)
  } else if (base == "fspline"){
    betabasis  = create.fourier.basis(c(0,1), nbasis = fbasis)
    genobasis  = create.fourier.basis(c(0,1), nbasis = gbasis)
  }else { }

  B = eval.basis(pos, genobasis)

  to_mul = ginv(t(B) %*% B) %*% t(B)
  #to_mul = solve( nearPD(t(B)%*%B)$mat ) %*% t(B)  ###nearPD requires library Matrix###
  U      = geno %*% t( to_mul )

  J      = inprod(genobasis, betabasis)    ### added by Fan ###

  UJ = matrix( U %*% J, ncol = ncol(J) )

  ### Make sure UJ has full rank of bbasis or fbasis ###
  UJdqr   = qr(UJ)
  UJindex = UJdqr$pivot[1:UJdqr$rank]
  UJ      = UJ[, UJindex]
  ###

  #fitNull   =  glm (pheno ~ covariate, family = "binomial")
  pval = list()

  gtfit    = gt(pheno ~ covariate, pheno ~ covariate + UJ, model = "logistic")

  if (interaction == FALSE)
  {
    fit        = glm (pheno ~ covariate + UJ, family = "binomial")
    pval$LRT   = anova(fit, test = "LRT")[3,5]
    pval$Chisq = anova(fit, test = "Chisq")[3,5]
    pval$Rao   = anova(fit, test = "Rao")[3,6]
    pval$gt    =  p.value(gtfit)
  }else
  {
    Interaction = matrix (0, ncol =  ncol(UJ), nrow =nrow(UJ) )
    for (i in 1:nrow(Interaction) ) { Interaction[i,] = covariate[i,2] * UJ[i,]}
    gtfitInt = gt(pheno ~ covariate + UJ, pheno ~ covariate + UJ + Interaction, model = "logistic")
    fitInt     = glm (pheno ~ covariate + UJ + Interaction, family = "binomial")
    pval$LRT   = anova(fitInt, test = "LRT")[4,5]
    pval$Chisq = anova(fitInt, test = "Chisq")[4,5]
    pval$Rao   = anova(fitInt, test = "Rao")[4,6]
    pval$gt    =  p.value(gtfitInt)
  }

  pval
}

