library(fda)
library(MASS)
library(Matrix)
library(survival)

### This Function is modified by Ruzong Fan, March 20, 2013 ###


#' Title
#'
#' @param time a numeric vector of follow up time for right censored data.
#' @param delta a numeric vector of status indicator for right censored data, normally 0=no event, 1=event.
#' @param mode a character value specifying the discrete realization of genetic variant function. The default of “Additive” assumes additive effect of minor alleles, “Dom” assume dominant effect, and “Rec” assumes recessive effect.
#' @param geno a matrix of genotype data only. The number of rows should be equal to the number of individuals in pheno. The number of columns is the total number of SNPs.
#' @param pos a vector of physical positions of the SNP used in geno. The length of pos should be equal to the column number of geno.
#' @param order an integer specifying the order of b-splines basis function, which is one higher than their degree. The default of 4 gives cubic splines. If base=”fspline”, order is ignored.
#' @param bbasis an integer variable specifying the number of b-spine basis functions used to expand genetic effect function.
#' @param fbasis an odd integer variable specifying the number of Fourier spline basis functions used to expand genetic effect function.
#' @param gbasis an integer variable specifying the number of basis functions used to expand genetic variant function. If base=”fspline”, basis should be an odd integer.
#' @param covariate a data frame of covariate data. This data frame should contrain at least three columns: pedigree id, person id, and at least one covariate. The first two colums should be named as ped and person, respectively.
#' @param base a character value to specify basis function system to expand the genetic effect function and genetic variant function. Options are “bspline” for b-splines and “fspline” for Fourier splines.
#'
#' @return a numeric value of p values of testing the genetic effect given by likelihood ratio test.
#' @export

flm_surv_fixed_model=function(time, delta, mode = "Additive", geno, pos, order, bbasis, fbasis, gbasis, covariate, base = "bspline")
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

  pval = list()

  fit        = coxph(Surv(time, delta) ~ as.matrix( covariate ) + as.matrix( UJ ) )
  #pval$LRT   = anova(fit, test = "LRT")[3, 5]
  pval$Chisq = anova(fit, test = "Chisq")[3, 4]
  #pval$Rao   = anova(fit, test = "Rao")[3, 6]

  if (pval$Chisq == 0)
  {
    Stat = anova(fit, test = "Chisq")[3, 2]
    df = anova(fit, test = "Chisq")[3, 3]
    pval$Chisq = pchisq(Stat, df, lower.tail=FALSE)
  }

  pval
}

