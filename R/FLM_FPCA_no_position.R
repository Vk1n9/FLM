library(globaltest)
library(fda)
library(MASS)
library(Matrix)
library(mgcv)
### This Function is created by Ruzong Fan, Dec 18, 2012 ###

#' Title
#'
#' @param x describe
#'
#' @return a value

Trace.Matrix = function(x)
{
  sum(diag(x))
}

#' Title
#'
#' @param pheno @param pheno a data frame of phenotype data. This data frame should contain three columns: pedigree id, person id, and phenotype value. These three columns should be named as ped, person, and trait, respectively.
#' @param mode "Additive"
#' @param geno a matrix of genotype data only. The number of rows should be equal to the number of individuals in pheno. The number of columns is the total number of SNPs.
#' @param covariates a data frame of covariate data. This data frame should contrain at least three columns: pedigree id, person id, and at least one covariate. The first two colums should be named as ped and person, respectively.
#' @param kz 30
#' @param kb 30
#' @param smooth.cov FALSE
#' @param family "gaussian"
#'
#' @importFrom stats predict
#' @importFrom stats anova
#' @importFrom stats quantile
#' @importFrom stats glm
#' @importFrom globaltest gt
#' @importFrom globaltest p.value
#' @importFrom mgcv gam
#' @import rainbow
#' @import MASS
#'
#' @return a list of p values of testing the genetic effect, given by likelihood ratio test, chi-squared test, and F-distributed test, respectively.
#' @export
#'
flm_fpca_no_position <- function(pheno, mode = "Additive", geno, covariates, kz = 30, kb = 30, smooth.cov = FALSE, family = "gaussian") {

  ####
  geno[is.na(geno)] = 0
  covariates[is.na(covariates)] = 0

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

  if ( mode == "Rec" || mode == "Dom")
    geno = geno_X

  idx        = is.na(pheno)
  pheno      = pheno[!idx]
  geno       = geno[!idx,]
  if (is.vector(covariates)==FALSE )
  {covariates = covariates[!idx,]
  dimcov = dim(covariates)[2]
  }
  else if (is.vector(covariates))
  {covariates[!idx]
    dimcov = 1
  }

  dqr     = qr(geno)
  index   = dqr$pivot[1:dqr$rank]
  geno    = geno[, index]

  kb = min(kz, kb)
  n  = length(pheno)
  p  = ifelse(is.null(covariates), 0, dimcov)

  if (is.matrix(geno)){
    Funcs = list(length=1)
    Funcs[[1]] = geno
  } else {
    Funcs = geno
  }

  # functional predictors
  N.Pred = length(Funcs)

  t = phi = psi = list(length=N.Pred)
  for (i in 1:N.Pred){
    t[[i]] = seq(0, 1, length = dim(Funcs[[i]])[2])
    N_obs  = length(t[[i]])

    # de-mean the functions
    meanFunc   = apply(Funcs[[i]], 2, mean, na.rm=TRUE)
    resd       = sapply(1:length(t[[i]]), function(u) Funcs[[i]][,u]-meanFunc[u])
    Funcs[[i]] = resd

    # construct and smooth covariance matrices
    G.sum   <- matrix(0, N_obs, N_obs)
    G.count <- matrix(0, N_obs, N_obs)

    for (j in 1:dim(resd)[1]){
      row.ind = j
      temp    = resd[row.ind, ] %*% t( resd[row.ind, ])
      G.sum   <- G.sum + replace(temp, which(is.na(temp)), 0)
      G.count <- G.count + as.numeric(!is.na(temp))
    }
    G <- ifelse(G.count==0, NA,  G.sum/G.count)

    ## get the eigen decomposition of the smoothed variance matrix
    if (smooth.cov){
      G2 <- G
      M  <- length(t[[i]])
      diag(G2)= rep(NA, M)
      g2 <- as.vector(G2)
      ## define a N*N knots for bivariate smoothing
      N <- 10

      ## bivariate smoothing using the gamm function
      t1 <- rep(t[[i]], each=M)
      t2 <- rep(t[[i]], M)
      newdata <- data.frame(t1 = t1, t2 = t2)
      K.0  <- matrix(predict(gam(as.vector(g2) ~ te(t1, t2, k = N)), newdata), M, M) # smooth K.0
      K.0 <- (K.0 + t(K.0)) / 2

      eigenDecomp <- eigen(K.0)
    } else {
      eigenDecomp <- eigen(G)
    }

    ### make sure no "Error in eigenDecomp$vectors[, 1:kz] : subscript out of bounds"
    kz1 = min(kz, length(eigenDecomp$values))
    psi[[i]] = eigenDecomp$vectors[,1:kz1]

    # set the basis to be used for beta(t)
    num = kb-2
    qtiles   <- seq(0, 1, length = num + 2)[-c(1, num + 2)]
    knots    <- quantile(t[[i]], qtiles)
    phi[[i]] = cbind(1, t[[i]], sapply(knots, function(k) ((t[[i]] - k > 0) * (t[[i]] - k))))

    C = matrix(0, nrow=dim(resd)[1], ncol=kz1)

    for (j in 1:dim(resd)[1]){
      C[j,] <-  replace(resd[j,], which(is.na(resd[j,])), 0) %*% psi[[i]][ ,1:kz1 ]
    }

    J = t(psi[[i]]) %*% phi[[i]]
    CJ = C %*% J
  }

  pval = list() #anova(fit0, fit, test = "LRT")[2, 5]

  fit   = glm (pheno ~ covariates + CJ, family )

  if (family == "gaussian"){
    gtmodel = "linear"
  } else if (family == "binomial"){
    gtmodel = "logistic"
  }else { }

  gtfit = gt(pheno ~ covariates, pheno ~ covariates + CJ, model = gtmodel )

  if (family == "binomial")
  {
    pval$LRT   = anova(fit, test = "LRT")[3,5]
    pval$Chisq = anova(fit, test = "Chisq")[3,5]
    pval$Rao   = anova(fit, test = "Rao")[3,6]
    pval$gt    = p.value(gtfit)
  }
  else if (family == "gaussian")
  {
    pval$LRT = anova(fit, test = "LRT")[3,5]
    pval$Chisq = anova(fit, test = "Chisq")[3,5]
    pval$F   = anova(fit, test = "F")[3,6]
    #pval$gt    = p.value(gtfit)
  }

  pval
}

