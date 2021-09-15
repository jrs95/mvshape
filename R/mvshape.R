#####################################################
##### Fractional polynomials #####
#####################################################

#' fracpoly
#'
#' fracpoly fits the best fitting fractional polynomial of degree 1 and 2.
#' @param y outcome.
#' @param x exposure.
#' @param covar data.frame with covariates.
#' @param family the glm family (options: gaussian and binomial).
#' @return List of best-fitting polynomials of degrees 1 and 2 as well as associated statistics.
#' @return \item{power_d1}{power of the best-fitting fractional polynomial of degree 1}
#' @return \item{fp1}{model of the best-fitting fractional polynomial of degree 1}
#' @return \item{power_d2}{powers of the best-fitting fractional polynomial of degree 2}
#' @return \item{fp2}{model of the best-fitting fractional polynomial of degree 2}
#' @return \item{p_d1}{p-value testing the best-fitting fractional polynomial of degree 1 against the linear model}
#' @return \item{p_d2}{p-value testing the best-fitting fractional polynomial of degree 2 against the best-fitting fractional polynomial of degree 2}
#' @return \item{xmin}{miniumum value of the exposure}
#' @return \item{xmax}{maximum value of the exposure}
#' @return \item{family}{family used in the analysis}
#' @examples
#' ### Data
#' y <- rnorm(5000)
#' x <- rnorm(5000,10,1)
#' c1 <- rbinom(5000,1,0.5)
#' c2 <- rnorm(5000)
#' covar <- data.frame(c1=c1, c2=c2)
#'
#' ### Analyses
#' res <- fracpoly(y=y, x=x, covar=covar, family="gaussian")
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export

fracpoly <- function(y=y, x=x, covar=NULL, family="gaussian"){
  
  # Errors
  if(length(y)!=length(x) | (if(!is.null(covar)){(nrow(covar)!=length(y))}else{FALSE})) stop("the size of the outcome is not equal to the size of the exposure or covariates")
  if(!is.null(covar)){
    if(!is.data.frame(covar)) stop("the covar object is not a data.frame")
    if(any(names(covar) %in% c("y", "x"))) stop("do not call covariates 'y' or 'x'")
  }
  if(any(x<=0)) stop("The values of x are not all positive")
  if(!(family=="gaussian" | family=="binomial")) stop("The generalised linear model has to be either binomial or gaussian")
  
  # Remove missing values
  if(is.null(covar)){missing <- is.na(y) | is.na(x)}else{missing <- is.na(y) | is.na(x) | apply(is.na(covar),1,any)}
  if(length(y)!=length(missing)) stop("the size of the outcome is not equal to the size of the missing variable")
  y <- y[!missing]
  x <- x[!missing]
  if(!is.null(covar)){covar_names <- names(covar); covar <- as.data.frame(covar[!missing,,drop=F]); names(covar) <- covar_names}
  
  # FP degree 1
  likelihood_d1 <- NULL
  powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  
  for(pi in powers){
    if(pi==0){xfp <- log(x)}else{xfp <- x^pi}
    if(length(covar)==0){model <- glm(y~xfp, family=family)}else{model <- glm(y~xfp + ., data=covar, family=family)}
    if(pi==-2){fp1 <- model; power_bfp1 <- pi}
    else{if(logLik(model)>=max(likelihood_d1)){fp1 <- model; power_bfp1 <- pi}}
    likelihood_d1 <- c(likelihood_d1, logLik(model)) 
  }
  
  maxlik_d1 <- max(likelihood_d1)
  # power_bfp1 <- powers[which.max(rank(likelihood_d1, ties.method = "first", na.last=FALSE))]
  chi2 <- (-2*likelihood_d1[6]) - (-2*maxlik_d1)
  p_d1 <- 1 - pchisq(chi2, df=1)
  
  # FP degree 2
  likelihood_d2 <- NULL
  powers1 <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  powers2 <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  # power_d2 <- data.frame(p1=NULL, p2=NULL)
  
  for(pi1 in powers1){
    if(pi1==0){xfp1 <- log(x)}else{xfp1 <- x^pi1}
    for(pi2 in powers2){
      if(pi1==pi2){if(pi2==0){xfp2 <- log(x)*log(x)}else{xfp2 <- x^pi2*log(x)}}
      else{if(pi2==0){xfp2 <- log(x)}else{xfp2 <- x^pi2}}
      if(length(covar)==0){model <- glm(y~xfp1 + xfp2, family=family)}else{model <- glm(y~xfp1 + xfp2 + ., data=covar, family=family)}
      if(pi1==-2 & pi2==-2){fp2 <- model; powers_bfp2 <- c(pi1, pi2)}
      else{if(logLik(model)>=max(likelihood_d2)){fp2 <- model; powers_bfp2 <- c(pi1, pi2)}}
      likelihood_d2 <- c(likelihood_d2, logLik(model))
      # power_d2 <- rbind(power_d2, data.frame(p1=pi1, p2=pi2)) 
    }
    powers2 <- powers2[-1]
  }
  
  maxlik_d2 <- max(likelihood_d2)
  # powers_bfp2 <- power_d2[which.max(rank(likelihood_d2, ties.method = "first", na.last=FALSE)),]
  chi2 <- (-2*maxlik_d1) - (-2*maxlik_d2)
  p_d2 <- 1 - pchisq(chi2,df=2)
  
  # Results
  results <- list(power_d1=power_bfp1, fp1=fp1, power_d2=powers_bfp2, fp2=fp2, p_d1=p_d1, p_d2=p_d2, xmin=min(x), xmax=max(x), family=family)
  class(results) <- "fracpoly"
  
  # Return
  return(results)
  
}

#' Print fracpoly
#'
#' print method for class "fracpoly".
#' @param x an object of class "fracpoly".
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
print.fracpoly <- function(x, ...){
  cat("Call: \nfracpoly")
  cat("\n\nBest-fitting fractional polynomial of degree 1", sep="")
  cat("\nPower: ",x$power_d1,"\n", sep="")
  cat("\nCoefficients:\n", sep="")
  cat(x$fp1$coef)
  cat("\n\nBest-fitting fractional polynomial of degree 2", sep="")
  cat("\nPowers:",as.vector(unlist(x$power_d2)),"\n", sep=" ")
  cat("\nCoefficients:\n", sep="")
  cat(x$fp2$coef)
  cat("\n\n")
}

#' Summary fracpoly
#'
#' summary method for class "fracpoly".
#' @param x an object of class "fracpoly".
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
summary.fracpoly <- function(x, ...){
  cat("Call: \nfracpoly")
  cat("\n\nBest-fitting fractional polynomial of degree 1", sep="")
  cat("\nPower: ",x$power_d1,"\n", sep="")
  print(summary(x$fp1), row.names = FALSE)
  cat("\nBest-fitting fractional polynomial of degree 2", sep="")
  cat("\nPowers:",as.vector(unlist(x$power_d2)),"\n", sep=" ")
  print(summary(x$fp2), row.names = FALSE)
}

#####################################################
##### Mvshape #####
#####################################################

#' mvshape
#'
#' mvshape fits a multivariate meta-analysis for groups of the exposure (e.g. deciles).
#' @import mvmeta
#' @param y outcome.
#' @param x exposure.
#' @param covar data.frame with covariates.
#' @param study study variable.
#' @param ngrp number of quantiles of the exposure.
#' @param refgrp reference group.
#' @param family the glm family (options: gaussian and binomial).
#' @param float floating point variances.
#' @param method meta-analysis method.
#' @return List of multivariate meta-analysis results for each group.
#' @return \item{results}{data.frame of results: q is the quantile group, xbeta is the mean of x in each quantile, xse is the standard error of the mean of x in each quantile, beta is the regression coefficient of association between y and each quantile of x, se is the standard error of the regression coefficient of association between y and each quantile of x}
#' @return \item{varcor}{variance-covariance matrix}
#' @return \item{xmin}{miniumum value of the exposure}
#' @return \item{xmax}{maximum value of the exposure}
#' @return \item{family}{family used in the analysis}
#' @examples
#' ### Data
#' y <- rnorm(5000)
#' x <- rnorm(5000,10,1)
#' c1 <- rbinom(5000,1,0.5)
#' c2 <- rnorm(5000)
#' covar <- data.frame(c1=c1, c2=c2)
#' study <- c(rep("study1",1000),rep("study2",1000),rep("study3",1000),rep("study4",1000),rep("study5",1000))
#'
#' ### Analyses
#' res <- mvshape(y=y, x=x, covar=covar, study=study, family="gaussian")
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export

mvshape <- function(y=y, x=x, covar=NULL, study=NULL, ngrp=10, refgrp=1, family="gaussian", float=FALSE, method="reml"){

  # Reference group
  refgrp <- as.integer(refgrp)
  
  # Errors
  if(length(y)!=length(x) | (if(!is.null(covar)){(nrow(covar)!=length(y))}else{FALSE})) stop("the size of the outcome is not equal to the size of the exposure or covariates")
  if(!is.null(covar)){
    if(!is.data.frame(covar)) stop("the covar object is not a data.frame")
    if(any(names(covar) %in% c("y", "x", "xq"))) stop("do not call covariates 'y', 'x' or 'xq'")
  }
  if(!is.null(study)){if(length(y)!=length(study)) stop("the size of the outcome is not equal to the size of the study variable")}
  if(!is.null(study)){if(any(is.na(study))) stop("there are NAs in the study variable")}
  if(refgrp > ngrp) stop("the refgrp is larger than the number of groups")
  if(refgrp < 1) stop("the refgrp has to be at least one")
  if(!(family=="gaussian" | family=="binomial")) stop("the generalised linear model has to be either binomial or gaussian")

  # Remove missing values
  if(length(covar)==0){missing <- is.na(y) | is.na(x)}else{missing <- is.na(y) | is.na(x) | apply(is.na(covar),1,any)}
  if(length(y)!=length(missing)) stop("the size of the outcome is not equal to the size of the missing variable")
  y <- y[!missing]
  x <- x[!missing]
  if(!is.null(covar)){covar_names <- names(covar); covar <- as.data.frame(covar[!missing,,drop=F]); names(covar) <- covar_names}
  if(!is.null(study)){study <- study[!missing]; study <- as.factor(study)}else{study <- as.factor(rep("study",length(y)))}
  
  # Quantiles
  prob <- 1/ngrp
  quantiles <- quantile(x, probs=seq(0,1, prob))
  xq <- cut(x, quantiles, include.lowest=T)
  xq <- factor(as.numeric(xq))
  xq <- relevel(xq, ref = refgrp)

  # Study-specific estimates
  est <- data.frame()
  est_var <- data.frame()
  x_mean <- data.frame()
  x_var <- data.frame()

  for(coh in levels(study)){
    
    # Outcome estimates
    if(!is.null(covar)){model <- glm(y[study==coh] ~ xq[study==coh], family=family)}else{covar_coh <- as.data.frame(covar[study==coh,,drop=F]); names(covar_coh) <- covar_names; model <- glm(y[study==coh] ~ xq[study==coh] + ., data=covar_coh, family=family)}
    if(family=="gaussian"){
      if(any(is.na(model$coef))) stop("there are missing regression coefficients") 
      b <- model$coef[1:ngrp]
      varcov <- vcov(model)[1:ngrp,1:ngrp]
      v <- varcov[lower.tri(varcov, diag = TRUE)]
    }else{
      if(any(is.na(model$coef))) stop("there are missing regression coefficients") 
      b <- model$coef[2:ngrp]
      varcov <- vcov(model)[2:ngrp,2:ngrp]
      v <- varcov[lower.tri(varcov, diag = TRUE)]
    }
    est <- rbind(est, b)
    names(est)<- paste0("beta_",1:length(b))
    est_var <- rbind(est_var, v)
    names(est_var)<- paste0("cov_",1:length(v))
    
    # Exposure means
    model_mean <- lm(x[study==coh]~xq[study==coh] - 1)
    if(any(is.na(model_mean$coef))) stop("there are missing mean exposure estimates")     
    x_mean <- rbind(x_mean, model_mean$coef)
    names(x_mean)<- paste0("mean_",1:length(x_mean))
    varcov_mean <- vcov(model)
    v_mean <- varcov_mean[lower.tri(varcov_mean, diag = TRUE)]
    x_var <- rbind(x_var, v_mean)
    names(x_var)<- paste0("cov_",1:length(x_var))
    
  }
  
  if(length(levels(study))>1){
    
    # Outcome mvmeta
    est <- as.matrix(est); class(est) <- "numeric"
    est_var <- as.matrix(est_var); class(est_var) <- "numeric"
    mvmodel <- suppressWarnings(mvmeta(est,est_var,method=method))
    
    if(family=="gaussian"){
      if(length(mvmodel$coef)!=ngrp | any(is.na(mvmodel$coef))) stop("there are missing mvmeta coefficients")
      beta <- c(mvmodel$coef)
      se <- c(summary(mvmodel)$coefficients[,2]); names(se) <- NULL
      varcov <- vcov(mvmodel); row.names(varcov) <- paste0("q_", levels(xq)); colnames(varcov) <- paste0("q_", levels(xq))
    }else{
      if(length(mvmodel$coef)!=(1-ngrp) | any(is.na(mvmodel$coef))) stop("there are missing mvmeta coefficients")
      beta <- c(0,mvmodel$coef)
      se <- c(0,summary(mvmodel)$coefficients[,2]); names(se) <- NULL
      if(float==TRUE){se <- sqrt(floatvar(vcov(mvmodel))$variance); names(se) <- NULL}
      varcov <- vcov(mvmodel); row.names(varcov) <- paste0("q_", as.character(levels(xq))[-1]); colnames(varcov) <- paste0("q_", as.character(levels(xq))[-1])
    }
    
    # Exposure mvmeta
    x_mean <- as.matrix(x_mean); class(x_mean) <- "numeric"
    x_var <- as.matrix(x_var); class(x_var) <- "numeric"
    mvmodel_mean <- suppressWarnings(mvmeta(x_mean,x_var,method=method))
    xbeta <- c(mvmodel_mean$coef); names(xbeta) <- NULL
    xse <- c(summary(mvmodel_mean)$coefficients[,2]); names(xse) <- NULL
    xbeta[xse==0] <- x_mean[1,xse==0]
    
  }
  else{
    
    # Outcome model
    if(family=="gaussian"){
      beta <- c(model$coef[1:ngrp])
      se <- c(summary(model)$coefficients[1:ngrp,2]); names(se) <- NULL
      varcov <- varcov <- vcov(model)[1:ngrp,1:ngrp]; row.names(varcov) <- paste0("q_", levels(xq)); colnames(varcov) <- paste0("q_", levels(xq))
    }else{
      beta <- c(0,model$coef[2:ngrp])
      se <- c(0,summary(model)$coefficients[2:ngrp,2]); names(se) <- NULL
      if(float==TRUE){se <- sqrt(floatvar(vcov(model)[2:ngrp,2:ngrp])$variance); names(se) <- NULL}
      varcov <- vcov(model)[2:ngrp,2:ngrp]; row.names(varcov) <- paste0("q_", as.character(levels(xq))[-1]); colnames(varcov) <- paste0("q_", as.character(levels(xq))[-1])
    }
    
    # Exposure mean
    xbeta <- c(model_mean$coef); names(xbeta) <- NULL
    xse <- c(summary(model_mean)$coefficients[,2]); names(xse) <- NULL
    
  }
  
  # Results
  results <- list(results=data.frame(q=as.character(levels(xq)), xbeta=xbeta, xse=xse, beta=beta, se=se, stringsAsFactors = F), varcov=varcov, xmin=min(x), xmax=max(x), family=family)
  class(results) <- "mvshape"
  return(results)
}

#' Print mvshape
#'
#' print method for class "mvshape".
#' @param x an object of class "mvshape".
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
print.mvshape <- function(x, ...){
  cat("Call: \nmvshape\n")
  cat("\nCoefficients:\n", sep="")
  x$results$q[1] <- paste0(x$results$q[1], " (intercept)")
  print(x$results[,c("q", "xbeta", "beta", "se")], row.names = FALSE)
}
  
#' Summary mvshape
#'
#' summary method for class "mvshape".
#' @param x an object of class "mvshape".
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
summary.mvshape <- function(x, ...){
  cat("Call: \nmvshape\n")
  cat("\nCoefficients:\n", sep="")
  x$results$q[1] <- paste0(x$results$q[1], " (intercept)")
  print(x$results[,c("q", "xbeta", "beta", "se")], row.names = FALSE)
}
  
#####################################################
##### Figures #####
#####################################################

#' fracpoly_plot
#'
#' fracpoly_plot plots the best fitting fractional polynomial.
#' @import ggplot2
#' @param fracpoly a fracpoly object.
#' @param degree the fractonal polyniomal degree to be plotted (options: 1, 2 and both).
#' @param xref the reference point for the figures of binary outcomes.
#' @param logx plot the x-axis on the log-scale.
#' @param logy plot the y-axis on the log-scale.
#' @param pref_x the prefix for the x-axis.
#' @param pref_y the prefix for the y-axis.
#' @param xbreaks breaks in the x-axis.
#' @param ybreaks breaks in the y-axis.
#' @param cicolour colour of the 95\% confidence intervals.
#' @param betacolour colour of the regression line.
#' @param intcolour colour of the intercept line.
#' @param xlim x-axis limits.
#' @return fractional polynomial plot.
#' @examples
#' ### Data
#' y <- rnorm(5000)
#' x <- rnorm(5000,10,1)
#' c1 <- rbinom(5000,1,0.5)
#' c2 <- rnorm(5000)
#' covar <- data.frame(c1=c1, c2=c2)
#'
#' ### Analyses
#' res <- fracpoly(y=y, x=x, covar=covar, family="gaussian")
#' fracploy_plot(res)
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export

fracpoly_plot <- function(fracpoly, degree="both", xref=NULL, logx=FALSE, logy=FALSE, pref_x="x", pref_y="y", xbreaks=NULL, ybreaks=NULL, cicolour="grey", betacolour="black", intcolour="grey", xlim=NULL){

  # Family
  family <- fracpoly$family
  
  # Extract fractional polynomial estimates
  if(degree=="both"){if(fracpoly[[6]]>=0.05){degree<-1}; if(fracpoly[[6]]<0.05){degree<-2}}
  if(degree==1){fp <- fracpoly[[2]]; powers <- as.numeric(fracpoly[[1]])}
  if(degree==2){fp <- fracpoly[[4]]; powers <- as.numeric(fracpoly[[3]])}
  x <- runif(10000, fracpoly[[7]], fracpoly[[8]])
  if(!is.null(xlim)){x <- runif(10000, xlim[1], xlim[2])}
  if(is.null(xref)){xref <- mean(x)}

  # Continuous outcome
  if(family=="gaussian"){
    if(degree==1 & powers[1]==0){
      yest <- fp$coef[1] + fp$coef[2]*log(x)
      yse <- sqrt(vcov(fp)[1,1] + (log(x))^2*vcov(fp)[2,2] + 2*(log(x))*vcov(fp)[1,2])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==1 & powers[1]!=0){
      yest <- fp$coef[1] + fp$coef[2]*x^powers[1]
      yse <- sqrt(vcov(fp)[1,1] + (x^powers[1])^2*vcov(fp)[2,2] + 2*(x^powers[1])*vcov(fp)[1,2])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==2 & powers[1]==0 & powers[2]==0){
      yest <- fp$coef[1] + fp$coef[2]*log(x) + fp$coef[3]*log(x)*log(x)
      yse <- sqrt(vcov(fp)[1,1] + (log(x))^2*vcov(fp)[2,2] + (log(x)*log(x))^2*vcov(fp)[3,3] + 2*(log(x))*vcov(fp)[1,2] + 2*(log(x)*log(x))*vcov(fp)[1,3] + 2*(log(x))*(log(x)*log(x))*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==2 & powers[1]==0 & powers[2]!=0){
      yest <- fp$coef[1] + fp$coef[2]*log(x) + fp$coef[3]*x^powers[2]
      yse <- sqrt(vcov(fp)[1,1] + (log(x))^2*vcov(fp)[2,2] + (x^powers[2])^2*vcov(fp)[3,3] + 2*(log(x))*vcov(fp)[1,2] + 2*(x^powers[2])*vcov(fp)[1,3] + 2*(log(x))*(x^powers[2])*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==2 & powers[1]!=0 & powers[2]==0){
      yest <- fp$coef[1] + fp$coef[2]*x^powers[1] + fp$coef[3]*log(x)
      yse <- sqrt(vcov(fp)[1,1] + (x^powers[1])^2*vcov(fp)[2,2] + (log(x))^2*vcov(fp)[3,3] + 2*(x^powers[1])*vcov(fp)[1,2] + 2*(log(x))*vcov(fp)[1,3] + 2*(x^powers[1])*(log(x))*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==2 & powers[1]!=0 & powers[2]!=0 & powers[1]==powers[2]){
      yest <- fp$coef[1] + fp$coef[2]*x^powers[1] + fp$coef[3]*x^powers[2]*log(x)
      yse <- sqrt(vcov(fp)[1,1] + (x^powers[1])^2*vcov(fp)[2,2] + (x^powers[2]*log(x))^2*vcov(fp)[3,3] + 2*(x^powers[1])*vcov(fp)[1,2] + 2*(x^powers[2]*log(x))*vcov(fp)[1,3] + 2*(x^powers[1])*(x^powers[2]*log(x))*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    }
    if(degree==2 & powers[1]!=0 & powers[2]!=0 & powers[1]!=powers[2]){
      yest <- fp$coef[1] + fp$coef[2]*x^powers[1] + fp$coef[3]*x^powers[2]
      yse <- sqrt(vcov(fp)[1,1] + (x^powers[1])^2*vcov(fp)[2,2] + (x^powers[2])^2*vcov(fp)[3,3] + 2*(x^powers[1])*vcov(fp)[1,2] + 2*(x^powers[2])*vcov(fp)[1,3] + 2*(x^powers[1])*(x^powers[2])*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    }
    
    # Data for figure
    data <- data.frame(x=x, yest=yest, yse=yse, ylci=ylci, yuci=yuci)
    xminp <- fracpoly[[7]]; xmaxp <- fracpoly[[8]]; yintr <- 0
    if(!is.null(xlim)){xminp <- xlim[1]; xmaxp <- xlim[2]}
    if(logy==T){data$yest <- exp(data$yest); data$yuci <- exp(data$yuci); data$ylci <- exp(data$ylci); yintr <- 1}
    if(logx==T){data$xest <- exp(data$xest); xminp <- exp(xminp); xmaxp <- exp(xmaxp)}
    data$xminp <- xminp; data$xmaxp <- xmaxp
    hline <- data.frame(x=c((xminp-0.1*xminp), (xmaxp+0.1*xmaxp)), y=c(yintr, yintr))
    if(findInterval(0, c(min(ylci),max(yuci)))==1){hlinep <- TRUE}else{hlinep <- FALSE}
    
    # Figure
    figure <- ggplot(data, aes(x=x))
    if(hlinep==TRUE){figure <- figure + geom_line(aes(x=x, y=y), color=intcolour, size=0.35, data=hline)}
    figure <- figure + geom_line(aes(y=ylci), color=cicolour) + geom_line(aes(y=yuci), color=cicolour) + geom_line(aes(y=yest), color=betacolour) + theme_bw() + labs(x=pref_x,y=pref_y) + theme(axis.title.x = element_text(vjust=-0.5, size=14), axis.title.y = element_text(vjust=1.5, size=14), axis.text=element_text(size=12)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(expand = c(0,0), limits=c(hline$x[1], hline$x[2]))
    if(!is.null(ybreaks)==TRUE){figure <- figure + scale_y_continuous(breaks=ybreaks)}
    suppressMessages(if(!is.null(xbreaks)==TRUE){figure <- figure + scale_x_continuous(expand = c(0,0), breaks=xbreaks)})
    if(logx==FALSE & logy==TRUE){ybreaks <- ggplot_build(figure)$layout$panel_ranges[[1]]$y.major_source; figure <- figure + scale_y_continuous(trans="log", breaks=ybreaks)}
    if(logx==TRUE & logy==FALSE){xbreaks <- ggplot_build(figure)$layout$panel_ranges[[1]]$x.major_source; figure <- figure + scale_x_continuous(trans="log", breaks=xbreaks)}
    if(logx==TRUE & logy==TRUE){ybreaks <- ggplot_build(figure)$layout$panel_ranges[[1]]$y.major_source; figure <- figure + scale_y_continuous(trans="log", breaks=ybreaks); xbreaks <- ggplot_build(figure)$layout$panel_ranges[[1]]$x.major_source; figure <- figure + scale_x_continuous(trans="log", breaks=xbreaks)}
    
  }
  
  # Binary outcome
  if(family=="binomial"){
    if(degree==1 & powers[1]==0){
      yest <- fp$coef[2]*log(x) - fp$coef[2]*log(xref)
      yse <- sqrt((log(x)-log(xref))^2*vcov(fp)[2,2])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==1 & powers[1]!=0){
      yest <- fp$coef[2]*x^powers[1] - fp$coef[2]*xref^powers[1]
      yse <- sqrt((x^powers[1]-xref^powers[1])^2*vcov(fp)[2,2])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==2 & powers[1]==0 & powers[2]==0){
      yest <- fp$coef[2]*log(x) + fp$coef[3]*log(x)*log(x) - (fp$coef[2]*log(xref) + fp$coef[3]*log(xref)*log(xref))
      yse <- sqrt((log(x)-log(xref))^2*vcov(fp)[2,2] + (log(x)*log(x)-log(xref)*log(xref))^2*vcov(fp)[3,3] + 2*(log(x)-log(xref))*(log(x)*log(x)-log(xref)*log(xref))*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==2 & powers[1]==0 & powers[2]!=0){
      yest <- fp$coef[2]*log(x) + fp$coef[3]*x^powers[2] - (fp$coef[2]*log(xref) + fp$coef[3]*xref^powers[2])
      yse <- sqrt((log(x)-log(xref))^2*vcov(fp)[2,2] + (x^powers[2]-xref^powers[2])^2*vcov(fp)[3,3] + 2*(log(x)-log(xref))*(x^powers[2]-xref^powers[2])*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==2 & powers[1]!=0 & powers[2]==0){
      yest <- fp$coef[2]*x^powers[1] + fp$coef[3]*log(x) - (fp$coef[2]*xref^powers[2] + fp$coef[3]*log(xref))
      yse <- sqrt((x^powers[1]-xref^powers[1])^2*vcov(fp)[2,2] + (log(x)-log(xref))^2*vcov(fp)[3,3] + 2*(x^powers[1]-xref^powers[1])*(log(x)-log(xref))*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==2 & powers[1]!=0 & powers[2]!=0 & powers[1]==powers[2]){
      yest <- fp$coef[2]*x^powers[1] + fp$coef[3]*x^powers[2]*log(x) - (fp$coef[2]*xref^powers[1] + fp$coef[3]*xref^powers[2]*log(xref))
      yse <- sqrt((x^powers[1]-xref^powers[1])^2*vcov(fp)[2,2] + (x^powers[2]*log(x)-xref^powers[2]*log(xref))^2*vcov(fp)[3,3] + 2*(x^powers[1]-xref^powers[1])*(x^powers[2]*log(x)-xref^powers[2]*log(xref))*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    }
    if(degree==2 & powers[1]!=0 & powers[2]!=0 & powers[1]!=powers[2]){
      yest <- fp$coef[2]*x^powers[1] + fp$coef[3]*x^powers[2] - (fp$coef[2]*xref^powers[1] + fp$coef[3]*xref^powers[2])
      yse <- sqrt((x^powers[1]-xref^powers[1])^2*vcov(fp)[2,2] + (x^powers[2]-xref^powers[2])^2*vcov(fp)[3,3] + 2*(x^powers[1]-xref^powers[1])*(x^powers[2]-xref^powers[2])*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    }
    
    # Data for figure
    data <- data.frame(x=x, yest=yest, yse=yse, ylci=ylci, yuci=yuci)
    xminp <- fracpoly[[7]]; xmaxp <- fracpoly[[8]]; yintr <- 0
    if(!is.null(xlim)){xminp <- xlim[1]; xmaxp <- xlim[2]}
    logy <- T
    data$yest <- exp(data$yest); data$yuci <- exp(data$yuci); data$ylci <- exp(data$ylci); yintr <- 1
    if(logx==T){data$x <- exp(data$x); xminp <- exp(xminp); xmaxp <- exp(xmaxp)}
    data$xminp <- xminp; data$xmaxp <- xmaxp
    hline <- data.frame(x=c((xminp-0.1*xminp), (xmaxp+0.1*xmaxp)), y=c(yintr, yintr))
    pref_y <- paste0("Odds ratio of ", pref_y)

    # Figure
    figure <- ggplot(data, aes(x=x))
    figure <- figure + geom_line(aes(x=x, y=y), color=intcolour, size=0.35, data=hline) + geom_line(aes(y=ylci), color=cicolour) + geom_line(aes(y=yuci), color=cicolour) + geom_line(aes(y=yest), color=betacolour)  + theme_bw() + labs(x=pref_x,y=pref_y) + theme(axis.title.x = element_text(vjust=-0.5, size=14), axis.title.y = element_text(vjust=1.5, size=14), axis.text=element_text(size=12)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(expand = c(0,0), limits=c(hline$x[1], hline$x[2]))
    if(!is.null(ybreaks)==TRUE){figure <- figure + scale_y_continuous(breaks=ybreaks)}
    suppressMessages(if(!is.null(xbreaks)==TRUE){figure <- figure + scale_x_continuous(expand = c(0,0), breaks=xbreaks)})
    if(logx==FALSE & logy==TRUE){ybreaks <- ggplot_build(figure)$layout$panel_ranges[[1]]$y.major_source; figure <- figure + scale_y_continuous(trans="log", breaks=ybreaks)}
    if(logx==TRUE & logy==FALSE){xbreaks <- ggplot_build(figure)$layout$panel_ranges[[1]]$x.major_source; figure <- figure + scale_x_continuous(trans="log", breaks=xbreaks)}
    if(logx==TRUE & logy==TRUE){ybreaks <- ggplot_build(figure)$layout$panel_ranges[[1]]$y.major_source; figure <- figure + scale_y_continuous(trans="log", breaks=ybreaks); xbreaks <- ggplot_build(figure)$layout$panel_ranges[[1]]$x.major_source; figure <- figure + scale_x_continuous(trans="log", breaks=xbreaks)}
  }
  
  # Return
  return(figure)   
}

#' mvshape_plot
#'
#' mvshape_plot plots the meta-analysis results against the mean values of the exposure in each group.
#' @param mvshape a mvshape object.
#' @param logx plot the x-axis on the log-scale.
#' @param logy plot the y-axis on the log-scale.
#' @param pref_x the prefix for the x-axis.
#' @param pref_y the prefix for the y-axis.
#' @param xbreaks breaks in the x-axis.
#' @param ybreaks breaks in the y-axis.
#' @param cicolour colour of the 95\% confidence intervals.
#' @param betacolour colour of the regression line.
#' @param intcolour colour of the intercept line.
#' @param xlim x-axis limits.
#' @return mvshape plot.
#' @examples
#' ### Data
#' y <- rnorm(5000)
#' x <- rnorm(5000,10,1)
#' c1 <- rbinom(5000,1,0.5)
#' c2 <- rnorm(5000)
#' covar <- data.frame(c1=c1, c2=c2)
#' study <- c(rep("study1",1000),rep("study2",1000),rep("study3",1000),rep("study4",1000),rep("study5",1000))
#'
#' ### Analyses
#' res <- mvshape(y=y, x=x, covar=covar, study=study, family="gaussian")
#' mvshape_plot(res)
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export

mvshape_plot <- function(mvshape, logx=FALSE, logy=FALSE, pref_x="x", pref_y="y", xbreaks=NULL, ybreaks=NULL, betacolour="red", cicolour="red", intcolour="grey",  xlim=NULL){

  # Family
  family <- mvshape$family
  
  # Continuous outcome
  if(family=="gaussian"){
    xest <- as.numeric(mvshape$results$xbeta)
    yest <- as.numeric(mvshape$results$beta) + mvshape$results$beta[1]
    yest[1] <- yest[1] - mvshape$results$beta[1]
    yvar <- mvshape$varcov[1,1]
    for(i in 1:length(xest)){
      yvar <- yvar + mvshape$varcov[i,i] + 2*mvshape$varcov[1,i]
    }
    yse <- sqrt(yvar)
    ylci <- yest - 1.96*yse 
    yuci <- yest + 1.96*yse
  }

  # Binary outcome
  if(family=="binomial"){
    xest <- as.numeric(mvshape$results$xbeta)
    yest <- as.numeric(mvshape$results$beta)
    ylci <- mvshape$results$beta - 1.96*as.numeric(mvshape$results$se) 
    yuci <- mvshape$results$beta + 1.96*as.numeric(mvshape$results$se)
  }  
  
  # Data for figure
  xminp <- mvshape[[3]]; xmaxp <- mvshape[[4]]; yintr <- 0
  if(!is.null(xlim)){xminp <- xlim[1]; xmaxp <- xlim[2]}
  data <- data.frame(xest=xest, yest=yest, ylci=ylci, yuci=yuci) 
  if(family=="binomial"){logy <- T; pref_y <- paste0("Odds ratio of ", pref_y)}
  if(logy==T){data$yest <- exp(data$yest); data$yuci <- exp(data$yuci); data$ylci <- exp(data$ylci); yintr <- 1}
  if(logx==T){data$xest <- exp(data$xest); xminp <- exp(xminp); xmaxp <- exp(xmaxp)}
  hline <- data.frame(x=c((xminp-0.1*xminp), (xmaxp+0.1*xmaxp)), y=c(yintr, yintr))
  if(findInterval(yintr, c(min(data$ylci),max(data$yuci)))==1){hlinep <- TRUE}else{hlinep <- FALSE}
  
  # Figure
  figure <- ggplot(data, aes(x=xest))
  if(hlinep==TRUE){figure <- figure + geom_line(aes(x=x, y=y), color=intcolour, size=0.35, data=hline)}
  figure <- figure + geom_errorbar(mapping=aes(x=xest, ymin=ylci, ymax=yuci), color=cicolour, width=0.025) + geom_point(aes(x=xest, y=yest), color=betacolour, data=data, size=4, shape=15) + theme_bw() + labs(x=pref_x,y=pref_y) + theme(axis.title.x = element_text(vjust=-0.5, size=14), axis.title.y = element_text(vjust=1.5, size=14), axis.text=element_text(size=12)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(expand = c(0,0), limits=c(hline$x[1], hline$x[2]))
  if(!is.null(ybreaks)==TRUE){figure <- figure + scale_y_continuous(breaks=ybreaks)}
  suppressMessages(if(!is.null(xbreaks)==TRUE){figure <- figure + scale_x_continuous(expand = c(0,0), breaks=xbreaks)})
  if(logx==FALSE & logy==TRUE){ybreaks <- ggplot_build(figure)$layout$panel_ranges[[1]]$y.major_source; figure <- figure + scale_y_continuous(trans="log", breaks=ybreaks)}
  if(logx==TRUE & logy==FALSE){xbreaks <- ggplot_build(figure)$layout$panel_ranges[[1]]$x.major_source; figure <- figure + scale_x_continuous(trans="log", breaks=xbreaks)}
  if(logx==TRUE & logy==TRUE){ybreaks <- ggplot_build(figure)$layout$panel_ranges[[1]]$y.major_source; figure <- figure + scale_y_continuous(trans="log", breaks=ybreaks); xbreaks <- ggplot_build(figure)$layout$panel_ranges[[1]]$x.major_source; figure <- figure + scale_x_continuous(trans="log", breaks=xbreaks)}
  
  # Return
  return(figure)
}

#' fracpoly_mvshape_plot
#'
#' fracpoly_mvshape_plot overlays the fractional polynomial plot with the mvshape plot.
#' @param fracpoly a fracpoly object.
#' @param degree the fractonal polyniomal degree to be plotted (options: 1, 2 and both).
#' @param mvshape a mvshape object.
#' @param xref the reference point for the figures of binary outcomes.
#' @param logx plot the x-axis on the log-scale.
#' @param logy plot the y-axis on the log-scale.
#' @param pref_x the prefix for the x-axis.
#' @param pref_y the prefix for the y-axis.
#' @param xbreaks breaks in the x-axis.
#' @param ybreaks breaks in the y-axis.
#' @param cicolour colour of the 95\% confidence intervals.
#' @param betacolour colour of the regression line.
#' @param intcolour colour of the intercept line.
#' @param xlim x-axis limits.
#' @return fractional polynomial and mvshape plot.
#' @examples
#' ### Data
#' y <- rnorm(5000)
#' x <- rnorm(5000,10,1)
#' c1 <- rbinom(5000,1,0.5)
#' c2 <- rnorm(5000)
#' covar <- data.frame(c1=c1, c2=c2)
#' study <- c(rep("study1",1000),rep("study2",1000),rep("study3",1000),rep("study4",1000),rep("study5",1000))
#'
#' ### Analyses
#' res <- mvshape(y=y, x=x, covar=covar, study=study, family="gaussian")
#' fracpoly_mvshape_plot(res)
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export

fracpoly_mvshape_plot <- function(fracpoly, degree="both", mvshape, xref=NULL, logx=FALSE, logy=FALSE, pref_x="x", pref_y="y", xbreaks=NULL, ybreaks=NULL, cicolour="grey", betacolour="black", intcolour="grey", mvbetacolour="red", mvcicolour="red", xlim=NULL, ylim=NULL){
  
  # Family
  family <- fracpoly$family
  family_mvshape <- mvshape$family
  if(family!=family_mvshape) stop("the family used in the fractional polynomial analysis is not the same as the family used in the mvshape analysis")
  
  # Fractional polynomial plot
  if(degree=="both"){if(fracpoly[[6]]>=0.05){degree<-1}; if(fracpoly[[6]]<0.05){degree<-2}}
  if(degree==1){fp <- fracpoly[[2]]; powers <- as.numeric(fracpoly[[1]])}
  if(degree==2){fp <- fracpoly[[4]]; powers <- as.numeric(fracpoly[[3]])}
  x <- runif(10000, fracpoly[[7]], fracpoly[[8]])
  if(!is.null(xlim)){x <- runif(10000, xlim[1], xlim[2])}
  if(is.null(xref)){xref <- as.numeric(mvshape$results$xbeta[1])}
  
  if(family=="gaussian"){
    if(degree==1 & powers[1]==0){
      yest <- fp$coef[1] + fp$coef[2]*log(x)
      yse <- sqrt(vcov(fp)[1,1] + (log(x))^2*vcov(fp)[2,2] + 2*(log(x))*vcov(fp)[1,2])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==1 & powers[1]!=0){
      yest <- fp$coef[1] + fp$coef[2]*x^powers[1]
      yse <- sqrt(vcov(fp)[1,1] + (x^powers[1])^2*vcov(fp)[2,2] + 2*(x^powers[1])*vcov(fp)[1,2])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==2 & powers[1]==0 & powers[2]==0){
      yest <- fp$coef[1] + fp$coef[2]*log(x) + fp$coef[3]*log(x)*log(x)
      yse <- sqrt(vcov(fp)[1,1] + (log(x))^2*vcov(fp)[2,2] + (log(x)*log(x))^2*vcov(fp)[3,3] + 2*(log(x))*vcov(fp)[1,2] + 2*(log(x)*log(x))*vcov(fp)[1,3] + 2*(log(x))*(log(x)*log(x))*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==2 & powers[1]==0 & powers[2]!=0){
      yest <- fp$coef[1] + fp$coef[2]*log(x) + fp$coef[3]*x^powers[2]
      yse <- sqrt(vcov(fp)[1,1] + (log(x))^2*vcov(fp)[2,2] + (x^powers[2])^2*vcov(fp)[3,3] + 2*(log(x))*vcov(fp)[1,2] + 2*(x^powers[2])*vcov(fp)[1,3] + 2*(log(x))*(x^powers[2])*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==2 & powers[1]!=0 & powers[2]==0){
      yest <- fp$coef[1] + fp$coef[2]*x^powers[1] + fp$coef[3]*log(x)
      yse <- sqrt(vcov(fp)[1,1] + (x^powers[1])^2*vcov(fp)[2,2] + (log(x))^2*vcov(fp)[3,3] + 2*(x^powers[1])*vcov(fp)[1,2] + 2*(log(x))*vcov(fp)[1,3] + 2*(x^powers[1])*(log(x))*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==2 & powers[1]!=0 & powers[2]!=0 & powers[1]==powers[2]){
      yest <- fp$coef[1] + fp$coef[2]*x^powers[1] + fp$coef[3]*x^powers[2]*log(x)
      yse <- sqrt(vcov(fp)[1,1] + (x^powers[1])^2*vcov(fp)[2,2] + (x^powers[2]*log(x))^2*vcov(fp)[3,3] + 2*(x^powers[1])*vcov(fp)[1,2] + 2*(x^powers[2]*log(x))*vcov(fp)[1,3] + 2*(x^powers[1])*(x^powers[2]*log(x))*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    }
    if(degree==2 & powers[1]!=0 & powers[2]!=0 & powers[1]!=powers[2]){
      yest <- fp$coef[1] + fp$coef[2]*x^powers[1] + fp$coef[3]*x^powers[2]
      yse <- sqrt(vcov(fp)[1,1] + (x^powers[1])^2*vcov(fp)[2,2] + (x^powers[2])^2*vcov(fp)[3,3] + 2*(x^powers[1])*vcov(fp)[1,2] + 2*(x^powers[2])*vcov(fp)[1,3] + 2*(x^powers[1])*(x^powers[2])*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    }
    if(findInterval(0, c(min(ylci),max(yuci)))==1){hlinep <- TRUE}else{hlinep <- FALSE}
  }
  if(family=="binomial"){
    if(degree==1 & powers[1]==0){
      yest <- fp$coef[2]*log(x) - fp$coef[2]*log(xref)
      yse <- sqrt((log(x)-log(xref))^2*vcov(fp)[2,2])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==1 & powers[1]!=0){
      yest <- fp$coef[2]*x^powers[1] - fp$coef[2]*xref^powers[1]
      yse <- sqrt((x^powers[1]-xref^powers[1])^2*vcov(fp)[2,2])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==2 & powers[1]==0 & powers[2]==0){
      yest <- fp$coef[2]*log(x) + fp$coef[3]*log(x)*log(x) - (fp$coef[2]*log(xref) + fp$coef[3]*log(xref)*log(xref))
      yse <- sqrt((log(x)-log(xref))^2*vcov(fp)[2,2] + (log(x)*log(x)-log(xref)*log(xref))^2*vcov(fp)[3,3] + 2*(log(x)-log(xref))*(log(x)*log(x)-log(xref)*log(xref))*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==2 & powers[1]==0 & powers[2]!=0){
      yest <- fp$coef[2]*log(x) + fp$coef[3]*x^powers[2] - (fp$coef[2]*log(xref) + fp$coef[3]*xref^powers[2])
      yse <- sqrt((log(x)-log(xref))^2*vcov(fp)[2,2] + (x^powers[2]-xref^powers[2])^2*vcov(fp)[3,3] + 2*(log(x)-log(xref))*(x^powers[2]-xref^powers[2])*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==2 & powers[1]!=0 & powers[2]==0){
      yest <- fp$coef[2]*x^powers[1] + fp$coef[3]*log(x) - (fp$coef[2]*xref^powers[2] + fp$coef[3]*log(xref))
      yse <- sqrt((x^powers[1]-xref^powers[1])^2*vcov(fp)[2,2] + (log(x)-log(xref))^2*vcov(fp)[3,3] + 2*(x^powers[1]-xref^powers[1])*(log(x)-log(xref))*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    } 
    if(degree==2 & powers[1]!=0 & powers[2]!=0 & powers[1]==powers[2]){
      yest <- fp$coef[2]*x^powers[1] + fp$coef[3]*x^powers[2]*log(x) - (fp$coef[2]*xref^powers[1] + fp$coef[3]*xref^powers[2]*log(xref))
      yse <- sqrt((x^powers[1]-xref^powers[1])^2*vcov(fp)[2,2] + (x^powers[2]*log(x)-xref^powers[2]*log(xref))^2*vcov(fp)[3,3] + 2*(x^powers[1]-xref^powers[1])*(x^powers[2]*log(x)-xref^powers[2]*log(xref))*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    }
    if(degree==2 & powers[1]!=0 & powers[2]!=0 & powers[1]!=powers[2]){
      yest <- fp$coef[2]*x^powers[1] + fp$coef[3]*x^powers[2] - (fp$coef[2]*xref^powers[1] + fp$coef[3]*xref^powers[2])
      yse <- sqrt((x^powers[1]-xref^powers[1])^2*vcov(fp)[2,2] + (x^powers[2]-xref^powers[2])^2*vcov(fp)[3,3] + 2*(x^powers[1]-xref^powers[1])*(x^powers[2]-xref^powers[2])*vcov(fp)[2,3])
      ylci <- yest - 1.96*yse
      yuci <- yest + 1.96*yse
    }
    hlinep <- TRUE; pref_y <- paste0("Odds ratio of ", pref_y)
  }
  data <- data.frame(x=x, yest=yest, yse=yse, ylci=ylci, yuci=yuci)
  xminp <- fracpoly[[7]]; xmaxp <- fracpoly[[8]]; yintr <- 0
  if(!is.null(xlim)){xminp <- xlim[1]; xmaxp <- xlim[2]}
  if(family=="binomial"){logy <- T}
  if(logy==T){data$yest <- exp(data$yest); data$yuci <- exp(data$yuci); data$ylci <- exp(data$ylci); yintr <- 1}
  if(logx==T){data$x <- exp(data$x); xminp <- exp(xminp); xmaxp <- exp(xmaxp)}
  data$xminp <- xminp; data$xmaxp <- xmaxp
  
  # Mvshape plot
  if(family=="gaussian"){
    xest <- as.numeric(mvshape$results$xbeta)
    yest <- as.numeric(mvshape$results$beta) + mvshape$results$beta[1]
    yest[1] <- yest[1] - mvshape$results$beta[1]
    yvar <- mvshape$varcov[1,1]
    for(i in 1:length(xest)){
      yvar <- yvar + mvshape$varcov[i,i] + 2*mvshape$varcov[1,i]
    }
    yse <- sqrt(yvar)
    ylci <- yest - 1.96*yse 
    yuci <- yest + 1.96*yse
  }
  if(family=="binomial"){
    xest <- as.numeric(mvshape$results$xbeta)
    yest <- as.numeric(mvshape$results$beta)
    ylci <- mvshape$results$beta - 1.96*as.numeric(mvshape$results$se) 
    yuci <- mvshape$results$beta + 1.96*as.numeric(mvshape$results$se)
  }  
  datamv <- data.frame(xest=xest, yest=yest, ylci=ylci, yuci=yuci) 
  if(logy==T){datamv$yest <- exp(datamv$yest); datamv$yuci <- exp(datamv$yuci); datamv$ylci <- exp(datamv$ylci)}
  if(logx==T){datamv$xest <- exp(datamv$xest)}

  # Intercept
  hline <- data.frame(x=c((xminp-0.1*xminp), (xmaxp+0.1*xmaxp)), y=c(yintr, yintr))
  
  # Figure
  figure <- ggplot(data, aes(x=x))
  if(hlinep==TRUE){figure <- figure + geom_line(aes(x=x, y=y), color=intcolour, size=0.35, data=hline)}
  figure <- figure + geom_line(aes(y=ylci), color=cicolour) + geom_line(aes(y=yuci), color=cicolour) + geom_line(aes(y=yest), color=betacolour) + theme_bw() + labs(x=pref_x,y=pref_y) + geom_errorbar(mapping=aes(x=xest, ymin=ylci, ymax=yuci), color=mvcicolour, width=0.025, data=datamv) + geom_point(aes(x=xest, y=yest), color=mvbetacolour, data=datamv, size=4, shape=15) + theme(axis.title.x = element_text(vjust=-0.5, size=14), axis.title.y = element_text(vjust=1.5, size=14), axis.text=element_text(size=12)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(expand = c(0,0), limits=c(hline$x[1], hline$x[2]))
  suppressMessages(if(!is.null(ybreaks)==TRUE & !is.null(ylim)==FALSE){figure <- figure + scale_y_continuous(breaks=ybreaks)})
  suppressMessages(if(!is.null(ybreaks)==TRUE & !is.null(ylim)==TRUE){figure <- figure + scale_y_continuous(breaks=ybreaks, limits=ylim)})
  suppressMessages(if(!is.null(ybreaks)==FALSE & !is.null(ylim)==TRUE){figure <- figure + scale_y_continuous(limits=ylim)})
  suppressMessages(if(!is.null(xbreaks)==TRUE){figure <- figure + scale_x_continuous(expand = c(0,0), breaks=xbreaks)})
  if(logx==FALSE & logy==TRUE){ybreaks <- ggplot_build(figure)$layout$panel_ranges[[1]]$y.major_source; figure <- figure + scale_y_continuous(trans="log", breaks=ybreaks)}
  if(logx==TRUE & logy==FALSE){xbreaks <- ggplot_build(figure)$layout$panel_ranges[[1]]$x.major_source; figure <- figure + scale_x_continuous(trans="log", breaks=xbreaks)}
  if(logx==TRUE & logy==TRUE){ybreaks <- ggplot_build(figure)$layout$panel_ranges[[1]]$y.major_source; figure <- figure + scale_y_continuous(trans="log", breaks=ybreaks); xbreaks <- ggplot_build(figure)$layout$panel_ranges[[1]]$x.major_source; figure <- figure + scale_x_continuous(trans="log", breaks=xbreaks)}
  
  # Return 
  return(figure)  
  
}
  
#####################################################
##### Floating point variance #####
#####################################################
  
#' floatvar
#'
#' floatvar fits the best fitting fractional polynomial of degree 1 and 2.
#' @param v variance matrix.
#' @param tol tolerance parameter.
#' @param iter.max maximum number of iterations.
#' @return List of floating point variances, error limits and divergence parameters.
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
floatvar <- function(V, tol = 0.001, iter.max = 50) {
  m <- nrow(V)
  if (!is.matrix(V) || ncol(V) != m || m == 1) 
    stop("V must be a square matrix of size 2 x 2 or more")
  evals <- eigen(V, only.values = TRUE)$values
  if (any(evals < 0)) 
    stop("V not positive definite")
  R <- V - diag(diag(V))
  V00 <- sum(R)/(m * (m - 1))
  V10 <- apply(R, 1, sum)/(m - 1)
  fv <- c(V00, V00 - 2 * V10 + diag(V))
  for (iter in 1:iter.max) {
    w <- 1/fv
    S <- sum(w)
    w1 <- w[-1]/S
    V10 <- as.vector(V %*% w1)
    V00 <- as.vector(1/S + t(w1) %*% V %*% w1)
    fv.old <- fv
    fv <- c(V00, V00 - 2 * V10 + diag(V))
    if (max(abs(fv.old - fv)/fv) < tol) 
      break
  }
  if (iter == iter.max) 
    warning("Floated variance estimates did not converge")
  Vmodel.inv <- S * (diag(w1) - w1 %*% t(w1))
  evals <- 1/(eigen(V %*% Vmodel.inv, only.values = TRUE)$values)
  divergence <- sum(1/evals - 1 + log(evals))/2
  return(list(variance = fv, error.limits = sqrt(range(evals)), divergence = divergence))
}
