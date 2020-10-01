

rm(list=ls())     # Clean memory
graphics.off()    # Close graphs
cat("\014")       # Clear Console

#install.packages("foreign")       # Package to read .dta files.
library(foreign)


######################################################################
################ Import data and define the variables ################
######################################################################

df <- read.dta("Mroz1987.dta")    # Read the dataframe.

df <- df[1:428,]    # Keep the first 428 observations in the dataframe.

n <- nrow(df)   # n is the sample size.

y <- df$hours    # Define y as our dependent variable (wife's hours worked). It is an nx1 vector.

X <- cbind(rep(1,n), df$lwage, df$nwifeinc, df$kidslt6, df$kidsge6, df$age, df$educ)   # Define X as our nxk matrix of independent variables.
colnames(X) <- c("INTERCEPT", "LWAGE", "NONWIFEINC", "KIDSLT6", "KIDSGE6", "AGE", "EDUC")

k <- ncol(X)  # k is the number of independent variables (including the intercept).

### equation: hoursi = b1+b2lwagei+b2NonWifeInci+b3kidslt6i+b4kidsge6i+b5agei+b6educi+ei:

#################################################
################ 2SLS estimation ################
#################################################

IV <- cbind(df$age^2, df$educ^2, df$age*df$educ, df$age^3, df$educ^3, df$age^2*df$educ, df$age*df$educ^2, df$unem, df$city, df$fatheduc, df$motheduc)   # IV is our nx11 matrix of external instruments.
colnames(IV) <- c("AGE^2", "EDUC^2", "AGE*EDUC", "AGE^3", "EDUC^3", "AGE^2*EDUC", "AGE*EDUC^2", "UNEM", "CITY", "FATHEDUC", "MOTHEDUC")

Z <- cbind(X[,-2], IV)   # Z is the nxm matrix of exogenous variables (all columns of X except for lwage which corresponds to its second column) and the external instruments.

m <- ncol(Z)    # m is the number of all instruments (that is the exogenous variables in equation 1 and the external instruments).

lwage <- df$lwage    # lwage is the endogenous variable.

gamma.hat <- solve(t(Z)%*%Z)%*%(t(Z)%*%lwage)     # Regress lwage on Z and get the OLS estimates of the gamma parameters. This is the first stage.

lwage.hat <- Z%*%gamma.hat     # Get the fitted values of lwage.

X.hat <- cbind(X[,1],lwage.hat,X[,3:ncol(X)])    # Regressors for the second stage of the 2SLS: lwage is substituted with lwage.hat.

beta.hat.2sls <- solve(t(X.hat)%*%X.hat)%*%t(X.hat)%*%y  # Regress y on the X.hat and get the OLS estimates of the parameters. This is the second stage.

res.2sls <- y-(X.hat%*%beta.hat.2sls)        # Get the residuals from the 2SLS estimation. 
res.2sls <- as.vector(res.2sls)      # Convert res.2sls into a vector if R treats it as a nx1 matrix.
res.squared.2sls <- (res.2sls%*%t(res.2sls))  # Get the squared residuals

aux1 <- diag(res.squared.2sls)   # Take the diagonal entries -> these are the disturbance variances
aux2 <- matrix(0, nrow = n, ncol =n)
diag(aux2) <- aux1  # Create a nxn matrix with diagonal entries equal to aux1 and the rest equal to 0

var.hat.2sls <- sum(aux2)/(n-k)    # Get the variance of the 2SLS residuals.
var.hat.2sls <- as.numeric(var.hat.2sls)     # Convert var.hat.2sls into a number if R treats it as a 1x1 matrix.

var.cov.2sls <- var.hat.2sls*solve(t(X.hat)%*%X.hat)
se.beta.2sls <- sqrt(diag(var.cov.2sls))     # Get the standard errors of the 2SLS estimates, under homoskedasticity assumption.

white.var.cov.2sls <- solve((t(X.hat)%*%X.hat))%*%(t(X.hat)%*%aux2%*%X.hat)%*%solve((t(X.hat)%*%X.hat))
white.se.beta.2sls <- sqrt(diag(white.var.cov.2sls))   # Get the White standard errors of the estimated parameters.


########################################################################
################ Sargan (J) test for overidentification ################
########################################################################

P.Z <- Z%*%solve(t(Z)%*%Z)%*%t(Z)       # Get the projection matrix of Z.

J <- (t(res.2sls)%*%P.Z%*%res.2sls)/(t(res.2sls)%*%res.2sls/(n-k))      # Get the J statistic.(L3S47)

pchisq(J, df = ncol(Z)-k, lower.tail = F)      # compare with a Chi-squared distribution with m-k degrees of freedom, to get the p-value. 

## p-value = 0.2379 -> do not reject the null of overidentifying restrictions -> leave all the IV's in


################################################
################ GMM estimation ################
################################################

W.hat.2sls <- solve(((var.hat.2sls*t(Z))%*%(Z))/n)   # Get the two-stages least squares estimate of the weighting matrix.(L4S11)

beta.hat.gmm.2sls <- solve(t(X)%*%Z%*%W.hat.2sls%*%t(Z)%*%X)%*%(t(X)%*%Z%*%W.hat.2sls%*%t(Z)%*%y)  # Get the GMM parameter estimates.(L4S9).

g.bar.hat.2sls <- (t(Z)%*%(y-X%*%beta.hat.gmm.2sls))/n      # Get the estimated (based on beta.hat.gmm.2sls) sample average of the elementary zero function.(L4S16)

J.gmm.2sls <- n*(t(g.bar.hat.2sls)%*%W.hat.2sls%*%g.bar.hat.2sls)  # this is exactly equal to the J above (= 12.75116)

pchisq(J.gmm.2sls, df = m-k, lower.tail = F)      # compare with a Chi-squared distribution with ? degrees of freedom, to get the p-value

## the p-value equals 0.2379, so we do not reject the null of overidentifying restrictions -> leave all the IV's in

W.hat.white <- solve((t(Z)%*%aux2%*%Z)/n)  # Get the White estimate of the weighting matrix.(L4S12)

beta.hat.gmm.white <- solve(t(X)%*%Z%*%W.hat.white%*%t(Z)%*%X)%*%(t(X)%*%Z%*%W.hat.white%*%t(Z)%*%y)    # Get the GMM parameter estimates.(L4S12)

g.bar.hat.white <- (t(Z)%*%(y-X%*%beta.hat.gmm.white))/n      # Get the estimated (based on beta.hat.gmm.white) sample average of the elementary zero function.(L4S16)

J.gmm.white <- n*(t(g.bar.hat.white)%*%W.hat.white%*%g.bar.hat.white)     # Get the J statistic for overidentification, based on g.bar.hat.white.(L4S33andL4S41) 
  
pchisq(J.gmm.white, df = m-k, lower.tail = F)       # compare with a Chi-squared distribution with m-k degrees of freedom, to get the p-value.(L4S30) 

## the p-value equals 0.12269, so we do not reject the null of overidentifying restrictions -> leave all the IV's in

 