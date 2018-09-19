rm(list=ls())
########################NEW YORK STATE########################
### 1. Data Analysis
library(fUnitRoots)
library(quantmod)

### 1.1 Reading in Housing Data
ny=read.csv(file=file.choose(),header = T) # original NYHPI data
ny=ts(ny,start=c(1995,01),frequency = 12)

tail(ny)

#Removing last row that contains NA values
ny=ny[-265,]

#Checking for Negative Values in each column
for (i in 1:dim(ny)[2])
{
  print(sum(ny[,i]<0))
}
colnames(ny)

# "NYSLIND" column has 35 negative values. I'm removing this column 
ny=ny[,-7]

### 1.2 Calculating monthly returns of each variable
mlr=matrix(rep(0,2104),nrow=263,ncol=8)
for (i in 1:dim(ny)[2])
{
  nyr=diff(log(ny[,i]))
  mlr[,i]=nyr
}
colnames(mlr)<-c("HPr","UPr","OPr","M2r","GPr","FRr","STr","CPr")

# We have mlr, which contains our monthly return data
mlr=as.data.frame(mlr)



### 1.3 Variable Selection Process
### We are going to be using 4 different criteria for doing a Variable Selection Process

## 1.3.1 BIC, CP
# The first step is to run a "Best Subsets Regression".
# From the candidate models that this regression gives us, 
# we will use Mallow's Cp and BIC to select two candidate models.

library(ISLR)
library(leaps)
regfit.full = regsubsets(HPr ~., data=mlr, nvmax=19)
reg.summary = summary(regfit.full)
reg.summary ##
#which.max(reg.summary$adjr2) # excpet FRr
which.min(reg.summary$cp) # UPr, M2r,GPr,STr,CPr
which.min(reg.summary$bic) ## UPr and M2r

# Using Mallow's Cp, we see that the candidate model has the variables: 
# UPr, OPr, M2r, GPr, CPr

# Using BIC, we see that the candidate model has the variables:
# UPr and M2r

coef(regfit.full,5)
coef(regfit.full,2)

## 1.3.2 LASSO
# The next criteria we will use is a Lasso Regression

library(glmnet)
mlr1=data.frame(scale(mlr)) # scale data
x1=as.matrix(mlr1[,2:8])
y1=as.matrix(mlr1[,1])
lamb=seq(from=0, to=0.2, by=0.01)
fit2.lasso=glmnet(x1,y1, lambda=lamb)
plot(fit2.lasso, xvar="lambda",label=T)
set.seed(3)
#Doing a Lasso Regression with Cross-Validation
cv.lasso2=cv.glmnet(x1,y1,intercep=F,standardize=FALSE,nfolds=10)
plot(cv.lasso2)

## Finding the optimal lambda values, and the corresponding coefficients
opCoef= cv.lasso2$lambda.min
myCoefs=coef(cv.lasso2, s=opCoef)
#coef(cv.lasso2)

## Results of cv.lasso change all the time, next, we then select the best lamba with min mean error
lambdas = NULL
temp=NULL
for (i in 1:100)
{
  fit <- cv.glmnet(x1,y1)
  errors = data.frame(fit$lambda,fit$cvm)
  lambdas <- rbind(lambdas,errors)
}
# take mean cvm for each lambda
lambdas <- aggregate(lambdas[, 2], list(lambdas$fit.lambda), mean)
lambdas
# select the best one
bestindex = which(lambdas[2]==min(lambdas[2]))
bestlambda = lambdas[bestindex,1]
bestlambda
# and now run glmnet once more with it
fit <- glmnet(x1,y1,lambda=bestlambda)
coef(fit) 

# We see that the variables chosen by Lasso Regression are:
# UPr, OPr, M2r, GPr, STr, CPr

library(lars)
lc=lars(x1,y1,type="lasso")
summary(lc)
coef(lc)# all the variables are chosen


## 1.3.3 Adaptive LASSO
# This is the final criteria that we're using to select variables
# It is an adaptaion of Lasso Regression, and has the same variable selection property
library(parcor)
set.seed(1)
adalasso(x1, y1, k = 10) 

# The results are same as the Lasso Regression
# i.e, the variables chosen by Adaptive Lasso are :
# UPr, OPr, M2r, GPr, STr, CPr


## 1.4
# Here, we figure out what order of ARIMA to apply before we begin comparing the different
# candidate models
library(forecast)

# 1.4.1
# Mallow's Cp candidate model
# UPr, OPr, M2r, GPr, CPr

m1=lm(HPr~UPr+OPr+M2r+GPr+CPr,data=mlr)
acf(m1$residuals, lag.max=25)
pacf(m1$residuals) # ar 1
acf(diff(m1$residuals,12))
pacf(diff(m1$residuals,12))
auto.arima(diff(m1$residuals,12))

mts1=arima(x=mlr[,1],order=c(2,0,3),seasonal=list(order=c(1,0,1), period=12)
           ,xreg=mlr[,c(2,3,4,5,8)],method = "CSS" )
mts1
tsdiag(mts1) ##  adequate-> The residueals are NOT correlated 

# 1.4.2
# for BIC selections
#UPr and M2r
m2=lm(HPr~UPr+M2r,data=mlr)
acf(m2$residuals)
pacf(m2$residuals)
mts2=arima(x=mlr[,1],order=c(4,0,2), xreg=mlr[,c(2,4)],seasonal=list(order=c(1,0,1),period=12))
mts2
tsdiag(mts2)#

# 1.4.3
# for LASSO and Adatpive LASSO selections
# UPr, OPr, M2r, GPr, STr, CPr
m3=lm(HPr~UPr+OPr+M2r+GPr+STr+CPr,data=mlr)
acf(m3$residuals)
pacf(m3$residuals)
auto.arima(m3$residuals)

mts3=arima(x=mlr[,1],order=c(3,0,2), xreg=mlr[,c(2,3,4,5,7,8)]
           ,seasonal=list(order=c(1,0,1),period=12))
mts3
tsdiag(mts3)

## 1.5 
## Model Comparison
# 1.5.1: Mallow's CP
# UPr, OPr, M2r, GPr, CPr
pred=rep(0,62)
err=rep(0,62)
for(i in 200:262)
{
  mm=arima(x=mlr[1:i,1],order=c(2,0,3),seasonal=list(order=c(1,0,1), period=12)
            ,xreg=mlr[1:i,c(2,3,4,5,8)],method = "CSS")
  nextRow=mlr[i+1,c(2,3,4,5,8)]
  fore=forecast(mm1,1,xreg=nextRow)
  pred[i-200+1]=fore$mean
  err[i-200+1]=mlr[i+1,1]-pred[i-200+1]
} ## get the predictions and errors
sum_mabso=0
sum_rmse=0
for(i in 1:length(err))
{
  sum_mabso=sum_mabso+abs(err[i])
  sum_rmse=sum_rmse+(err[i]^2)
}

mabso=sum_mabso/length(err)
rmse=sqrt(sum_rmse/length(err))
#mabso= 0.001098447
#rmse= 0.001447683


## 1.5.2: BIC
## UPr and M2r
pred=rep(0,62)
err=rep(0,62)
for(i in 200:262)
{
  mm=arima(x=mlr[1:i,1],order=c(4,0,2),seasonal=list(order=c(1,0,1), period=12)
           ,xreg=mlr[1:i,c(2,4)],method = "CSS")
  nextRow=mlr[i+1,c(2,4)]
  fore=forecast(mm,1,xreg=nextRow)
  pred[i-200+1]=fore$mean
  err[i-200+1]=mlr[i+1,1]-pred[i-200+1]
} ## get the predictions and errors
sum_mabso=0
sum_rmse=0
for(i in 1:length(err))
{
  sum_mabso=sum_mabso+abs(err[i])
  sum_rmse=sum_rmse+(err[i]^2)
}

mabso=sum_mabso/length(err)
rmse=sqrt(sum_rmse/length(err))
#mbaso= 0.001073746
#rmse= 0.001454871

## 1.5.3: Lasso and Adaptive Lasso
# UPr, OPr, M2r, GPr, STr, CPr
pred=rep(0,62)
err=rep(0,62)
for(i in 200:262)
{
  mm=arima(x=mlr[1:i,1],order=c(3,0,2),seasonal=list(order=c(1,0,1), period=12)
           ,xreg=mlr[1:i,c(2,3,4,5,7,8)],method = "CSS")
  nextRow=mlr[i+1,c(2,3,4,5,7,8)]
  fore=forecast(mm,1,xreg=nextRow)
  pred[i-200+1]=fore$mean
  err[i-200+1]=mlr[i+1,1]-pred[i-200+1]
} 
# get the predictions and errors
sum_mabso=0
sum_rmse=0
for(i in 1:length(err))
{
  sum_mabso=sum_mabso+abs(err[i])
  sum_rmse=sum_rmse+(err[i]^2)
}

mabso=sum_mabso/length(err)
rmse=sqrt(sum_rmse/length(err))
#mbaso= 0.001128367
#rmse= 0.001479593

## We see that the candidate model from the BIC measure has the lowest Mean Absolute Error
## and Root Mean Square Error!
## Therefore the variables we will be moving ahead with is: UPr and M2r
## UPr= Unemployment Rate, M2r= M2 Money Supply(U.S.A)

### 1.6
## Fluctuation tests/ Chow Test
# examine the stability of the parameters by rolling OLS regressions ##UPr and M2r
library(zoo)

roll <- function(x) coef(lm(HPr ~ UPr+M2r,data = as.data.frame(x)))
rollregression <- rollapplyr(lny,120,roll,by.column = FALSE) #size of rolling = half of the length of time series
for (i in 1:3){
  plot(rollregression[,i])
}
# Fluctuation tests
library(strucchange)
lny=as.data.frame(lny)
HP_ocus <- efp(HPr ~ UPr+M2r, data = lny, type="OLS-CUSUM")
sctest(HP_ocus)                                       # Significance testing
bound1 <- boundary(HP_ocus, alpha=0.05)               # Boundaries and plotting
plot(HP_ocus, main = "HP_ocus Fluctuation tests")
# Tests based on F statistics (Chow Test) 
HP_fs <- Fstats(HPr ~ UPr+M2r, data = lny, from = 0.1)
sctest(HP_fs)
plot(HP_fs, main = "HP F tests")
# Dating structural changes
HP_bp <- breakpoints(HPr ~ UPr+M2r, data = lny, h = 0.1)
coef(HP_bp)
plot(lny[,1])
lines(fitted(HP_bp), col = 4)
lines(confint(HP_bp))

#From our structural change analysis, we see that a breakpoint occurs approximately
#every 3 years. We shall use this window in our Rolling Window model to make predictions

###1.7
###Rolling Window
#3-Year-Period
pred_rolln=rep(0,227)
err_rolln=rep(0,227)
for(i in 1:227)
{
  mm=arima(mlr[i:(i+35),1],order=c(4,0,2), seasonal=list(order=c(1,0,0),period=12)
           , xreg=mlr[i:(i+35),c(2,4)],method = "CSS")
  nextRow=mlr[i+35+1,c(2,4)]
  fore=forecast(mm,1,xreg=nextRow)
  pred_rolln[i]=fore$mean
  err_rolln[i]= mlr[i+35+1,1]-pred_rolln[i]
}
#Plotting Log-return Housing prices NY
plot(pred_rolln, type="l", main="LogReturn-Housing Prices NY")

#Let's plot actual Housing Index values and the Predicted values on the same plot
plot(mlr[37:263,1], type="l", main="LogReturn-Housing Prices NY", ylab="LogReturn", 
     xlab="Time in Months")
lines(pred_rolln, col="blue")

legend("topright", c("Actual", "Prediction"), lty=c(1,1), lwd=c(2,2), 
       col=c("black","blue"))

#Pretty good fit!

#Calculating RMSE
sum_rmse_rolln=0
for(i in 1:(length(err_rolln)-1))
{
  sum_rmse_rolln=sum_rmse_rolln+err_rolln[i]^2
}
rmse_roll=sqrt(sum_rmse_rolln/(length(err_rolln)-1))
paste("RMSE is ",rmse_roll)





########################CALIFORNIA STATE########################
###2.
###2.1 Reading in California Data

ca=read.csv(file=file.choose(),header = T) 
ca=ts(ca,start=c(1995,01),frequency = 12)

# The last row has NA values, therefore I'm removing it
ca=ca[-265,]
tail(ca)

#Checking for Negative Values in each column
for (i in 1:dim(ca)[2])
{
  print(sum(ca[,i]<0))
}
colnames(ca)

#CASLIND has 31 negative values, I'm removing this column
ca=ca[,-7]

### 2.2 Calculating returns of the california dataset
lca=log(ca) # log transformation
cmlr=matrix(rep(0,2104),nrow=263,ncol=8) # log return of original dataset

for (i in 1:dim(ca)[2])
{
  car=diff(log(ca[,i]))
  cmlr[,i]=car
}
colnames(cmlr)<-c("HPr","UPr","OPr","M2r","GPr","FRr","STr","CPr")
cmlr=as.data.frame(cmlr)


### 2.3 Variable Selection
## 2.3.1 BIC, CP
library(ISLR)
library(leaps)
c.regfit.full= regsubsets(HPr ~., data=cmlr1, nvmax=7)
c.reg.summary = summary(c.regfit.full)
c.reg.summary ##
which.max(c.reg.summary$adjr2)# UPr, OPr,GPr,FRr
which.min(c.reg.summary$cp)# UPr, OPr,GPr,FRr
which.min(c.reg.summary$bic)## FRr

# Using Mallow's Cp as a measure, we get the candidate model with the variables:
# UPr,OPr,GPr,FRr
# Using BIC as a measure, we get the candidate model with the variables:
# FRr

## 2.3.2 LASSO
library(glmnet)
cmlr1=data.frame(scale(cmlr))
cx1=as.matrix(cmlr1[,2:8])
cy1=as.matrix(cmlr1[,1])
lamb=seq(from=0, to=0.2, by=0.01)
c.fit.lasso=glmnet(cx1,cy1, lambda=lamb, standardize = F)
plot(c.fit.lasso,xvar="lambda",label=T)
set.seed(1)
c.cv.lasso=cv.glmnet(cx1,cy1,intercep=F,standardize=FALSE)
plot(c.cv.lasso)
opCoef=c.cv.lasso$lambda.min
opCoef
coef(c.cv.lasso,s = opCoef) # UPr,OPr,M2r,GPr,FRr are chosen


# Results of cv.lasso change all the time, 
# next, we then select the best lamba with min average error
c.lambdas = NULL
for (i in 1:100)
{
  fit <- cv.glmnet(cx1,cy1, lambda=lamb, standardize=F)
  errors = data.frame(fit$lambda,fit$cvm)
  c.lambdas <- rbind(c.lambdas,errors)
}
# take mean cvm for each lambda
c.lambdas <- aggregate(c.lambdas[, 2], list(c.lambdas$fit.lambda), mean)
# select the best one
which(c.lambdas[2]==min(c.lambdas[2]))
c.bestlambda = c.lambdas[3,1]
c.bestlambda
# and now run glmnet once more with it
c.fit2 <- glmnet(cx1,cy1,lambda=c.bestlambda)
coef(c.fit2) 

# The variables chosen from Lasso Regression are:
# UPr,OPr,M2r,GPr,FRr,CPr

#library(lars)
#lc=lars(cx1,cy1,type="lasso")
#summary(lc)
#coef(lc) #except STr


## 2.3.3 Adaptive Lasso
library(parcor)
adalasso(cx1, cy1, k = 10) # UPr,OPr,M2r,GPr,FRr

## The variables chosen by an Adaptive Lasso regression are:
## UPr,OPr,M2r,GPr,FRr

### 2.4 
## Now, we will figure out what model to apply for each of our candidate models

## 2.4.1 
## Mallow's Cp
## UPr,OPr,GPr,FRr
c.m1=lm(HPr~UPr+OPr+GPr+FRr,data=cmlr1)
acf(c.m1$residuals, lag.max=36)
pacf(c.m1$residuals, lag.max = 36) # ar 1
auto.arima(c.m1$residuals)

c.mts1=arima(x=cmlr1[,1],order=c(1,0,2),xreg= cmlr1[,c(2,3,5,6)]
             ,seasonal=list(order=c(1,0,1),period=12))
c.mts1
tsdiag(c.mts1) ##  adequate 

## 2.4.2
## BIC selections
## FRr

c.m2=lm(HPr~FRr,data=cmlr1)
acf(c.m2$residuals, lag.max=36)
pacf(c.m2$residuals, lag.max=36)
auto.arima(c.m2$residuals)

c.mts2=arima(x=cmlr1[,1],order=c(1,0,2), xreg=cmlr1[,6]
           ,seasonal=list(order=c(1,0,1),period=12))
c.mts2
tsdiag(c.mts2)#

## 2.4.3
## LASSO selections 
## UPr,OPr,M2r,GPr,FRr,CPr

c.m3=lm(HPr~UPr+OPr+M2r+GPr+FRr+CPr,data=cmlr1)
acf(c.m3$residuals, lag.max=36)
pacf(c.m3$residuals, lag.max=36)
c.mts3=arima(x=cmlr1[,1],order=c(1,0,2), xreg=cmlr1[,c(2,3,4,5,6,8)]
           ,seasonal=list(order=c(1,0,1),period=12)) ## 
c.mts3
tsdiag(c.mts3)

## 2.4.4
## Adaptive LASSO selections 
## UPr,OPr,M2r,GPr,FRr

c.m4=lm(HPr~UPr+OPr+M2r+GPr+FRr,data=cmlr1)
acf(c.m4$residuals, lag.max=36)
pacf(c.m4$residuals, lag.max=36)
c.mts4=arima(x=cmlr1[,1],order=c(1,0,2), xreg=cmlr1[,c(2,3,4,5,6)]
             ,seasonal=list(order=c(1,0,1),period=12)) ## 
c.mts4
tsdiag(c.mts4)


### 2.5
## Now, we compare the different models based on the RMSE values they give

## 2.5.1
## Mallow's Cp
## UPr,OPr,GPr,FRr

pred=rep(0,63)
err=rep(0,63)
for(i in 200:262)
{
  mm=arima(x=cmlr1[1:i,1],order=c(1,0,2), seasonal=list(order=c(1,0,1),period=12)
           , xreg=cmlr1[1:i, c(2,3,5,6)])
  nextRow=cmlr1[i+1,c(2,3,5,6)]
  fore=forecast(mm,1,xreg=nextRow)
  pred[i-200+1]=fore$mean
  err[i-200+1]=cmlr1[i+1,1]-pred[i-200+1]
}

sum_mabso=0
sum_rmse=0
for(i in 1:(length(err)))
{
  sum_mabso=sum_mabso+abs(err[i])
  sum_rmse=sum_rmse+err[i]^2
}

mabso=sum_mabso/(length(err))
rmse=sqrt(sum_rmse/(length(err)))
# mabso= 0.1443725
# rmse= 0.1859875


## 2.5.2
## BIC
## Variables: FRr

pred=rep(0,63)
err=rep(0,63)
for(i in 200:262)
{
  dat= as.data.frame(cmlr1[1:i,6])
  colnames(dat)="FRr"
  mm= arima(x=cmlr1[1:i,1],order=c(1,0,2), seasonal=list(order=c(1,0,1),period=12)
           ,xreg=dat)
  nextRow=as.data.frame(cmlr1[i+1,6])
  colnames(nextRow)="FRr"
  rownames(nextRow)=i+1
  fore=forecast(mm,1,xreg=nextRow)
  pred[i-200+1]=fore$mean
  err[i-200+1]=cmlr1[i+1,1]-pred[i-200+1]
}


sum_mabso=0
sum_rmse=0
for(i in 1:(length(err)))
{
  sum_mabso=sum_mabso+abs(err[i])
  sum_rmse=sum_rmse+err[i]^2
}

mabso=sum_mabso/(length(err))
rmse=sqrt(sum_rmse/(length(err)))
#mabso= 0.1474149
#rmse= 0.1859264


## 2.5.3
## Lasso
## UPr,OPr,M2r,GPr,FRr,CPr

pred=rep(0,63)
err=rep(0,63)
for(i in 200:262)
{
  mm=arima(x=cmlr1[1:i,1],order=c(1,0,2), seasonal=list(order=c(1,0,1),period=12)
           , xreg=cmlr1[1:i, c(2,3,4,5,6,8)])
  nextRow=cmlr1[i+1,c(2,3,4,5,6,8)]
  fore=forecast(mm,1,xreg=nextRow)
  pred[i-200+1]=fore$mean
  err[i-200+1]=cmlr1[i+1,1]-pred[i-200+1]
}

sum_mabso=0
sum_rmse=0
for(i in 1:(length(err)))
{
  sum_mabso=sum_mabso+abs(err[i])
  sum_rmse=sum_rmse+err[i]^2
}

mabso=sum_mabso/(length(err))
rmse=sqrt(sum_rmse/(length(err)))
# mabso= 0.1440821
# rmse= 0.1854916


## 2.5.4
## Adaptive Lasso
## UPr,OPr,M2r,GPr,FRr

pred=rep(0,63)
err=rep(0,63)
for(i in 200:262)
{
  mm=arima(x=cmlr1[1:i,1],order=c(1,0,2), seasonal=list(order=c(1,0,1),period=12)
           , xreg=cmlr1[1:i, c(2,3,4,5,6)])
  nextRow=cmlr1[i+1,c(2,3,4,5,6)]
  fore=forecast(mm,1,xreg=nextRow)
  pred[i-200+1]=fore$mean
  err[i-200+1]=cmlr1[i+1,1]-pred[i-200+1]
}

sum_mabso=0
sum_rmse=0
for(i in 1:(length(err)))
{
  sum_mabso=sum_mabso+abs(err[i])
  sum_rmse=sum_rmse+err[i]^2
}

mabso=sum_mabso/(length(err))
rmse=sqrt(sum_rmse/(length(err)))
# mabso= 0.1420605
# rmse= 0.1828836


## The candidate model that gives is the least values of both Mean Absolute Error and
## Root Mean Square Error is the Adaptive Lasso.

### 2.6 Structure change
### We now examine the stability of the parameters by rolling OLS regressions 
lca=log(ca)
colnames(lca)<-c("HPr","UPr","OPr","M2r","GPr","FRr","STr","CPr")
roll <- function(x) coef(lm(HPr ~ UPr+OPr+M2r+GPr+FRr,data = as.data.frame(x)))
c.rollregression <- rollapplyr(lca,120,roll,by.column = FALSE) #size of rolling = half of the length of time series
for (i in 1:5){
  plot(c.rollregression[,i])
}
# Fluctuation tests
library(strucchange)
c.HP_ocus <- efp(HPr ~ UPr+OPr+M2r+GPr+FRr, data = lca, type="OLS-CUSUM")
sctest(c.HP_ocus)                                       # Significance testing
bound1 <- boundary(c.HP_ocus, alpha=0.05)               # Boundaries and plotting
par(mfrow=c(1,1))
plot(c.HP_ocus, main = "HP_ocus Fluctuation tests")
# Tests based on F statistics (Chow Test) 
c.HP_fs <- Fstats(HPr ~ UPr+OPr+M2r+GPr+FRr, data = lca, from = 0.1)
sctest(c.HP_fs)
plot(c.HP_fs, main = "HP F tests")
# Dating structural changes
c.HP_bp <- breakpoints(HPr ~ UPr+OPr+GPr+FRr, data = lca, h = 0.1)
coef(c.HP_bp)
plot(lca[,1])
lines(fitted(c.HP_bp), col = 4)
lines(confint(c.HP_bp))

#From our structural change analysis, we see that a breakpoint occurs approximately
#every 2 and a half years. We shall use this window in our Rolling Window model 
#to make predictions

###2.7
### Now we will use the window of 2 and a half years, that we got from our
### structural change analysis to make predictions using our final candidate model
#Rolling Window
#California
#2-Year-Period

pred_rollc=rep(0,233)
err_rollc=rep(0,233)
for(i in 1:233)
{
  mm=arima(cmlr1[i:(i+29),1],order=c(1,0,2), seasonal=list(order=c(1,0,1),period=12)
           , xreg=cmlr1[i:(i+29),c(2,3,4,5,6)],method="CSS")
  nextRow=cmlr1[i+29+1,c(2,3,4,5,6)]
  fore=forecast(mm,1,xreg=nextRow)
  pred_rollc[i]=fore$mean
  err_rollc[i]=cmlr1[i+29+1,1]-pred_rollc[i]
}
# Starting predictions from the 31st to 263rd values


## Now, let's plot our predictions against the actual values of the CA index to see
## how we've done
plot(cmlr1[31:263,1], type="l", main="LogReturn-Housing Prices CA", ylim=c(-5,5), 
     ylab="(LogReturn)HPI California", xlab="Time in Months")
lines(pred_rollc,type="l",col="red")

legend("topright", c("Actual", "Prediction"), lty=c(1,1), lwd=c(2,2), col=c("black","red"))



#Calculating RMSE 
sum_rmse_rollc=0
for(i in 1:length(err_rollc))
{
  sum_rmse_rollc=sum_rmse_rollc+err_rollc[i]^2
}
rmse_rollc=sqrt(sum_rmse_rollc/length(err_rollc))
paste("RMSE is ",rmse_rollc)


