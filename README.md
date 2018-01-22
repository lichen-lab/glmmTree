# glmmTree Tutorial

## Installation
```R
install.packages('rrBLUP')
```
```R
install.packages('Matrix')
```
```R
nlme: https://cran.r-project.org/src/contrib/Archive/nlme/nlme_3.1-122.tar.gz
R CMD INSTALL nlme_3.1-122.tar.gz
```

## example data
### data1.rda
* z.tr: 778 OTU proportions for 100 samples in training set
* y.tr: Continous outcome for 100 samples in training set
* z.te: 778 OTU proportions for 200 samples in testing set
* y.te: Continous outcome for 200 samples in training set
* D: Patristic distance matrix among 778 OTUs

### data2.rda
* z.tr: 778 OTU proportions for 100 samples in training set
* y.tr: Binary outcome for 100 samples in training set
* z.te: 778 OTU proportions for 200 samples in testing set
* y.te: Binary outcome for 200 samples in training set
* D: Patristic distance matrix among 778 OTUs


## example for continous outcome
```R
library(rrBLUP)
library(Matrix)
source('lib.R')
source('glmmTreeg.R')

load('data1.rda')
lambda1=c(0.1,0.5)
lambda2=c(0.1,10)
obj.cv=cv.glmmTreeg(y=y.tr,Z=z.tr,X=NULL,D,lambda1=lambda1,lambda2=lambda2)
yhat=predict.cv.glmmTreeg(obj.cv,X=NULL,z.te)
plot(yhat,y.te,main='Continuous outcome')
```

## example for binary outcome
```R
library(Matrix)
source('lib.R')
source('glmmTreeb.R')
load("data2.rda")
obj.cv=cv.glmmTreeb(y=y.tr,Z=z.tr,X=NULL,D,lambda1=c(1,2),lambda2=c(1,2))
yhat=predict.cv.glmmTreeb(obj.cv,X=NULL,z.te)
boxplot(yhat ~ y.te, main='Binary outcome')
```



