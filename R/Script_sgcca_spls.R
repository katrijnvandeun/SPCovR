############################################
#             NEEDED PACKAGES              #
############################################

packages <- c("Biobase", "annotate", "hgu133plus2.db")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(setdiff(packages, rownames(installed.packages())))  
}
packages <- c("pls", "spcr", "RGCCA", "spls")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

library(pls)
library(spls)
library(spcr)
library(RGCCA)
library(Biobase)
library(annotate)
library("hgu133plus2.db")


#####################################
#
#      Extract data
#
#####################################

y<-read.table("../DATA/TIVtiter.txt", header = FALSE, sep = "", 
           na.strings = "NA") #TITER2008_centered.txt
y2007<-read.table("../DATA/TIVtiter2007.txt", header = FALSE, sep = "", 
           na.strings = "NA") #TITER2007_centered.txt
X<-read.table("../DATA/TIVD3_rev.txt", header = FALSE, sep = "", 
           na.strings = "NA")
X_2007<-read.table("../DATA/TIVD3_2007_rev.txt", header = FALSE, sep = "", 
           na.strings = "NA")
X <- as.matrix(X)
X<-t(X)
X<-X[,1:54675]#selection of probe sets shared between 2008 and 2007 season
X2007<-t(as.matrix(X_2007))
#pre-processing titers as described by Nakaya and centering (supposed by (s)pls and RGCCA)
Y <- as.matrix(log2(y))
Y <- Y - mean(Y)
Y2007 <- as.matrix(log2(y2007))
Y2007 <- Y2007 - mean(Y2007)
 
 #Get official gene symbols: TAKES SOME TIME!
 x <- hgu133plus2SYMBOL 
 symbolall<-Lkeys(x)
 for (i in 1:length(symbolall))
 {
   symbolall[i]<-get(symbolall[i],env=hgu133plus2SYMBOL)
 }

#####################################
#
#      SPLS (Chun & Keles, 2010)
#
#####################################

#1. different values for the sparseness penalty
penv<-seq(0,.9999,length=30)
RsqYv<-c()
RsqY2007v<-c()
nrnonzerov<-c()

for (i in 1:30){
  eta<-penv[i]
  result<-spls(X,Y,K=2,eta) 
  nrnonzerov<-cbind(nrnonzerov,length(result$A))#nr of non-zero coefficients
  n1<-predict(result,X,type="fit")
  RsqY <- cor(n1,Y)[1]^2
  RsqYv<-cbind(RsqYv,RsqY)#r²(yhat,yobs) for training data
  n3<-predict(result,X2007,type="fit")
  RsqY2007 <- cor(n3,Y2007)[1]^2
  RsqY2007v<-cbind(RsqY2007v,RsqY2007)#r²(yhat,yobs) for test data
}

RESULTspls<-rbind(nrnonzerov,RsqYv,RsqY2007v)
write.table(RESULTspls, "TABLEspls.txt", row.names=F, col.names=F,sep=" ")

#2. penalty fixed such that 418 non-zero weights result
result<-spls(X,Y,K=2,0.7875) #eta=.7875 results in 418 weights not equal to 0, closest to result SPCovR
print(result,type="fit")
summary(result)

n1<-predict(result,X,type="fit")
RsqY <- cor(n1,Y)[1]^2
RsqY#value reported in Table 1

n3<-predict(result,X2007,type="fit")
RsqY2007 <- cor(n3,Y2007)[1]^2
RsqY2007#value reporred in Table 1
expvarY2007 <- 1-(sum((n3-Y2007)^2)/sum(Y2007^2))
expvarY2007

#store official gene symbols for selected probe sets
write.table(symbolall[result$A], "../GENEIDS_spls.txt", row.names=F, col.names=F,sep=" ")


#####################################
#
#      SPCR (Kawano, 2015)
#
#####################################

spcr(X,as.vector(Y),k=2,lambda.B=6,lambda.gamma=2,scale=TRUE)
#does not work!! out-of-memory

#######################################
#
#      RGCCA (Tenenhaus et al., 2014)
#
#######################################

A = list(X,Y)#work with same scaling of data as for spcovr: X center + unit var;Y: centered
#build design matrix to obtain (sparse) pls regression of X on y
C = matrix(c(0, 1, 1, 0), 2, 2)

ssq1<-sum(X^2)
ssq2<-sum(Y^2)

#c1=c(0.038,1) results in 418 non-zero coefs
result.rgcca = sgcca(A,C, c1=c(0.038,1),ncomp=c(2,1), scheme="factorial",verbose=TRUE)#note that RGCCA scales all variables to unit var
result.rgcca$AVE
sum(colSums(result.rgcca$astar[[1]]!=0))#check nr of non-zeroes

s1<-result.rgcca$AVE[[1]]
sum(s1[[1]])#VAF by the two components in pred. data reported in Table 1
sum(s1[[2]])#VAF by the two components in outcome reported in Table 1

w1<-result.rgcca$astar[[1]]
tx<-X%*%w1#construct component scores (these are not the same as in rgcca$Y)
b2<-lm(Y~tx)#find optimal regression weights to use for out-of-sample prediction

t2007<-X2007%*%w1#construct comp. scores for external data
ypred2<-t2007%*%(b2$coefficients[2:3])#created predicted titers using regr. coefs. of training data
cor(Y2007,ypred2)^2#sq. correlation (yobs, ypred) reported in Table 1

#store official gene symbols for selected probe sets
write.table(symbolall[w1[,1]!=0], "../GENEIDS_sgcca_pc1.txt", row.names=F, col.names=F,sep=" ")
write.table(symbolall[w1[,2]!=0], "../GENEIDS_sgcca_pc2.txt", row.names=F, col.names=F,sep=" ")