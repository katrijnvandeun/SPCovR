#sparse PCA with missing data using PMA package of Witten et al.

setwd('C:/Users/Katrijn/OneDrive/Documents/Manuscripts/InProgress/SPCovR/REVISION')
#library(PMA)  #not, is penalty on loadings
library(elasticnet)
library(RGCCA)
library(psych)

RESULT<-c()
#RESULT=read.table("../result_spcr_spls.txt",sep="\t")
#for (teller in ((dim(result)[1]+1):(20*27)))
for (teller in (1:(20*27)))
{
  #read data
  X<-read.table(paste('SIMDATA/X',teller,'.txt',sep=""),sep="\t")
	X<-scale(X)
	Xout<-read.table(paste('SIMDATA/Xout',teller,'.txt',sep=""),sep="\t")
	Xout<-scale(Xout)
	y<-read.table(paste('SIMDATA/Y',teller,'.txt',sep=""),sep="\t")
	y<-scale(y)
	yout<-read.table(paste('SIMDATA/Yout',teller,'.txt',sep=""),sep="\t")
	yout<-scale(yout)
	WTRUE<-read.table(paste('SIMDATA/WTRUE',teller,'.txt',sep=""),sep="\t")
	nrnonzero = colSums(WTRUE!=0)
	PTRUE<-read.table(paste('SIMDATA/PTRUE',teller,'.txt',sep=""),sep="\t")
	d<-dim(PTRUE)
	R<-d[2]
	
	###SPARSE PCR analysis using elastic net package
	resspca=spca(X,R,para=nrnonzero,sparse="varnum")
	cbind(resspca$loadings,WTRUE)
	sum(colSums(abs((WTRUE==0)-(resspca$loadings==0))))
	
	##PERFORMANCE MEASURES
	#nr false positives
	spcr_nfalsepos=min(sum(resspca$loadings[which(WTRUE==0,arr.ind=TRUE)]!=0),sum(resspca$loadings[which(WTRUE[,c(1,2)]==0,arr.ind=TRUE)]!=0))
	#nr false negatives
	spcr_nfalseneg=min(sum(resspca$loadings[which(WTRUE!=0,arr.ind=TRUE)]==0),sum(resspca$loadings[which(WTRUE[,c(1,2)]!=0,arr.ind=TRUE)]==0))
	#Tucker congruence
	T=as.matrix(X)%*%resspca$loadings
	TTRUE<-as.matrix(read.table(paste('SIMDATA/TTRUE',teller,'.txt',sep=""),sep="\t"))
	congrmat=factor.congruence(TTRUE,T)
	spcr_tucongr=max(abs(congrmat[1,1])+abs(congrmat[2,2]),abs(congrmat[2,1])+abs(congrmat[1,2]))/2
	spcr_tucongr
	#predictive accuracy
	b=lm(y~T+0)
	Tnew=as.matrix(Xout)%*%resspca$loadings
	yhat=Tnew%*%b$coefficients
	dev=yhat-yout
	spcr_prederr=sum(dev^2)/sum(yout^2)

	###SPARSE PLS
	A = list(X,y)#work with same scaling of data as for spcovr: X center + unit var;Y: centered
	#build design matrix to obtain (sparse) pls regression of X on y
	C = matrix(c(0, 1, 1, 0), 2, 2)
	
	ssq1<-sum(X^2)
	ssq2<-sum(y^2)
	
	#NOTE: a are outer weights; astars basically to - difference is wrt deflation so on 1st comp a=astar
	result.rgcca = sgcca(A,C, c1=c(sqrt(sum(colSums(WTRUE!=0))/(2*prod(dim(WTRUE)))),1),ncomp=c(2,1), scheme="factorial",verbose=TRUE)#note that RGCCA scales all variables to unit var
	result.rgcca$AVE#AVE inner: cor(result.rgcca$Y[[1]],result.rgcca$Y[[2]])^2; outer is average sq. correl of vars with components
	sum(colSums(result.rgcca$astar[[1]]!=0))#check nr of non-zeroes
	w1<-result.rgcca$astar[[1]]
	#nr false positives
	spls_nfalsepos=min(sum(w1[which(WTRUE==0,arr.ind=TRUE)]!=0),sum(w1[which(WTRUE[,c(1,2)]==0,arr.ind=TRUE)]!=0))
	#nr false negatives
	spls_nfalseneg=min(sum(w1[which(WTRUE!=0,arr.ind=TRUE)]==0),sum(w1[which(WTRUE[,c(1,2)]!=0,arr.ind=TRUE)]==0))
	tx<-X%*%w1#construct component scores (these are not the same as in rgcca$Y)
	b2<-lm(y~tx+0)#find optimal regression weights to use for out-of-sample prediction
	congrmat=factor.congruence(TTRUE,tx)
	spls_tucongr=max(abs(congrmat[1,1])+abs(congrmat[2,2]),abs(congrmat[2,1])+abs(congrmat[1,2]))/2
	tout<-Xout%*%w1#construct comp. scores for external data
	ypred<-tout%*%(b2$coefficients)#created predicted titers using regr. coefs. of training data
	spls_prederr=sum((ypred-yout)^2)/sum(yout^2)#sq. correlation (yobs, ypred) reported in Table 1
	RESULT<-rbind(RESULT,c(spcr_nfalsepos,spcr_nfalseneg,spcr_tucongr,spcr_prederr,spls_nfalseneg,spls_nfalsepos,spls_tucongr,spls_prederr))
	write.table(RESULT,file="../result_spcr_spls.txt",sep="\t",eol="\n")
}

#load spcovr results
a=read.table("PERFORMANCE_a.txt",sep="\t")
b=read.table("PERFORMANCE_b.txt",sep="\t")
c=read.table("PERFORMANCE_c.txt",sep="\t")
spcovr=rbind(a,b,c)

###CHECK
a2=read.table("SIMDATA/PERFORMANCE_corrected_a.txt",sep="\t")


#make box plot; first define factors
errX<-rep(c("0.01", "0.40", "0.70"),c(180,180,180))
errY<-rep(c(rep(c("0.80", "0.50", "0.02"),c(60,60,60))),3)
dominance<-rep(rep(c("0.10/0.90","0.50/0.50","0.90/0.10"),c(20,20,20)),9)

#make data for spca, spcovr (a=.99) and spls, as required by the boxplot function
tucker=c(RESULT[,3],RESULT[,7],spcovr[,3])
#tucker[(540+1):(540+180)]=a2[,3]
prederr=c(RESULT[,4],RESULT[,8],spcovr[,4])
#prederr[(540+1):(540+180)]=a2[,4]

method=rep(c("spcr","spls","spcovr"),c(540,540,540))
errX=rep(errX,3)
dominance=rep(dominance,3)
errY=rep(errY,3)

pdf("boxplot.pdf",width=24,height = 8,paper="special")
par(mfrow=c(1,3))

a=boxplot(tucker[errY=="0.02"]~method[errY=="0.02"]+dominance[errY=="0.02"]+errX[errY=="0.02"],at=c(0.1,0.7,1.3, 3.1,3.7,4.3, 6.1,6.7,7.3,  9.1,9.7,10.3, 12.1,12.7,13.3, 15.1,15.7,16.3,  18.1,18.7,19.3, 21.1,21.7,22.3, 24.1,24.7,25.3),
          col=c(rep(c(1,2,3),9)),xaxt="n",ylab="Tucker Congruence",main="VAFY=0.02",ylim=c(0,1))
axis(1, at = c(3.6 , 12.6 , 21.6), labels = c("VAFX=0.01","VAFX=0.40","VAFX=0.70") , tick=FALSE , cex=0.3)
abline(v=c(8,17),lty=1, col="grey")
# Add a legend
legend("bottomright", legend = c("SPCOVR","SPCR", "SPLS"), col=c(1 , 2, 3),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = 0.01)

a=boxplot(tucker[errY=="0.50"]~method[errY=="0.50"]+dominance[errY=="0.50"]+errX[errY=="0.50"],at=c(0.1,0.7,1.3, 3.1,3.7,4.3, 6.1,6.7,7.3,  9.1,9.7,10.3, 12.1,12.7,13.3, 15.1,15.7,16.3,  18.1,18.7,19.3, 21.1,21.7,22.3, 24.1,24.7,25.3),
          col=c(rep(c(1,2,3),9)),xaxt="n",ylab="Tucker Congruence",main="VAFY=0.50",ylim=c(0,1))
axis(1, at = c(3.6 , 12.6 , 21.6), labels = c("VAFX=0.01","VAFX=0.40","VAFX=0.70") , tick=FALSE , cex=0.3)
abline(v=c(8,17),lty=1, col="grey")
# Add a legend
legend("bottomright", legend = c("SPCOVR","SPCR", "SPLS"), col=c(1 , 2,3),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = 0.01)

a=boxplot(tucker[errY=="0.80"]~method[errY=="0.80"]+dominance[errY=="0.80"]+errX[errY=="0.80"],at=c(0.1,0.7,1.3, 3.1,3.7,4.3, 6.1,6.7,7.3,  9.1,9.7,10.3, 12.1,12.7,13.3, 15.1,15.7,16.3,  18.1,18.7,19.3, 21.1,21.7,22.3, 24.1,24.7,25.3),
          col=c(rep(c(1,2,3),9)),xaxt="n",ylab="Tucker Congruence",main="VAFY=0.80",ylim=c(0,1))
axis(1, at = c(3.6 , 12.6 , 21.6), labels = c("VAFX=0.01","VAFX=0.40","VAFX=0.70") , tick=FALSE , cex=0.3)
abline(v=c(8,17),lty=1, col="grey")
# Add a legend
legend("bottomright", legend = c("SPCOVR","SPCR", "SPLS"), col=c(1 , 2,3),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = 0.01)
dev.off()


#######################
pdf("boxplotMSE.pdf",width=24,height = 8,paper="special")
par(mfrow=c(1,3))

a=boxplot(prederr[errY=="0.02"]~method[errY=="0.02"]+dominance[errY=="0.02"]+errX[errY=="0.02"],at=c(0.1,0.7,1.3, 3.1,3.7,4.3, 6.1,6.7,7.3,  9.1,9.7,10.3, 12.1,12.7,13.3, 15.1,15.7,16.3,  18.1,18.7,19.3, 21.1,21.7,22.3, 24.1,24.7,25.3),
          col=c(rep(c(1,2,3),9)),xaxt="n",ylab="Prediction error",main="VAFY=0.02",ylim=c(0,2))
axis(1, at = c(3.6 , 12.6 , 21.6), labels = c("VAFX=0.01","VAFX=0.40","VAFX=0.70") , tick=FALSE , cex=0.3)
abline(v=c(8,17),lty=1, col="grey")
# Add a legend
legend("topright", legend = c("SPCOVR","SPCR", "SPLS"), col=c(1 , 2, 3),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = 0.01)

a=boxplot(prederr[errY=="0.50"]~method[errY=="0.50"]+dominance[errY=="0.50"]+errX[errY=="0.50"],at=c(0.1,0.7,1.3, 3.1,3.7,4.3, 6.1,6.7,7.3,  9.1,9.7,10.3, 12.1,12.7,13.3, 15.1,15.7,16.3,  18.1,18.7,19.3, 21.1,21.7,22.3, 24.1,24.7,25.3),
          col=c(rep(c(1,2,3),9)),xaxt="n",ylab="Prediction error",main="VAFY=0.50",ylim=c(0,2))
axis(1, at = c(3.6 , 12.6 , 21.6), labels = c("VAFX=0.01","VAFX=0.40","VAFX=0.70") , tick=FALSE , cex=0.3)
abline(v=c(8,17),lty=1, col="grey")
# Add a legend
legend("topright", legend = c("SPCOVR","SPCR", "SPLS"), col=c(1 , 2,3),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = 0.01)

a=boxplot(prederr[errY=="0.80"]~method[errY=="0.80"]+dominance[errY=="0.80"]+errX[errY=="0.80"],at=c(0.1,0.7,1.3, 3.1,3.7,4.3, 6.1,6.7,7.3,  9.1,9.7,10.3, 12.1,12.7,13.3, 15.1,15.7,16.3,  18.1,18.7,19.3, 21.1,21.7,22.3, 24.1,24.7,25.3),
          col=c(rep(c(1,2,3),9)),xaxt="n",ylab="Prediction error",main="VAFY=0.80",ylim=c(0,2))
axis(1, at = c(3.6 , 12.6 , 21.6), labels = c("VAFX=0.01","VAFX=0.40","VAFX=0.70") , tick=FALSE , cex=0.3)
abline(v=c(8,17),lty=1, col="grey")
# Add a legend
legend("topright", legend = c("SPCOVR","SPCR", "SPLS"), col=c(1 , 2,3),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = 0.01)
dev.off()

