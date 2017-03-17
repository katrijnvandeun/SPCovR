
############################################
#             NEEDED PACKAGES              #
############################################

packages <- c("affy", "GEOquery")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(setdiff(packages, rownames(installed.packages())))  
}
library(affy)
library(GEOquery)


############################################
#             GET DATA FROM GEO            #
############################################


###2007 sample
accnr<-c("GSE29614")#enter here the accession number: GSE29614 for TIV 2007

# load series and platform data from GEO
getGEOSuppFiles(accnr)
setwd(paste("./",accnr,sep=""))
untar(paste(accnr,"_RAW.tar",sep=""), exdir="rawdata")
setwd("./rawdata")
cels <- list.celfiles()
data<-ReadAffy(filenames=cels)
#show(data): 27 samples = 9 subjects x 3 times (D0,D3,D7);


eset <- rma(data)#create RMA pre-processed expression matrix

setwd("../../")
write.exprs(eset,sep = "\t", file=("DATA/TIVRMA_2007.txt"),quote=FALSE,row.names = FALSE,col.names = FALSE)
#affyids may be checked to be the same for 2008 and 2007 because two different expression arrays are used
#c<-row.names(exprs(eset))
#write.table(c, sep = "\t", file=("GSE29614_affyID.txt"),quote=FALSE,row.names = FALSE,col.names = FALSE)


###2008 sample
accnr<-c("GSE29617")#enter here the accession number: GSE29617 for TIV 2008

# load series and platform data from GEO
getGEOSuppFiles(accnr)
setwd(paste("./",accnr,sep=""))
untar(paste(accnr,"_RAW.tar",sep=""), exdir="rawdata")
setwd("./rawdata")
cels <- list.celfiles()
data<-ReadAffy(filenames=cels)
#show(data): 80 samples = 28 subjects X 3 times (D0,D3,D7) yet not all vaccinees showed up at each of the planned days;


eset <- rma(data)

setwd("../../")
write.exprs(eset,sep = "\t", file=("DATA/TIVRMA.txt"),quote=FALSE,row.names = FALSE,col.names = FALSE)
#c<-row.names(exprs(eset))
#write.table(c, sep = "\t", file=("GSE29617_affyID.txt"),quote=FALSE,row.names = FALSE,col.names = FALSE)
