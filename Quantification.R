#########################################################################
# Author: S.M. Vidanagamachchi
# Function: Quantification()
#           Quantify patient and control peptides based on isobaric tags
#
# Args:
#   pathFiles: Directory name/file path raw spectrum data in mzXML format
#   
# Returns:
#   Intensities of patient and control sample precursor ions
#
# Prerequisites:
#   Install MSnbase package 
#   [source("http://bioconductor.org/biocLite.R");biocLite("MSnbase")]
##########################################################################
Quantification <- function(pathFiles) {
  t1<-Sys.time()
  require("MSnbase")
  mainDir <- getwd()#the directory you run the main program
  subDir <- "out"
  
  #Get names of the files in the given directory
  quantFile <- dir(file.path(mainDir,pathFiles),full.names = TRUE, pattern = ".mzXML")
  
  for (p in 1:(length(quantFile))){
    
    #Set path and name of each file
    filename<-file.path(mainDir,pathFiles,paste("ChkGnya_2plex_iTRAQ_",p,".mzXML",sep=""))
    
    #read mzXML formatted ms data
    msexp <- readMSData(filename,verbose = FALSE)
    
    #Quantifying an experiment using iTRAQ4-plex tagging (quantify reporter peaks)-quantification based on sum of intensities
    msset <- quantify(msexp, method = "sum", reporters = iTRAQ4, verbose=FALSE) 
    
    #Use impurity matrix of iTRAQ 4-plex as a correction after reporter peak quantitation (as provided by manufacturer)
    impurities <- matrix(c(0.929,0.059,0.002,0.000, 
                           0.020,0.923,0.056,0.001,
                           0.000,0.030,0.924,0.045,
                           0.000,0.001,0.040,0.923),
                         nrow=4, byrow = TRUE)
    
    #Purity correction based on above impurity matrix 
    qnt.crct <- purityCorrect(msset, impurities)
    mssetx <- qnt.crct
    
    #Remove rows having non-applicable values in all columns
    mssety<-mssetx[rowSums(is.na(mssetx[,1:2]))!=2, ]
    
    #Bayesian missing value imputation/"bpca" for values still has NA positions(without the rows that has NA values for all positions)
    qnt.imp<- impute(mssety[,1:2], method="bpca")  
    #qnt.imp<- impute(mssety[,1:2], method="knn")  
    
    #Write data of MSnSet object to a file
    dir.create(file.path(mainDir, subDir))
    write.table(t(qnt.imp), file.path(mainDir,subDir,paste(p,"data.txt",sep="")), sep="\t",col.names=NA,row.names = TRUE)
    file <- read.csv(file.path(mainDir,subDir,paste(p,"data.txt",sep="")),header=TRUE,sep="")
    m1 <- as.matrix(file)

    #Calculate normal ratios of patient/control and intensity ratios (for callculating precursor intensity)
    #of patient/control that is used in biomarker discovery
    ratio1<-as.numeric(m1[,2])/(as.numeric(m1[,2])+as.numeric(m1[,3]))
    ratio2<-as.numeric(m1[,3])/(as.numeric(m1[,2])+as.numeric(m1[,3]))
    inten1 <-matrix(nrow=dim(m1)[1], ncol=1)
    inten2 <-matrix(nrow=dim(m1)[1], ncol=1)

    inten1[,1] <-ratio1*as.numeric(m1[,8])
    inten2[,1] <-ratio2*as.numeric(m1[,8])
    
    dataOut <- data.frame (m1,inten1, inten2, as.numeric(m1[,3])/as.numeric(m1[,2]))#outputs calculated precursor ion intensities and reporter ion intensities
    
    #Write the ratios to the file 
    write.table(dataOut, file.path(mainDir,subDir,paste(p,".txt",sep="")), sep="\t", col.names=FALSE, row.names = FALSE)
  }
  t2<-Sys.time()
  t3<-t2-t1
  t3
}
