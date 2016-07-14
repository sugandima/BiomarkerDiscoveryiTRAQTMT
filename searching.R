###############################################################################
# Author: S.M. Vidanagamachchi
# Function: searching() 
# Match experimental precursor masses to the theoretical masses of the database
# (Peptide Mass Fingerprinting)
#
# Args:
#   dir: directory where the output of quantification keeps
#   experimentNo: number of experiment files (fraction output)
#   tol : mass tolerance
#   #tol<- 20
##############################################################################
searching <- function(dir,experimentNo,tol){
  
  t1<-Sys.time()
  mainDir<-getwd()
  subDir<- dir
  out1<-'matching_ratios'
  
  file_pepList <- read.csv(file.path(mainDir,paste("peptideMassModified.txt",sep="")),header=FALSE,sep="\t")#RefSeq db Digested peptides
  matrixMass <- as.matrix(file_pepList)#
  sortedmatrixMZ <- matrixMass[sort.list(as.numeric(matrixMass[,5])),]#sort the mass list
  #a1 <- sortedmatrixMZ[as.numeric(sortedmatrixMZ[,5])>350,]#3030207
  #a2 <- a1[as.numeric(a1[,5])<1795,]#2949710#matrix having peptide mz in the precursor range
  a2<-sortedmatrixMZ
  a4<-a2[duplicated(a2[,5]),]#2592503
  #remove all duplicate entries(mz) including the entry of duplicate
  a6<- a2[!(duplicated(as.numeric(a2[,5]),fromLast=TRUE) | duplicated(as.numeric(a2[,5]))),]#86143x5
  a7<- as.matrix(a6)
  dir.create(file.path(mainDir, out1))
  
  
  for(p in 1:experimentNo){
    ex1 <- read.csv(file.path(mainDir,subDir,paste(p,".txt",sep="")),header=FALSE,sep="\t")
    m_ex1 <- as.matrix(ex1) 
    low <- 1
    high <- dim(a6)[1]
    for(i in 1:((dim(m_ex1)[1]))){
      tolerance <- as.numeric(m_ex1[i,7]) *(as.numeric(tol)/1000000)
      r<-as.numeric(a7[,5])-as.numeric(m_ex1[i,7])
      a8<- as.matrix(a7[which.min(abs(r)),])#min value (selected raw outputs as one column)
      if (r<=tolerance){
        df1 <- data.frame(a8[1,1],a8[2,1],a8[3,1],a8[4,1],a8[5,1],m_ex1[i,1],m_ex1[i,16],m_ex1[i,17],m_ex1[i,18])
        write.table(df1,file.path(mainDir,out1,paste("ex",p,"_matched_peptides_",tol,"_Noise.txt",sep="")), sep="\t",append=TRUE, col.names=FALSE, row.names=FALSE) 
      }
    }
  }
  
  t2<-Sys.time()
  t3<-t2-t1
  t3
}
#1.630466 hours(chaththa machine)-bpca
#1.518879 hours-knn
