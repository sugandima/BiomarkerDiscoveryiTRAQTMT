###################################################################################
# Author: S.M. Vidanagamachchi
# Function: proteinIdentification()
#           Find number of unique proteins and calculate the protein ratio   
#           as the average of peptide abundance ratio
#
# Args:
#   experimentNo: number of experiments
#   Tol: mass tolerance of the precursor ions
# Returns:
#   Identified proteins sequences, their ratios in fractions seperately
###################################################################################
#experimentNo:22
#tol:20
###################################################################################

proteinIdentificationSeperateRuns <- function(experimentNo,tol){
#protein identification from identified peptides
t1<-Sys.time()
mainDir <- getwd()
out1<-"matching_ratios"
out2<-"ratio"

for(p in 1:experimentNo){
  
file1 <- read.csv(file.path(mainDir,out1,paste("ex",p,"_matched_peptides_",tol,"_Noise.txt",sep="")),header=FALSE,sep="\t")#RefSeq db Digested peptides
matrix_pepProteinList <- as.matrix(file1)
uniqueProteins <- unique(matrix_pepProteinList[,2])
m1<-as.matrix(uniqueProteins)
proteinRatio <-vector("list",dim(m1)[1])
df1 <-vector("list",dim(m1)[1])
peptideProteinIntensity114 <- vector("list",dim(m1)[1])#keep sum of peptide intensities of controls (114) 
peptideProteinIntensity115 <- vector("list",dim(m1)[1])#keep sum of peptide intensities of patients (115) 

#Protein ratio calculations from intensity based method
for(j in 1:(dim(m1)[1])){
  count<-0
  proteinRatio[j] <-0
  peptideProteinIntensity114[j] <-0
  peptideProteinIntensity115[j] <-0
  for(i in 1:((dim(matrix_pepProteinList)[1]))){
    if (m1[j,1]==matrix_pepProteinList[i,2]){
      peptideProteinIntensity114[j]<- as.numeric(peptideProteinIntensity114[j])+as.numeric(matrix_pepProteinList[i,7])
      peptideProteinIntensity115[j]<- as.numeric(peptideProteinIntensity115[j])+as.numeric(matrix_pepProteinList[i,8])
    }    
  }
  peptideFinal114 <- as.numeric(peptideProteinIntensity114[j])
  peptideFinal115 <- as.numeric(peptideProteinIntensity115[j])
  df1[[j]]<-data.frame(m1[j,1],peptideFinal114,peptideFinal115,peptideFinal115/peptideFinal114)
}

dir.create(file.path(mainDir, out2))

for(r in 1:(dim(m1)[1])){
  write.table(df1[[r]],file.path(mainDir,out2,paste("ex",p,"_noise",tol,"_protein_ratio.txt",sep="")), sep="\t",append=TRUE, col.names=FALSE, row.names=FALSE) 
}

}

t2<-Sys.time()
t3<-t2-t1
t3

}
