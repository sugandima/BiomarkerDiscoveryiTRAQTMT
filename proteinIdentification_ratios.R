###################################################################################
# Author: S.M. Vidanagamachchi
# Function: proteinIdentification()
#           Find number of unique proteins and calculate the protein ratio   
#           as the average of peptide abundance ratio
#
# Args:
#   fractionNo: number of fractions in the experiments
#   tol: mass tolerance of the precursor ions
#   
# Returns:
#   Identified proteins sequences, their ratios 
###################################################################################
#fractionNo:22
#Tol:20
###################################################################################

proteinIdentification_ratios <- function(fractionNo,tol){
  t1<-Sys.time()
  #protein identification from identified peptides
  mainDir <- getwd()
  out1<-"matching_ratios"
  out2<-"allProteins_ratios"
  path2<-"allPeptides_ratios"
  df1 <-vector("list",fractionNo)
  
  #Read files of matched peptides
  for(p in 1:fractionNo){
    
    file1 <- read.csv(file.path(mainDir,out1,paste("ex",p,"_matched_peptides_",tol,"_Noise.txt",sep="")),header=FALSE,sep="\t")#RefSeq db Digested peptides
    matrix_pepProteinList <- as.matrix(file1)
    df1[[p]] <-  matrix_pepProteinList
    
  }
  
  df<- df1[[1]]
  
  #row binding
  for (p in 2:fractionNo){
    df<- rbind(df,df1[[p]])    
  }
  
  dir.create(file.path(mainDir, path2))
  
  write.table(df,file.path(mainDir,path2,paste("peptidesAll",tol,".txt",sep="")), sep="\t",append=TRUE, col.names=FALSE, row.names=FALSE) 
  
  file2 <- read.csv(file.path(mainDir,path2,paste("peptidesAll",tol,".txt",sep="")),header=FALSE,sep="\t")#read the file with all peptides
  matrix_pepProteinList <- as.matrix(file2)
  
  uniquePep<-matrix_pepProteinList[!duplicated(matrix_pepProteinList[,1]),]#13174 Unique peptides
  s1<-dim(as.matrix(uniquePep))
  
  uniquePep1 <-matrix(nrow=s1[1],ncol=4)
  
  #calculate 3 ratios of peptides (intensity based, average and median) for identical peptides from different/same proteins
  for(x in 1:s1[1]){
    count<-0
    in1<-0
    in2<-0
    avg_in1<-0
    avg_in2<-0
    pep_ratio1 <- 0
    pep_list1 <-vector("list",s1[1])
    avg_pep_ratio1 <-0
    median_pep_ratio1<- 0

    for(y in 1:dim(matrix_pepProteinList)[1]){
      if(matrix_pepProteinList[y,1]==uniquePep[x,1]){
        in1 <- as.numeric(matrix_pepProteinList[y,7]) + in1
        in2 <- as.numeric(matrix_pepProteinList[y,8]) + in2
        count <- count+1
        pep_ratio1 <- as.numeric(matrix_pepProteinList[y,9]) + pep_ratio1
        pep_list1[[count]] <- as.numeric(matrix_pepProteinList[y,9]) 

      }
      
    }
    avg_in1<- in1/count #average intensity of peptide  (control)
    avg_in2<- in2/count #average intensity of peptide  (patient)
    
    
    avg_pep_ratio1 <- pep_ratio1/count #average of the peptide ratios
    #avg_pep_ratio2 <- pep_ratio2/count #average of the peptide ratios
    
    median_pep_ratio1 <- median(unlist(pep_list1))
    #median_pep_ratio2 <- median(unlist(pep_list2))
  
    uniquePep1[x,1]<- avg_in1
    uniquePep1[x,2]<- avg_in2
    uniquePep1[x,3]<- avg_pep_ratio1
    uniquePep1[x,4]<- median_pep_ratio1

  }
  
  write.table(as.matrix(uniquePep1),file.path(mainDir,path2,paste("peptidesUnique",tol,".txt",sep="")), sep="\t",append=TRUE, col.names=FALSE, row.names=FALSE) 
  
  uniqueProteins <- unique(matrix_pepProteinList[,2])#7234 unique proteins
  m1<-as.matrix(uniqueProteins)
  
  df1 <-vector("list",dim(m1)[1])
  peptideProteinIntensity114 <- vector("list",dim(m1)[1])#keep sum of peptide intensities of patients (114) 
  peptideProteinIntensity115 <- vector("list",dim(m1)[1])#keep sum of peptide intensities of patients (115) 
  
  #calculate 3 ratios of proteins (intensity based, average and median) 
  for(j in 1:(dim(m1)[1])){
    count1<-0
    pep_ratio1 <-0
    
    avg_pro_ratio1 <- 0
    
    peptideProteinIntensity114[j] <-0
    peptideProteinIntensity115[j] <-0
    
    median_pro_ratio1<- 0
  
    pep_list1 <-vector("list",dim(m1)[1])
    
    for(i in 1:((dim(uniquePep)[1]))){
      if (m1[j,1]== uniquePep[i,2]){
        peptideProteinIntensity114[j]<- as.numeric(peptideProteinIntensity114[j])+as.numeric(uniquePep1[i,1])
        peptideProteinIntensity115[j]<- as.numeric(peptideProteinIntensity115[j])+as.numeric(uniquePep1[i,2])
        pep_ratio1 <- as.numeric(uniquePep1[i,3]) + pep_ratio1
        count1 <- count1+1
        pep_list1[[count1]] <- as.numeric(uniquePep1[i,4]) 
        
      }    
    }
    
    peptideFinal114 <- as.numeric(peptideProteinIntensity114[j])
    peptideFinal115 <- as.numeric(peptideProteinIntensity115[j])
    avg_pro_ratio1 <- pep_ratio1/count1
    median_pro_ratio1 <- median(unlist(pep_list1))
  
    df1[[j]]<-data.frame(m1[j,1],abs(as.numeric(peptideFinal115/peptideFinal114)), abs(avg_pro_ratio1), abs(median_pro_ratio1))#115/114
  }
  
  dir.create(file.path(mainDir, out2))
  
  for(r in 1:(dim(m1)[1])){
    write.table(df1[[r]],file.path(mainDir,out2,paste("all_protein_ratioUniquePep",tol,".txt",sep="")), sep="\t",append=TRUE, col.names=FALSE, row.names=FALSE) 
  }
  
  t2<-Sys.time()
  t3<-t2-t1
  t3
}
