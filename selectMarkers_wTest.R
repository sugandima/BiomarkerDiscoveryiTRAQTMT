######################################################################
# Author: S.M. Vidanagamachchi
# Function: extractCommonProteins()
#           Extracting common proteins in all experiments
#
# Args:
#   matrixCombination: matrix of experiment numbers with best combination
#   tol: mass tolerance
# Returns:
#   Protein ratios and p-values 
###################################################################
#
#matrixCombination <- number of fractions
#tol<-20
###################################################################

selectMarkers_wTest <- function(matrixCombination,tol) {#parameters=22,20
  
  t1<-Sys.time()
  mainDir <- getwd()
  path1 <- "ratio"
  path2<-"allProteins_ratios"
  files<-list.files(file.path(mainDir,path1), full.names = TRUE, pattern="_protein_ratio.txt")
  df1 <-vector("list",length(files))
  
  for (p in 1:length(files)){
    filename <- paste(p, sep="")
    file<-assign(filename, file.path(mainDir,path1,paste("ex",p,"_noise",tol,"_protein_ratio.txt",sep="")))
    file1<-read.csv(file,header=FALSE,sep="\t")
    m_file<-as.matrix(file1[,1])
    df1[[p]] <- m_file
  }
  
  m_csv<-read.csv(file.path(mainDir,path2,paste("all_protein_ratioUniquePep20.txt")),header=FALSE,sep="\t")
  commonList<-as.matrix(m_csv)

  
  m_commonList <- as.matrix(commonList)[,1]
  fc1 <- as.matrix(commonList)[,2]

  length (m_commonList)
  
  output<-matrix(nrow = length(m_commonList), ncol = matrixCombination)
  output_intensity<-matrix(nrow = length(m_commonList), ncol = matrixCombination)
  output_intenP1 <- matrix(nrow = length(m_commonList), ncol = matrixCombination)
  output_intenC1 <- matrix(nrow = length(m_commonList), ncol = matrixCombination)
    
  df2 <-vector("list",length(m_commonList))
  df3 <-vector("list",length(files))
  
  #Write the fold change (patient/control) ratios of proteins in to a file
  for(k in 1:(length(m_commonList))){#
    for (p in 1:matrixCombination){#search throgh all the experiment files
      filename <- paste(p, sep="")
      file<-assign(filename, file.path(path1,paste("ex",p,"_noise",tol,"_protein_ratio.txt",sep="")))
      file1<-read.csv(file,header=FALSE,sep="\t")
      m_file1<-as.matrix(file1)
      m_file<-as.matrix(file1[,1])
      m_intenP1 <- as.matrix(file1[,2])#
      m_intenC1 <- as.matrix(file1[,3])#
      
      #m_ratio_weighted_intensity<-as.matrix(file1[,5])
      for(r in 1: length(m_file)){
        if (m_file[r]==m_commonList[k]){#Get ratio of the common element
          #output_intensity[k,p] <- as.numeric(m_ratio_weighted_intensity[r])
          output_intenP1[k,p] <- as.numeric(m_intenC1[r])
          output_intenC1[k,p] <- as.numeric(m_intenP1[r])
                    
        }
      }
      
    }
  }
  
  dx2 <- data.frame(m_commonList,output_intenP1[,1:matrixCombination],output_intenC1[,1:matrixCombination],fc1)

  write.table(dx2,file.path(mainDir,path1,paste("22ExpTestIntensitiesAllProteins_Chikungunya_",tol,".txt",sep="")), sep="\t",append=TRUE, col.names=FALSE, row.names=FALSE) 
   
  m_in1<-read.csv(file.path(mainDir,path1,paste("22ExpTestIntensitiesAllProteins_Chikungunya_",tol,".txt",sep="")),header=FALSE,sep="\t")
  m_in2<-as.matrix(m_in1)#207x46
  m_in2[is.na(m_in2)]<-0
  
  write.table(m_in2,file.path(mainDir,path1,paste("22ExpTestIntensitiesAllProteins_Chikungunya_",tol,"_1.txt",sep="")), sep="\t",append=TRUE, col.names=FALSE, row.names=FALSE) 
  m_in3<-read.csv(file.path(mainDir,path1,paste("22ExpTestIntensitiesAllProteins_Chikungunya_",tol,"_1.txt",sep="")),header=FALSE,sep="\t")
  m_in4<-as.matrix(m_in3)#
  pvalue<- matrix(nrow=length(m_in4[,1]), ncol=1)
  
  for(y in 1:length(m_commonList)){
    ta<-wilcox.test(as.numeric(m_in4[y,2:11]),as.numeric(m_in4[y,12:21]), paired = TRUE, alternative = "greater")
    pvalue[y,1]<-ta$p.value
  }
  
  dxx <- data.frame(m_commonList,fc1,pvalue)
  write.table(dxx,file.path(mainDir,path1,paste("22ExpTestIntensitiesAllProteins_Chikungunya_w_test",tol,".txt",sep="")), sep="\t",append=TRUE, col.names=FALSE, row.names=FALSE) 
  
  
 t2<-Sys.time()
 t3<-t2-t1
 t3
  
}
#4.104143 hours
#4.193821 hours-1 modification