######################################################################
# Author: S.M. Vidanagamachchi
# Function: trypsinDigestion()
#           Digest a given set of proteins in a FASTA file
#
# Args:
#   fastaFile: name of the FASTA file containing protein sequences
#
# Output:
#   Digested peptide set of each protein (peptideListAll.txt)
######################################################################
trypsinDigestion <- function (fastaFile){
  
  t1<-Sys.time()
  require("Biostrings")
  mainDir <- getwd()#the directory you run the main program
  
  #Read FASTA sequences in UniRef database
  s <- readAAStringSet(file.path(mainDir,fastaFile))
  
  #Read sequence data as characters
  proteinSeq <- as.character(s)
  proteinID <- names(s)
  mPID <- as.matrix(proteinID)
  
  #split the long protein id to extract the short id
  proteinIDShort <- strsplit(mPID, "|", fixed=TRUE)
  shortID<-matrix(unlist(proteinIDShort), ncol=5, byrow=TRUE)
  proteinlist1 <- as.matrix(proteinSeq)

  #Perform trypsin digestion
  for (p in 1: length(proteinlist1)){#Go through the protein list
  indexList <- list()
  indexList[[1]]<- 0 # Start position of digestion
  z <- 0# keep last index (length of protein)
  c <- 2# keep next index
  
  #Process proteins one by one
  for (q in 1: nchar(proteinlist1[p])){
    charCurrent <- substr(proteinlist1[p],q,q)
    if (q==1){charPrevious <- ''}else{charPrevious <- substr(proteinlist1[p],q-1,q-1)}
    charNext <- substr(proteinlist1[p],q+1,q+1)
    
    #Check for digestion rules
    if (charCurrent== 'K' && charNext != 'P'){#!KP
      if ((charNext != 'Y' && charPrevious != 'C') || (charNext != 'D' && charPrevious != 'D') || (charNext != 'H' && charPrevious != 'C') || (charNext != 'D' && charPrevious != 'C')){#!CKY, !DKD, !CKH, !CKD
        indexList [[c]] <- sapply(q,as.numeric)
        c<- c+1
      }
    }else if(charCurrent == 'R' && charNext != 'P'){#!RP
      if ((charNext != 'H' && charPrevious != 'R') || (charNext != 'R' && charPrevious != 'R') || (charNext != 'K' && charPrevious != 'C') || (charNext != 'D' && charPrevious != 'D') || (charNext != 'F' && charPrevious != 'R') || (charNext != 'R' || charPrevious != 'K')){
        indexList [[c]] <- sapply(q,as.numeric)
        c<- c+1
      }   
    }else{
    }
    z<-c
  }
  
  indexList [[z]] <- nchar (proteinlist1[p])#get last index of peptide digestion index array
  
  #Write peptides into a file with the protein ID
  filename<-file.path(mainDir,paste("peptideListAll.txt",sep=""))#Set path of the peptide file
  y <-length(indexList)-1
  for (k in 1:y){
    starti<-indexList[[k]]
    stopi<-indexList[[k+1]]
    peptide <- substring(proteinlist1[p],starti+1,stopi)#get peptide sequence
    if (nchar(peptide)>2){
      dframe <- data.frame(peptide,shortID[p,4]);
      write.table(dframe, filename, sep="\t", append= TRUE,col.names=FALSE, row.names=FALSE)
    }
  }
  
  }
  peptidelist <- read.csv(filename,head=FALSE,sep="\t")#Read peptide list
  
  #return(peptidelist)
  t2<-Sys.time()
  t3<-t2-t1
  t3
}
