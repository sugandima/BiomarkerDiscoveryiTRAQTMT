######################################################################
# Author: S.M. Vidanagamachchi
# Function: massToCharge() 
# Computes the mass and charge of a set of peptide sequences from a database
#
# Args:
#
#   pH: pH value of the medium
#   modification1: 'cystein'
#   modification2: 'lysine'
# Masses Reference:http://www.weddslist.com/ms/tables.html
#
# Returns:
#   Mass to charge ratio of the peptide for the given pH value
#   
####################################################################
#eg:
#peptideSequence <- "ACDDR"
#pH <- 2.7
####################################################################

massToCharge <- function(pH,modification1,modification2){
#Peptides from RefSeq database on 31st March 2015
t1<-Sys.time()
mainDir <- getwd()
filedup <- read.csv(file.path(mainDir,paste("peptideListAll.txt",sep="")),header=FALSE,sep="")
#calculate masses for all peptides because duplicates remove based on m/z (in next function)
file <- as.matrix(filedup)
#Monoisotopic masses
aa_G<-57.021464 #Glycine
aa_A<-71.037114 #Alanine
aa_S<-87.032028 #Serine
aa_P<-97.052764 #Proline
aa_V<-99.068414 #Valine
aa_T<-101.047678 #Threonine
aa_C<-103.009184 #Cysteine
aa_I<-113.084064 #Isoleucine
aa_L<-113.084064 #Leucine
aa_N<-114.042927 #Aspagine
aa_D<-115.026943 #Aspartic acid
aa_Q<-128.058578 #Glutamine
aa_K<-128.094963 #Lysine
aa_E<-129.042593 #Glutamic acid
aa_M<-131.040485 #Methionine
aa_H<-137.058912 #Histidine
aa_F<-147.068414 #Phenylalanine
aa_R<-156.101111 #Arginine
aa_Y<-163.063329 #Tyrosine
aa_W<-186.079313 #Tryptophan
c_terminus <- 17.00734 #c-terminus
n_terminus <- 1.00794 #n-terminus
MH<- 1#A charged peptide, corresponding to the mass of a hydrogen atom
#Pk values of c terminal, n terminal, negatively charged, positively charged and polar uncharged amini acid groups
n_pk <- 9.0 #N-terminal
c_pk <- 3.5 #C-terminal
aa_D_pk <-3.9 #Aspartic acid
aa_E_pk <- 4.1 #Glutamic acid
aa_H_pk <- 6.0 #Histidine
aa_C_pk <- 8.4 #Cysteine
aa_Y_pk <- 10.5 #Tyrosine
aa_K_pk <- 10.5 #Lysine
aa_R_pk <- 12.5 #Arginine

#####Static Modifications#####
cystein <- 57.02146 #Cystein residue modification-carbamidomethylation (modification1)
n_term_mod <- 144.102063 #N- terminal modification (modification2)
lysine <- 144.102063 #Lysine modification (modification3)
  
seq_matrix <- as.matrix(file[,1])

#modification1 <- 'cystein'
#modification2 <- 'lysine'
sum_weight <- 0
sum_charge <- 0
for(i in 1:length(seq_matrix)){
  sum_weight <- 0
  sum_charge <- 0
  for (j in 1:nchar(seq_matrix[i,1])){#iterate each character
    charCurrent <- substr(seq_matrix[i,1],j,j)
    switch(charCurrent, 
           'G'={
             sum_weight <- sum_weight+aa_G
           },
           'A'={
             sum_weight <- sum_weight+aa_A
           },
           'S'={
             sum_weight <- sum_weight+aa_S
           },
           'P'={
             sum_weight <- sum_weight+aa_P
           },
           'V'={
             sum_weight <- sum_weight+aa_V
           },
           'T'={
             sum_weight <- sum_weight+aa_T
           },
           'C'={
             if (modification1 == 'cystein'){
               sum_weight <- sum_weight+aa_C+cystein
             }
             else{
               sum_weight <- sum_weight+aa_C
             }
             if (pH>aa_C_pk){
               sum_charge <- sum_charge-1
             }
           },
           'I'={
             sum_weight <- sum_weight+aa_I
           },
           'L'={
             sum_weight <- sum_weight+aa_L
           },
           'N'={
             sum_weight <- sum_weight+aa_N
           },
           'D'={
             sum_weight <- sum_weight+aa_D
             if (pH>aa_D_pk){
               sum_charge <- sum_charge-1
             }
           },
           'Q'={
             sum_weight <- sum_weight+aa_Q
           },
           'K'={
             if (modification2 == 'lysine'){
               sum_weight <- sum_weight+aa_K+lysine
             }
             else{
               sum_weight <- sum_weight+aa_K
             }
             
             if(pH<aa_K_pk){
               sum_charge <- sum_charge+1
             }
           },
           'E'={
             sum_weight <- sum_weight+aa_E
             if (pH>aa_E_pk){
               sum_charge <- sum_charge-1
             }
           },
           'M'={
             sum_weight <- sum_weight+aa_M
           },
           'H'={
             sum_weight <- sum_weight+aa_H
             if (pH<aa_H_pk){
               sum_charge <- sum_charge+1
             }
           },
           'F'={
             sum_weight <- sum_weight+aa_F
           },
           'R'={
             sum_weight <- sum_weight+aa_R
             if (pH<aa_R_pk){
               sum_charge <- sum_charge+1
             }
           },
           'Y'={
             sum_weight <- sum_weight+aa_Y
             if (pH>aa_C_pk){
               sum_charge <- sum_charge-1
             }
           },
           'W'={
             sum_weight <- sum_weight+aa_W
           },
{#default
}         
    )  
    }
    
  #}
  
  #charge n-terminus
  if (pH>n_pk){
    sum_charge <- sum_charge+0
  }else{
    sum_charge <- sum_charge+1
  }
  
  #charge c-terminus
  if (pH>c_pk){
    sum_charge <- sum_charge-1
  }else{
    sum_charge <- sum_charge+0
  }

  sum_weight <- sum_weight+n_terminus+c_terminus+(MH*sum_charge)+n_term_mod;#n terminus weight, c terminus weight and MH+ weight
  modifiedMZ = sum_weight/sum_charge;
  dfrm1 <- data.frame(seq_matrix[i,1],file[i,2],sum_weight,sum_charge,modifiedMZ)
  filename<-file.path(mainDir,paste("peptideMassModified.txt",sep=""))#Set name of the output file
  write.table(dfrm1, filename, sep="\t",append=TRUE,  col.names=FALSE, row.names=FALSE)#
    
}
t2<-Sys.time()
t3<- t2-t1
t3
}
#5.915355 hours

