primer_f="GGTGCTCAAGGCGGGACATTCGTT"
primer_r="TCGGGATTGGCCACAGCGTTGAC"
name_primer_f='SP6'
name_primer_r='T7'
source="ftp://140.109.56.5/104DATA/0504/"
username="rm208"
password="167cm"
desdir="/home/mclab/R/git/M.C.Lab"

setwd(desdir)
source('./Vec_Aln_function.R') 

VA(primer_f, primer_r, name_primer_f, name_primer_r, source, username, password, desdir)


#pure too short?
#too many pri NOT FOUND?
data_seq
data_seq[,4]%>%summary()
#READ FASTA FILES
seq= readDNAStringSet(list.files())
#LOADING DATABASE
db_nt <- blast(db="../../db/nt") 
#VECTOR SCREEN
cl <- predict(v1, seq[1]) 
cl[order(-cl$Bits),] %>% (function(x){x[1:50,]})
