primer_f="GGTGCTCAAGGCGGGACATTCGTT"
primer_r="TCGGGATTGGCCACAGCGTTGAC"
setwd("/home/mclab/R/git/M.C.Lab/seq/ch11/")

#INSTALLATION
 library(annotate)
 library(Biostrings) #DNAStringSet Object
 library(rBLAST)
 library(rMSA)
 library(devtools)
 library(magrittr)
 library(seqinr)
 library(ape) #read.dna
 library(data.table)
 library(lubridate)
 library(RCurl)
 library(magrittr)
 library(R.utils)
 library(downloader)
 library(ggplot2)
 library(gridExtra)
 primer_f2=primer_f%>%DNAString()%>%reverseComplement()%>%as.character()
 primer_r2=primer_r%>%DNAString()%>%reverseComplement()%>%as.character()
 primer=c(primer_f,primer_f2,primer_r,primer_r2)
 
# filenames= getURL("ftp://140.109.56.5/104DATA/0504/",userpwd="rm208:167cm",verbose=TRUE,ftp.use.epsv=TRUE, dirlistonly = TRUE) %>%
#   strsplit("[\\\\]|[^[:print:]]",fixed = FALSE) %>%
#   unlist() %>% (function(x){x[grep('seq',x)]})
# filepath=sprintf('ftp://rm208:167cm@140.109.56.5/104DATA/0504/%s',filenames)
# 
# for (i in 1:(length(filenames))){
#   download.file(filepath[i],
#                 paste('./',filenames[i],sep=''))
#                 }
# names=list.files()
# new_names=list.files()%>% gsub("^[0-9][0-9]\\.","",.) %>% gsub("(T7)","_R",.) %>% gsub("(SP6)","_F",.)
# 
# for (j in 1:length(names)){
#   text=readLines(names[j])
#   for (i in length(text):1){
#     text[i+1]=text[i]
#   }
#   text[1]=paste(">",new_names[j],sep='')
#   writeLines(text,new_names[j])
# }

#11114_R, 32, 2
 
db_vec= blast(db="../../db/UniVec") 
seq=readDNAStringSet(new_names)
data_seq=seq@ranges %>% data.frame()
index_r=new_names%>%grep('R',.)
seq[index_r]=seq[index_r]%>% reverseComplement()


l.vec=array(dim=c(2,2,length(new_names)),
        dimnames = list(c(1:2),c('Start.pt','End.pt'),new_names))
for (k in 1: length(new_names)){
  if(nrow(predict(db_vec,seq[k]))!=0){
    num=predict(db_vec,seq[k])%>% (function(x){x[x$Mismatches<10,]})
    qs= num$Q.start
    qe= num$Q.end
    d3= data.frame(qs=qs, qe=qe)
    d3=d3[order(d3$qs),]
    #d4=data.table(y=c(d3$qs, d3$qe), line= seq(1, length(qs),1)%>% rep(2))
    #ggplot(d4, aes(x=factor(line),y=y,group=line))+geom_point()+geom_line(aes(group=line))+scale_y_continuous(limits=c(0,1500))
    vec=matrix(ncol=2,nrow=2)
    vec[1,1]=d3[1,1]
    vec[1,2]=d3[1,2]
    for (i in 1:(nrow(d3)-1)){
      if(d3[i+1,1]<=vec[1,2]){
        if(d3[i+1,2]>vec[1,2]){
          vec[1,2]=d3[i+1,2]
        }}
      else if(vec[2,1]%>%is.na()){
        vec[2,1]=d3[i+1,1]
        vec[2,2]=d3[i+1,2]
      }else if(d3[i+1,2]>vec[2,2]){
        vec[2,2]=d3[i+1,2]
      }
    }
    l.vec[,,k]=vec
  }
}

l.seq=array(dim=c(3,2,length(new_names)),
            dimnames = list(c(1:3),c('Start.pt','End.pt'),new_names))
for (k in 1: length(new_names)){
  if(l.vec[2,1,k]%>%is.na()){
    l.seq[1,,k]=c(1,l.vec[1,1,k]-1)
    l.seq[2,,k]=c(l.vec[1,2,k]+1,data_seq[k,3])
  }else{
    l.seq[1,,k]=c(1,l.vec[1,1,k]-1)
    l.seq[2,,k]=c(l.vec[1,2,k]+1,l.vec[2,1,k]-1)
    if(l.vec[2,2,k]!=data_seq[k,3]){
      l.seq[3,,k]=c(l.vec[2,2,k]+1,data_seq[k,3])  
    }
  }
}

seq.pure=seq
for (k in 1: length(new_names)){
  if(nrow(predict(db_vec,seq[k]))!=0){
    s.seq=seq[k]%>%unlist()
    c=matrix(ncol=3,nrow=4)
    if(l.seq[3,1,k]%>%is.na()){
      for (i in 1:2){
        for (j in 1:4){
          c[j,i]=s.seq[l.seq[i,1,k]:l.seq[i,2,k]]%>%as.character()%>%grepl(primer[j],.)
        }
      }
    }else{
      for (i in 1:3){
        for (j in 1:4){
          c[j,i]=s.seq[l.seq[i,1,k]:l.seq[i,2,k]]%>%as.character()%>%grepl(primer[j],.)
        }
      }
    }
    if(!(grep('TRUE',c)/4)%>%isEmpty()){
      seq.pure[k]=seq[k]%>%
        unlist()%>%
        (function(x){x[l.seq[(grep('TRUE',c)/4) %>% ceiling(),1,k]:
                         l.seq[(grep('TRUE',c)/4) %>% ceiling(),2,k]]})%>%
        as("DNAStringSet")
    }
  }
}

set=DNAStringSet(c(seq[1:10]))
set

msa=muscle(set)


msa%>%unlist()
detail(msa)
print(msa)
sink(file = "aln.il", type = "output")
print(msa, from = 1, to = msa$length)
sink()
plot(mas,1,20)
mas

distCV(msa)


#pure too short?
#too many pri NOT FOUND?

#
testd=cbind(
  seq@ranges%>%data.frame()%>%(function(x){x[,c(4,3)]}),
  seq.pure@ranges%>%data.frame()%>%(function(x){x[,c(3)]})
)
colnames(testd)=c('names','ori','pure')
testd2=data.frame(testd,pri=(testd$ori!=testd$pure))
testd2
testd2[,4]%>%summary()
#

#READ FASTA FILES
seq= readDNAStringSet(list.files())

#LOADING DATABASE
db_nt <- blast(db="../../db/nt") 

#VECTOR SCREEN
cl <- predict(v1, seq[1]) 
cl[order(-cl$Bits),] %>% (function(x){x[1:50,]})
mas=muscle(seq)
