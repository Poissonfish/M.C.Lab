#INSTALLATION
  # install.packages('devtools')
  # source("https://bioconductor.org/biocLite.R")
  # biocLite("annotate")
  # biocLite("Biostrings")
  # install_github("mhahsler/rBLAST")
  library(annotate)
  library(Biostrings) #DNAStringSet Object
  library(rBLAST)
  #install_github("mhahsler/rMSA")
  library(rMSA)
  library(devtools)
  library(magrittr)
  library(seqinr)
  library(ape) #read.dna
  library(ggplot2)
  library(data.table)
  # install.packages('ghostscript')
  library(lubridate)
  
  wd="/home/mclab/R/git/M.C.Lab"

  #READ FASTA FILES
  setwd(paste(wd,'/seq',sep=''))
  seq= readDNAStringSet('2A6F.txt')
  seqr= readDNAStringSet('2A6R.txt')

#LOADING DATABASE
  setwd(wd)
  v1 <- blast(db="./db/nt") 

#VECTOR SCREEN
  cl <- predict(v1, seq) %>% (function(x){x[x$Perc.Ident>=97.5,]})

#PLOTTING
  qs= cl$Q.start
  qe= cl$Q.end
  data=data.table(y= c(qs,qe), line= seq(1,length(qs),1) %>% rep(2))
  ggplot(data, aes(x=factor(line),y=y))+geom_point()+geom_line(aes(group=factor(line)))

#ALIGNMENT
  setwd(paste(wd,'/seq',sep=''))
  seqset=c(
    '2A6F.txt',
    '2A6R.txt',
    '1A1F.txt',
    '1A6F.txt'
  )
  set= readDNAStringSet('seq.txt')
  mas=muscle(set)
  clustal_help()
  
  detail(mas)
  plot(mas,1,40)
  distCV(mas)

  distCV(mas) %>% 
    hclust() %>% 
    as.dendrogram() %>% 
    plot(horiz=TRUE, type='triangle')
  