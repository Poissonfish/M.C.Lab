#INSTALLATION
library(annotate) 
library(Biostrings) #DNAStringSet Object
library(rBLAST)
library(rMSA)
library(devtools)
library(magrittr)
library(seqinr)
library(ape) #read.dna, write.fasta
library(data.table)
library(lubridate)
library(RCurl)
library(magrittr)
library(R.utils)
library(downloader)
library(ggplot2)
library(gridExtra)
library(plyr)
library(taxize)
library(rentrez)

MCLab=function(primer_f, primer_r, name_primer_f, name_primer_r, source, username, password, desdir,folder,local_path, local, nt_search){
  packageStartupMessage("Creating Folders...", appendLF = FALSE)
  setwd(desdir)
  if (!file.exists(file.path(desdir,'seq',folder,'raw'))){
    dir.create(file.path(desdir,'seq',folder,'raw'),recursive = TRUE)
  }
  if (!file.exists(file.path(desdir,'seq',folder,'fasta'))){
    dir.create(file.path(desdir,'seq',folder,'fasta'),recursive = TRUE)
  }
  if (!file.exists(file.path(desdir,'seq',folder,'fasta.aln'))){
    dir.create(file.path(desdir,'seq',folder,'fasta.aln'),recursive = TRUE)
  }
  if(local){system(paste('cp -r',file.path(local_path,'.'),file.path(desdir,'seq',folder,'raw')))}
  packageStartupMessage(" Done!")

  packageStartupMessage("Primer Analyzing...", appendLF = FALSE)
  primer_f2=primer_f%>%DNAString()%>%reverseComplement()%>%as.character()
  primer_r2=primer_r%>%DNAString()%>%reverseComplement()%>%as.character()
  primerall=c(primer_f,primer_f2,primer_r,primer_r2)
  
  pro=c('N','R','Y','K','M','W')
  sp_N=c('A','T','C','G')
  sp_R=c('A','G')
  sp_Y=c('T','C')
  sp_K=c('T','G')
  sp_M=c('A','C')
  sp_W=c('A','T')
  sp_list=list(sp_N,sp_R,sp_Y,sp_K,sp_M,sp_W)
  names(sp_list)=pro
  
  primer=c()
  for (pri in 1:4){
    listP=list()
    for (i in 1:6){
      listP[[i]]=strsplit(primerall[pri],'')%>%unlist() %>%grep(pro[i],.)
    }
    seqn=c()
    for (i in 1:6){
      seqn[listP[[i]]]=pro[i]
    }
    seqn=seqn%>%na.omit()%>%as.character()
    grid=  expand.grid(if(!is.na(seqn[1])){sp_list[seqn[1]]%>%unlist()}else{NA},
                       if(!is.na(seqn[2])){sp_list[seqn[2]]%>%unlist()}else{NA},
                       if(!is.na(seqn[3])){sp_list[seqn[3]]%>%unlist()}else{NA},
                       if(!is.na(seqn[4])){sp_list[seqn[4]]%>%unlist()}else{NA},
                       if(!is.na(seqn[5])){sp_list[seqn[5]]%>%unlist()}else{NA},
                       if(!is.na(seqn[6])){sp_list[seqn[6]]%>%unlist()}else{NA},
                       if(!is.na(seqn[7])){sp_list[seqn[7]]%>%unlist()}else{NA})
    primer=c(primer, primerall[pri]%>% gsub('[NRYKMW]','%s',.) %>% 
               sprintf(.,grid[,1],grid[,2],grid[,3],grid[,4],grid[,5],grid[,6],grid[,7]))
  }
  primerraw=primer
  primern=c()
  for (i in 1:length(primer)){
    primern=c(primern,primer[i] %>% substr(.,1,nchar(.)-5),primer[i] %>% substr(.,6,nchar(.))) 
  }
  primer=primern
  packageStartupMessage(" Done!")
  
  setwd(file.path(desdir,'seq',folder,'raw'))
  if(!local){
    packageStartupMessage("File Downloading...", appendLF = FALSE)
    filenames= getURL(source,userpwd=paste0(username,':',password),
                      verbose=TRUE,ftp.use.epsv=TRUE, dirlistonly = TRUE) %>%
               strsplit("[\\\\]|[^[:print:]]",fixed = FALSE) %>%
               unlist() %>% (function(x){x[grep('seq',x)]})
    filepath= sprintf(paste0('ftp://',
                            paste0(username,':',password),'@',
                            (source%>%gsub('ftp://','',.)),'%s'),filenames)
    
    for (i in 1:(length(filenames))){
      download.file(filepath[i],
                    file.path(getwd(),filenames[i]))
    }
    packageStartupMessage(" Done!") 
  }
  
  packageStartupMessage("Renaming...", appendLF = FALSE)  
  names=list.files()
  new_names=names
  new_names[grep(name_primer_r,names)]= names[grep(name_primer_r,names)] %>% 
                                        gsub("^[0-9][0-9]\\.","",.) %>%
                                        gsub(paste0('(',name_primer_r,')'),"_R",.)
  new_names[grep(name_primer_f,names)]= names[grep(name_primer_f,names)] %>% 
                                        gsub("^[0-9][0-9]\\.","",.) %>%
                                        gsub(paste0('(',name_primer_f,')'),"_F",.)
   
  new_names_split= new_names %>% strsplit('_') %>% unlist() 
  SN= new_names_split[seq(1,length(new_names_split),2)]  
  FR= new_names_split[seq(2,length(new_names_split),2)]  
  nchr= SN %>%  nchar() %>% (function(x){c(min(x),max(x))})
  
  if (nchr[1]!=nchr[2]){
    s_index= SN%>%nchar() == nchr[1]
    l_index= SN%>%nchar() == nchr[2]
    
    new_names[l_index]=paste0(SN[l_index] %>% substr(1,nchr[2]-3),'-',
                              SN[l_index] %>% substr(nchr[2]-2,nchr[2]-2),'-',
                              SN[l_index] %>% substr(nchr[2]-1,nchr[2]),'_', FR[l_index])
    new_names[s_index]=paste0(SN[s_index] %>% substr(1,nchr[1]-2),'-',
                              SN[s_index] %>% substr(nchr[1]-1,nchr[1]-1),'-',
                              SN[s_index] %>% substr(nchr[1],nchr[1]),'_', FR[s_index])
    new_names=new_names%>%gsub('.seq','.fasta',.)
    
    for (j in 1:length(names)){
      text=readLines(names[j])
      for (i in length(text):1){
        text[i+1]=text[i]
      }
      text[1]=paste0(">",new_names[j])
      writeLines(text,file.path('../fasta',new_names[j]%>%gsub('seq','fasta',.)))
    }
  }else{
    new_names=paste0(SN %>% substr(1,nchr[2]-3),'-',
                     SN %>% substr(nchr[2]-2,nchr[2]-2),'-',
                     SN %>% substr(nchr[2]-1,nchr[2]),'_',FR)
    new_names=new_names%>%gsub('.seq','.fasta',.)
    for (j in 1:length(names)){
      text=readLines(names[j])
      for (i in length(text):1){
        text[i+1]=text[i]
      }
      text[1]=paste0(">",new_names[j])
      writeLines(text,file.path('../fasta',new_names[j]%>%gsub('seq','fasta',.)))
    }
  }
  packageStartupMessage(" Done!")
  
  packageStartupMessage("Vector Screening...", appendLF = FALSE)
  setwd('../fasta')
  db_vec= blast(db="../../../db/UniVec") 
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
  packageStartupMessage(" Done!")
  
  packageStartupMessage("Candidate Sequences Evaluating...", appendLF = FALSE)
  l.seq=array(dim=c(3,2,length(new_names)),
              dimnames = list(c(1:3),c('Start.pt','End.pt'),new_names))
  for (k in 1: length(new_names)){
    if(l.vec[2,1,k]%>%is.na()){
      l.seq[1,,k]=c(1,l.vec[1,1,k]-1)
      l.seq[2,,k]=c(l.vec[1,2,k]+1,data_seq$width[k])
    }else{
      if(l.vec[1,1,k]==1){
        l.seq[1,,k]=c(1,1)
      }else{
        l.seq[1,,k]=c(1,l.vec[1,1,k]-1)
      }
      l.seq[2,,k]=c(l.vec[1,2,k]+1,l.vec[2,1,k]-1)
      if(l.vec[2,2,k]!=data_seq$width[k]){
        l.seq[3,,k]=c(l.vec[2,2,k]+1,data_seq$width[k])  
      }
    }
  }

  seq.pure=seq
  for (k in 1: length(new_names)){
    if(nrow(predict(db_vec,seq[k]))!=0){
      s.seq=seq[k]%>%unlist()
      c=matrix(ncol=3,nrow=length(primer))
      if(l.seq[3,1,k]%>%is.na()){
        for (i in 1:2){
          for (j in 1:length(primer)){
            c[j,i]=s.seq[l.seq[i,1,k]:l.seq[i,2,k]]%>%as.character()%>%grepl(primer[j],.)
          }
        }
      }else{
        for (i in 1:3){
          for (j in 1:length(primer)){
            c[j,i]=s.seq[l.seq[i,1,k]:l.seq[i,2,k]]%>%as.character()%>%grepl(primer[j],.)
          }
        }
      }
      if(!(grep('TRUE',c)/length(primer))%>%isEmpty()){
        seq.pure[k]=seq[k]%>%
                    unlist()%>%
                    (function(x){x[l.seq[(grep('TRUE',c)/length(primer)) %>% ceiling()%>%mean(),1,k]:
                    l.seq[(grep('TRUE',c)/length(primer)) %>% ceiling()%>%mean(),2,k]]})%>%
                    as("DNAStringSet")
      }
    }
  }
  for (i in 1:length(seq.pure)){
    write.fasta(seq.pure[i]%>%as.DNAbin(), names(seq.pure[i])%>%gsub('.seq','',.), names(seq.pure[i])%>%gsub('.seq','.fasta',.))
  }
  packageStartupMessage(" Done!")
  
  packageStartupMessage("Sequences Alignment...", appendLF = FALSE)
  data_seq=cbind(seq@ranges%>%data.frame()%>%(function(x){cbind(x$names, x$width)}),
                 seq.pure@ranges%>%data.frame()%>%(function(x){x$width})) %>% data.frame()
  colnames(data_seq)=c('names','original','selected')
  data_seq$original=data_seq$original %>% as.character()%>% as.numeric()
  data_seq$selected=data_seq$selected %>% as.character()%>% as.numeric()
  data_seq=data.frame(data_seq,success=(data_seq$original!=data_seq$selected))
  
  seq_names=data_seq[data_seq$success==TRUE,]$names %>%
            (function(x){x[grep('_F',x)]}) %>%
            gsub('_F.fasta','',.)
  table_length=data.table()
  for (sn in 1:length(seq_names)){
    setwd(file.path(desdir,'seq',folder,'fasta'))
    set=seq.pure[grep(paste0(seq_names[sn],'_'),new_names)]
    if(length(set)==2){
      set=DNAStringSet(c(set[names(set)%>%grep('_F',.)],set[names(set)%>%grep('_R',.)]))
      aln.r=system(paste('blastn','-subject',set[2]%>%names,'-query',set[1]%>%names),intern=TRUE)
      aln.r2=aln.r[(aln.r %>% grep('Score',.)%>%min()):
                     ((aln.r %>% grep('Lambda',.)%>%min())-4)]
      
      aln.in=aln.r2 %>% grep('Score',.)
      aln=matrix(nrow=length(aln.in), ncol=2)
      for (i in 1:length(aln.in)){
        aln[i,1]=aln.in[i]
        if (i==length(aln.in)){
          aln[i,2]=length(aln.r2)
        }else{
          aln[i,2]=aln.in[i+1]-3 
        }
      }
      aln.len=aln.r2%>% grep('Identities',.)
      aln.bigind= aln.r2[aln.len] %>% 
                  strsplit(' ')%>% unlist() %>% 
                  (function(x){x[seq(4,length(x),9)]}) %>%
                  strsplit('/')%>% unlist() %>%
                  (function(x){x[seq(2,length(x),2)]}) %>%
                  (function(x){x==max(x)})
      aln.truind=aln[aln.bigind,]
      aln.t=aln.r2[aln.truind[1]:aln.truind[2]]
      aln.q=aln.t%>% grep('Query',.,value=TRUE)%>%gsub('[A-Z][a-z]*','',.)%>% 
            strsplit(' ')%>%unlist()%>%as.numeric()%>%na.omit()
      aln.s=aln.t%>% grep('Sbjct',.,value=TRUE)%>%gsub('[A-Z][a-z]*','',.)%>% 
            strsplit(' ')%>%unlist()%>%as.numeric()%>%na.omit()
      
      if(mean(aln.q)>mean(aln.s)){
        seq.aln=paste0(set[1]%>%as.character()%>% substr(1,mean(c(max(aln.q),min(aln.q)))),
                       set[2]%>%as.character()%>% substr(mean(c(max(aln.s),min(aln.s)))+1,nchar(.))) %>% 
                DNAStringSet()
      }else{
        seq.aln=paste0(set[2]%>%as.character()%>% substr(1,mean(c(max(aln.s),min(aln.s)))),
                       set[1]%>%as.character()%>% substr(mean(c(max(aln.q),min(aln.q)))+1,nchar(.))) %>% 
                DNAStringSet()
      }
      setwd(file.path(desdir,'seq',folder,'fasta.aln'))
      seq.aln%>% as.DNAbin()%>%write.fasta(names=seq_names[sn],file.out=paste(seq_names[sn],'.fasta',sep=''))
      table_length=rbind(table_length,data.table(seq_names[sn],seq.aln@ranges@width))
    }
  } 
  seq.aln= readDNAStringSet(list.files())
  packageStartupMessage(" Done!")
######################################

for (sss in 1:length(seq.aln)){
  text=seq.aln[sss]%>% as.character()%>%strsplit('')%>%unlist()
  l.pri=matrix(nrow=length(primer), ncol=2)
  for (i in 1:length(primer)){
    l.pri[i,1]=text%>%paste(collapse='')%>%gregexpr(pattern =primer[i],.)%>%as.numeric
    l.pri[i,2]=nchar(primer[i])
  }
  data_pri=data.frame(pri.pos=l.pri[,1], pri.len=l.pri[,2], NO.pri=NA, head_tail=NA, dup=NA, pair=NA)
  data_pri=data_pri[data_pri$pri.pos>0,]
  data_pri$head_tail=data_pri %>% rownames %>% as.numeric %>% (function(x){x%%2==1})
  data_pri$NO.pri=data_pri %>% rownames %>% as.numeric %>% (function(x){(x/2)%>%ceiling})
  data_pri$dup=data_pri$NO.pri %>% duplicated()
  for (i in 1: nrow(data_pri)){
    boolean=data_pri[(data_pri$NO.pri== data_pri$NO.pri[i]),]$dup
    data_pri[i,]$pair=boolean[1]|boolean[2]
  }

  data_pri
  data.frame(pri.pos=c)
    
  if(data_pri[(data_pri$pri.pos)==(data_pri$pri.pos %>%min),]$pair){   # if head_primer is paired
    if(data_pri[(data_pri$pri.pos)==(data_pri$pri.pos %>%max),]$pair){ # if tail_primer is paired given head_primer is paired
      seq.aln[sss]=text[(data_pri$pri.pos%>%min):                      # 1. both tail_primer and head_primer are paired
                       ((data_pri$pri.pos%>%max)+(data_pri$pri.len[data_pri[,1]==(data_pri[,1]%>%max)])-1)]  %>% 
                   paste(collapse='')       
      print(1)
    }else if(data_pri[data_pri$pri.pos==(data_pri$pri.pos%>%max),]$head_tail){
      seq.aln[sss]=text[(data_pri$pri.pos%>%min):                      # 2. only head_primer is paired, with head of tail_primer
                       ((data_pri$pri.pos%>%max)-1)] %>% 
                   paste(collapse='') %>% paste0(.,primerraw[data_pri$NO.pri[data_pri[,1]==(data_pri[,1]%>%max)]],collapse='')
      print(2)
    }else{
      seq.aln[sss]=text[(data_pri$pri.pos%>%min):                      # 3. only head_primer is paired, with tail of tail_primer
                       ((data_pri$pri.pos%>%max)+(data_pri$pri.len[data_pri[,1]==(data_pri[,1]%>%max)])-1)]  %>% 
                   paste(collapse='')
      print(3)
    }
  }else if(data_pri[(data_pri$pri.pos)==(data_pri$pri.pos %>%max),]$pair){ # if tail_primer is paired given head_primer isn't paired
    if(data_pri[data_pri$pri.pos==(data_pri$pri.pos%>%min),]$head_tail){
      seq.aln[sss]=text[(data_pri$pri.pos%>%min):                      # 4. only tail_primer is paired, with head of head_primer
                       ((data_pri$pri.pos%>%max)+(data_pri$pri.len[data_pri[,1]==(data_pri[,1]%>%max)])-1)]  %>% 
                   paste(collapse='')
      print(4)
    }else{                                                             # 5. only tail_primer is paired, with tail of head_primer
      seq.aln[sss]=text[((data_pri$pri.pos%>%min)+(data_pri$pri.len[data_pri[,1]==(data_pri[,1]%>%min)])): 
                        ((data_pri$pri.pos%>%max)+(data_pri$pri.len[data_pri[,1]==(data_pri[,1]%>%max)])-1)] %>%
                   paste(collapse='') %>% paste0(primerraw[data_pri$NO.pri[data_pri[,1]==(data_pri[,1]%>%min)]],., collapse='')
      print(5)
    }
  }else if(data_pri$head_tail[data_pri$pri.pos==(data_pri$pri.pos %>% min)]){
    if(data_pri$head_tail[data_pri$pri.pos==(data_pri$pri.pos %>% max)]){
      seq.aln[sss]=text[(data_pri$pri.pos%>%min):                       # 6. exist head of head_primer and tail_primer
                       ((data_pri$pri.pos%>%max)-1)] %>% 
                   paste(collapse='') %>% paste0(.,primerraw[data_pri$NO.pri[data_pri[,1]==(data_pri[,1]%>%max)]], collapse='')
    }else{
      seq.aln[sss]=text[(data_pri$pri.pos%>%min):                       # 7. exist head of head_primer and tail of tail_primer
                       ((data_pri$pri.pos%>%max)+(data_pri$pri.len[data_pri[,1]==(data_pri[,1]%>%max)])-1)]  %>% 
                   paste(collapse='')
    }
  }else if(data_pri$head_tail[data_pri$pri.pos==(data_pri$pri.pos %>% max)]){  # 8. exist tail of head_primer and head of tail_primer
    seq.aln[sss]=text[((data_pri$pri.pos%>%min)+(data_pri$pri.len[data_pri[,1]==(data_pri[,1]%>%min)])):
                      ((data_pri$pri.pos%>%max)-1)] %>%
                 paste(collapse='') %>% paste0 (primerraw[data_pri$NO.pri[data_pri[,1]==(data_pri[,1]%>%min)]],
                                                ., 
                                                primerraw[data_pri$NO.pri[data_pri[,1]==(data_pri[,1]%>%max)]], collapse='')
  }else{                                                                # 9. exist tail of head_primer and tail_primer
    seq.aln[sss]=text[((data_pri$pri.pos%>%min)+(data_pri$pri.len[data_pri[,1]==(data_pri[,1]%>%min)])): 
                      ((data_pri$pri.pos%>%max)+(data_pri$pri.len[data_pri[,1]==(data_pri[,1]%>%max)])-1)] %>%
                 paste(collapse='') %>% paste0(primerraw[data_pri$NO.pri[data_pri[,1]==(data_pri[,1]%>%min)]],.,collapse='')
  }
}

######################################

  if (nt_search){
    packageStartupMessage(paste0("Fetching sequence information...\nIt may take at least ",
                                 length(seq.aln)*3,
                                 " mins to complete."), appendLF = FALSE)
    db_nt <- blast(db="../../../db/nt") 
    report=data.frame()
    options(download.file.method = "wininet")
    for (i in 1: length(seq.aln)){
      if (nchar(seq.aln[i])>100){
        cl <- predict(db_nt, seq.aln[i]) 
        if(!cl[1]%>%isEmpty()){
          if (nrow(cl)<10){
            x=cl[order(-cl$Bits),] %>%
            (function(x){x[1:nrow(cl),2]})
          }else{
            x=cl[order(-cl$Bits),] %>%
            (function(x){x[1:10,2]})
          }
          x2=x %>%
             as.character() %>%
             strsplit('\\|') %>%
             unlist()
          y=x2[seq(4,length(x2),4)]%>%
            genbank2uid()%>%
            ncbi_get_taxon_summary()
          z=x2[seq(2,length(x2),4)]%>%
            entrez_summary(db='Nucleotide',id=.)%>%
            extract_from_esummary('title')
          if (nrow(cl)<10){
            report=rbind(report,
                         data.frame(Seq.Names=cl[1:nrow(cl),c(1)],
                                    description=z,y[,c(2,3)],
                                    cl[1:nrow(cl),c(3,4,5,6,7,8,9,10)]))
          }else{
            report=rbind(report,
                         data.frame(Seq.Names=cl[1:10,c(1)],
                                    description=z,y[,c(2,3)],
                                    cl[1:10,c(3,4,5,6,7,8,9,10)]))
          }
        }
      }
    }
    packageStartupMessage(" Done!")
    write.csv(report,'../summary_seq.csv',row.names = FALSE)
  }
  
  packageStartupMessage("Exporting Summary Information...", appendLF = FALSE)
  sum_names= gsub('_F.fasta','',data_seq$names) %>%
             grep('_R',.,value=TRUE,invert=TRUE)
  data_summary=matrix(nrow=length(sum_names),ncol=7)
  colnames(data_summary)=c('seq.names','F_length_raw','R_length_raw','F_length_vs','R_length_vs','Suc_or_Fal','aln_length')
  
  for (su in 1:length(sum_names)){
    index=data_seq$names %>% grep(paste(sum_names[su],'_',sep=''),.)    
    index_f= index[data_seq$names[index] %>% grep('F',.)]
    index_r= index[data_seq$names[index] %>% grep('R',.)]
    data_summary[su,1]=sum_names[su]
    if(index_f%>%isEmpty()){
      data_summary[su,2]=NA
      data_summary[su,3]=data_seq[index_r,2]
      data_summary[su,4]=NA
      data_summary[su,5]=data_seq[index_r,3]
      data_summary[su,6]='Failure'
      data_summary[su,7]=NA
    }else if(index_r%>%isEmpty()){
      data_summary[su,2]=data_seq[index_f,2]
      data_summary[su,3]=NA
      data_summary[su,4]=data_seq[index_f,3]
      data_summary[su,5]=NA
      data_summary[su,6]='Failure'
      data_summary[su,7]=NA
    }else{
      data_summary[su,2]=data_seq[index_f,2]
      data_summary[su,3]=data_seq[index_r,2]
      data_summary[su,4]=data_seq[index_f,3]
      data_summary[su,5]=data_seq[index_r,3]
      data_summary[su,6]=if(data_seq[index_f,4]&data_seq[index_r,4]){'Success'}else{'Failure'}
      if(!table_length[,V1]%>%grep(sum_names[su],.)%>%isEmpty()){
        data_summary[su,7]=table_length[table_length[,V1]%>%grep(paste0(sum_names[su],'$'),.),V2]
      }else{
        data_summary[su,7]=NA
      }
    }
  }
  data_summary=data_summary%>% as.data.frame()
  write.csv(data_summary,'../summary_aln.csv')
  packageStartupMessage(paste0('Done! Program End! \n\n\nFiles Location: ',
                               file.path(desdir,'seq',folder)))
}

DownloadFTP=function(source, username, password, des_folder){
  filenames= getURL(source,userpwd=paste0(username,':',password),
                    verbose=TRUE,ftp.use.epsv=TRUE, dirlistonly = TRUE) %>%
             strsplit("[\\\\]|[^[:print:]]",fixed = FALSE) %>%
             unlist() %>% (function(x){x[grep('seq',x)]})
  filepath=sprintf(paste0('ftp://',
                         paste0(username,':',password),'@',
                         (source%>%gsub('ftp://','',.)),'%s'),filenames)
  if (!file.exists(file.path(desdir,'Download',des_folder))){
    dir.create(file.path(desdir,'Download',des_folder),recursive = TRUE)
  }
  for (i in 1:(length(filenames))){
    download.file(filepath[i],
                  file.path(desdir,'Download',des_folder,filenames[i]))
  }
}

FetchingSeq=function(folder){
  setwd(filePath(desdir,'seq',folder,'fasta.aln'))
  seq.aln= readDNAStringSet(list.files())
  db_nt <- blast(db="../../../db/nt") 
  report=data.frame()
  options(download.file.method = "wininet")
  packageStartupMessage(paste0("Fetching sequence information...\nIt may take at least ",
                               length(seq.aln)*3,
                               " mins to complete."), appendLF = FALSE)
  for (i in 1: length(seq.aln)){
    if (nchar(seq.aln[i])>100){
      cl <- predict(db_nt, seq.aln[i]) 
      if(!cl[1]%>%isEmpty()){
        if (nrow(cl)<10){
          x=cl[order(-cl$Bits),] %>%
            (function(x){x[1:nrow(cl),2]})
        }else{
          x=cl[order(-cl$Bits),] %>%
            (function(x){x[1:10,2]})
        }
        x2=x %>%
          as.character() %>%
          strsplit('\\|') %>%
          unlist()
        y=x2[seq(4,length(x2),4)]%>%
          genbank2uid()%>%
          ncbi_get_taxon_summary()
        z=x2[seq(2,length(x2),4)]%>%
          entrez_summary(db='Nucleotide',id=.)%>%
          extract_from_esummary('title')
        if (nrow(cl)<10){
          report=rbind(report,
                       data.frame(Seq.Names=cl[1:nrow(cl),c(1)],
                                  description=z,y[,c(2,3)],
                                  cl[1:nrow(cl),c(3,4,5,6,7,8,9,10)]))
        }else{
          report=rbind(report,
                       data.frame(Seq.Names=cl[1:10,c(1)],
                                  description=z,y[,c(2,3)],
                                  cl[1:10,c(3,4,5,6,7,8,9,10)]))
        }
      }
    }
  }
}
