# #######################
# folder='20489527'
# primer_f="GGTGCTCAAGGCGGGACATTCGTT"
# primer_r="TCGGGATTGGCCACAGCGTTGAC"
# name_primer_f='SP6'
# name_primer_r='T7'
# #######################
# source="ftp://140.109.56.5/104DATA/0504/"
# username="rm208"
# password="167cm"
# #######################
# local=TRUE
# source_path='/home/mclab/R/git/M.C.Lab/seq/ch11-2048/raw'

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
  primern=c()
  for (i in 1:length(primer)){
    primern=c(primern,primer[i] %>% substr(.,1,nchar(.)-5),primer[i] %>% substr(.,6,nchar(.))) 
  }
  primer=primern
  packageStartupMessage(" Done!")
  
  setwd(file.path(desdir,'seq',folder,'raw'))
  if(!local){
    packageStartupMessage("File Downloading...", appendLF = FALSE)
    filenames= getURL(source,userpwd=paste(username,':',password,sep=''),
                      verbose=TRUE,ftp.use.epsv=TRUE, dirlistonly = TRUE) %>%
      strsplit("[\\\\]|[^[:print:]]",fixed = FALSE) %>%
      unlist() %>% (function(x){x[grep('seq',x)]})
    filepath=sprintf(paste('ftp://',
                           paste(username,':',password,sep=''),'@',
                           (source%>%gsub('ftp://','',.)),'%s',sep=''),filenames)
    
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
    gsub(paste('(',name_primer_r,')',sep=''),"_R",.)
  new_names[grep(name_primer_f,names)]= names[grep(name_primer_f,names)] %>% 
    gsub("^[0-9][0-9]\\.","",.) %>%
    gsub(paste('(',name_primer_f,')',sep=''),"_F",.)
   
  new_names_split= new_names %>% strsplit('_') %>% unlist() 
  SN= new_names_split[seq(1,length(new_names_split),2)]  
  FR= new_names_split[seq(2,length(new_names_split),2)]  
  nchr= SN %>%  nchar() %>% (function(x){c(min(x),max(x))})
  
  s_index= SN%>%nchar() == nchr[1]
  l_index= SN%>%nchar() == nchr[2]
  
  new_names[l_index]=paste(SN[l_index] %>% substr(1,nchr[2]-3),'-',
                           SN[l_index] %>% substr(nchr[2]-2,nchr[2]-2),'-',
                           SN[l_index] %>% substr(nchr[2]-1,nchr[2]),'_', FR[l_index],sep='')
  new_names[s_index]=paste(SN[s_index] %>% substr(1,nchr[1]-2),'-',
                           SN[s_index] %>% substr(nchr[1]-1,nchr[1]-1),'-',
                           SN[s_index] %>% substr(nchr[1],nchr[1]),'_', FR[s_index],sep='')
  new_names=new_names%>%gsub('.seq','.fasta',.)
  
  for (j in 1:length(names)){
    text=readLines(names[j])
    for (i in length(text):1){
      text[i+1]=text[i]
    }
    text[1]=paste(">",new_names[j],sep='')
    writeLines(text,file.path('../fasta',new_names[j]%>%gsub('seq','fasta',.)))
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
      }else{l.seq[1,,k]=c(1,l.vec[1,1,k]-1)}
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
  
  
  # for (sn in 1:length(seq_names)){
    set=DNAStringSet(seq.pure[grep(paste0(seq_names[sn],'_'),new_names)])
    aln.r=system(paste('blastn','-subject',set[1]%>%names,'-query',set[2]%>%names),intern=TRUE)
    
  #   msa=muscle(set)
  #   data_msa=msa%>% as.character()
  #   
  #   if(length(set)==2){
  #     ident=c()
  #     for (i in 1:ncol(data_msa)){
  #       ident=c(ident,(data_msa[1,i]==data_msa[2,i]))
  #     }
  #     idnum=(1:ncol(data_msa))[ident==FALSE]
  # 
  #     if(isEmpty(idnum)!=TRUE){
  #       if (length(idnum)!=1){
  #         gap=c()
  #         for(i in 2:length(idnum)){
  #           gap=c(gap,idnum[i]-idnum[i-1])
  #         }
  #         if (mean(gap)!=1){
  #           index_for=1:idnum[(1:length(gap))[gap==max(gap)][(((1:length(gap))[gap==max(gap)] -(length(idnum)/2)) %>% abs() %>% order()) ==1]]
  #           index_rev=idnum[(1:length(gap))[gap==max(gap)][(((1:length(gap))[gap==max(gap)] -(length(idnum)/2)) %>% abs() %>% order()) ==1]+1]:ncol(data_msa)
  #           index_con=(max(index_for)+1):(min(index_rev)-1)
  #           seq_aln = c(data_msa[row.names(data_msa)%>%grep('F',.),index_for],
  #                       data_msa[1,index_con],
  #                       data_msa[row.names(data_msa)%>%grep('R',.),index_rev])
  #           if(isEmpty(grep('-',seq_aln))!=TRUE){seq_aln=seq_aln[-grep('-',seq_aln)]}
  #         }else{
  #           if(mean(idnum)<mean(ncol(data_msa))){
  #             if(data_msa[grep('R',rownames(data_msa)),1]=='-'){
  #               seq_aln = c(data_msa[grep('F',rownames(data_msa)),1:((ncol(data_msa)/2)%>%ceiling())],
  #                           data_msa[grep('R',rownames(data_msa)),(((ncol(data_msa)/2)%>%ceiling())+1):ncol(data_msa)])
  #             }else{
  #               data_msa=data_msa[,-(data_msa[grep('F',rownames(data_msa)),]%>% grep('-',.))]
  #               seq_aln = c(data_msa[grep('F',rownames(data_msa)),1:((ncol(data_msa)/2)%>%ceiling())],
  #                           data_msa[grep('R',rownames(data_msa)),(((ncol(data_msa)/2)%>%ceiling())+1):ncol(data_msa)])
  #             }
  #           }else{
  #             if(data_msa[grep('F',rownames(data_msa)),1]=='-'){
  #               seq_aln = c(data_msa[grep('F',rownames(data_msa)),1:((ncol(data_msa)/2)%>%ceiling())],
  #                           data_msa[grep('R',rownames(data_msa)),(((ncol(data_msa)/2)%>%ceiling())+1):ncol(data_msa)])
  #             }else{
  #               data_msa=data_msa[,-(data_msa[grep('R',rownames(data_msa)),]%>% grep('-',.))]
  #               seq_aln = c(data_msa[grep('F',rownames(data_msa)),1:((ncol(data_msa)/2)%>%ceiling())],
  #                           data_msa[grep('R',rownames(data_msa)),(((ncol(data_msa)/2)%>%ceiling())+1):ncol(data_msa)])
  #             }
  #           }
  #         }
  #       }else{
  #         seq_aln = c(data_msa[row.names(data_msa)%>%grep('F',.),1:idnum],
  #                     data_msa[row.names(data_msa)%>%grep('R',.),(idnum+1):(ncol(data_msa))]) 
  #         if(isEmpty(grep('-',seq_aln))!=TRUE){seq_aln=seq_aln[-grep('-',seq_aln)]}
  #       } 
  #     }else{
  #       seq_aln = data_msa[1,]  
  #     }
  #     setwd(file.path(desdir,'seq',folder,'fasta.aln'))
  #     seq_aln%>%as.DNAbin()%>%write.fasta(names=seq_names[sn],file.out=paste(seq_names[sn],'.fasta',sep=''))
  #     table_length=rbind(table_length,data.table(seq_names[sn],length(seq_aln)))
  #   }
  # }
  # seq.aln= readDNAStringSet(list.files())
  # packageStartupMessage(" Done!")

  if (nt_search){
    packageStartupMessage(paste0("Fetching sequence information...\nIt may take at least ",length(seq.aln)*3," mins to complete."), appendLF = FALSE)
    db_nt <- blast(db="../../../db/nt") 
    report=data.frame()
    options(download.file.method = "wininet")
    for (i in 1: length(seq.aln)){
      if (nchar(seq.aln[i])>100){
        cl <- predict(db_nt, seq.aln[i]) 
        if(!cl[1]%>%isEmpty()){
          x=cl[order(-cl$Bits),] %>% 
            (function(x){x[1:10,2]}) %>% 
            as.character() %>%
            strsplit('\\|') %>%
            unlist()
          y=x[seq(4,length(x),4)]%>%
            genbank2uid()%>%
            ncbi_get_taxon_summary()
          z=x[seq(2,length(x),4)]%>%
            entrez_summary(db='Nucleotide',id=.)%>%
            extract_from_esummary('title')
          report=rbind(report,
                       data.frame(Seq.Names=cl[1:10,c(1)],
                                  description=z,y[,c(2,3)],
                                  cl[1:10,c(3,4,5,6,7,8,9,10)]))
        }
      }
    }
    packageStartupMessage(" Done!")
    write.csv(report,'../summary_seq.csv',row.names = FALSE)
  }
  
  packageStartupMessage("Exporting Summary Information...", appendLF = FALSE)
  sum_names= gsub('_F.seq','',data_seq$names) %>%
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
      }else{data_summary[su,7]=NA}
    }
  }
  data_summary=data_summary%>% as.data.frame()
  write.csv(data_summary,'../summary_aln.csv')
  packageStartupMessage(paste0('Done! Program End! \n\n\nFiles Location: ',file.path(desdir,'seq',folder)))
}

DownloadFTP=function(source, username, password, des_folder){
  filenames= getURL(source,userpwd=paste(username,':',password,sep=''),
                    verbose=TRUE,ftp.use.epsv=TRUE, dirlistonly = TRUE) %>%
    strsplit("[\\\\]|[^[:print:]]",fixed = FALSE) %>%
    unlist() %>% (function(x){x[grep('seq',x)]})
  filepath=sprintf(paste('ftp://',
                         paste(username,':',password,sep=''),'@',
                         (source%>%gsub('ftp://','',.)),'%s',sep=''),filenames)
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
  packageStartupMessage(paste0("Fetching sequence information...\nIt may take at least ",length(seq.aln)*3," mins to complete."), appendLF = FALSE)
  for (i in 1: length(seq.aln)){
    if (nchar(seq.aln[i])>100){
      cl <- predict(db_nt, seq.aln[i]) 
      x=cl[order(-cl$Bits),] %>% 
        (function(x){x[1:10,2]}) %>% 
        as.character() %>%
        strsplit('\\|') %>%
        unlist()
      y=x[seq(4,length(x),4)]%>%
        genbank2uid()%>%
        ncbi_get_taxon_summary()
      z=x[seq(2,length(x),4)]%>%
        entrez_summary(db='Nucleotide',id=.)%>%
        extract_from_esummary('title')
      report=rbind(report,
                   data.frame(Seq.Names=cl[1:10,c(1)],
                              description=z,y[,c(2,3)],
                              cl[1:10,c(3,4,5,6,7,8,9,10)]))
    }
  }
  packageStartupMessage(" Done!")
  write.csv(report,'../summary_seq.csv',row.names = FALSE)
}