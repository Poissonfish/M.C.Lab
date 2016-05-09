primer_f="GGTGCTCAAGGCGGGACATTCGTT"
primer_r="TCGGGATTGGCCACAGCGTTGAC"
source="ftp://140.109.56.5/104DATA/0504/"
username="rm208"
password="167cm"
desdir="/home/mclab/R/git/M.C.Lab/seq/"


VA=function(primer_f, primer_r, source, username, password, desdir){

  primer_f2=primer_f%>%DNAString()%>%reverseComplement()%>%as.character()
  primer_r2=primer_r%>%DNAString()%>%reverseComplement()%>%as.character()
  primer=c(primer_f,primer_f2,primer_r,primer_r2)
  
  filenames= getURL(source,userpwd=paste(username,':',password,sep=''),
                    verbose=TRUE,ftp.use.epsv=TRUE, dirlistonly = TRUE) %>%
             strsplit("[\\\\]|[^[:print:]]",fixed = FALSE) %>%
             unlist() %>% (function(x){x[grep('seq',x)]})
  filepath=sprintf(paste('ftp://',
                         paste(username,':',password,sep=''),'@',
                         (source%>%gsub('ftp://','',.)),'%s',sep=''),filenames)

  for (i in 1:(length(filenames))){
    download.file(filepath[i],
                  paste('./',filenames[i],sep=''))
                  }
  names=list.files()
  new_names=list.files()%>% gsub("^[0-9][0-9]\\.","",.) %>% gsub("(T7)","_R",.) %>% gsub("(SP6)","_F",.)

  for (j in 1:length(names)){
    text=readLines(names[j])
    for (i in length(text):1){
      text[i+1]=text[i]
    }
    text[1]=paste(">",new_names[j],sep='')
    writeLines(text,new_names[j])
  }
  
  db_vec= blast(db="../db/UniVec") 
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
  
  data_seq=cbind(seq@ranges%>%data.frame()%>%(function(x){x[,c(4,3)]}),
                 seq.pure@ranges%>%data.frame()%>%(function(x){x[,c(3)]}))
  colnames(data_seq)=c('names','original','selected')
  data_seq=data.frame(data_seq,success=(data_seq$original!=data_seq$selected))
  
  seq_names=data_seq[data_seq$success==TRUE,]$names %>%
    (function(x){x[grep('_F',x)]}) %>%
    gsub('_F.seq','',.)
  
  for (sn in 1:length(seq_names)){
    set=DNAStringSet(seq.pure[grep(seq_names[sn],new_names)])
    msa=muscle(set)
    data_msa=msa%>% as.character()
    
    if(length(set)==2){
      ident=c()
      for (i in 1:ncol(data_msa)){
        ident=c(ident,(data_msa[1,i]==data_msa[2,i]))
      }
      idnum=(1:ncol(data_msa))[ident==FALSE]
      
      if(isEmpty(idnum)!=TRUE){
        if (length(idnum)!=1){
          gap=c()
          for(i in 2:length(idnum)){
            gap=c(gap,idnum[i]-idnum[i-1])
          }
          if (mean(gap)!=1){
            index_for=1:idnum[(1:length(gap))[gap==max(gap)]]
            index_rev=idnum[(1:length(gap))[gap==max(gap)]+1]:ncol(data_msa)
            index_con=(max(index_for)+1):(min(index_rev)-1)
            seq_aln = c(data_msa[row.names(data_msa)%>%grep('F',.),index_for],
                        data_msa[1,index_con],
                        data_msa[row.names(data_msa)%>%grep('R',.),index_rev])
            if(isEmpty(grep('-',seq_aln))!=TRUE){seq_aln=seq_aln[-grep('-',seq_aln)]}
          }else{
            if(mean(idnum)<mean(ncol(data_msa))){
              if(data_msa[grep('R',rownames(data_msa)),1]=='-'){
                seq_aln = c(data_msa[grep('F',rownames(data_msa)),1:((ncol(data_msa)/2)%>%ceiling())],
                            data_msa[grep('R',rownames(data_msa)),(((ncol(data_msa)/2)%>%ceiling())+1):ncol(data_msa)])
              }else{
                data_msa=data_msa[,-(data_msa[grep('F',rownames(data_msa)),]%>% grep('-',.))]
                seq_aln = c(data_msa[grep('F',rownames(data_msa)),1:((ncol(data_msa)/2)%>%ceiling())],
                            data_msa[grep('R',rownames(data_msa)),(((ncol(data_msa)/2)%>%ceiling())+1):ncol(data_msa)])
              }
            }else{
              if(data_msa[grep('F',rownames(data_msa)),1]=='-'){
                seq_aln = c(data_msa[grep('F',rownames(data_msa)),1:((ncol(data_msa)/2)%>%ceiling())],
                            data_msa[grep('R',rownames(data_msa)),(((ncol(data_msa)/2)%>%ceiling())+1):ncol(data_msa)])
              }else{
                data_msa=data_msa[,-(data_msa[grep('R',rownames(data_msa)),]%>% grep('-',.))]
                seq_aln = c(data_msa[grep('F',rownames(data_msa)),1:((ncol(data_msa)/2)%>%ceiling())],
                            data_msa[grep('R',rownames(data_msa)),(((ncol(data_msa)/2)%>%ceiling())+1):ncol(data_msa)])
              }
            }
          }
        }else{
          seq_aln = c(data_msa[row.names(data_msa)%>%grep('F',.),1:idnum],
                      data_msa[row.names(data_msa)%>%grep('R',.),(idnum+1):(ncol(data_msa))]) 
          if(isEmpty(grep('-',seq_aln))!=TRUE){seq_aln=seq_aln[-grep('-',seq_aln)]}
        } 
      }else{
        seq_aln = data_msa[1,]  
      }
      seq_aln%>%as.DNAbin()%>%write.fasta(names=seq_names[sn],file.out=paste(seq_names[sn],'.fasta',sep=''))
    }
  }
}