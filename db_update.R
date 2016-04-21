library(RCurl)
library(magrittr)

setwd("/home/mclab/R/git/M.C.Lab/db/")
db_names=getURL("ftp://ftp.ncbi.nlm.nih.gov/blast/db/",verbose=TRUE,ftp.use.epsv=TRUE, dirlistonly = TRUE) %>%
  strsplit("[\\\\]|[^[:print:]]",fixed = FALSE) %>%
  unlist() %>%
  (function(x){x[grep('^nt\\.',x)]})

#UPADATING DABABASE FROM NCBI FTP
  filepath=sprintf('ftp://ftp.ncbi.nlm.nih.gov/blast/db/%s',db_names)
  for (i in 1:(length(db_names)/2)){
    download.file(filepath[i],
                  paste(getwd(),'/',db_names[i],sep=''),
                  mode='wb')
  }
  Sys.sleep(60)
  for (i in ((length(db_names)/2)+1):length(db_names)){
    download.file(filepath[i],
                  paste(getwd(),'/',db_names[i],sep=''),
                  mode='wb')
  }

#Extract files from compressed items
  for (i in grep("^(?!.*md5)",db_names, perl=TRUE)){
    untar(db_names[i])
  }
 