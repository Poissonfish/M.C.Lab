desdir="/home/mclab/R/git/M.C.Lab" 
setwd(desdir)
source('./MCLab_Function.R') 

#######################
folder='ch11-11'
primer_f="GGTGCTCAAGGCGGGACATTCGTT"
primer_r="TCGGGATTGGCCACAGCGTTGAC"
name_primer_f='SP6'
name_primer_r='T7'
#######################
source="ftp://140.109.56.5/104DATA/0504/"
username="rm208"
password="167cm"
#######################
local=FALSE
source_path='/home/mclab/R/git/M.C.Lab/11-4'
#######################
##############  #######
#############   #######
############   #######
MCLab(primer_f, primer_r, name_primer_f, name_primer_r, source, username, password, desdir,folder,source_path, local) ###
#############   #######
##############  #######
#######################





#######################
des_folder='0519'
source="ftp://140.109.56.5/104DATA/0509/"
username="rm208"
password="167cm"
#######################
##############  #######
#############   #######
############    #######
DownloadFTP(source, username, password, folder)
#############   #######
##############  #######
#######################

#11-6
primer_f="ATGGCGAGAACTAARCANAT"
primer_r="AGYTGNATRTCCTTYTGCAT"
#11-4
primer_f="ATGGCGAGAACTAARCANAT"
primer_r="TYAYGAAAAATGTCKWGMRCCA"
#11-11


