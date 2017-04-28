
library(dplyr)
setwd('/home/ark06eu/evolgen/data/1135_data_Rdata/')

source('chr_sort_func.r')
C_1 <- Chr_1()


colnames(C_1) <- as.numeric(matrix(unlist(strsplit(colnames(C_1), split='1- ')),
                                 ncol = 2, byrow = T)[,2])
C_1 <- rbind(as.numeric(colnames(C_1)),C_1)
Time <- Sys.time()
ws = 500
ma <- max(C_1[1,])

my_len <- length(seq(1,(ma-ws), by = ws))
Haplos <- vector(length=my_len, mode = 'numeric')
SNPs <- vector(length=my_len, mode = 'numeric')
Pos_start <- vector(length=my_len, mode = 'numeric' )
Pos_end <- vector(length=my_len, mode = 'numeric')
k = 1
is.numeric(Pos_start)

for(i in seq(1,(ma-ws), by = ws)){
  j = i  + ws
  mywin <- which(C_1[1,]>=i & C_1[1,]<j)
  if(length(mywin)>1){ 
    mma <- max(mywin)
    mmi <- min(mywin)
    A <- Chr_1[2:dim(C_1)[1],mmi:mma]
    uniq <- distinct(as.data.frame(A))
    Haplos[k] <- as.numeric(dim(uniq)[1])
    SNPs[k] <- as.numeric(dim(uniq)[2])
    Pos_start[k] <- as.numeric(colnames(A)[1])
    Pos_end[k] <-  as.numeric(colnames(A)[length(colnames(A))])
  }else if(length(mywin)==0){
    if(k==1){
      Haplos[k] <- 0
      SNPs[k]  <- 0
      Pos_start[k] <- 1
      Pos_end[k] <- 500
    }else{
      Haplos[k] <- 0
      SNPs[k]  <- 0
      Pos_start[k] <- as.numeric(Pos_end[k-1])
      Pos_end[k] <- as.numeric(Pos_start[k])+ws
    }
  }else if(length(mywin)==1){
    Haplos[k] <- 2
    SNPs[k]  <- 1
    Pos_start[k] <- as.numeric(colnames(C_1)[mywin])
    Pos_end[k] <-as.numeric(colnames(C_1)[mywin])
  }  
  k = k+1
  if(k%%100==0){
    print(paste(k, 'von', my_len),sep =' ')
  }
}


Haplos_kb <- as.data.frame(cbind(Pos_start,Pos_end,Haplos,SNPs))
colnames(Haplos_kb) <- c('Start','End','Haplos','SNPs')
Jane <- Time - Sys.time()


setwd('/home/ark06eu/evolgen/jan/Haplos_type_scripte/')


write.csv(Haplos_kb,"Chr_1_5HB.csv")



