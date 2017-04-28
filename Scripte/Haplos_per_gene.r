### R script fuer Haplotypes per gene und wo wir schon mal dabein sind SNPs MAC ..

setwd('/home/ark06eu/evolgen/data/1135_data_Rdata/')
library(dplyr)

Ara11 <- read.table('Araport11_protein_coding.201606.bed')
Ara11_Chr1 <- subset(Ara11, Ara11$V1=="Chr1")
Ara11_1 <- Ara11_Chr1[,1:5]
rm(Ara11_Chr1)
Ara11_1$V6 <- NA
Ara11_1$V7 <- NA
Ara11_1$V8 <- NA
colnames(Ara11_1) <- c('Chr','Start','End','Gene','kp', 'Haplos','SNPs','length')
u = 6
Time <- Sys.time()
for(u in 1:6){
filename <- paste('X_1135_',u,'.rda.gz',sep = '')
load(filename)
if(u==6){
  X <- X[,1:97826]
}
my_coln <- as.vector(as.integer(strsplit(colnames(X),'- ') %>% unlist))
my_rm <- which(my_coln == 1)
my_coln <- my_coln[-my_rm]

colnames(X) <- my_coln
X_ <- (rbind(my_coln,X))

my_max <- max(X_[1,])
my_min <- min(X_[1,])
my_start <- which.min(abs(Ara11_1$Start - my_min ))
my_end <- which.min(abs(Ara11_1$End - my_max )) 

for(i in my_start:my_end){
cols <- as.character(subset(X_[1,], (X_)[1,] <= Ara11_1$End[i] & X_[1,] >= Ara11_1$Start[i]))
window <- (which(colnames(X_)%in%cols))
if(!length(window)==0){
dw <- (distinct(as.data.frame(X_[2:dim(X_)[1],min(window):max(window)])))
Ara11_1$Haplos[i]<- dim(dw)[1]
Ara11_1$SNPs[i] <- dim(dw)[2]
Ara11_1$length[i] <- (Ara11_1$End[i] - Ara11_1$Start[i])
} else {
  Ara11_1$Haplos[i]<- 'NOSNP'
  Ara11_1$SNPs[i]  <- 'NOSNP'
  Ara11_1$length[i]<- (Ara11_1$End[i] - Ara11_1$Start[i])
}
if(i%%250 == 0){
  print(i)
}
}
print(u)
print(Time - Sys.time())
}

my_vec <- which(Ara11_1$Haplos == 'NOSNP')
Ara11_1$Haplos
write('Ara11_1_HPG.rda',Ara11_1)
write.csv(Ara11_1, 'Ara11_1_HPG.csv')
dir()
