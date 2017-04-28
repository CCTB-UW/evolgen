hallo<-Sys.time()
require(adegenet)
require(pegas)
setwd("/home2/ark06eu/evolgen/data/")
source('../jan/find_haplo.r')
Haplos_NAO_Gesamt <- c()
for (v in seq(0,99000,by=1000)){
Data <- read.vcf("../vcf1.vcf", quiet = T,from = (v+1), to = (v+1000))
#colnames(Data) <- c(paste("SNP_",seq(1:ncol(Data)),sep=""))
setwd("/home2/ark06eu/evolgen/jan/")
Data_OT <- Remove_Triall(Data)
Ref     <- Ref_All_2(Data_OT)
Data_OT_01 <- Transform_Data_Fast(Data_OT,Ref)
a <- Find_Haplos_2(100,Data_OT_01)
Haplos_NAO_Gesamt <- append(Haplos_NAO_Gesamt,a)
print(v)
}
save(Haplos_NAO_Gesamt,file='output_function.rda')

install.packages(c("rafalib","vioplot","dplyr"))

library(rafalib)
library(vioplot)
library(dplyr)
Hallo <- Sys.time()
setwd("/home2/ark06eu/evolgen/jan/Haplos_type_scripte/")
source('find_haplo.r')
setwd("/home2/ark06eu/evolgen/data/1135_data_Rdata/")
ws = 200
for (u in 1:2) {
  filename<-paste('X_1135_',u,'.rda.gz',sep='')
  load(filename)
  A <- as.matrix(X)
  if (u==1) {
    Results <- Find_Haplos_PLYR(ws,A)
    rm(X)
    print(u)
  } else{
    Mehr_Results <- Find_Haplos_PLYR(ws,X)
    Results <-rbind(Results,Mehr_Results)
    rm(X)
    print(u)
  }
}
  
  
End <- Sys.time()-Hallo

pdf("Plots_Haplos_R")
mypar(1,3)
plot(Results_2$Haplos)
hist(Results_2$Haplos)
vioplot(Results_2$Haplos)
dev.off()

End
