library(adegenet)
library(pegas)
getwd()
setwd("/home2/ark06eu/evolgen/data/")
Data <- read.vcf("f0b654b16c7000f8b6f1f216384def4e_snpeff.vcf", quiet = T)
colnames(Data) <- c(paste("SNP_",seq(1:ncol(Data)),sep=""))

Ref_Allel(Data)
Transform_Data(Data,Ref)
Find_Haplos_2(50,B_)
A <- B_
Ref_2 <- Ref_Allel(Data_2)
B_[45,8]



Ref_Allel <- function(Data){
Ref = vector(length = 0)
for (i in 1:ncol(Data)){
  a <- which(Data[,i]==levels(Data[,i])[1])
  b <- which(Data[,i]==levels(Data[,i])[2])
  c <- which(Data[,i]==levels(Data[,i])[3])
  if (max(length(a),length(b), length(c))==length(a)){
    Ref <- append(Ref,levels(Data[,i])[1])
  } else if (max(length(a),length(b), length(c))==length(b)){
    Ref <- append(Ref,levels(Data[,i])[2])
  } else {
    Ref <- append(Ref,levels(Data[,i])[3])
  }  
  if (Ref[i]=="./."){
    if (levels(Data[,i])[1]!= "./."){
      Ref[i] <- levels(Data[,i])[1] 
     } else if (levels(Data[,i])[1]== "./.") {
       Ref[i] <- levels(Data[,i])[2]
     }
  }
}
print(Ref)
}
B__ <- Transform_Data(Data_2,Ref_2)

Transform_Data <- function(Data,Ref){
B <- matrix(ncol=ncol(Data),nrow=nrow(Data))
for (i in 1:nrow(B)){
  for (j in 1:ncol(B)){
    if ((Data[i,j])==Ref[j]){
      B[i,j] <- 0
    } else if (Data[i,j]=="./."){
      B[i,j] <- 0.5
    } else {
      B[i,j] <- 1
    }
  } 
}
}

Find_Haplos <-function(ws,A){
  ws = ws -1
  Haplos <- vector(mode="numeric", length=0)
  Position <- vector(mode="numeric", length=0)
  for(i in 1:dim(A)[2]){
    A_ <- A[1:dim(A)[1], (ws*(i/2)):((ws*(i/2))+ws)]
    A__ <- c(which(!duplicated(A_)))
    Number_Haplos <- length(A__)
    #print(Number_Haplos)#
    #print(i)#
    Haplos <- append(Haplos, Number_Haplos)
    Position<- append(Position,i)
    if (i == (round(((dim(A)[2])/(ws/2)),0)-2)){
      #if ( i == 5){#
      break
    }
  }
  plot(Position,Haplos)
  print(Haplos)
  #shapiro.test(Haplos)#
}

Ref_Allel(Data)
Transform_Data(Data,Ref)
Find_Haplos(10,B)

NAS <- vector(length=0)
for (i in 1:nrow(B)){
  NAS[i] <- sum(is.na(B[,i]))
  
}
NAS
write
library(MASS)
getwd()
setwd("/home2/jaf81qa/Documents/")
write.matrix(B, file = "Matrix_B")
B_ <- read.table("Matrix_B")





