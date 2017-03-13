Remove_Triall <- function(Data){
  Triall = vector(mode ="numeric") 
  for (i in 1:dim(Data)[2]){
    if ((is.na(levels(Data[,i])[4]))==FALSE){
      Triall= append(Triall,i)  
    }
  }
  Data_OT <- Data[,-Triall]
}

Ref_Allel <- function(Data_OT){
  Ref = vector(length = 0)
  for (i in 1:ncol(Data_OT)){
    a <- which(Data_OT[,i]==levels(Data[,i])[1])
    b <- which(Data_OT[,i]==levels(Data_OT[,i])[2])
    c <- which(Data_OT[,i]==levels(Data_OT[,i])[3])
    if (max(length(a),length(b), length(c))==length(a)){
      Ref <- append(Ref,levels(Data_OT[,i])[1])
    } else if (max(length(a),length(b), length(c))==length(b)){
      Ref <- append(Ref,levels(Data_OT[,i])[2])
    } else {
      Ref <- append(Ref,levels(Data_OT[,i])[3])
    }  
    if (Ref[i]=="./."){
      if (levels(Data_OT[,i])[1]!= "./."){
        Ref[i] <- levels(Data_OT[,i])[1] 
      } else if (levels(Data_OT[,i])[1]== "./.") {
        Ref[i] <- levels(Data_OT[,i])[2]
      }
    }
  }
  Ref
}

Transform_Data_Fast <- function(Data_OT,Ref){
  for (i in 1:ncol(Data_OT)){
    Data_OT[,i] <- ifelse(Data_OT[,i]%in%Ref[i],0,ifelse(Data_OT[,i]=="./.",0.5,1))
  }
  Data_OT <- as.matrix(Data_OT)
  Data_OT
}

Find_Haplos_2 <-function(ws,A,raus=0.5){
  Haplos <- numeric(0)
  Position <- numeric(0)
  Position_NAO <- numeric(0)
  Haplos_NAO <- numeric(0)
  for(i in seq(0,(ncol(A)/ws)-1, 0.5)){ 
    A_ <- A[,(i*ws+1):((i+1)*ws)]
    A__ <- which(!duplicated(A_)) 
    Number_Haplos <-  length(A__) 
    AA <- matrix(ncol = ws, nrow = length(A__)) 
    for (k in 1:length(A__)){
      AA[k,] <-(A_[A__[k],])
    }
    rm <- c()
    for (r in 1:nrow(AA)){
     if ((length(which(AA[r,]==0.5))/ws) > raus){
      rm <- append(rm,r)  
      }
    }
    if (length(rm)>0){
    AA <- AA[-rm,]
    }
    NAS <- c()
    for (q in 1:nrow(AA)){
      NAS <- append(NAS, length(which(AA[q,]==0.5)))
    }
    AA <- cbind(matrix(NAS,ncol=1),AA)
    AA <- AA[order(AA[,1]),]
    AA <- AA[,-1]
    Haplos <- append(Haplos,Number_Haplos) 
    Position<- append(Position,i)
    AAA <- AA
    a = 1
    while(a<=(dim(AAA)[1])-1){
      b =1+a
      r= c()
      while(b<dim(AAA)[1]){
        if (all.equal(AAA[a,],AAA[b,],tolerance = 0.5)==TRUE && a!=b){
          r <- append(r,b)
        }
        b = b+1
      }
      if (length(r) >=1){
        AAA <- AAA[-r,]
      }
      a= a+1
    }
    #if (i%%5 ==0){#
    #print(i)#
    #}#  
    Haplos_NAO <- append(Haplos_NAO,nrow(AAA))
    Position_NAO <- append(Position_NAO, i)
  }
  #pdf(file='Haplotypes_Plot.pdf')#
  #plot(Position_NAO,Haplos_NAO)#
  #dev.off()
  print(Haplos_NAO)
  Name <- paste("Haplos_NAO", v, sep = "")
  assign(Name, Haplos_NAO)
  return(Haplos_NAO)
  #Haplos_NAO_Gesamt <- append(Haplos_NAO_Gesamt,Haplos_NAO)#
  #save(Haplos_NAO_Gesamt,file='output_function.rda')#
}

Find_Haplos_Imputed <- function(ws,A,raus=0.5){
  Haplos <- numeric(0)
  Position <- numeric(0)
  for(i in seq(0,(ncol(A)/ws)-1, 0.5)){ 
    A_ <- A[,(i*ws+1):((i+1)*ws)]
    A__ <- which(!duplicated(A_)) 
    Number_Haplos <-  length(A__) 
    Haplos <- append(Haplos,Number_Haplos) 
    Position<- append(Position,i)
  }
    return(Haplos)
}

Find_Haplos_PLYR <- function(ws,A){
  Position <- NULL
  Haplos <- numeric(0)
  for(i in seq(0,(ncol(A)/ws)-1, 0.5)){ 
    A_ <- A[,(i*ws+1):((i+1)*ws)]
    A__ <-distinct(as.data.frame(A_)) 
    Number_Haplos <-  (dim(A__)[1])
    Haplos <- append(Haplos,Number_Haplos) 
    Position<- append(Position,colnames(A)[i*ws+1])
  }
  Results <- data.frame(Position,Haplos)
  return(Results)
}


