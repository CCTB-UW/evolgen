nPheno = 1000
nSNP = 10000

A <- matrix(nrow = nPheno, ncol = nSNP)
for (i in 1:nPheno){
  A[i,] = round(rnorm(nSNP, mean=0.5, sd = 0.1),0)
}
dim(A)

Find_Haplos <- function(ws,A){
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

Find_Haplos(10,A)


