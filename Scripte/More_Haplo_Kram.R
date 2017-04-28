setwd("/home2/jaf81qa/Documents/Results Haplos/2029/Results_200/")


write.csv(Chr_5_200,"Chr_5_200.csv")

library(rafalib)

Chr_1_200 <- read.csv("Chr_1_200.csv")
Chr_2_200 <- read.csv("Chr_2_200.csv")
Chr_3_200 <- read.csv("Chr_3_200.csv")
Chr_4_200 <- read.csv("Chr_4_200.csv")
Chr_5_200 <- read.csv("Chr_5_200.csv")
### Ungerade und Gerade überlappen sich nicht mehr##
Gerade <- NULL
Ungerade <- NULL
for (i in 1:dim(Chr_5_200)[1]){
  if (i%%2 == 0){
    Gerade <- append(Gerade,(Chr_5_200[i,3]) )
  } else {
    Ungerade <- append(Ungerade,(Chr_5_200[i,3]))
  }  
}

if (!length(Gerade)== length(Ungerade){
Ungerade <- Ungerade[-length(Ungerade)]
}
### Zwei vecs aus gerade und ungerade machen... um zu gucken
### wie die benachbarten miteinanden cor.

Ung1 <- NULL
Ung2 <- NULL
for (i in 1:length(Ungerade)){
  if (i%%2 == 0){
    Ung1 <- append(Ung1,Ungerade[i])
  } else {
    Ung2 <- append(Ung2,Ungerade[i])
  }
}
if((length(Ung1)==length(Ung2))== FALSE){
Ung2 <- Ung2[-length(Ung2)]  
}


Ger1 <- NULL
Ger2 <- NULL
for (i in 1:length(Gerade)){
  if (i%%2 == 0){
    Ger1 <- append(Ger1,Gerade[i])
  } else {
    Ger2 <- append(Ger2,Gerade[i])
  }
}
if((length(Ger1)==length(Ger2))== FALSE){
  Ger2 <- Ger2[-length(Ger2)]  
}

cor(Ung1,Ung2)
cor(Ger1,Ger2)

mypar(1,3)
plot(Ger1,Ger2)
plot(Ung1,Ung2)
plot(Gerade,Ungerade)

pdf("Cor_Nachbarn_Chr_5_200")
mypar(1,3)
plot(Ungerade,Gerade, main = c("Cor_Überlappung",cor(Ungerade,Gerade)))
abline(lm(Ungerade~Gerade))

plot(Ung1,Ung2, main = c("Cor_Nachbarn_1",cor(Ung1,Ung2)))
abline(lm(Ung1~Ung2))

plot(Ger1,Ger2, main = c("Cor_Nachbarn_2",cor(Ger1,Ger2)))
abline(lm(Ger1~Ger2))

dev.off()

pdf("Haplos_200")
mypar(5,1)
plot(Chr_1_200$Haplos)
plot(Chr_2_200$Haplos)
plot(Chr_3_200$Haplos)
plot(Chr_4_200$Haplos)
plot(Chr_5_200$Haplos)
dev.off()
