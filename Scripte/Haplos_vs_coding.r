### The purpose of this is to analyze whether a region is part of a coding or a non-coding regeion Input of the following functions are the tables as produced bye the find_haplo.r function set.
##In the first step we need the actual numeric information of the start position of the specif window which by default is only save as a factor class..

### First things first we add a 3 col full of NAs to the data table


Chr_1 <- read.csv("Chr_1_20_1135.csv")

Chr_1 <- cbind(Chr_1,rep(NA,dim(Chr_1)[1]))
###assigin a Name to the new data table 
colnames(Chr_1)[4] <- "Pos"

#now we are going to use the results stored in the first column of the table as input... it should be in a format like: "1- 55" and we need the 55 which is the position in basepaars as a numerical value

for(i in 1:dim(Chr_1)[1]){
    Chr_1$Pos[i]<-as.integer(strsplit(as.character(Chr_1$Position[i]),"-")[[1]][2])
}

# Time < 1 min.
## Alles rauswerfen aus der Tabelle was wir nicht mehr brauch::

Chr_1 <- Chr_1[,3:4]

## Split Ara11 up into the 5 chromsomes shouldnt be too complicated i guess


Ara11_Chr1 <- subset(Ara11, Ara11$V1=="Chr1")
Ara11_Chr2 <- subset(Ara11, Ara11$V1=="Chr2")
Ara11_Chr3 <- subset(Ara11, Ara11$V1=="Chr3")
Ara11_Chr4 <- subset(Ara11, Ara11$V1=="Chr4")
Ara11_Chr5 <- subset(Ara11, Ara11$V1=="Chr5")




### So und nun noch irgendwo hinzufügen, ob es auf einer coding oder einer non coding region liegt... Super eays
##Dazu brauchen wir die Datei Ara11(AraELF), den Nachfolger von Tair10 wo alle bekannten Gene gekennzeichnet sind mit Position
## und ich will jetzt erstmal wissen ob das Anfangs SNP eine Windows welches in der Datei augenommen ist coding region liegt oder nicht... das wäre doch schon mal was :-)

##Die folgenden beiden Funktionen machen genau das gleiche....
##Die erste mach einen großen Vektor und braucht ewigkeiten, da der Vector irgendwan ins gigantische wächste dadurch braucht er ca. ein 12 Stunden

my_vector <- c()
for (i in 1:dim(Ara11_1)[1]){
  my_vector <- append(my_vector, ((Ara11_1$V2[i]):(Ara11_1$V3[i])))
if (i%%10 == 0){
my_vector <- unique(my_vector)
}}

##Die funktion speichert den Vektor alle 1000 iterations ab und löchst den Vektor und beginnt mit Vektorlänge 0 wieder von vorn und braucht dadurch nichteinal eine Mintue..


my_vector <- c()
for (i in 1:dim(Ara11_1)[1]){
  my_vector <- append(my_vector, ((Ara11_1$V2[i]):(Ara11_1$V3[i])))
  if (i%%100 == 0){
    my_vector <- unique(my_vector)
    print(i)
  }
  if (i%%1000 == 0){
    assign(paste('Vec',i,sep=''),my_vector)
    rm(my_vector)
    my_vector <- c()
    
  }
}

#jetz hat man natürlich ganz viel Vektoren die man wieder zusammenbasteln muss-...


Vec_Fin <- c(Vec1000,Vec2000,Vec3000,Vec4000,Vec5000,Vec6000,Vec7000,Vec8000,Vec9000,Vec10000,Vec11000,Vec12000,my_vector)



## achja wieder so unfassbar eleganter Code....

Chr_1 <- cbind(Chr_1,rep(0,dim(Chr_1)[1]))
colnames(Chr_1)[3] <- "Yay"
 

for (i in 1:dim(Chr_1)[1]){
  if((Chr_1$Pos[i]%in%Vec_Fin)==TRUE){
  Chr_1$Yay <- 1
  }
}



##Hier fehlt dann der Code mit der Liste der bissl schneller war,,,,...

my_list <- list()




My_list <- list(Vec1000,Vec2000,Vec3000,Vec4000,Vec5000,Vec6000,Vec7000,Vec8000,Vec9000,Vec10000,Vec11000,Vec12000,my_vector)



###auch hier stest Chr_1 in Chr_2 usw.. äandern
j =1
for (i in 1:dim(Chr_1)[1]){
    if((Chr_1$Pos[i]%in%My_list[[j]])==TRUE){
        Chr_1$Yay[i] <- 1
    }else{
        Chr_1$Yay[i] <- 0
    }
    if(max(My_list[j][[1]])<=Chr_1$Pos[i]){
        j = j+1
        print(j)
  }
}



### Plot den Kram mit ifelse für color auswahl

setwd("/home/ark06eu/evolgen/data/1135_data_Rdata/Ara11/")
dir()
for(u in 1:5){
  filename <- paste('Chr_',u,'_cod.csv',sep='')
  Chr_5 <- read.csv(filename)
  filename <- paste('plot_chr_',u,'.pdf',sep='') 
pdf(filename, width = 12)
plot(Chr_1$Pos_2,Chr_1$Haplos, col = ifelse(Chr_1$Yay == 0, 'red', 'blue'))
dev.off()
}
head(Chr_5)
Chr_5$Yay <- as.character(Chr_5$Yay)
colnames(Chr_5) <- c('X','Haplos','Pos_2', 'coding')

setwd('/home/jaf81qa/Pictures/')
jpeg('Chr_3_cod_gg.jpg', width = 800, height = 600)
g = ggplot(Chr_3, aes(Pos_2,Haplos, color=coding))
g+geom_point() +
  theme(legend.position = c(0.9,0.9))+
  xlab('position (Bp)')+
  ylab('number of haplotypes per 20 SNPs over 1135 accessions')

dev.off()

setwd('/home/ark06eu/evolgen/jan/Haplos_type_scripte/Results_Haplo_Rich_Regs/1035/')

g+geom_point() 


Chr_1 <- read.csv('Ara11_1_HPG.csv')
head(Chr_1)
