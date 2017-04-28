## IM ERSTEN SCHRITT WOLLEN WIR DIE MARKER DICHTE ANHAND DER PHYSICALISCHEN POSITION BESTIMMEN..
### ERSTMAL DIE DATEI EINLESEN BISSCHEN DOOF, DASS ICH DIE GANZE TABELLE EINLESEN MUSS OBWOHL ICH NUR DIE MARKERNAMEN BRAUCHE...
library(dplyr)
Hi <- Sys.time()
setwd("/home/ark06eu/evolgen/data/1135_data_Rdata/")

### collect all the column names from the data files at first##

for (u in 1:22){
    filename <- paste('X_1135_',u,'.rda.gz',sep = "")
    load(filename)
    if (u==1){
        Pos  <- colnames(X)
        Freq <- apply(X,2,sum)
        rm(X)
    }else {
        Pos   <- append(Pos,colnames(X))
        Freqs <- apply(X,2,sum)
        Freq  <- append(Freq,Freqs)
        rm(X)
        print(u)
    }
}

#### Nun sollten zwei Vectoren vorhanden sein Freq und Pos der eine enthält die Position der Marker und der andere enthält
### die freqenzen der einezelnen Marker diese sollte gleichlang sein ansosnten wird es eine Fehlermeldung geben...,
### um das zu testen habe ich unten einen Test hinzugefügt. Es gibt noch ein paar mehr Tests, wenn man End der Vector
### Tests nur TRUE enthält ist alles soweit gut geangen

Tests <- c()
Testfunc  <- function() {
    if (length(Freq)==length(Pos)){
        Tests <- append(Tests,TRUE)
    }else{
        Tests <- append(Tests,FALSE)
    }
}
Tests <- Testfunc()

#### Okay nun bauen wird die o.g. Vectoren zu Tabellen mit dim(2,X) um für jedes Xsomen eine Tabelle logischerweise
#### Hier sollte alles gut gehen ansonsten wird er sich bescheren, dass die zu rbind()en vectoren nicht gleich lang sind
#### wenn er einen Fehler hat kann das dran liegen, dass im Chr1,2,3,4,5 ein Marker die Pos 1,2,3,4,5
### Hier wäre natürlich noch ein loop besser. In dem man oben einfach die Anzahl der Chromosomen eingibt.


Chr_1 <-subset(Pos,startsWith(Pos,'1')==TRUE)
A1    <- (strsplit(Chr_1,"-") %>% unlist)
A_1   <- as.integer((A1[A1!="1"]))
Posfreq_1 <- rbind(Freq[1:length(A_1)],A_1)
Freq  <- Freq[-c(1:length(A_1))]
rm(Chr_1,A1,A_1)

Tests <- append(Tests,Testfunc())

Chr_2 <-subset(Pos,startsWith(Pos,'2')==TRUE)
A2    <- (strsplit(Chr_2,"-") %>% unlist)
A_2   <- as.integer((A2[A2!="2"]))
Posfreq_2 <- rbind(Freq[1:length(A_2)],A_2)
Freq  <- Freq[-c(1:length(A_2))]
rm(Chr_2,A2,A_2)

Tests <- append(Tests,Testfunc())

Chr_3 <-subset(Pos,startsWith(Pos,'3')==TRUE)
A3    <- (strsplit(Chr_3,"-") %>% unlist)
A_3   <- as.integer((A3[A3!="3"]))
Posfreq_3 <- rbind(Freq[1:length(A_3)],A_3)
Freq  <- Freq[-c(1:length(A_3))]
rm(Chr_3,A3,A_3)

Tests <- append(Tests,Testfunc())

Chr_4 <-subset(Pos,startsWith(Pos,'4')==TRUE)
A4    <- (strsplit(Chr_4,"-") %>% unlist)
A_4   <- as.integer((A4[A4!="4"]))
Posfreq_4 <- rbind(Freq[1:length(A_4)],A_4)
Freq  <- Freq[-c(1:length(A_4))]
rm(Chr_4,A4,A_4)

Tests <- append(Tests,Testfunc())

Chr_5 <-subset(Pos,startsWith(Pos,'5')==TRUE)
A5    <- (strsplit(Chr_5,"-") %>% unlist)
A_5   <- as.integer((A5[A5!="5"]))
Posfreq_5 <- rbind(Freq[1:length(A_5)],A_5)
Freq  <- Freq[-c(1:length(A_5))]
rm(Chr_5,A5,A_5)

Tests <- append(Tests,Testfunc())




Time <- Sys.time()-Hi



### well now there should be five vectores with the postions 
### first things first let's do a simple density plot

Dens_1 <-  density(Posfreq_1[2,])
Dens_2 <-  density(Posfreq_2[2,])
Dens_3 <-  density(Posfreq_3[2,])
Dens_4 <-  density(Posfreq_4[2,])
Dens_5 <-  denstiy(Posfreq_5[2,])

### und nun nur noch die Density plotten und wir haben einen ersten Entwurf
### nun wissen wir was über die absolute Dichte des Vorkommens von SNPs aber die SNP freq muss noch mit einbezogen werden, wenn wir nützuliche Aussagen treffen wollen
Dens_1 <-  density(Posfreq_1[2,])
Dens_2 <-  density(Posfreq_2[2,])
Dens_3 <-  density(Posfreq_3[2,])
Dens_4 <-  density(Posfreq_4[2,])
Dens_5 <-  density(Posfreq_5[2,])

### und nun nur noch die Density plotten und wir haben einen ersten Entwurf
### nun wissen wir was über die absolute Dichte des Vorkommens von SNPs aber die SNP freq muss noch mit einbezogen werden, wenn wir nützuliche Aussagen treffen wollen
library(rafalib)

mypar(5,1)
plot(Dens_1)
plot(Dens_2)
plot(Dens_3)
plot(Dens_4)
plot(Dens_5)
