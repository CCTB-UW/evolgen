library(dplyr)
library(rafalib)


##### Wie verändert sich die Markerdensity, wenn man bestimmte Marker excluded z.B. alle die geringer freqs als 50 haben oder mehr als 50 haben oder wie auch immer...





Posfreq_1_f <- t(subset(t(Posfreq_1),t(Posfreq_1)[,1]>1000))
Posfreq_2_f <- t(subset(t(Posfreq_2),t(Posfreq_2)[,1]>1000))
Posfreq_3_f <- t(subset(t(Posfreq_3),t(Posfreq_3)[,1]>1000))
Posfreq_4_f <- t(subset(t(Posfreq_4),t(Posfreq_4)[,1]>1000))
Posfreq_5_f <- t(subset(t(Posfreq_5),t(Posfreq_5)[,1]>1000))

Dens_1_f <-  density(Posfreq_1_f[2,])
Dens_2_f <-  density(Posfreq_2_f[2,])
Dens_3_f <-  density(Posfreq_3_f[2,])
Dens_4_f <-  density(Posfreq_4_f[2,])
Dens_5_f <-  density(Posfreq_5_f[2,])

Dens_1 <-  density(Posfreq_1[2,])
Dens_2 <-  density(Posfreq_2[2,])
Dens_3 <-  density(Posfreq_3[2,])
Dens_4 <-  density(Posfreq_4[2,])
Dens_5 <-  density(Posfreq_5[2,])

pdf("plot_all_vs_S1000.pdf")
mypar(5,2)
plot(Dens_1)
plot(Dens_1_f)
plot(Dens_2)
plot(Dens_2_f)
plot(Dens_3)
plot(Dens_3_f)
plot(Dens_4)
plot(Dens_4_f)
plot(Dens_5)
plot(Dens_5_f)
dev.off()








### Häufigkeiten der Marker freqs

len <- max(Posfreq_1[1,])
Haufig_1 <-(rbind(1:len,vector(length=len)))
for (i in 1:(max(Posfreq_1[1,]))){
  Haufig_1[2,i]<- length(which(Posfreq_1[1,]==i))
}

len <- max(Posfreq_2[1,])
Haufig_2 <-(rbind(1:len,vector(length=len)))
for (i in 1:(max(Posfreq_2[1,]))){
  Haufig_2[2,i]<- length(which(Posfreq_2[1,]==i))
}

len <- max(Posfreq_3[1,])
Haufig_3 <- (rbind(1:len,vector(length=len)))
for (i in 1:(max(Posfreq_3[1,]))){
  Haufig_3[2,i]<- length(which(Posfreq_3[1,]==i))
}

len <- max(Posfreq_4[1,])
Haufig_4<- (rbind(1:len,vector(length=len)))
for (i in 1:(max(Posfreq_4[1,]))){
  Haufig_4[2,i]<- length(which(Posfreq_4[1,]==i))
}

len <- max(Posfreq_5[1,])
Haufig_5<- (rbind(1:len,vector(length=len)))
for (i in 1:(max(Posfreq_5[1,]))){
  Haufig_5[2,i]<- length(which(Posfreq_5[1,]==i))
}

mypar(5,1)
plot(Haufig_1[2,50:1134])
plot(Haufig_2[2,50:1134])
plot(Haufig_3[2,50:1134])
plot(Haufig_4[2,50:1134])
plot(Haufig_5[2,50:1134])
