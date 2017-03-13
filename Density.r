Marker_Position
dim(Marker_Position)
length(Marker_Position)
Mark_Pos_Sub <- Marker_Position[1:1000]

Marker_Density <- function(Marker_Position,ws){
  Dens <- c() 
for (i in 0:(length(Marker_Position)/ws)){
  Lang <- length(which(Marker_Position>(i*ws+1) & Marker_Position<(i*ws+ws)))
  Dens <- append(Dens,Lang)
}
Dens
plot(seq(1,length(Dens)),Dens)
cor(seq(1,length(Dens)),Dens)
}

Mark_Pos_Sub

ws = 100
i = 3
length(which(Mark_Pos_Sub>(i*ws+1) & Mark_Pos_Sub<(i*ws+ws)))

Marker_Density(Marker_Position,200)
rm(Dens)
.1 == .3/3
unique(c(.3, .4 - .1, .5 - .2, .6 - .3, .7 - .4))
Quad <- function (a, b, c){
  rad <- b^2 - 4 * a * c
  if(is.complex(rad) || all(rad >= 0)) {
    rad <- sqrt(rad)
  } else {
    rad <- sqrt(as.complex(rad))
  }
  cbind(-b - rad, -b + rad) / (2 * a)
}
print(Quad(1/3,-5/3,6/3),digits= 19)

