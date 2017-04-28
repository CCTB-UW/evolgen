

Chr_1
X_1
X_2
X_3
X_4
X_5
X_6[,1:97826]

Chr2

X6[,97827:dim(X)[2]]
X7
X8
X9[,1:466695]

Chr3

X9[,466696:dim(X)[2]]
X10
X11
X12
X13
X14[,1:161065]

Chr4

X14[,161066:dim(X)[2]]
X15
X16
X17[,1:428148]

Chr5

X17[,428149:dim(X)[2]]
X18
X19
X20
X21
X22


### generate 5 empty matricies for each Chroomsome
for(i in 1:5){
  assign(paste('Chr_',i,sep=''),matrix())
}


for(u in 1:22){
    filename <- paste('X_1135_', u, '.rda.gz', sep = '')
    load(filename)
    ## if not otherwise specified the load function will the desired file
    ## as 'X' in the R environment
    if(u == 1){
        Chr_1 <- X
    }else if(u %in% 2:5){
        Chr_1 <- cbind(Chr_1, X)
    }else if(u == 6){
        Chr_1 <- cbind(Chr_1,X[,1:97826])
        Chr_2 <- X[,97827:dim(X)[2]]
    }else if(u %in% 7:8){
        Chr_2 <- cbind(Chr_2,X)
    }else if(u == 9){
        Chr_2 <- cbind(Chr_2,X[,1:466695])
        Chr_3 <- X[,466696:dim(X)[2]]
    }else if (u %in% 10:13){
        Chr_3 <- cbind(Chr_3,X)
    }else if (u == 14){
        Chr_3 <- cbind(Chr_3,X[,1:161061])
        Chr_4 <- X[,161061:dim(X)[2]]
    }else if (u %in% 15:16){
        Chr_4 <- cbind(Chr_4, X)
    }else if (u == 17){
        Chr_4 <- cbind(Chr_4,X[,1:428148])
        Chr_5 <- X[,428149:dim(X)[2]]
    }else {
        Chr_5 <- cbind(Chr_5,X)
    }
}
                       
    
        
        
