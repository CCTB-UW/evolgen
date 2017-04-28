
Chr_1 <- function(){
    setwd('/home/ark06eu/evolgen/data/1135_data_Rdata/')
    for(u in 1:6){
        filename <- paste('X_1135_',u,'.rda.gz',sep='')
        load(filename)
        assign(paste('X',u,sep=''),X)
        }
        C_1 <- cbind(X1,X2,X3,X4,X5,X6[,1:97826])
        return(C_1)
    }

    
Chr_2 <- function(){
    setwd('/home/ark06eu/evolgen/data/1135_data_Rdata/')
    for(u in 6:9){
        filename <- paste('X_1135_',u,'.rda.gz',sep='')
        load(filename)
        assign(paste('X',u,sep=''),X)
        }
        C_2 <- cbind(X6[,97827:dim(X)[2]],X7,X8,X9[,1:466695])

        return(C_2)
    }

    
Chr_3 <- function(){
    setwd('/home/ark06eu/evolgen/data/1135_data_Rdata/')
    for(u in 9:14){
        filename <- paste('X_1135_',u,'.rda.gz',sep='')
        load(filename)
        assign(paste('X',u,sep=''),X)
        }
        C_3 <- cbind(X9[,466696:dim(X)[2]],X10,X11,X12,X13,X14[,1:161060])
        reutrn(C_3)
    }

    
Chr_4 <- function(){
    setwd('/home/ark06eu/evolgen/data/1135_data_Rdata/')
    for(u in 14:17){
        filename <- paste('X_1135_',u,'.rda.gz',sep='')
        load(filename)
        assign(paste('X',u,sep=''),X)
        }
        C_4 <- cbind(X14[,161061:dim(X)[2]], X15, X16, X7[,1:428148])
        return(C_4)
    }

    
Chr_5 <- function(){
    setwd('/home/ark06eu/evolgen/data/1135_data_Rdata/')
    for(u in 17:22){
        filename <- paste('X_1135_',u,'.rda.gz', )
        load(filename)
        assign(paste('X',u,sep=''),X)
        }
        C_5 <- cbind(X17[,428149:dim(X)[2]],X18,X19,X20,X21,X22)
        return(C_5)
    }

    
