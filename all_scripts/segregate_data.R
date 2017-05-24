
load('Sub_tair10.Rda')
load('all_chr_net_data.Rda')

library(plyr)
library(dplyr)
library(tidyr)

data_seg<-function(){
  
  file_name<-subset(all_chr_net_data,(all_chr_net_data$p_value_over < 10^-9))  # pval_col #### ex p_value_over #### change -log10 value
   
  
  sb_tair_g1<-sub_tair10[,1:4]
  sb_tair_g2<-sub_tair10[,5:7]


  suba<-data.frame ( table ( file_name$gene1 ) )
  subb<-data.frame ( table ( file_name$gene2 ) )
  s_a<-suba[!(suba$Freq == 0),]
  s_b<-subb[!(subb$Freq == 0),]
  a<-rbind(s_a, s_b)
  x<-unique(a$Var1)
  
  tair_g1<-filter(sb_tair_g1, gene1 %in% x )
  tair_g2<-filter(sb_tair_g2, gene2 %in% x )
  
  incl_dist<-merge(file_name,tair_g2,by='gene2')
  S_dist<-merge(incl_dist,tair_g1,by='gene1')
  
  S_dist$abs_dist<- abs(S_dist$Stop - S_dist$Start2)

  ### excluding very close same chr data 
  W<-S_dist
  W$samechr<-rep(0)
  W$samechr<-as.numeric(W$Chr)-as.numeric(W$Chr2)
  W$samechr[which(W$samechr<0)]<-1
  remove<-which(W$abs_dist<100000&W$samechr==0)
  W_<-W[-remove,]
  
  over_kb<-filter(W_, abs_dist  > '20000') ## change threshold value at over_....kb and abs_dist > '.....'
  
  suba<-data.frame ( table (over_kb$gene1 ) )
  subb<-data.frame ( table (over_kb$gene2 ) )
  s_a<-suba[!(suba$Freq == 0),]
  s_b<-subb[!(subb$Freq == 0),]
  a<-rbind(s_a, s_b)
  result<-unique(a$Var1)
  result_<-length(result)

  return(result_) 
 }
























