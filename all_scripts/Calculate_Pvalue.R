


library(plyr)
library(dplyr)
library(DSR)
library(tidyr)
library(data.table)


pvalue<-function(){


a<-separate_rows(tab, IDs, sep = ", ")
###to get which gene is knocked out in which IDs?
gene_ID<-ddply(a, "Gene", summarize, IDs = paste(IDs, collapse = ", "))


###which genes have been knocked out in how many IDs?
accession_centric<-read.csv('accession_centric.csv',sep = ',')
gen_data<-accession_centric
gen_data$no_of_stop<-NULL
gen_data$X<-NULL

x<-separate_rows(gen_data, Gene, sep = ",")
x<-cSplit(gen_data, 'Gene', "," , "long")
## z<-count(x$Gene)
z<-data.frame ( table ( x$Gene ) )
           ### OR ###


#### to extract bigger gene knockouts
Z<-subset(z,z$freq>100)


####Creating a dataframe to collect all values and calculating p value


net<-read.csv('Net_data_com.csv', sep = ',')
net$X<-NULL
test1<-Q
test1[,4]<-NA
test1[,5]<-NA
test1[,6]<-NA
test1[,7]<-NA
colnames(test1)<-c('gene1','gene2','connection','gene1_count','gene2_count','p_value_under','p_value_over')

test1[,8]<-NA
test1[,8]<-test1[,2]
test1[,2]<-test1[,3]
test1[,3]<-test1[,8]
test1[,8]<-NULL


z<-read.csv('gene_in_acc.csv', sep = ',')
z$X<-NULL
z1<-z
colnames(z)<-c('gene1', 'freq')
colnames(z1)<-c('gene2', 'freq')
test2<-merge(test1,z1,by='gene2')
test3<-merge(test2,z,by='gene1')


gene_1<-test1$gene1
gene_2<-test1$gene2

fil_1 <- filter(z, x %in% gene_1)
fil_2 <- filter(z, x %in% gene_2)
colnames(fil_1)<-c('gene1','num1')
colnames(fil_2)<-c('gene2','num2')
test2<-merge(test1,fil_2,by='gene2')
test3<-merge(test2,fil_1,by='gene1')
test3$gene1_count<-test3$num1
test3$gene2_count<-test3$num2
test3$num1<-NULL
test3$num2<-NULL
write.csv(test3, file = paste('Pre_pvalue.csv', sep = ''))


#### to calculate p value via hypergeometric distribution

connection<-test3$connection
gene1_count<-test3$gene1_count
gene2_count<-test3$gene2_count

test3$p_value_over<-phyper(connection,gene1_count,1135-gene1_count,gene2_count)
test3$p_value_under<-phyper(connection,gene1_count,1135-gene1_count,gene2_count,lower.tail=FALSE)

test4<-test3
write.csv(test4, file = paste('Post_pvalue.csv', sep = ''))

return()
}


### Key to understand ###

#pop_size<-1135
#sample size : gene2_count
#Number of items in the pop that are classified as successes : gene1_count
#Number of items in the sample that are classified as successes : connection





#### Other Methods to calculate Pvalue ### (take time)



#pval<-read.csv("Pre_pvalue.csv", sep=",")
#A<-pval
#A$small<-apply(A[,4:5],1,min)
#A$big<-apply(A[,4:5],1,max)
#A$pattern<-paste(A$connection,A$gene1_count,A$gene2_count,sep='-')
#A$pattern2<-paste(A$connection,A$small,A$big,sep='-')
#D<-unique(A$pattern2)
#sub_pval<-as.data.frame(matrix(nrow=2039579,ncol=6))
#colnames(sub_pval)=c("pattern2",'connection','gene1_count','gene2_count','p_value_over','p_value_under')
#sub_pval$pattern2<-D
#library(stringr)
#test1<-str_split_fixed(sub_pval$pattern2, "-", 3)




#a<-1:1135
#for (t in 1:nrow(sub_pval)) {
#  d<-list()
#  for ( i in 1:1000) {
#    b<-sample(a,sub_pval[t,3])
#    c<-sample(a,sub_pval[t,4])
#    d[i]<-length(which(b%in%c))
#  }
#  sub_pval[t,5]<-length(which(d>=sub_pval[t,2]))/1000
#  sub_pval[t,6]<-length(which(d<=sub_pval[t,2]))/1000
#  }
#  if(t%%10000==0){print(t)}
#}

#test4<-test3
write.csv(test4, file = paste('Post_pvalue.csv', sep = ''))




### Key to understand ###

#pop_size<-49985001
#sample size : 131
#Number of items in the pop that are classified as successes : 1998
#Number of items in the sample that are classified as successes : 62






