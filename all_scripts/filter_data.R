###to extract allele frequency information from vcf.snpeff.gz version v4.1
###vcftools --gzvcf 1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff.gz --snp --out out_test.frequencies

load('Ann_gene.rda')
load('tair10.rda')

library(plyr)
library(dplyr)
library(DSR)
library(tidyr)
library(splitstackshape)

fil_data<-function(Summary_1001Genome=read.csv('1001genomes_vcf_summary.csv')){
colnames(Summary_1001Genome)=c('GENE_ID','GENE_NAME','BIO_TYPE','SYNONYMOUS_VARIANT','NON_SYNONYMOUS_VARIANT','START_GAINED','START_LOST','STOP_GAINED','STOP_LOST')
freq_geneid<-read.csv('Allele_frq_all.csv')

############## STOP_GAINED #############

f1=filter(Summary_1001Genome, STOP_GAINED > "0")
gene_list_SG<-f1$GENE_ID
fil1_SG <- filter(tair10, Gene %in% gene_list_SG)
fil2_SG <- filter(anno, gene %in% gene_list_SG)


### filter data according to gene_list_SG

anno_stop<-anno[grep('STOP_GAIN',anno[,6]),]

SNP_list_SG<-anno_stop$SNP    ### take SNP list from anno
freq_stop_gained <- filter(freq_geneid, SNP %in% SNP_list_SG)   #### filter stop gain allele freq by snp list
colnames(freq_stop_gained)[27] <- "Gene"    ### Rename GENE to Gene
d<-merge(tair10,freq_stop_gained,by='Gene')    ### merge start and stop from tair 10
d$ref_start_SG = d$POS - d$Start
d$ref_stop_SG = d$POS - d$Stop
d$total<-d$Stop-d$Start
d$rel_start<-d$ref_start_SG/d$total
d$rel_stop<-d$ref_stop_SG/d$total

AC_ID<-read.delim('c.txt',header =FALSE, sep = '\t')
colnames(AC_ID)=c('SNP','IDs')
d_<-merge(d,AC_ID1,by='SNP')

d<-read.delim('stop_gain_last.txt', sep= '') #### extracts data for all snps
d_<-read.delim('stop_gain_del.txt', sep= '')
write.txt(d_, file = paste('stop_gain_last','.txt',sep = ''))


### gene centric

### how many stop gained SNPS each gene has??
## gene_centric<-count(d_, "Gene")
gene_centric<-as.data.frame(table(d_$Gene))
colnames(gene_centric)<-c('Gene','no_of_stop')
gene_centric_count<-as.data.frame(unique(table(gene_centric[,2])))
colnames(gene_centric_count)<-c('no_of_stop','no_of_genes')
hist(gene_centric_count$no_of_stop)


ID_num<-count(na.omit(a$Gene))
colnames(ID_num)<-c('Gene','count')
gene_table<-merge(gene_centric,ID_num,by="Gene")



##### accession_centric

stop_table<-read.delim('stop_gain_del.txt', sep= '')
gene_data<-as.data.frame(matrix(nrow=28148,ncol=2))
colnames(gene_data)=c('Gene','IDs') 
gene_data$Gene<-stop_table$Gene
gene_data$IDs<-stop_table$IDs
tab<-ddply(gene_data, "Gene", summarize, IDs = paste(IDs, collapse = ", "))



#tab<-as.data.frame(matrix(nrow=28452,ncol=2))
#tab[,1]<-d_$Gene
#tab[,2]<-d_$IDs
#colnames(tab)<-c('Gene','IDs')
#a<-separate_rows(tab, IDs, sep = ", ")

df3 <- cSplit(tab, "IDs", sep = ", ", direction = "long")
pre_idcount<-unique(df3)
ID_Count<-data.frame ( table ( pre_idcount$IDs ) )
colnames(ID_count)<-c('IDs','no_of_stop')
e<-ddply(pre_idcount,.(IDs),summarise,  Gene = paste(unique(Gene),collapse = ', '))
df4<-ddply(pre_idcount,.(Gene),summarise,  IDs = paste(unique(IDs),collapse = ', '))
accession_centric<-merge(ID_Count,e,by="IDs")



#b<-as.numeric(a$IDs)
#ID_Count<-(count(na.omit(b)))
#colnames(ID_Count)<-c('IDs','no_of_stop')


#### 

SNP_Count_SG<-count(freq_stop_gained$GENE)
a_SG<-as.numeric(SNP_Count_SG$freq)

if(length(a_SG)==length(fil1_SG)){
  start_pos_SG<-rep(fil1_SG[,2],a_SG)
  stop_pos_SG<-rep(fil1_SG[,3],a_SG)
}
  else{
    y_SG=unique(freq_stop_gained$GENE)
    z_SG<-fil1_SG$Gene
    missing_genes_SG<-z_SG[-match(y_SG,z_SG)]
    write.table(missing_genes_SG, file=paste('missing_genes_SG','.txt',sep = ''))
    fil1_subset_SG <- fil1_SG[ ! fil1_SG$Gene %in% missing_genes_SG, ]
    start_pos_SG<-rep(fil1_subset_SG[,2],a_SG)
    stop_pos_SG<-rep(fil1_subset_SG[,3],a_SG)
    }
  
stop_gain_combine_table<-cbind(freq_stop_gained,start_pos_SG,stop_pos_SG)
d<-stop_gain_combine_table
d$ref_start_SG = d$POS - d$Start
d$ref_stop_SG = d$POS - d$Stop
head(d)


write.csv(d, file = paste('stop_gain_final','.csv',sep = ''))
write.csv(f1,file=paste('stop_gained','.csv',sep=''))
write.csv(fil2, file = paste('stop_gained_snpeff','.csv',sep = ''))


############## STOP LOST #############


f2=filter(Summary_1001Genome, STOP_LOST > "0")
gene_list_SL<-f2$GENE_ID
freq_stop_lost <- filter(freq_geneid, GENE %in% gene_list_SL)
write.csv(f2,file=paste('stop_lost','.csv',sep=''))


############## SYNONYMOUS_VARIANT #############




f3=filter(Summary_1001Genome, SYNONYMOUS_VARIANT > "0")
gene_list_SV<-f3$GENE_ID
freq_syn_var <- filter(freq_geneid, GENE %in% gene_list_SV)
write.csv(f3,file=paste('synonymous_variant','.csv',sep=''))


############## NON_SYNONYMOUS_VARIANT #############




f4=filter(Summary_1001Genome, NON_SYNONYMOUS_VARIANT > "0")
gene_list_NSV<-f4$GENE_ID
freq_nonsyn_var <- filter(freq_geneid, GENE %in% gene_list_NSV)
write.csv(f4,file=paste('non_synonymous_variant','.csv',sep=''))

############## START_GAINED #############




f5=filter(Summary_1001Genome, START_GAINED > "0")
gene_list_STG<-f5$GENE_ID
freq_start_gained <- filter(freq_geneid, GENE %in% gene_list_STG)
write.csv(f5,file=paste('start_gained','.csv',sep=''))

############## START_LOST #############



f6=filter(Summary_1001Genome, START_LOST > "0")
gene_list_STL<-f6$GENE_ID
freq_start_lost <- filter(freq_geneid, GENE %in% gene_list_STL)
write.csv(f6,file=paste('start_lost','.csv',sep=''))

return(Summary_1001Genome)
}




################# END #################


## Anno<-merge(anno,SNPs[,c(4,6)],by='SNP',all.x=TRUE)

###adding gene list from tair 10 to anno in loop

anno$gene<-NA
for ( i in 1:nrow(tair10)) {
  anno[which(anno[,1]==tair10[i,1]&anno[,2]>tair10[i,2]&anno[,2]<tair10[i,3]),7]<-tair10[i,5]
  cat(i,'\n')
  flush.console()
  save(anno, file = (anno_upd.rda))
}
Anno[which(Anno[,1]==tair10[1,1]&Anno[,2]>tair10[1,2]&Anno[,2]<tair10[1,3]),7]<-tair10[1,5]
tail(Anno)



#### read frequency data from .csv file
freq<-read.csv('out_test.frequencies.csv')
### split whole data into separate columns
freq_SPLITTED<-separate(freq, CHROM.POS.N_ALLELES.N_CHR..ALLELE.FREQ. , into = c('CHR', 'POS', 'N_ALLELES', 'N_CHR','ALLELE.FREQ_1','ALLELE.FREQ_2','ALLELE.FREQ_3', 'ALLELE.FREQ_4', 'ALLELE.FREQ_5','ALLELE.FREQ_6','ALLELE.FREQ_7', 'ALLELE.FREQ_8', 'ALLELE.FREQ_9', 'ALLELE.FREQ_10', 'ALLELE.FREQ_11', 'ALLELE.FREQ_12','ALLELE.FREQ_13','ALLELE.FREQ_14', 'ALLELE.FREQ_15', 'ALLELE.FREQ_16', 'ALLELE.FREQ_17','ALLELE.FREQ_18', 'ALLELE.FREQ_19', 'ALLELE.FREQ_20'),sep='\\t', convert = TRUE)
write.csv(freq_SPLITTED, file = paste('Allele_frq_all','.csv',sep = ''))


###adding gene_list in allele frequency table to filter data 

freq_SPLITTED$SNP<-anno$SNP
freq_SPLITTED$Gene<-anno$gene
#rename(freq_SPLITTED,c('V25'='SNP','V26'='GENE'))





###### COMBINE ALL DATA ACCORDING TO GENE ######


a<-ddply(freq_stop_gained, "gene", summarize, SNP = paste(SNP, collapse = ", "))
b<-ddply(freq_stop_gained, "gene", summarize, CHR = paste(CHR, collapse = ", "))
c<-ddply(freq_stop_gained, "gene", summarize, POS = paste(POS, collapse = ", "))
d<-ddply(freq_stop_gained, "gene", summarize, N_ALLELES = paste(N_ALLELES, collapse = ", "))
e<-ddply(freq_stop_gained, "gene", summarize, N_CHR = paste(N_CHR, collapse = ", "))
f<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_1 = paste(ALLELE.FREQ_1, collapse = ", "))
g<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_2 = paste(ALLELE.FREQ_2, collapse = ", "))
h<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_3 = paste(ALLELE.FREQ_3, collapse = ", "))
i<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_4 = paste(ALLELE.FREQ_4, collapse = ", "))
j<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_5 = paste(ALLELE.FREQ_5, collapse = ", "))
k<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_6 = paste(ALLELE.FREQ_6, collapse = ", "))
l<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_7 = paste(ALLELE.FREQ_7, collapse = ", "))
m<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_8 = paste(ALLELE.FREQ_8, collapse = ", "))
n<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_9 = paste(ALLELE.FREQ_9, collapse = ", "))
o<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_10 = paste(ALLELE.FREQ_10, collapse = ", "))
p<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_11 = paste(ALLELE.FREQ_11, collapse = ", "))
q<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_12 = paste(ALLELE.FREQ_12, collapse = ", "))
r<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_13 = paste(ALLELE.FREQ_13, collapse = ", "))
s<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_14 = paste(ALLELE.FREQ_14, collapse = ", "))
t<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_15 = paste(ALLELE.FREQ_15, collapse = ", "))
u<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_16 = paste(ALLELE.FREQ_16, collapse = ", "))
v<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_17 = paste(ALLELE.FREQ_17, collapse = ", "))
w<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_18 = paste(ALLELE.FREQ_18, collapse = ", "))
x<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_19 = paste(ALLELE.FREQ_19, collapse = ", "))
y<-ddply(freq_stop_gained, "gene", summarize, ALLELE.FREQ_20 = paste(ALLELE.FREQ_20, collapse = ", "))

frq_stop_G<-as.data.frame(cbind(a[,2],b[,2],c[,2],d[,2],e[,2],f[,2],g[,2],h[,2],i[,2],j[,2],k[,2],l[,2],m[,2],n[,2],
                                o[,2],p[,2],q[,2],r[,2],s[,2],t[,2],u[,2],v[,2],w[,2],x[,2],y[,2]))

colnames(frq_stop_G)=c('CHR', 'POS', 'N_ALLELES', 'N_CHR','ALLELE.FREQ_1','ALLELE.FREQ_2','ALLELE.FREQ_3', 'ALLELE.FREQ_4', 'ALLELE.FREQ_5','ALLELE.FREQ_6','ALLELE.FREQ_7', 'ALLELE.FREQ_8', 'ALLELE.FREQ_9', 'ALLELE.FREQ_10', 'ALLELE.FREQ_11', 'ALLELE.FREQ_12','ALLELE.FREQ_13','ALLELE.FREQ_14', 'ALLELE.FREQ_15', 'ALLELE.FREQ_16', 'ALLELE.FREQ_17','ALLELE.FREQ_18', 'ALLELE.FREQ_19', 'ALLELE.FREQ_20')

COMBINE_STOP_GAINED<- cbind(f1,fil1,frq_stop_G)




## to convert all columns from data frame to character
acc<-paste(dQuote(t(a)), collapse=',')

m<-as.numeric(gsub(",", "", d_$IDs)) ### it also adds NAs 
ID_Count<-(count(na.omit(m)))

acc<-read.delim('newfile1.txt',header = FALSE)
t(acc)
o<-unique(na.omit(m))
p<-acc[,1]
miss<-as.numeric(p[-match(o,p)])







