

load('Sub_tair10.Rda')
load('all_chr_net_data.Rda')

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

add_GO_SLIM_term<-function(){
   
  GO_SLIM=read.delim(file = 'ATH_GO_GOSLIM.txt (2)', header = FALSE)
  colnames(GO_SLIM)<-c('locus_name', 'TAIR_accession','object_name', 'relationship_type', 'GO_term', 'GO_ID','TAIR_keyword_ID', 'aspect', 'GO_slim_term', 'evidence_code', 'evidence_description', 'evidence_with', 'reference', 'annotator', 'data_annotated' )
  subset_GO_SLIM<-GO_SLIM
  subset_GO_SLIM<-subset_GO_SLIM[,1:6]
  subset_GO_SLIM[,2:5]<-NULL
  subset_GO_SLIM<-subset_GO_SLIM[1:277537,1:2]
  
  mol_fun<-c('GO:0016787','GO:0016301','GO:0016740','GO:0003824','GO:0003700','GO:0003677','GO:0003723','GO:0003676','GO:0000166','GO:0005515','GO:0005102','GO:0004872','GO:0005488','GO:0005198','GO:0005215','GO:0005554','GO:0003674')
  
  t1<-filter(subset_GO_SLIM, GO_ID %in% mol_fun)
  
  data_subset<-subset(all_chr_net_data,(all_chr_net_data$p_value_under < 10^-9))
  
  suba<-data.frame ( table ( data_subset$gene1 ) )
  subb<-data.frame ( table ( data_subset$gene2 ) )
  s_a<-suba[!(suba$Freq == 0),]
  s_b<-subb[!(subb$Freq == 0),]
  a<-rbind(s_a, s_b)
  x<-unique(a$Var1)
  
  sb_tair_g1<-sub_tair10[,1:4]
  sb_tair_g2<-sub_tair10[,5:7]
  tair_g1<-filter(sb_tair_g1, gene1 %in% x )
  tair_g2<-filter(sb_tair_g2, gene2 %in% x )
  
  incl_dist<-merge(data_subset,tair_g2,by='gene2')
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
  unique_gene_list<-unique(a$Var1)
  d<-filter(t1, locus_name %in% unique_gene_list)
  
  graph_table<-as.data.frame(matrix(nrow=15, ncol=6))
  rownames(graph_table)<-c('hydrolase_activity', "kinase_activity" , "transferase_activity" , "other_enzyme_activity" , "transcription_factor_activity", "DNA_or_RNA_binding", "nucleic_acid_binding" , "nucleotide_binding",  "protein_binding" , "receptor_binding_or_activity", "other_binding", "structural_molecule_activity", "transporter_activity", "unknown_molecular_functions", "other_molecular_functions")

  graph_table[1,2]<-count(filter(d, GO_ID == 'GO:0016787'))
  graph_table[2,2]<-count(filter(d, GO_ID == 'GO:0016301'))
  graph_table[3,2]<-count(filter(d, GO_ID == 'GO:0016740'))
  graph_table[4,2]<-count(filter(d, GO_ID == 'GO:0003824'))
  graph_table[5,2]<-count(filter(d, GO_ID == 'GO:0003700'))
  graph_table[6,2]<-count(filter(d, GO_ID == 'GO:0003677')) + count(filter(d, GO_ID == 'GO:0003723'))
  graph_table[7,2]<-count(filter(d, GO_ID == 'GO:0003676'))
  graph_table[8,2]<-count(filter(d, GO_ID == 'GO:0000166'))
  graph_table[9,2]<-count(filter(d, GO_ID == 'GO:0005515'))
  graph_table[10,2]<-count(filter(d, GO_ID == 'GO:0005102')) + count(filter(d, GO_ID == 'GO:0004872'))
  graph_table[11,2]<-count(filter(d, GO_ID == 'GO:0005488'))                  
  graph_table[12,2]<-count(filter(d, GO_ID == 'GO:0005198'))
  graph_table[13,2]<-count(filter(d, GO_ID == 'GO:0005215'))
  graph_table[14,2]<-count(filter(d, GO_ID == 'GO:0005554'))
  graph_table[15,2]<-count(filter(d, GO_ID == 'GO:0003674'))

  dat2 = data.frame(count=graph_table[,2], category=rownames(graph_table))
  dat2$percentage = dat2$count / sum(dat2$count)* 100
  dat2$fraction = dat2$count / sum(dat2$count)
  dat2 = dat2[order(dat2$fraction), ]
  dat2$ymax = cumsum(dat2$fraction)
  dat2$ymin = c(0, head(dat2$ymax, n=-1))
  

  donut = ggplot(dat2, aes(fill = category, ymax = ymax, ymin = ymin, xmax = 100, xmin = 80)) +
  geom_rect(colour = "black") +
  coord_polar(theta = "y") + 
  xlim(c(0, 100)) +
  geom_label_repel(aes(label = paste(round(percentage,2),"%"), x = 100, y = (ymin + ymax)/2),inherit.aes = F, show.legend = F, size = 5)+
  theme(legend.title = element_text(colour = "black", size = 16, face = "bold"), 
          legend.text = element_text(colour = "black", size = 15), 
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank()) +
  annotate("text", x = 0, y = 0, size = 6, label = "seg_under9_20kb") ####### change label of graph according to requirement
  donut
  return(donut)  
}  
  ## to save plot as pdf
  pdf(file='seg_over20_50kb.pdf') ####### change name of graph according to requirement
  donut
  dev.off()

