setwd("F:/工作/3.诸葛菜/1.Plant_diversity/0.data")

library("DESeq2")
#library(dplyr)
library(tidyverse)
###read counts import to matrix
#all_counts<- read.csv("/usr_storage/forR/jcf/input/all_gene_count.csv",header=T, sep = ",")
library(pheatmap)

meta<-read.csv("ALL_sample_data.csv", header = T)
flav_KEGG<-read_tsv("flavor.txt", col_names = F)
c<-right_join(meta,flav_KEGG,by=c("cpd_ID"="X1"))


meta<-meta %>%
  select(Index, Class.I,Class.II,grep("^T",colnames(.), value = T)) %>%
  arrange(Class.I, Class.II)

rownames(meta)<-meta$Index  
meta<-meta[,-1]

meta %>% select(grep("^T",colnames(.), value = T))

dist.obs.tis<-as.dist(1-cor(t( meta %>% select(grep("^T",colnames(.), value = T)))))

dist.obs.tis.tre<- hclust(dist.obs.tis, method = "ward.D")


ann_colors<-list(Class.I=c('Alkaloids'='#a663cc',  "Amino acids and derivatives"= '#023e8a',
                           "Flavonoids"='#c72c2c',
                           "Lignans and Coumarins"='#fdc500',
                           "Lipids"='#7FFF7F',
                           "Nucleotides and derivatives"= '#ff7b00',
                           "Organic acids"='#ce4257',
                           "Phenolic acids"= '#a2d6f9',
                           "Tannins"='#b08968',
                           "Terpenoids"='#227c9d', 
                           "Others"='#fdf0d5'))


pheatmap(mat = meta %>% select(grep("^T",colnames(.), value = T)) , 
         annotation_row =  meta %>% select(Class.I),
         show_rownames = F, cluster_rows = dist.obs.tis.tre, cluster_cols = F, scale = "row",
         #color  = colorRampPalette(c("#0000FF", "white", "#FF1010"))(50),
         color  = colorRampPalette(c("#020763", "white", "#ad0508"))(50),
         #color  = colorRampPalette(c("#8a95a7", "white", "#ae545e"))(50),
         annotation_colors = ann_colors,
         cutree_rows = 9, treeheight_row = 25,
         gaps_col = c(4,8,12)
)  


dist.obs.tis.sam<-as.dist(1-cor( meta %>% select(grep("^T",colnames(.), value = T))))

dist.obs.tis.tre.sam<- hclust(dist.obs.tis.sam, method = "ward.D")

plot(dist.obs.tis.tre.sam)


meta.pca<-prcomp(t(meta %>% select(grep("^T",colnames(.), value = T))),scale. = T)
    
meta.sum=summary(meta.pca)

library(factoextra)
group=c(rep("T1",4),rep("T2",4),rep("T3",4),rep("T4",4))

fviz_pca_ind(meta.pca,  pointshape = 21,
             col.ind=group,
             pointsize = 2,
             fill=group,
             geom = "point",
             mean.point=F,
             addEllipses = T, 
             legend.title="Groups",
             ellipse.type="confidence",
             ellipse.level=0.95,
             #palette = c("#CC3333", "#339999"))+ #Cell配色哦
             #palette = c("#b390b1", "#acc294" ,"#a1bfbf", "#b9a395")
             palette = "Dark2"
             )+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) 




#all_counts<- read_tsv("Ovio.counts")

all_tpm<-as.data.frame(read_tsv("Ovio_tpm.txt"))
rownames(all_tpm)<-all_tpm$geneid
all_tpm<-all_tpm[,-1]
all_tpm<-all_tpm%>%select(grep("P0",colnames(.), value = T, invert = T)) %>% filter_all(any_vars(.>1)) 


dist.obs.tis<-as.dist(1-cor(t( all_tpm )))
dist.obs.tis.tre<- readRDS("dist.tre.RDS")
#saveRDS(dist.obs.tis.tre, file = "dist.tre.RDS" )


ann_colors<-list(class=c('#c72c2c','#023e8a','#a663cc','#fdc500','#7FFF7F','#ff7b00','#ce4257','#a2d6f9'))

pheatmap(mat = all_tpm , 
         annotation_col = data.frame( row.names = colnames(all_tpm), class=  colnames(all_tpm) %>% gsub("OV_","",.) %>% gsub("\\..*","",.) %>% gsub("A","",.)   )    ,
         show_rownames = F, cluster_rows = dist.obs.tis.tre, cluster_cols = F, scale = "row",
         #color  = colorRampPalette(c("#0000FF", "white", "#FF1010"))(50),
         color  = colorRampPalette(c("#020763", "white", "#ad0508"))(50),
         #color  = colorRampPalette(c("#8a95a7", "white", "#ae545e"))(50),
         annotation_colors = ann_colors,
        
         #cutree_rows = 9, treeheight_row = 25,
         gaps_col = c(3,6,9,12,16,20,24)
)  




dist.obs.tis.sam<-as.dist(1-cor( all_tpm   ))

dist.obs.tis.tre.sam<- hclust(dist.obs.tis.sam, method = "ward.D")

plot(dist.obs.tis.tre.sam)


#all_tpm<-all_tpm%>%select(grep("P0",colnames(.), value = T, invert = T)) %>% filter_all(any_vars(.>1)) 

tpm.pca<-prcomp(t(all_tpm ),scale. = T)
colnames(all_tpm)
tpm.pca.sum=summary(tpm.pca)

#library(factoextra)
#tpm.group=c(  rep("F",4),  rep("Fr",1),  rep("L",4),  rep("R",4),  rep("S",4), rep("T1",4),rep("T2",4),rep("T3",4),rep("T4",4))

tpm.group=c(  rep("F",3),  rep("L",3),  rep("R",3),  rep("S",3), rep("T1",4),rep("T2",4),rep("T3",4),rep("T4",4))

fviz_pca_ind(tpm.pca,  pointshape = 21, geom = c("point","text"),
             #col.ind=tpm.group,
             pointsize = 2,
             fill=tpm.group,
             #geom = "point",
             mean.point=F,
             addEllipses = T, 
             legend.title="Groups",
             ellipse.type="confidence",
             ellipse.level=0.90,
             #palette = c("#CC3333", "#339999"))+ #Cell配色哦
             #palette = c("#b390b1", "#acc294" ,"#a1bfbf", "#b9a395")
            # palette = "Dark2"
  )  +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) 


