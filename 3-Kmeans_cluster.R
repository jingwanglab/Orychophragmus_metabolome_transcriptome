

###
###Step1: meta and tpm data pre-processing 
###

setwd("C:/Users/fjj/Desktop/jcf/PD/1.kmeans")
library(tidyverse)
#library("Mfuzz")
#devtools::install_github("junjunlab/ClusterGVis")

library(ClusterGVis)
seed<-13

meta<-read.csv("ALL_sample_data.csv", header = T)
#flav_KEGG<-read_tsv("flavor.txt", col_names = F)
#c<-right_join(meta,flav_KEGG,by=c("cpd_ID"="X1"))

#read data

meta<-meta %>%
  select(Index, Class.I,Class.II,grep("^T",colnames(.), value = T)) %>%
  arrange(Class.I, Class.II)

rownames(meta)<-meta$Index  
meta<-meta[,-1]

meta<-meta %>% select(grep("^T",colnames(.), value = T))
meta<-meta%>%
  mutate(T1=(T1.1+T1.2+T1.3+T1.4)/4 ,T2=(T2.1+T2.2+T2.3+T2.4)/4 , T3=(T3.1+T3.2+T3.3+T3.4)/4, T4=(T4.1+T4.2+T4.3+T4.4)/4  ) %>%
  select(T1,T2,T3,T4)


all_tpm<-as.data.frame(read_tsv("Ovio_tpm.txt"))
rownames(all_tpm)<-all_tpm$geneid
all_tpm<-all_tpm[,-1]
all_tpm<-all_tpm%>%select(grep("P0",colnames(.), value = T, invert = T)) %>% 
  select( grep("OV_T",colnames(.), value = T)  )  %>%
  mutate(T1=(OV_T1A.1+OV_T1A.2+OV_T1A.3+OV_T1A.4)/4 ,T2=(OV_T2A.1+OV_T2A.2+OV_T2A.3+OV_T2A.4)/4 , 
         T3=(OV_T3A.1+OV_T3A.2+OV_T3A.3+OV_T3A.4)/4, T4=(OV_T4A.1+OV_T4A.2+OV_T4A.3+OV_T4A.4)/4  ) %>%
  select(T1,T2,T3,T4) %>% filter_all(any_vars(.>3))


###
###Step2 : Kmeans or Cmeans cluster construction based on visCluster using a merged scale
###

#merge
seed<-25

cm<-clusterData(exp=rbind(meta,all_tpm), cluster.method = "mfuzz", cluster.num = 10)
visCluster(object = cm,
           plot.type = "line")
ck<-clusterData(exp=rbind(meta,all_tpm), cluster.method = "kmeans", cluster.num = 4, seed = 33)
visCluster(object = ck,
           plot.type = "line")

ckk<-ck
###
###Step 3: Evaluation for the quality of cluster classification using a specific seed
###

FA<-read_tsv("Ov_FA.id")
FLA<-read_tsv("Ov_FLA.id")

ck[["wide.res"]]
meta<-read.csv("ALL_sample_data.csv", header = T)
cc<-meta
cc<- cc%>%filter( Class.I=="Flavonoids"|Class.I=='Lipids' )

FLA_id<-read_tsv("flavor.txt", col_names = F) %>%
  left_join(meta,by=c("X1"="cpd_ID")) %>%
  select(Index, Compounds) %>% 
  na.omit()


#left_join( FLA_id, ck[["wide.res"]], by=c("Index"="gene") )
#c<-left_join( FLA, ck[["wide.res"]], by=c("geneID"="gene") )

FLA_all<-rbind(
FLA %>%
  mutate(X1=paste(V1,"_",V2,sep="")) %>%
  select(geneID,X1),
FLA_id %>%
  rename("X1"="Compounds", "geneID"="Index" ) 
)


FLA_all2<-rbind( FA %>% select(geneID, X1) %>% mutate(Class.II="FA",Class.I="FA")  ,
  rbind(
  FLA %>%
    mutate(X1=paste(V1,"_",V2,sep="")) %>%
    select(geneID,X1) %>%
    mutate( Class.II="FLA", Class.I="FLA" ),
  cc %>% select(Index,Compounds,Class.II, Class.I)  %>%
    rename("X1"="Compounds", "geneID"="Index" ) 
))

c<-left_join( FLA_all2 , ck[["wide.res"]], by=c("geneID"="gene") ) %>%
  group_by(cluster,Class.I) %>% add_count() %>%
  group_by(cluster,Class.II) %>% add_count




###
###Step4: dividing tpm-meta cluster into tpm and meta with specific gene name mark changes 
###

FA_metaId<-read_tsv("coreFA.txt",col_names = F)

meta_id_change<-rbind(
meta[which(meta$Index%in%FA_metaId$X1),] %>%
  select(Index,Compounds),
FLA_id)

tpm_id_change<- FLA_all2 %>%
  group_by(X1,Class.II) %>%
  mutate(n= seq(1, length.out=n() ) ) %>% ungroup() %>%
  mutate(Compounds= paste(Class.II, gsub( "_.*","", X1 ),n, sep="_")  ) %>%
  select(geneID, Compounds) %>%
  rename("Index"="geneID") %>%
  filter(Index%in%grep("^Ovio" ,Index, value = T ))

id_change<-rbind(meta_id_change, tpm_id_change) #%>% mutate_at(.vars = vars(Compounds), .funs = function(x){gsub("\\*","", gsub(" ","_",x))} )
rownames(id_change)<-1:length(id_change$Index)

ck<-ckk

#t<-ck$wide.res %>% group_by(gene) %>% add_count %>% filter(n>1)
  
#meta
#FLA_all2



###
ck_meta_ext<-ck
ck_tpm_ext<-ck

ck_meta_ext$wide.res<- ck_meta_ext$wide.res %>%
  filter(gene%in%base::grep("Ovio", gene, value = T, invert = T)) 
ck_meta_ext$long.res<- ck_meta_ext$long.res %>%
  filter(gene%in%base::grep("Ovio", gene, value = T, invert = T)) %>%
  filter(gene%in%base::grep("Ovio", gene, value = T, invert = T)) %>%
  separate(col=cluster_name,into=c("c","n","num")) %>%
  group_by(n, cell_type) %>%
  add_count() %>%
  ungroup() %>%
  mutate(cluster_name=paste(c," ", n, " (" ,nn, ")", sep="") )


ck_meta_ext$long.res<-left_join(ck_meta_ext$long.res, id_change, by=c("gene"= "Index") ) %>%
  mutate(new=ifelse(is.na(Compounds), gene, Compounds)) %>%
  rename( "old"="gene", "gene"="new" ) %>%
  select(cluster,gene,cell_type ,norm_value ,cluster_name )

ck_meta_ext$wide.res<-left_join(ck_meta_ext$wide.res, id_change, by=c("gene"= "Index") )%>%
  mutate(new=ifelse(is.na(Compounds), gene, Compounds)) %>%
  rename( "old"="gene", "gene"="new" ) %>%
  select(T1,T2,T3,T4,gene,cluster)

#rownames(ck_meta_ext$wide.res)<-ck_meta_ext$wide.res$gene



FLA_all2


ck_tpm_ext$wide.res<- ck_tpm_ext$wide.res %>%
  filter(gene%in%base::grep("Ovio", gene, value = T, invert = F)) 
ck_tpm_ext$long.res<- 
  ck_tpm_ext$long.res %>%
  filter(gene%in%base::grep("Ovio", gene, value = T, invert = F)) %>%
  separate(col=cluster_name,into=c("c","n","num")) %>%
    group_by(n, cell_type) %>%
    add_count() %>%
    ungroup() %>%
    mutate(cluster_name=paste(c," ", n, " (" ,nn, ")", sep="") )
  


ck_tpm_ext$long.res<-left_join(ck_tpm_ext$long.res, id_change, by=c("gene"= "Index") ) %>%
  mutate(new=ifelse(is.na(Compounds), gene, Compounds)) %>%
  rename( "old"="gene", "gene"="new" ) %>%
  select(cluster,gene,cell_type ,norm_value ,cluster_name )

ck_tpm_ext$wide.res<-left_join(ck_tpm_ext$wide.res, id_change, by=c("gene"= "Index") )%>%
  mutate(new=ifelse(is.na(Compounds), gene, Compounds)) %>%
  rename( "old"="gene", "gene"="new" ) %>%
  select(T1,T2,T3,T4,gene,cluster)

pdf('tpm_line44.pdf',height = 4,width = 8,onefile = F)
visCluster(object = ck_tpm_ext,
           plot.type = "line") 
dev.off()
pdf('meta_line44.pdf',height = 4,width = 8,onefile = F)
visCluster(object = ck_meta_ext,
           plot.type = "line")
dev.off()


markGenes = meta_id_change$Compounds

pdf('meta4.pdf',height = 10,width = 10,onefile = F)
visCluster(object = ck_meta_ext,
           plot.type = "both",column_names_rot = 45,line.side = "right",markGenes = markGenes,
           markGenes.side = "right"
           )
dev.off()


#markGenes=tpm_id_change$Compounds
#c("GPAT", "LACS", "MYB", "HLH")

markGenes<-tpm_id_change[grep("GPAT|LACS|MYB|HLH",tpm_id_change$Compounds, invert = T  ),]$Compounds


pdf('tpm4.pdf',height = 15,width = 15,onefile = F)
visCluster(object = ck_tpm_ext,
           plot.type = "both",column_names_rot = 45,line.side = "right",markGenes=markGenes,
           markGenes.side = "right")
dev.off()


saveRDS(ck_tpm_ext,file = "ck_tpm_ext.RDS")
saveRDS(ck_meta_ext, file="ck_meta_ext.RDS" )
saveRDS(ck, file="ck.RDS")



