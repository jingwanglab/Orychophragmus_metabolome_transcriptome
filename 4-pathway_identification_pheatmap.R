
####
####STEP1:  filter AT FA pathway gene id 
####


setwd("F:/工作/3.诸葛菜/1.Plant_diversity/3.FA_pathway")

library(tidyverse)
#install.packages("tidyverse")

ara_lip<-read.csv("aralip_data.csv", header = F)
FA_id<-read_tsv("FA_geneid.txt", col_names = F)

join<-left_join(FA_id,ara_lip, by=c("X1"="V3")) %>%
  filter(V8!="Seq. Similarity", V7!="")

#manual advise
#write.csv(join, file="1.FA_id_filter_first.csv", row.names = FALSE)

ara_filter<-read.csv("1.FA_id_filter_first.csv") %>%
  mutate_at(vars(V5), ~ toupper(.))



####
####
#### STEP2: read syntenic realtion ship between Bra species and then fitler out Ov FA gene id 
#### 
####




setwd("F:/工作/3.诸葛菜/00.diploid/14.genefamily/00.analysis/ortho")
all_tree<-readRDS("allTrees.RDS")
FAE<- c("Ovio08648", "Ovio10684", "Ovio26944")


OGtbl <- read_tsv("OGtbl.tsv",col_types = cols())

ScaffoldGene <- read.table("gene_on_scaffold.txt")

OGtbl<-OGtbl[which(!OGtbl$geneID %in% ScaffoldGene$V1),]


#how many N11s per N1
Ovtbl<-
  OGtbl %>% 
  select(OG, N1, N11, spc, geneID) %>%
  filter(!is.na(N11) & spc == "Ovio") 

#AtperN1<- 
AtperN1<-
  OGtbl %>%
  select(OG, N1, spc, geneID) %>%
  filter(spc == "Atha") %>%
  count(N1) 





AtperN0<-
  OGtbl %>%
  select(OG, N0, spc, geneID) %>%
  filter(spc == "Atha") %>%
  count(N0) 


OvperN1 <- 
  OGtbl %>% 
  select(OG, N1, N11, spc, geneID) %>%
  filter(!is.na(N11) & spc == "Ovio") %>%
  count(N1)

Ov_resolved_gene<-
  Ovtbl %>%
  left_join(AtperN1,by = "N1") %>%
  left_join(OvperN1, by = "N1") %>%
  na.omit() %>%
  left_join(  OGtbl %>%
                select(OG, N1, spc, geneID) %>%
                filter(spc == "Atha") , by = "N1") %>%
  mutate_at( vars(geneID.y), ~gsub("\\..*","",. ))


####
#### STEP3: get Atha and Ovio relationship across FA pathway
####

Ov_resolved_gene
ara_filter

ara_filter_orth<-left_join(ara_filter,Ov_resolved_gene,by=c("V5"= "geneID.y" ) )
ara_filter_orth_first<- ara_filter_orth %>%
  filter(!is.na( spc.y))

ara_filter_orth_un<- ara_filter_orth %>% 
  filter(is.na( spc.y) ) %>%
  left_join( ara_filter_orth %>% filter(!is.na(spc.y)), by=c("X1")  ) %>%
  filter(is.na(OG.y.y))

ara_filter_orth_un_sign<-left_join(ara_filter_orth_un[,1:11], OGtbl %>% mutate_at( vars(geneID), ~gsub("\\..*","",. )), by=c("V5.x"="geneID") ) %>%
  left_join(OGtbl,by=c("OG")) %>%
  filter(spc.y=="Ovio") 

###all result
ara_filter_orth_un_sign
ara_filter_orth_first

FA_id<-rbind(
  ara_filter_orth_un_sign %>% select(X1,geneID),
  ara_filter_orth_first%>% select(X1,geneID.x) %>% rename("geneID"="geneID.x")
)

write_tsv(FA_id,file="Ov_FA.id", col_names = T)


###
### STEP 3.2 : Flav pathway ID match 
###
setwd("F:/工作/3.诸葛菜/1.Plant_diversity/2.fla_pathway")

fla_id<-read.csv("flavor_gene.csv", header = F)


ara_filter_orth<-left_join(fla_id,Ov_resolved_gene,by=c("V3"= "geneID.y" ) )
ara_filter_orth_first<- ara_filter_orth %>%
  filter(!is.na( spc.y))


ara_filter_orth_un<- ara_filter_orth %>% 
  filter(is.na( spc.y) ) %>%
  left_join( ara_filter_orth %>% filter(!is.na(spc.y)), by=c("V2")  ) %>%
  filter(is.na(OG.y.y))


ara_filter_orth_un_sign<-left_join(ara_filter_orth_un[,1:3], OGtbl %>% mutate_at( vars(geneID), ~gsub("\\..*","",. )), by=c("V3.x"="geneID") ) %>%
  left_join(OGtbl,by=c("OG")) %>%
  filter(spc.y=="Ovio") 


###all result
ara_filter_orth_un_sign
ara_filter_orth_first

FLA_id<-rbind(
  ara_filter_orth_un_sign %>% select(V1.x,V2,geneID)  %>% rename("V1"="V1.x") ,
  ara_filter_orth_first%>% select(V1,V2,geneID.x) %>% rename("geneID"="geneID.x")
)

write_tsv(FLA_id,file="Ov_FLA.id", col_names = T)




setwd("F:/工作/3.诸葛菜/1.Plant_diversity")

library(tidyverse)
library(pheatmap)

all_tpm<-as.data.frame(read_tsv("0.data/Ovio_tpm.txt"))
rownames(all_tpm)<-all_tpm$geneid
all_tpm<-all_tpm[,-1]
all_tpm<-all_tpm%>%select(grep("P0",colnames(.), value = T, invert = T)) %>% 
  select( grep("OV_T",colnames(.), value = T)  )  %>%
  mutate(T1=(OV_T1A.1+OV_T1A.2+OV_T1A.3+OV_T1A.4)/4 ,T2=(OV_T2A.1+OV_T2A.2+OV_T2A.3+OV_T2A.4)/4 , 
         T3=(OV_T3A.1+OV_T3A.2+OV_T3A.3+OV_T3A.4)/4, T4=(OV_T4A.1+OV_T4A.2+OV_T4A.3+OV_T4A.4)/4  ) %>%
  select(T1,T2,T3,T4) #%>% filter_all(any_vars(.>3))




FA<-read_tsv("Ov_FA.id")
FLA<-read_tsv("Ov_FLA.id")



idlist<-rbind(FA, FLA %>% mutate( X1=paste(V1,V2,sep="_") ) %>% select(X1,geneID)) %>%
  group_by(X1) %>%
  mutate(n=seq(1, length.out = n() ))  %>%
  mutate( sym=paste(X1,n,sep = "_") )

plot_tpm<-all_tpm[idlist$geneID, ]

rownames(plot_tpm) <- idlist$sym
plot_tpm<-plot_tpm%>% na.omit()

pdf(file="pathway_pheat.pdf", width = 5, height = 12)
pheatmap(mat = plot_tpm, scale = "row", cluster_rows = F, cluster_cols = F  )
dev.off()



###all tissue 
all_tpm<-as.data.frame(read_tsv("0.data/Ovio_tpm.txt"))
rownames(all_tpm)<-all_tpm$geneid
all_tpm<-all_tpm[,-1]
all_tpm<-all_tpm%>%select(grep("P0",colnames(.), value = T, invert = T)) %>% 
  #select( grep("OV_T",colnames(.), value = T)  )  %>%
  mutate(flower=(OV_F.P1+OV_F.P2+OV_F.P3)/3,leaf=(OV_L.P1+OV_L.P2+OV_L.P3)/3,
         root=(OV_R.P1+OV_R.P2+OV_R.P3)/3,stem=(OV_S.P1+OV_S.P2+OV_S.P3)/3,
         T1=(OV_T1A.1+OV_T1A.2+OV_T1A.3+OV_T1A.4)/4 ,T2=(OV_T2A.1+OV_T2A.2+OV_T2A.3+OV_T2A.4)/4 , 
         T3=(OV_T3A.1+OV_T3A.2+OV_T3A.3+OV_T3A.4)/4, T4=(OV_T4A.1+OV_T4A.2+OV_T4A.3+OV_T4A.4)/4  ) %>%
  select(flower, leaf, root, stem,T1,T2,T3,T4) #%>% filter_all(any_vars(.>3))



plot_tpm<-all_tpm[idlist$geneID, ]

rownames(plot_tpm) <- idlist$sym
plot_tpm<-plot_tpm%>% na.omit()

pdf(file="pathway_pheat_all.pdf", width = 6, height = 12)
pheatmap(mat = plot_tpm, scale = "row", cluster_rows = F, cluster_cols = F  )
dev.off()
