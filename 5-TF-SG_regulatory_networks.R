

setwd("F:/工作/3.诸葛菜/1.Plant_diversity/4.network")

library(tidyverse)

ck<-readRDS("ck.RDS")
ck_meta_ext<-readRDS("ck_meta_ext.RDS")
ck_tpm_ext<-readRDS("ck_tpm_ext.RDS")

meta<-read.csv("ALL_sample_data.csv", header = T)

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
  select(T1,T2,T3,T4) %>% 
  filter_all(any_vars(.>3))

###



FA<-read_tsv("Ov_FA.id")
FLA<-read_tsv("Ov_FLA.id")
meta<-read.csv("ALL_sample_data.csv", header = T)
cc<-meta
cc<- cc%>%filter( Class.I=="Flavonoids"|Class.I=='Lipids' )


FLA_id<-read_tsv("flavor.txt", col_names = F) %>%
  left_join(meta,by=c("X1"="cpd_ID")) %>%
  select(Index, Compounds) %>% 
  na.omit()


FLA_all2<-rbind( FA %>% select(geneID, X1) %>% mutate(Class.II="FA",Class.I="FA")  ,
                 rbind(
                   FLA %>%
                     mutate(X1=paste(V1,"_",V2,sep="")) %>%
                     select(geneID,X1) %>%
                     mutate( Class.II="FLA", Class.I="FLA" ),
                   cc %>% select(Index,Compounds,Class.II, Class.I)  %>%
                     rename("X1"="Compounds", "geneID"="Index" ) 
                 ))

SG_net<-FLA_all2[FLA_all2$geneID %>% grep("^Ovio",.),]# %>% filter(Class.II=="FA")
META_all2 <-FLA_all2 %>% filter(geneID%in%grep("Ovio",.$geneID,value = T,invert = T)) 
META_all2$geneID
out<-list()
for (i in 1:3){
  a<-(ck$wide.res%>%filter(cluster==i, gene%in%SG_net$geneID  ))[,1:4]
  aa<-(ck$wide.res%>%filter(cluster==i, gene%in%META_all2$geneID ))[,1:4]
  aaa<-as.data.frame(cor(t(aa),t(a))) %>%
    mutate(meta=rownames(.)) %>%
    gather(key=geneID, value=PCC, 1:(ncol(.)-1) ) %>%
    mutate(cluster=i )
  #aaa[,"gene"]<-rownames(aaa) 
  #aaa <- 
  out[[paste("cluster",i,sep = "_")]]<-aaa #as_tibble(aaa)
}

out_dat<-do.call("rbind", lapply(out, function(x){ x%>% left_join(FLA_all2,by=c("geneID")) %>% left_join(FLA_all2,by=c("meta"="geneID"))  } ))
out_dat<-out_dat %>% left_join(cc %>% select(`物质`,  Compounds), by=c("X1.y"="Compounds") )
PC_ind<-(out_dat %>% filter(Class.II.y=="Glycerol ester"))$meta %>% unique(.)

ck$wide.res[PC_ind,]

out_dat[which(out_dat$meta %in% rownames(ck$wide.res[PC_ind,] %>% filter(cluster==1))),]
cc<-out_dat[which(out_dat$meta %in% rownames(ck$wide.res[PC_ind,] %>% filter(cluster==3))),]



fimo<-read_tsv("TFBS_prediction/fimo.txt", col_names = F) %>%
  mutate_at(vars(X3),function(x){gsub("::.*","",x)})

TF<-read.table("TFBS_prediction/TF_and_best1_in_Ath.list.final.ov",header=F,sep = "\t" , quote = "")

TFBS<-fimo %>% select(X3,X2) %>% unique()
TF$V2

out2<-list()

#non-pvalue person
#for (i in 1:3){
#  a<-(ck$wide.res%>%filter(cluster==i, gene%in%SG_net$geneID  ))[,1:4]
#  aa<-(ck$wide.res%>%filter(cluster==i, gene%in%TF$V2 ))[,1:4]
#    aaa<-as.data.frame(cor(t(aa),t(a), method="pearson")) %>%
#    mutate(meta=rownames(.)) %>%
#    gather(key=geneID, value=PCC, 1:(ncol(.)-1) ) %>%
#    mutate(cluster=i )
#  out2[[paste("cluster",i,sep = "_")]]<-aaa #as_tibble(aaa)
#}

#pearson with pvalue
for (i in 1:3){
  a<-(ck$wide.res%>%filter(cluster==i, gene%in%SG_net$geneID  ))[,1:4]
  aa<-(ck$wide.res%>%filter(cluster==i, gene%in%TF$V2 ))[,1:4]
  out2[[paste("cluster",i,sep = "_")]]<-data.frame( meta=rep(0, nrow(a)*nrow(aa)), geneID= rep(0, nrow(a)*nrow(aa)), PCC=rep(0, nrow(a)*nrow(aa)),pvalue=rep(0, nrow(a)*nrow(aa)) , cluster=rep(0, nrow(a)*nrow(aa)) )
  for (j in 1:nrow(a)){
    xx<-a[j,]
    for (g in 1:nrow(aa)){
      yy<-aa[g,]
      corr<-cor.test( as.numeric(xx), as.numeric(yy), method = "pearson" )
      
      out2[[paste("cluster",i,sep = "_")]][(j-1)*g+g,]<-c(rownames(yy), rownames(xx),   corr$estimate, corr$p.value,i )
    }
  }
  #aaa<-as.data.frame(cor(t(aa),t(a), method="pearson")) %>%
  # mutate(meta=rownames(.)) %>%
  #gather(key=geneID, value=PCC, 1:(ncol(.)-1) ) %>%
  #mutate(cluster=i )
  #out2[[paste("cluster",i,sep = "_")]]<-aaa #as_tibble(aaa)
}



out_dat2<-do.call("rbind", lapply(out2, function(x){ x%>% left_join(FLA_all2,by=c("geneID")) %>% left_join(TF,by=c("meta"="V2"))  } ))

out_dat2<- out_dat2 %>% left_join(TFBS,by=c("geneID"="X3")) %>%
  filter(V3==X2) %>%
  #filter(V3==X2, PCC>0.95|PCC< -0.95 ) %>%
  mutate( TF = ifelse(geneID==meta, X1,"") ) %>%
  group_by(meta) %>%
  mutate(TF_c= paste(TF,collapse = "")  )  %>%
  ungroup() %>%
  #arrange(desc(TF_c)) %>%
  filter( meta!="geneID")




out_dat2_FA<-out_dat2 %>% filter( Class.I=="FA" , cluster==3) %>% filter(pvalue<=0.05)

out_dat2_FLA<-out_dat2 %>% filter( Class.I=="FLA", cluster==3 )%>% filter(pvalue<=0.05)


GetNode<-function(x){
  rbind(
    x%>% mutate(X4=X1)  %>%  select(geneID, X1,Class.I,X4) %>% unique(.) %>% rename(X1="geneID", X2="X1", X3="Class.I" ),
    x%>%  group_by(cluster, V3) %>% mutate(n=seq(1,length.out = n()), X4=ifelse(n==1,V3,""))  %>% ungroup() %>%  select(meta, V1,V3,X4) %>% unique(.) %>% rename(X1="meta", X2="V1", X3="V3")
  ) #%>%
  #group_by(X3) %>%
  #mutate(n=seq(1,length.out = n())) %>%
  #mutate(X4=ifelse(X3=="FA", X2, ifelse(n==1, X3, "")  )  )
}


###file for network 

out_dat2_FA


write_tsv(GetNode(out_dat2_FA), file ="FA_C1_node.txt"  , col_names = T)
write_tsv(GetNode(out_dat2_FLA),file ="FLA_C1_node.txt", col_names = T)
write_tsv(out_dat2_FA, file="FA_C1_edge.txt", col_names = T)
write_tsv(out_dat2_FLA, file="FLA_C1_edge.txt", col_names = T)


#cor pheat map for ANR_ban
all_tpm[(out_dat2_FLA %>% filter(X1=="ANR_ban") %>% select(geneID) %>% unique(.))$geneID,]
#pheatmap( all_tpm[(out_dat2_FLA %>% filter(X1=="ANR_ban") %>% select(geneID) %>% unique(.))$geneID,]  )
#cor pheat map for SAD
all_tpm[(out_dat2_FA %>% filter(X1=="SAD") %>% select(geneID) %>% unique(.))$geneID,]


library(PerformanceAnalytics)
#install.packages("PerformanceAnalytics")
chart.Correlation(t(all_tpm[(out_dat2_FA %>% filter(X1=="SAD") %>% select(geneID) %>% unique(.))$geneID,]),  histogram=TRUE)
chart.Correlation( t(all_tpm[(out_dat2_FLA %>% filter(X1=="ANR_ban") %>% select(geneID) %>% unique(.))$geneID,]),  histogram=TRUE  )


pheatmap( all_tpm[c("Ovio46619", (out_dat2_FA %>% filter( geneID=="Ovio46619") %>% filter( X2=="bZIP"|X2=="bHLH"|X2=="NAC"|X2=="HD-ZIP"|X2=="G2-like" ) %>% arrange(X2) )$meta),], 
          cluster_cols = F ,cluster_rows = F, scale = "row" )
out_dat2_FA %>% filter( geneID=="Ovio46619") %>% filter( X2=="bZIP"|X2=="bHLH"|X2=="NAC"|X2=="HD-ZIP"|X2=="G2-like" ) %>% arrange(X2)
#pheatmap( all_tpm[c("Ovio43049", (out_dat2_FLA %>% filter( geneID=="Ovio43049") %>% filter(X2=="B3"|X2=="bHLH"|X2=="MYB_related"|X2=="HD-ZIP"|X2=="MYB") %>% arrange(X2) )$meta),], 
#          cluster_cols = F ,cluster_rows = F, scale = "row" )
out_dat2_FLA %>% filter( geneID=="Ovio43049") %>% filter(X2=="B3"|X2=="bHLH"|X2=="MYB_related"|X2=="HD-ZIP"|X2=="MYB") %>% arrange(desc(TF))

pheatmap(all_tpm[c("Ovio43049", (head(out_dat2_FLA %>% filter( geneID=="Ovio43049") %>% filter(X2=="B3"|X2=="bHLH"|X2=="MYB_related"|X2=="HD-ZIP"|X2=="MYB") %>% arrange(desc(PCC)),n = 20 ) %>% arrange(X2))$meta ),],
         cluster_cols = F ,cluster_rows = F, scale = "row")
head(out_dat2_FLA %>% filter( geneID=="Ovio43049") %>% filter(X2=="B3"|X2=="bHLH"|X2=="MYB_related"|X2=="HD-ZIP"|X2=="MYB") %>% arrange(desc(PCC)),n = 20 ) %>% arrange(X2)
#WGD analysis

write_csv( out_dat2_FA %>% filter( geneID=="Ovio46619") %>%
             arrange(X2) %>% select(meta,geneID,PCC,pvalue,cluster,X1,V1,V3,V4,X2) %>% left_join(all_tpm%>%mutate(meta=rownames(.)), by=c("meta")), 
           file = "SAD_sup.csv" )
write_csv(out_dat2_FLA %>% filter( geneID=="Ovio43049") %>% 
            arrange(desc(PCC)) %>% select(meta,geneID,PCC,pvalue,cluster,X1,V1,V3,V4,X2)%>% left_join(all_tpm%>%mutate(meta=rownames(.)), by=c("meta")) , file="ANR_sup.csv")

write_csv(out_dat2 %>% arrange(desc(PCC)) %>% select(meta,geneID,PCC,pvalue,cluster,X1,V1,V3,V4,X2), file="all_cor_sup.csv"  )
write_csv(FLA_all2, file = "FLA_FA_gene_meta.csv" )
write_csv(ck$wide.res, file = "kmeans_cluster_classification.csv")
write_csv(all_tpm, file="gene_tpm_used.csv")

length(all_tpm$T1)


dup<-read_tsv("dup_gene.tsv") %>% filter(class=="wgd")
FLA_all3<- FLA_all2 %>% 
  filter( geneID%in%grep("Ovio",.$geneID, value = T) ) %>%
  select(X1,geneID) %>%
  group_by(X1) %>%
  mutate( c= paste("copy", seq(1,length.out = n()) , sep="_") ) %>%
  spread( key=c, value = geneID) %>% ungroup() %>%
  mutate(WGD1="none", WGD2="none")


for(i in 1:nrow(FLA_all3)){
  b<-c("none","none")
  for (j in 1:nrow(dup)){
    if(length(which(as.character(FLA_all3[i,][,2:8]) %in% as.character(dup[j,2:3]))) == 2){
      b<- as.character(FLA_all3[i,][,2:8])[which(as.character(FLA_all3[i,][,2:8]) %in% as.character(dup[j,2:3]) )]
    }
  }
  FLA_all3[i,9]<-b[1]
  FLA_all3[i,10]<-b[2]
  
}

FLA_all3_WGD<-FLA_all3 %>% filter(WGD1!="none") %>% mutate(pcc=0,pvalue=0) %>% mutate(WGD1_C=0, WGD2_C=0) %>%  filter(X1!="MYB5_MYB5")

for(i in 1:nrow(FLA_all3_WGD)){
  cco<-cor.test( as.numeric(all_tpm[as.character(FLA_all3_WGD[i,"WGD1"]),]) , as.numeric(all_tpm[as.character(FLA_all3_WGD[i,"WGD2"]),]), method = "pearson")
  FLA_all3_WGD[i,"pcc"]<-cco$estimate
  FLA_all3_WGD[i,"pvalue"]<-cco$p.value
  FLA_all3_WGD[i,"WGD2_C"] <-ck$wide.res[as.character(FLA_all3_WGD[i,"WGD2"]),"cluster"]
  FLA_all3_WGD[i,"WGD1_C"] <-ck$wide.res[as.character(FLA_all3_WGD[i,"WGD1"]),"cluster"]
}
#corr$estimate, corr$p.value

ggplot(FLA_all3_WGD, aes(y=pcc))+
  # geom_histogram() +
  geom_density( fill="#7fd6dd", color="#7fd6dd" ) +
  geom_hline(yintercept = 0.6, linetype="dashed" ) +
  theme_bw() +
  #scale_color_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD" ))+
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black"),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",hjust=0.5))  +
  ylab("Pearson correlation between WGD duplicates\n among SOG-related structrual genes") +
  xlab("Relative density") +
  annotate("text", y=0.45,x=1.3, label="PCC=0.6" ,size =5 ) 

ggsave("SOG_related_WGD_pearson.pdf")

all_tpm_all<-as.data.frame(read_tsv("Ovio_tpm.txt"))
rownames(all_tpm_all)<-all_tpm_all$geneid
all_tpm_all<-all_tpm_all[,-1]
all_tpm_all<-all_tpm_all%>%select(grep("P0",colnames(.), value = T, invert = T)) %>%
  mutate(F=(OV_F.P1+OV_F.P2+OV_F.P3)/3, L=(OV_L.P1+OV_L.P2+OV_L.P3)/3, R=(OV_R.P1+OV_R.P2+OV_R.P3)/3, S=(OV_S.P1+OV_S.P2+OV_S.P3)/3,
        T1=(OV_T1A.1+OV_T1A.2+OV_T1A.3+OV_T1A.4)/4 ,T2=(OV_T2A.1+OV_T2A.2+OV_T2A.3+OV_T2A.4)/4 , 
        T3=(OV_T3A.1+OV_T3A.2+OV_T3A.3+OV_T3A.4)/4, T4=(OV_T4A.1+OV_T4A.2+OV_T4A.3+OV_T4A.4)/4  ) %>%
  select(F,L,R,S,T1,T2,T3,T4) 
  #filter_all(any_vars(.>3))
  
  ###
  
#library(pheatmap)

pdf("DGAT_all_tis_cor_pheat.pdf")
pheatmap( cor( t(all_tpm_all[FLA_all2[ which( FLA_all2$X1 %in% c("DGAT")),]$geneID,]) ), display_numbers = TRUE )
dev.off()
color = colorRampPalette(c("#6D9EC1", "white", "#E46726") )(50)
pdf("DGAT_all_tis_exp_pheat.pdf")
pheatmap(log2(t(all_tpm_all[FLA_all2[ which( FLA_all2$X1 %in% c("DGAT")),]$geneID,])+0.01), scale = "column", cluster_cols = T, cluster_rows = F  )
dev.off()
color = colorRampPalette(c("#495086", "white", "#9a2722"))(50)
pdf("DGAT_seed_cor_pheat.pdf")
pheatmap( cor( t(all_tpm[FLA_all2[ which( FLA_all2$X1 %in% c("DGAT")),]$geneID,]) ) )
dev.off()




library(ggcorrplot)
round(cor( t(all_tpm_all[FLA_all2[ which( FLA_all2$X1 %in% c("DGAT")),]$geneID,])),3 )
ggcorrplot::ggcorrplot( round(cor( t(all_tpm_all[FLA_all2[ which( FLA_all2$X1 %in% c("DGAT")),]$geneID,])),3 ),hc.order=T,outline.color="white",
                        lab = TRUE,    colors = c("#6D9EC1", "white", "#E46726")     )
ggsave("DGAT_all_tis_cor_pheat.pdf")
ggcorrplot::ggcorrplot( round(cor( t(all_tpm[FLA_all2[ which( FLA_all2$X1 %in% c("DGAT")),]$geneID,]) ),3 ),hc.order=T,outline.color="white",
                        lab = TRUE,    colors = c("#6D9EC1", "white", "#E46726")     )

ggsave("DGAT_seed_cor_pheat.pdf")


write.table(all_tpm_all[FLA_all2[ which( FLA_all2$X1 %in% c("DGAT")),]$geneID,], file ="F:/工作/3.诸葛菜/1.Plant_diversity/4.network/DGAT/DGAT_tpm.tsv", row.names = TRUE, quote=F  )



FLA_all3_WGD %>%
  mutate(cor_class=ifelse(pcc>0.6,"high-cor","low-cor") ) %>%
  ggplot(aes(x=cor_class))+
  geom_bar() 

ck$wide.res %>%
  filter( gene%in%(grep("Ovio",.$gene, value = T, invert = T)) ) %>%
  left_join(meta%>%select(Index, Class.I, Class.II, kegg_map) ,by=c("gene"="Index") ) %>%
  group_by(cluster, Class.I) %>%
  add_count() %>%
  select(Class.I, n) %>%
  unique(.) %>%
  ggplot(aes(y=as.character(cluster),x=n,fill=factor(Class.I)))+
  geom_col() +
  theme_classic() +
  scale_fill_manual(values = c('Alkaloids'='#a663cc',  "Amino acids and derivatives"= '#023e8a',
                               "Flavonoids"='#c72c2c',
                               "Lignans and Coumarins"='#fdc500',
                               "Lipids"='#7FFF7F',
                               "Nucleotides and derivatives"= '#ff7b00',
                               "Organic acids"='#ce4257',
                               "Phenolic acids"= '#a2d6f9',
                               "Tannins"='#b08968',
                               "Terpenoids"='#227c9d', 
                               "Others"='#fdf0d5')) +
  xlab("") +
  xlim(0,500)



ck$wide.res %>%
  filter( gene%in%(grep("Ovio",.$gene, value = T, invert = T)) ) %>%
  left_join(meta%>%select(Index, Class.I, Class.II, kegg_map) ,by=c("gene"="Index") ) %>%
  mutate(Class=ifelse(Class.II=="Free fatty acids", Class.II, "Others")) %>% ungroup()%>%
  filter(Class.I=="Lipids")%>%
  group_by(cluster, Class) %>%
  add_count() %>%
  select(Class, n,cluster) %>%
  unique(.) %>%
  #group_by(Class) %>%
  ungroup()%>%
  mutate(nn=sum(n))%>%
  mutate(label=paste(round(n/nn,2)*100, "%", sep="" )) %>%
  ggplot(aes(y=as.character(cluster),x=n,fill=factor(Class,levels=rev(c("Free fatty acids", "Others")))))+
  geom_col() +
  geom_text( aes(y=as.character(cluster),x=n, label=label  ), color="white", size=5 )+
  theme_classic() +
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5,  vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black"),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",hjust=0.5)) 


  

FLA_all2 %>%filter(geneID%in%grep("Ovio",.$geneID, value = T, invert = T )) %>% 
  left_join(ck$wide.res, by=c("geneID"="gene")) %>%
  group_by(Class.I,cluster) %>%
  add_count() %>%
  select(Class.I, cluster,n) %>% unique(.) %>%
  group_by(Class.I) %>% mutate(sum=sum(n),prop=n/sum(n), percen=paste(round(prop,2)*100 , "%", sep="" )) %>%
  
  filter(Class.I=="Flavonoids") %>%
  
  ggplot(aes(y=as.character(cluster),x=prop,fill=Class.I))+
  geom_col(position = "dodge")+
  geom_text(aes(y=as.character(cluster),x=prop/2,label=paste(n,"(",percen, ")",sep="")), color="white", size=5)+
  #facet_wrap(~Class.I,nrow = 2)+
  theme_classic() +
  #scale_color_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD" ))+
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, face="italic", vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black"),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",face="italic",hjust=0.5)) 


ggsave("cluster_distribution.pdf")
  
