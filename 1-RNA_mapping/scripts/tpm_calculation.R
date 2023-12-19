library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

counts<-list()
counts[["in"]]<-as.data.frame(read.table(args[1], header = T ,sep="\t"))


fpkm_cal<-function(x){
  rownames(x)<-x$geneid
  x<-x[,-1]
  oneB<-10^9
  sum<-apply(x[1:(ncol(x)-1)],2,function(y){sum(y)})
  fpkm_self<-sapply(1:ncol(x[1:(ncol(x)-1)]),function(y){ (oneB*x[,y])/(sum[y]*as.double(x[,"gene_length"]))  })
  dimnames(fpkm_self)<-dimnames(x[,1:ncol(x)-1])
  return(as.data.frame(fpkm_self))
}

fpkm<-lapply(counts, fpkm_cal)

fpkmToTpm <- function(fpkm) {
exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

tpm<-lapply(fpkm, function(x){ apply(x,2, fpkmToTpm)} )

write.table(fpkm[["in"]], file=args[2], sep="\t", row.names=T, col.names=TRUE, quote= F   )
write.table(tpm[["in"]], file=args[3], sep="\t", row.names=T, col.names=TRUE, quote= F   )

