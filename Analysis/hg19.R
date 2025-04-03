source("http://bioconductor.org/biocLite.R")

library(BiocManager)
BiocManager::install("StructuralVariantAnnotation")
library(rtracklayer)
library(StructuralVariantAnnotation)
help(package="StructuralVariantAnnotation")

library(SummarizedExperiment)
library(DelayedArray)
library(matrixStats)
library(BiocParallel)
library(DelayedArray)
library(Rsamtools)
library(Biostrings)
library(XVector)
library(Biostrings)
library(VariantAnnotation)


setwd("E:/zxl/11_WGS(9.16-)/WGS/data/SV_result/")
vcf<-readVcf('E11TD_E11ND/somatic.filter.final.variants.vcf')
head(rowData(vcf))
head(fixed(vcf),3)
alternate <- alt(vcf)
geno(vcf)


hg19<-read.table("E:/hg19.txt",header=T,sep="\t")
gene<-as.data.frame(table(hg19$Gene))
g1<-gene[gene[,2]==1,]
g2<-gene[gene[,2]!=1,]
colnames(g1)[1]<-"Gene"
colnames(g2)[1]<-"Gene"
g1_data<-merge(g1,hg19,by="Gene")
g2_data<-merge(g2,hg19,by="Gene")
g3<-g2[-c(1:11),]
g3_data<-g2_data[-c(1:36),]
H2<-data.frame()
for(i in 1:dim(g3)[1]){
a<-g3_data[g3_data$Gene==g3[i,1],]
start<-min(c(as.vector(a$Start),as.vector(a$Stop)))
stop<-max(c(as.vector(a$Start),as.vector(a$Stop)))
b<-cbind(a[1,1:3],start,stop)
H2<-rbind(H2,b)
}



shuju<-H[,-7]
hg19<-read.table("E:/hg19.txt",header=T,sep="\t")
for(i in 1:dim(shuju)[1]){
a<-hg19[hg19$Chr==shuju[1,2],]
begin<-as.numeric(as.vector(shuju[1,3]))
s1<-a[which(a[,3]>begin & a[,3]<end),]
end<-as.numeric(as.vector(shuju[1,4]))
s2<-a[a[,4]<end,]


}

