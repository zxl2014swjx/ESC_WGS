CNV<-read.delim("E:/E11ND_E11TD.seg.txt",header=TRUE,sep="\t")
CNV<-subset(CNV,select=c("SampleID","chrom","Corrected_Call"))
CNV<-CNV[!duplicated(CNV),]

row<-as.vector(unique(as.vector(CNV$chrom)))
col<-as.vector(unique(as.vector(CNV$SampleID)))
pinpu<-matrix(rep(0,length(row)*length(col)),nrow=length(row))
rownames(pinpu)<-row
colnames(pinpu)<-col

for(i in 1:length(row)){
for(j in 1:length(col)){
a<-which(CNV$chrom==rownames(pinpu)[i])
b<-which(CNV$SampleID==colnames(pinpu)[j])
c<-intersect(a,b)
if(length(c)!=0){
if(CNV[c,]$Corrected_Call =="GAIN"){
pinpu[i,j]=1
}
if(CNV[c,]$Corrected_Call =="LOSS"){
pinpu[i,j]=-1
}
}
}
}
library(pheatmap)
pheatmap(pinpu,
cluster_cols = TRUE,cluster_rows = FALSE,
legend=TRUE,fontsize_row=10,fontsize_col=10,
show_rownames=T,show_colnames=T,
col=colorRampPalette(c("blue","white","red"))(100))

