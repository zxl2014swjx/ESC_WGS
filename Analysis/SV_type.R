setwd("E:/项目-食管癌基因组不稳定性的机制研究_WGS/organised_data")

dir()

sv_indel<-read.table("sv_indel_gene_anno.txt",header=TRUE,sep="\t")
sv_transl<-read.table("sv_transl_gene_anno.txt",header=TRUE,sep="\t")

head(sv_indel)
head(sv_transl)

setwd("E:/项目-食管癌基因组不稳定性的机制研究_WGS/项目-食管癌基因组不稳定性的机制研究_WGS_09.26/SV")

temp<-dir()
sample<-unlist(strsplit(temp,"[_]"))
sample_name<-matrix(sample,byrow=TRUE,ncol=2)
sample<-as.vector(sample_name[,1])
H<-data.frame()
for(i in 1:length(temp)){
a<-read.delim(temp[i],header=F,sep=" ")
b<-cbind(sample[i],a)
colnames(b)[1]<-"sample"
H<-rbind(H,b)
}

SV<-H[,-c(4:6)]
colnames(SV)<-c("sample","type","mechine","chr","start","end")
write.csv(SV,"SV.csv",row.names=FALSE)

a<-as.vector(SV[,2])
a[grep(paste0("^","tandem_dup"),a,perl=TRUE)]<-"tandem_dup"
a[grep(paste0("^","inver*"),a,perl=TRUE)]<-"invers"
a[grep(paste0("^","del*"),a,perl=TRUE)]<-"del"
a[grep(paste0("^","ins+"),a,perl=TRUE)]<-"ins"
a[grep(paste0("^","trans*"),a,perl=TRUE)]<-"transl"
SV[,2]<-a


library(ggplot2)
ggplot(SV,aes(x=sample,fill=type))+
geom_bar(stat="count",position="fill")+
theme(axis.text.x=element_text(size=10))+Sv
ggtitle("SV_type")+
theme(plot.title = element_text(hjust = 0.5))+
coord_flip()


