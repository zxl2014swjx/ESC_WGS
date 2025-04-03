setwd("E:/项目-食管癌基因组不稳定性的机制研究_WGS")

dim(all_cnv)
dim(hot)
dim(hot_cnv_gene)


cnv<-read.table("cytoband_anno_organised.txt",header=TRUE,sep="\t")

sample<-as.data.frame(table(cnv$sample))
library(ggplot2)
lable<-factor(c(1:22))

lis<-list()

for(i in 1:dim(sample)[1]){
data<-cnv[which(cnv$sample==as.vector(sample[i,1])),]
num<-c(1:dim(data)[1])
data<-cbind(data,factor(num))
#pdf(paste0(as.vector(sample[i,1])," cna position.pdf"))
lis[[i]]<-ggplot(data,aes(x=chr,y=log2_copy_ratio,
          fill=factor(num),col=type))+
          geom_bar(stat="identity",position="dodge")+
          guides(fill=FALSE)+
          scale_fill_manual(values=rep("grey",length(num)))+
          scale_x_discrete(limits=lable)+
          ggtitle(paste0(as.vector(sample[i,1])," cna position"))+
          theme(plot.title = element_text(hjust = 0.5))
#print(s)
#dev.off()
}


library(magrittr)
library(ggpubr)
res <- ggarrange(plotlist=lis,ncol = 1, nrow = 15)
ggsave("bestPlot.pdf",res,width = 12,height = 16)


