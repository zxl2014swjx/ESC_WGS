setwd("E:/项目-食管癌基因组不稳定性的机制研究_WGS")

cnv<-read.table("cytoband_anno_organised.txt",header=TRUE,sep="\t")
cnv_ratio<-cnv[which(abs(cnv$log2_copy_ratio)>0.04),]
write.csv(cnv_ratio,"cnv_ratio.csv",row.names=FALSE)



row<-as.vector(as.data.frame(table(cnv_ratio$cytoband_id))[,1])
col<-as.vector(as.data.frame(table(cnv_ratio$sample))[,1])
pinpu<-matrix(rep(0,length(row)*length(col)),nrow=length(row))
rownames(pinpu)<-row
colnames(pinpu)<-col


for(i in 1:length(row)){
for(j in 1:length(col)){
a<-which(cnv_ratio$cytoband_id==rownames(pinpu)[i])
b<-which(cnv_ratio$sample==colnames(pinpu)[j])
c<-intersect(a,b)
d<-length(table(as.vector(cnv_ratio[c,]$type)))
if(d==1){
if(cnv_ratio[c,7]=="gain"){pinpu[i,j]=length(c)}
else{pinpu[i,j]=-length(c)}
}
}
}



sample <-as.vector(as.data.frame(table(cnv_ratio$sample))[,1])
lable<-as.data.frame(c(rep("early",7),rep("late",8)))
rownames(lable)<-sample
colnames(lable)<-"sample"


library(pheatmap)
pdf(file="cnv_position.pdf",height=10,width=10)
pheatmap(pinpu,display_numbers = F,
cluster_cols = FALSE,cluster_rows = TRUE,
annotation_col = lable,border_color='black',
col=colorRampPalette(c("#45FF45",'white',"#B95CF4"))(100),
show_rownames=F,fontsize_col=10,
cutree_row = 2,
gaps_col = 7
)
dev.off()
colorRampPalette(c("white","purple"))(256)[187]
colorRampPalette(c("white","green"))(256)[187]
pinpu[1:4,1:4]


#######################################

#hot_gain
gain<-pinpu
for(i in 1:dim(gain)[1]){
for(j in 1:dim(gain)[2]){
if(gain[i,j]>=1){}
else{gain[i,j]=0}
}
}
gain_chr<-apply(ga)



gain_chr<-gain_chr[order(gain_chr[,2],decreasing=TRUE),]
gain<-as.vector(gain_chr[1:20,1])
a<-as.data.frame(rep("gain",length(gain)))
rownames(a)<-gain
colnames(a)<-"type"

loss<-cnv_ratio[which(cnv_ratio$type=="loss"),]
loss_chr<-as.data.frame(table(loss$cytoband_id))
loss_chr<-loss_chr[order(loss_chr[,2],decreasing=TRUE),]
loss<-as.vector(loss_chr[1:20,1])
b<-as.data.frame(rep("loss",length(loss)))
rownames(b)<-loss
colnames(b)<-"type"
anno_row<-rbind(a,b)

length(unique(rownames(pinpu)))

b<-vector()
for(i in 1:dim(anno_row)[1]){
a<-grep(paste0("^",rownames(anno_row)[i],"$"),rownames(pinpu),perl=TRUE)
b<-c(b,a)
}
hot<-pinpu[b,]

pdf(file="hot_cnv_position.pdf",height=10,width=10)
pheatmap(hot,display_numbers = F,
cluster_cols = FALSE,cluster_rows = FALSE,
annotation_col = lable,
annotation_row = anno_row,
border_color=NA,
col=colorRampPalette(c("#45FF45",'white',"#B95CF4"))(100),
fontsize_col=10,fontsize_row=10,
cutree_row = 2,
gaps_col = 7,
gaps_row=20
)
dev.off()


