setwd("E:/项目-食管癌基因组不稳定性的机制研究_WGS")

cnv<-read.table("cytoband_anno_organised.txt",header=TRUE,sep="\t")
cnv_ratio<-cnv[which(abs(cnv$log2_copy_ratio)>0.04),]
head(cnv_ratio)


gain<-cnv_ratio[which(cnv_ratio$type=="gain"),]
chr<-as.data.frame(table(gain$cytoband_id))
all_late<-length(grep("EL+",as.vector(gain[,1]),perl=TRUE))
all_early<-dim(gain)[1]-all_late

H<-data.frame()
for(i in 1:dim(chr)[1]){
a<-gain[which(gain$cytoband_id==as.vector(chr[i,1])),]
b<-as.vector(a[,1])
late<-length(grep("EL+",b,perl=TRUE))
early<-length(b)-late
not_late<-all_late-late
not_early<-all_early-early
OR1<-(early*not_late)/(late*not_early)

x1 <- matrix(c(early,late,not_early,not_late),ncol = 2,byrow=TRUE)
p1<-chisq.test(x)$p.value
q1<-p.adjust(p,"fdr")

OR2<-(late*not_early)/(early*not_late)

c<-cbind(as.vector(chr[i,1]),early,late,not_early,not_late,OR,p,q)
colnames(c)[1]<-"chr"
H<-rbind(H,c)
}

write.csv(H,"diff_hot_chr.csv",row.names=FALSE)
diff<-H[which(as.numeric(as.vector(as.matrix(H$p)))<0.05),]



