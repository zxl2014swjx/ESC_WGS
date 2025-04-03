CNV<-read.csv("E:/merge_cytoband_anno.csv")

H<-data.frame()
for(i in 1:dim(CNV)[1]){
a<-as.vector(as.matrix(CNV[i,]))
a<-a[which(a!="NA")]
sample<-a[1]
b<-a[-1]
c<-cbind(sample,b)
H<-rbind(H,c)
}


ID<-unlist(strsplit(as.vector(H[,2]),"_"))

chr<-matrix(ID,byrow=TRUE,ncol=2)
data<-cbind(as.vector(H[,1]),chr)
colnames(data)<-c("sample","chr","type")
write.csv(data,"E:/data.csv",row.names=FALSE)

data<-data[!duplicated(data),]
data<-data[which(data[,3]!=0),]

H<-data.frame()
a<-as.data.frame(table(data[,1]))
for(i in 1:dim(a)[1]){
b<-data[which(data[,1]==a[i,1]),]
c<-as.data.frame(table(b[,2]))
d<-cbind(as.vector(b[1,1]),c)
H<-rbind(H,d)
}

fa<-H[which(H[,3]!=1),]
colnames(fa)<-c("sample","chr","freq")
false<-merge(data,fa,by=c("sample","chr"))
write.csv(false,"E:/false.csv",row.names=FALSE)