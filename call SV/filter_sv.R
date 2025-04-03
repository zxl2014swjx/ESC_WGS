

temp<-dir(pattern="*callsv.anno.filter.csv")
H1<-data.frame()
for(i in 1:length(temp)){
a<-read.csv(temp[i],header=T,sep=",")
sample<-unlist(strsplit(temp[i],"[.]"))[1]
sample<-unlist(strsplit(sample,"_"))[2]
b<-cbind(sample,a)
H1<-rbind(H1,b)
}

sv<-H1
data<-sv[which(as.vector(sv$AutoInterpretation)!="N"),]

data<-data[which(as.vector(data$Gene1)!="Unknown" & 
as.vector(data$Gene2)!="Unknown"),]

write.table(data,"sv_result.txt",row.names=F,quote=F,sep="\t")

need<-data[which(data$IsCase=="Y"),]
need<-data[which(data$Freq>0.01),]
need<-data[which(data$IsFunc=="Y"),]
need<-data[which(data$geneCross=="Y"),]
need<-data[which(data$NumOfexa!=0),]



sv_num<-as.data.frame(table(need[,1]))
colnames(sv_num)[1]<-"sample"
sv_num<-merge(clinical,sv_num,by="sample")
tapply(sv_num$Freq,sv_num$Type,mean,na.rm=TRUE)


data[which(as.vector(data$Gene1)!="Unknown" & 
as.vector(data$Gene2)!="Unknown"),]


need<-subset(data,select=c(1,4,18,20,22,24))
write.table(need,"sv_need.txt",row.names=F,quote=F,sep="\t")


exa<-read.table("sample_gene.txt",header=T,sep="\t")
row<-as.vector(unique(as.vector(exa[,3])))
col<-as.vector(unique(as.vector(exa[,1])))
pinpu<-matrix(rep("",length(row)*length(col)),nrow=length(row))
rownames(pinpu)<-row
colnames(pinpu)<-col

for(i in 1:length(row)){
for(j in 1:length(col)){
a<-which(exa$Gene==rownames(pinpu)[i])
b<-which(exa$sample==colnames(pinpu)[j])
c<-intersect(a,b)
if(length(c)!=0){
if(length(unique(as.vector(exa[c,2])))!=1){
pinpu[i,j]<-paste0(unique(as.vector(exa[c,2])),collapse=",")
}else{
pinpu[i,j]=unique(as.vector(exa[c,2]))}
}
}
}
pinpu<-unique(pinpu)
