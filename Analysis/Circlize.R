
library(stringr)	
library(circlize)	
library(grid)	
library(ComplexHeatmap)	

#rm(list=ls())
sample<-read.table("E:/zxl/11_WGS(9.16-)/WGS/data/sample.txt",header=T)

####################################

for(i in i:dim(sample)[1]){
setwd("E:/zxl/11_WGS(9.16-)/WGS/data/bed")
CNV<-list.files(pattern="*.seg") 
cnv<-read.delim(CNV[i],header=TRUE)
cnv<-subset(cnv,select=c(1,2,3,5))
colnames(cnv)<-c("chr","start","end","log")
cnv[,1]<-paste0("chr",cnv[,1])


setwd("E:/zxl/11_WGS(9.16-)/WGS/data/SV")
SV<-list.files(pattern="*.rp") 
sv<-read.delim(SV[i],header=F,sep=" ")
sv<-subset(sv,select=c(1,6,7,8))
colnames(sv)<-c("event","chr","start","end")
sv<-sv[sv[,3]<sv[,4],]
sv[,2]<-paste0("chr",sv[,2])

a<-as.vector(sv[,1])
a[grep("tandem_dup",a,perl=TRUE)]<-4
a[grep("inver*",a,perl=TRUE)]<-3
a[grep("del*",a,perl=TRUE)]<-1
a[grep("ins*",a,perl=TRUE)]<-2
sv<-cbind(sv[,2:4],as.data.frame(a))
sv[[4]] <- as.numeric(sv[[4]])


setwd("E:/zxl/11_WGS(9.16-)/WGS/data/fusion_result")
FUSION<-list.files(pattern="*.fusions") 
fusion<-read.delim(FUSION[i],header=F)
fusion<-subset(fusion,select=c(4,5,9,10))
ICT<-fusion[as.matrix(fusion[,1])!=as.matrix(fusion[,3]),]
if(dim(ICT)[1]!=0){
ICT[,1]<-paste0("chr",ICT[,1])
ICT[,3]<-paste0("chr",ICT[,3])
}

setwd("E:/zxl/11_WGS(9.16-)/WGS/data")
pdf(file=paste0(sample[i,1],".pdf"))
circos.initializeWithIdeogram()
circos.genomicPosTransformLines(cnv)
circos.genomicTrackPlotRegion(
	cnv, ylim = c(-1, 1), panel.fun = function(region, value, ...) {
    col = ifelse(value[[1]] > 0, "red", "green")
    circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
    cell.xlim = get.cell.meta.data("cell.xlim")
    for(h in c(-1, -0.5, 0, 0.5, 1)) {
        circos.lines(cell.xlim, c(h, h), col = "#00000040")
    }
}, track.height = 0.1)


cnv_legend <- Legend(
	at = c(-1, 1), labels = c('del', 'amp'), labels_gp = gpar(fontsize = 8), title = 'CNV type', title_gp = gpar(fontsize = 9),  
	grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'), type = 'points', pch = NA, background = c('#FB564A', '#74C476'))


circos.genomicPosTransformLines(sv)
color_assign <- colorRamp2(breaks = c(1:4), colors = c('red','green','blue','purple'))
circos.genomicTrackPlotRegion(sv, stack = TRUE, panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), 
        border = color_assign(value[[1]]), lwd = 0.1, posTransform = posTransform.default, ...)
}, track.height = 0.1)

sv_legend <- Legend(
	at = c(1:4), labels = c('del','ins','inver','tandem_dup'),
	labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'), 
	type = rep('line',4), pch = NA, background = c('red','green','blue','purple'),
	legend_gp = gpar(col = c('red','green','blue','purple'), lwd = 1),
	title = 'SV type', title_gp = gpar(fontsize = 9))


circos.genomicLink(ICT[c(1, 2, 2)], ICT[c(3, 4, 4)],col = '#FF46346E',lwd = 3)



pushViewport(viewport(x = 0.9, y = 0.85))
grid.draw(cnv_legend)
y_coord <- y_coord
upViewport()

pushViewport(viewport(x = 0.1, y = 0.85))
grid.draw(sv_legend)
y_coord <- y_coord
upViewport()


circos.clear()
dev.off()
}

