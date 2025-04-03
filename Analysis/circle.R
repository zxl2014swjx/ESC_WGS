library(circlize)
set.seed(999)
n <- 1000
a <- data.frame(factors = sample(letters[1:8], n, replace = TRUE), x = rnorm(n), y = runif(n))



par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.6)

circos.par(track.height = 0.1)

circos.initialize(factors = a$factors, x = a$x)


circos.trackPlotRegion(factors = a$factors, y = a$y,
panel.fun = function(x, y) {
circos.axis()})


col <- rep(c("#FF0000", "#00FF00"), 4) 
circos.trackPoints(a$factors, a$x, a$y, col = col, pch = 16, cex = 0.5) 
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)

circos.text(1, 0.5, "right", sector.index = "a")
circos.clear()

par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.6)
circos.par(track.height = 0.1)

circos.initialize(factors = a$factors, x = a$x)
circos.trackPlotRegion(factors = a$factors, y = a$y,
	panel.fun = function(x, y) {
		circos.axis()})
col <- rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(a$factors, a$x, a$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
circos.text(1, 0.5, "right", sector.index = "a")
bg.col <- rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(a$factors, a$x, bg.col = bg.col, col = NA)
circos.clear()

par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.6)
circos.par(track.height = 0.1)
circos.initialize(factors = a$factors, x = a$x)
circos.trackPlotRegion(factors = a$factors, y = a$y,
panel.fun = function(x, y) {
circos.axis()})
col <- rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(a$factors, a$x, a$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
circos.text(1, 0.5, "right", sector.index = "a")
bg.col <- rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(a$factors, a$x, bg.col = bg.col, col = NA)
circos.trackPlotRegion(factors = a$factors, x = a$x, y = a$y,
panel.fun = function(x, y) {
grey = c("#FFFFFF", "#CCCCCC", "#999999")
sector.index = get.cell.meta.data("sector.index") 
xlim = get.cell.meta.data("xlim")
ylim = get.cell.meta.data("ylim")
circos.text(mean(xlim), mean(ylim), sector.index)
circos.points(x[1:10], y[1:10], col = "red", pch = 16, cex = 0.6)
circos.points(x[11:20], y[11:20], col = "blue", cex = 0.6)})
circos.clear()


par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.6)
circos.par(track.height = 0.1)
circos.initialize(factors = a$factors, x = a$x)
circos.trackPlotRegion(factors = a$factors, y = a$y,panel.fun = function(x, y) {
circos.axis()})
col <- rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(a$factors, a$x, a$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
circos.text(1, 0.5, "right", sector.index = "a")
circos.trackPoints(a$factors, a$x, a$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
circos.text(1, 0.5, "right", sector.index = "a")
bg.col <- rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(a$factors, a$x, bg.col = bg.col, col = NA)
circos.trackPlotRegion(factors = a$factors, x = a$x, y = a$y,panel.fun = function(x, y) {
grey = c("#FFFFFF", "#CCCCCC", "#999999")
sector.index = get.cell.meta.data("sector.index")
xlim = get.cell.meta.data("xlim")
ylim = get.cell.meta.data("ylim")
circos.text(mean(xlim), mean(ylim), sector.index)
circos.points(x[1:10], y[1:10], col = "red", pch = 16, cex = 0.6)
circos.points(x[11:20], y[11:20], col = "blue", cex = 0.6)})
#update第2个track中标记为d的sector
circos.updatePlotRegion(sector.index = "d", track.index = 2)
circos.points(x = -2:2, y = rep(0, 5))
xlim <- get.cell.meta.data("xlim")
ylim <- get.cell.meta.data("ylim")
circos.text(mean(xlim), mean(ylim), "updated")
circos.clear()


par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.6)
circos.par(track.height = 0.1)
circos.initialize(factors = a$factors, x = a$x)
circos.trackPlotRegion(factors = a$factors, y = a$y,panel.fun = function(x, y) { 
circos.axis()})
col <- rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(a$factors, a$x, a$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
circos.text(1, 0.5, "right", sector.index = "a")
circos.trackPoints(a$factors, a$x, a$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
circos.text(1, 0.5, "right", sector.index = "a")
bg.col <- rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(a$factors, a$x, bg.col = bg.col, col = NA)
circos.trackPlotRegion(factors = a$factors, x = a$x, y = a$y,panel.fun = function(x, y) {
	grey = c("#FFFFFF", "#CCCCCC", "#999999")
sector.index = get.cell.meta.data("sector.index")
xlim = get.cell.meta.data("xlim")
ylim = get.cell.meta.data("ylim")
circos.text(mean(xlim), mean(ylim), sector.index)
circos.points(x[1:10], y[1:10], col = "red", pch = 16, cex = 0.6)
circos.points(x[11:20], y[11:20], col = "blue", cex = 0.6)})
# update第2个track中标记为d的sector
circos.updatePlotRegion(sector.index = "d", track.index = 2)
circos.points(x = -2:2, y = rep(0, 5))
xlim <- get.cell.meta.data("xlim")
ylim <- get.cell.meta.data("ylim")
circos.text(mean(xlim), mean(ylim), "updated")
circos.clear()




par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.6)
circos.par(track.height = 0.1)
circos.initialize(factors = a$factors, x = a$x)
circos.trackPlotRegion(factors = a$factors, y = a$y,panel.fun = function(x, y) { 
circos.axis()})
col <- rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(a$factors, a$x, a$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
circos.text(1, 0.5, "right", sector.index = "a")
circos.trackPoints(a$factors, a$x, a$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
circos.text(1, 0.5, "right", sector.index = "a")
bg.col <- rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(a$factors, a$x, bg.col = bg.col, col = NA)
circos.trackPlotRegion(factors = a$factors, x = a$x, y = a$y,panel.fun = function(x, y) {
grey = c("#FFFFFF", "#CCCCCC", "#999999")
sector.index = get.cell.meta.data("sector.index")
xlim = get.cell.meta.data("xlim")
ylim = get.cell.meta.data("ylim")
circos.text(mean(xlim), mean(ylim), sector.index)
circos.points(x[1:10], y[1:10], col = "red", pch = 16, cex = 0.6)
circos.points(x[11:20], y[11:20], col = "blue", cex = 0.6)})



set.seed(1234)
data <- matrix(rnorm(100 * 10), nrow = 10, ncol = 100)
col <- colorRamp2(c(-2, 0, 2), c("green", "black", "red"))
factors <- rep(letters[1:2], times = c(30, 70))
data_list <- list(a = data[, factors == "a"], b = data[, factors == "b"])
dend_list <- list(a = as.dendrogram(hclust(dist(t(data_list[["a"]])))),
	b = as.dendrogram(hclust(dist(t(data_list[["b"]])))))
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 5)
circos.initialize(factors = factors, xlim = cbind(c(0, 0), table(factors)))
circos.track(ylim = c(0, 10), bg.border = NA,
	panel.fun = function(x, y) {
	sector.index = get.cell.meta.data("sector.index")
	d = data_list[[sector.index]]
	dend = dend_list[[sector.index]]
	d2 = d[, order.dendrogram(dend)]
	col_data = col(d2)
	nr = nrow(d2)
	nc = ncol(d2)
	for (i in 1:nr) {
		circos.rect(1:nc - 1, rep(nr - i, nc), 1:nc, rep(nr - i + 1, nc),border = col_data[i, ], col = col_data[i, ]) }})
max_height <- max(sapply(dend_list, function(x) attr(x, "height")))
circos.track(ylim = c(0, max_height),bg.border = NA, track.height = 0.3,panel.fun = function(x, y) {sector.index = get.cell.meta.data("sector.index")
dend = dend_list[[sector.index]]
	circos.dendrogram(dend, max_height = max_height)})
circos.clear()


library(circlize) 

circos.initializeWithIdeogram()

bed = generateRandomBed(nr = 200, nc = 4)
circos.genomicPosTransformLines(bed, posTransform = posTransform.default, horizontalLine = "top")
om = circos.par("track.margin")
oc = circos.par("cell.padding")
circos.par(track.margin = c(om[1], 0), cell.padding = c(0, 0, 0, 0))
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("blue", "white", "red"))
circos.genomicTrackPlotRegion(bed, stack = TRUE, panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = f(value[[1]]), 
        border = f(value[[1]]), lwd = 0.1, posTransform = posTransform.default, ...)
}, bg.border = NA, track.height = 0.1)
circos.par(track.margin = om, cell.padding = oc)

bed = generateRandomBed(nr = 500, fun = function(k) runif(k)*sample(c(-1, 1), k, replace = TRUE))
circos.genomicTrackPlotRegion(bed, ylim = c(-1, 1), panel.fun = function(region, value, ...) {
    col = ifelse(value[[1]] > 0, "red", "green")
    circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
    cell.xlim = get.cell.meta.data("cell.xlim")
    for(h in c(-1, -0.5, 0, 0.5, 1)) {
        circos.lines(cell.xlim, c(h, h), col = "#00000040")
    }
}, track.height = 0.1)

bed = generateRandomBed(nr = 500, fun = function(k) rnorm(k, 0, 50))
circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
    x = (region[[2]] + region[[1]]) / 2
    y = value[[1]]
    loess.fit = loess(y ~ x)
    loess.predict = predict(loess.fit, x, se = TRUE)
    d1 = c(x, rev(x))
    d2 = c(loess.predict$fit + loess.predict$se.fit, rev(loess.predict$fit - loess.predict$se.fit))
    circos.polygon(d1, d2, col = "#CCCCCC", border = NA)
    circos.points(x, y, pch = 16, cex = 0.5)
    circos.lines(x, loess.predict$fit)
}, track.height = 0.1)


bed_list = list(generateRandomBed(nr = 500, fun = function(k) runif(k)),
                generateRandomBed(nr = 500, fun = function(k) runif(k)))
col = c("#FF000040", "#0000FF40")
circos.genomicTrackPlotRegion(bed_list, ylim = c(-1, 1), panel.fun = function(region, value, ...) {
    i = getI(...)
    if(i == 1) {
        circos.genomicLines(region, value, area = TRUE, baseline = 0, col = "orange", border = NA, ...)
    } else {
        circos.genomicLines(region, -value, area = TRUE, baseline = 0, col = "yellow", border = NA, ...)
    }
}, track.height = 0.1)

region1 = generateRandomBed(nr = 1000); region1 = region1[sample(nrow(region1), 20), ]
region2 = generateRandomBed(nr = 1000); region2 = region2[sample(nrow(region2), 20), ]
circos.genomicLink(region1, region2, col = sample(10, 20, replace = TRUE))

highlight.chromosome("chr1")

circos.clear()