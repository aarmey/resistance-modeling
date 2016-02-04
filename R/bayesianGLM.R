library("ggplot2")
require("MASS")
library("Hmisc")
source("multiplot.R")
source("runbGLM.R")

results1 <-  as.data.frame(EGFR$Sol)
results2 <-  as.data.frame(HER2$Sol)

p1 <- qplot(results1$JNK, results1$Erk, geom = "density2d")
p2 <- qplot(results2$Akt, results2$Erk, 
            color = cut(results2$Erk, breaks=c(-3, 0, 3)), geom = "density2d")
p3 <- qplot(results1$Erk, results1$cJun, geom = "density2d")
p4 <- qplot(results2$P38, results2$JNK, geom = "density2d")

multiplot(p1, p2, p3, p4, cols = 2)

mattx <- as.matrix(results1)
mattx <- mattx[order(mattx[,3]),]

viol <- reshape2::melt(results2, id.vars = NULL)

#gg <- ggplot(viol, aes(x = variable, y = value)) + geom_violin(scale = 'width')
#plot(gg)
#nba_heatmap <- heatmap(mattx, Rowv=NA, Colv=NA, col = terrain.colors(256), scale="none", margins=c(5,10))

#pairs(mattx, pch='.')

ss <- summary(mattx)

