library(ggplot2)
library(reshape2)




setwd("/home/je714/Troponin/IAN_Troponin/completehowarthcut/phospho/hmr_runs/try_readdata")




pc1 <- read.table("./KL-PC1.dat", quote="\"")
pc2 <- read.table("./KL-PC3.dat", quote="\"")
pc3 <- read.table("./KL-PC3.dat", quote="\"")


colnames(pc1) <- colnames(pc2) <- colnames(pc3) <- c("Frame", "1-2", "1-3", "2-3")


d1 <- melt(pc1, id.vars = "Frame")
d2 <- melt(pc2, id.vars = "Frame")
d3 <- melt(pc3, id.vars = "Frame")

p1 <- ggplot(d1, aes(x = Frame, y = value, colour = variable)) + geom_line() + ylab("Principal Component 1") + scale_color_discrete("") + theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(fill = "white", colour = "black"))

p2 <- ggplot(d2, aes(x = Frame, y = value, colour = variable)) + geom_line() + ylab("Principal Component 1") + scale_color_discrete("") + theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(fill = "white", colour = "black"))

p3 <- ggplot(d3, aes(x = Frame, y = value, colour = variable)) + geom_line() + ylab("Principal Component 1") + scale_color_discrete("") + theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(fill = "white", colour = "black"))

ppi=900
png("KL-PC.png", width = 8*ppi, height = 8*ppi, res = ppi)
multiplot(p1,p2,p3)
dev.off()