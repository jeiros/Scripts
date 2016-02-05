library(ggplot2)
library(reshape2)


pc1 <- read.table("./KL-PC1_time.dat", quote="\"", skip = 1)
pc2 <- read.table("./KL-PC2_time.dat", quote="\"", skip = 1)
pc3 <- read.table("./KL-PC3_time.dat", quote="\"", skip = 1)


pc1[pc1 == 0] <- NA
pc2[pc2 == 0] <- NA
pc3[pc3 == 0] <- NA

colnames(pc1) <- colnames(pc2) <- colnames(pc3) <- c("Time", "1-2", "1-3", "1-4", "1-5", "1-6", "1-7")


d1 <- melt(pc1, id.vars = "Time")
d2 <- melt(pc2, id.vars = "Time")
d3 <- melt(pc3, id.vars = "Time")

p1 <- ggplot(d1, aes(x = Time, y = value, colour = variable)) + geom_line(size=1) + ylab("KLD") +
  scale_color_discrete("") + scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(breaks = seq(0, 750, by = 150)) + ggtitle("PC1") + xlab("Time (ns)") +
  theme_bw(15) + theme(legend.position="none") 

p2 <- ggplot(d2, aes(x = Time, y = value, colour = variable)) + geom_line(size=1) + ylab("KLD") +
  scale_color_discrete("WT runs") +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(breaks = seq(0, 750, by = 150)) + ggtitle("PC2") + xlab("Time (ns)") +
  theme_bw(15) + theme(legend.position="none") 

p3 <- ggplot(d3, aes(x = Time, y = value, colour = variable)) + geom_line(size=1) + ylab("KLD") +
  scale_color_discrete("WT runs") + scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(breaks = seq(0, 750, by = 150)) + ggtitle("PC3") + xlab("Time (ns)") +
  theme_bw(15) + theme(legend.justification=c(1,0), legend.position="bottom")

cairo_ps("KLD.eps", width = 10, height = 10)
source("/Users/je714/Scripts/Rplots/Multiple_plot_function.R")
multiplot(p1,p2,p3)
dev.off()