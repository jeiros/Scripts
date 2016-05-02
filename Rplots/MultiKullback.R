library(ggplot2)
library(reshape2)

setwd("~/wt_data/Analysis/KLD/")
pc1 <- read.table("./KL-PC1_time.dat", quote="\"", skip = 1)
pc2 <- read.table("./KL-PC2_time.dat", quote="\"", skip = 1)
pc3 <- read.table("./KL-PC3_time.dat", quote="\"", skip = 1)
pc4 <- read.table("./KL-PC4_time.dat", quote="\"", skip = 1)
pc5 <- read.table("./KL-PC5_time.dat", quote="\"", skip = 1)
pc6 <- read.table("./KL-PC6_time.dat", quote="\"", skip = 1)
pc7 <- read.table("./KL-PC7_time.dat", quote="\"", skip = 1)
pc8 <- read.table("./KL-PC8_time.dat", quote="\"", skip = 1)
pc9 <- read.table("./KL-PC9_time.dat", quote="\"", skip = 1)




pc1[pc1 == 0] <- NA
pc2[pc2 == 0] <- NA
pc3[pc3 == 0] <- NA
pc4[pc4 == 0] <- NA
pc5[pc5 == 0] <- NA
pc6[pc6 == 0] <- NA
pc7[pc7 == 0] <- NA
pc8[pc8 == 0] <- NA
pc9[pc9 == 0] <- NA


colnames(pc1) <- colnames(pc2) <- colnames(pc3) <- colnames(pc4) <- colnames(pc5) <- colnames(pc6) <-colnames(pc7) <- colnames(pc8) <- colnames(pc9) <- c("Time", "1-2", "1-3", "1-4", "1-5", "1-6", "1-7", "1-8", "1-9")







d1 <- melt(pc1, id.vars = "Time")
d2 <- melt(pc2, id.vars = "Time")
d3 <- melt(pc3, id.vars = "Time")
d4 <- melt(pc4, id.vars = "Time")
d5 <- melt(pc5, id.vars = "Time")
d6 <- melt(pc6, id.vars = "Time")
d7 <- melt(pc7, id.vars = "Time")
d8 <- melt(pc8, id.vars = "Time")
d9 <- melt(pc9, id.vars = "Time")

p1 <- ggplot(d1, aes(x = Time, y = value, colour = variable)) +
    geom_line(size=1.5) +
    ylab("KLD") +
    scale_color_discrete("") +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(breaks = seq(0, 1500, by = 150)) +
    ggtitle("PC1") + xlab("Time (ns)") +
    theme_bw() + theme(legend.position="none")

p2 <- ggplot(d2, aes(x = Time, y = value, colour = variable)) +
    geom_line(size=1.5) +
    ylab("KLD") +
    scale_color_discrete("") +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(breaks = seq(0, 1500, by = 150)) +
    ggtitle("PC2") + xlab("Time (ns)") +
    theme_bw() + theme(legend.position="none")

p3 <- ggplot(d3, aes(x = Time, y = value, colour = variable)) +
    geom_line(size=1.5) +
    ylab("KLD") +
    scale_color_discrete("") +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(breaks = seq(0, 1500, by = 150)) +
    ggtitle("PC3") + xlab("Time (ns)") +
    theme_bw() + theme(legend.position="none")

p4 <- ggplot(d4, aes(x = Time, y = value, colour = variable)) +
    geom_line(size=1.5) +
    ylab("KLD") +
    scale_color_discrete("") +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(breaks = seq(0, 1500, by = 150)) +
    ggtitle("PC4") + xlab("Time (ns)") +
    theme_bw() + theme(legend.position="none")

p5 <- ggplot(d5, aes(x = Time, y = value, colour = variable)) +
    geom_line(size=1.5) +
    ylab("KLD") +
    scale_color_discrete("") +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(breaks = seq(0, 1500, by = 150)) +
    ggtitle("PC5") + xlab("Time (ns)") +
    theme_bw() + theme(legend.position="none")

p6 <- ggplot(d6, aes(x = Time, y = value, colour = variable)) +
    geom_line(size=1.5) +
    ylab("KLD") +
    scale_color_discrete("") +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(breaks = seq(0, 1500, by = 150)) +
    ggtitle("PC6") + xlab("Time (ns)") +
    theme_bw() + theme(legend.position="none")

p7 <- ggplot(d7, aes(x = Time, y = value, colour = variable)) +
    geom_line(size=1.5) +
    ylab("KLD") +
    scale_color_discrete("") +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(breaks = seq(0, 1500, by = 150)) +
    ggtitle("PC7") + xlab("Time (ns)") +
    theme_bw() + theme(legend.position="none")

p8 <- ggplot(d8, aes(x = Time, y = value, colour = variable)) +
    geom_line(size=1.5) +
    ylab("KLD") +
    scale_color_discrete("") +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(breaks = seq(0, 1500, by = 150)) +
    ggtitle("PC8") + xlab("Time (ns)") +
    theme_bw() + theme(legend.position="none")

p9 <- ggplot(d9, aes(x = Time, y = value, colour = variable)) +
    geom_line(size=1.5) +
    ylab("KLD") +
    scale_color_discrete("") +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(breaks = seq(0, 1500, by = 150)) +
    ggtitle("PC9") + xlab("Time (ns)") +
    theme_bw() + theme(legend.position="none")

source("/Users/je714/Scripts/Rplots/Multiple_plot_function.R")

multiplot(p1, p4, p7, p2, p5, p8, p3, p6, p9, cols=3)

cairo_ps("~/Dropbox (Imperial)/Poster/Pictures/KLD.eps", width = 20, height = 20)
source("/Users/je714/Scripts/Rplots/Multiple_plot_function.R")
multiplot(p1,p2,p3)
dev.off()