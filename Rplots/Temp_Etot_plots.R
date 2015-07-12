
library(ggplot2)
library(plyr)
setwd("~/Troponin/IAN_Troponin/completehowarthcut/phospho/hmr_runs/CTnI_hmr/run1/S1P/")


d1 <- read.table("~/Troponin/IAN_Troponin/completehowarthcut/phospho/hmr_runs/CTnI_hmr/run1/S1P/000-050/summary.TEMP", quote="\"")

d2 <- read.table("~/Troponin/IAN_Troponin/completehowarthcut/phospho/hmr_runs/CTnI_hmr/run1/S1P/000-050rep/summary.TEMP", quote="\"")

colnames(d1) = colnames(d2) = c("step", "Temp")
dfNew <- rbind(data.frame(d1, Group = "d1"),
               data.frame(d2, Group = "d2"))
cdat1 <- ddply(dfNew, "Group", summarise, rating.mean = mean(Temp))

# Density plot
p1 <- ggplot(dfNew, aes(x = Temp, color = Group, fill = Group, alpha = Group)) + geom_density() + ylab("Density") + xlab("Temperature (K)") + theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(fill = "white", colour = "black")) + ggtitle("Temperature Distribution") + geom_vline(data = cdat1, aes(xintercept = rating.mean, color = Group), linetype = "dashed", size = 1) + scale_colour_discrete("MD run",labels = c("4000 frames", "2500 frames")) + scale_alpha_manual("MD run", values = c("d1" = 0.8, "d2" = 0.4), labels = c("4000 frames", "2500 frames")) + scale_fill_discrete("MD run",labels = c("4000 frames", "2500 frames"))

# Timeseries plot
p2 <- ggplot(data = dfNew, aes(x = step, y = Temp, color = Group)) + geom_line() + scale_color_discrete("MD run",labels = c("4000 frames", "2500 frames")) + xlab("Time (ps)") + theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(fill = "white", colour = "black")) + ggtitle("Temperature") + ylab("Temperature (K)")




etot1 <- read.table("~/Troponin/IAN_Troponin/completehowarthcut/phospho/hmr_runs/CTnI_hmr/run1/S1P/000-050/summary.ETOT", quote="\"")
etot2 <- read.table("~/Troponin/IAN_Troponin/completehowarthcut/phospho/hmr_runs/CTnI_hmr/run1/S1P/000-050rep/summary.ETOT", quote="\"")
colnames(etot1) = colnames(etot2) = c("step", "Etot")

dfNew_etot <- rbind(data.frame(etot1, Group = "d1"),
               data.frame(etot2, Group = "d2"))
cdat2 <- ddply(dfNew_etot, "Group", summarise, rating.mean = mean(Etot))

# Density plot
p3 <- ggplot(dfNew_etot, aes(x = Etot, color = Group, fill = Group, alpha = Group)) + geom_density() + ylab("Density") + xlab("Energy (kcal/mol)") + theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(fill = "white", colour = "black")) + ggtitle("Total Energy Distribution") + geom_vline(data = cdat2, aes(xintercept = rating.mean, color = Group), linetype = "dashed", size = 1) + scale_alpha_manual("MD run", values = c("d1" = 0.8, "d2" = 0.4), labels = c("4000 frames", "2500 frames")) + scale_fill_discrete("MD run", labels = c("4000 frames", "2500 frames")) + scale_colour_discrete("MD run", labels = c("4000 frames", "2500 frames"))


# Timeseries plot
p4 <- ggplot(data = dfNew_etot, aes(x = step, y = Etot, color = Group)) + geom_line() + scale_color_discrete("MD run",labels = c("4000 frames", "2500 frames")) + xlab("Time (ps)") + theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(fill = "white", colour = "black")) + ggtitle("Total Energy") + ylab("Total Energy (kcal/    mol)")