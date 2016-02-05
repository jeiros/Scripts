#!/usr/bin/env Rscript


# run as ./script.R distance_file1.dat distance_file2.dat plot.png


library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)


args

d1 <- read.table(args[1])
d2 <- read.table(args[2])


dfNew <- rbind(data.frame(d1, Group = "d1"),
               data.frame(d2, Group = "d2"))

substr_1 <- substr(args[1], nchar(args[1]) - 16, nchar(args[1]) - 14)
substr_2 <- substr(args[2], nchar(args[2]) - 16, nchar(args[2]) - 14)

p1 <- ggplot(data=dfNew, aes(x=V2, color=Group)) +
    geom_line(stat='density', size = 1.5) +
    labs(y='Density', x='Distance to Ca (Å)') +
    scale_x_continuous(breaks = seq(0, 10)) +
    expand_limits(x=0) +
    theme_bw(15) +
    theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_rect(fill = 'white',
        colour = 'black')) +
    scale_colour_discrete("Atom",
        labels = c(paste(substr_1, "WT", collapse=""),
        paste(substr_2, "SP23/SP24", collapse="")))



ggsave(paste(toString(args[3]),".eps",sep=''), plot=p1, dpi=900, width=10, height=10, units='in')

# colnames(d1) = colnames(d2) = c("step", "Temp")
# dfNew <- rbind(data.frame(d1, Group = "d1"),
#                data.frame(d2, Group = "d2"))
# cdat1 <- ddply(dfNew, "Group", summarise, rating.mean = mean(Temp))

# p1 <- ggplot(dfNew, aes(x = Temp, color = Group, fill = Group, alpha = Group)) 
# + geom_density() + ylab("Density") + xlab("Temperature (K)") + theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(fill = "white", colour = "black")) 
# + ggtitle("Temperature Distribution") + geom_vline(data = cdat1, aes(xintercept = rating.mean, color = Group), linetype = "dashed", size = 1)
# + scale_colour_discrete("MD run",labels = c("4000 frames", "2500 frames")) + scale_alpha_manual("MD run", values = c("d1" = 0.8, "d2" = 0.4), labels = c("4000 frames", "2500 frames"))
# + scale_fill_discrete("MD run",labels = c("4000 frames", "2500 frames"))


# pplot <- ggplot(data=df, aes(x=dist, color = Group)) + geom_line(stat='density') + theme_classic(15) + ylab("Density") +
#   theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(fill = "white", colour = "black")) +
#   scale_colour_discrete("cTn complex",labels = c("WT", "SP23/SP24")) + labs(y='Density', x='Distance to Ca (Å)', title='Thr 71') + scale_x_continuous(breaks = seq(0,10)) + expand_limits(x=0)


# colnames(d1) = colnames(d2) = c("step", "Temp")
# dfNew <- rbind(data.frame(d1, Group = "d1"),
#                data.frame(d2, Group = "d2"))
# cdat1 <- ddply(dfNew, "Group", summarise, rating.mean = mean(Temp))