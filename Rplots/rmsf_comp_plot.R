library(ggplot2)





d1 <- read.delim("~/Troponin/IAN_Troponin/completehowarthcut/salted/Analysis/RMSF/rmsf_std.dat", header=FALSE, stringsAsFactors=FALSE, skip = 1)
d2 <- read.delim("~/Troponin/IAN_Troponin/completehowarthcut/phospho/hmr_runs/Analysis/RMSF/mean_std_ci.dat", header=FALSE, stringsAsFactors=FALSE, skip = 1)

colnames(d1) <- colnames(d2) <- c("res", "mean", "std", "ci")

dfNew <- rbind(data.frame(d1, Group = "d1"),
               data.frame(d2, Group = "d2"))

p <- ggplot(data = dfNew, aes(x = res, y = mean, color = Group)) + geom_line() + scale_color_manual("MD run", labels = c("Normal", "Phosphorylation"), values = c("red", "blue")) + geom_ribbon(aes(ymin = mean - std, ymax = mean + std, colour = NULL, fill = Group, alpha = 0.1)) + guides(alpha = F) + scale_fill_discrete("MD run", labels = c("Normal", "Phosphorylation")) + ylab("RMSF (Å)") + xlab("Residue number") + scale_x_continuous(breaks = c(0,100,200,300,400,161,249,419)) + ggtitle("RMSF comparison")



