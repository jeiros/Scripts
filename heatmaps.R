
library(ggplot2)
library(scales)


colnames(cont) <- c("Res1", "Res2", "Contact", "NContacts")

cont$Contact <- rescale(cont$Contact)

ggplot(cont, aes(x = Res1, y = Res2)) + geom_tile(aes(fill = Contact)) + scale_x_continuous(breaks = seq.int(min(cont$Res1),max(cont$Res1), by = 2)) + scale_y_continuous(breaks = seq.int(min(cont$Res2),max(cont$Res2),by=1)) + theme_bw() + ylab("cTnI residue number") + xlab("cTnC residue number") + ggtitle("Phosphorylated system CTnI_runs") + scale_fill_gradient(low = "white", high = "steelblue")

ggsave(args[1] + ".png")