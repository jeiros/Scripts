errorbars_gnuplot <- read.delim("~/errorbars_gnuplot", header=FALSE)

colnames(errorbars_gnuplot) <- c("Residue", "Mean", "Std")

ggplot(errorbars_gnuplot, aes(x = Residue, y = Mean)) + ylab("RMSF (Å)") + xlim(1,420) + 
  geom_ribbon(aes(ymin = Mean - Std, ymax = Mean + Std), alpha = 0.2) + geom_line() +
  scale_x_continuous(breaks = c(0,50,100,150,200,250,300,350,400,419))

