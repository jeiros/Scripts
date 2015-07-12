library(ggplot2)





`d1` <- read.table("~/Troponin/IAN_Troponin/completehowarthcut/phospho/hmr_runs/rmsf_CTnI_hmr.000-400ns.dat", quote="\"")
`d2` <- read.table("~/Troponin/IAN_Troponin/completehowarthcut/phospho/hmr_runs/rmsf_CTnI_runs.000-400ns.dat", quote="\"")
`d3` <- read.table("~/Troponin/IAN_Troponin/completehowarthcut/phospho/hmr_runs/rmsf_CTnT_hmr.000-450ns.dat", quote="\"")
`d4` <- read.table("~/Troponin/IAN_Troponin/completehowarthcut/phospho/hmr_runs/rmsf_CTnT_runs.000-450ns.dat", quote="\"")

colnames(d1) <- colnames(d2) <- colnames(d3) <- colnames(d4) <- c("res", "rmsf")


dfNew <- rbind(data.frame(d1, Group = "d1"),
               data.frame(d2, Group = "d2"),
			   data.frame(d3, Group = "d3"),
               data.frame(d4, Group = "d4"))


TnC <- rbind(data.frame(d1[1:161,], Group = "d1"),
               data.frame(d2[1:161,], Group = "d2"),
			   data.frame(d3[1:161,], Group = "d3"),
               data.frame(d4[1:161,], Group = "d4"))

TnT <- rbind(data.frame(d1[162:248,], Group = "d1"),
               data.frame(d2[162:248,], Group = "d2"),
			   data.frame(d3[162:248,], Group = "d3"),
               data.frame(d4[162:248,], Group = "d4"))
TnT$res <- TnT$res - 161

TnI <- rbind(data.frame(d1[249:419,], Group = "d1"),
               data.frame(d2[249:419,], Group = "d2"),
			   data.frame(d3[249:419,], Group = "d3"),
               data.frame(d4[249:419,], Group = "d4"))
TnI$res <- TnI$res - 248


doplot <- function(df) {
	library(ggplot2)
	#legend <- theme(legend.position = c(1,1), legend.justification = c(1,1),
					 #legend.background = element_rect(fill = "white", colour = "black"))
  last_res <- max(df$res)
  
  
	legend_labels <- scale_color_discrete("Run", labels = c("CTnI_hmr", "CTnI_runs", "CTnT_hmr", "CTnT_runs"))

  axis_breaks <- scale_x_continuous(breaks = round(seq(from = 1, to = last_res, length.out = 10)))
  
	plot <- ggplot(df, aes(x = res, y = rmsf, color = Group))+ geom_line() + legend_labels +
		ylab("RMSF (Å)") + xlab("Residue number") + axis_breaks

	

	return(plot)
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


tnc <- doplot(TnC) + ggtitle("TnC")
tni <- doplot(TnI) + ggtitle("TnI")
tnt <- doplot(TnT) + ggtitle("TnT")
ppi = 900
png("myplot.png", width = 8*ppi, height = 8*ppi, res = ppi)
multiplot(tnc, tni, tnt)
dev.off()