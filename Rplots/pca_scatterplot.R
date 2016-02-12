wt1 <- read.table("juan30mhr1/myevecs_WT_1.000-750ns.dat_time", header = T)
wt2 <- read.table("juan30mhr2/myevecs_WT_2.000-750ns.dat_time", header = T)
wt3 <- read.table("juan30mhr3/myevecs_WT_3.000-750ns.dat_time", header = T)
wt4 <- read.table("juan30mhr4/myevecs_WT_4.000-750ns.dat_time", header = T)

xmin <- min(wt1$Mode1, wt2$Mode1, wt3$Mode1, wt4$Mode1)
xmax <- max(wt1$Mode1, wt2$Mode1, wt3$Mode1, wt4$Mode1)
ymin <- min(wt1$Mode2, wt2$Mode2, wt3$Mode2, wt4$Mode2)
ymax <- max(wt1$Mode2, wt2$Mode2, wt3$Mode2, wt4$Mode2)



pca_plots <- function(data, title, x_max, x_min, y_max, y_min) {
  plot <- ggplot(data, aes(x=Mode1, y=Mode2, color = X0)) + geom_jitter() +
    labs(title=title, x="PC1 (Å)", y="PC2 (Å)") + theme_bw(15) +
    coord_cartesian(xlim=c(x_min-20,x_max+20), ylim = c(y_min-20, y_max+20)) +
    scale_colour_continuous("Time (ns)") +
    theme(legend.position="none")
  return(plot)
}

# PCA projection plots
p_wt1 <- pca_plots(wt1, "WT run1", x_max=xmax, x_min=xmin, y_max=ymax, y_min=ymin)
p_wt2 <- pca_plots(wt2, "WT run2", x_max=xmax, x_min=xmin, y_max=ymax, y_min=ymin)
p_wt3 <- pca_plots(wt3, "WT run3", x_max=xmax, x_min=xmin, y_max=ymax, y_min=ymin)
p_wt4 <- pca_plots(wt4, "WT run4", x_max=xmax, x_min=xmin, y_max=ymax, y_min=ymin)

source("~/Scripts/Rplots/Multiple_plot_function.R")
multiplot(p_wt1, p_wt3, p_wt2, p_wt4, cols = 2)


# PC1 distribution plot

wt1 <- wt1[,1:2]
wt2 <- wt2[,1:2]
wt3 <- wt3[,1:2]
wt4 <- wt4[,1:2]


ggplot() + theme_bw(15) +
  geom_line(stat='density', data=wt1, aes(x=Mode1, color = 'Run1'), size = 1.5) +
  geom_line(stat='density', data=wt2, aes(x=Mode1, color = 'Run2'), size = 1.5) +
  geom_line(stat='density', data=wt3, aes(x=Mode1, color = 'Run3'), size = 1.5) +
  geom_line(stat='density', data=wt4, aes(x=Mode1, color = 'Run4'), size = 1.5) +
  labs(title='PC1 distribution', y='Density', x='PC1 (Å)') +
  scale_color_discrete("WT run") +
  theme(legend.justification=c(1,0), legend.position=c(.95,.8))