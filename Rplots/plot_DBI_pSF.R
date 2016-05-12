library(ggplot2)


reduce_matrix <- function(data) {
  
  unique_clusterlist <- unique(data[,1])
  new_data_frame <- matrix(ncol = 3, nrow = length(unique_clusterlist))
  new_data_frame[,1] <- unique_clusterlist
  iter_no <- 1
  
  for (i in unique_clusterlist) {
    average_DBI <- mean(data[data[,1]==i,2])
    average_pSF <- mean(data[data[,1]==i,3])
    
    new_data_frame[iter_no,2] <- average_DBI
    new_data_frame[iter_no,3] <- average_pSF
    
    iter_no <- iter_no + 1
  }
  return(data.frame(new_data_frame))
}

plotdbipsf <- function(data){

  colnames(data) <- c("Clusters", "DBI", "pSF")

  x_axis <- scale_x_continuous(breaks=round(seq(from=min(data$Clusters),
                                                to=max(data$Clusters),
                                                length.out = 20)))

  g_top <- ggplot(data, aes(x = Clusters, y = pSF)) +
           geom_line(size=1) +
           geom_point(size=2.5) +
           theme_bw() +
           labs(y = "pSF", x='Clusters') +
           x_axis
  g_bot <- ggplot(data, aes(x = Clusters, y = DBI)) +
           geom_line(size=1) +
           geom_point(size=2.5) +
           theme_bw() +
           labs(y = "DBI", x='Clusters') +
           x_axis
  source("~/Scripts/Rplots/Multiple_plot_function.R")
  ## plot graphs and set relative heights
  pdf("DBI_PSF.pdf", width = 15, height = 12)
  multiplot(g_top, g_bot)
  dev.off()
}