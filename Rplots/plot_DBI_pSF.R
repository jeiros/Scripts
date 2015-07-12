library(ggplot2)
library(gridExtra)

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
  g.top <- ggplot(data, aes(x = Clusters, y = pSF)) +
    geom_line() +
    theme_bw() +
    theme(plot.margin = unit(c(1,5,-30,-2),units="points"),
          axis.title.y = element_text(vjust =0.25)) +
    labs(y = "pSF") + scale_x_continuous(breaks = seq(2,data[nrow(data),1],by=2))

  g.bottom <- ggplot(data, aes(x = Clusters, y = DBI)) +
    geom_line() +
    theme_bw() +
    theme(plot.margin = unit(c(0,5,1,6),units="points")) +
    labs(x = "Clusters", y = "DBI") + scale_x_continuous(breaks = seq(0,data[nrow(data),1],by=2))

## plot graphs and set relative heights
  grid.arrange(g.top,g.bottom, heights = c(45/100,55/100))
  #g <- arrangeGrob(g.top, g.bottom, nrow = 2)
  #ggsave("dbipsf_prueba.png", g, width = 12, height = 10, dpi = 800)

}