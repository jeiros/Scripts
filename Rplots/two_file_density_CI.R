#!/usr/bin/env Rscript

# Plot densities with 95% CI using bootstrap
# http://stats.stackexchange.com/a/207455/76846

# run as ./script.R distance_file1.dat distance_file2.dat title_plot red_label blue_label file_name

get_CI_shade  <- function(data, iterations, x_0=0, x_f=10) {
  evalpoints <- nrow(data)
  distances <- data[,2] # The distances are in the second column
  x <- na.omit(distances) # Clean any NA values 
  xax <- seq(x_0, x_f, length.out=evalpoints)
  d <- density(x, n=evalpoints, from=x_0, to=x_f)
  estimates <- matrix(NA, nrow=evalpoints, ncol=iterations)
  
  for (b in 1:iterations){
    print("Iteration number:")
    print(b)
    xstar <- sample(x, replace=TRUE)
    dstar <- density(xstar, n=evalpoints, from=x_0, to=x_f)
    estimates[,b] <- dstar$y
  }
  CI_95 <- apply(estimates, 1, quantile, probs=c(0.05, 0.95))
  xshade <- c(xax, rev(xax))
  yshade <- c(CI_95[2,], rev(CI_95[1,]))
  newList <- list("density"=d, "CI_intervals"=CI_95, "xshade"=xshade,
                  "yshade"=yshade)
  return(newList)
}
args <- commandArgs(trailingOnly = TRUE)

df1 <- read.table(args[1])
df2 <- read.table(args[2])
title <- args[3]
red_label <- args[4]
blue_label <- args[5]


file_name <- args[6]

head(df1)
head(df2)

last_x <- floor(max(na.omit(df1[,2]), na.omit(df2[,2])))

outputs_df1 <- get_CI_shade(df1, iterations=10, x_0=0, x_f=last_x)

outputs_df2 <- get_CI_shade(df2, iterations=10, x_0=0, x_f=last_x)



png(paste(toString(file_name),".png",sep=''), width=20, height=20)
  plot(x=c(0,last_x), y=c(0,max(outputs_df1$density$y, outputs_df2$density$y)),
       type='n', xlab="Distance to Ca (Å)", ylab="Density", main = title)
  lines(x=outputs_df1$density$x, y=outputs_df1$density$y, col='blue', lw=2)
  lines(x=outputs_df2$density$x, y=outputs_df2$density$y, col='red', lw=2)
  legend(x="topright", c(red_label, blue_label), title="System",
         lty=c(1,1), lwd=c(2.5, 2.5), col = c('blue', 'red'))
  polygon(outputs_df1$xshade, outputs_df1$yshade, border=NA,
          col=adjustcolor("blue", 0.4))
  polygon(outputs_df2$xshade, outputs_df2$yshade, border=NA,
          col=adjustcolor("red", 0.4))
  axis(1, at=seq(0, last_x, 1))
dev.off()