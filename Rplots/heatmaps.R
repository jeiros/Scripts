#!/usr/bin/env Rscript


# run as ./script.R cmapfile1.dat cmapfile2.dat cmap1_label cmap2_label "Region To Plot"
library(ggplot2)
library(scales)



avg_contacts <- function(data){
  new_array <- c()
  for(column in 1:ncol(data)) {
    if(column%%4 == 3) {
      new_array <- c(new_array, (data[, column]/data[, column+ 1]))
    }
  }
  cosa <- matrix(new_array, nrow=nrow(data))
  medias <- rowMeans(cosa)
  new_data <- cbind(data[,1:2], medias)
  return(new_data)
}


heat_map_single <- function(data,Title="", xlab="", ylab="", stride_x=3, stride_y=1, min_contact, max_contact) {
  
  library(ggplot2)
  
  colnames(data) <- c("Res1", "Res2", "Contact") 

  if ((max(data$Res1) > 161) & (max(data$Res1) < 249)) {
    data$Res1 <- data$Res1 + 50
  } else if (max(data$Res1) >= 249) {
    data$Res1 <- data$Res1 - 248
  }
  if ((max(data$Res2) > 161) & (max(data$Res2) < 249)) {
    data$Res2 <- data$Res2 + 50
  } else if (max(data$Res2) >= 249) {
    data$Res2 <- data$Res2 - 248
  }

  p <- ggplot(data, aes(Res1, Res2)) + geom_tile(aes(fill=Contact)) +
              scale_x_continuous(breaks=seq.int(min(data$Res1),
                                 max(data$Res1), by=stride_x)) +
              scale_y_continuous(breaks=seq.int(min(data$Res2),
                                 max(data$Res2), by=stride_y)) +
              scale_fill_gradientn(colours=c("white", "red"),
                                   limits=c(min_contact, max_contact)) + 
              theme_classic(15) + labs(title=Title, x=xlab, y=ylab)
  return(p)
}


heat_map_diff <- function(data,Title="", xlab="", ylab="", stride_x=3, stride_y=1, min_contact, max_contact) {
  
  library(ggplot2)
  
  colnames(data) <- c("Res1", "Res2", "Contact") 

  if ((max(data$Res1) > 161) & (max(data$Res1) < 249)) {
    data$Res1 <- data$Res1 + 50
  } else if (max(data$Res1) >= 249) {
    data$Res1 <- data$Res1 - 248
  }
  if ((max(data$Res2) > 161) & (max(data$Res2) < 249)) {
    data$Res2 <- data$Res2 + 50
  } else if (max(data$Res2) >= 249) {
    data$Res2 <- data$Res2 - 248
  }

  plot <- ggplot(data, aes(Res1, Res2)) + geom_tile(aes(fill=Contact)) +
                 scale_x_continuous(breaks=seq.int(min(data$Res1),
                                    max(data$Res1), by=stride_x)) +
                 scale_y_continuous(breaks=seq.int(min(data$Res2),
                                    max(data$Res2), by=stride_y)) +
                 scale_fill_gradientn(colours=c("blue", "white", "red"),
                                      limits=c(min_contact, max_contact)) +
                 # gradientn is NOT a typo, don't change it!!
                 theme_classic(15) + labs(title=Title, x=xlab, y=ylab)
  return(plot)
}


difference <- function(df1, df2) {
  
  df1_avg <- avg_contacts(df1)
  df2_avg <- avg_contacts(df2)
  
  new_df <- cbind(df1_avg[,1:2], df2_avg[,3] - df1_avg[,3])
  colnames(new_df) <- c("Res1", "Res2", "Contact_diff")
  return(new_df)
}


args <- commandArgs(trailingOnly=TRUE)

df1 <- read.table(args[1])
df2 <- read.table(args[2])

# Get maps
df1_avg <- avg_contacts(df1)
df2_avg <- avg_contacts(df2)
difference_map <- difference(df1=df1,df2=df2)

# Get mins/maxs
df1_max <- max(df1_avg$medias)
df1_min <- min(df1_avg$medias)
df2_max <- max(df2_avg$medias)
df2_min <- min(df2_avg$medias)
dff_max <- max(difference_map$Contact_diff)
dff_min <- min(difference_map$Contact_diff)
abs_diff=max(abs(dff_max), abs(dff_min))

# Set the min/max values for the two types of plots
min_single <- min(df1_min, df2_min)
max_single <- max(df1_max, df2_max)
min_diff <- 0 - abs_diff
max_diff <- 0 + abs_diff

# Plot labels
df1_label <- args[3]
df2_label <- args[4]
region <- args[5]


if (toString(region) == "CcTnT-NcTnI") {
  xstride <- 1
  ystride <- 2
  xlabel <- "cTnT residue"
  ylabel <- "cTnI residue"
}
if (toString(region) == "CcTnT-inhibitorypeptide") {
  xstride <- 1
  ystride <- 1
  xlabel <- "cTnT residue"
  ylabel <- "cTnI residue (inhibitory peptide)"
}
if (toString(region) == "NcTnC-NcTnI") {
  xstride <- 3
  ystride <- 2
  xlabel <- "cTnC residue"
  ylabel <- "cTnI residue"
}
if (toString(region) == "NcTnC-switchpeptide") {
  xstride <- 3
  ystride <- 1
  xlabel <- "cTnC residue"
  ylabel <- "cTnI residue (switch peptide)"
}
if (toString(region) == "NcTnI-inhibitorypeptide") {
  xstride <- 3
  ystride <- 1
  xlabel <- "cTnI residue"
  ylabel <- "cTnI residue (inhibitory peptide)"
}
if (toString(region) == "cTnC-inhibitorypeptide") {
  xstride <- 10
  ystride <- 1
  xlabel <- "cTnC residue"
  ylabel <- "cTnI residue (inhibitory peptide)"
}


df1_p <- heat_map_single(data=df1_avg,
                         Title=paste(df1_label, toString(region), sep=" "),
                         xlab=xlabel, ylab=ylabel, stride_x=xstride,
                         stride_y=ystride, min_contact=min_single,
                         max_contact=max_single)

df2_p <- heat_map_single(data=df2_avg,
                         Title=paste(df2_label, toString(region), sep=" "),
                         xlab=xlabel, ylab=ylabel, stride_x=xstride,
                         stride_y=ystride, min_contact=min_single,
                         max_contact=max_single)

diff_p <-  heat_map_diff(data=difference(df1=df1,df2=df2),
                         Title=paste(paste(df1_label, df2_label, sep='-'),
                                     toString(region), sep=' '),
                         xlab=xlabel, ylab=ylabel, stride_x=xstride,
                         stride_y=ystride, min_contact=min_diff,
                         max_contact=max_diff)

ggsave(paste("df2_",toString(region),".eps",sep=''), plot=df2_p, dpi=900,
       width=10, height=10)
ggsave(paste("df1_",toString(region),".eps",sep=''), plot=df1_p, dpi=900,
       width=10, height=10)
ggsave(paste("diff_",toString(region),".eps",sep=''), plot=diff_p, dpi=900,
       width=10, height=10)










