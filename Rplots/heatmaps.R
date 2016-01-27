#!/usr/bin/env Rscript


# run as ./script.R cmapfile1.dat cmapfile2.dat "Region To Plot"


library(ggplot2)
library(scales)
args <- commandArgs(trailingOnly = TRUE)


if (toString(args[3]) == "CcTnT-NcTnI") {
  xstride = 1
  ystride = 2
  xlabel = "cTnT residue"
  ylabel = "cTnI residue"
}
if (toString(args[3]) == "CcTnT-inhibitory peptide") {
  xstride = 1
  ystride = 1
  xlabel = "cTnT residue"
  ylabel = "cTnI residue (inhibitory peptide)"
}
if (toString(args[3]) == "NcTnC-NcTnI") {
  xstride = 3
  ystride = 2
  xlabel = "cTnC residue"
  ylabel = "cTnI residue"
}
if (toString(args[3]) == "NcTnC-switch peptide") {
  xstride = 3
  ystride = 1
  xlabel = "cTnC residue"
  ylabel = "cTnI residue (switch peptide)"
}
if (toString(args[3]) == "NcTnI-inhibitory peptide") {
  xstride = 2
  ystride = 1
  xlabel = "cTnI residue"
  ylabel = "cTnI residue (inhibitory peptide)"
}
if (toString(args[3]) == "cTnC-inhibitory peptide") {
  xstride = 5
  ystride = 1
  xlabel = "cTnC residue"
  ylabel = "cTnI residue (inhibitory peptide)"
}





avg_contacts <- function(data){
  new_array <- c()
  for(column in 1:ncol(data)) {
    if(column%%4 == 3) {
      new_array <- c(new_array, (data[, column]/data[, column+1]))
    }
  }
  cosa <- matrix(new_array, nrow=nrow(data))
  medias <- rowMeans(cosa)
  new_data <- cbind(data[,1:2], medias)
  return(new_data)
}


heat_map_single <- function(data,Title = "", xlab = "", ylab = "", stride_x = 3, stride_y = 1) {
  
  library(ggplot2)
  
  colnames(data) <- c("Res1", "Res2", "Contact") 

  if ((max(data$Res1) > 161) & (max(data$Res1) < 249)) {
    data$Res1 <- data$Res1 - 161
  } else if (max(data$Res1) >= 249) {
    data$Res1 <- data$Res1 - 248
  }
  if ((max(data$Res2) > 161) & (max(data$Res2) < 249)) {
    data$Res2 <- data$Res2 - 161
  } else if (max(data$Res2) >= 249) {
    data$Res2 <- data$Res2 - 248
  }



  p <- ggplot(data, aes(Res1, Res2)) + geom_tile(aes(fill = Contact)) +
    scale_x_continuous(breaks = seq.int(min(data$Res1),max(data$Res1), by = stride_x)) +
    scale_y_continuous(breaks = seq.int(min(data$Res2),max(data$Res2), by = stride_y)) +
    scale_fill_gradientn(colours = c("white", "red"), limits=c(0,1)) + 
    theme_classic(15) +
    labs(title = Title, x = xlab, y = ylab)
  return(p)
}


heat_map_diff <- function(data,Title = "", xlab = "", ylab = "", stride_x = 3, stride_y = 1) {
  
  library(ggplot2)
  
  colnames(data) <- c("Res1", "Res2", "Contact") 

  if ((max(data$Res1) > 161) & (max(data$Res1) < 249)) {
    data$Res1 <- data$Res1 - 161
  } else if (max(data$Res1) >= 249) {
    data$Res1 <- data$Res1 - 248
  }
  if ((max(data$Res2) > 161) & (max(data$Res2) < 249)) {
    data$Res2 <- data$Res2 - 161
  } else if (max(data$Res2) >= 249) {
    data$Res2 <- data$Res2 - 248
  }



  p <- ggplot(data, aes(Res1, Res2)) + geom_tile(aes(fill = Contact)) +
    scale_x_continuous(breaks = seq.int(min(data$Res1),max(data$Res1), by = stride_x)) +
    scale_y_continuous(breaks = seq.int(min(data$Res2),max(data$Res2), by = stride_y)) +
    scale_fill_gradientn(colours = c("blue","white", "red"), limits=c(-1,1)) + 
    theme_classic(15) +
    labs(title = Title, x = xlab, y = ylab)
  return(p)
}


difference <- function(df1, df2) {
  
  df_WT <- avg_contacts(df1)
  df_P <- avg_contacts(df2)
  
  new_df <- cbind(df_WT[,1:2], df_P[,3] - df_WT[,3])
  colnames(new_df) <- c("Res1", "Res2", "Contact_diff")
  return(new_df)
}



wt <- read.table(args[1])
p <- read.table(args[2])

wt_avg <- avg_contacts(wt)
p_avg <- avg_contacts(p)


p_p <- heat_map_single(data=p_avg, Title=paste("SP23/SP24", toString(args[3]), sep=" "), xlab=xlabel, ylab=ylabel,
 stride_x=xstride, stride_y=ystride)
wt_p <- heat_map_single(data=wt_avg, Title=paste("WT", toString(args[3]), sep=" "), xlab=xlabel, ylab=ylabel,
 stride_x=xstride, stride_y=ystride)

diff_p <-  heat_map_diff(data=difference(df1=wt,df2=p), Title=paste("(SP23/SP24-WT)", " ", toString(args[3]),sep=""),
  xlab=xlabel, ylab=ylabel,
  stride_x=xstride, stride_y=ystride)

ggsave(paste("p_",toString(args[3]),".eps",sep=''), plot=p_p, dpi=900, width=10, height=10)
ggsave(paste("wt_",toString(args[3]),".eps",sep=''), plot=wt_p, dpi=900, width=10, height=10)
ggsave(paste("diff_",toString(args[3]),".eps",sep=''), plot=diff_p, dpi=900, width=10, height=10)










