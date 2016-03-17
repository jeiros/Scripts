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
  min_single = 0
  max_single = 0.17
  min_diff = -0.12
  max_diff = 0.12
}
if (toString(args[3]) == "CcTnT-inhibitorypeptide") {
  xstride = 1
  ystride = 1
  xlabel = "cTnT residue"
  ylabel = "cTnI residue (inhibitory peptide)"
  min_single = 0
  max_single = 0.14
  min_diff = -0.13
  max_diff = 0.13
}
if (toString(args[3]) == "NcTnC-NcTnI") {
  xstride = 3
  ystride = 2
  xlabel = "cTnC residue"
  ylabel = "cTnI residue"
  min_single = 0
  max_single = 0.7
  min_diff = -0.32
  max_diff = 0.32
}
if (toString(args[3]) == "NcTnC-switchpeptide") {
  xstride = 3
  ystride = 1
  xlabel = "cTnC residue"
  ylabel = "cTnI residue (switch peptide)"
  min_single = 0
  max_single = 0.83
  min_diff = -0.20
  max_diff =    0.20
}
if (toString(args[3]) == "NcTnI-inhibitorypeptide") {
  xstride = 3
  ystride = 1
  xlabel = "cTnI residue"
  ylabel = "cTnI residue (inhibitory peptide)"
  min_single = 0
  max_single = 0.06
  min_diff = -0.06
  max_diff = 0.06
}
if (toString(args[3]) == "cTnC-inhibitorypeptide") {
  xstride = 10
  ystride = 1
  xlabel = "cTnC residue"
  ylabel = "cTnI residue (inhibitory peptide)"
  min_single = 0
  max_single = 0.59
  min_diff = -0.3
  max_diff = +0.3
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


heat_map_single <- function(data,Title = "", xlab = "", ylab = "", stride_x = 3, stride_y = 1, min_contact, max_contact) {
  
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

  p <- ggplot(data, aes(Res1, Res2)) + geom_tile(aes(fill = Contact)) +
    scale_x_continuous(breaks = seq.int(min(data$Res1),max(data$Res1), by = stride_x)) +
    scale_y_continuous(breaks = seq.int(min(data$Res2),max(data$Res2), by = stride_y)) +
    scale_fill_gradientn(colours = c("white", "red"), limits=c(min_contact, max_contact)) + 
    theme_classic(15) +
    labs(title = Title, x = xlab, y = ylab)
  return(p)
}


heat_map_diff <- function(data,Title = "", xlab = "", ylab = "", stride_x = 3, stride_y = 1, min_contact, max_contact) {
  
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

  p <- ggplot(data, aes(Res1, Res2)) + geom_tile(aes(fill = Contact)) +
    scale_x_continuous(breaks = seq.int(min(data$Res1),max(data$Res1), by = stride_x)) +
    scale_y_continuous(breaks = seq.int(min(data$Res2),max(data$Res2), by = stride_y)) +
    scale_fill_gradientn(colours = c("blue","white", "red"), limits=c(min_contact,max_contact)) + 
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
 stride_x=xstride, stride_y=ystride, min_contact=min_single, max_contact=max_single)
wt_p <- heat_map_single(data=wt_avg, Title=paste("WT", toString(args[3]), sep=" "), xlab=xlabel, ylab=ylabel,
 stride_x=xstride, stride_y=ystride, min_contact=min_single, max_contact=max_single)

diff_p <-  heat_map_diff(data=difference(df1=wt,df2=p), Title=paste("(SP23/SP24-WT)", " ", toString(args[3]),sep=""),
  xlab=xlabel, ylab=ylabel,
  stride_x=xstride, stride_y=ystride, min_contact=min_diff, max_contact=max_diff)

ggsave(paste("p_",toString(args[3]),".eps",sep=''), plot=p_p, dpi=900, width=10, height=10)
ggsave(paste("wt_",toString(args[3]),".eps",sep=''), plot=wt_p, dpi=900, width=10, height=10)
ggsave(paste("diff_",toString(args[3]),".eps",sep=''), plot=diff_p, dpi=900, width=10, height=10)










