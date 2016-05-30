#!/usr/bin/env Rscript


# run as ./script.R file1.dat file2.dat title_of_files red_label blue_label

# red_label should match file1.dat
# blue_label should match file2.dat

library(ggplot2)

rmsfs_averaging <- function(data){
  # Input is the rmsf raw data
  even <- seq_len(ncol(data)) %% 2
  resids <- data.frame(data[even])  # Select the columns with the residue numbers
  rmsfs <- data.frame(data[!even])  # Select the columns with the RMSFs values
  resids_avg <- rowMeans(resids)    # Residue mean just gives the same number
  
  rmsfs_avg <- rowMeans(rmsfs)
  rmsfs_sd <- apply(rmsfs, 1, sd)
  
  error <- qt(0.975, df = ncol(rmsfs)-1) * rmsfs_sd/sqrt(ncol(rmsfs))
  
  left <- rmsfs_avg - error
  right <- rmsfs_avg + error
  
  
  clean <- data.frame(cbind(resids_avg, rmsfs_avg, rmsfs_sd, left, right))
  colnames(clean) <- c("Residue", "Mean", "Stdev", "CI_left", "CI_right")
  return(clean)
}

rmsfs_return <- function(data) {
  # Input is the rmsf raw data
  even <- seq_len(ncol(data)) %% 2
  rmsfs <- data.frame(data[!even])  # Select the columns with the RMSFs values
  return(rmsfs)
}

rmsfs_ttests <- function(df1, df2) {
  # Input data-frames need to be only RMSF values (output of rmsfs_return)
  output <- data.frame(c(1:419))

  for (row in 1:nrow(df1)) {
    rowdf1 <- df1[row,]
    rowdf2 <- df2[row,]
    t_test <- t.test(rowdf1, rowdf2)  # Unpaired T-test, assuming unequal variances
    left <- t_test$conf.int[1]        # Lower bound of the 95% CI of the mean difference
    right <- t_test$conf.int[2]       # Upper bound of the 95% CI of the mean difference
    output[row,2] <- left
    output[row,3] <- right
    output[row,4] <- t_test$p.value
  }
  colnames(output) <- c("Residue", "Mean_Difference_Left", "Mean_Difference_Right", "pvalue")
  return(output)
}

pvalue_plotter <- function(df) {
  library(ggplot2)
  
  pvalues_count <- df$pvalue < 0.05  # Generates a TRUE or FALSE sequence
  trues <- table(pvalues_count)["TRUE"]
  falses <- table(pvalues_count)["FALSE"]
  total <- trues + falses
  proportion <- trues * 100 / total # Count the % of residues that have a p-value lower than 0.05
  if (is.na(proportion)) {
    proportion <- 0
  }

  specify_decimal <- function(x, k) format(round(x, k), nsmall=k)
  proportion <- specify_decimal(proportion,1)
  proportion_string <- toString(proportion)
  x_third_q <- quantile(df$Residue)["100%"]
  y_third_q <- quantile(df$pvalue)["100%"]
  
  last_res <- max(df$Residue)
  first_res <- min(df$Residue)
  
  axis_breaks <- scale_x_continuous(breaks = round(seq(from = first_res, to = last_res, length.out = 10)))
  
  plot <- ggplot(df, aes(x=Residue, y = pvalue)) +
          geom_point() +
          geom_line(aes(x=Residue, y = 0.05, color = 'red')) +
                        axis_breaks +
          theme_bw(15) +
          theme(legend.position="none") +
          annotate("text", x=x_third_q , y=0.75,
                   label = paste(proportion_string, "%")) +
          ylab("p-value") +
          scale_y_continuous(limits=c(0, 1))
  return(plot)
}

doplot <- function(data, title, red_label, blue_label){
  
  last_res <- max(data$Residue)
  first_res <- min(data$Residue)

  axis_breaks <- scale_x_continuous(breaks = round(seq(from = first_res, to = last_res, length.out = 10)))
  
  if ( last_res == 171 ) {
    p <- ggplot(data, aes(x = Residue, y = Mean, color = Group)) +
         geom_line(size=1.5) +
         scale_color_manual("System", labels = c(red_label, blue_label),
                            values = c("red", "blue")) +
         geom_ribbon(aes(ymin = CI_left, ymax = CI_right, colour = NULL,
                         fill = Group, alpha = 0.1)) +
         guides(alpha = F) +
         scale_fill_discrete("System", labels = c(red_label, blue_label)) +
         ylab("RMSF (Å)") + xlab("Residue number") +
         axis_breaks +
         theme_bw(15) +
         ggtitle(title) +
         theme(legend.justification=c(1,0), legend.position="bottom") +
         scale_y_continuous(limits=c(0,20))
  } else {
    p <- ggplot(data, aes(x = Residue, y = Mean, color = Group)) +
         geom_line(size=1.5) +
         scale_color_manual("System", labels = c(red_label, blue_label),
                            values = c("red", "blue")) +
         geom_ribbon(aes(ymin = CI_left, ymax = CI_right, colour = NULL,
                         fill = Group, alpha = 0.1)) +
         guides(alpha = F) +
         scale_fill_discrete("System", labels = c(red_label, blue_label)) +
         ylab("RMSF (Å)") + xlab("Residue number") +
         axis_breaks +
         theme_bw(15) +
         ggtitle(title) +
         theme(legend.position="none") +
         scale_y_continuous(limits=c(0,20))
  }
  return(p)
}


rmsf_comparison_plot <- function(d1, d2, red_label, blue_label) {
  # Columns of d1 and d2 must be
  # Residue Mean Stdev CI_left CI_right 
  TnC <- rbind(data.frame(d1[1:161,], Group = "d1"),
               data.frame(d2[1:161,], Group = "d2"))

  TnT <- rbind(data.frame(d1[162:248,], Group = "d1"),
               data.frame(d2[162:248,], Group = "d2"))
  TnT$Residue <- TnT$Residue + 50

  TnI <- rbind(data.frame(d1[249:419,], Group = "d1"),
               data.frame(d2[249:419,], Group = "d2"))
  TnI$Residue <- TnI$Residue - 248

  tnc_plot <- doplot(TnC, "cTnC", red_label, blue_label)
  tnt_plot <- doplot(TnT, "cTnT", red_label, blue_label)
  tni_plot <- doplot(TnI, "cTnI", red_label, blue_label)

  source("/Users/je714/Scripts/Rplots/Multiple_plot_function.R")
  return(multiplot(tnc_plot, tnt_plot, tni_plot))
}

pvalue_comparison_plot <- function(df) {
  # df must be the output of the T tests of two
  # rmsf data frames (rmsfs_ttests output)
  TnC <- df[1:161,]
  TnT <- df[162:248,]
  TnI <- df[249:419,]

  TnT$Residue <- TnT$Residue + 50
  TnI$Residue <- TnI$Residue - 248

  tnc_plot <- pvalue_plotter(TnC) + ggtitle("cTnC")
  tnt_plot <- pvalue_plotter(TnT) + ggtitle("cTnT")
  tni_plot <- pvalue_plotter(TnI) + ggtitle("cTnI")

  source("/Users/je714/Scripts/Rplots/Multiple_plot_function.R")
  return(multiplot(tnc_plot, tnt_plot, tni_plot))
}

args <- commandArgs(trailingOnly = TRUE)


print(args[1])
print(args[2])
print(args[3])
print(args[4])
print(args[5])

d1 <- read.table(args[1])
d2 <- read.table(args[2])

d1_clean <- rmsfs_averaging(d1)
d2_clean <- rmsfs_averaging(d2)

d1_rmsfs <- rmsfs_return(d1)
d2_rmsfs <- rmsfs_return(d2)

# Save RMSF comparison picture
png(paste(toString(args[3]),".png",sep=''), width=20, height=20*1.6180,
    units='cm', res=300)
rmsf_comparison_plot(d1_clean, d2_clean, red_label=args[4], blue_label=args[5])
dev.off()

# Save p value analysis picture
png(paste(toString(args[3]),"_pValues.png",sep=''), width=20, height=20*1.6180,
    units='cm', res=300)
pvalue_comparison_plot(rmsfs_ttests(d1_rmsfs, d2_rmsfs))
dev.off()

