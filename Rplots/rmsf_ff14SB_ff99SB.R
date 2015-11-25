rmsfs_averaging <- function(data){

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

doplot <- function(data){
  
  last_res <- max(data$Residue)
  
  axis_breaks <- scale_x_continuous(breaks = round(seq(from = 1, to = last_res, length.out = 10)))
  
  plot <- ggplot(data, aes(x = Residue, y = Mean, color = Group)) + ylab("RMSF (Å)") + axis_breaks + 
    geom_ribbon(aes(ymin = CI_left, ymax = CI_right), alpha = 0.2) + geom_line() + xlab("Residue number")
  #return(plot)
  p <- ggplot(data, aes(x = Residue, y = Mean, color = Group)) + geom_line() +
    scale_color_manual("System", labels = c("S1P phosphorylation", "SEP phosphorylation"), values = c("red", "blue")) +
    geom_ribbon(aes(ymin = CI_left, ymax = CI_right, colour = NULL, fill = Group, alpha = 0.1)) +
    guides(alpha = F) + scale_fill_discrete("System", labels = c("S1P phosphorylation", "SEP phosphorylation")) + ylab("RMSF (Å)") + xlab("Residue number") + axis_breaks + theme_bw(15)
  return(p)
}

rmsfs_return <- function(data) {
  even <- seq_len(ncol(data)) %% 2
  resids <- data.frame(data[even])  # Select the columns with the residue numbers
  rmsfs <- data.frame(data[!even])  # Select the columns with the RMSFs values
  return(rmsfs)
}

rmsfs_ttests <- function(df1, df2) {
  # Input data-frames need to be only RMSF values
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
  specify_decimal <- function(x, k) format(round(x, k), nsmall=k)
  proportion <- specify_decimal(proportion,1)
  proportion_string <- toString(proportion)
  x_third_q <- quantile(df$Residue)["100%"]
  y_third_q <- quantile(df$pvalue)["100%"]
  
  last_res <- max(df$Residue)
  axis_breaks <- scale_x_continuous(breaks = round(seq(from = 1, to = last_res, length.out = 10)))
  plot <- ggplot(df, aes(x=Residue, y = pvalue)) + geom_point() + geom_line(aes(x=Residue, y = 0.05, color = 'red')) +
    axis_breaks + theme_classic(15) + theme(legend.position="none") + annotate("text", x=x_third_q , y=y_third_q, label = paste(proportion_string, "%")) + ylab("p-value")
  return(plot)
}

# ff99SB <- read.table("rmsf_99SB.dat")
# S1P <- read.table("all_rmsf_S1P.dat")
# SEP <- read.table("all_rmsf_SEP.dat")
# all_phos <- read.table("all_rmsf_anyphos.dat")
# View(all_phos)
# View(ff99SB)
# View(S1P)
# View(SEP)
# View(all_phos)
# source('~/Scripts/Rplots/rmsf_ff14SB_ff99SB.R')
# ff99SB_clean <- rmsfs_averaging(ff99SB)
# View(ff99SB_clean)
# S1P_clean <- rmsfs_averaging(S1P)
# SEP_clean <- rmsfs_averaging(SEP)
# all_phos_clean <- rmsfs_averaging(all_phos)
# View(all_phos_clean)
# ff99SB_rmsfs <- rmsfs_return(ff99SB)
# S1P_rmsfs <- rmsfs_return(S1P)
# SEP_rmsfs <- rmsfs_return(SEP)
# all_phos_rmsfs <- rmsfs_return(all_phos)