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


heat_map <- function(data, Title = "", xlab = "", ylab = "", stride_x = 3, stride_y = 1) {
  
  library(ggplot2)
  
  colnames(data) <- c("Res1", "Res2", "Contact") 

  if ((max(data$Res1) > 161) & (max(data$Res1) < 249)) {
    data$Res1 <- data$Res1 - 162
  } else if (max(data$Res1) >= 249) {
    data$Res1 <- data$Res1 - 248
  }
  if ((max(data$Res2) > 161) & (max(data$Res2) < 249)) {
    data$Res2 <- data$Res2 - 162
  } else if (max(data$Res2) >= 249) {
    data$Res2 <- data$Res2 - 248
  }



  p <- ggplot(data, aes(Res1, Res2)) + geom_tile(aes(fill = Contact)) +
    scale_x_continuous(breaks = seq.int(min(data$Res1),max(data$Res1), by = stride_x)) +
    scale_y_continuous(breaks = seq.int(min(data$Res2),max(data$Res2), by = stride_y)) +
    scale_fill_gradient(low = "white", high = "steelblue") + theme_classic(15) +
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

