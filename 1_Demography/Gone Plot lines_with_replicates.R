#libraries
library(ggplot2)
library(doBy)
library(gridExtra)
library("ggpubr")


# Change directory
setwd("~/my_dir")

# Define a function to estimate mean and 95% CI for a vector
mean_ci <- function(x) {
  n <- length(x)
  mean_x <- mean(x)
  se_x <- sd(x)/sqrt(n)
  ci <- qt(0.975, n-1) * se_x
  return(c(mean_x, mean_x-ci, mean_x+ci))
}
####

# Load data

#the table contains 2 + N_replicates columns such as

#Gen years Pop1  Pop1  Pop1  Pop1  Pop1  Pop1  Pop1  Pop1  Pop1  Pop1
#1 5 2.50181 2.50035 2.50732 2.50219 2.50045 2.51021 2.50296 2.50072 2.50373 2.50199
#.
#.
#200 1000  9753.32 5700.14 15800.3 8810.43 23111.5 3153.46 12543.5 9242.67 14455.7 9731.95


res <- read.table("my_multiple_GONE_runs.txt",header=T)
tail(res)
head(res)

#make a df
df_res<-data.frame(res)

#compute means and CIs and make plot

pdf("my_multiple_GONE_runs.pdf", width=7, height=7)
# Subset the data frame to only include columns starting with 'MTR'
crater_cols <- grep("^Pop1", names(df_res))
crater_df <- df_res[, crater_cols]
# Apply the function to each row of the subsetted data frame
crater_row_means_ci <- t(apply(crater_df, 1, mean_ci))
head(crater_row_means_ci)
# Add row names and column names to the resulting matrix
rownames(crater_row_means_ci) <- rownames(crater_df)
crater_row_means_ci<-cbind(res$years,crater_row_means_ci)
colnames(crater_row_means_ci) <- c("years","mean", "lower", "upper")

plot(crater_row_means_ci[,'years'], crater_row_means_ci[,'lower'], main="Pop1", type = "l",col="darkred",lwd=0,xlim=c(10,1200),xlab="Time (Years BP)",ylab="Ne",xaxt='n',ylim=c(0,20000))
lines(crater_row_means_ci[,'years'], crater_row_means_ci[,'upper'], type = "l",col="darkred",lwd=0,xlim=c(10,1100),xlab="", xaxt='n')
# Shade the area between the curves
#my_color <- rgb(0.545, 0, 0, alpha = 0.5)
my_color <- rgb(1, 0.65, 0, alpha = 0.5)
polygon(c(crater_row_means_ci[,'years'], rev(crater_row_means_ci[,'years'])), c(crater_row_means_ci[,'lower'], rev(crater_row_means_ci[,'upper'])), col = my_color, border = NA)
lines(crater_row_means_ci[,'years'], crater_row_means_ci[,'mean'],  type = "l",col="black",lwd=2,xlim=c(10,1100),xlab="", xaxt='n')
ticks = c(0, 100, 200, 400, 600, 800, 1000, 1200)
axis(side = 1, at = ticks)
#abline(v=100, col="black",lwd=2, lty=2)
#abline(v=300, col="black",lwd=2, lty=2)
dev.off()

########### inverted

pdf("my_multiple_GONE_runs.pdf", width=7, height=7)
# Subset the data frame to only include columns starting with 'MTR'
crater_cols <- grep("^Pop1", names(df_res))
crater_df <- df_res[, crater_cols]
# Apply the function to each row of the subsetted data frame
crater_row_means_ci <- t(apply(crater_df, 1, mean_ci))
head(crater_row_means_ci)
# Add row names and column names to the resulting matrix
rownames(crater_row_means_ci) <- rownames(crater_df)
crater_row_means_ci<-cbind(res$years,crater_row_means_ci)
colnames(crater_row_means_ci) <- c("years","mean", "lower", "upper")

plot(crater_row_means_ci[,'years'], crater_row_means_ci[,'lower'], main="Pop1 ", type = "l",col="darkred",lwd=0,xlim=c(1200,10),xlab="Time (Years BP)",ylab="Ne",xaxt='n',ylim=c(0,20000))
lines(crater_row_means_ci[,'years'], crater_row_means_ci[,'upper'], type = "l",col="darkred",lwd=0,xlim=c(1100,10),xlab="", xaxt='n')
# Shade the area between the curves
#my_color <- rgb(0.545, 0, 0, alpha = 0.5)
my_color <- rgb(1, 0.65, 0, alpha = 0.5)
polygon(c(crater_row_means_ci[,'years'], rev(crater_row_means_ci[,'years'])), c(crater_row_means_ci[,'lower'], rev(crater_row_means_ci[,'upper'])), col = my_color, border = NA)
lines(crater_row_means_ci[,'years'], crater_row_means_ci[,'mean'],  type = "l",col="black",lwd=2,xlim=c(1100,10),xlab="", xaxt='n')
ticks = c(0, 100, 200, 400, 600, 800, 1000, 1200)
axis(side = 1, at = ticks)
#abline(v=100, col="black",lwd=2, lty=2)
#abline(v=300, col="black",lwd=2, lty=2)
dev.off()
