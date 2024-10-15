library(ggplot2)
library(doBy)
library(rio)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(boot)
library(bootstrap)
library(dplyr)


#### colours ###################################

#install.packages("wesanderson")
#library(wesanderson)
#library(hrbrthemes) # for plot themes
#library(gapminder) # for data
#library(ggbump) # for the bump plot

# Change directory
setwd("~/my_dir_with_csv_files")
colors = brewer.pal(n = 8, name = "BrBG")[1:8]

#Version 3 (several bootstraps)
#main run (all chromosomes)
CV_Comp1<- read.csv("Data/results_Pop1_Pop2_boot_1.csv", sep=',',header=T)


#chatgtp code

# Get the list of files
file_list <- c(list.files(path = "Data", pattern = "results_Pop1_Pop2_boot_"))

# Create an empty list to store the data frames
file_data <- list()

# Read each file and store the data frame in the list
for (file in file_list) {
  file_data[[file]] <- read.csv(file.path("Data", file))
}
####

# Define a function to estimate the confidence interval for the last value
estimate_ci <- function(data) {
  ci <- boot::boot(data, statistic = function(x, i) tail(x[i], n = 1), R = 1000)$t
  quantile(ci, c(0.05, 0.95))
}

#### stats

# Create empty vectors to store the values
last_values <- vector()
mean_values <- vector()
median_values <- vector()
percentile_values <- vector()
min_values <- vector()
max_values <- vector()
ci_values <- matrix(ncol = 2, nrow = length(file_data))

# Iterate over the file_data list and calculate the statistics for each file
for (i in seq_along(file_data)) {
  last_value <- tail(file_data[[i]]$x, n = 1)
  last_values[i] <- last_value
  mean_values[i] <- mean(last_value)
  median_values[i] <- median(last_value)
  percentile_values[i] <- quantile(last_value, probs = 0.95)
  min_values[i] <- tail(file_data[[i]]$x, n = 1)
  max_values[i] <- tail(file_data[[i]]$x, n = 1)
  
  ci <- estimate_ci(file_data[[i]]$x)  # Estimate the confidence interval #seems wrong
  ci_values[i, ] <- ci
  
  
}


# Calculate overall statistics (values for all runs, even weird ones)
overall_mean <- mean(mean_values) 
overall_median <- median(median_values) 
percentile_up <- quantile(percentile_values, probs = 0.95) 
percentile_down <- quantile(percentile_values, probs = 0.05)
min <- min(min_values)
max <- max(max_values)
estimate_ci(last_values)
estimate_ci(median_values)

# Print the results
cat("Overall Mean:", overall_mean, "\n") #18381.1 
cat("Overall Median:", overall_median, "\n") # 
cat("95% Percentile Upper Bound:", percentile_up, "\n") # 
cat("95% Percentile Lower Bound:", percentile_down, "\n") #
cat("Overall Minimum:", min, "\n") #
cat("Overall Maximum:", max, "\n") #


# Create the plot
Myplot_boot <- ggplot()
colors = brewer.pal(n = 8, name = "BrBG")[1:8]
# Iterate over the file_data list and add each line to the plot
for (i in seq_along(file_data)) {
  Comp1_boot <- Comp1_boot + 
    geom_line(data=CV_Comp1, aes(x=x, y=y, group=label,color=label), size=1.5)+
    #geom_line(data=CV_Comp1_boot_3, aes(x=x, y=y, group=label,color=label), size=0.1)+
    geom_line(data = file_data[[i]], aes(x = x, y = y, group = label, color = label), size = 0.1)+
    labs(x = "Time (Years BP)", y =  "Ne") + 
    #scale_x_log10(breaks=c(1,10,100,1000,10000,100000), limits=c(10,100000),
                  #minor_breaks = c(1000,2000,3000,4000,5000,6000,7000,8000,9000),labels = c("1", "10","100", "1k", "10k", "100k")) +
    #scale_y_log10(breaks=c(1,10,100,1000,10000,100000), limits=c(10,100000)) +
    
    scale_x_log10(
      breaks = c(1, 10, 100, 1000, 10000, 100000),
      limits = c(10, 100000),
      minor_breaks = c(10^(1:5) * 2, 10^(1:5) * 5),
      labels = c("1", "10", "100", "1k", "10k", "100k")
    ) +
    scale_y_log10(
      breaks = c(1, 10, 100, 1000, 10000, 100000),
      limits = c(10, 100000),
      minor_breaks = c(10^(1:5) * 2, 10^(1:5) * 5),
      labels = c("1", "10", "100", "1k", "10k", "100k")
    ) +
    annotation_logticks(size = 0.3) +
    
    
    theme(axis.text.x = element_text(size=10),axis.text.y = element_text(size=10)) +
    scale_color_manual(values = c(colors[1],colors[3]),labels=c('Pop1', 'Pop2'))+
    theme(panel.background = element_blank(),axis.text=element_text(size=20),
          axis.title=element_text(size=14),
          legend.position = c(0.28, 0.18),
          legend.background = element_rect(fill="white"),
          legend.title = element_blank(), axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    
    geom_vline(xintercept=percentile_up ,linetype='dashed', color='black', size=0.2)+
    geom_vline(xintercept=percentile_down,linetype='dashed', color='black', size=0.2) +
  geom_vline(xintercept=overall_mean,linetype='dashed', color='black', size=0.5)
}

Myplot_boot


### Fig        #############################################################
pdf("My_Divcergence_plot.pdf", width=5, height=5)
grid.arrange(Myplot_boot)
dev.off()
############################################################################
