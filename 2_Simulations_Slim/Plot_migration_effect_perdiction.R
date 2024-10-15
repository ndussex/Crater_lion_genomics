library(ggplot2)
library(doBy)
library(rio)
library(ggpubr)
library("RColorBrewer")
library("gridExtra")
library("grid")


# get names of all simulations
#allFiles <- Sys.glob(file.path(getwd(), "Migration_M*_K2=*.txt"))

setwd("~/Slim/Prediction/output")
allFiles <- list.files(pattern="Lion_M=*")

####made a mistake in the headers (BP1 and BP2 outputed twice) - so add extra header columns in bash

#(for i in Migration_*
#    do sed '1s/\avgNeutP1/avgNeutP1,xxx/g' $i > $i.temp 
#  mv $i.temp $i
#  done)
#
#(for i in Migration_*
#    do sed '1s/\avgNeutP2/avgNeutP2,yyy/g' $i > $i.temp 
#  mv $i.temp $i
#  done)

# loop over all simulations and save in a single data frame
# make an extra column for the pattern (e.g "M=0_KF=50")
df_all <- data.frame()
for (x in allFiles) {
  df_x <- read.table(x, sep=",", header=T)
  pattern <- unlist(strsplit(x, split="_"))
  df_x$pattern <- paste(pattern[10], sep="_") #, pattern[3],
  df_all <- rbind(df_all, df_x)
}

head(df_all)
tail(df_all)
#df_test <- read.table('Migration_M0.1_K2=50_9623.txt', sep=",", header=T)
#pattern<-unlist(strsplit('Migration_M0.1_K2=50_9623.txt', split="_"))
#df_test$pattern <- paste(pattern[2], pattern[3], sep="_")

## timing (call once)
start <- 410000  # end of 2020 (i.e. now-ish - last migrant from 2021)
end <- 410100 # until 2120

# loop over each unique pattern
# creates one data frame with pattern and het

results_all <- data.frame()
for (i in unique(df_all$pattern)){
  print(i)
  df_pattern <- df_all[which(df_all$pattern==i),]
  data_start <- df_pattern[which(df_pattern$Year==start),]
  data_end <- df_pattern[which(df_pattern$Year==end),] 
  #het <- (data_end$meanHetP2)/(data_start$meanHetP2)
  Pi<-(((data_end$PiP2)-(data_start$PiP2))/data_start$PiP2)*100
  het <- (((data_end$meanHetP2)-(data_start$meanHetP2))/data_start$meanHetP2)*100
  froh<- (((data_end$FROHP2)-(data_start$FROHP2))/data_start$FROHP2)*100
  avgvStrDel<-(((data_end$avgvStrDelP2)-(data_start$avgvStrDelP2))/data_start$avgvStrDelP2)*100
  avgStrDel <-(((data_end$avgStrDelP2)-(data_start$avgStrDelP2))/data_start$avgStrDelP2)
  avgModDel <-(((data_end$avgModDelP2)-(data_start$avgModDelP2))/data_start$avgModDelP2)
  avgWkDel  <-(((data_end$avgWkDelP2)-(data_start$avgWkDelP2))/data_start$avgWkDelP2)
  avgNeut  <-(((data_end$avgNeutP2)-(data_start$avgNeutP2))/data_start$avgNeutP2)*100
  Genload  <-(((1-data_end$FitnessP2)-(1-data_start$FitnessP2))/(1-data_start$FitnessP2))*100
  B         <-(((data_end$BP2)-(data_start$BP2))/data_start$BP2)*100
  K<- data_end$KCrater
  # create data frame with number of rows equal to number of elements in het
  results <- data.frame(matrix(nrow = length(het), ncol = 0))
  # fill data frame with het and pattern
  results$Pi <- Pi
  results$het <- het
  results$froh <- froh
  results$avgvStrDel <-avgvStrDel
  results$avgStrDel <- avgStrDel
  results$avgModDel <-avgModDel
  results$avgWkDel <-avgWkDel
  results$avgNeut <-avgNeut
  results$Genload <- Genload
  results$B <- B
  results$K <-K
  results$pattern <- i
  results_all <- rbind(results_all, results)
}


#Summary table for plotting
###########################

#1. Het and FROH
################

#Heterozygosity

#see below for summary function
tgc <- summarySE(results_all, measurevar="het", groupvars=c("K","pattern"))
tgc

### plot Works the best for now - version 3
p1<- ggplot(tgc, aes(x=factor(K), y=het,color=pattern)) + 
    geom_errorbar(aes(ymin=het-se, ymax=het+se), width=.2) +
    #geom_errorbar(aes(ymin=het-sd, ymax=het+sd), colour="red",width=.1) +
    geom_point(aes(colour = factor(pattern))) +
    geom_line(aes(group = pattern)) +
    geom_hline(yintercept=c(-5,0,5), linetype=c('dotted','dashed','dotted'), color='black') +
    scale_colour_manual(breaks = c("M=0", "M=1", "M=5", "M=10"), 
                        values = unique(as.character(color_codes))) +
  #labs(title = "", x = expression(K[Crater]), y = expression("% change in nucleotide diversity ("*pi*")")   ) +
  labs(title = "", x = expression(K[Crater]), y = "% change in Heterozygosity") +
  scale_y_continuous(limits = c(-25, 15), breaks = seq(-25, 15, 10))
 
p1b <-p1 + theme(axis.line = element_line(colour = "black"),
                 legend.title=element_blank(),
                 legend.position = "none",
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) 

#FROH
tgc2 <- summarySE(results_all, measurevar="froh", groupvars=c("K","pattern"))
tgc2

tgc2$pattern <- factor(tgc2$pattern, levels = c("M=0", "M=1", "M=5", "M=10"))

p2<- ggplot(tgc2, aes(x=factor(K), y=froh,color=pattern)) + 
  geom_errorbar(aes(ymin=froh-se, ymax=froh+se), width=.2) +
  geom_point(aes(colour = factor(pattern))) +
  geom_line(aes(group = pattern)) +
  geom_hline(yintercept=c(-5,0,5), linetype=c('dotted','dashed','dotted'), color='black') +
  scale_colour_manual(breaks = c("M=0", "M=1", "M=5", "M=10"), 
                      values = unique(as.character(color_codes))) +
  labs(title = "", x = expression(K[Crater]), y =  expression("% change in Inbreeding (F"["ROH"]*")"))  +
  scale_y_continuous(limits = c(-30, 70), breaks = seq(-30, 70, 10))
  
p2b<-p2 + theme(axis.line = element_line(colour = "black"),
                legend.title=element_blank(),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) 

#2. load (Inverse of fitness) vs Inbreeding Load (B)

#Genload
tgc3 <- summarySE(results_all, measurevar="Genload", groupvars=c("K","pattern"))
tgc3


### plot Works the best for now - version 3
p3<- ggplot(tgc3, aes(x=factor(K), y=Genload,color=pattern)) + 
  geom_errorbar(aes(ymin=Genload-se, ymax=Genload+se), width=.2) +
  #geom_errorbar(aes(ymin=het-sd, ymax=het+sd), colour="red",width=.1) +
  geom_point(aes(colour = factor(pattern))) +
  geom_line(aes(group = pattern)) +
  geom_hline(yintercept=c(-5,0,5), linetype=c('dotted','dashed','dotted'), color='black') +
  scale_colour_manual(breaks = c("M=0", "M=1", "M=5", "M=10"), 
                      values = unique(as.character(color_codes))) +
  #labs(title = "", x = "K", y = "Proportion Heterozygosity retained") 
  labs(title = "", x = expression(K[Crater]), y = "% change in Realised load") +
  scale_y_continuous(limits = c(-20, 35), breaks = seq(-20, 35, 10))

p3b <-p3 + theme(axis.line = element_line(colour = "black"),
                 legend.title=element_blank(),
                 legend.position = "none",
                 #panel.grid.major = element_blank(),
                 #panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) 

#Masked load
tgc4 <- summarySE(results_all, measurevar="B", groupvars=c("K","pattern"))
tgc4

tgc4$pattern <- factor(tgc4$pattern, levels = c("M=0", "M=1", "M=5", "M=10"))

p4<- ggplot(tgc4, aes(x=factor(K), y=B,color=pattern)) + 
  geom_errorbar(aes(ymin=B-se, ymax=B+se), width=.2) +
  #geom_errorbar(aes(ymin=het-sd, ymax=het+sd), colour="red",width=.1) +
  geom_point(aes(colour = factor(pattern))) +
  geom_line(aes(group = pattern)) +
  geom_hline(yintercept=c(-5,0,5), linetype=c('dotted','dashed','dotted'), color='black') +
  scale_colour_manual(breaks = c("M=0", "M=1", "M=5", "M=10"), 
                      values = unique(as.character(color_codes))) +
  labs(title = "", x = expression(K[Crater]), y = "% change in Masked load") +
  scale_y_continuous(limits = c(-25, 45), breaks = seq(-25, 45, 10))

p4b<-p4 + theme(axis.line = element_line(colour = "black"),
                legend.title=element_blank(),
                #panel.grid.major = element_blank(),
                #panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) 

#pdf("../Plots_mig_effect/Proportion_change_Pi_FROH_Loads.pdf", width=9, height=7)
#grid.arrange(p1b, p2b, p3b,p4b, nrow = 2, ncol = 2, widths=c(2,2.5))
#dev.off()


# Create a list of the plots
plots <- list(p1b, p2b, p3b, p4b)

# Create a list of labels
labels <- c("(a)", "(b)", "(c)", "(d)")

# Create grobs for the labels
label_grobs <- lapply(labels, function(label) {
  textGrob(label, x = unit(0.008, "npc"), y = unit(0.1, "npc"), 
           just = c("left", "top"), gp = gpar(fontsize = 12, fontface = "bold"))
})

pdf("../Plots_mig_effect/Proportion_change_Pi_FROH_Loads_5pct.pdf", width=9, height=8)
# Arrange the plots with labels outside the panels
grid.arrange(
  arrangeGrob(plots[[1]], top = label_grobs[[1]]),
  arrangeGrob(plots[[2]], top = label_grobs[[2]]),
  arrangeGrob(plots[[3]], top = label_grobs[[3]]),
  arrangeGrob(plots[[4]], top = label_grobs[[4]]),
  nrow = 2, ncol = 2, widths = c(2, 2.5)
)
dev.off()




#############################
### summary function(tgc) ###
#############################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which South handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

