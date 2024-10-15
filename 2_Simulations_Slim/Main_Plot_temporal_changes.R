library("RColorBrewer")


#Temporal change in gen. variation through time: 1820 to 2120

# idea
# 1. import all runs (e.g. M0_K50, M0_K100,... M1_K400) as separate files

setwd("~/Slim/Hist_recap/output")
#loading files for each simulation set 
#only one dataset because the first part of the simulation is the same
#it only changes after 2020 with K and M

Lion_new_repro_M0_KF50  <- list.files(pattern="Lion_new_repro_5e09_murate_KGSE_5KSplit200yBP_long_sampling2yearsNEW_M=0_KF=50")

### !!!! note
# for 5000 genes, need to change FitnessP2 anf BP2 by a factor of 4

#FitnessP2 -> 1- (X^4)
#BP2 -> X * 4

#2. make dataframes for each dataset
df_Lion_new_repro_M0_KF50 <- data.frame(matrix(nrow = 0, ncol = 12))
for(rep in seq_along(Lion_new_repro_M0_KF50)){
  data <- read.table(Lion_new_repro_M0_KF50[rep], sep=",", header=T)
  df_Lion_new_repro_M0_KF50 <- rbind(df_Lion_new_repro_M0_KF50,cbind(data,rep))
}


#convert loads to whole genome for main figure
##############################################
#realised load
df_Lion_new_repro_M0_KF50$FitnessP1 <-1 - df_Lion_new_repro_M0_KF50$FitnessP1^4
df_Lion_new_repro_M0_KF50$FitnessP2 <-1 - df_Lion_new_repro_M0_KF50$FitnessP2^4

#masked load
df_Lion_new_repro_M0_KF50$BP1 <- (df_Lion_new_repro_M0_KF50$BP1)*4
df_Lion_new_repro_M0_KF50$BP2 <- (df_Lion_new_repro_M0_KF50$BP2)*4

## get mean of all reps(i.e. file) for each time point
# for temporal figure below (changes in indices through time)

#KF50
#M0
Year <- unique(df_Lion_new_repro_M0_KF50$Year)
means_Lion_new_repro_M0_KF50 <- data.frame(matrix(nrow = 0, ncol = 12))

for(Years in Year){
  data_this_Year <- df_Lion_new_repro_M0_KF50[which(df_Lion_new_repro_M0_KF50$Year==Years),]
  means_Lion_new_repro_M0_KF50 <- rbind(means_Lion_new_repro_M0_KF50,colMeans(data_this_Year))
  
}
colnames(means_Lion_new_repro_M0_KF50) <- colnames(df_Lion_new_repro_M0_KF50)

#stats
df_Lion_new_repro_M0_KF50[which(means_Lion_new_repro_M0_KF50$Year %in% c(309802, 310000)), ]
write.table(df_Lion_new_repro_M0_KF50[which(means_Lion_new_repro_M0_KF50$Year %in% c(309802, 310000)), ], "../Summary_time.txt", sep="\t", row.names=FALSE)

## timing
#300 burnin
start <- 309800 # (~1820, 200 y BP) = split
end <- 310000 #until 2020
####################################################################################################
# 1. plot temporal changes of indices 
####################################################################################################

######
xmin <- start
xmax <- end

#pdf(paste("../Plots/MAIN_Lion_new_repro_M0_KF50_N_FROH_load.pdf", sep=""), width=8, height=9)

pdf(paste("../Plots/MAIN_Lion_new_repro_M0_KF50_N_FROH_load_300Kburnin_samp2years_NEW.pdf", sep=""), width=8, height=9)

colors = brewer.pal(n = 8, name = "BrBG")[1:8]
lwd = 1
lab = 1.5
axis=1.2

par(mfrow=c(5,1), mar = c(2,5,0.5,0.5), oma=c(0,0,0,0)) #, bty = "n"
#first plot the average for N
plot(means_Lion_new_repro_M0_KF50$Year, means_Lion_new_repro_M0_KF50$popSizeP2, type = "l", ylab = "Population size", lwd = 1, cex.lab=lab, xlab = "",xlim = c(xmin+10,xmax-5),
     ylim=c((min(means_Lion_new_repro_M0_KF50$popSizeP2) - 20),max(means_Lion_new_repro_M0_KF50$popSizeP2) + 20), 
     col="black", yaxt = "n",xaxt = "n", cex.axis=axis) 

axis(2, at = seq(from = floor(0), to = ceiling(max(means_Lion_new_repro_M0_KF50$popSizeP2) + 10), by = 20), las = 1, cex.axis = axis)


abline(v = 309942, col = "black", lty = 3) 

#plot all runs
for(file in seq_along(Lion_new_repro_M0_KF50)){
  
  data <- df_Lion_new_repro_M0_KF50[which(df_Lion_new_repro_M0_KF50$rep==file),]
  lines(data$Year,data$popSizeP2, col="#FFA60080")
}
lines(means_Lion_new_repro_M0_KF50$Year,means_Lion_new_repro_M0_KF50$popSizeP2, lwd=2, col="black")

#first plot the average for het
plot(means_Lion_new_repro_M0_KF50$Year, means_Lion_new_repro_M0_KF50$meanHetP2, type = "l", ylab = "Heterozygosity", lwd = 1, cex.lab=lab, xlab = "",xlim = c(xmin+10,xmax-5),
     ylim=c((min(means_Lion_new_repro_M0_KF50$meanHetP2)),max(means_Lion_new_repro_M0_KF50$meanHetP2)+0.000005), 
     col="black", xaxt = "n", yaxt="n", cex.axis=axis) 
axis(2, at = axTicks(2), labels = formatC(axTicks(2), format = "e", digits = 1), cex.axis = axis)
abline(v = 309942, col = "black", lty = 3) 

#plot all runs for P2
for(file in seq_along(Lion_new_repro_M0_KF50)){
  
  data <- df_Lion_new_repro_M0_KF50[which(df_Lion_new_repro_M0_KF50$rep==file),]
  lines(data$Year,data$meanHetP2, col="#FFA60080")
}
lines(means_Lion_new_repro_M0_KF50$Year,means_Lion_new_repro_M0_KF50$meanHetP2, lwd=2, col="black")

# plot all mean for P1
lines(means_Lion_new_repro_M0_KF50$Year,means_Lion_new_repro_M0_KF50$meanHetP1, lwd=1, lty=2, col="black")

#first plot the average for Froh
plot(means_Lion_new_repro_M0_KF50$Year, means_Lion_new_repro_M0_KF50$FROHP2, type = "l", ylab = expression("Inbreeding (F"["ROH"]*")"), lwd = 1, cex.lab=lab, xlab = "",xlim = c(xmin+10,xmax-5),
     ylim=c(0,0.44), 
     col="black", xaxt = "n", cex.axis=axis) #ylim=c(0,550)
abline(v = 309942, col = "black", lty = 3) 
#plot all runs for P2
for(file in seq_along(Lion_new_repro_M0_KF50)){
  
  data <- df_Lion_new_repro_M0_KF50[which(df_Lion_new_repro_M0_KF50$rep==file),]
  lines(data$Year,data$FROHP2, col="#FFA60080")
}
lines(means_Lion_new_repro_M0_KF50$Year,means_Lion_new_repro_M0_KF50$FROHP2, lwd=2, col="black")

# plot all mean for P1
lines(means_Lion_new_repro_M0_KF50$Year,means_Lion_new_repro_M0_KF50$FROHP1, lwd=1, lty=2, col="black")


###FitnessP1 = mean genetic load (calculated multiplicatively across sites) from Kyriazis
plot(means_Lion_new_repro_M0_KF50$Year, means_Lion_new_repro_M0_KF50$FitnessP2, type = "l", ylab = "Realised load", lwd = 1, cex.lab=lab, xlab = "",xlim = c(xmin+10,xmax-5),
     ylim=c((min(means_Lion_new_repro_M0_KF50$FitnessP2-0.01)),max(means_Lion_new_repro_M0_KF50$FitnessP2+0.07)), 
     col="black", xaxt = "n", cex.axis=axis)
abline(v = 309942, col = "black", lty = 3) 
#plot all runs
for(file in seq_along(Lion_new_repro_M0_KF50)){
  
  data <- df_Lion_new_repro_M0_KF50[which(df_Lion_new_repro_M0_KF50$rep==file),]
  lines(data$Year,data$FitnessP2, col="#FFA60080")
}
lines(means_Lion_new_repro_M0_KF50$Year,means_Lion_new_repro_M0_KF50$FitnessP2, lwd=2, col="black")

# plot all mean for P1
lines(means_Lion_new_repro_M0_KF50$Year,means_Lion_new_repro_M0_KF50$FitnessP1, lwd=1, lty=2, col="black")


## masked load
plot(means_Lion_new_repro_M0_KF50$Year, means_Lion_new_repro_M0_KF50$BP2, type = "l", ylab = "Masked load", lwd = 1, cex.lab=lab, xlab = "Years (BP)",xlim = c(xmin+10,xmax-5),
     ylim=c((min(means_Lion_new_repro_M0_KF50$BP2-0.35)),max(means_Lion_new_repro_M0_KF50$BP2+2.5)), 
     col="black", xaxt = "n", cex.axis=axis)
abline(v = 309942, col = "black", lty = 3) 
#plot all runs
for(file in seq_along(Lion_new_repro_M0_KF50)){
  
  data <- df_Lion_new_repro_M0_KF50[which(df_Lion_new_repro_M0_KF50$rep==file),]
  lines(data$Year,data$BP2, col="#FFA60080")
}
lines(means_Lion_new_repro_M0_KF50$Year,means_Lion_new_repro_M0_KF50$BP2, lwd=2, col="black")

# plot mean for P1
lines(means_Lion_new_repro_M0_KF50$Year,means_Lion_new_repro_M0_KF50$BP1, lwd=1, lty=2, col="black")

axis(1, at=seq(xmin,xmax, by=10), labels=seq(1820, 2020, by=10), cex.axis=1)

dev.off()

################################################################
### 2. plot comparing allele counts before and after bottleneck
################################################################
#colors = brewer.pal(n = 8, name = "BrBG")[1:8]

## timing
t1950 <- 309930 #only need to call once
t1964 <- 309944 #only need to call once #post-epizootic 
t1972 <- 309950 #post epizootic
t2020 <- 310000 #only need to call once

###
#Lion_new_repro_KF50
Lion_new_repro_KF50_t1950 <- df_Lion_new_repro_M0_KF50[which(df_Lion_new_repro_M0_KF50$Year ==t1950) ,]
Lion_new_repro_KF50_t1964 <- df_Lion_new_repro_M0_KF50[which(df_Lion_new_repro_M0_KF50$Year ==t1964) ,]
Lion_new_repro_KF50_t1972 <- df_Lion_new_repro_M0_KF50[which(df_Lion_new_repro_M0_KF50$Year == t1972) ,]
Lion_new_repro_KF50_t2020 <- df_Lion_new_repro_M0_KF50[which(df_Lion_new_repro_M0_KF50$Year == t2020) ,]

pdf("../Plots/MAIN_Lion.pdf", width = 10, height = 14)
par(mfrow=c(2,2))

lab = 1.5
axis = 1
Title = 1.5

boxplot(Lion_new_repro_KF50_t2020$avgWkDelP1,Lion_new_repro_KF50_t1950$avgWkDelP2,Lion_new_repro_KF50_t1964$avgWkDelP2,Lion_new_repro_KF50_t1972$avgWkDelP2,Lion_new_repro_KF50_t2020$avgWkDelP2, col=c("#EEE8CD","#FFA60080","#FFA60080","#FFA60080","#FFA60080"), names = c("GSE 2020","C. 1950","C. 1964","C. 1972","C. 2020"), ylab="avg. N. alleles / ind.", main = "Weakly deleterious", cex.lab=lab, cex.axis=axis, cex.main=Title)
boxplot(Lion_new_repro_KF50_t2020$avgModDelP1,Lion_new_repro_KF50_t1950$avgModDelP2,Lion_new_repro_KF50_t1964$avgModDelP2,Lion_new_repro_KF50_t1972$avgModDelP2,Lion_new_repro_KF50_t2020$avgModDelP2, col=c("#EEE8CD","#FFA60080","#FFA60080","#FFA60080","#FFA60080"), names = c("GSE 2020","C. 1950","C. 1964","C. 1972","C. 2020"), ylab="", main = "Moderately deleterious", cex.lab=lab, cex.axis=axis, cex.main=Title)
boxplot(Lion_new_repro_KF50_t2020$avgStrDelP1,Lion_new_repro_KF50_t1950$avgStrDelP2,Lion_new_repro_KF50_t1964$avgStrDelP2,Lion_new_repro_KF50_t1972$avgStrDelP2,Lion_new_repro_KF50_t2020$avgStrDelP2, col=c("#EEE8CD","#FFA60080","#FFA60080","#FFA60080","#FFA60080"), names = c("GSE 2020","C. 1950","C. 1964","C. 1972","C. 2020"), ylab="avg. N. alleles / ind.", main = "Strongly deleterious deleterious", cex.lab=lab, cex.axis=axis, cex.main=Title)
boxplot(Lion_new_repro_KF50_t2020$avgvStrDelP1,Lion_new_repro_KF50_t1950$avgvStrDelP2,Lion_new_repro_KF50_t1964$avgvStrDelP2,Lion_new_repro_KF50_t1972$avgvStrDelP2,Lion_new_repro_KF50_t2020$avgvStrDelP2, col=c("#EEE8CD","#FFA60080","#FFA60080","#FFA60080","#FFA60080"), names = c("GSE 2020","C. 1950","C. 1964","C. 1972","C. 2020"), ylab="", main = "Very Strongly deleterious", cex.lab=lab, cex.axis=axis, cex.main=Title)
dev.off()
