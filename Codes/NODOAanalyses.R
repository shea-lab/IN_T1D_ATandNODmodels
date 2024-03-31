# Gene expression analysis of OpenArray data for NOD model
# Base code written by Russell Urie and adapted by Jessica King
#Adapted from Aaron's OA analysis matlab code and Russ' RNAseq pipeline

  sample_set <- "T1D Open Array"
  currentDate <- Sys.Date()

#Load packages
library(pacman)
  pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, 
               colorspace, ggplot2, biomaRt, fmsb, car, caret, mixOmics, DESeq2, apeglm, rgl, 
               devtools, qpcR, reshape2, gridExtra, factoextra, edgeR, cowplot,
               pheatmap, randomForest, ROCR, genefilter, Hmisc, rdist, factoextra,
               ggforce, glmnet, ggpubr, grid, ggvenn, ggridges)

# Part A: load and prepare data 
# 1 Data Tables ----
 r_file = read.csv(".\\Data\\OAdata.csv", header=T)
 genes = r_file[,2]
 r_counts = r_file[,-c(1:3)] #remove labeling columns
 samples = c("t1_MZ4", "t1_DS3", "t1_DS7", "t1_DS8", "t1_DS9", "t2_MZ4", "t2_DS3",
		"t2_DS7", "t2_DS8", "t2_DS9", "t3_MZ4", "t3_DS3", "t3_DS7",
		 "t3_DS8", "t3_DS9")
 sample_num <- length(samples)
 colnames(r_counts) = c(samples) # Add sample names as column names
 genes = make.names(genes, unique=T)
 rownames(r_counts) = genes # Add gene names as row names
 r_counts[is.na(r_counts)] <- 0 # Convert NA to 0
 r_counts = r_counts[rowSums(r_counts)>0,] # Remove 0 expression genes
 genes = rownames(r_counts) # Non-zero genes
 n_counts = mapply('/',r_counts,r_counts["Actb", ]) #Normalize to houskeeping gene
 rownames(n_counts) = genes

# 2 Define Factors ----
time <- c(rep("t1",5), rep("t2",5), rep("t3",5))
mouse <- c(rep(c("MZ4", "DS3", "DS7", "DS8", "DS9"),3))
simple <- c(rep("Healthy", 10), rep("Diabetic", 5))
grouped = c(rep("t1",5), rep("t2",5), rep("t3",5))

# 3 Pseudo-normalize ----
ps_counts <- log(n_counts + 1, base = 2) # Transform to pseudo-normalized

# 4 Expression Filters ----
filt_low <- 0.5 # These filters are per sample per gene
  filt_high <- 16
keep1 <- rowSums(ps_counts) > (filt_low*sample_num) # Low count filter, required
ps_filt <- ps_counts[keep1,] # Apply filter to pseudonorm counts
dim(ps_filt) # Dimensions after applying low filter
 keep2 <- rowSums(ps_filt) < (filt_high*sample_num) # High count filter, optional
 ps_filt <- ps_filt[keep2,] # Apply filter to pseudonorm counts
 dim(ps_filt) # Dimensions after applying high filter
 ps_high <- ps_counts[rowSums(ps_counts) >= (filt_high*sample_num),]
 genes_high <- rownames(ps_high) # To check high filter genes for relevance
g_filt <- rownames(ps_filt) # Vector of filtered genes
n_filt <- n_counts[g_filt,] 

# C) Analyze --------------------------------------------------------
 #5 Define groups
   t1 = n_filt[,c(1:5)]
   t2 = n_filt[,c(6:10)]
   t3 = n_filt[,c(11:15)]
   Healthy = n_filt[ , simple == "Healthy"]
   Diabetic = n_filt[ , simple == "Diabetic"]

 #6 T-tests 
  t1t2 = pvalgenes( t1, t2)
  t1t3 = pvalgenes( t1, t3)
  t2t3 = pvalgenes( t2, t3) #essentially constant
  HallD = pvalgenes( Healthy, Diabetic) #essentially constant
   sig_g <- sig_pvaltable(list(t1, t2, t3, Healthy, Diabetic), c("t1", "t2", 
	"t3", "Healthy", "Diabetic"))
   write.csv(sig_g, file=paste(currentDate, sample_set, filt_low, ",", filt_high, "filt",
                              "_all_t-test_p-value.csv", sep=""))

 #7 Specify genesets
  include <- unique(c(t1t3, t2t3, HallD))
  exclude <- t1t2
  n_addfilt = n_counts[setdiff(include, exclude),]

# 8 Elastic net: Simple EN by diseased v healthy----
#Factors for EN
simplef <-factor(simple, labels = c(1,0))

#Colors for simple
colSide <- simple
for (x in 1:sample_num) {
  if (simple[x] == "Healthy") {colSide[x] <- "skyblue" 
  } else if (simple[x] == "Diabetic") {colSide[x] <- "red4"  
  }}

names(colSide) = c(rep("Healthy", 10), rep("Diabetic", 5))
ecolor = c("red4", "skyblue")

#After t-tests
  xfactors <- cbind.data.frame(simplef) # SELECT THE FACTOR YOU WANT!!!!!
samplesbc = colnames(n_addfilt)
  rownames(xfactors) <- samplesbc
  dataset <- t(n_addfilt) # Load the data. Should it be raw or normalized!?!?!?!
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first colmn
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$simplef %>% createDataPartition(p = 1.0, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(simplef~., train.data)[,-1] # Predictor variables; need reduced gene set to run, runs w/ 9814
  y <- train.data$simplef # Outcome variable
  genes_EN <- binom_EN( x, y, 100) # number following multinom is number of groups
genes_EN_simpleOApostt = genes_EN
datainput = n_addfilt
  g_EN1 <- unique(rownames(datainput)[genes_EN$ENgenes$`1`])
  g_EN.95 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.95`])
  g_EN.9 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.9`])
  g_EN.85 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.85`])
  g_EN.8 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.8`])
  g_EN.75 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.75`])
  g_EN.7 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.7`])
  g_EN.65 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.65`])
  g_EN.6 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.6`])
  g_EN.55 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.55`])
  g_EN.5 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.5`])
  g_EN.45 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.45`])
  g_EN.4 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.4`])
  g_EN.35 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.35`])
  g_EN.3 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.3`])
  g_EN.25 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.25`])
  g_EN.2 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.2`])
  g_EN.15 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.15`])
  g_EN.1 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.1`])
  g_EN.05 <- unique(rownames(datainput)[genes_EN$ENgenes$`0.0499999999999999`])
  g_EN0 <- unique(rownames(datainput)[genes_EN$ENgenes$`0`])
  n_EN1 <- na.omit(datainput[g_EN1, ])
  n_EN.95 <- na.omit(datainput[g_EN.95, ])
  n_EN.9 <- na.omit(datainput[g_EN.9, ])
  n_EN.85 <- na.omit(datainput[g_EN.85, ])
  n_EN.8 <- na.omit(datainput[g_EN.8, ])
  n_EN.75 <- na.omit(datainput[g_EN.75, ])
  n_EN.7 <- na.omit(datainput[g_EN.7, ])
  n_EN.65 <- na.omit(datainput[g_EN.65, ])
  n_EN.6 <- na.omit(datainput[g_EN.6, ])
  n_EN.55 <- na.omit(datainput[g_EN.55, ])
  n_EN.5 <- na.omit(datainput[g_EN.5, ])
  n_EN.45 <- na.omit(datainput[g_EN.45, ])
  n_EN.4 <- na.omit(datainput[g_EN.4, ])
  n_EN.35 <- na.omit(datainput[g_EN.35, ])
  n_EN.3 <- na.omit(datainput[g_EN.3, ])
  n_EN.25 <- na.omit(datainput[g_EN.25, ])
  n_EN.2 <- na.omit(datainput[g_EN.2, ])
  n_EN.15 <- na.omit(datainput[g_EN.15, ])
  n_EN.1 <- na.omit(datainput[g_EN.1, ])
  n_EN.05 <- na.omit(datainput[g_EN.05, ])
  n_EN0 <- na.omit(datainput[g_EN0, ])

#Figure S1D
  heatmap.2(as.matrix(n_EN.9), scale = "row", col = coul, key = T, 
            xlab = "", ylab="", labCol=FALSE, # density.info="none",
            margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")

#Figure S1C
  ggbiplot(prcomp(t(n_EN.9), scale.=T), ellipse=T, groups=names(colSide), var.axes=F, 
           var.scale=1, circle=T) + theme_classic(base_size=15) + 
    theme(legend.position= c(0.75, 0.85)) +
    #ggtitle("Post t-tests EN alpha=.9 Geneset") + 
    geom_point(size=3, color=colSide) +
    scale_color_manual(name="Group", values = ecolor)

#9 SVD on EN data
n_cent <- n_EN.9 - rowMeans(n_EN.9) # First, center data on genes
svd2 <- svd(t(n_cent)) # Apply SVD on transposed data
centr_healthy <- colMeans(rbind(svd2$u[simple == "Healthy",]),dims=1)
euclid <- rbind(centr_healthy, svd2$u) # Euclidean Distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:15]
euc <- cbind.data.frame(samples, euclid_dist)
rownames(euc) <- c(samples)
euc_SVD <- cbind(euc, mouse, time, simple, grouped)

#Colors for simple
colSide <- simple
for (x in 1:sample_num) {
  if (simple[x] == "Healthy") {colSide[x] <- "skyblue" 
  } else if (simple[x] == "Diabetic") {colSide[x] <- "red4"  
  }}

# 10 Random Forest ---- 
predictor_data <- t(n_EN.9) # Transpose data & assign genes as predictors.
target <- simple # Set variable to predict target (reject status)
target[target=="Healthy"] <- "Healthy"
target[target=="Diabetic"] <- "Diabetic"
target <- as.factor(target)
tmp <- as.vector(table(target))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=F)[1]]
sampsizes <- rep(min_size,num_classes)
rf_output <- randomForest(x=predictor_data, y=target, importance=T, ntree=10001, 
                          proximity=T, sampsize=sampsizes, na.action=na.omit)
rf_importances <- importance(rf_output, scale=F) # Importance for each variable
confusion <- rf_output$confusion # Calculates sensitivity, specificity, accuracy
sensitivity <- (confusion[2,2]/(confusion[2,2]+confusion[2,1]))*100
specificity <- (confusion[1,1]/(confusion[1,1]+confusion[1,2]))*100
overall_error <- rf_output$err.rate[length(rf_output$err.rate[,1]),1]*100
overall_accuracy <- 1-overall_error
class1_error <- paste(rownames(confusion)[1]," error rate= ",confusion[1,3], sep="")
class2_error <- paste(rownames(confusion)[2]," error rate= ",confusion[2,3], sep="")
overall_accuracy <- 100-overall_error
sens_out <- paste("sensitivity=",sensitivity, sep="")
spec_out <- paste("specificity=",specificity, sep="")
err_out <- paste("overall error rate=",overall_error,sep="")
acc_out <- paste("overall accuracy=",overall_accuracy,sep="")
misclass_1 <- paste(confusion[1,2], rownames(confusion)[1],"misclassified as", 
                    colnames(confusion)[2], sep=" ")
misclass_2 <- paste(confusion[2,1], rownames(confusion)[2],"misclassified as", 
                    colnames(confusion)[1], sep=" ")
confusion_out <- confusion[1:2,1:2]
confusion_out <- cbind(rownames(confusion_out), confusion_out)
p_predictors <- varImpPlot(rf_output, type=2, n.var=12, scale=F, # Top variables
                           main="Variable Importance (Gini) for EN predictors")
target_labels <- as.vector(target) # MDS Class Separation
target_labels[target_labels=="Healthy"] <- "H"
target_labels[target_labels=="Diabetic"] <- "D"
plot_MDS <- MDSplot(rf_output, target,k=2,xlab="",ylab="",pch=target_labels, 
                    palette=c("red", "blue"), main="MDS plot")
plot_MDS <- plot_MDS$points
plot_MDS <- as.data.frame(plot_MDS)
p_mds <- ggplot(plot_MDS, aes(x=`Dim 1`,y=`Dim 2`, color=target_labels)) + 
  geom_point(aes(shape=donor),size=3) + 
  geom_text(aes(label = time),nudge_x=0.03, nudge_y=-0.01)
centr_healthy <- colMeans(rbind(plot_MDS[clin_info$state == 0,]),dims=1)
euclid <- rbind(centr_healthy, plot_MDS) # Euclidean distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:15]
euc <- cbind.data.frame(samples, euclid_dist)
rownames(euc) <- c(samples)
euc_RF <- cbind(euc, mouse, time, simple, grouped)
predictions <- as.vector(rf_output$votes[,2]) # ROC Curve, D:H votes as prediction
pred <- prediction(predictions,target)
perf_AUC <- performance(pred,"auc") # First calculate the AUC value
AUC <- perf_AUC@y.values[[1]]
perf_ROC <- performance(pred,"tpr","fpr") # Plot the actual ROC curve

#Figure S1F
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=4, scientific=F)))
options(digits=2) # Vote Distributions

# 11 Scoring System ----
scores <- cbind.data.frame(euc_SVD$euclid_dist, euc_RF$euclid_dist)
colnames(scores) <- c("SVD", "RF")
rownames(scores) <- samples

scores$SVD <- (scores$SVD - min(scores$SVD))/(max(scores$SVD) - min(scores$SVD))
scores$RF <- (scores$RF - min(scores$RF))/(max(scores$RF) - min(scores$RF))

#Figure S1E
par(mfrow=c(1,1), mar = c(4,4,1,1), oma = c(0.5,0.5,0.5,0.5)) # Adjust margins
plot(scores$SVD, scores$RF, col=colSide, pch=19, cex=2, 
     mgp=c(2.5,1,0), cex.lab=1.25,
     xaxs="i", yaxs="i",
     xlim=c(0,1.1), ylim=c(0,1.1),
     ylab="Random Forest prediction (prob.)", 
     xlab="Singular Value Decomposition (arb. units)") 
#Add Confidence Interval Ellipses
dataEllipse(scores$SVD[simple == "Diabetic"], 
            scores$RF[simple == "Diabetic"], fill=T, fill.alpha=0.075,
            levels=c(0.7), center.pch=FALSE, draw=TRUE, add=TRUE, segments=51, 
            robust=FALSE, plot.points=FALSE, col="red4", pch=1, lwd=1, lty=1)
dataEllipse(scores$SVD[simple == "Healthy"], 
            scores$RF[simple == "Healthy"], levels=c(0.7), 
            fill=T, fill.alpha=0.075, center.pch=FALSE, draw=TRUE, add=TRUE, segments=51, 
            robust=FALSE, plot.points = FALSE, col="skyblue", pch=1, lwd=1, lty=1)
legend(.1, .95, legend=c("Diabetic", "Healthy"),
       fill=c("red4", "skyblue"), cex=1.3)


