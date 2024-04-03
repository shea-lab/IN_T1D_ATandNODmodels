# Gene expression analysis of RNAseq for T1D adoptive transfer (AT) cohort
# Analysis pipeline code written by Russell Urie and adapted by Jessica King

  sample_set <- "T1D Adoptive Transfer"
  currentDate <- Sys.Date()

# Set-up packages
library(pacman)
pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, 
               colorspace, ggplot2, biomaRt, fmsb, car, caret, mixOmics, DESeq2, apeglm, rgl, 
               devtools, qpcR, reshape2, gridExtra, factoextra, edgeR, cowplot,
               pheatmap, randomForest, ROCR, genefilter, Hmisc, rdist, factoextra,
               ggforce, glmnet, ggpubr, grid, ggvenn, ggridges, clusterProfiler,
		   GOSemSim, AnnotationHub, org.Mm.eg.db) 

# Part A: load and prepare data 
# (Sections 1-3) Raw and CPM data are read in & organized. Raw data is pseudo-
# normalized for subsequent analyses. All later sections rely on Part A.

# 1 Data Tables ----
  # * 1.1 Sample Names ----
  samples <- c("OVAd00_M1", "OVAd00_M2", "OVAd00_M4", "OVAd00_M5", 
               "BDCd00_M15", "BDCd00_M20", "BDCd00_M11", "BDCd00_M12", "BDCd00_M13", 
               "OVAd10_M1", "OVAd10_M2", "OVAd10_M4", "OVAd10_M5", 
               "BDCd10_M15", "BDCd10_M20", "BDCd10_M11", "BDCd10_M12", "BDCd10_M13")
  sample_num <- length(samples)

  # * 1.2 Raw Data ----
  r_file <- read.csv(".\\Data\\ATrawcounts.csv", header=T)
  genes <- r_file[,1] # Vector of gene names (rows: genes, columns: samples)
  r_counts <- r_file[,-1] # Remove first column of gene names temporarily
  colnames(r_counts) <- c(samples) # Add sample names as column names
  genes <- make.names(genes, unique = T) # Cant start w/ numbers as row names
  rownames(r_counts) <- genes # Add gene names as row names
  r_counts <- r_counts[rowSums(r_counts)>0,] # Remove 0 expression genes
  genes <- rownames(r_counts) # Non-zero genes
  
  # * 1.3 Normalized Data ----
  n_file <- read.csv(".\\Data\\ATCPM.csv", header=T)
  n_counts <- n_file[,-1] # Remove first column of gene names temporarily
  colnames(n_counts) <- c(samples) # Add sample names as column names
  n_counts <- n_counts[rowSums(n_counts)>0,] # Remove 0 expression genes
  rownames(n_counts) <- genes # Add gene names back as row names

# 2 Define Factors ----
cohort <- c(rep("Day 0",9), rep("OVA Day 10",4), rep("BDC Day 10",5))
peptide <- c(rep("Before",9), rep("OVA",4), rep("BDC",5))
time <- c(rep("Day00",9), rep("Day10",9))
day <- c(rep( 0, 9), rep( 10, 9))
mouse <- c(rep(c("M1", "M2", "M4", "M5", "M15", "M20", "M11", "M12", "M13"),2))
# Two ways to categorize samples
simple <- c(rep("Healthy", 13), rep("Diabetic", 5))
grouped <- c(rep("Healthy", 4), rep("Before Diabetes", 5), 
            rep("Healthy", 4), rep("Diabetic", 5))
full <- c(rep("OVADay00", 4), rep("BDCDay00", 5), rep("OVADay10", 4), 
          rep("BDCDay10 Pre-Diabetic", 2), rep("BDCDay10 Diabetic", 3))
clin_info <- cbind.data.frame( peptide, cohort, time, mouse, simple, grouped)
rownames(clin_info) <- samples
clin_info <- cbind(clin_info, c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1)) 
# 0=healthy, 1=diabetic
colnames(clin_info)[7] <- "state"
 
# B) REDUCE NOISE --------------------------------------------------------
# 3 Pseudo-normalize ----
ps_counts <- log(n_counts + 1, base = 2) # Transform to pseudo-normalized

# 4 Expression Filters ----
filt_low <- 2 # filters are per sample per gene
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
r_filt <- r_counts[g_filt,] # Used for DESeq2
n_filt <- n_counts[g_filt,] # Used for EN


# C) Analyze --------------------------------------------------------
 #5 Define groups
   OVA0 = n_filt[,c(1:4)]
   BDC0 = n_filt[,c(5:9)]
   OVA10 = n_filt[,c(10:13)]
   BDC10 = n_filt[,c(14:18)]
   OVAall = cbind(OVA0, OVA10)
   BDCall = cbind(BDC0, BDC10)
   Bothd0 = n_filt[ , time == "Day00"]
   Bothd10 = n_filt[ , time == "Day10"]
   Healthy = n_filt[ , simple == "Healthy"]
   dOVA = OVA10 - OVA0
   dBDC = BDC10 - BDC0

 #6 T-tests
   OVA10OVA0 <- pvalgenes( OVA10, OVA0)
   BDC10BDC0 <- pvalgenes( BDC10, BDC0)
   OVA0BDC0 <- pvalgenes( OVA0, BDC0)
   OVA10BDC10 <- pvalgenes( OVA10, BDC10)
   OVAallBDC10 <- pvalgenes( OVAall, BDC10)
   OVAallBDCall <- pvalgenes( OVAall, BDCall)
   both0both10 = pvalgenes( Bothd0, Bothd10)
   both0BDC10 = pvalgenes( Bothd0, BDC10)
   both0OVA10 = pvalgenes( Bothd0, OVA10)
   HealthyBDC10 = pvalgenes( Healthy, BDC10)
   delt <- pvalgenes( dBDC, dOVA)
   sig_g <- sig_pvaltable(list(OVA0, OVA10, BDC0, BDC10, OVAall, BDCall, Bothd0, 
	Bothd10, Healthy, dOVA, dBDC), c("OVA0", "OVA10", "BDC0", 
	"BDC10", "OVAall", "BDCall", "Bothd0", "Bothd10", "Healthy", "dOVA", "dBDC"))
   write.csv(sig_g, file=paste(currentDate, sample_set, filt_low, ",", filt_high, "filt",
                              "_all_t-test_p-value.csv", sep=""))

 #7 Specify genesets
  include <- unique(c(BDC10BDC0, OVAallBDC10, OVAallBDCall, both0BDC10, HealthyBDC10, delt))
  exclude <- unique(c(OVA10OVA0, OVA0BDC0, both0OVA10))
  n_addfilt = n_counts[setdiff(include, exclude),]
  g_addfilt <- rownames(n_addfilt) # Vector of filtered genes
  r_addfilt <- r_counts[g_addfilt,] # Used for DESeq2

# 8 Elastic net ----
#Factors for EN
groupedf <-factor(cohort, labels = c(1,2,0)) # Factor(s) for Elastic Net (EN) later
simplef <-factor(simple, labels = c(1,0))

# * 8.1 Grouped EN by cohort----

#Colors for grouped
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
colSide <- clin_info$grouped
for (x in 1:sample_num) {
  if (cohort[x] == "Day 0") {colSide[x] <- "skyblue" 
  } else if (cohort[x] == "OVA Day 10") {colSide[x] <- "dodgerblue"  
  } else if (cohort[x] == "BDC Day 10") {colSide[x] <- "red4"
  }}

names(colSide) = cohort
ecolor = c("red4", "skyblue", "dodgerblue")

#Before t-tests
  xfactors <- cbind.data.frame(groupedf) 
samplesbc = colnames(n_filt)
  rownames(xfactors) <- samplesbc
  dataset <- t(n_filt) # Load the data. 
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$groupedf %>% createDataPartition(p = 1.0, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(groupedf~., train.data)[,-1] # Predictor variables
  y <- train.data$groupedf # Outcome variable
  genes_EN <- multinom3_EN( x, y, 100) # number following multinom is number of groups
genes_EN_cohortAT = genes_EN
datainput = n_filt
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

#Figure 3A
  ggbiplot(prcomp(t(n_EN.95), scale.=T), ellipse=T, groups=names(colSide), var.axes=F, 
           var.scale=1, circle=T) + theme_classic(base_size=15) +
    theme(legend.position= c(0.25, 0.85)) +
    #ggtitle("Pre t-tests EN alpha=.95 Geneset") + 
    geom_point(size=3, color=colSide) +
    scale_color_manual(name="Group", values = ecolor)

# * 8.2 Simple EN by diseased v healthy----

#Colors for grouped
colSide <- clin_info$simple
for (x in 1:sample_num) {
  if (simple[x] == "Healthy") {colSide[x] <- "skyblue" 
  } else if (simple[x] == "Diabetic") {colSide[x] <- "red4"  
  }}

names(colSide) = c(rep("Non-diabetic", 13), rep("Diabetic", 5))
ecolor = c("red4", "skyblue")

#Before t-tests
  xfactors <- cbind.data.frame(simplef) 
samplesbc = colnames(n_filt)
  rownames(xfactors) <- samplesbc
  dataset <- t(n_filt) 
  dataset <- as.data.frame(cbind(xfactors, dataset)) 
  dataset <- dataset[complete.cases(dataset), ] 
  set.seed(123) 
  training.samples <- dataset$simplef %>% createDataPartition(p = 1.0, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(simplef~., train.data)[,-1] 
  y <- train.data$simplef 
  genes_EN <- binom_EN( x, y, 100) 
genes_EN_simple = genes_EN
datainput = n_filt
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

#Figure 3C
 heatmap.2(as.matrix(n_EN.95), scale = "row", col = coul, key = T, 
            xlab = "", ylab="", labCol=FALSE, # density.info="none",
            margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")

#Figure 3B
  ggbiplot(prcomp(t(n_EN.95), scale.=T), ellipse=T, groups=names(colSide), var.axes=F, 
           var.scale=1, circle=T) + theme_classic(base_size=15) +  
    theme(legend.position= c(0.25, 0.85)) +
    geom_point(size=3, color=colSide) +
    scale_color_manual(name="Group", values = ecolor)

# 9 SVD ----
 
n_cent <- n_EN.95 - rowMeans(n_EN.95) # First, center data on genes
svd2 <- svd(t(n_cent)) # Apply SVD on transposed data
centr_healthy <- colMeans(rbind(svd2$u[clin_info$state == 0,]),dims=1)
euclid <- rbind(centr_healthy, svd2$u) # Euclidean Distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:18]
euc <- cbind.data.frame(samples, euclid_dist)
rownames(euc) <- c(samples)
euc_SVD <- cbind(euc, mouse, time, cohort, simple, grouped)

colSide <- clin_info$grouped
for (x in 1:sample_num) {
  if (cohort[x] == "Day 0") {colSide[x] <- "skyblue" 
  } else if (cohort[x] == "OVA Day 10") {colSide[x] <- "dodgerblue"  
  } else if (cohort[x] == "BDC Day 10") {colSide[x] <- "red4"
  }}

names(colSide) <- cohort

# 10 Random Forest ---- 
predictor_data <- t(n_EN.95) # Transpose data & assign genes as predictors.
target <- clin_info[,"state"] # Set variable to predict target (reject status)
target[target==0] <- "Healthy"
target[target==1] <- "Diabetic"
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
euclid_dist <- euclid_dist[1:18]
euc <- cbind.data.frame(samples, euclid_dist)
rownames(euc) <- c(samples)
euc_RF <- cbind(euc, mouse, time, cohort, simple, grouped)
predictions <- as.vector(rf_output$votes[,2]) # ROC Curve, D:H votes as prediction
pred <- prediction(predictions,target)
perf_AUC <- performance(pred,"auc") # First calculate the AUC value
AUC <- perf_AUC@y.values[[1]]
perf_ROC <- performance(pred,"tpr","fpr") # Plot the actual ROC curve
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=4, scientific=F)))
options(digits=2) # Vote Distributions

# 11  Scoring System ----
scores <- cbind.data.frame(euc_SVD$euclid_dist, euc_RF$euclid_dist)
colnames(scores) <- c("SVD", "RF")
rownames(scores) <- samples

scores$SVD <- (scores$SVD - min(scores$SVD))/(max(scores$SVD) - min(scores$SVD))
scores$RF <- (scores$RF - min(scores$RF))/(max(scores$RF) - min(scores$RF))
 
#For figure 4D
ATRF = scores$RF

#Figure 3D
par(mfrow=c(1,1), mar = c(4,4,1,1), oma = c(0.5,0.5,0.5,0.5)) # Adjust margins
plot(scores$SVD, scores$RF, col=colSide, pch=19, cex=2, 
     mgp=c(2.5,1,0), cex.lab=1.25,
     xaxs="i", yaxs="i",
     xlim=c(0,1.1), ylim=c(0,1.1),
     ylab="Random Forest prediction (prob.)", 
     xlab="Singular Value Decomposition (arb. units)") 
#Add Confidence Interval Ellipses
dataEllipse(scores$SVD[clin_info$cohort == "BDC Day 10"], 
            scores$RF[clin_info$cohort == "BDC Day 10"], fill=T, fill.alpha=0.075,
            levels=c(0.7), center.pch=FALSE, draw=TRUE, add=TRUE, segments=51, 
            robust=FALSE, plot.points=FALSE, col="red4", pch=1, lwd=1, lty=1)
dataEllipse(scores$SVD[clin_info$cohort == "OVA Day 10"], 
            scores$RF[clin_info$cohort == "OVA Day 10"], levels=c(0.7), 
            fill=T, fill.alpha=0.075, center.pch=FALSE, draw=TRUE, add=TRUE, segments=51, 
            robust=FALSE, plot.points = FALSE, col="dodgerblue", pch=1, lwd=1, lty=1)
dataEllipse(scores$SVD[clin_info$cohort == "Day 0"], 
            scores$RF[clin_info$cohort == "Day 0"], levels=c(0.7), 
            fill=T, fill.alpha=0.075, center.pch=FALSE, draw=TRUE, add=TRUE, segments=51, 
            robust=FALSE, plot.points = FALSE, col="skyblue", pch=1, lwd=1, lty=1)
legend(.1, .95, legend=c("BDC Day 10", "OVA Day 10", "Day 0"),
       fill=c("red4", "dodgerblue", "skyblue"), cex=1.3)

# D Diff expression ----------------------------------------------------------
# 12 DESeq2 ----
# DESeq2 uses raw counts; factors defined previously; rows:genes & cols:samples

  # * 12.1 DE Analysis ----
  factor_DEseq <- clin_info$simple
  coldata <- data.frame(cbind(factor_DEseq))
  row.names(coldata) <- samples
  dds <- DESeqDataSetFromMatrix(countData=r_addfilt, colData=coldata, 
                                design = ~ factor_DEseq)
  paste(nrow(dds), " genes input into DESeq2 from these filters", sep="")
  dds <- DESeq(dds)
  res <- results(dds) # Table of log2 fold changes, p-values, & p-adj values

  # * 12.2 Organize Results ----
  resultsordered <- res[order(res$padj),]   # Order results by p-adj
  sum(res$padj < 0.1, na.rm=T)   # How many genes have p-adj < 0.1
  resSig <- subset(resultsordered, padj < 0.1) # Subset by p-adj < 0.1
  
  dds_genes <- rownames(resSig) # Differentialy-expressed genes
  r_DEG <- r_counts[dds_genes,] 
  n_DEG <- n_counts[dds_genes,]

  # 12.3 GSEA ----

egoAT <- enrichGO(gene         = dds_genes,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                readable  = TRUE)
#Figure 5A
plot(barplot(egoAT ,showCategory = 10))

#Figure 5C
cnetplot(egoAT, node_label="category")

#Supplemental table 1
ATgo = as.data.frame(egoAT)
write.csv(ATgo,file="ATgoterms.csv", row.names = FALSE)
