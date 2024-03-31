# Gene expression analysis of RNAseq from 6wk and diabetic timepoints in NOD cohort
# Base code written by Russell Urie and adapted by Jessica King

  sample_set <- "NOD 6wk vs diabetic"
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
  samples <- c("NODt0_M1", "NODt0_M2", "NODt0_M3", "NODt0_M4", "NODt0_M5", 
               "NODt4_M1", "NODt4_M2", "NODt4_M3", "NODt4_M4", "NODt4_M5")
  sample_num <- length(samples)

  # * 1.2 Raw Data ----
  r_file = read.csv(".\\Data\\NODseq.csv", header=T)
  r_file <- r_file[ grep("Gm", r_file[,3], invert = T) , ] # Remove Gm genes
  genes <- r_file[,3] # Vector of gene names (rows: genes, columns: samples)
  r_counts <- r_file[,-(1:4)] # Remove first column of gene names temporarily

#First explant
 t0 = r_counts[,c('X6872.JK.1','X6872.JK.2','X6872.JK.3','X6872.JK.4','X6872.JK.77')]
#Disease onset
 t4 = r_counts[,c('X6872.JK.21','X6872.JK.22','X6872.JK.23','X6872.JK.24','X6872.JK.25')]
#Combine groups
 r_counts = cbind(t0,t4)

  colnames(r_counts) <- c(samples) # Add sample names as column names
  genes <- make.names(genes, unique = T) 
  rownames(r_counts) <- genes # Add gene names as row names
  r_counts <- r_counts[rowSums(r_counts)>0,] # Remove 0 expression genes
  genes <- rownames(r_counts) # Non-zero genes

# 2 Define Factors ----
time <- c(rep("t0",5), rep("td",5))
mouse <- c(rep(c("M1", "M2", "M3", "M4", "M5"),2))
simple <- c(rep("Healthy", 5), rep("Diabetic", 5))
grouped = c(rep("t0",5), rep("td",5))

# 3 Pseudo-normalize ----
ps_counts <- log(r_counts + 1, base = 2) # Transform to pseudo-normalized

# 4 Expression Filters ----
filt_low <- 2 # These filters are per sample per gene
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
r_filt <- r_counts[g_filt,] 

# C) Analyze --------------------------------------------------------
 #5 Define groups
   t0 = r_filt[,c(1:5)]
   td = r_filt[,c(6:10)]
   Healthy = r_filt[ , simple == "Healthy"]
   Diabetic = r_filt[ , simple == "Diabetic"]

 #6 T-tests 
  t0td = pvalgenes( t0, td)
   sig_g <- sig_pvaltable(list(t0, td), c("t0", "td"))
   write.csv(sig_g, file=paste(currentDate, sample_set, filt_low, ",", filt_high, "filt",
                              "_all_t-test_p-value.csv", sep=""))

 #7 Specify genesets
  r_addfilt = r_counts[t0td,]

# 8 Elastic net: Simple EN by diseased v healthy ----
#Factors for EN
simplef <-factor(simple, labels = c(1,0))

#Colors for simple
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
colSide <- simple
for (x in 1:sample_num) {
  if (simple[x] == "Healthy") {colSide[x] <- "skyblue" 
  } else if (simple[x] == "Diabetic") {colSide[x] <- "red4"  
  }}

names(colSide) = c(rep("Healthy", 5), rep("Diabetic", 5))
ecolor = c("red4", "skyblue")

#Without t-tests
  xfactors <- cbind.data.frame(simplef) 
samplesbc = colnames(r_filt)
  rownames(xfactors) <- samplesbc
  dataset <- t(r_filt) # Load the data. 
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$simplef %>% createDataPartition(p = 1.0, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(simplef~., train.data)[,-1] # Predictor variables
  y <- train.data$simplef # Outcome variable
  genes_EN <- binom_EN( x, y, 100) 
genes_EN_simplet0tdpretnoGM = genes_EN
datainput = r_filt
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

#Figure 4B
 heatmap.2(as.matrix(n_EN.95), scale = "row", col = coul, key = T, 
            xlab = "", ylab="", labCol=FALSE, # density.info="none",
            margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")

#Figure 4A
  ggbiplot(prcomp(t(n_EN.95), scale.=T), ellipse=T, groups=names(colSide), var.axes=F, 
           var.scale=1, circle=T) + theme_classic(base_size=15) + 
    theme(legend.position= c(0.35, 0.85)) +
    geom_point(size=3, color=colSide) +
    scale_color_manual(name="Group", values = ecolor)


#9 SVD 
n_cent <- n_EN.95 - rowMeans(n_EN.95) # First, center data on genes
svd2 <- svd(t(n_cent)) # Apply SVD on transposed data
centr_healthy <- colMeans(rbind(svd2$u[simple == "Healthy",]),dims=1)
euclid <- rbind(centr_healthy, svd2$u) # Euclidean Distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:10]
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
predictor_data <- t(n_EN.95) # Transpose data & assign genes as predictors.
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
centr_healthy <- colMeans(rbind(plot_MDS[simple == "Healthy",]),dims=1)
euclid <- rbind(centr_healthy, plot_MDS) # Euclidean distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:10]
euc <- cbind.data.frame(samples, euclid_dist)
rownames(euc) <- c(samples)
euc_RF <- cbind(euc, mouse, time, simple, grouped)
predictions <- as.vector(rf_output$votes[,2]) # ROC Curve, D:H votes as prediction
pred <- prediction(predictions,target)
perf_AUC <- performance(pred,"auc") # First calculate the AUC value
AUC <- perf_AUC@y.values[[1]]
perf_ROC <- performance(pred,"tpr","fpr") # Plot the actual ROC curve
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=4, scientific=F)))
options(digits=2) # Vote Distributions

# 11 Scoring System ----
scores <- cbind.data.frame(euc_SVD$euclid_dist, euc_RF$euclid_dist)
colnames(scores) <- c("SVD", "RF")
rownames(scores) <- samples

scores$SVD <- (scores$SVD - min(scores$SVD))/(max(scores$SVD) - min(scores$SVD))
scores$RF <- (scores$RF - min(scores$RF))/(max(scores$RF) - min(scores$RF))

#Figure 4C
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

#For figure 4D
NODRF = scores$RF

# D Diff expression ----------------------------------------------------------
# 12 DESeq2 ----
# DESeq2 uses raw counts; factors defined previously; rows:genes & cols:samples

  # * 12.1 DE Analysis ----
  factor_DEseq <- c(simple)
  coldata <- data.frame(cbind(factor_DEseq))
  row.names(coldata) <- samples
  dds <- DESeqDataSetFromMatrix(countData=r_filt, colData=coldata, 
                                design = ~ factor_DEseq)
  paste(nrow(dds), " genes input into DESeq2 from these filters", sep="")
  dds <- DESeq(dds)
  res <- results(dds) # Table of log2 fold changes, p-values, & p-adj values

  # * 12.2 Organize Results ----
  resultsordered <- res[order(res$padj),]   # Order results by p-adj
  sum(res$padj < 0.1, na.rm=T)   # How many genes have p-adj < 0.1
  resSig <- subset(resultsordered, padj < 0.1) # Subset by p-adj < 0.1, 4365 genes
  
  dds_genes <- rownames(resSig) # Differentialy-expressed genes
  r_DEG <- r_counts[dds_genes,] 

  # 12.3 GSEA ----

egoNOD <- enrichGO(gene         = dds_genes,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                readable  = TRUE)

#Figure 5B
barplot(egoNOD, showCategory = 10)

#Figure 5D
cnetplot(egoNOD, node_label="category")

#Supplemental table 2
NODgo = as.data.frame(egoNOD)
write.csv(NODgo,file="NODgoterms.csv", row.names = FALSE)

#Overlapping pathways with adoptive transfer cohort.
#Run GSEA section of adoptive transfer analysis first before running this section.
bg = intersect(egoAT$ID, egoNOD$ID)

ATboth = egoAT[egoAT$ID %in% bg, ]
NODboth = egoNOD[egoNOD$ID %in% bg, ]

Supplemental tables 3 and 4
write.csv(ATboth,file="ATbothgoterms.csv", row.names = FALSE)
write.csv(NODboth,file="NODbothgoterms.csv", row.names = FALSE)

