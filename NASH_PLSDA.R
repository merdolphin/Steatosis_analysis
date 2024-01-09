# Load necessary libraries
library(caret)        # for data preprocessing
library(pls)          # for PLSDA
library(ggplot2)      # for visualization
library(ROCR)         # for ROC curve
library("guidedPLS")
library("mixOmics")
library(doMC)
library(foreach)
library(VennDiagram)
library(doParallel)
registerDoParallel(cores = 96)  #

set.seed(1209)
setwd("/home/lina/scratch/18.NASH/Rscripts/")
source("NASH_preprocessing_data.R")

metabolomics_data <-
  lapply(liver_meta.S, as.numeric) %>% as.data.frame()

metabolomics_data <-
  lapply(plasma_meta.S, as.numeric) %>% as.data.frame()

metabolomics_data <-
  lapply(mRNA_count.S, as.numeric) %>% as.data.frame()

metabolomics_data %>% head()
data_for_plsda <- metabolomics_data
# Convert the response variable to a factor


# Perform PLSDA
res.plsda <-
  plsda(
    X = as.matrix(metabolomics_data[, -1]),
    Y = metabolomics_data$Histo_STEATOSIS,
    ncomp = 3
  )

plotIndiv(
  res.plsda,
  ind.names = TRUE,
  comp = c(1, 2),
  ellipse = TRUE,
  style = "ggplot2",
  cex = 2,
  legend = TRUE
)


mRNA_count.S %>% dim()
# run pca method on data
pca.mRNA = pca(mRNA_count.S[,-1], ncomp = 10, center = TRUE)

pca.mRNA
plot(pca.mRNA)  # barplot of the eigenvalues (explained variance per component)
plotIndiv(
  pca.mRNA,
  group = mRNA_count.S$Histo_STEATOSIS,
  ind.names = FALSE,
  # plot the samples projected
  legend = TRUE,
  title = 'PCA on SRBCT, comp 1 - 2'
) # onto the PCA subspace


run_sPLSDA <- function(data) {
  
  X <- data[,-1]
  Y <- data$Histo_STEATOSIS
  srbct.splsda <-
    splsda(X, Y, ncomp = 15)  # set ncomp to 10 for performance assessment later
  # plot the samples projected onto the first two components of the PLS-DA subspace
  plotIndiv(
    srbct.splsda ,
    comp = 1:2,
    group = data$Histo_STEATOSIS,
    ind.names = FALSE,
    # colour points by class
    ellipse = TRUE,
    # include 95% confidence ellipse for each class
    legend = TRUE,
    title = '(a) PLSDA with confidence ellipses'
  )
  
  # use the max.dist measure to form decision boundaries between classes based on PLS-DA data
  background = background.predict(srbct.splsda, comp.predicted = 2, dist = "max.dist")

  # plot the samples projected onto the first two components of the PLS-DA subspace
  plotIndiv(
    srbct.splsda,
    comp = 1:2,
    group = data$Histo_STEATOSIS,
    ind.names = FALSE,
    # colour points by class
    background = background,
    # include prediction background for each class
    legend = TRUE,
    title = " (b) PLSDA with prediction background"
  )

  
  perf.splsda.srbct <-perf(
    srbct.splsda,
    validation = "Mfold",
    folds = 10, ## 10
    nrepeat = 100, ## 50-100
    # use repeated cross-validation
    progressBar = FALSE,
    auc = TRUE
  ) # include AUC values


  # plot the outcome of performance evaluation across all ten components
  plot(
    perf.splsda.srbct,
    col = color.mixo(5:7),
    sd = TRUE,
    legend.position = "horizontal"
  )

  perf.splsda.srbct$choice.ncomp # what is the optimal value of components according to perf()
  perf.splsda.srbct$choice.ncomp[1,1]
  # grid of possible keepX values that will be tested for each component
  list.keepX <- c(1:10,  seq(20, 300, 10))
  dim(plasma.meta.count)
 list.keepX
  perf.splsda.srbct$choice.ncomp[1,1] -> ncomp.1
  # undergo the tuning process to determine the optimal number of variables
  tune.splsda.srbct <-
    tune.splsda(
      X,
      Y,
      ncomp =    perf.splsda.srbct$choice.ncomp[1,1],
      # calculate for first 4 components
      validation = 'Mfold',
      folds = 10, ## 10
      nrepeat = 100, ## 50-100
      # use repeated cross-validation
      dist = 'max.dist',
      # use max.dist measure
      measure = "BER",
      # use balanced error rate of dist measure
      test.keepX = list.keepX, 
      cpus=96
    ) 

  tune.splsda.srbct -> plasma.tune.splsda.srbct
  plot(tune.splsda.srbct, col = color.jet(ncomp.1)) # plot output of variable number tuning
  tune.splsda.srbct$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda()
  tune.splsda.srbct$choice.keepX # what are the optimal values of variables according to tune.splsda()
  optimal.ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
  optimal.keepX <- tune.splsda.srbct$choice.keepX[1:optimal.ncomp]
  return(tune.splsda.srbct)
  }

## data with last column as factor
data <- liver_meta.S
data <- lapply(data, as.numeric) %>% as.data.frame()

pdf("liver_meta.pdf")

liver_meta.splsda <- run_sPLSDA(data)

dev.off()

data <- plasma_meta.S
data <- lapply(data, as.numeric) %>% as.data.frame()
pdf("plasma_meta.pdf")
plasma_plasma.splsda <- run_sPLSDA(data)
save(tune.splsda.srbct, file="plasma_meta_splisda.tunned.rda")
dev.off()

data <- mRNA_count.S
data <- lapply(data, as.numeric) %>% as.data.frame()
pdf("mRNA_meta.pdf")
save(tune.splsda.srbct, file = "liver_meta_splisda.tunned.rda")
mRNA_plasma.splsda <- run_sPLSDA(data)
dev.off()

final.splsda$loadings

final.splsda <- splsda(X, Y, 
                       ncomp = 10, 
                       keepX = optimal.keepX)

plotIndiv(final.splsda, comp = c(2,7), # plot samples from final model
          group = liver_meta.S$Histo_STEATOSIS, ind.names = FALSE, # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = ' (a) sPLS-DA on SRBCT, comp 1 & 2')

final.splsda

plotIndiv(final.splsda, comp = c(1,3), # plot samples from final model
          group = liver_meta.S$Histo_STEATOSIS, ind.names = FALSE,  # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = '(b) sPLS-DA on SRBCT, comp 1 & 3')


# set the styling of the legend to be homogeneous with previous plots
legend=list(legend = levels(Y), # set of classes
            col = unique(color.mixo(Y)), # set of colours
            title = "Steatosis stages", # legend title
            cex = 0.7) # legend size

# generate the CIM, using the legend and colouring rows by each sample's class
cim <- cim(final.splsda, row.sideColors = color.mixo(Y), 
           legend = legend)

final.splsda$X 


data_opls <- opls(t(mRNA.gene.count.1), Histo_STEATOSIS, predI = 1, orthoI = 2)
# Set the number of permutations
n_permutations <- 1000



q2_permutations <- numeric(n_permutations)

# Perform permutation test
for (i in 1:n_permutations) {
  # Permute the response variable
  permuted_response <- sample(response)
  
  # Fit OPLS-DA model with permuted response
  permuted_opls <- opls(permuted_response ~., data = data, predI = 1, orthoI = 2)
  
  # Calculate Q2 for the permuted model
  q2_permutations[i] <- summary(permuted_opls)$qa$q2[2]
}

# Calculate the Q2 value for the original model
q2_original <- summary(data_opls)$qa$q2[2]

# Calculate the p-value for the permutation test
p_value <- sum(q2_permutations >= q2_original) / n_permutations

# Print the p-value
cat("Permutation test p-value:", p_value, "\n")

library(VIM)
library(rgl) # Interactive 3D plot will load the rgl library.
library(doParallel)
library(future)
getDoParWorkers()
plan(multisession, workers = 124)
registerDoParallel(cores = 96)  # or the number of cores you want to use
# For example, using z-score standardization
set.seed(1209)

X <- t(mRNA.gene.count.1)
colnames(X) <- mRNA.gene.dict$Gene.stable.ID
colnames(X) %>% unique() %>% length()
mRNA.pca <- pca(X)
plotIndiv(mRNA.pca,
          group=mRNA_count.S$Histo_STEATOSIS, 
          ind.names = FALSE,
          legend=TRUE)
plotVar(mRNA.pca)

Y=mRNA_count.S$Histo_STEATOSIS
help(plsda)
plsda.srbct <- splsda(X, Y, ncomp=10)
set.seed(1209)
registerDoParallel(cores = 98)
perf.splsda.srbct <- perf(plsda.srbct, validation = 'Mfold', folds=10,
                          progressBar=TRUE,
                          nrepeat=50,
                          BPPARAM = MulticoreParam(verbose = TRUE)
)
plot(perf.splsda.srbct, sd = TRUE, legend.position = 'horizontal')

perf.splsda.srbct$choice.ncomp

register(MulticoreParam(workers = 98))
warnings()
grid.keepX <- c(1:10, seq(20,300,10))
grid.kset.seed(30) # For reproducibility with this handbook, remove otherwise
tune.spca.result <- tune.splsda(X, Y, ncomp = 8, 
                                validation = 'Mfold',
                                folds = 10, 
                                test.keepX = grid.keepX, 
                                nrepeat = 50,
                                progressBar = TRUE,
                                cpus = 118,
                                measure="BER") 
tune.splsda.srbct <- tune.spca.result
plot(tune.spca.result, col = color.jet(8)) # plot output of variable number tuning

tune.spca.result$choice.keepX
tune.splsda.srbct$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda()

save(tune.splsda.srbct, file="tune.splsda.mRNA.rda")
tune.splsda.srbct$choice.keepX # what are the optimal values of variables according to tune.splsda()
optimal.ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
optimal.ncomp <- 4
optimal.keepX <- tune.splsda.srbct$choice.keepX[1:optimal.ncomp]
optimal.keepX

# form final model with optimised values for component and variable count
final.splsda <- splsda(X, Y, 
                       ncomp = optimal.ncomp, 
                       keepX = optimal.keepX)

plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = mRNA_count.S$Histo_STEATOSIS, ind.names = FALSE, # colour by class label
          comp.legend = 1,  # specify the component for the legend
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = ' (a) sPLS-DA on SRBCT, comp 1 & 2')

plotIndiv(final.splsda, comp = c(1,4), # plot samples from final model
          group = mRNA_count.S$Histo_STEATOSIS, ind.names = FALSE,  # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = '(b) sPLS-DA on SRBCT, comp 1 & 3')
# set the styling of the legend to be homogeneous with previous plots

#set the styling of the legend to be homogeneous with previous plots
legend=list(legend = levels(Y), # set of classes
            col = unique(color.mixo(Y)), # set of colours
            title = "Steatosis", # legend title
            cex = 1) # legend size

# generate the CIM, using the legend and colouring rows by each sample's class
cim <- cim(final.splsda, row.sideColors = color.mixo(Y), 
           legend = legend)

dim(mRNA_count.S)
