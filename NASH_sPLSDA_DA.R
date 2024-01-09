source("/home/lina/scratch/18.NASH/Rscripts/NASH_preprocessing_data.R")
liver.meta.count <- read.delim("3.metabolimic/FromChoi/Normalized_data/liver_data_290523.txt",
                               header=T, as.is=T)



######################## ##########
names(liver.meta.count) <- gsub(x=names(liver.meta.count), pattern="X", replacement = "")
liver.meta.count
liver.meta.count %>% filter(Batch == "HILIC_Neg") -> liver.HILIC_Neg

liver.HILIC_Neg %>% dplyr::select(all_of(RNA_plasma_liver.90)) -> liver.HILIC_Neg.90
rownames(liver.HILIC_Neg.90)

samples_information <- read_xlsx("1.sample/HumanPatients_230216.xlsx")
samples_information %>% dplyr::select(ID, `Histo_STEATOSIS (0-3)`) %>% 
  setNames(c("ID", "Histo_STEATOSIS")) %>%
  mutate(Histo_STEATOSIS = str_replace_all(Histo_STEATOSIS, "3","2")) %>% 
  mutate(ID=as.character(ID)) %>%
  mutate(Histo_STEATOSIS = as.factor(Histo_STEATOSIS)) -> sampled.ID.Steatosis

sampled.ID.Steatosis

liver.HILIC_Neg.90 %>% t() %>% as.data.frame() %>% rownames_to_column(var="ID.1") %>% 
  left_join(sampled.ID.Steatosis, by=c("ID.1" = "ID")) %>% dplyr::select(-ID.1) -> liver.HILIC_Neg.90.s


liver.HILIC_Neg.90.s -> data

X = data %>% dplyr::select(-Histo_STEATOSIS)
Y = data$Histo_STEATOSIS
pca.data <- pca(X, ncomp=10, center=TRUE)
plot(pca.data)

plotIndiv(pca.data, group = Y, ind.names = FALSE, # plot the samples projected
          legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2') 

srbct.splsda <- splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda , comp = 1:2, 
          group = Y, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')

# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(srbct.splsda, comp.predicted=2, dist = "max.dist")

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda, comp = 1:2,
          group = Y, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = " (b) PLSDA with prediction background")

# undergo performance evaluation in order to tune the number of components to use
perf.splsda.srbct <- perf(srbct.splsda, validation = "Mfold", 
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          progressBar = TRUE, 
                          auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.srbct, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")

perf.splsda.srbct$choice.ncomp # what is the optimal value of components according to perf()

# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 270, 10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 4, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 10, nrepeat = 50, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
                                 progressBar = TRUE,
                                 cpus = 124) # allow for paralleliation to decrease runtime


plot(tune.splsda.srbct, col = color.jet(4)) # plot output of variable number tuning

tune.splsda.srbct$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda()

tune.splsda.srbct$choice.keepX # what are the optimal values of variables according to tune.splsda()

optimal.ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
optimal.keepX <- tune.splsda.srbct$choice.keepX[1:optimal.ncomp]

save(optimal.ncomp, file = "liver_meta_HILIC_neg.optimal.ncomp.rda")
save(optimal.keepX, file = "liver_meta_HILIC_opttimal.keepX.rda")
# form final model with optimised values for component and variable count
final.splsda <- splsda(X, Y, 
                       ncomp = optimal.ncomp, 
                       keepX = optimal.keepX)

save(final.splsda, file = "liver.final.splsda.model.rda")

plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = Y, ind.names = FALSE, # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = ' (a) sPLS-DA on SRBCT, comp 1 & 2')

plotIndiv(final.splsda, comp = c(1,3), # plot samples from final model
          group = Y, ind.names = FALSE,  # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = '(b) sPLS-DA on SRBCT, comp 1 & 3')




#########################
###########
###########  HILIC postive
###########
####################
liver.meta.count %>% filter(Batch == "HILIC_Pos") -> liver.HILIC_Pos

liver.HILIC_Pos %>% dplyr::select(all_of(RNA_plasma_liver.90)) -> liver.HILIC_Pos.90


liver.HILIC_Pos.90 %>% t() %>% as.data.frame() %>% rownames_to_column(var="ID.1") %>% 
  left_join(sampled.ID.Steatosis, by=c("ID.1" = "ID")) %>% dplyr::select(-ID.1) -> liver.HILIC_Pos.90.s


liver.HILIC_Pos.90.s -> data

X = data %>% dplyr::select(-Histo_STEATOSIS)
Y = data$Histo_STEATOSIS
pca.data <- pca(X, ncomp=10, center=TRUE)
plot(pca.data)

plotIndiv(pca.data, group = Y, ind.names = FALSE, # plot the samples projected
          legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2') 

srbct.splsda <- splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda , comp = 1:2, 
          group = Y, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')

# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(srbct.splsda, comp.predicted=2, dist = "max.dist")

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda, comp = 1:2,
          group = Y, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = " (b) PLSDA with prediction background")

# undergo performance evaluation in order to tune the number of components to use
perf.splsda.srbct <- perf(srbct.splsda, validation = "Mfold", 
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          progressBar = TRUE, 
                          auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.srbct, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")

perf.splsda.srbct$choice.ncomp # what is the optimal value of components according to perf()

# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 270, 10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 4, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 10, nrepeat = 50, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
                                 progressBar = TRUE,
                                 cpus = 124) # allow for paralleliation to decrease runtime


plot(tune.splsda.srbct, col = color.jet(4)) # plot output of variable number tuning

tune.splsda.srbct$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda()

tune.splsda.srbct$choice.keepX # what are the optimal values of variables according to tune.splsda()
optimal.ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
optimal.keepX <- tune.splsda.srbct$choice.keepX[1:optimal.ncomp]

save(optimal.ncomp, file = "liver_meta_HILIC_pos.optimal.ncomp.rda")
save(optimal.keepX, file = "liver_meta_HILIC_pos_opttimal.keepX.rda")
# form final model with optimised values for component and variable count
final.splsda <- splsda(X, Y, 
                       ncomp = optimal.ncomp, 
                       keepX = optimal.keepX)

save(final.splsda, file = "liver.final.pos.splsda.model.rda")

plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = Y, ind.names = FALSE, # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = ' (a) sPLS-DA on SRBCT, comp 1 & 2')

plotIndiv(final.splsda, comp = c(1,3), # plot samples from final model
          group = Y, ind.names = FALSE,  # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = '(b) sPLS-DA on SRBCT, comp 1 & 3')



################## only lipid and lipid-like molecules #############3

liver.HILIC_Neg.90 %>% t() %>% as.data.frame() %>% rownames_to_column(var="ID.1") %>% 
  left_join(sampled.ID.Steatosis, by=c("ID.1" = "ID")) %>% dplyr::select(-ID.1) -> liver.HILIC_Neg.90.s


liver.HILIC_Neg.90.s -> data

X = data %>% dplyr::select(-Histo_STEATOSIS)
Y = data$Histo_STEATOSIS
pca.data <- pca(X, ncomp=10, center=TRUE)
plot(pca.data)

plotIndiv(pca.data, group = Y, ind.names = FALSE, # plot the samples projected
          legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2') 

srbct.splsda <- splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda , comp = 1:2, 
          group = Y, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')

# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(srbct.splsda, comp.predicted=2, dist = "max.dist")

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda, comp = 1:2,
          group = Y, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = " (b) PLSDA with prediction background")

# undergo performance evaluation in order to tune the number of components to use
perf.splsda.srbct <- perf(srbct.splsda, validation = "Mfold", 
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          progressBar = TRUE, 
                          auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.srbct, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")

perf.splsda.srbct$choice.ncomp # what is the optimal value of components according to perf()

# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 270, 10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 4, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 10, nrepeat = 50, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
                                 progressBar = TRUE,
                                 cpus = 124) # allow for paralleliation to decrease runtime


plot(tune.splsda.srbct, col = color.jet(4)) # plot output of variable number tuning

tune.splsda.srbct$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda()

tune.splsda.srbct$choice.keepX # what are the optimal values of variables according to tune.splsda()

optimal.ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
optimal.keepX <- tune.splsda.srbct$choice.keepX[1:optimal.ncomp]

save(optimal.ncomp, file = "liver_meta_HILIC_neg.optimal.ncomp.rda")
save(optimal.keepX, file = "liver_meta_HILIC_opttimal.keepX.rda")
# form final model with optimised values for component and variable count
final.splsda <- splsda(X, Y, 
                       ncomp = optimal.ncomp, 
                       keepX = optimal.keepX)

save(final.splsda, file = "liver.final.splsda.model.rda")

plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = Y, ind.names = FALSE, # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = ' (a) sPLS-DA on SRBCT, comp 1 & 2')

plotIndiv(final.splsda, comp = c(1,3), # plot samples from final model
          group = Y, ind.names = FALSE,  # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = '(b) sPLS-DA on SRBCT, comp 1 & 3')



################## only krust.test < 0.05 ones  #############3
setwd("/home/lina/scratch/18.NASH/Rscripts/")
load(file="liver.meta.AD0.05.rda")
liver.meta.AD0.05$key

liver_meta.S %>% dplyr::select(liver.meta.AD0.05$key) -> data

X = data 
Y = liver_meta.S$Histo_STEATOSIS
pca.data <- pca(X, ncomp=10, center=TRUE)
plot(pca.data)

plotIndiv(pca.data, group = Y, ind.names = FALSE, # plot the samples projected
          legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2') 

srbct.splsda <- splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda , comp = 1:2, 
          group = Y, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')

# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(srbct.splsda, comp.predicted=2, dist = "max.dist")

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda, comp = 1:2,
          group = Y, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = " (b) PLSDA with prediction background")

# undergo performance evaluation in order to tune the number of components to use
perf.splsda.srbct <- perf(srbct.splsda, validation = "Mfold", 
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          progressBar = TRUE, 
                          auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.srbct, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")

perf.splsda.srbct$choice.ncomp # what is the optimal value of components according to perf()

# grid of possible keepX values that will be tested for each component
list.keepX <- c(seq(80, 350, 10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 2, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 10, nrepeat = 50, # use repeated cross-validation
                                 #dist = 'max.dist', # use max.dist measure
                                 dist = "mahalanobis.dist",
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
                                 progressBar = TRUE,
                                 cpus = 124) # allow for paralleliation to decrease runtime


plot(tune.splsda.srbct, col = color.jet(2)) # plot output of variable number tuning

tune.splsda.srbct$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda()

tune.splsda.srbct$choice.keepX # what are the optimal values of variables according to tune.splsda()

optimal.ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
optimal.keepX <- tune.splsda.srbct$choice.keepX[1:optimal.ncomp]

save(optimal.ncomp, file = "liver_meta_DA0.05.optimal.ncomp.rda")
save(optimal.keepX, file = "liver_meta_DA0.05_opttimal.keepX.rda")
# form final model with optimised values for component and variable count
final.splsda <- splsda(X, Y, 
                       ncomp = optimal.ncomp, 
                       keepX = optimal.keepX)

save(final.splsda, file = "liver.final.splsda.model.rda")

plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = Y, ind.names = FALSE, # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = ' (a) sPLS-DA on SRBCT, comp 1 & 2')

BiocManager::install("ggdist")
library(ggdist)  # For the ellipse_geoms function

# Access the data for the plot
X_scores <- final.splsda$variates$X[, c(1, 2)]
Y_labels <- final.splsda$Y

# Create a data frame for ggplot2
plot_data <- data.frame(X1 = X_scores[, 1], X2 = X_scores[, 2], Group = as.factor(Y_labels))

# Plot using ggplot2
plot_plsda <- ggplot(plot_data, aes(x = X1, y = X2, color = Group, shape = Group)) +
  geom_point(size = 4) + stat_ellipse(aes(group = Group)) +
  scale_color_manual(values = c("orange", "blue", "gray"), name = "Steatosis") + 
  scale_shape_manual(values = c(17, 19,21), name = "Steatosis")+
  labs(title = "sPLS-DA model based on\n significantly different liver metabolites",
    x = "Comp 1: 29% expl. var",
       y = "Comp 2: 7% expl. var") +
  clean_background +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        legend.background = element_blank()) 
#  scale_color_discrete(labels=c("S0","S1","S2")) 

plot_plsda

pdf("try.pdf")
plot_plsda
dev.off()

# set the styling of the legend to be homogeneous with previous plots
legend=list(legend = levels(Y), # set of classes
            col = unique(color.mixo(Y)), # set of colours
            title = "Tumour Type", # legend title
            cex = 0.7) # legend size

# generate the CIM, using the legend and colouring rows by each sample's class
cim <- cim(final.splsda, row.sideColors = color.mixo(Y), 
           legend = legend)

# form new perf() object which utilises the final model
perf.splsda.srbct <- perf(final.splsda, 
                          folds = 10, nrepeat = 50, # use repeated cross-validation
                          validation = "Mfold", dist = "mahalanobis.dist",  # use max.dist measure
                          progressBar = TRUE)

# Extract classification rate
classification_rate <- 1 - perf.splsda.srbct$error.rate.all$overall$mahalanobis.dist

# Print classification rate
print(paste("Classification Rate:", classification_rate))

# plot the stability of each feature for the first three components, 'h' type refers to histogram
par(mfrow=c(1,2))
plot(perf.splsda.srbct$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.srbct$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)

perf.splsda.srbct$error.rate.all$overall$mahalanobis.dist


train <- sample(1:nrow(X), 60) # randomly select 50 samples in training
test <- setdiff(1:nrow(X), train) # rest is part of the test set
dim(X) 
90*0.7
# store matrices into training and test set:
X.train <- X[train, ]
X.test <- X[test,]
Y.train <- Y[train]
Y.test <- Y[test]
train.splsda.srbct <- splsda(X.train, Y.train, ncomp = optimal.ncomp, keepX = optimal.keepX)
# use the model on the Xtest set
predict.splsda.srbct <- predict(train.splsda.srbct, X.test, 
                                dist = "mahalanobis.dist")

# evaluate the prediction accuracy for the first two components
predict.comp2 <- predict.splsda.srbct$class$mahalanobis.dist[,2]

table(factor(predict.comp2, levels = levels(Y)), Y.test)
auc.splsda.1 = auroc(final.splsda, roc.comp = 1, print = TRUE) # AUROC for the first component
auc.splsda.2 = auroc(final.splsda, roc.comp = 2, print = FALSE) # AUROC for


auc.comp1 <- auc.splsda.1$graph.Comp1$data %>% ggplot(aes(x=Specificity, y=Sensitivity, color=Outcome)) +
  geom_line(size=1.4) +
  geom_abline(slope = 1, intercept = 0, linetype="dashed", color = "black", size=1.2) +
  labs(title ="sPLS-DA model predicted steatosis classification", x = "False positive rate (%)", y = "Sensitivity (%)") +
  clean_background + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        legend.position = c(0.75,0.25),       
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  scale_color_discrete(name="ROC curve using Comp. 1",
                       labels=c("S0 vs. Others: 0.95",
                                "S1 vs. Others: 0.49",
                                "S2 vs. Others: 0.98"))
auc.comp1
  

auc.comp2 <- auc.splsda.2$graph.Comp$data %>% ggplot(aes(x=Specificity, y=Sensitivity, color=Outcome)) +
  geom_line(size=1.4) +
  geom_abline(slope = 1, intercept = 0, linetype="dashed", color = "black", size=1.2) +
  labs(title = "sPLS-DA model predicted steatosis classification", x = "False positive rate (%)", y = "Sensitivity (%)") +
  clean_background +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        legend.position = c(0.75,0.25),       
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  scale_color_discrete(name="ROC curve using Comp. 2",
                       labels=c("S0 vs. Others: 0.97",
                                "S1 vs. Others: 0.86",
                                "S2 vs. Others: 0.98"))
auc.comp2

clean_background <- theme(plot.background = element_rect("white"),
                                 panel.background = element_rect("white"),
                                 panel.grid = element_line("white"),
                                 axis.line = element_line("black"),
                                 axis.text = element_text(size = 12, color = "black", face="bold"),
                                 axis.title = element_text(color = "black"),
                                 axis.text.x = element_text(size=12, face="bold"),
                                 axis.text.y = element_text(size=12, face="bold"),
                                 legend.text = element_text(size = 12),
                                 legend.key = element_rect("white"))


pdf("splsda_liver_meta_DA_auc.pdf", width=16, height=4)
gridExtra::grid.arrange(plot_plsda, auc.comp1, auc.comp2, ncol = 3)
dev.off()


