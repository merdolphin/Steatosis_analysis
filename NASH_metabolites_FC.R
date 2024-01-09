setwd("/home/lina/scratch/18.NASH/Rscripts/")
source("NASH_preprocessing_data.R")
setwd("/home/lina/scratch/18.NASH/Rscripts/")
source("NASH_preprocessing_data.R")
library(limma)
library(edgeR)
library(DESeq2)


# Create a DESeqDataSet object

samples_information <- read_xlsx("1.sample/HumanPatients_230216.xlsx")
samples_information %>% dplyr::select(ID, `Histo_STEATOSIS (0-3)`, Age, `Gender (0=female, 1=male, 2=see Key)`) %>% 
  setNames(c("ID", "Histo_STEATOSIS", "Age","Gender")) %>%
  mutate(Histo_STEATOSIS = str_replace_all(Histo_STEATOSIS, "3","2")) %>% 
  mutate(ID=as.character(ID)) %>%
  mutate(Gender = as.factor(Gender)) %>% 
  mutate(Histo_STEATOSIS = factor(Histo_STEATOSIS,  levels = c(0, 1, 2), labels = c("S0", "S1", "S2"))) -> sampled.ID.Steatosis.age.gender





############### liver metabolites #######
###########################################



sampled.ID.Steatosis.age.gender[sampled.ID.Steatosis.age.gender$ID %in% colnames(liver.meta.count),] -> sam94

sam94[match(colnames(liver.meta.count), sam94$ID),] -> sam94.1

liver.meta.count %>% colnames() == sam94.1$ID

aggr(sam94.1, plot=F)

sam94.1$Histo_STEATOSIS
liver.meta.count %>% colSums()
dge <- DGEList(counts = liver.meta.count, group = sam94.1$Histo_STEATOSIS)

bcv <- 0.2
et <- exactTest(dge, dispersion = bcv^2)
et$table %>% as.data.frame() %>% arrange(PValue)

dge
liver.meta.dict[1329,] 
liver.meta.dict[1320,]
liver.meta.dict[390,]
sam94.1
design <- model.matrix(~  Age + Gender + Histo_STEATOSIS, sam94.2)

sam94.1$Histo_STEATOSIS %>% summary()

# Perform common dispersion estimation
dge <- estimateCommonDisp(dge)

# Perform tagwise dispersion estimation
dge <- estimateTagwiseDisp(dge)

dge <- estimateDisp(dge, design)

dge
# Perform differential expression analysis
fit <- glmFit(dge, design)
fit <- glmQLFit(dge, design)
qlf.2vs0 <- glmQLFTest(fit, coef = 5) 
qlf.1vs0 <- glmQLFTest(fit, coef = 4) 
topTags(qlf.1vs0)
topTags(qlf.2vs0) %>% as.data.frame() %>% dplyr::filter(PValue < 0.05)
qlf.2vs1 <- glmQLFTest(fit, contrast = c(0,0,0,-1,1))
topTags(qlf.2vs1)

results <- topTags(fit, n=100)
dge 

help(makeContrasts)
# Contrast of interest, e.g., comparing 'Normal' to 'Steatosis'
contrast <- makeContrasts(S1-S0, S2-S1, S2-S0, levels = sam94.1$Histo_STEATOSIS)
# Set up the correct contrast matrix

# Fit the contrast
fit_contrast <- glmQLFit(dge, design)

# Extract results
results <- topTags(fit_contrast, n=1000)
results %>% tail()

# Filter for significantly different metabolites
significant_results <- results$table[results$table$FDR < 0.05, ]

significant_results
# View the top significant results
significant_results -> SigRes_S1vsS0

# Save the results to a CSV file
write.csv(significant_results, file = "significant_metabolites_edgeR.csv")

# Create MA plot
plotMA(fit_contrast)

# Create a volcano plot
plotSmear(fit_contrast, de.tags = 100)