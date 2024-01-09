source("~/scratch/18.NASH/Rscripts/NASH_preprocessing_data.R")
library(VIM)
library(mediation)
library(car)
library(mice)
library(ordinal)


samples_information <- read_xlsx("1.sample/HumanPatients_230216.xlsx")
samples_information %>% dplyr::select(ID,
                                      `Histo_STEATOSIS (0-3)`,  `Trig (<1.5)`, `LDL (<3.5 mmol/L)`,`HDL (>1.0 mmol/L = normal)`, 
                                      `TotalChol (3.5-5.5 mM)`, `NonHDL-Chol (mmol/L)`, `ALT(5-30`, `AST (10-35)`,
                                      `GGT (5-35 U/L)`, `Total Bili (3-15)`, `Ferritin (20-500 ng/mL)`,
                                      `HOMA2%B`,`HOMA2%S`,HOMA2IR,
                                      Age, `Gender (0=female, 1=male, 2=see Key)`,`DOS weight (kg)`,`DOS BMI`,
                                      `Diabetes (0=No, 1=Yes, 2= Ex diabetic)`
                                      ) %>% 
  setNames(c("ID", "Histo_STEATOSIS","Trig", "LDL", "HDL", "TotalChol", "NonHDLChol", "ALT", "AST", "GGT", "TotalBili", "Ferritin", 
             "HOMA2B", "HOMA2S", "HOMA2IR",
             "Age","Gender","Weight","BMI",
             "Diabetes")) %>%
  mutate(Histo_STEATOSIS = str_replace_all(Histo_STEATOSIS, "3","2")) %>% 
  mutate(ID=as.character(ID)) %>%
  mutate(Histo_STEATOSIS = as.factor(Histo_STEATOSIS)) -> sampled.ID.Steatosis



plasma.meta.1 %>% t() %>% as.data.frame() %>% rownames_to_column(var="ID.1") %>% 
  left_join(sampled.ID.Steatosis, by=c("ID.1" = "ID")) %>% dplyr::select(-ID.1) %>% dplyr::select(c(1088:1105)) %>%  
  mutate(HOMA2B = as.numeric(HOMA2B))  -> sam90.14


tempData <- mice(sam90.14, m=5, maxit = 50, meth="pmm", seed=500)

sam90.14.1 <- complete(tempData,1)

liver_meta.S %>% dim()
liver.meta.dict[c(37,59,89,90),] %>% dplyr::select(Compound)

names(mRNA_count.S) <- gsub(x=names(mRNA_count.S), pattern="V", replacement = "RNA")
names(liver_meta.S) <- gsub(x=names(liver_meta.S), pattern="V", replacement = "Liver")
names(plasma_meta.S) <- gsub(x=names(plasma_meta.S), pattern="V", replacement = "Plasma")


cbind(mRNA_count.S, liver_meta.S, plasma_meta.S) %>% dplyr::select(!contains("Histo_STEATOSIS")) %>% cbind(sam90.14.1) -> data

cbind(liver_meta.S, plasma_meta.S) %>% dplyr::select(!contains("Histo_STEATOSIS")) -> metabolites.data

regression_models <- list()

######################################################
###########
########## top association of mRNA with metabolites
##########
######################################################

mRNA_count.S[,colSums(mRNA_count.S[,1:60362]==0) < 18] -> mRNA_count.S.20per 



# Initialize an empty data frame to store results
pcor.test.res.mRNA.metabolites <- data.frame(mRNAs = character(),
                                 metabolites = character(),
                                 estimate = numeric(),
                                 pvalue = numeric(),
                                 stringsAsFactors = FALSE)



for (i in 10235:18159) {
  mRNA <- colnames(mRNA_count.S.20per)[i]
  for (j in 1:nrow(metabolites.data)) {
    metabolite <- colnames(metabolites.data)[j]
    
    # Perform Spearman correlation test
    tmp <- cor.test(
      mRNA_count.S.20per[[mRNA]] %>% as.numeric(),
      metabolites.data[[metabolite]],
      method = "spearman",
      exact = FALSE,
      na.action = na.omit
    )
    
    # Create a data frame for the current correlation test
    tmp.df <- data.frame(
       mRNAs=mRNA,
       metabolites = metabolite,
      estimate = as.numeric(tmp$estimate),
      pvalue = tmp$p.value,
      stringsAsFactors = FALSE
    )
    
    # Append the results to pcor.test.res
    pcor.test.res.mRNA.metabolites <- rbind(pcor.test.res.mRNA.metabolites, tmp.df)
  }
}


pcor.test.res.mRNA.metabolites %>% tail()

save(pcor.test.res.mRNA.metabolites, file="pcor.test.res.mRNA.metabolites.rda")

pcor.test.res.mRNA.metabolites %>% arrange(desc(abs(estimate)))
pcor.test.res.mRNA.metabolites %>% arrange(desc(abs(estimate))) %>% dplyr::filter(pvalue < 0.05 & abs(estimate) > 0.5) -> mRNA.trig.0.5

# Perform FDR correction on p-values
pcor.test.res.mRNA.metabolites$pvalue.adjusted <- p.adjust(pcor.test.res.mRNA.metabolites$pvalue, method = "BH")


pcor.test.res.mRNA.metabolites

# Filter significant results based on adjusted p-values
significant_results <- pcor.test.res.mRNA.metabolites %>% filter(pvalue.adjusted < 0.1)

significant_results %>% arrange(desc(abs(estimate))) %>% filter(abs(estimate) > 0.5) %>% filter(metabolites %in% c("Liver89","Liver90")) %>%
  dplyr::select(mRNAs) %>% distinct() -> liver.mRNAs.9

significant_results %>% arrange(desc(abs(estimate))) %>% filter(abs(estimate) > 0.5) %>% filter(metabolites %in% c("Liver37","Liver59","Liver89","Liver90")) %>%
  dplyr::select(mRNAs) %>% distinct() 


significant_results %>% arrange(desc(abs(estimate))) %>% filter(abs(estimate) > 0.5) %>% filter(metabolites %in% c("Liver37","Liver59","Liver89","Liver90")) 

# Display significant results
significant_results %>% arrange(desc(abs(estimate))) %>% filter(abs(estimate) > 0.5) %>% 
  dplyr::select(mRNAs) %>% unique() -> mRNAs.49

mRNAs.49

mRNA_count.S.20per %>% dplyr::select(mRNAs.49$mRNAs) -> mRNA_count.S.20per.49


significant_results %>% arrange(desc(abs(estimate))) %>% filter(abs(estimate) > 0.5) %>% 
  dplyr::select(metabolites) %>% unique() -> metabolites.11

metabolites.data %>% dplyr::select(metabolites.11$metabolites) -> metabolites.data.11

metabolites.data.11 %>% dim()



########################################
##############
############## metabolites vs. Clinical output
##############
###################

pcor.test.res.metabolite.11.ClinicalProfile <- data.frame(phenotypes = character(),
                                                 metabolites = character(),
                                                 estimate = numeric(),
                                                 pvalue = numeric(),
                                                 stringsAsFactors = FALSE)



for (i in 1:14){
 phenotype <- colnames(sam90.14)[i]
 for(j in 1:11){
   metabolite = colnames(metabolites.data.11)[j]
   # Perform Spearman correlation test
   tmp <- cor.test(
     sam90.14[[phenotype]] %>% as.numeric(),
     metabolites.data.11[[metabolite]],
     method = "spearman",
     exact = FALSE,
     na.action = na.omit
   )
   
   # Create a data frame for the current correlation test
   tmp.df <- data.frame(
    phenotypes = phenotype,
    metabolites = metabolite,
     estimate = as.numeric(tmp$estimate),
     pvalue = tmp$p.value,
     stringsAsFactors = FALSE
   )
   
   # Append the results to pcor.test.res
   pcor.test.res.metabolite.11.ClinicalProfile <- bind_rows(pcor.test.res.metabolite.11.ClinicalProfile, tmp.df)
 }
}

pcor.test.res.metabolite.11.ClinicalProfile %>% arrange(desc(abs(estimate))) %>% filter(abs(estimate) > 0.4) 
pcor.test.res.metabolite.11.ClinicalProfile %>% arrange(desc(abs(estimate))) %>% filter(abs(estimate) > 0.4) %>% dplyr::select(metabolites) %>%
  distinct() -> metabolites.4
 
metabolites.data %>% dplyr::select(metabolites.4$metabolites) -> metabolites.data.4



aggr(data, plot=F)
sam90.14.1 %>% dplyr::select(-Histo_STEATOSIS) -> sam90.14.1

sam90.14



sam90.14 %>% dplyr::select(Histo_STEATOSIS, ALT, AST) -> sam90.14.ALT.AST


# Assume your data frame is named 'data'
data <- data.frame(mRNA_count.S.20per.49, metabolites.data.4, sam90.14.1)

for (mRNA_col in colnames(mRNA_count.S.20per.49)){
  mRNA <- get(mRNA_col, data)
  model.test <- glm(Liver90 ~ mRNA + Age + Gender + Weight + BMI,  data=data)
  model.test %>% summary() %>% print()
}

# Initialize an empty data frame to store results
mediation_results <- data.frame(
  mRNA = character(),
  metabolite = character(),
  sam90_column = character(),
  indirect_effect = numeric(),
  total_effect = numeric(),
  stringsAsFactors = FALSE
)

liver.meta.dict[c(37,59,89,90),] %>% dplyr::select(Compound)

model.liver90 <- polr(Histo_STEATOSIS ~ Liver90, data=data)
model.liver89 <- clm(Histo_STEATOSIS ~ Liver89, data=data)
model.liver37 <- clm(Histo_STEATOSIS ~ Liver37, data=data)
model.liver59 <- clm(Histo_STEATOSIS ~ Liver59, data=data)

model.liver90 %>% summary()
sam90.14.1[,1:14]
# Create empty lists to store results
significant_mRNAs <- list()
significant_metabolites <- list()
significant_sam <- list()
significant_med_results <- list()
metabolites.4

for (mRNA_col in colnames(mRNA_count.S.20per.49)) {
  
    mRNA <- get(mRNA_col, data)
  
    
    model.mRNA.metabolite <-
      glm(metabolite ~ mRNA + Age + Gender + Weight + BMI, data = data)
    
    model.liver90 <-
      polr(Histo_STEATOSIS ~ Liver90 + mRNA + Age + Gender + Weight + BMI,
           data = data)

    med.result <- mediate(
      model.mRNA.metabolite,
      model.liver90,
      treat = "Gender",
      mediator = "Liver90",
      boot = TRUE
    )
    
    med.result %>% summary()
    med.result$d0.p
    
    # Check if p-values meet the criteria
    
    # Check if p-values meet the criteria
    if (any(med.result$d0.p < 0.05) &&
        any(med.result$d1.p < 0.05) &&
        any(med.result$z0.p < 0.05) && 
        any(med.result$z1.p < 0.05) ) {
      significant_mRNAs[[mRNA_col]] <- mRNA
      significant_metabolites[[metabolite_col]] <- metabolite 
      significant_sam[[sam_col]] <- sam
      significant_med_results[[mRNA_col]] <- med.result
    }
}

significant_mRNAs 

################################3

for (mRNA_col in colnames(mRNA_count.S.20per.49)) {
  
  mRNA <- get(mRNA_col, data)
  
  
  model.mRNA.metabolite <-
    glm(metabolite ~ mRNA + Age + Gender + Weight + BMI, data = data)
  
  model.liver89 <-
    polr(Histo_STEATOSIS ~ Liver89 + mRNA + Age + Gender + Weight + BMI,
         data = data)
  metabolite_col
  med.result <- mediate(
    model.mRNA.metabolite,
    model.liver89,
    treat = "Gender",
    mediator = "Liver89",
    boot = TRUE
  )
  
  med.result %>% summary()
  med.result
  
  # Check if p-values meet the criteria
  if (med.result$d0.p < 0.05 &&
      med.result$d1.p < 0.05 &&
      med.result$z0.p < 0.05 && med.result$z1.p < 0.05) {
    significant_mRNAs[[mRNA_col]] <- mRNA
    significant_metabolites[[metabolite_col]] <- metabolite 
    significant_med_results[[mRNA_col]] <- med.result
  }
}

############################
for (mRNA_col in colnames(mRNA_count.S.20per.49)) {
  
  mRNA <- get(mRNA_col, data)
  
  
  model.mRNA.metabolite <-
    glm(metabolite ~ mRNA + Age + Gender + Weight + BMI, data = data)
  
  model.liver59 <-
    polr(Histo_STEATOSIS ~ Liver59 + mRNA + Age + Gender + Weight + BMI,
         data = data)
  metabolite_col
  med.result <- mediate(
    model.mRNA.metabolite,
    model.liver59,
    treat = "Gender",
    mediator = "Liver59",
    boot = TRUE
  )
  
  med.result %>% summary()
  med.result
  
  # Check if p-values meet the criteria
  if (med.result$d0.p < 0.05 &&
      med.result$d1.p < 0.05 &&
      med.result$z0.p < 0.05 && med.result$z1.p < 0.05) {
    significant_mRNAs[[mRNA_col]] <- mRNA
    significant_metabolites[[metabolite_col]] <- metabolite 
    significant_med_results[[mRNA_col]] <- med.result
  }
}

##############################################
for (mRNA_col in colnames(mRNA_count.S.20per.49)) {
  
  mRNA <- get(mRNA_col, data)
  
  
  model.mRNA.metabolite <-
    glm(metabolite ~ mRNA + Age + Gender + Weight + BMI, data = data)
  
  model.liver37 <-
    polr(Histo_STEATOSIS ~ Liver37 + mRNA + Age + Gender + Weight + BMI,
         data = data)
  metabolite_col
  med.result <- mediate(
    model.mRNA.metabolite,
    model.liver37,
    treat = "Gender",
    mediator = "Liver37",
    boot = TRUE
  )
  
  med.result %>% summary()
  med.result
 
   
  # Check if p-values meet the criteria
  if (med.result$d0.p < 0.05 &&
      med.result$d1.p < 0.05 &&
      med.result$z0.p < 0.05 && med.result$z1.p < 0.05) {
    significant_mRNAs[[mRNA_col]] <- mRNA
    significant_metabolites[[metabolite_col]] <- metabolite 
    significant_med_results[[mRNA_col]] <- med.result
  }
}
########################################################


significant_mRNAs
# Display the results for significant mRNA variables
cat("Significant mRNA variables:\n", names(significant_mRNAs), "\n")

# Access the corresponding mediate results
significant_med_results[["RNA1693"]] %>% summary()
significant_med_results[["RNA12257"]] %>% summary()  
significant_med_results[["RNA7794"]] %>% summary()

mRNA.gene.dict[c(1693,12257,1197),]




########################################
##############
############## mRNA.49 not associatied with Clinical output
##############
###################

pcor.test.res.mRNA.49.ClinicalProfile <- data.frame(phenotypes = character(),
                                                          mRNAs = character(),
                                                          estimate = numeric(),
                                                          pvalue = numeric(),
                                                          stringsAsFactors = FALSE)



for (i in 1:14){
  phenotype <- colnames(sam90.14)[i]
  for(j in 1:49){
    mRNA = colnames(mRNA_count.S.20per.49)[j]
    # Perform Spearman correlation test
    tmp <- cor.test(
      sam90.14[[phenotype]] %>% as.numeric(),
      mRNA_count.S.20per.49[[mRNA]] %>% as.numeric(),
      method = "spearman",
      exact = FALSE,
      na.action = na.omit
    )
    
    # Create a data frame for the current correlation test
    tmp.df <- data.frame(
      phenotypes = phenotype,
      mRNAs = mRNA,
      estimate = as.numeric(tmp$estimate),
      pvalue = tmp$p.value,
      stringsAsFactors = FALSE
    )
    
    # Append the results to pcor.test.res
    pcor.test.res.mRNA.49.ClinicalProfile <- bind_rows(pcor.test.res.mRNA.49.ClinicalProfile, tmp.df)
  }
}

pcor.test.res.mRNA.49.ClinicalProfile  %>% filter(mRNAs %in% liver.mRNAs.9$mRNAs)




########################################
##############
############## mRNA vs. Clinical output
##############
###################
for (i in 1:18159) {
  mRNA <- colnames(mRNA.gene.count.trig.20per)[i]
  for (j in 18160:18171) {
    trig <- colnames(mRNA.gene.count.trig.20per)[j]
    
    # Perform Spearman correlation test
    tmp <- cor.test(
      mRNA.gene.count.trig.20per[[trig]] %>% as.numeric(),
      mRNA.gene.count.trig.20per[[mRNA]],
      method = "spearman",
      exact = FALSE,
      na.action = na.omit
    )
    
    # Create a data frame for the current correlation test
    tmp.df <- data.frame(
      trig = trig,
      mRNAs = mRNA,
      estimate = as.numeric(tmp$estimate),
      pvalue = tmp$p.value,
      stringsAsFactors = FALSE
    )
    
    # Append the results to pcor.test.res
    pcor.test.res.mRNA.ClinicalProfile <- bind_rows(pcor.test.res.mRNA.ClinicalProfile, tmp.df)
  }
}


pcor.test.res.mRNA.ClinicalProfile %>% arrange(desc(abs(estimate)))
pcor.test.res.mRNA.ClinicalProfile %>% arrange(desc(abs(estimate))) %>% dplyr::filter(pvalue < 0.05 & abs(estimate) > 0.5) -> mRNA.trig.0.5

# Loop through each mRNA column
for (i in 1:ncol(mRNA_count.S)-1) {
  for(j in 1:ncol(metabolites.data)){
   mRNA_name <- colnames(mRNA_count.S)[i]
   metabolite_name <- colnames(metabolites.data)[j] 

  paired_data <- data.frame(mRNA = mRNA_count.S[, i],
                          metabolite = plasma_meta.S[, j])  # Include other relevant variables
  
  # Remove missing values
  paired_data <- na.omit(paired_data)
  
 
  model <- glm(mRNA ~ metabolite, data = paired_data)
  
  # Store the model in the list
  regression_models[[mRNA_name]] <- model
  }
}


