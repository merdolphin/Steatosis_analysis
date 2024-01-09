setwd("/home/lina/scratch/18.NASH/Rscripts")


ClinicalData <- read_xlsx("/home/lina/scratch/18.NASH/1.sample/HumanPatients_230216.xlsx",
                          sheet = "Sheet1")

colnames(ClinicalData) -> tmp1
write.csv(tmp1,"colnames.1.csv")
colnames.1 <- read.table("colnames.3.csv", sep=",")
colnames(ClinicalData) <- colnames.1

ClinicalData %>% dplyr::select(ID, Trig, LDL,HDL, TotalChol, `NonHDL-Chol`, ALT, AST,GGT, TotalBili, Ferritin,
                               `HOMA2%B`,`HOMA2%S`,HOMA2IR) %>%
  mutate_if(is.character, as.numeric)   %>% 
  mutate(ID=as.character(ID))  -> C.trig

norm_mRNA_count <- read.delim("../4.NAFLD_Transcriptomics/NAFLD_Transcriptomics/CleanData/gene_level.txt",
                              header=T, as.is=T)

norm_mRNA_count[,1:4] -> mRNA.gene.dict
norm_mRNA_count[,-c(1:4)] -> mRNA.gene.count
names(mRNA.gene.count)  <- gsub(x=names(mRNA.gene.count),pattern="S", replacement = "")

C.trig %>% na.omit() %>% str()

mRNA.gene.count %>% t() %>% as.data.frame() %>% rownames_to_column(var="ID") %>% left_join(C.trig) %>% 
  dplyr::select(-ID)  -> mRNA.gene.count.trig

mRNA.gene.count.trig 
0.2*nrow(mRNA.gene.count.trig)


mRNA.gene.count.trig[,colSums(mRNA.gene.count.trig[,1:60362]==0) < 18] -> mRNA.gene.count.trig.20per 

mRNA.gene.count.trig.20per %>% dim()


# Initialize an empty data frame to store results
pcor.test.res.mRNA <- data.frame(mRNAs = character(),
                            trig = character(),
                            estimate = numeric(),
                            pvalue = numeric(),
                            stringsAsFactors = FALSE)

C.trig %>% dim()
mRNA.gene.count.trig.20per
colnames(mRNA.gene.count.trig.20per)[18160]
# Loop through pairs of trig and mRNA
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
    pcor.test.res.mRNA <- bind_rows(pcor.test.res.mRNA, tmp.df)
  }
}


pcor.test.res.mRNA %>% arrange(desc(abs(estimate)))
pcor.test.res.mRNA %>% arrange(desc(abs(estimate))) %>% dplyr::filter(pvalue < 0.05 & abs(estimate) > 0.5) -> mRNA.trig.0.5

pcor.test.res.mRNA %>% arrange(desc(abs(estimate))) %>% dplyr::filter(pvalue < 0.05 & abs(estimate) > 0.5) %>%
  dplyr::select(mRNAs)  -> mRNAtrig.8

mRNAtrig.8
pcor.test.res.mRNA %>% filter(mRNAs %in% mRNAtrig.8$mRNAs) -> mRNA.trig.0.5
# Convert the modified string into a vector of individual values
selected_row_ids <- unlist(strsplit(gsub("V", "", mRNA.trig.0.5$mRNAs), ", "))
selected_row_ids

result.mRNA <- mRNA.gene.dict[selected_row_ids, ]
result.mRNA
result.mRNA %>% dplyr::select(HGNC.symbol) -> result.mRNA.1
cbind(mRNA.trig.0.5, result.mRNA.1) %>%
  filter(str_detect(HGNC.symbol, "\\S")) -> mRNA.nonHDL.Chol

mRNA.nonHDL.Chol %>% dplyr::select(HGNC.symbol) %>% distinct()



mRNA.nonHDL.Chol %>% filter(pvalue > 0.05)

mRNA.nonHDL.Chol %>%
  mutate(trig = factor(trig, levels=c("Trig", "LDL", "HDL", "TotalChol", "NonHDL-Chol", 
                                      "ALT", "AST", "GGT", "TotalBili", "Ferritin",
                                      "HOMA2%B", "HOMA2%S", "HOMA2IR"))) %>%
  mutate(shape=ifelse(pvalue<0.05,1,NA)) %>%
  ggplot(aes(x = HGNC.symbol, y = trig, fill = estimate)) +
  geom_tile() +
  geom_point(aes(size=shape),color="white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value = "gray") +
  clean_background +
  labs(
    title = "Heatmap of Gene-Triglyceride Relationships",
    x = "Genes",
    y = "Triglyceride",
    fill = "Estimate"
  )+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate gene labels by 45 degrees
  ) +
  scale_y_discrete(labels = c("Trig" = "Triglycerides", 
                              "LDL" = "LDL",
                              "HDL" = "HDL",
                              "TotalChol" = "Total Chol.",
                              "NonHDL-Chol" = "Non-HDL Chol.",
                              "ALT" = "ALT",
                              "AST" = "AST",
                              "GGT" = "GGT",
                              "TotalBili" = "TBIL",
                              "Ferritin" = "Ferritin",
                              "HOMA2%B" = "HOMA2%B",
                              "HOMA2%S" = "HOMA2%S",
                              "HOMA2IR" = "HOMA2-IR")) 


liver.meta.count <- read.delim("../3.metabolimic/FromChoi/Normalized_data/liver_data_290523.txt",
                               header=T, as.is=T)
liver.meta.count[,1:30] -> liver.meta.dict
liver.meta.count[,-c(1:30)] -> liver.meta.count
names(liver.meta.count) <- gsub(x=names(liver.meta.count), pattern="X", replacement = "")

liver.meta.count %>% t() %>% as.data.frame() %>% rownames_to_column(var="ID") %>% 
  left_join(C.trig, by=c("ID" = "ID")) %>%
  dplyr::select(-ID) -> liver.meta.count.trig
C.trig

liver.meta.count.trig %>% dim()
# Initialize an empty data frame to store results
pcor.test.res <- data.frame(mRNAs = character(),
                            trig = character(),
                            estimate = numeric(),
                            pvalue = numeric(),
                            stringsAsFactors = FALSE)

# Loop through pairs of trig and mRNA
for (i in 1:1352) {
  for (j in 1353:1365) {
    trig <- colnames(liver.meta.count.trig)[j]
    mRNA <- colnames(liver.meta.count.trig)[i]
    
    # Perform Spearman correlation test
    tmp <- cor.test(
      liver.meta.count.trig[[trig]] %>% as.numeric(),
      liver.meta.count.trig[[mRNA]],
      method = "spearman",
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
    pcor.test.res <- bind_rows(pcor.test.res, tmp.df)
  }
}

pcor.test.res %>% arrange(desc(abs(estimate))) %>% dplyr::filter(pvalue < 0.05 & abs(estimate) > 0.5) 
pcor.test.res %>% arrange(desc(abs(estimate))) %>% dplyr::filter(pvalue < 0.05 & abs(estimate) > 0.5) -> livertrig.0.5

pcor.test.res %>% arrange(desc(abs(estimate))) %>% dplyr::filter(pvalue < 0.05 & abs(estimate) > 0.5) %>%
  dplyr::select(mRNAs)  -> livertrig.4

pcor.test.res %>% filter(mRNAs %in% livertrig.0.5$mRNAs) -> liver.trig.0.5

# Convert the modified string into a vector of individual values
selected_row_ids <- unlist(strsplit(gsub("V", "", liver.trig.0.5$mRNAs), ", "))


result <- liver.meta.dict[selected_row_ids, ]
result
result %>% dplyr::select(Compound)
result %>% dplyr::select(abbrev..LMSD.) -> result.LMSD
result %>% dplyr::select(accession..HMDB.)


result %>% dplyr::select(abbrev..LMSD., InChIKey, accession..HMDB.) %>% filter(str_detect(abbrev..LMSD., "\\S")) %>% 
  separate_rows(InChIKey, sep = ",") %>% distinct() %>% as.data.frame()


cbind(liver.trig.0.5, result.LMSD) %>%
  filter(str_detect(abbrev..LMSD., "\\S")) -> liver.nonHDL.Chol

liver.nonHDL.Chol  %>% mutate(cat = str_extract(abbrev..LMSD., "^[^ ]+")) %>% arrange(estimate) %>%
  mutate(trig = factor(trig, levels = c("Trig", "LDL", "HDL", "TotalChol", "NonHDL-Chol", 
                                        "ALT", "AST", "GGT", "TotalBili", "Ferritin",
                                        "HOMA2%B", "HOMA2%S", "HOMA2IR"))) %>%
  mutate(shape = ifelse(pvalue < 0.05, 1, NA)) -> liver.data



liver.data %>% dplyr::select(abbrev..LMSD.) %>% distinct() 
plasma.data %>% dplyr::select(abbrev..LMSD.) %>% distinct() %>% dim()

p1 <- liver.data %>% mutate(cat = factor(cat, 
                                          levels=c( "Cer", "PC", "PS",   "PE" , "SM"  ))) %>%
  ggplot(aes(x = abbrev..LMSD., y = trig, fill = estimate)) +
  geom_tile() +
  geom_point(aes(size = shape), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value = "gray") +
  #theme_void() +
  clean_background +
  labs(
    title = "Heatmap of Gene-Triglyceride Relationships",
    x = "Metabolite Categories",
    y = "Triglyceride",
    fill = "Estimate"
  ) +
  scale_y_discrete(labels = c("Trig" = "Triglycerides", 
                              "LDL" = "LDL",
                              "HDL" = "HDL",
                              "TotalChol" = "Total Chol.",
                              "NonHDL-Chol" = "Non-HDL Chol.",
                              "ALT" = "ALT",
                              "AST" = "AST",
                              "GGT" = "GGT",
                              "TotalBili" = "TBIL",
                              "Ferritin" = "Ferritin",
                              "HOMA2%B" = "HOMA2%B",
                              "HOMA2%S" = "HOMA2%S",
                              "HOMA2IR" = "HOMA2-IR")) +
  theme( axis.text.x = element_text(angle = 90, hjust = 1,size = 12, face = "bold"),
         axis.line = element_line("gray"),
         panel.spacing.x = unit(-0.2, "lines"),
         strip.text.x = element_text(size = 10, face = "bold"),
         panel.border = element_rect(colour = "gray", fill = NA, size = 3)
  )

pdf("liver_clinical.pdf",width=11)

p1+ ggforce::facet_row(vars(cat), scales='free', space='free')

dev.off()


plasma.meta.count <- read.delim("../3.metabolimic/FromChoi/Normalized_data/plasma_data_290523.txt",
                                header=T, as.is=T)
plasma.meta.count[,1:30] -> plasma.meta.dict
plasma.meta.count[,-c(1:30)] -> plasma.meta.count
names(plasma.meta.count) <- gsub(x=names(plasma.meta.count), pattern="X", replacement = "")

plasma.meta.count %>% t() %>% as.data.frame() %>% rownames_to_column(var="ID") %>% left_join(C.trig) %>% 
  dplyr::select(-ID) -> plasma.meta.count.trig



# Initialize an empty data frame to store results
pcor.test.res <- data.frame(trig = character(),
                            mRNAs = character(),
                            estimate = numeric(),
                            pvalue = numeric(),
                            stringsAsFactors = FALSE)

# Loop through pairs of trig and mRNA
for (i in 1:1087) {
  for (j in 1088:1100) {
    trig <- colnames(plasma.meta.count.trig)[j]
    mRNA <- colnames(plasma.meta.count.trig)[i]
    
    # Perform Spearman correlation test
    tmp <- cor.test(
      plasma.meta.count.trig[[trig]] %>% as.numeric(),
      plasma.meta.count.trig[[mRNA]],
      method = "spearman",
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
    pcor.test.res <- bind_rows(pcor.test.res, tmp.df)
  }
}


pcor.test.res %>% arrange(desc(abs(estimate))) %>% dplyr::filter(pvalue < 0.05 & abs(estimate) > 0.5)
pcor.test.res %>% arrange(desc(abs(estimate))) %>% dplyr::filter(pvalue < 0.05 & abs(estimate) > 0.5) -> plasma.trig.0.5


pcor.test.res %>% arrange(desc(abs(estimate))) %>% dplyr::filter(pvalue < 0.05 & abs(estimate) > 0.5) %>%
  dplyr::select(mRNAs)  -> plasma.trig.0.5

pcor.test.res %>% filter(mRNAs %in% plasma.trig.0.5$mRNAs) -> plasma.trig.0.5

plasma.trig.0.5 %>% dplyr::select(trig) %>% distinct()

# Convert the modified string into a vector of individual values
selected_row_ids <- unlist(strsplit(gsub("V", "", plasma.trig.0.5$mRNAs), ", "))

result.plasma <- liver.meta.dict[selected_row_ids, ]
result.plasma  %>% str()
result.plasma %>% dplyr::select(Compound)
result.plasma %>% dplyr::select(abbrev..LMSD.) -> plasma.LMSD

result.plasma %>% dplyr::select(abbrev..LMSD., formula)

cbind(plasma.trig.0.5, plasma.LMSD) %>% filter(str_detect(abbrev..LMSD., "\\S")) -> plasma.noHDL.Chol

result.plasma %>% dplyr::select(abbrev..LMSD., InChIKey) %>% filter(str_detect(abbrev..LMSD., "\\S")) %>% 
  separate_rows(InChIKey, sep = ",") %>% distinct() -> result.plasma.1


result.plasma %>% dplyr::select(abbrev..LMSD., InChIKey, accession..HMDB., abbrev_chains..LMSD.) %>% filter(str_detect(abbrev..LMSD., "\\S")) %>% 
  separate_rows(InChIKey, sep = ",") %>% filter(grepl("^PC", abbrev..LMSD.)) %>% distinct() %>% as.data.frame()


result.plasma.1 %>% filter(grepl("LPC", abbrev..LMSD.)) %>% distinct() %>% as.data.frame() 

result.plasma.1 %>% filter(grepl("^PC", abbrev..LMSD.)) %>% distinct() %>% as.data.frame() %>% dplyr::select(abbrev..LMSD.) %>% distinct()

result.plasma %>% dplyr::select(Batch, abbrev..LMSD.,formula,abbrev_chains..LMSD.,Compound) %>% filter(grepl("LPC", abbrev..LMSD.))

plasma.noHDL.Chol %>% dplyr::select(estimate, trig, abbrev..LMSD.) %>% pivot_wider(names_from = trig, values_from = estimate) %>% 
  filter(grepl("^PC", abbrev..LMSD.)) %>% dplyr::select(abbrev..LMSD.)


result.plasma %>%  filter(grepl("LPC 18:0", abbrev..LMSD.)) %>% 


spearman_bond <- read_xlsx("Analysis_results/LPC/LPC_smiles.xlsx",sheet="Sheet6")
spearman_bond %>% str()
spearman_bond[[i]]
bond.cor.test.res <- data.frame(bondInform = character(),
                            clinicalParameters = character(),
                            estimate = numeric(),
                            pvalue = numeric(),
                            stringsAsFactors = FALSE)

for(i in 1:5){
  for(j in 7:19){
    cor.test(spearman_bond[[i]], 
             spearman_bond[[j]],
             method="spearman") -> tmp
    tmp.df <- data.frame(
      bondInform = colnames(spearman_bond)[i] ,
      clinicalParameters = colnames(spearman_bond)[j],
      estimate = as.numeric(tmp$estimate),
      pvalue = tmp$p.value,
      stringsAsFactors = FALSE
    )
    
    # Append the results to pcor.test.res
    bond.cor.test.res <- bind_rows(bond.cor.test.res, tmp.df)
  }
}


pdf("try.pdf")
bond.cor.test.res %>% arrange(desc(abs(estimate))) %>% 
  mutate(clinicalParameters = factor(clinicalParameters, levels = c("Trig", "LDL", "HDL", "TotalChol", "NonHDL-Chol", 
                                        "ALT", "AST", "GGT", "TotalBili", "Ferritin",
                                        "HOMA2%B", "HOMA2%S", "HOMA2IR"))) %>%
  mutate(shape1 = ifelse(pvalue < 0.05, 1, NA)) %>% 
  ggplot(aes(x=bondInform, y=clinicalParameters, fill=estimate, size=shape1)) +
  geom_tile() +  scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value = "gray", 
                                      name = "Correlation Estimate", breaks = seq(-1, 1, by = 0.25)) +
  geom_point() +
  labs(title = "Correlation Estimates for Bond Information",
       x = "Bond Information",
       y = "Clinical Parameters") +
  #theme(legend.position = "bottom") +
  clean_background

dev.off()

help("pivot_wider")

library(ggplot2)
library(ggtext)


smiles <- c(
  "CCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)O   LPC 15:0",
  "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)O   LPC 16:0",
  "CCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)O  LPC 17:0",
  "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)O LPC 18:0",
  "CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)O LPC 18:1"
)



# Create a data frame with metabolite categories
metabolite_categories <- data.frame(
  abbrev..LMSD. = c("Cer", "FAL", "FA", "LPC", "LPE", "LPG", "LPI", "LPS", "NAE", "PA", "PC", "PE", "PG", "PI"),
  category = c("Ceramides", "Fatty Acyls", "Fatty Acids", "Lysophosphatidylcholines", "Lysophosphatidylethanolamines", "Lysophosphatidylglycerols", 
               "Lysophosphatidylinositols", "Lysophosphatidylserines", "N-Acyl Ethanolamines", "Phosphatidic Acids", "Phosphatidylcholines", 
               "Phosphatidylethanolamines", "Phosphatidylglycerols", "Phosphatidylinositols")
)



plasma.noHDL.Chol %>% mutate(cat = str_extract(abbrev..LMSD., "^[^ ]+")) %>% arrange(estimate) %>%
  mutate(trig = factor(trig, levels = c("Trig", "LDL", "HDL", "TotalChol", "NonHDL-Chol", 
                                        "ALT", "AST", "GGT", "TotalBili", "Ferritin",
                                        "HOMA2%B", "HOMA2%S", "HOMA2IR"))) %>%
  mutate(shape = ifelse(pvalue < 0.05, 1, NA)) -> plasma.data
plasma.data %>% dplyr::select(trig) %>% distinct()

plasma.data$abbrev..LMSD. <- str_extract(plasma.data$abbrev..LMSD., "^[^,]+")
plasma.data %>% dplyr::select(abbrev..LMSD.) 
plasma.data %>% str()
library(ggforce)

p <- plasma.data %>% mutate(cat = factor(cat, 
                                         levels=c( "LPE" ,"PC", "FAL" , "LPC","LPS", "LPG" ,  "FA" ,
                                                   "Cer","NAE" ,"PI", "LPI", "PA" , "PE", "PG"))) %>%
  ggplot(aes(x = abbrev..LMSD., y = trig, fill = estimate)) +
  geom_tile() +
  geom_point(aes(size = shape), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value = "gray") +
  theme_void() +
  #clean_background +
  labs(
    title = "Heatmap of Gene-Triglyceride Relationships",
    x = "Metabolite Categories",
    y = "Triglyceride",
    fill = "Estimate"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1) # Rotate category labels by 45 degrees
  ) +
  scale_y_discrete(labels = c("Trig" = "Triglycerides", 
                              "LDL" = "LDL",
                              "HDL" = "HDL",
                              "TotalChol" = "Total Chol.",
                              "NonHDL-Chol" = "Non-HDL Chol.",
                              "ALT" = "ALT",
                              "AST" = "AST",
                              "GGT" = "GGT",
                              "TotalBili" = "TBIL",
                              "Ferritin" = "Ferritin",
                              "HOMA2%B" = "HOMA2%B",
                              "HOMA2%S" = "HOMA2%S",
                              "HOMA2IR" = "HOMA2-IR")) +
 theme( axis.text.x = element_text(angle = 90, hjust = 1,size = 12, face = "bold"),
        axis.line = element_line("gray"),
        panel.spacing.x = unit(-0.2, "lines"),
        strip.text.x = element_text(size = 10, face = "bold"),
        panel.border = element_rect(colour = "gray", fill = NA, size = 3)
        )


pdf("plasma_clinical.pdf",width=12)
plasma.data$cat %>% unique()
p+ ggforce::facet_row(vars(cat), scales='free', space='free')

dev.off()

liver.nonHDL.Chol %>% mutate(source="liver") -> liver.nonHDL.Chol
plasma.noHDL.Chol %>% mutate(source="plasma") -> plasma.noHDL.Chol
plasma.LMSD %>% dim()
liver.nonHDL.Chol
plasma.noHDL.Chol
rbind(liver.nonHDL.Chol, plasma.noHDL.Chol) %>% dplyr::select(1:4,abbrev..LMSD., source) -> Liver.Plasma.Clinical.association.0.5

rbind
mRNA.nonHDL.Chol %>% setNames(c( "mRNAs", "trig","estimate","pvalue","abbrev..LMSD.")) %>% 
  mutate(source="mRNA") %>%
  rbind(Liver.Plasma.Clinical.association.0.5) -> mat.cor.0.5

mat.cor.0.5 %>% dplyr::select(trig) 


# Convert abbrev..LMSD. to a factor for proper ordering in the plot
data$abbrev..LMSD. <- factor(data$abbrev..LMSD., levels = unique(data$abbrev..LMSD.))

# Create a scatter plot
ggplot(data, aes(y = abbrev..LMSD., x = 0, fill = estimate, label=abbrev..LMSD.)) +
  geom_tile() +
  xlim(-1,1) +
  labs(title = "Correlation between variables",
       y = "Abbreviation (LMSD)",
       x = "Estimate") + clean_background +
  scale_fill_gradient2(low="red",high = "blue",mid="white")


geom_text_repel(
  data = data %>% filter(estimate > 0 ),
  aes(
    segment.square  = FALSE,
    segment.inflect = TRUE
  ),
  box.padding = 0.5,
  max.overlaps = Inf,
  nudge_x = -1,
  direction = "y",
  vjust = "right",
  segment.size      = 0.5,
  segment.curvature = -0.2,
  size=5
  # color="red"
) +
  geom_text_repel(
    data = data %>% filter(estimate < 0 ),
    aes(
      segment.square  = FALSE,
      segment.inflect = TRUE
    ),
    box.padding = 0.5,
    max.overlaps = Inf,
    nudge_x = 1,
    direction = "y",
    vjust = "right",
    segment.size      = 0.5,
    segment.curvature = -0,
    size=5
    # color="red"
  ) +



library(ChemmineR)
library(PubChemR)
library(webchem)
# Function to get SMILES from InChIKey
getSmilesFromInChIKey <- function(inchikey) {
  compound <- pc_info(inchikey, type = "inchikey", quiet = FALSE)
  smiles <- compound$info$canonical_smiles
  return(smiles)
}

# Function to calculate molecular properties
calculateMolecularProperties <- function(smiles) {
  mol <- parse.smiles(smiles)
  numDoubleBonds <- num.doublebonds(mol)
  numRotatableBonds <- num.rotatablebonds(mol)
  maxSingleBondLength <- max.singlebondlength(mol)
  return(c(NumDoubleBonds = numDoubleBonds, NumRotatableBonds = numRotatableBonds, MaxSingleBondLength = maxSingleBondLength))
}


chemspider_api_key <- Sys.getenv("CHEMSPIDER_KEY")
chemspider_api_key
mat.nonHDL.chol.cor %>% separate_rows(InChIKey, sep=",") -> mat.nonHDL.chol.cor.1

mat.nonHDL.chol.cor.1$InChIKey 

# Assuming your data frame is named 'df'
library(ChemmineR)
library("ChemmineOB")

# Select the row of interest, e.g., the first row
row_index <- 1
inchi_key <- mat.nonHDL.chol.cor.1$InChIKey[row_index]
inchi_key
# Use ChemmineR to convert InChIKey to SMILES
smiles <- inchi2smiles(inchi_key)

ChemmineR::pubchemInchikey2sdf(inchikeys = inchi_key) 

# Replace 'your_inchi_key' with the actual InChIKey you have
inchi_key <- "HXFPPRPLRSPNIB-VARSQMIESA-N"

# Convert InChIKey to SMILES
 ChemmineR::pubchemInchikey2sdf(inchikeys = inchi_key)$sdf_set -> tmp

ChemmineR::sdf2smiles(tmp)%>% str()


# Print the result


?cs_check_key()

# Apply the functions to the data frame
mat.nonHDL.chol.cor <- mat.nonHDL.chol.cor %>%
  rowwise() %>%
  mutate(SMILES = getSmilesFromInChIKey(InChIKey),
         Properties = list(calculateMolecularProperties(SMILES))) %>%
  separate(Properties, into = c("NumDoubleBonds", "NumRotatableBonds", "MaxSingleBondLength"), sep = ",")

# Display the result
print(mat.nonHDL.chol.cor)

