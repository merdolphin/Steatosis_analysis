  
library(ggplot2)
library(ggpubr)

Chronos.Score <- read.csv("/home/lina/scratch/18.NASH/Analysis_results/metaboanalysis_Deseq2FDR0.05_association_metabolite/total309metabolites/CRISPR_(DepMap_Public_23Q2+Score,_Chronos)_subsetted.csv",
                          header = T)


Chronos.Score[,1:8] %>% str()
Chronos.Score %>% dim()

selected_columns <- Chronos.Score[, 9:79]

# Use logical indexing to get values less than -0.5
selected_columns_less_than_minus_0.5 <- selected_columns < -0.5

# Use colSums to get the sum for each column
sums_less_than_minus_0.5 <- colSums(selected_columns_less_than_minus_0.5, na.rm = TRUE)

# Display or use the result as needed
print(sums_less_than_minus_0.5) %>% enframe() %>% arrange(desc(value)) %>% head(n=10) 
# Check if the values in Chronos.Score are numeric and less than -0.5

values_less_than_minus_0.5 <- Chronos.Score[sapply(Chronos.Score[, 9:79], is.numeric) & Chronos.Score[, 9:79] < -0.5,]



genes <- c("PPIL1", "CFL1",  "CAP1",  "TKT" ,  "RSRC2", "B3GAT3", "RRAD", "ELOVL1", "CD151","SPI1")

# Create an empty data frame to store results
noLiverDependentLines <- data.frame(Gene = character(), Liver.number = numeric(), Liver.number = numeric(), stringsAsFactors = FALSE)

Chronos.Score %>%
  dplyr::select(c(1:8, c(PPIL1, CFL1))) %>%
  filter(lineage_1 == "Liver", PPIL1 < -0.5)

for (gene in genes) {
  # Filter data for liver and genes with scores less than -0.5
  liverDependentLines <- Chronos.Score %>%
    dplyr::select(c(1:8, all_of(gene))) %>%
    filter(lineage_1 == "Liver", .data[[gene]] < -0.5)
  
  # Count the number of lines
  Liver.number <- nrow(liverDependentLines)
  Total.number <- liverDependentLines <- Chronos.Score %>%
    dplyr::select(c(1:8, all_of(gene))) %>%
    filter(.data[[gene]] < -0.5) %>% nrow()
  
  # Create a temporary data frame
  tmp <- data.frame(Gene = gene, Liver.number =  Liver.number, TotalNumber = Total.number, stringsAsFactors = FALSE)
  
  # Append the temporary data frame to the main results
  noLiverDependentLines <- rbind(noLiverDependentLines, tmp)
}

# Print or use the results as needed
noLiverDependentLines %>% arrange(desc(Liver.number))

load(file= "mRNA_kruskal_005.rda")


gsub("V","",mRNA_kruskal_005) -> mRNA_kruskal_005.key 
mRNA_kruskal_005.key
mRNA.gene.dict[mRNA_kruskal_005.key,] %>% dplyr::select(HGNC.symbol) -> AD.gene.list.5485


PPIL1interactionGenes <- c("AQR","BUD31","CDC40","CDC5L","EFTUD2","PPIL1","PRPF19","PRPF8","RBM22","SNW1","SYF2")
CFL1interactionGenes <- c("ACTB","ACTG1","CFL1","CTTN","HCLS1","LIMK1","LIMK2","PFN1","SSH1","WASL","WDR1")

intersect(AD.gene.list.5485$HGNC.symbol, PPIL1interactionGenes) 
intersect(AD.gene.list.5485$HGNC.symbol, CFL1interactionGenes)

PPIL1interactionGenes <- c("AQR","BUD31","CDC40","CDC5L","EFTUD2","PPIL1","PRPF19","PRPF8","RBM22","SNW1","SYF2")
CFL1interactionGenes <- c("ACTB","ACTG1","CFL1","CTTN","HCLS1","LIMK1","LIMK2","PFN1","SSH1","WASL","WDR1")

intersect(AD.gene.list.5485$HGNC.symbol, PPIL1interactionGenes) 
intersect(AD.gene.list.5485$HGNC.symbol, CFL1interactionGenes)



essential_genes <- c("PPIL1","CFL1")

mRNA.gene.dict %>% rownames_to_column(var = "rowid") %>%
  filter(HGNC.symbol %in% essential_genes) %>% dplyr::select(rowid, HGNC.symbol) %>%
  mutate(rowid = paste0("V", rowid)) -> mRNA.gene.dict.2

mRNA.gene.dict.2


mRNA_count.S %>% dplyr::select(any_of(mRNA.gene.dict.2$rowid),Histo_STEATOSIS) %>% 
  gather(key, value, -Histo_STEATOSIS) -> df2

df2

df2 %>% left_join(mRNA.gene.dict.2, by=c("key" = "rowid")) -> df2.1

subset(S2vsS1vsS0FDR0.05, mm_symbol %in% importantGenes) %>% as.data.frame() %>% arrange(log2FoldChange) %>% 
  dplyr::select(mm_symbol)



subset(S2vsS1vsS0FDR0.05, mm_symbol %in% essential_genes) %>% as.data.frame() %>% arrange(log2FoldChange) -> df2.2

my_comparisons <- list(c("2","1"),c("0","2"))


pdf("try.pdf", width=8, height=4)
ggboxplot(df2.1, x="Histo_STEATOSIS", y="value", fill="Histo_STEATOSIS", color="Histo_STEATOSIS",
          alpha=0.2, outlier.shape = NA, add = "jitter",
          add.params = list(alpha = 0.5),width=0.8) +
  stat_compare_means(label.y = (1.5*min(df2.1$value)))+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", 
                     method = "wilcox.test", paired = F, hide.ns = T, show.legend = F, size = 8,
                     vjust = 0.5, tip.length = 0.03,
                     bracket.size = 0.3)+
  scale_fill_npg() +
  scale_color_npg() +
  facet_wrap(~HGNC.symbol, scales="free_y", ncol=2)+
  clean_background +
  scale_x_discrete(labels=c("S0","S1","S2+S3")) +
  # Set y-axis limits # Replace lower_limit and upper_limit with your desired values
  #+ coord_cartesian(ylim = ylim_values) +
  scale_fill_discrete(name = "Steatosis", labels  = c("S0 (n = 37)", "S1 (n = 41)", "S2+S3: (n = 36)")) +
  labs(x = "Steatosis level",
       y = "Normalized level") + # guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 16)
  ) 

dev.off()

#################### Association with liver biomarkers ########
ClinicalData <- read_xlsx("/home/lina/scratch/18.NASH/1.sample/HumanPatients_230216.xlsx",
                          sheet = "Sheet1")

colnames.1 <- read.table("/home/lina/scratch/18.NASH/Rscripts/colnames.3.csv", sep=",")
colnames(ClinicalData) <- colnames.1
colnames(ClinicalData)
colSums(is.na(ClinicalData)) 


ClinicalData %>% dplyr::select(contains("Histo")) %>% mutate_if(is.numeric, as.factor) %>%
  summary()

ClinicalData$Diabetes[ClinicalData$Diabetes == "2"] <- "1"
ClinicalData$T2DM[ClinicalData$T2DM == "2"] <- "1"


ClinicalData[,colSums(is.na(ClinicalData))<25] -> CD.1


CD.1 %>% mutate(Gender = factor(Gender, levels=c(0,1), labels=c("Female","Male")), 
                Diabetes = as.factor(Diabetes),
                Previous.gestational.DM = as.factor(Previous.gestational.DM),
                Metformin_treatment = factor(Metformin_treatment, levels=c(0,1), labels=c("No","Yes")), 
                T2DM = as.factor(T2DM) 
                #Chol_Rx_name = as.factor(Chol_Rx_name)
) -> CD.2


CD.2 %>% dplyr::select(ID, Histo_STEATOSIS, ALT, AST, `ALT/AST`, GGT, ALP, TotalBili, Albumin, Ferritin,
                       Urea,SerumCreatine, eGFR ) -> CD.2


pdf("try.pdf", width=6, height=4)
CD.2 %>% dplyr::select(Histo_STEATOSIS, `HOMA2%B`, `HbA1c(%)`, `HOMA2%S`, HOMA2IR) %>% na.omit() %>% 
  gather(key = "Variable", value = Value, -Histo_STEATOSIS) %>%
  ggplot(aes(x = Histo_STEATOSIS, y = as.numeric(Value), color = Variable)) +
  #geom_point() +
  geom_smooth(method = "lm", se=FALSE, linetype="dashed", linewidth=2, span=1, level=0.7) +  # You can change the method if needed
  labs(title = "Smoothed Lines for Variables by Histo_STEATOSIS",
       x = "Histo_STEATOSIS",
       y = "Value") +
  scale_y_log10() +
  ylim(0,150) +
  clean_background
dev.off()

help("geom_smooth")


mRNA.gene.count[substr(mRNA.gene.dict.2$rowid, 2, nchar(mRNA.gene.dict.2)),] -> essentialGeneCounts.1
essentialGeneCounts.1 %>% t() -> essentialGeneCounts.2


CD.2 %>% filter(ID %in% rownames(essentialGeneCounts.2)) -> CD.2.2

common_ids <- intersect(CD.2$ID, rownames(essentialGeneCounts.2))
CD.2.2 <- CD.2[CD.2$ID %in% common_ids, ]
essentialGeneCounts.2[rownames(essentialGeneCounts.2) %in% common_ids,] -> essentialGeneCounts.3
essentialGeneCounts.3[match(rownames(essentialGeneCounts.3),CD.2.2$ID),] -> essentialGeneCounts.4
essentialGeneCounts.4 %>% as.data.frame() %>% rownames_to_column(var="ID") %>% mutate(ID=as.numeric(ID)) -> essentialGeneCounts.5

CD.2.2 %>% left_join(essentialGeneCounts.5) %>% mutate(across(c(3:13), as.numeric)) %>%
  na.omit()-> essential2genes.sam

essential2genes.sam %>% aggr(plot=F)

# Create an empty result data frame
pcor.test.res <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(pcor.test.res) <- c("gene", "biomarkers", "estimate", "pvalue")
colnames(essential2genes.sam)


for (i in 2:55) {
  biomarker <- colnames(essential2genes.sam)[i]
  tmp1 <- cor.test(essential2genes.sam$`13357`, essential2genes.sam[[biomarker]], method="spearman")
  tmp1.df <- data.frame(
     gene = "CFL1",
     biomarkers = biomarker,
    estimate = as.numeric(tmp1$estimate),
    pvalue = tmp1$p.value
  )
  tmp2 <- cor.test(essential2genes.sam$`7418`, essential2genes.sam[[biomarker]], method="spearman")
  tmp2.df <- data.frame(
    gene = "PPIL1",
    biomarkers = biomarker,
    estimate = as.numeric(tmp2$estimate),
    pvalue = tmp2$p.value
  )
  tmp.df <- rbind(tmp1.df, tmp2.df)
  pcor.test.res <- rbind(pcor.test.res, tmp.df) 
}


pcor.test.res %>% arrange(desc(abs(estimate)))


essential2genes.sam %>%
  dplyr::select(AST, `7418`, `13357`) %>%
  setNames(c("AST", "PPIL1", "CFL1")) -> essential2Genes.AST
essential2Genes.AST %>% str()
ggscatter(essential2Genes.AST, x="AST", y = "CFL1",  color = "AST", size = "AST",
            conf.int = TRUE) +
  stat_cor(method = "spearman") +
  scale_color_gradient2(low = "red", high = "blue", midpoint = 75) +
  labs(title = "Association of Genes with AST",
       x = "AST U/L",
       y = "CFL1 normalized level",
       color = "AST Levels") +
  clean_background



# Create plot for CFL1
plot_cfl1 <- ggplot(essential2Genes.AST, aes(x = AST, y = CFL1, color = AST, size = AST)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x = 80, label.y = 370, size=5) +
  scale_color_gradient2(low = "red", high = "blue", midpoint = 75) +
  labs(title = "Association of CFL1 with AST",
       x = "AST (U/L)",
       y = "CFL1 Normalized Level",
       color = "AST Levels") +
  clean_background +
  theme(axis.title.y.right = element_text(color = "blue"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, vjust = 2))

plot_cfl1
# Create plot for PPIL1
plot_ppil1 <- ggplot(essential2Genes.AST, aes(x = AST, y = PPIL1, color = AST, size = AST)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE,  color = "black") +
  stat_cor(method = "spearman", label.x = 80, label.y = 17, size=5) +
  scale_color_gradient2(low = "red", high = "blue", midpoint = 75) +
  labs(title = "Association of PPIL1 with AST",
       x = "AST (U/L)",
       y = "PPIL1 Normalized Level",
       color = "AST Levels") +
  clean_background +
  theme(axis.title.y.right = element_text(color = "blue"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, vjust = 2)) 

plot_ppil1

pdf("try.pdf", width=9, height=5.5)
ggarrange(plot_ppil1, plot_cfl1)
dev.off()




######### association_with_liver_metabolites ########
mRNA.gene.dict.2

mRNA_count.S %>% dplyr::select(any_of(mRNA.gene.dict.2$rowid),Histo_STEATOSIS) -> essentialGeneCounts

mRNA_count.S %>% rownames() == rownames(liver_meta.S)

essentialGeneCounts

# Create an empty result data frame
pcor.test.res <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(pcor.test.res) <- c("mRNA", "liver", "estimate", "pvalue")

# Perform parallelized loop
pcor.test.res

for (i in 1:ncol(essentialGeneCounts[,-1]))
{
  for (j in 1:ncol(liver_meta.S[,-1])) {
    tmp <-
      cor.test(essentialGeneCounts[, i],
               liver_meta.S[, j],
               method = "spearman")
    tmp.df <- data.frame(
      mRNA = colnames(essentialGeneCounts)[i],
      liver = j,
      estimate = as.numeric(tmp$estimate),
      pvalue = tmp$p.value
    )
    pcor.test.res <- rbind(pcor.test.res, tmp.df)
  }
}

# Print the result
print(pcor.test.res)
pcor.test.res %>% filter(pvalue < 0.05 & abs(estimate) > 0.3)
pcor.test.res %>% filter(pvalue < 0.05 & abs(estimate) > 0.3)  %>% dplyr::select(liver) -> rho0.3liverID
pcor.test.res %>% filter(pvalue < 0.05 & abs(estimate) > 0.4 ) %>% dplyr::select(liver, estimate, pvalue) -> liver.ID
liver.ID

liver.meta.dict[liver.ID$liver,] 

liver.meta.dict[rho0.3liverID$liver,] %>% dplyr::select(accession..HMDB.) %>% 
  mutate(accession..HMDB. = str_replace_all(accession..HMDB., "[\" ]", "")) %>%
  separate_rows(accession..HMDB., sep = ",") %>% na.omit %>% as.data.frame()

liver.meta.dict[liver.ID$liver,] %>% dplyr::select( abbrev..LMSD., sub.class..HMDB.) %>% cbind(liver.ID) %>%
  dplyr::select(-liver) %>% setNames(c("Name", "Class", "Rho", "pvalue"))
dplyr::select(accession..HMDB.) %>% as.vector()

# Assuming your data frame is named result_data

result_data <- liver.meta.dict[liver.ID$liver,] %>% 
  dplyr::select(abbrev..LMSD., sub.class..HMDB.) %>%
  cbind(liver.ID) %>%
  dplyr::select(-liver) %>%
  setNames(c("Name", "Class", "Rho", "pvalue"))

result_data
# Write the data frame to a text file with tab as the delimiter
write.table(result_data, file = "output_file.txt", sep = "\t", quote = FALSE, row.names = FALSE)


liver_meta.S %>% dplyr::select(V940, V952, Histo_STEATOSIS)  -> data


# Assuming data is your dataframe
my_comparisons <- list(c("2","1"),c("0","1"), c("0","2"))

data %>% 
  gather(key, value, -Histo_STEATOSIS) -> data.1

pdf("try.pdf", width = 8, height = 4)
ggviolin(data.1, x="Histo_STEATOSIS", y="value", fill="Histo_STEATOSIS", color="black",
                        alpha=0.2, outlier.shape = NA, add = "jitter",
                        add.params = list(alpha = 0.5),width=0.8) +
  stat_compare_means(label.y = (0.8*min(data.1$value)))+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", 
                     method = "wilcox.test", paired = F, hide.ns = T, show.legend = F, size = 8,
                     vjust = 0.5, tip.length = 0.03,
                     bracket.size = 0.3)+
  scale_fill_npg() +
 # scale_color_npg() +
  facet_wrap(~key, scales="free_y", ncol=2, labeller = labeller(key = c("V940" = "PE 34:2", "V952" = "PC 35:5"))) +
  clean_background +
  scale_x_discrete(labels=c("S0","S1","S2+S3")) +
  # Set y-axis limits # Replace lower_limit and upper_limit with your desired values
  #+ coord_cartesian(ylim = ylim_values) +
  scale_fill_discrete(name = "Steatosis", labels  = c("S0 (n = 37)", "S1 (n = 41)", "S2+S3: (n = 36)")) +
  labs(x = "Steatosis level",
       y = "Normalized level") + # guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 16)
  ) 
dev.off()
######### association_with_plasma_metabolites ########
##### no association was found

mRNA_count.S %>% rownames() == rownames(plasma_meta.S) 

essentialGeneCounts

# Create an empty result data frame
pcor.test.res <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(pcor.test.res) <- c("mRNA", "plasma", "estimate", "pvalue")

# Perform parallelized loop
pcor.test.res

for (i in 1:ncol(essentialGeneCounts[,-1]))
{
  for (j in 1:ncol(plasma_meta.S[,-1])) {
    tmp <-
      cor.test(essentialGeneCounts[, i],
               plasma_meta.S[, j],
               method = "spearman")
    tmp.df <- data.frame(
      mRNA = colnames(essentialGeneCounts)[i],
      plasma = j,
      estimate = as.numeric(tmp$estimate),
      pvalue = tmp$p.value
    )
    pcor.test.res <- rbind(pcor.test.res, tmp.df)
  }
}

# Print the result
print(pcor.test.res)

pcor.test.res %>% filter(pvalue < 0.05 & abs(estimate) > 0.4) 
p
