library(tidyverse)
library(ggplot2)
library(ggrepel)



setwd("/home/lina/scratch/18.NASH/")
norm_mRNA_count <- read.delim("4.NAFLD_Transcriptomics/NAFLD_Transcriptomics/CleanData/gene_level.txt",
                              header=T, as.is=T)

isomer_mRNA_count <- read.delim("4.NAFLD_Transcriptomics/NAFLD_Transcriptomics/CleanData/isoform_level.txt",
                                header=T, as.is=T)

norm_mRNA_count[,1:4] -> mRNA.gene.dict
norm_mRNA_count[,-c(1:4)] -> mRNA.gene.count
names(mRNA.gene.count)  <- gsub(x=names(mRNA.gene.count),pattern="S", replacement = "")


isomer_mRNA_count[,1:5] -> isomer.gene.dict
isomer_mRNA_count[,-c(1:5)] -> isomer.gene.count
names(isomer.gene.count)  <- gsub(x=names(isomer.gene.count),pattern="S", replacement = "")

mRNA.gene.count %>% dim()
isomer_mRNA_count %>% dim()


mRNA.gene.count[rowSums(mRNA.gene.count == 0) < 98,] %>% dim()
isomer.gene.count[rowSums(isomer.gene.count == 0) < 98, ] %>% dim()


mRNA.gene.count.filtered <- norm_mRNA_count[rowSums(norm_mRNA_count[,-c(1:4)]) != 0, ] 
mRNA.gene.count.filtered %>% dplyr::select(Gene.stable.ID) %>% dim()

isomer_mRNA_count.filter <- isomer_mRNA_count[rowSums(isomer_mRNA_count[,-c(1:5)]) !=0, ]

isomer_mRNA_count.filter %>% dplyr::select(Gene.stable.ID) %>% distinct() %>% dim()
isomer_mRNA_count.filter %>% dplyr::select(Gene.stable.ID) %>% count(Gene.stable.ID) %>% 
  filter(n>1) %>% dim()


total <- 35927
slice <- 19716

(total - slice)/total
# Create a vector with the values
data <- c(slice, total - slice)

# Labels for the pie chart
labels <- c()

# Colors for each section
colors <- c("yellowgreen", "gray")

pdf("try.pdf")
# Create the pie chart
pie(data, label=NA,  col = colors, main = "Pie Chart")
dev.off()


load(file="S2vsS1vsS0FDR0.05.rda")

S2vsS1vsS0FDR0.05 

rbind(S2vsS1res0.05, S1vsS0res0.05) -> S2vsS1vsS0FDR0.05
S2vsS1vsS0FDR0.05 %>% rownames() -> DEGs.91

mRNA.gene.dict %>% filter(Gene.stable.ID  %in% DEGs.91) %>% dplyr::select(HGNC.symbol)

isomer_mRNA_count.filter %>% filter(Gene.stable.ID %in% DEGs.91) %>% dplyr::select(Gene.stable.ID) %>% 
  table() %>% enframe() %>% filter(value !=1) %>% dplyr::select(name) -> ENSG.multi.isomers.70

ENSG.multi.isomers.70

S2vsS1vsS0FDR0.05 %>% 
  as.data.frame() %>% dplyr::select(log2FoldChange)

result <- S2vsS1vsS0FDR0.05 %>% as.data.frame()

result

# Provided gene lists
metabolite_related_genes <- c(
  "ASAH1", "CYB5R3", "HAAO", "SLC27A4", "ELOVL1", 
  "CYP2C19", "OAT", "ACOX2", "AKR1B10", "FADS1", 
  "GLYCTK", "M6PR", "PLIN3", "CYB5R3", "HAAO", "SAT2", "SMOX"
)

inflammation_related_genes <- c(
  "ASAH1", "CYB5R3", "HAAO", "C2CD4A", "CAP1", 
  "CD68", "CHI3L1", "CYB561A3", "DYNLT1", "GRN", 
  "IL10RB", "LYZ", "MCUB", "MVP", "PLA2G7", "RTN4", 
  "S100A9", "SPI1", "TCIRG1"
)

# Integrated list
all_related_genes <- unique(c(metabolite_related_genes, inflammation_related_genes))

result %>%
  arrange(log2FoldChange) %>%
  mutate(log10FC = log2FoldChange / log10(2)) %>%
  mutate(rownames_factor = factor(rownames(.), levels = rownames(.))) %>%
  mutate(col = case_when(
    rownames(.) %in% all_related_genes ~ "purple",  # Change the color as needed
    TRUE ~ "gray"
  )) %>%
  ggplot(aes(x = rownames_factor, y = log2FoldChange, color = col)) +
  geom_point() +
  geom_text(aes(label = ifelse(col != "gray", as.character(rownames_factor), "")),
            vjust = -0.5, hjust = 1, size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Gene", y = "log2FoldChange") +
  clean_background


Inflammation_related_genes <- c("MVP", "TYMP", "GRN", "CYB5R3", "CD68", "CHI3L1", "S100A9", "ANXA2", "TGM2")
metabolite_related_genes <- c("ADH1C", "GNMT", "MT-ATP6", "MT-CYB", "GLYCTK", "HAO2", "OAT", "OAT", "GNMT", "RPL22L")

result %>%
  mutate(rownames_factor = factor(rownames(.), levels = rownames(.))) %>%
  mutate(col1 = case_when(
    rownames(.) %in% ENSG.multi.isomers.70$name ~ "1",
    TRUE ~ "Other"
  )) %>%
  mutate(col2 = case_when(
    mm_symbol %in% Inflammation_related_genes ~ "Inflammation",
    mm_symbol %in% metabolite_related_genes ~ "Metabolite",
    TRUE ~ "other"
  ))
  ggplot(aes(x = rownames_factor, y = log2FoldChange, color = col, label = mm_symbol)) +
  geom_point() +
  scale_color_manual(values = c("NA" = "gray", 
                                "1" = "darkgreen", 
                                "Inflammation" = "red", 
                                "Metabolite" = "blue", 
                                "Other" = "black")) +
  geom_text_repel(box.padding = 0.3, max.overlaps = Inf) +
  #  geom_text_repel(
   #   nudge_x = -2, direction = "x", vjust = "right"
  #  ) +
 #   geom_text_repel(
 #   nudge_x = 2, direction = "x", vjust = "left"
#    )+
 # coord_flip() +
  labs(x = "Gene", y = "log2FoldChange") +
  ylim(-3, 3) +
  geom_hline(yintercept = 0, color = "gray") +
  clean_background

  
