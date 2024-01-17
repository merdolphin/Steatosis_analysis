# Load necessary libraries
library(vegan)  # for the Mantel test
library(proxy)
library(dplyr)
library(tibble)
library(readxl)
library(tidyverse)
library(stringr)
library(mixOmics)
library(VIM)
setwd("/home/lina/scratch/18.NASH/")
## Data preparation and normalization




liver.meta.count.0 <- read.delim("3.metabolimic/FromChoi/Normalized_data/liver_data_290523.txt",
                               header=T, as.is=T)

plasma.meta.count.0 <- read.delim("3.metabolimic/FromChoi/Normalized_data/plasma_data_290523.txt",
                                header=T, as.is=T)

CM_metabolomics_data.0 <- read.delim("3.metabolimic/FromChoi/Normalized_data/cm_data_290523.txt",
                                   header=T, as.is = T)


norm_mRNA_count.0 <- read.delim("4.NAFLD_Transcriptomics/NAFLD_Transcriptomics/CleanData/gene_level.txt",
                              header=T, as.is=T)

norm_mRNA_count.0[,1:4] -> mRNA.gene.dict
norm_mRNA_count.0[,-c(1:4)] -> mRNA.gene.count
names(mRNA.gene.count)  <- gsub(x=names(mRNA.gene.count),pattern="S", replacement = "")
mRNA.gene.count

liver.meta.count.0[,1:30] -> liver.meta.dict
liver.meta.count.0[,-c(1:30)] -> liver.meta.count
names(liver.meta.count) <- gsub(x=names(liver.meta.count), pattern="X", replacement = "")

plasma.meta.count.0[,1:30] -> plasma.meta.dict
plasma.meta.count.0[,-c(1:30)] -> plasma.meta.count
names(plasma.meta.count) <- gsub(x=names(plasma.meta.count), pattern="X", replacement = "")

identical(liver.meta.dict, plasma.meta.dict)
liver.meta.dict

CM_metabolomics_data[,1:30] -> CM.meta.dict
CM_metabolomics_data[,-c(1:30)] -> CM.meta.count
names(CM.meta.count) <- gsub(x=names(CM.meta.count), pattern="X", replacement = "")
mRNA.gene.count %>% dim()
liver.meta.count[, colSums(is.na(liver.meta.count)) < 10 ] -> liver.meta.count
plasma.meta.count[, colSums(is.na(plasma.meta.count)) < 10] -> plasma.meta.count
CM.meta.count[, colSums(is.na(CM.meta.count)) < 10] -> CM.meta.count

intersect(colnames(mRNA.gene.count), colnames(plasma.meta.count)) -> RNA_plasma.80
intersect(colnames(mRNA.gene.count), colnames(liver.meta.count)) -> RNA_liver.78
intersect(intersect(colnames(mRNA.gene.count), colnames(liver.meta.count)), 
          intersect(colnames(CM.meta.count),  colnames(plasma.meta.count))) -> RNA_plasma_liver.90

liver.meta.count %>% dplyr::select(all_of(RNA_plasma_liver.90)) %>% scale() -> liver.meta.1
plasma.meta.count %>% dplyr::select(all_of(RNA_plasma_liver.90)) %>% scale() -> plasma.meta.1
CM.meta.count %>% dplyr::select(all_of(RNA_plasma_liver.90)) %>% scale() -> CM.meta.1
mRNA.gene.count %>% as.data.frame() %>% dplyr::select(all_of(RNA_plasma_liver.90)) %>% scale()-> mRNA.gene.count.1

colnames(mRNA.gene.count.1) == colnames(liver.meta.1)
colnames(mRNA.gene.count.1) == colnames(plasma.meta.1)

samples_information <- read_xlsx("1.sample/HumanPatients_230216.xlsx")
samples_information %>% dplyr::select(ID, `Histo_STEATOSIS (0-3)`,  `Trig (<1.5)`, `LDL (<3.5 mmol/L)`,`HDL (>1.0 mmol/L = normal)`, 
                                      `TotalChol (3.5-5.5 mM)`, `NonHDL-Chol (mmol/L)`, `ALT(5-30`, `AST (10-35)`,
                                      `GGT (5-35 U/L)`, `Total Bili (3-15)`, `Ferritin (20-500 ng/mL)`,
                                      `HOMA2%B`,`HOMA2%S`,HOMA2IR) %>% 
  setNames(c("ID", "Histo_STEATOSIS","Trig", "LDL", "HDL", "TotalChol", "NonHDL-Chol", "ALT", "AST", "GGT", "TotalBili", "Ferritin", "HOMA2B", "HOMA2S", "HOMA2IR")) %>%
  mutate(Histo_STEATOSIS = str_replace_all(Histo_STEATOSIS, "3","2")) %>% 
  mutate(ID=as.character(ID)) %>%
  mutate(Histo_STEATOSIS = as.factor(Histo_STEATOSIS)) -> sampled.ID.Steatosis

samples_information %>% dplyr::select(ID, `Histo_INFLAM(0-3)`) %>% 
  setNames(c("ID", "Histo_INFLAM")) %>%
  mutate(Histo_INFLAM = str_replace_all(Histo_INFLAM, "3","2")) %>% 
  mutate(ID=as.character(ID)) %>%
  mutate(Histo_INFLAM = as.factor(Histo_INFLAM)) -> sampled.ID.INFLAM

liver.meta.1 %>% t() %>% as.data.frame() %>% rownames_to_column(var="ID.1") %>% 
  left_join(sampled.ID.Steatosis, by=c("ID.1" = "ID")) %>% dplyr::select(-ID.1) -> liver_meta.S

liver.meta.1 %>% t() %>% as.data.frame() %>% rownames_to_column(var="ID.1") %>% 
  left_join(sampled.ID.INFLAM, by=c("ID.1" = "ID")) %>% dplyr::select(-ID.1) -> liver_meta.Inflammation

plasma.meta.1 %>% t() %>% as.data.frame() %>% rownames_to_column(var="ID.1") %>% 
  left_join(sampled.ID.Steatosis, by=c("ID.1" = "ID")) %>% dplyr::select(-ID.1) -> plasma_meta.S

CM.meta.1 %>% t() %>% as.data.frame() %>% rownames_to_column(var="ID.1") %>% 
  left_join(sampled.ID.Steatosis, by=c("ID.1" = "ID")) %>% dplyr::select(-ID.1) -> CM_meta.S


mRNA.gene.count.1 %>% t() %>% as.data.frame() %>% rownames_to_column(var="ID.1") %>% 
  left_join(sampled.ID.Steatosis, by=c("ID.1" = "ID")) %>% dplyr::select(-ID.1) -> mRNA_count.S

clean_background.strong <- theme(plot.background = element_rect("white"),
                                panel.background = element_rect("white"),
                                panel.grid = element_line("white"),
                                axis.line = element_line("black"),
                                axis.text = element_text(size = 14, color = "black"),
                                axis.title = element_text(color = "black"),
                                axis.text.x = element_text(size=14, face="bold"),
                                axis.text.y = element_text(size=14, face="bold"),
                                legend.text = element_text(size = 14),
                                legend.key = element_rect("white"))
pal <- c("lightsalmon1", "gold1", "palegreen4")

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

