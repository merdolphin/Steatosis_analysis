
source("NASH_preprocessing_data.R")
library(org.Hs.eg.db)
library(DESeq2)
library(ggplot2)
library(superb)
library(ggpubr)
library(ggsci)

samples_information <- read_xlsx("1.sample/HumanPatients_230216.xlsx")
samples_information %>% dplyr::select(ID, `Histo_STEATOSIS (0-3)`,  `Trig (<1.5)`, `LDL (<3.5 mmol/L)`,`HDL (>1.0 mmol/L = normal)`, 
                                      `TotalChol (3.5-5.5 mM)`, `NonHDL-Chol (mmol/L)`, `ALT(5-30`, `AST (10-35)`,
                                      `GGT (5-35 U/L)`, `Total Bili (3-15)`, `Ferritin (20-500 ng/mL)`,
                                      `HOMA2%B`,`HOMA2%S`,HOMA2IR, HISTO_NAS, 
                                      `Histo_INFLAM(0-3)`,`Histo_BALLOON (0-2)`,`Histo_FIBROSIS (0-4)`, 
                                      `Liver Fibrosis grouping`, 
                                      `Liver Gross grouping`, `Diabetes (0=No, 1=Yes, 2= Ex diabetic)`) %>% 
  setNames(c("ID", "Histo_STEATOSIS","Trig", "LDL", "HDL", "TotalChol", "NonHDL-Chol", 
             "ALT", "AST", "GGT", "TotalBili", "Ferritin", "HOMA2B", "HOMA2S", "HOMA2IR", "HISTO_NAS",
             "Histo_INFLAM","Histo_BALLOON","Histo_FIBROSIS",
             "Liver_Fibrosis_grouping","Liver_Gross_grouping", "Diabetes")) %>%
  mutate(ID=as.character(ID)) %>%
  mutate(Histo_STEATOSIS = as.factor(Histo_STEATOSIS)) -> samples.18



norm_mRNA_count <- read.delim("4.NAFLD_Transcriptomics/NAFLD_Transcriptomics/CleanData/gene_level.txt",
                              header=T, as.is=T)

norm_mRNA_count[,1:4] -> mRNA.gene.dict
norm_mRNA_count[,-c(1:4)] -> mRNA.gene.count
names(mRNA.gene.count)  <- gsub(x=names(mRNA.gene.count),pattern="S", replacement = "")
mRNA.gene.count
rownames(mRNA.gene.count) <- mRNA.gene.dict$Gene.stable.ID

mRNA.gene.count[,colnames(mRNA.gene.count) %in% samples.18$ID] -> mRNA.gene.count.94
mRNA.gene.count.94 %>% rownames()

save(mRNA.gene.count.94, file="Rscripts/SLC19A1/mRNA.gene.count.94.rda")

samples.18[samples.18$ID %in% colnames(mRNA.gene.count.94),] -> sam94

sam94[match(colnames(mRNA.gene.count.94), sam94$ID),] %>%
  mutate(HOMA2B = as.numeric(HOMA2B)) %>%
  mutate(Liver_Gross_grouping = factor(Liver_Gross_grouping, ordered=TRUE))%>%
  mutate(Liver_Fibrosis_grouping = factor(Liver_Fibrosis_grouping, ordered = TRUE)) %>%
  mutate(Histo_INFLAM = factor(Histo_INFLAM, ordered = TRUE)) %>%
  mutate(Histo_BALLOON = factor(Histo_BALLOON)) %>%
  mutate(HISTO_NAS = factor(HISTO_NAS, ordered = TRUE)) %>%
  mutate(Histo_FIBROSIS = factor(Histo_FIBROSIS, ordered = TRUE)) %>%
  mutate(Diabetes=as.factor(Diabetes)) -> sam94.1

sam94.1 %>% str()

mRNA.gene.count.94 %>% colnames() == sam94.1$ID

aggr(sam94.1, plot=F)

mRNA.gene.dict %>% filter(HGNC.symbol == "SLC19A1") 

mRNA.gene.count.94["ENSG00000173638",] %>% t() %>% as.data.frame() %>% rownames_to_column(var="ID") -> SLC19A1.count
SLC19A1.count$ID == sam94.1$ID
SLC19A1.count %>% str()
left_join(SLC19A1.count, sam94.1, by=c("ID"="ID")) %>% dplyr::select(-ID) -> SLC19A1.sam94

samples_information %>% mutate(ID = as.character(ID)) %>% 
  right_join(SLC19A1.count, by=c("ID"="ID")) %>% dplyr::select(-ID) -> SLC19A1.all
SLC19A1.all %>% colnames()
SLC19A1.all[, colSums(is.na(SLC19A1.all)) < 5] -> SLC19A1.all.1
SLC19A1.all.1 %>% colnames()
correlationAll <- cor(SLC19A1.all.1 %>% mutate_if(is.factor, as.numeric) %>% mutate_if(is.character, as.numeric),
    method="spearman",
    use="pairwise.complete.obs")

correlationAll %>% dim()

correlationAll[,44] %>% enframe() %>% arrange(desc(abs(value)))

SLC19A1.sam94 %>% str()

# Basic summary statistics for SLC19A1 levels
summary(SLC19A1.sam94$ENSG00000173638)



# Boxplot of SLC19A1 levels 
ggplot(SLC19A1.sam94, aes(x = Histo_BALLOON, y = ENSG00000173638)) +
  geom_boxplot() +
  labs(title = "Boxplot of SLC19A1 Levels by BALLOON Level",
       x = "BALLOONING staging",
       y = "SLC19A1 Levels")

SLC19A1.sam94 %>% dplyr::select(ENSG00000173638, Histo_BALLOON) %>%
  kruskal.test(ENSG00000173638 ~ Histo_BALLOON) -> kruskal_result

library(dunn.test)
# Post hoc Dunn test
posthoc_dunn <- dunn.test(SLC19A1.sam94$ENSG00000173638, g = SLC19A1.sam94$Histo_BALLOON, method = "bonferroni")

print(kruskal_result)
print(posthoc_dunn)


pdf("Rscripts/SLC19A1/SLC19A1_Ballooning.pdf")
my_comparisons <- list(c("1","0"))
ggboxplot(SLC19A1.sam94, x="Histo_BALLOON", y="ENSG00000173638", fill="Histo_BALLOON", color="Histo_BALLOON",
          alpha=0.5, outlier.shape = NA, add = "jitter",
          add.params = list(alpha = 0.5),width=0.8) +
  stat_compare_means(label.y = 0.5) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", 
                     method = "wilcox.test", paired = F, hide.ns = T, show.legend = F, size = 8)+
  scale_fill_nejm() +
  ylim(0,7) +
  scale_color_nejm() +
  labs(x = "Ballooning",
       y = "SLC19A1 normalized level") +  guides(color=FALSE, alpha = FALSE)+
  clean_background +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))
dev.off()


SLC19A1.sam94 %>% ggplot(aes(x = HISTO_NAS, y = ENSG00000173638)) +
  geom_boxplot() +
  labs(title = "Boxplot of SLC19A1 Levels by NAS Level",
       x = "NAS",
       y = "SLC19A1 Levels")

SLC19A1.sam94 %>% ggplot(aes(x = HOMA2S, y = ENSG00000173638)) +
  geom_col() +
  labs(title = "Boxplot of SLC19A1 Levels by ",
       x = "HOMA2IR",
       y = "SLC19A1 Levels")



# Correlation matrix between SLC19A1 and other numeric variables
correlation_matrix <- cor(SLC19A1.sam94 %>% mutate_if(is.factor, as.numeric),
                          method = "spearman", 
                          use = "pairwise.complete.obs")

correlation_matrix[,1] %>% enframe() %>% arrange(desc(abs(value)))

cor(SLC19A1.sam94$ENSG00000173638, SLC19A1.sam94$Trig)


SLC19A1.sam94$ENSG00000173638_log <- log(SLC19A1.sam94$ENSG00000173638 + 1)


model <- glm(Histo_BALLOON ~ ENSG00000173638, data = SLC19A1.sam94)
model <- glm(ENSG00000173638 ~ Histo_BALLOON + HOMA2IR + HOMA2S + Ferritin, data = SLC19A1.sam94)
summary(model)

library(car)
vif(model) 
model <- lm(ENSG00000173638 ~ poly(Histo_BALLOON, degree = 2) + poly(HISTO_NAS, degree = 2), data = SLC19A1.sam94)

summary(model)

library(glmnet)



model <- step(lm(ENSG00000173638 ~ Histo_BALLOON + HISTO_NAS, data = SLC19A1.sam94 ), direction = "both")

SLC19A1.sam94$Histo_BALLOON


# Load the randomForest package
library(randomForest)
correlationAll
correlationAll[,44] %>% enframe() %>% arrange(desc(abs(value))) %>% head(n=10) %>% 
  dplyr::select(name) -> topSpearman.names

topSpearman.names
SLC19A1.all.1[,topSpearman.names$name] %>% na.omit() -> subset_data

subset_data
# Fit a random forest
rf_model <- randomForest(ENSG00000173638 ~ ., 
                         data = subset_data, importance=TRUE)

print(rf_model)

############ BALLOONING ########

correlationAll %>% colnames()
correlationAll[,38] %>% enframe() %>% arrange(desc(abs(value))) %>% print(n=15)

###########################



sam94.1

dds <- DESeqDataSetFromMatrix(countData = round(mRNA.gene.count.94), colData = sam94.1, 
                              design = ~ Histo_BALLOON)



# Run DESeq2 analysis
dds <- DESeq(dds)
resultsNames(dds)
# Differential expression analysis


convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}


res <- results(dds)
res <- results(dds, contrast = c("Histo_BALLOON","2","1"))
res <- results(dds, contrast = c("Histo_BALLOON","1","0"))
res$mm_symbol <- convertIDs(row.names(res), "ENSEMBL","SYMBOL", org.Hs.eg.db)
res$mm_entrezgene <- convertIDs(row.names(res), "ENSEMBL","ENTREZID", org.Hs.eg.db)

DEGs <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

DEGs <- subset(res, padj < 0.05)
res
DEGs$mm_symbol
