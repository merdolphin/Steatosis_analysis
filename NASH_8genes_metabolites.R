source("/home/lina/scratch/18.NASH/Rscripts/NASH_preprocessing_data.R")

library(ggplot2)
library(superb)
library(ggpubr)
library(ggsci)
library(pathview)
library(ggrepel)

importantGenes <- c( "AKR1B10", "GLYCTK",  "SAT2",    "SMOX",    "OAT",
                     "TKT",     "GLYCTK",  "AK1",     "CYP2C19", "ACOX2")



mRNA.gene.dict %>% rownames_to_column(var = "rowid") %>%
  filter(HGNC.symbol %in% importantGenes) %>% dplyr::select(rowid, HGNC.symbol) %>%
  mutate(rowid = paste0("V", rowid)) -> mRNA.gene.dict.8




mRNA_count.S %>% dplyr::select(any_of(mRNA.gene.dict.8$rowid),Histo_STEATOSIS) %>% 
  gather(key, value, -Histo_STEATOSIS) -> df8


df8 %>% left_join(mRNA.gene.dict.8, by=c("key" = "rowid")) -> df8.1




my_comparisons <- list(c("2","1"))

subset(S2vsS1vsS0FDR0.05, mm_symbol %in% importantGenes) %>% as.data.frame() %>% arrange(log2FoldChange) %>% 
  dplyr::select(mm_symbol)

mm_symbol <- c("CYP2C19","ACOX2", "GLYCTK","SAT2","OAT","TKT","SMOX","AK1","AKR1B10")

subset(S2vsS1vsS0FDR0.05, mm_symbol %in% importantGenes) %>% as.data.frame() %>% arrange(log2FoldChange) -> df8.2


df8.1$HGNC.symbol <- factor(df8.1$HGNC.symbol, levels = mm_symbol)

pdf("try.pdf", width=20, height=4)
ggboxplot(df8.1, x="Histo_STEATOSIS", y="value", fill="Histo_STEATOSIS", color="Histo_STEATOSIS",
          alpha=0.2, outlier.shape = NA, add = "jitter",
          add.params = list(alpha = 0.5),width=0.8) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", 
                     method = "wilcox.test", paired = F, hide.ns = T, show.legend = F, size = 8,
                     vjust = 5, 
                     bracket.size = 0)+
  scale_fill_nejm() +
  scale_color_nejm() +
  facet_wrap(~HGNC.symbol, scales="free_y", ncol=9)+
  clean_background +
  scale_x_discrete(labels=c("S0","S1","S2+S3")) +
  scale_fill_discrete(name = "Steatosis", labels  = c("S0 (n = 37)", "S1 (n = 41)", "S2+S3: (n = 36)")) +
  labs(x = "Steatosis level",
       y = "Normalized level") +  guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 16)
  )

dev.off()




# Define the list of compounds to highlight
compounds <- c(
  "C00183", "C00180", "C00258", "C00001", "C01693", "C00005", "C00002", 
  "C00249", "C00122", "C00003", "C00006", "C00016", "C00009", "C00004", 
  "C00008", "C00345", "C00197", "C00007", "C00162", "C00010", "C00011", 
  "C00074", "C00392", "C00165", "C00024"
)

# Load the pathway data for hsa00561
pathway_data <- pathview("hsa00561", species = "hsa", kegg.dir = NULL, pathway.info = NULL)

# Create a pathway view with highlighted compounds
pv.out <- pathview(
  pathway.id = "hsa00561",
  gene.data = importantGenes,
  metabolite.data = NULL,
  compound.highlight = compounds
)


# View the pathway
print(pv.out)
BiocManager::install("clusterProfiler")
library(clusterProfiler)
dev.off()

load(file="/home/lina/scratch/18.NASH/S2vsS1vsS0FDR0.05.rda")

S2vsS1vsS0FDR0.05 %>% as.data.frame() %>% dplyr::select(mm_symbol)


subset(S2vsS1vsS0FDR0.05, mm_symbol %in% importantGenes) %>% as.data.frame() %>% arrange(baseMean)

subset(S2vsS1vsS0FDR0.05, mm_symbol %in% importantGenes) %>% as.data.frame() %>%
  dplyr::select(log2FoldChange, pvalue, padj, mm_symbol) -> FC.8genes

save(FC.13genes, file="FC.13genes.rda")
FC.8genes %>%
  arrange(desc(log2FoldChange)) %>%
  mutate(mm_symbol = factor(mm_symbol, levels = mm_symbol)) %>%
  mutate(size1 = ifelse(pvalue < 0.05, 0.05, NA)) -> FC.8genes.1

FC.8genes.1 %>% arrange((log2FoldChange)) %>% dplyr::select(mm_symbol)


pdf("try.pdf")
dev.off()
FC.8genes.1%>% ggplot(aes(x = mm_symbol, y = log2FoldChange, fill = log2FoldChange, 
                         label = mm_symbol,
                         size=size1)) +
  geom_col(alpha=0.3) +

  geom_smooth(aes(x=1:8, y=FC.8genes.1$log2FoldChange), method = "loess",
              se = FALSE, color = "red", span=.3, levels=999) +  # Add smooth line
  geom_point(size=7, color="orange") +
  geom_point(color="black") +
  #              size = 5) +  # Add repel labels
  geom_text(aes(x=mm_symbol, y=log2FoldChange),
            hjust =  ifelse(FC.8genes.1$log2FoldChange > 0, -0.2, 1.2), size = 5) +
  labs(y = "log2 Fold Change", x = "") +
  clean_background +
  scale_fill_gradient2("Set1") +
  ylim(-2.7,2.7) +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16)) +
  coord_flip()


