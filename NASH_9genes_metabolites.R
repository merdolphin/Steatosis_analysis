source("/home/lina/scratch/18.NASH/Rscripts/NASH_preprocessing_data.R")

library(ggplot2)
library(superb)
library(ggpubr)
library(ggsci)
library(KEGGREST)
library(tidyverse)
library(org.Hs.eg.db)

library("dgof")



integratedPathway <- read.csv(file = "Analysis_results/metaboanalysis_Deseq2FDR0.05_association_metabolite/MetaboAnalyst_result_pathway.csv")

integratedPathway$X.log10.p. %>% summary()

mRNA_liver_pathway.p <- integratedPathway %>% ggplot(aes(x=Impact, y=X.log10.p., label=X))+
  geom_point(aes(size=Impact*20, fill=X.log10.p.), shape=21) +
  scale_fill_gradient2(midpoint=0.46, 
                       low="orange", 
                       mid="yellow",
                       high="red", 
                       space ="Lab",
                       breaks = c(0.5,1.0,1.5,2,2.5))+
  geom_text_repel(
    data = subset(integratedPathway, Impact > 0.10 & X.log10.p. > 0.88),  # Filter points with Impact > 4
    box.padding = 1, 
    point.padding = 0.5,
    segment.color = "grey",
    segment.size = 0.3,
    aes(size=6)
  ) +
  labs(title = "Integration of Transcriptomic and Metabolic Pathway Analysis",
       x = "Pathway Impact",
       y = "-log10(p-value)") +
  scale_y_continuous(breaks = seq(0.5,4, 0.5)) +  # Exclude 0 from y-axis ticks
  clean_background   

mRNA_liver_pathway.p 

pdf("mRNA_liver_pathway.pdf", width=10, height=7)
mRNA_liver_pathway.p
dev.off()

# Gene symbols for the provided gene IDs
gene_ids <- c("57016","132158","1557","427","112483","54498","4942","7086","132158","203")

# Convert Entrez Gene IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, keytype = "ENTREZID", column = "SYMBOL")

# Display the results
result_df <- data.frame(Entrez_Gene_ID = gene_ids, Gene_Symbol = gene_symbols)
print(result_df)


importantGenes <- result_df$Gene_Symbol



mRNA.gene.dict %>% rownames_to_column(var = "rowid") %>%
  filter(HGNC.symbol %in% importantGenes) %>% dplyr::select(rowid, HGNC.symbol) %>%
  mutate(rowid = paste0("V", rowid)) -> mRNA.gene.dict.8



mRNA_count.S %>% dplyr::select(any_of(mRNA.gene.dict.8$rowid),Histo_STEATOSIS) %>% 
  gather(key, value, -Histo_STEATOSIS) -> df8



df8 %>% left_join(mRNA.gene.dict.8, by=c("key" = "rowid")) -> df8.1



my_comparisons <- list(c("2","1"))

subset(S2vsS1vsS0FDR0.05, mm_symbol %in% importantGenes) %>% as.data.frame() %>% arrange(log2FoldChange) %>% 
  dplyr::select(mm_symbol)

mm_symbol <- c("CYP2C19", "GLYCTK","SAT2","OAT","TKT", "ASAH1", "SMOX","AK1","AKR1B10")

subset(S2vsS1vsS0FDR0.05, mm_symbol %in% importantGenes) %>% as.data.frame() %>% arrange(log2FoldChange) -> df8.2

df8.1 %>% str()
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



load(file="/home/lina/scratch/18.NASH/S2vsS1vsS0FDR0.05.rda")

S2vsS1vsS0FDR0.05 %>% as.data.frame() %>% dplyr::select(mm_symbol)


subset(S2vsS1vsS0FDR0.05, mm_symbol %in% importantGenes) %>% as.data.frame() %>% arrange(baseMean)

subset(S2vsS1vsS0FDR0.05, mm_symbol %in% importantGenes) %>% as.data.frame() %>%
  dplyr::select(log2FoldChange, pvalue, padj, mm_symbol) -> FC.8genes

FC.8genes %>%
  arrange(desc(log2FoldChange)) %>%
  mutate(mm_symbol = factor(mm_symbol, levels = mm_symbol)) %>%
  mutate(size1 = ifelse(pvalue < 0.05, 0.05, NA)) -> FC.8genes.1

FC.8genes.1 %>% arrange((log2FoldChange)) %>% dplyr::select(mm_symbol)


FC.8genes.1%>% ggplot(aes(x = mm_symbol, y = log2FoldChange, fill = log2FoldChange, 
                          label = mm_symbol,
                          size=size1)) +
  geom_col(alpha=0.3) +
  #geom_line(aes(group = 1), color = "blue") +
  geom_smooth(aes(x=1:9, y=FC.8genes.1$log2FoldChange), method = "loess",
              se = FALSE, color = "red", span=.3, levels=999) +  # Add smooth line
  #geom_text_repel(box.padding = 0.8, segment.size = 0.2, 
  #               hjust = ifelse(FC.13genes.1$log2FoldChange > 0, 0.1, -0.1),
  geom_point(size=7, color="orange") +
  geom_point(color="black") +
  #              size = 5) +  # Add repel labels
  geom_text(aes(x=mm_symbol, y=log2FoldChange),
            hjust =  ifelse(FC.8genes.1$log2FoldChange > 0, -0.2, 1.2), size = 5) +
  labs(y = "log2 Fold Change", x = "") +
  #guides(color = FALSE, alpha = FALSE) +
  clean_background +
  scale_fill_gradient2("Set1") +
  ylim(-2.7,2.7) +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16)) +
  coord_flip()



