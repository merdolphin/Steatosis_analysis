source("/home/lina/scratch/18.NASH/Rscripts/NASH_preprocessing_data.R")

library(ggplot2)
library(superb)
library(ggpubr)
library(ggsci)
library("dgof")

importantGenes <- c("LPCAT1","GPCPD1","CDS2","LPGAT1","CYP2C8","CYP3A4","VNN1",
                    "PPCS","BAAT","GPX1","GPX2","MGST1","RRM2B")



mRNA.gene.dict %>% rownames_to_column(var = "rowid") %>%
  filter(HGNC.symbol %in% importantGenes) %>% dplyr::select(rowid, HGNC.symbol) %>%
  mutate(rowid = paste0("V", rowid)) -> mRNA.gene.dict.13



mRNA_count.S %>% dplyr::select(any_of(mRNA.gene.dict.13$rowid),Histo_STEATOSIS) %>% 
  gather(key, value, -Histo_STEATOSIS) -> df13


df13 %>% left_join(mRNA.gene.dict.13, by=c("key" = "rowid")) -> df13.1




df13.1 %>% ggplot(aes(x=HGNC.symbol, y=value, fill=as.factor(Histo_STEATOSIS), group=Histo_STEATOSIS))+
  geom_col(position = "dodge") +
  clean_background

list(c("3\",\"3")) %>% unlist()
my_comparisons <- list(c("1","2"), c("1","0"),c("2","0"))

df13.1 %>% filter(! HGNC.symbol == "GPCPD1") -> df12.1

df12.1$HGNC.symbol <- factor(df12.1$HGNC.symbol, 
                                   levels = c("CYP3A4","CYP2C8","BAAT","GPCPD1","PPCS","CDS2","MGST1","LPCAT1","RRM2B","VNN1","GPX1","LPGAT1","GPX2"))

pdf("try.pdf", width=10, height=3)
ggboxplot(df12.1, x="Histo_STEATOSIS", y="value", fill="Histo_STEATOSIS", color="Histo_STEATOSIS",
          alpha=0.2, outlier.shape = NA, add = "jitter",
          add.params = list(alpha = 0.5),width=0.8) +
  #stat_compare_means(label.y=min(df13$value))+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", 
                     method = "wilcox.test", paired = F, hide.ns = T, show.legend = F, size = 8)+
  scale_fill_nejm() +
  scale_color_nejm() +
  facet_wrap(~HGNC.symbol, scales="free_y", ncol=6)+
  theme_minimal() +
  scale_x_discrete(labels=c("S0","S1","S2+S3")) +
  scale_fill_discrete(name = "Steatosis", labels  = c("S0 (n = 37)", "S1 (n = 41)", "S2+S3: (n = 36)")) +
  labs(x = "Steatosis level",
       y = "Normalized level") +  guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6),
        legend.position = "none",
        strip.background = element_blank()
  )

dev.off()




ggboxplot(df13, x = "Histo_STEATOSIS", y = "value", fill = "Histo_STEATOSIS", color = "Histo_STEATOSIS",
          alpha = 0.2, outlier.shape = NA, add = "jitter",
          add.params = list(alpha = 0.5), width = 0.8) +
  scale_fill_nejm() +
  scale_color_nejm() +
  facet_wrap(~ key , scales = "free_y") +
  clean_background.strong +
  scale_x_discrete(labels = c("S0", "S1", "S2+S3")) +
  scale_fill_discrete(name = "Steatosis", labels = c("S0 (n = 37)", "S1 (n = 41)", "S2+S3: (n = 36)")) +
  labs(x = "Steatosis level", y = "Normalized level") +
  guides(color = FALSE, alpha = FALSE) +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom")


ggplot(mpg, aes(class, hwy)) +
  geom_boxplot() +
  geom_signif(
    comparisons = list(
      c("compact", "pickup"),
      c("subcompact", "suv")
    ),
    map_signif_level = function(p) sprintf("p = %.2g", p)
  )



df13 %>% gather(key, value, -Histo_STEATOSIS) %>% 
  ggplot(aes(x=Histo_STEATOSIS, y=value, color=Histo_STEATOSIS)) +
  geom_boxplot() + 
  facet_wrap(~key, scales="free_y")+
  clean_background



