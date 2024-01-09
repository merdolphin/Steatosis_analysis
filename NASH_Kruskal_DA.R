
source("/opt2/scratch/18.NASH/Rscripts/NASH_preprocessing_data.R")
library(tidyr)
library(broom)


###############
##########
############ inflammation
###########
############################


liver_meta.Inflammation %>% gather(key, value, -Histo_INFLAM) %>% group_by(key) %>%
  do(tidy(kruskal.test(x= .$value, g = .$Histo_INFLAM)))  -> liver.meta.Kruskal.res.inflammation

cor(as.numeric(liver_meta.S$Histo_STEATOSIS), as.numeric(liver_meta.Inflammation$Histo_INFLAM), method="spearman")

liver.meta.Kruskal.res.inflammation %>% dplyr::select(key, statistic, p.value) %>% filter(p.value < 0.05) -> liver.meta.AD0.05

liver.meta.AD0.05 %>% dim()

gsub("V","",liver.meta.AD0.05$key) -> liverMetaAD005

liver.meta.dict  %>% dplyr::select(Batch, super.class..HMDB.) %>% rownames_to_column(var="key") %>%
  filter(key %in% liverMetaAD005) %>% dplyr::select(super.class..HMDB.)  %>% 
  dplyr::select(super.class..HMDB.) %>%
  mutate(super.class..HMDB. = strsplit(super.class..HMDB., " --- ")) %>%
  unnest(super.class..HMDB.) %>%
  group_by(super.class..HMDB.) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) 

liver.meta.dict$Batch %>% as.factor() %>% summary()

liver.meta.dict  %>% dplyr::select(Batch, super.class..HMDB.) %>% rownames_to_column(var="key")  %>%
  filter(key %in% liverMetaAD005) %>% dplyr::select(Batch)  %>% 
  # mutate(Batch = strsplit(Batch, ",")) %>%
  #  unnest(Batch) %>%
  group_by(Batch) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) 



############
##########
############ Steatosis
###########
############################


shapiro.test.res <- lapply(liver.meta.count, shapiro.test)
shapiro.test.res %>% str()

res1 <- sapply(shapiro.test.res, `[`,c("stasticic","p.value"))

t(res1)
liver_meta.S$Histo_STEATOSIS 

mRNA_count.S %>% filter(Histo_STEATOSIS != 0) %>% gather(key, value, -Histo_STEATOSIS) %>% group_by(key) %>%
  do(tidy(wilcox.test(x= .$value, g = .$Histo_STEATOSIS)))  -> mRNA.count.wilcox.res

mRNA.count.wilcox.res %>% dplyr::select(key, statistic, p.value) %>% filter(p.value < 0.05) -> mRNA.wilcox.AD0.05

save(mRNA.wilcox.AD0.05,file="mRNA.wilcox.AD0.05.rda")


liver_meta.S %>% filter(Histo_STEATOSIS != 0) %>% gather(key, value, -Histo_STEATOSIS) %>% group_by(key) %>%
  do(tidy(wilcox.test(x= .$value, g = .$Histo_STEATOSIS)))  -> liver.meta.wilcox.res

liver.meta.wilcox.res %>% dplyr::select(key, statistic, p.value) %>% filter(p.value < 0.05) -> liver.wilcox.AD0.05

save(liver.wilcox.AD0.05,file="liver.wilcox.AD0.05.rda")

gsub("V","",liver.meta.AD0.05$key) -> liverMetaAD005



liver.wilcox.AD0.05 %>% arrange(desc(statistic))


mRNA_count.S %>% gather(key, value, -Histo_STEATOSIS) %>% group_by(key) %>%
  do(tidy(kruskal.test(x= .$value, g = .$Histo_STEATOSIS)))  -> mRNA.count.Kruskal.res
  

mRNA.count.Kruskal.res %>% dplyr::select(key, statistic, p.value) %>% filter(p.value < 0.05) -> mRNA.AD0.05

#gsub("V", "", mRNA.AD0.05$key) -> mRNA_kruskal_005
mRNA.AD0.05$key -> mRNA_kruskal_005
save(mRNA_kruskal_005, file= "mRNA_kruskal_005.rda")


gsub("V","",mRNA_kruskal_005) -> mRNA_kruskal_005.key 
mRNA_kruskal_005.key
mRNA.gene.dict[mRNA_kruskal_005.key,] %>% dplyr::select(HGNC.symbol) -> AD.gene.list.5485

write.table(AD.gene.list.5485, file = "AD.gene.list.5485.txt")


PPIL1interactionGenes <- c("AQR","BUD31","CDC40","CDC5L","EFTUD2","PPIL1","PRPF19","PRPF8","RBM22","SNW1","SYF2")
CFL1interactionGenes <- c("ACTB","ACTG1","CFL1","CTTN","HCLS1","LIMK1","LIMK2","PFN1","SSH1","WASL","WDR1")

intersect(AD.gene.list.5485$HGNC.symbol, PPIL1interactionGenes) 
intersect(AD.gene.list.5485$HGNC.symbol, CFL1interactionGenes)

#####3 liver 

liver_meta.S %>% gather(key, value, -Histo_STEATOSIS) %>% group_by(key) %>%
  do(tidy(kruskal.test(x= .$value, g = .$Histo_STEATOSIS)))  -> liver.meta.Kruskal.res

liver.meta.Kruskal.res %>% dplyr::select(key, statistic, p.value) %>% 
  filter(p.value < 0.05) -> liver.meta.AD0.05

liver.meta.AD0.05 %>% arrange(desc(statistic)) %>% as.data.frame()
mRNA.AD0.05 %>% arrange(desc(statistic)) %>% head(n=200) %>% as.data.frame()  -> top200
gsub("V","", top200$key) -> top200key
top200key %>% head()
mRNA.gene.dict[top200key,] %>% dplyr::select(HGNC.symbol)


gsub("V","",liver.meta.AD0.05$key) -> liverMetaAD005

save(liver.meta.AD0.05, file="liver.meta.AD0.05.rda")
liver.meta.AD0.05$key

liver.meta.dict  %>% dplyr::select(Batch, super.class..HMDB.) %>% rownames_to_column(var="key")  %>%
  filter(key %in% liverMetaAD005) %>%dplyr::select(super.class..HMDB.)  %>% 
  dplyr::select(super.class..HMDB.) %>%
  mutate(super.class..HMDB. = strsplit(super.class..HMDB., " --- ")) %>%
  unnest(super.class..HMDB.) %>%
  group_by(super.class..HMDB.) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) 

liver.meta.dict  %>% dplyr::select(Batch, super.class..HMDB.) %>% rownames_to_column(var="key")  %>%
  filter(key %in% liverMetaAD005) %>% dplyr::select(Batch)  %>% 
 # mutate(Batch = strsplit(Batch, ",")) %>%
#  unnest(Batch) %>%
  group_by(Batch) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) 

