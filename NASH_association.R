source("/home/lina/scratch/18.NASH/Rscripts/NASH_preprocessing_data.R")

library(reshape2)
library(Hmisc)
library(doParallel)
library(foreach)
library(stringr)
library(alluvial)
library(ggalluvial)
library(ggsci)


##############
######### 
############
mRNA_count.S %>% dim()

liver_meta.S %>% dim()

load(file = "/home/lina/scratch/18.NASH/mRNA_kruskal_005.rda")
load(file = "/home/lina/scratch/18.NASH/Rscripts/liver.meta.AD0.05.rda")
liver.meta.AD0.05 
mRNA_kruskal_005 %>% length()

mRNA_count.S %>% dplyr::select(all_of(mRNA_kruskal_005)) -> mRNA.DA
liver_meta.S %>% dplyr::select(all_of(liver.meta.AD0.05$key)) -> liver.DA


cor(mRNA.DA, liver.DA, method="spearman") %>% melt() %>% setNames(c("mRNA", "liver", "value")) -> mRNA.liver.sperman.cor


gsub("V","",liver.meta.AD0.05$key) -> liver.meta.AD0.05.key
liver.meta.dict %>% colnames()
liver.meta.dict   %>% rownames_to_column(var="key") %>% 
  filter(key %in% liver.meta.AD0.05.key) %>% 
  dplyr::select(key, Compound, super.class..HMDB.,accession..HMDB.) ->liver.AD.dict



liver.meta.dict   %>% rownames_to_column(var="key") %>% 
  filter(key %in% liver.meta.AD0.05.key) %>% str()


gsub("V","",mRNA_kruskal_005) -> mRNA_kruskal_005.key 

mRNA.gene.dict %>% rownames_to_column(var="key") %>% 
  filter(key %in% mRNA_kruskal_005.key) %>% dplyr::select(key, HGNC.symbol) -> mRNA.AD.dict
mRNA.AD.dict %>% head()
liver.AD.dict %>% head()
mRNA.liver.sperman.cor %>% filter(abs(value)>0.5) %>% mutate(liver=str_replace(liver, "V","")) %>%
  mutate(mRNA=str_replace(mRNA,"V","")) %>% left_join(mRNA.AD.dict, by=c("mRNA"="key")) %>%
  left_join(liver.AD.dict, by=c("liver" = "key")) %>% dplyr::select(HGNC.symbol) %>% distinct()

mRNA.liver.sperman.cor %>% filter(abs(value) > 0.5) %>% dim()



mRNA.liver.sperman.cor %>% filter(abs(value)>0.5) %>% mutate(liver=str_replace(liver, "V","")) %>%
  mutate(mRNA=str_replace(mRNA,"V","")) %>% left_join(mRNA.AD.dict, by=c("mRNA"="key")) %>%
  left_join(liver.AD.dict, by=c("liver" = "key"))%>%
  dplyr::select(super.class..HMDB.) %>% 
  mutate(super.class..HMDB. = strsplit(super.class..HMDB., " --- ")) %>%
  unnest(super.class..HMDB.) %>%
  group_by(super.class..HMDB.) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) 


  
  mRNA.liver.sperman.cor %>% filter(abs(value)>0.5) %>% mutate(liver=str_replace(liver, "V","")) %>%
    mutate(mRNA=str_replace(mRNA,"V","")) %>% left_join(mRNA.AD.dict, by=c("mRNA"="key")) %>%
    left_join(liver.AD.dict, by=c("liver" = "key")) %>% dplyr::select(HGNC.symbol) %>% distinct() %>% as.vector()
  
  
  mRNA.liver.sperman.cor %>% filter(abs(value)>0.5) %>% mutate(liver=str_replace(liver, "V","")) %>%
    mutate(mRNA=str_replace(mRNA,"V","")) %>% left_join(mRNA.AD.dict, by=c("mRNA"="key")) %>%
    left_join(liver.AD.dict, by=c("liver" = "key")) %>% dplyr::select(HGNC.symbol, accession..HMDB.) %>%
    dplyr::select(accession..HMDB.) %>% distinct() %>% as.vector()
  
  
  
  gene.pathway <- read_xlsx("Analysis_results/metascape_AD_association/metascape_result.t6iw70hgv.xlsx",
                            sheet="Enrichment")
  
  gene.pathway %>% dplyr::select(Term, Symbols) -> gene.pathway.1
  geneEnrichment <- read_xlsx("Analysis_results/metascape_AD_association/metascape_result.t6iw70hgv.xlsx", 
                              col_types = "text")
  geneEnrichment %>% dplyr::select(c(1, contains("GO:"))) %>% pivot_longer(-MyList) %>% 
    filter(value == 1) -> geneEnrichment.1

  geneEnrichment.1 %>% head()
  help(read_xlsx)
  geneEnrichment.1$name
  
  mRNA.liver.sperman.cor %>% filter(abs(value)>0.5) %>% mutate(liver=str_replace(liver, "V","")) %>%
    mutate(mRNA=str_replace(mRNA,"V","")) %>% left_join(mRNA.AD.dict, by=c("mRNA"="key")) %>%
    left_join(liver.AD.dict, by=c("liver" = "key")) %>% 
    left_join(geneEnrichment.1, by=c("HGNC.symbol" = "MyList")) %>%
    dplyr::select(HGNC.symbol, name, Compound, super.class..HMDB.) %>% na.omit() %>% 
    mutate(super.class..HMDB. = strsplit(super.class..HMDB., " --- "))  %>%
    unnest(super.class..HMDB.) -> AD.enrichedPathway.super.class

  AD.enrichedPathway.super.class %>% dplyr::select(HGNC.symbol) %>% distinct() %>% as.data.frame()
  AD.enrichedPathway.super.class %>% dplyr::select() %>% distinct() %>% as.data.frame()
  
  
  AD.enrichedPathway.super.class %>% group_by(name, super.class..HMDB.,HGNC.symbol, Compound) %>% 
    summarise(Freq = n()) %>% mutate(GO_number = str_extract(name, "GO:\\d+"))  -> 
    AD.enrichedPathway.super.class.2
  
  
  AD.enrichedPathway.super.class %>% group_by(name, super.class..HMDB.,HGNC.symbol) %>% 
    summarise(Freq = n()) %>% mutate(GO_number = str_extract(name, "GO:\\d+"))  -> 
    AD.enrichedPathway.super.class.1
  
  AD.enrichedPathway.super.class.1 %>% filter(Freq > 4) %>% ungroup %>% dplyr::select(name) %>% unique()
  
  
  pdf("RepresentDAgenesAlluviumGram.pdf", width=15, height=10) 
  AD.enrichedPathway.super.class.1 %>% filter(Freq > 4) %>% ggplot(
    aes(y = Freq, axis1=HGNC.symbol, axis2 = GO_number, 
        axis3 = super.class..HMDB.)) +
    geom_alluvium(aes(fill = name), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Gene","Enriched Pathways", "HMDB super class" ), expand = c(.05, 0.5)) +
    xlim(c(0.9,3.2)) +
    scale_fill_brewer(palette = "Set3") +
    clean_background.strong +
    ggtitle("")
 dev.off()
   
 ######## None of these mutation are significantly different among three Steatosis
  geneMutation <-c("PNPLA3","TM6SF2","APOE","GCKR","GPAM","PNPLA2","TMC4",
                   "MTARC1","APOH","ADH1B","HFE","SERPINA1","TM6SF2","NCAN","APOC3","FGD5","CITED2")
  
  
  AD.enrichedPathway.super.class.2 %>% filter(HGNC.symbol %in% geneMutation)
    
###################################################################

DEseq2FromChoiS2vsS1FDR0.05 <- c("ACOX2","AGTRAP","AK1","AKR1B10","APOLD1","ARPC1B","ASAH1","ATP6","B3GAT3","C2CD4A","CAP1","CD151","CD68","CFL1","CHI3L1","CYB561A3","CYB5R3","CYP2C19","CYTB","DENND2D","DYNLT1","ELOVL1","FADS1","FST","GALNT1","GIPC1","GLYCTK","GOLM1","GPRC5A","GRN","HAAO","HES1","HTRA1","IFI30","IFT43","IL10RB","IL32","JPT1","KCNN2","KIAA0930","LAPTM4B","LAPTM5","LOC105372432","LYZ","M6PR","MCUB","MFSD12","MSN","MVP","ND2","ND3","NFIL3","OAT","OR2I1P","PLA2G7","PLIN3","PPIL1","RASL11A","RRAD","RSRC2","RTN4","S100A9","SAT2","SFN","SLC25A4","SLC27A4","SMOX","SNX5","SPI1","SULF2","TCIRG1","TGFBI","TKT","TMEM87B","TRIM8","TWF2","TYMP")  
  
  AD.enrichedPathway.super.class.2 %>% filter(HGNC.symbol %in% DEseq2FromChoiS2vsS1FDR0.05) -> AD.enrichedPathway.super.class.3
  
  AD.enrichedPathway.super.class.3 %>% ungroup %>% dplyr::select(HGNC.symbol) %>% unique()
  AD.enrichedPathway.super.class.3$name
  

  AD.enrichedPathway.super.class.3 %>% filter(Freq >= 1) %>% ggplot(
                                            aes(y = Freq, axis1=HGNC.symbol, axis2 = GO_number, 
                                                axis3 = super.class..HMDB.,
                                                axis4 = Compound)) +
  geom_alluvium(aes(fill = name), width = 1/12) +
    scale_fill_brewer(palette = "Set3") +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  #scale_x_discrete(name="", limits = c("Gene","Enriched Pathways", "HMDB super class", "Compounds"), 
  #               expand = c(.05, 0.5)) +
    xlim(c(0.7,4.5)) +
 # scale_fill_npg() +
      clean_background.strong +
  ggtitle("")
  dev.off()


  ############ cor.test 
  
  # Set up parallel processing
  num_cores <- detectCores()
  cl <- makeCluster(124)
  registerDoParallel(cl)
  
  # Create an empty result data frame
  pcor.test.res <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(pcor.test.res) <- c("mRNA", "liver", "estimate", "pvalue")
  
  # Perform parallelized loop
  pcor.test.res <- foreach(i = 1:ncol(mRNA.DA), .combine = rbind) %dopar% {
    result_list <- list()
    for (j in 1:ncol(liver.DA)) {
      tmp <- cor.test(mRNA.DA[, i], liver.DA[, j], method = "spearman")
      tmp.df <- data.frame(
        mRNA = colnames(mRNA.DA)[i],
        liver = colnames(liver.DA)[j],
        estimate = tmp$estimate,
        pvalue = tmp$p.value
      )
      result_list[[j]] <- tmp.df
    }
    return(do.call(rbind, result_list))
  }
  
  # Stop the parallel processing cluster
  stopCluster(cl)
  
  # Print the result
  print(pcor.test.res)

  save(pcor.test.res, file="pcor.test.mRNAvs.Liver.res.rda")

  pcor.test.res %>% filter(pvalue < 0.01 & abs(estimate) > 0.5 ) %>% dim()
  

  pcor.test.res %>% arrange(desc(abs(estimate))) -> pcor.test.result.1
  
  pcor.test.result.1 %>% 
    dplyr::select(pvalue) %>% as.matrix() %>% p.adjust(method = "fdr") %>% as.data.frame() -> fdr.result
  
  fdr.result %>% tail()
  
  
  fdr.result %>% 
    mutate_if(is.numeric, ~ 1 * (. < 0.01))  -> fdr.result.1
  
  fdr.result.1[fdr.result.1 == 0] <- NA
  
  colnames(fdr.result.1) <- "fdr.res"
  
  cbind(pcor.test.result.1, fdr.result.1) -> pcor.test.result.all
  
  pcor.test.result.all %>% filter(fdr.res == 1 & abs(estimate) > 0.5) %>% dim()
  
  
