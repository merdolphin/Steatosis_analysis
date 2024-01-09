setwd("/home/lina/scratch/18.NASH/Rscripts")


ClinicalData <- read_xlsx("/home/lina/scratch/18.NASH/1.sample/HumanPatients_230216.xlsx",
                          sheet = "Sheet1")

colnames(ClinicalData) -> tmp1
write.csv(tmp1,"colnames.1.csv")
colnames.1 <- read.table("colnames.3.csv", sep=",")
colnames(ClinicalData) <- colnames.1

ClinicalData %>% dplyr::select(ID, Insulin, `HbA1c(%)`, `HOMA2%B`, `HOMA2%S`, HOMA2IR) %>% mutate(ID=as.character(ID)) -> C.insulin

norm_mRNA_count <- read.delim("../4.NAFLD_Transcriptomics/NAFLD_Transcriptomics/CleanData/gene_level.txt",
                              header=T, as.is=T)

norm_mRNA_count[,1:4] -> mRNA.gene.dict
norm_mRNA_count[,-c(1:4)] -> mRNA.gene.count
names(mRNA.gene.count)  <- gsub(x=names(mRNA.gene.count),pattern="S", replacement = "")



mRNA.gene.count %>% t() %>% as.data.frame() %>% rownames_to_column(var="ID") %>% left_join(C.insulin) %>% 
  dplyr::select(-ID) -> mRNA.gene.count.insulin


mRNA.gene.count.insulin[,colSums(mRNA.gene.count.insulin==0) < 13] -> mRNA.gene.count.insulin.20per 


pcor.test.res <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(pcor.test.res) <- c("mRNAs", "insulin", "estimate", "pvalue")




for (i in 1:18212) {
  mRNA <- colnames(mRNA.gene.count.insulin.20per)[i]
  for (j in 18213:18217) {
    insulin <- colnames(mRNA.gene.count.insulin.20per[j])
  tmp <-
      cor.test(
        mRNA.gene.count.insulin.20per[[insulin]] %>% as.numeric(),
        mRNA.gene.count.insulin.20per[[mRNA]],
        method = "spearman",
        exact = FALSE,
        na.action = na.omit
      )
    tmp.df <- data.frame(
      insulin  = "insulin",
      mRNAs = mRNA,
      estimate = as.numeric(tmp$estimate),
      pvalue = tmp$p.value
    )
  }
  pcor.test.res <- rbind(pcor.test.res, tmp.df)
}

pcor.test.res %>% arrange(desc(abs(estimate)))
pcor.test.res %>% arrange(desc(abs(estimate))) %>% dplyr::filter(pvalue < 0.05 & abs(estimate) > 0.4)

pcor.test.res %>% arrange(desc(abs(estimate))) %>% dplyr::filter(pvalue < 0.05 & abs(estimate) > 0.4) %>%
  dplyr::select(mRNAs)  -> mRNAinsulin.8


# Convert the modified string into a vector of individual values
selected_row_ids <- unlist(strsplit(gsub("V", "", mRNAinsulin.8$mRNAs), ", "))


result.mRNA <- mRNA.gene.dict[selected_row_ids, ]

result.mRNA %>% dplyr::select(HGNC.symbol) %>% na.omit() %>% as.vector() 


liver.meta.count <- read.delim("../3.metabolimic/FromChoi/Normalized_data/liver_data_290523.txt",
                               header=T, as.is=T)
liver.meta.count[,1:30] -> liver.meta.dict
liver.meta.count[,-c(1:30)] -> liver.meta.count
names(liver.meta.count) <- gsub(x=names(liver.meta.count), pattern="X", replacement = "")

liver.meta.count %>% t() %>% as.data.frame() %>% rownames_to_column(var="ID") %>% left_join(C.insulin) %>% 
  dplyr::select(-ID)  -> liver.meta.count.insulin

pcor.test.res <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(pcor.test.res) <- c("mRNAs", "insulin", "estimate", "pvalue")



for (i in 1:1352) {
  for (j in 1353:1357) {
    insulin <- colnames(liver.meta.count.insulin[j])
    mRNA <- colnames(liver.meta.count.insulin)[i]
    tmp <-
      cor.test(liver.meta.count.insulin[[insulin]] %>% as.numeric(),
               liver.meta.count.insulin[[mRNA]],
               method = "spearman",
               na.action = na.omit)
    tmp.df <- data.frame(
      insulin  = "insulin",
      mRNAs = mRNA,
      estimate = as.numeric(tmp$estimate),
      pvalue = tmp$p.value
    )
  }
  pcor.test.res <- rbind(pcor.test.res, tmp.df)
}

pcor.test.res %>% arrange(desc(abs(estimate))) %>% dplyr::filter(pvalue < 0.05 & abs(estimate) > 0.4)

pcor.test.res %>% arrange(desc(abs(estimate))) %>% dplyr::filter(pvalue < 0.05 & abs(estimate) > 0.3) %>%
  dplyr::select(mRNAs)  -> liverInsulin.4

liverInsulin.4
# Convert the modified string into a vector of individual values
selected_row_ids <- unlist(strsplit(gsub("V", "", liverInsulin.4$mRNAs), ", "))

result <- liver.meta.dict[selected_row_ids, ]

result %>% dplyr::select(Compound)
result %>% dplyr::select(abbrev..LMSD.)
result %>% dplyr::select(accession..HMDB.)

plasma.meta.count <- read.delim("../3.metabolimic/FromChoi/Normalized_data/plasma_data_290523.txt",
                                header=T, as.is=T)
plasma.meta.count[,1:30] -> plasma.meta.dict
plasma.meta.count[,-c(1:30)] -> plasma.meta.count
names(plasma.meta.count) <- gsub(x=names(plasma.meta.count), pattern="X", replacement = "")

plasma.meta.count %>% t() %>% as.data.frame() %>% rownames_to_column(var="ID") %>% left_join(C.insulin) %>% 
  dplyr::select(-ID) -> plasma.meta.count.insulin

pcor.test.res <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(pcor.test.res) <- c("mRNAs", "insulin", "estimate", "pvalue")


for (i in 1:1087) {
  for (j in 1088:1092) {
    insulin <- colnames(plasma.meta.count.insulin[j])
    mRNA <- colnames(plasma.meta.count.insulin)[i]
    tmp <-
      cor.test(
        plasma.meta.count.insulin[[insulin]] %>% as.numeric(),
        plasma.meta.count.insulin[[mRNA]],
        method = "spearman",
        na.action = na.omit
      )
    tmp.df <- data.frame(
      insulin  = "insulin",
      mRNAs = mRNA,
      estimate = as.numeric(tmp$estimate),
      pvalue = tmp$p.value
    )
  }
  pcor.test.res <- rbind(pcor.test.res, tmp.df)
}

pcor.test.res %>% arrange(desc(abs(estimate))) %>% dplyr::filter(pvalue < 0.05 & abs(estimate) > 0.4)

pcor.test.res %>% arrange(desc(abs(estimate))) %>% dplyr::filter(pvalue < 0.05 & abs(estimate) > 0.3) %>%
  dplyr::select(mRNAs)  -> plasma.Insulin.4

# Convert the modified string into a vector of individual values
selected_row_ids <- unlist(strsplit(gsub("V", "", plasma.Insulin.4$mRNAs), ", "))


result.plasma <- liver.meta.dict[selected_row_ids, ]

result.plasma %>% dplyr::select(Compound)
result.plasma %>% dplyr::select(abbrev..LMSD.)

