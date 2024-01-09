library(doParallel)
library(foreach)

load(file = "S2S1S0_FDR0.05.gene.list.rda")
S2S1S0_FDR0.05.gene.list

load( file="S2vsS1vsS0FDR0.05.rda")
load(file="mRNA.gene.count.94.rda")

S2vsS1vsS0FDR0.05 %>% rownames()

mRNA.gene.count.94[rownames(mRNA.gene.count.94) %in% rownames(S2vsS1vsS0FDR0.05),] -> selected.genesFDR0.05.count


# Find common column names
common_colnames <- intersect(colnames(liver.meta.count), colnames(selected.genesFDR0.05.count))
common_colnames %>% length()
# Subset based on common column names


liver.meta.count.1 <- liver.meta.count[,common_colnames]


selected.genesFDR0.05.count.1 <- selected.genesFDR0.05.count[, common_colnames] %>% t()

liver.meta.count.2 <- liver.meta.count.1[, match(rownames(selected.genesFDR0.05.count.1), colnames(liver.meta.count.1))] %>% t()

# Initialize vectors to store correlation coefficients and p-values

############ cor.test 


cl <- makeCluster(124)
registerDoParallel(cl)


# Create an empty result data frame
pcor.test.res <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(pcor.test.res) <- c("mRNA", "liver", "estimate", "pvalue")

# Perform parallelized loop
pcor.test.res

for (i in 1:ncol(selected.genesFDR0.05.count.1))
{
  for (j in 1:ncol(liver.meta.count.2)) {
    tmp <-
      cor.test(selected.genesFDR0.05.count.1[, i],
               liver.meta.count.2[, j],
               method = "spearman")
    tmp.df <- data.frame(
      mRNA = colnames(selected.genesFDR0.05.count.1)[i],
      liver = j,
      estimate = as.numeric(tmp$estimate),
      pvalue = tmp$p.value
    )
    pcor.test.res <- rbind(pcor.test.res, tmp.df)
  }
}

# Print the result
print(pcor.test.res)

for (i in seq(1,100,5)){
  print(i)
}

save(pcor.test.res, file="pcor.test.mRNAFDR0.05vs.Liver.rda")

pcor.test.res %>% filter(pvalue < 0.01 & abs(estimate) > 0.5 ) %>% dim()

pcor.test.res %>% filter(pvalue < 0.05 & abs(estimate) > 0.5 ) %>% dplyr::select(mRNA) %>% distinct() -> mRNA.transcript

S2vsS1vsS0FDR0.05[rownames(S2vsS1vsS0FDR0.05) %in% mRNA.transcript$mRNA,] %>% as.data.frame() %>%
  dplyr::select(mm_symbol)


pcor.test.res %>% filter(pvalue < 0.05 & abs(estimate) > 0.5 ) %>% dplyr::select(liver) %>% distinct() -> liver.ID
  liver.ID
liver.meta.dict[liver.ID$liver,] %>% dplyr::select(accession..HMDB.) %>% as.vector()


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


