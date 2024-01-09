source("~/scratch/18.NASH/Rscripts/NASH_preprocessing_data.R")


mantle_test <- function(){
  spearman_cor_mRNA <- cor(mRNA.gene.count.1, method = "spearman")
  spearman_cor_liver <- cor(liver.meta.1, method="spearman")
  spearman_cor_plasma <- cor(plasma.meta.1, method="spearman")
  spearman_cor_CM <- cor(CM.meta.1, method="spearman")
  
  
  mantel_mRNA_liver <- mantel(spearman_cor_mRNA, spearman_cor_liver, method = "spearman", permutations = 999)
  mantel_mRNA_liver
  
  mantel_mRNA_plasma <- mantel(spearman_cor_mRNA, spearman_cor_plasma, method = "spearman", permutations = 999)
  mantel_mRNA_plasma
  
  mantel_mRNA_CM <- mantel(spearman_cor_mRNA, spearman_cor_CM, method = "spearman", permutations = 999)
  mantel_mRNA_CM
  
  mantel_liver_plasma <- mantel(spearman_cor_liver, spearman_cor_plasma, method = "spearman", permutations = 999)
  mantel_liver_plasma
  
  mantel_liver_CM <- mantel(spearman_cor_liver, spearman_cor_CM, method = "spearman", permutations = 999)
  mantel_liver_CM
  
  mantel_plasma_CM <- mantel(spearman_cor_plasma, spearman_cor_CM, method = "spearman", permutations = 999)
  mantel_plasma_CM
  
}

liver.meta.1 %>% dim()
CM.meta.1 %>% dim() 

# Calculate pearson correlation
pearson_cor_mRNA <- cor(mRNA.gene.count.1, method = "pearson")
pearson_cor_liver <- cor(liver.meta.1, method="pearson")
pearson_cor_plasma <- cor(plasma.meta.1,  method="pearson")


mantel_mRNA_liver <- mantel(pearson_cor_mRNA, pearson_cor_liver, method = "pearson", permutations = 999)
mantel_mRNA_liver

mantel_mRNA_plasma <- mantel(pearson_cor_mRNA, pearson_cor_plasma, method = "pearson", permutations = 999)
mantel_mRNA_plasma

mantel_liver_plasma <- mantel(pearson_cor_liver, pearson_cor_plasma, method = "pearson", permutations = 999)
mantel_liver_plasma

