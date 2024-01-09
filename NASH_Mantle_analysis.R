source("~/scratch/18.NASH/Rscripts/NASH_preprocessing_data.R")

# Calculate spearman correlation
mantle_test <- function(){
  spearman_cor_mRNA <- cor(mRNA.gene.count.1, method = "spearman")
  spearman_cor_liver <- cor(liver_meta.1, method="spearman")
  spearman_cor_plasma <- cor(plasma_meta.1, method="spearman")
  spearman_cor_CM <- cor(CM_meta.1, method="spearman")
  
  
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

########### draw venna grame


metabolites_CM <- c("Metabolite1", "Metabolite2", "Metabolite3", ...)  # Replace with CM metabolites
metabolites_liver <- c("Metabolite2", "Metabolite4", "Metabolite5", ...)  # Replace with liver metabolites
metabolites_plasma <- c("Metabolite1", "Metabolite3", "Metabolite5", ...)  # Replace with plasma metabolites

# Create a list of sets
venn_list <- list(
  CM = metabolites_CM,
  Liver = metabolites_liver,
  Plasma = metabolites_plasma
)

# Create Venn diagram
venn.plot <- venn.plot <- venn.diagram(
  x = venn_list,
  category.names = c("CM", "Liver", "Plasma"),
  filename = NULL,
  output = TRUE
)

# Display the diagram
grid.draw(venn.plot)

# Calculate pearson correlation
pearson_cor_mRNA <- cor(mRNA.gene.count.1, method = "pearson")
pearson_cor_liver <- cor(liver_meta.1, method="pearson")
pearson_cor_plasma <- cor(plasma_meta.1,  method="pearson")


mantel_mRNA_liver <- mantel(pearson_cor_mRNA, pearson_cor_liver, method = "pearson", permutations = 999)
mantel_mRNA_liver

mantel_mRNA_plasma <- mantel(pearson_cor_mRNA, pearson_cor_plasma, method = "pearson", permutations = 999)
mantel_mRNA_plasma

mantel_liver_plasma <- mantel(pearson_cor_liver, pearson_cor_plasma, method = "pearson", permutations = 999)
mantel_liver_plasma


# PERMANOVA analysis
#permanova_result <- adonis2(mRNA.gene.count ~ liver.meta.count,  permutations = 999)


# Calculate Cosine similarity

cosine_similarity_mRNA <- similarity(norm_mRNA_count.1, norm_mRNA_count.1, method = "cosine")


####### 

# Standardize your data using Z-score standardization
norm_mRNA_count.1_zscore <- scale(norm_mRNA_count.1)
liver_meta.1_zscore <- scale(liver_meta.1)
plasma_meta.1_zscore <- scale(plasma_meta.1)




# Calculate Euclidean distances using standardized data
euclidean_distance_norm_mRNA <- dist(norm_mRNA_count.1_zscore, norm_mRNA_count.1_zscore, method = "euclidean")
euclidean_distance_norm_plasma <- dist(plasma_meta.1_zscore,plasma_meta.1_zscore, method = "euclidean")
euclidean_distance_norm_liver <- dist(liver_meta.1_zscore,liver_meta.1_zscore, method = "euclidean")


mantel_mRNA_liver <- mantel(euclidean_distance_norm_mRNA, euclidean_distance_norm_liver, method = "pearson", permutations = 999)
mantel_mRNA_liver

mantel_mRNA_plasma <- mantel(euclidean_distance_norm_mRNA, euclidean_distance_norm_plasma, method = "pearson", permutations = 999)
mantel_mRNA_plasma

mantel_liver_plasma <- mantel(euclidean_distance_norm_liver, euclidean_distance_norm_plasma, method = "pearson", permutations = 999)
mantel_liver_plasma



