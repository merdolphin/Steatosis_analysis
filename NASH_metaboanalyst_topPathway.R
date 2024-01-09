library(KEGGREST)
library(tidyverse)
library(httr)
library(jsonlite)

jointpa_matched_feature <- read.csv(file = "/home/lina/scratch/18.NASH/Analysis_results/metaboanalyst_AD_0.5/jointpa_matched_features.csv")


compound_data <- data.frame(
  CPD_ID = c(),
  Name = c(""),
  HMDB_ID = c("")
)

# List of compound IDs
compound_ids <- c("C00350", "C00157", "C04230", "C00346", "C00157", "C00819", "C00064", "C00099", "C00245", "C00550", "C00319", "C00195", "C00346")



gene_hsa <- c("hsa:79888", "hsa:56261", "hsa:8760", "hsa:9926", "hsa:1558", "hsa:1576", "hsa:8876", "hsa:79717", "hsa:570", "hsa:2876", "hsa:2877", "hsa:4257", "hsa:50484")
gene_symbol <-c ("ACSL5 (Acyl-CoA Synthetase Long Chain Family Member 5)", " MOGAT2 (Monoacylglycerol O-Acyltransferase 2)", " DGAT1 (Diacylglycerol O-Acyltransferase 1)", " LPCAT1 (Lysophosphatidylcholine Acyltransferase 1)", " CYP2E1 (Cytochrome P450 Family 2 Subfamily E Member 1)", " CYP2C9 (Cytochrome P450 Family 2 Subfamily C Member 9)", " SQLE (Squalene Epoxidase)", " SCD5 (Stearoyl-CoA Desaturase 5)", " BCAT1 (Branched-Chain Amino Acid Transaminase 1)", " GPX1 (Glutathione Peroxidase 1)", " GPX2 (Glutathione Peroxidase 2)", " MGMT (O-6-Methylguanine-DNA Methyltransferase)", " CPT2 (Carnitine Palmitoyltransferase 2)")

gene_hsa_symbol <- data.frame(hsa=gene_hsa, symbol=gene_symbol)
gene_hsa_symbol
library(org.Hs.eg.db)

# Gene symbols for the provided gene IDs
gene_ids <- c("79888", "56261", "8760", "9926", "1558", "1576", "8876", "79717", "570", "2876", "2877", "4257", "50484")

# Convert Entrez Gene IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, keytype = "ENTREZID", column = "SYMBOL")

# Display the results
result_df <- data.frame(Entrez_Gene_ID = gene_ids, Gene_Symbol = gene_symbols)
print(result_df)


compound_ids.1 <- c("C00550","C00319","C00195","C00346","C00350","C00157","C04230","C00346","C00819","C00064","C00157","C00245","C00099")
compound_ids.1 

liver.meta.dict %>% str()

get_kegg_name <- function(id){
  name <- keggGet(id)[[1]]$NAME %>% paste(collapse = "---")
  return(name)
}

keggGet("C00350")[[1]]$NAME 

liver.meta.dict %>% colnames()

get_kegg_name("C00350")
Names <- sapply(compound_ids.1, get_kegg_name) 

Names %>% enframe() -> CID.name
CID.name
CID.name %>% mutate(value=strsplit(value, "---")) %>%
  unnest(value) %>% 
  mutate(value = str_replace_all(value, ";", "")) -> CID.name.1
CID.name.1 

liver.meta.dict %>% dplyr::select(Compound) %>% rownames_to_column(var="key") %>% 
  mutate(Compound= strsplit(Compound, " --- ")) %>%
  unnest(Compound) %>% filter(Compound %in% CID.name.1$value)



# Function to query HMDB API for a given compound ID
get_hmdb_accession <- function(compound_id) {
  url <- paste0("http://www.hmdb.ca/metabolites/", compound_id, ".json")
  response <- GET(url)
  GET(url)
  if (http_status(response)$status == 200) {
    data <- fromJSON(content(response, "text"))
    accession <- data$metabolite$accession
    return(accession)
  } else {
    warning("Unable to retrieve data for ", compound_id)
    return(NULL)
  }
}

# Compound IDs
compound_ids <- c("C00550", "C00319", "C00195", "C00346", "C00350", "C00157", "C04230", "C00346", "C00819", "C00064")

# Get HMDB accessions for each compound ID
accessions <- lapply(compound_ids, get_hmdb_accession)

# Combine results into a data frame
result_df <- data.frame(Compound = compound_ids, HMDB_Accession = unlist(accessions))

# View the result
View(result_df)

