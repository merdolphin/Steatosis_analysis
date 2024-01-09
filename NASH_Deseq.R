library(tidyverse)
library(ggplot2)
library(tximport)
library(GenomicFeatures)
library(DESeq2)
library(caret)
library(party)
library(randomForest)
library(mice)
library(org.Hs.eg.db)
library(ggpubr)
library(biomaRt)
library(tximport)
library(VIM)
library(readxl)
library(clusterProfiler)
library(BiocParallel)

setwd("/home/lina/scratch/18.NASH/Rscripts/")


sampled_ID <- read.csv("../1.sample/sampleName_clientId.txt", sep = "\t") 

sampleInfo <- read_xlsx("/home/lina/scratch/18.NASH/1.sample/HumanPatients_230216.xlsx",
                          sheet = "Sheet1")


colnames.1 <- read.table("/home/lina/scratch/18.NASH/Rscripts/colnames.3.csv", sep=",")
colnames(sampleInfo) <- colnames.1
 
### Ex-Diabetes as 1
sampleInfo$Diabetes[sampleInfo$Diabetes == "2"] <- "1"
sampleInfo$T2DM[sampleInfo$T2DM == "2"] <- "1"


sampleInfo[sampleInfo$ID %in% sampled_ID$clientId,] -> samples94
samples94 %>% 
  mutate_if(is.character, as.factor) %>%
  mutate_if(is.integer, as.factor) %>%
  mutate(Histo_STEATOSIS=str_replace(Histo_STEATOSIS, "3","2")) %>%
  mutate(Histo_STEATOSIS = as.factor(Histo_STEATOSIS)) -> sampled

sampled_ID[sampled_ID$clientId %in% sampled$ID,] -> samples_filtered

aggr(sampled, plot=F)
samples_filtered 
samples_filtered$clientId == sampled$ID


left_join(samples_filtered, sampled, by=c("clientId" = "ID")) -> sampled.1

sampled.1  %>% dim()

files_import_cds <- paste0("2.cds_count_data/",
                           sampled.1$X.sampleName,"_good/quant.sf")
files_import_cds
names(files_import_cds) <- samples_filtered$clientId
file.exists(files_import_cds)

files_import_cdna <- paste0("/home/lina/scratch/16.R.singleCellRNAseq/1.Kaldis.Sequence/2.cdna_count_data/",sampled.1$X.sampleName,"_good/quant.sf")
files_import_cdna
names(files_import_cdna) <- samples_filtered$clientId
file.exists(files_import_cdna)


txdb <- makeTxDbFromGFF("/home/lina/scratch/16.R.singleCellRNAseq/0.ref_homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz")
columns(txdb)
keytypes(txdb)
keys(txdb)
k <- keys(txdb, keytype = "TXNAME")
k
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene


columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

gene2sym <- select(org.Hs.eg.db,
                   keys = tx2gene[,2],
                   column = "SYMBOL",
                   keytype = "ENSEMBL")

gene2sym[1:5,]


mat_gse <- tximport(files_import_cdna,
                    type = "salmon",
                    tx2gene = tx2gene,
                    ignoreAfterBar = TRUE,
                    ignoreTxVersion = TRUE)

save(mat_gse, file="mat_gse_cdna.rda")
load(file="/home/lina/scratch/16.R.singleCellRNAseq/1.Kaldis.Sequence/mat_gse_cdna.rda")



################################
############################### from me
mat_gse$counts %>% colnames()
dds <- DESeqDataSetFromTximport(mat_gse,
                                sampled.1
                              # ~ Histo_STEATOSIS)
                                ~ Gender + Weight + Histo_STEATOSIS) ## the condition at the end

colnames(dds) == sampled.1$clientId

dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(64))



res <- results(dds)
res
resultsNames(dds)
#res <- results(dds, contrast = c("Histo_STEATOSIS","3","1"))
res <- results(dds, contrast = c("Histo_STEATOSIS","2","1"))
res <- results(dds, contrast = c("Histo_STEATOSIS","1","0"))

convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}
res$mm_symbol <- convertIDs(row.names(res), "ENSEMBL","SYMBOL", org.Hs.eg.db)
res$mm_entrezgene <- convertIDs(row.names(res), "ENSEMBL","ENTREZID", org.Hs.eg.db)

res05 <- res[ which(res$padj < 0.05), ]

res05 <- res[ which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ] 

res05 <- res[ which(res$padj < 0.05 & res$log2FoldChange < - 1), ] 


res05 %>% as.data.frame() %>% dplyr::select(mm_symbol) 
res05 -> res05FC2.S1vsS0
res05FC2.S1vsS0 %>% na.omit() 
save(res05FC2.S1vsS0, file="res05FC2.S1vsS0.rda")
load(file="/home/lina/scratch/16.R.singleCellRNAseq/1.Kaldis.Sequence/res05FC2.S1vsS0.rda")
res05FC2.S1vsS0

plotCounts(dds, gene="ENSG00000085999", intgroup = "Histo_STEATOSIS", returnData = T) %>%
  group_by(Histo_STEATOSIS) %>% summarize(mean(count))

### only based on padj <- 0.05
res05 %>% rownames()-> res05_2_0
res05 %>% rownames()-> res05_2_1
res05 %>% rownames()-> res05_1_0

res05
res05$mm_symbol -> res05_2_0.516
res05$mm_symbol-> res05_2_1.2
res05$mm_symbol-> res05_1_0.581

c(res05_2_0, res05_1_0) %>% unique() %>% length()

c(res05_2_0, res05_2_1, res05_1_0) %>% unique() %>% as.data.frame() %>% setNames("listRes05")-> res05_adjusted

res05_adjusted %>% head()


intersect(up10, down10)

save(up10,file="up10.rda")
save(down10,file="dow10.rda")


resSig %>% as.data.frame() %>% filter(mm_symbol %in% variants)

variants <- c("PNPLA3", "TM6SF2","MBOAT7", "GCKR","HSD17B13",
              "SOD2","NCAN", "APOC3", "FGD5","CITED2")

res %>% as.data.frame() 

res[ which(res$mm_symbol != "HLA-DMA"),] -> res

pdf("try.pdf")
ggmaplot(res, main = expression("Steatosis (S1+S2+S3) vs. no steatosis (S0)"),
         fdr = 0.01, fc = 2, size = 2, label.rectangle = TRUE,
         palette = c("red", "blue", "darkgray"),
         genenames = as.vector(res$mm_symbol),
         legend = "top", top = 25,
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_classic())

dev.off()

help(ggmaplot)


rld <- rlog(dds)
plotPCA(rld, intgroup="Histo_STEATOSIS", ntop=500)

help("plotPCA")

# also possible to perform custom transformation:

# shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                           colData=colData(dds))
# the call to DESeqTransform() is needed to
# trigger our plotPCA method.
plotPCA( DESeqTransform( se ) )

vsd <- vst(dds.diabetes, blind=FALSE)

plotPCA(vsd, intgroup=c("Diabetes", "Liver.Gross.grouping"))

#############

convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

sampled$Liver.Gross.grouping

sampled %>% str()

resultsNames(dds)

resLFC <- lfcShrink(dds, coef="Liver.Gross.grouping_2_vs_0", type = "apeglm")

res$mm_symbol <- convertIDs(row.names(res), "ENSEMBL","SYMBOL", org.Hs.eg.db)
res$mm_entrezgene <- convertIDs(row.names(res), "ENSEMBL","ENTREZID", org.Hs.eg.db)

summary(res)

res05 <- results(dds, alpha=0.05)

summary(res05)

head( res[order(res$padj),], n=10 )

resSig <- res[ which(res$padj < 0.05 & abs(res$log2FoldChange) > 2), ] 

resSig_FC2 <- resSig <- res[ which(res$padj < 0.05 & abs(res$log2FoldChange) > 2), ] 

resSig_FC2
summary(resSig)


resSig  %>% data.frame() 
resSig
head(res)
res %>% data.frame() %>% dplyr::filter(mm_symbol=="CTSB")

mat_gse$abundance

rownames(res)

library("biomaRt")

ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

listAttributes(mart = ensembl) %>% filter(str_detect(name,'transcript'))



annot <- getBM(attributes = c('ensembl_gene_id',
                              'hgnc_symbol',
                              'uniprot_gn_symbol',
                              'ensembl_transcript_id'),
               filters = 'ensembl_gene_id',
               values = rownames(res),
               mart = ensembl)


res %>% data.frame() %>% rownames_to_column("ensembl_gene_id") %>% left_join(annot) -> res.genes

res.genes
res.genes %>% filter(hgnc_symbol == "KRT18")

res.genes %>%  filter(str_detect(hgnc_symbol,"KRT")) %>% arrange(desc(padj))



res.genes %>%  filter(str_detect(hgnc_symbol,"KRT")) %>% arrange(desc(baseMean)) %>% 
  filter(baseMean >0) %>% dplyr::select(mm_symbol,baseMean,log2FoldChange, pvalue) -> res.pie

res.pie %>% mutate(p=case_when(pvalue < 0.05 ~ "p < 0.05",
                               TRUE ~ "p > 0.05")) -> res.pie

res.pie

pdf("pieDonut.pdf")

PieDonut(res.pie, aes(p,mm_symbol, count=baseMean),
         explode=1, explodePos = 0.05,
         r0=0,
         r1=0.4,
         r2=1,
         showDonutName = TRUE,
         labelpositionThreshold = 0.01,
         donutLabelSize = 6
)

dev.off()



pdf("NASH_vs_noPathObeses.pdf")

ylim <- c(-15, 8)
ggmaplot(res, main = expression("NASH vs. noPath"),
         fdr = 0.05, fc = 4, size = 2,
         palette = c("red", "blue", "darkgray"),
         genenames = as.vector(res$mm_symbol),
         legend = "top", top = 20,
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_classic())

dev.off()


dds
sampled$Liver.Fibrosis.grouping

# variance stabilizing Transformation
vsd <- vst(dds, blind = FALSE, nsub=1000)

vsd


plotPCA(vsd, intgroup = c("Liver.Gross.grouping","Liver.Fibrosis.grouping"))

plotPCA(vsd, intgroup = "Liver.Fibrosis.grouping")

pdf("plotPCA_1.pdf")

plotPCA(vsd, intgroup = "Liver.Gross.grouping")

dev.off()


plotPCA_n(vsd, intgroup = "Liver.Gross.grouping", returnData = TRUE) -> data1

data1


par(mfrow=c(2,2))
a =0
for (i in colnames(sampled[60])){
  print(i)
  g <- i
  plotPCA_n(vsd, intgroup = g) 
}


########## cDNA analysis ##########


mat_gse <- tximport(files_import_cds,
                    type = "salmon",
                    tx2gene = tx2gene,
                    ignoreAfterBar = TRUE,
                    ignoreTxVersion = TRUE)


colnames(sampled)
colData(sampled)

## put NAFL in the same group as NASH
sampled$Liver.Gross.grouping
##sampled %>% mutate(Liver.Gross.grouping = str_replace_all(Liver.Gross.grouping,"1","0")) -> sampled

sampled$Liver.Gross.grouping

vignette("DESeq2")


sampled[, (colSums(is.na(sampled)) < 9)]

sampled$Liver.Fibrosis.grouping

help(mutate_if)

sampled %>% mutate(Liver.Fibrosis.grouping = str_replace(Liver.Fibrosis.grouping,"4","2")) -> sampled_4to0

colnames(sampled)
dds <- DESeqDataSetFromTximport(mat_gse,
                                sampled,
                                ~Liver.Gross.grouping )


dds <- DESeq(dds)
dds$Liver.Gross.grouping

res <- results(dds)

res

convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

sampled$Liver.Gross.grouping

resultsNames(dds)

resLFC <- lfcShrink(dds, coef="Liver.Gross.grouping_2_vs_0", type = "apeglm")

res$mm_symbol <- convertIDs(row.names(res), "ENSEMBL","SYMBOL", org.Hs.eg.db)
res$mm_entrezgene <- convertIDs(row.names(res), "ENSEMBL","ENTREZID", org.Hs.eg.db)


res

dds
sampled$Liver.Fibrosis.grouping

# variance stabilizing Transformation
vsd <- vst(dds, blind = FALSE, nsub=1000)

vsd


plotPCA(vsd, intgroup = c("Liver.Gross.grouping","Liver.Fibrosis.grouping"))

plotPCA(vsd, intgroup = "Liver.Fibrosis.grouping")

plotPCA(vsd, intgroup = "Liver.Gross.grouping")

plotPCA_n(vsd, intgroup = "Liver.Gross.grouping", returnData = TRUE) -> data1

data1


par(mfrow=c(2,2))
a =0
for (i in colnames(sampled[60])){
  print(i)
  g <- i
  plotPCA_n(vsd, intgroup = g) 
}






plotPCA(vsd, intgroup = "Total.Bili..3.15.")

plotPCA_n(vsd, intgroup = "Total.Bili..3.15.")


?plotPCA
getMethod("plotPCA","DESeqTransform")
data1

plotPCA_n = function (object, ...) 
{
  .local <- function (object, intgroup = "condition", ntop = 500, 
                      returnData = FALSE) 
  {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    print(round(percentVar[1] * 100)) 
    print(intgroup)
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }
    else {
      colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                    intgroup.df, name = colnames(object))
    if (returnData) {
      attr(d, "percentVar") <- percentVar[1:2]
      return(d)
    }
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
      geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 
                                                          100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
                                                                                                              100), "% variance")) + coord_fixed()
  }
  .local(object, ...)
}




#How to get PCA scree plot?

## calculate the variance for each gene
rv <- rowVars(assay(vsd))

## select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]

## perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsd)[select,]))

## the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

##plot the "percentVar"
scree_plot=data.frame(percentVar)
scree_plot[,2]<- c(1:94)

colnames(scree_plot)<-c("variance","component_number")
ggplot(scree_plot, mapping=aes(x=component_number, y=variance))+geom_bar(stat="identity")


summary(res)

plotMA(res)

resSig <- res[ which(res$padj < 0.05 & abs(res$log2FoldChange) > 2), ]
resSig$mm_symbol





vsd <- vst(dds, blind=FALSE)
corrs <- cor(assay(vsd), method="spearman")
corr.dists <- as.dist(1 - corrs)

corrs
library("pheatmap")
colors <- colorRampPalette(c("red","white","blue"))(99)
pheatmap(corrs,
         clustering_distance_rows=corr.dists,
         clustering_distance_cols=corr.dists,
         col=colors)

sampleDists <- dist(t(assay(vsd)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)


BiocManager::install("M3Drop")

library(edgeR)

files_import <- paste0("2.cds_count_data/",samples_filtered$X.sampleName,"_good/")
mat_edgeR <- catchSalmon(files_import)


d <- DGEList(counts = mat_edgeR$counts/mat_edgeR$annotation$Overdispersion,
             samples = sampled)

d.full <- d

apply(d$counts, 2, sum)
keep <- rowSums(cpm(d) >100) >=2
d <- d[keep,]
dim(d)

d$samples$Liver.Fibrosis.grouping



d <- calcNormFactors(d)
colnames(d) <- sampled$ID

plotMDS(d, method="bcv", col=as.numeric(d$samples$Liver.Fibrosis.grouping)) 




#plotMDS(d, method="bcv", col=as.numeric(d$samples$i), plot = FALSE) -> d0
#print(max(d0$var.explained))  


plotMDS(d, method="bcv", col=as.numeric(d$samples$Liver.Fibrosis.grouping))
max(d0$var.explained)

legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)

help(plotMDS)
