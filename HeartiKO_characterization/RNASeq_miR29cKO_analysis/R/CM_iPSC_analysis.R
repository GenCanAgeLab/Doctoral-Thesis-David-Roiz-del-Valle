#File paths and sample table  ------------------

directory <- "/home/droiva/miR29/20220725_RNASeq_CM_iPSCs_miR29ab1/Salmon_quants"

sampleID <- list.files(directory)
Genotypes <- c("WT","WT","WT","WT","KO","KO","KO","KO")

annotation_data <- readRDS("Data/hsap_annotation_data.rds")

filePaths <- file.path(directory, sampleID, "quant.sf")

sampleTable <- tibble(
  files = filePaths,
  names = sampleID,
  Genotype = Genotypes
)

# Import ------------------

# Import with tximeta

se <- tximeta(sampleTable)
gse <- summarizeToGene(se)

# Prefiltering low-count genes

gene_dds <- DESeqDataSet(gse, design = ~ Genotype)
gene_dds$Genotype <- relevel(gene_dds$Genotype, "WT") #Set WT as reference

smallestGroupSize <- 4
keep <- rowSums(counts(gene_dds) >= 10) >= smallestGroupSize
gene_dds <- gene_dds[keep,]

# DESeq analysis

gene_dds <- DESeq(gene_dds)

# Diagnostic plots ------------------

# PCA

gene_vsd <- vst(gene_dds,blind = TRUE)

plotPCA(gene_vsd,intgroup = "Genotype") + geom_label_repel(aes(label=name))

# Clustering

sampleDists <- dist(t(assay(gene_vsd)))


sampleDistMatrix <- as.matrix(sampleDists) # Creamos la matriz de distancias a partir de los datos de distancias
rownames(sampleDistMatrix) <- paste(colnames(gene_vsd), gene_vsd$Genotype, sep="-") # Denominamos a cada fila como muestra-genotipo-dieta
colnames(sampleDistMatrix) <- NULL # Dejamos a las columnas sin nombre
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) # Gama de colores negro


pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Results table ------------------

gene_res <- results(gene_dds, contrast = c("Genotype", "KO", "WT")) |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "gene_id") |>
  left_join(annotation_data,by="gene_id")

# Shrinked FC ------------------

shrinkedRes <- lfcShrink(gene_dds, coef = "Genotype_KO_vs_WT") |>
  as_tibble(rownames = "gene_id") |>
  left_join(annotation_data,by="gene_id")

# Vsd annotation ------------------

rowData(gene_vsd) <- merge(rowData(gene_vsd), annotation_data)

# Saving --------------------------------------------------------------

# Save results for metaanalysis

saveRDS(gene_vsd, file = "Data/vsd_iPSCs.rds")
saveRDS(shrinkedRes, file = "Data/dds_shrinked_iPSCs.rds")
saveRDS(gene_res, file = "Data/Results_iPSCs.rds")


#Remove environment
rm(list = ls())



