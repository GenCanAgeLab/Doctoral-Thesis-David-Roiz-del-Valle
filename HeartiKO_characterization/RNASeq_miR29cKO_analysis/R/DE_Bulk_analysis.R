
#File paths and sample table  ------------------

directory <- "/home/droiva/miR29/20220427_RNASeq_HeartBulk_Myh6Cre_miR29fl/salmon"
sampleID <- c("247","248","249","254","255","256","257","259")
Genotypes <- c("Heart-iKO","miR-29fl/fl","miR-29fl/fl","miR-29fl/fl",
               "Heart-iKO","Heart-iKO","Heart-iKO","miR-29fl/fl")

annotation_data <- readRDS("Data/mmus_annotation_data.rds")

filePaths <- file.path(directory, paste0("DavidRV_", sampleID), "quant.sf")

sampleTable <- tibble(
  files = filePaths,
  names = sampleID,
  Genotype = Genotypes
)

sampleTable <- sampleTable |>
  filter(names != "256") #Remove outlier

# Import  ------------------
# Import with tximeta

se <- tximeta(sampleTable)
gse <- summarizeToGene(se)

# Prefiltering low-count genes

gene_dds <- DESeqDataSet(gse, design = ~ Genotype)
gene_dds$Genotype <- relevel(gene_dds$Genotype, "miR-29fl/fl") #Set reference

smallestGroupSize <- 4
keep <- rowSums(counts(gene_dds) >= 10) >= smallestGroupSize
gene_dds <- gene_dds[keep,]

# DESeq analysis

gene_dds <- DESeq(gene_dds)

# Diagnostic plots  ------------------

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

# Results table  ------------------

gene_res <- results(gene_dds,
                    contrast = c("Genotype", "Heart-iKO", "miR-29fl/fl"),
                    alpha = 0.05) |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "gene_id") |>
  left_join(annotation_data,by="gene_id")

# Shrinked FC  ------------------

shrinkedRes <- lfcShrink(gene_dds, coef = "Genotype_Heart.iKO_vs_miR.29fl.fl") |>
  as_tibble(rownames = "gene_id") |>
  left_join(annotation_data,by="gene_id")

# Vsd annotation ------------------

rowData(gene_vsd) <- merge(rowData(gene_vsd), annotation_data)

# Saving  ------------------

# Save results for metaanalysis

saveRDS(gene_vsd, file = "Data/vsd_Bulk.rds")
saveRDS(shrinkedRes, file = "Data/dds_shrinked_Bulk.rds")
saveRDS(gene_res, file = "Data/Results_Bulk.rds")

#Remove environment
rm(list = ls())

