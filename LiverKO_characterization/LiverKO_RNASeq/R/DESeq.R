
source("R/init.R")

annotation_data <- readRDS("Data/annotation_data.rds")

#File paths and sample table  ------------------

samples <- list.files("Salmon_quants")
quantData <- file.path("Salmon_quants", samples, "quant.sf")

file.exists(quantData)

# Design

design <- read_tsv("Design.tsv")


sampleTable <- design |> 
  mutate(
    files = quantData,
    names = samples
  ) |> 
  mutate(Genotype = if_else(Genotype == "CRE/+",
                            "LiverKO",
                            "miR29fl")) |> 
  dplyr::select(files, names, Genotype)


# Import  ------------------
# Import with tximeta

se <- tximeta(sampleTable)
gse <- summarizeToGene(se)

# Prefiltering low-count genes

gene_dds <- DESeqDataSet(gse, design = ~ Genotype)
gene_dds$Genotype <- relevel(gene_dds$Genotype, "miR29fl") #Set reference

# smallestGroupSize <- 4
# keep <- rowSums(counts(gene_dds) >= 10) >= smallestGroupSize
# gene_dds <- gene_dds[keep,]

# DESeq analysis

gene_dds <- DESeq(gene_dds)

# Diagnostic plots  ------------------

# PCA

gene_vsd <- vst(gene_dds,blind = TRUE)

plotPCA(gene_vsd,intgroup = "Genotype") +
  geom_label_repel(aes(label=name)) +
  theme_classic()

ggsave("Out/PCA.pdf")

# Clustering

sampleDists <- dist(t(assay(gene_vsd)))


sampleDistMatrix <- as.matrix(sampleDists) # Creamos la matriz de distancias a partir de los datos de distancias
rownames(sampleDistMatrix) <- paste(colnames(gene_vsd), gene_vsd$Genotype, sep="-") # Denominamos a cada fila como muestra-genotipo-dieta
colnames(sampleDistMatrix) <- NULL # Dejamos a las columnas sin nombre
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) # Gama de colores negro

dev.off()
pdf("Out/heatmap_samples.pdf")

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

dev.off()

# Exclude outlier

# gene_dds <- gene_dds[,colnames(gene_dds) != "A991"]
# gene_dds <- DESeq(gene_dds)

# Results table  ------------------

gene_res <- results(gene_dds,
                    contrast = c("Genotype", "LiverKO", "miR29fl"),
                    alpha = 0.05) |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "gene_id") |>
  left_join(annotation_data,by="gene_id")

# Shrinked FC  ------------------

shrinkedRes <- lfcShrink(gene_dds, coef = "Genotype_LiverKO_vs_miR29fl") |>
  as_tibble(rownames = "gene_id") |>
  left_join(annotation_data,by="gene_id")

# Vsd annotation ------------------

rowData(gene_vsd) <- merge(rowData(gene_vsd), annotation_data)

# Saving  ------------------

# Save results for metaanalysis

saveRDS(gene_vsd, file = "Data/vsd.rds")
saveRDS(shrinkedRes, file = "Data/res_shrinked.rds")
saveRDS(gene_res, file = "Data/results.rds")

# Save results table

write_tsv(gene_res, file = "Out/results.tsv")

#Remove environment
# rm(list = ls())

