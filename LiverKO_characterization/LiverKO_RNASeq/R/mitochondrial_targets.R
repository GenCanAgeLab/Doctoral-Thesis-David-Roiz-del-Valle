
source("R/init.R")

# Load results and GSEA results from Hallmarks ----

gene_res <- readRDS("Data/results.rds")

GSEA_H <- read_tsv("Out/GSEA/GSEA_H.tsv")
GSEA_GO <- read_tsv("Out/GSEA/GSEA_M5_GO.tsv")

# Filter res for OxPhos genes ----

OxPhos_leading <- strsplit(GSEA_H |> filter(pathway == "OXIDATIVE_PHOSPHORYLATION") |> pull(leadingEdge),
                           ", ") |> 
  unlist()

ETC_leading <- strsplit(GSEA_GO |> filter(pathway == "GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY") |> pull(leadingEdge),
                           ", ") |> 
  unlist()



OxPhos <- gene_res |> 
  filter(gene_name %in% OxPhos_leading)

ETC <- gene_res |> 
  filter(gene_name %in% ETC_leading)

gene_res |>
  filter(padj < 0.01 & log2FoldChange > 0.5) |>
  filter(isTarget) |>
  pull(gene_name)

View(gene_res |>
       filter(padj < 0.01 & log2FoldChange > 0.5) |>
       filter(isTarget))
