
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(fgsea)

# Load consensus circRNAs and gene annotation ----

consensus_circDE <- read_tsv("Out/LAKI/circRNA/circRNA_results_consensus.tsv")
geneMetadata <- read_tsv("Data/geneMetadata.tsv")

# Generate Rank lists ---

SignTisGenes <- consensus_circDE |>
  filter(medianPadj < 0.05) |>
  mutate(Tissue_Change = paste(Tissue, Change, sep = "_")) |>
  group_by(Tissue_Change) |>
  summarise(
    genes = list(gene_name),
  ) |>
  ungroup() |> 
  tibble::deframe()


# Load MSigDB databases

mouseHallmarks <- gmtPathways(
  "/data/genomes/MsigDB/mh.all.v2024.1.Mm.symbols.gmt"
)

mouseM2 <- gmtPathways(
  "/data/genomes/MsigDB/m2.all.v2024.1.Mm.symbols.gmt"
)
mouseM5 <- gmtPathways(
  "/data/genomes/MsigDB/m5.all.v2024.1.Mm.symbols.gmt"
)

## Fgsea ----

ORA_H <- lapply(SignTisGenes, function(x){
  fora(mouseHallmarks, x, universe = unique(geneMetadata$gene_name))
})
names(ORA_H) <- names(SignTisGenes)

ORA_M2 <- lapply(SignTisGenes, function(x){
  fora(mouseM2, x, universe = unique(geneMetadata$gene_name))
})
names(ORA_M2) <- names(SignTisGenes)

ORA_GO <- lapply(SignTisGenes, function(x){
  fora(mouseM5, x, universe = unique(geneMetadata$gene_name))
})
names(ORA_GO) <- names(SignTisGenes)


allORA <- lapply(list(ORA_H, ORA_M2, ORA_GO),
                 function(x){
                   bind_rows(x, .id = "Set") |> 
                     rowwise() |> 
                     mutate(overlapGenes = paste(unlist(overlapGenes),
                                                 collapse = ", ")) |> 
                     ungroup()
                 })
names(allORA) <- c("Hallmarks","M2", "GO")

allORA <- bind_rows(allORA, .id = "Database")

write_tsv(allORA, file = "Out/LAKI/circRNA/ORA_genes.tsv")

#Remove environment
rm(list = setdiff(ls(), lsf.str()))
