# Load results ----

resultsList <- readRDS("Out/LAKI/DESeq_results_sexDiv_LAKI.rds")
tissueList <- names(resultsList)

# Generate Rank lists ---

getRank <- function(res) {
  res |>
    as_tibble(rownames = "gene_id") |>
    filter(gene_type == "protein_coding") |>
    drop_na(pvalue,log2FoldChange,gene_name) |>
    arrange(padj) |> 
    distinct(gene_name, .keep_all = TRUE) |> 
    mutate(rank = -log10(pvalue) * sign(log2FoldChange)) |>
    arrange(desc(rank)) |>
    pull(rank, name = gene_name)
}

rankList <- lapply(resultsList, getRank)

## Import GMTs ----

mouseHallmarks <- gmtPathways("/data/genomes/MsigDB/mh.all.v2024.1.Mm.symbols.gmt")

## Fgsea ----
set.seed(5)

listHGSEA <- lapply(rankList, function(x){
  fgsea(mouseHallmarks, x, nproc = 1)
})
names(listHGSEA) <- tissueList


GSEA_All <- bind_rows(listHGSEA,
                      .id = "Sex_Tissue") |>
  separate_wider_delim(Sex_Tissue, names = c("Sex","Tissue"), delim = " ") |> 
  mutate(pathway = str_replace(pathway,"HALLMARK_",""))

signGSEA <- GSEA_All |> 
  filter(padj < 0.05) |> 
  filter(abs(NES) > 1) |> 
  distinct(pathway) |> 
  left_join(GSEA_All,
            by = join_by(pathway)) |> 
  mutate(colorScale = if_else(padj < 0.05,
                              padj,
                              NA))

## Plot GSEA results ----

ggplot(signGSEA,aes(x = NES,
                    y = fct_reorder(pathway, NES),
                    size = -log10(padj),
                    color = Sex
)
) +
  geom_point() +
  facet_grid(cols = vars(Tissue),
             scales = "free") +
  labs(x = "Normalized enrichment score",
       y="Pathways") +
  geom_vline(xintercept = 0) +
  # scale_color_gradient(low = "lightblue", high = "darkblue", na.value = "gray") +
  theme_bw()

ggsave("Out/LAKI/GSEA.pdf",
       width = 18,
       height = 7)


#Remove environment
rm(list = setdiff(ls(), lsf.str()))
