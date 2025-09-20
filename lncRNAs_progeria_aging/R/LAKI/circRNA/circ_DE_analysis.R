library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

# Read files

circSexTisDds <- readRDS("Out/LAKI/circRNA/dds_circRNA.rds")
resultsList <- readRDS("Out/LAKI/circRNA/results_circRNA.rds")

resAll <- read_tsv("Out/LAKI/circRNA/circRNA_results_all.tsv")

# All circ annotation ----

circ_annotation <- read_tsv("/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/Results/circRNA_ID_across_tools.tsv")
circ_gene_annotation <- read_tsv("/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/Results/All_circRNAs_annotated_cross_with_circDBs.tsv")

geneMetadata <- read_tsv("Data/geneMetadata.tsv")

circ_annotation <- circ_annotation |>
  pivot_longer(
    cols = c(CIRI2, Terrace, circExplorer2, find_circ, starchip),
    names_to = "Tool",
    values_to = "circID_tool"
  ) |>
  left_join(circ_gene_annotation,
    by = join_by(circID, seqnames, start, end, width)
  ) |>
  left_join(
    geneMetadata |>
      dplyr::select(gene_id, gene_type, gene_name),
    by = join_by(gene_id)
  )

# Merge annotation with results ----

annotated_res <- resAll |>
  left_join(circ_annotation,
    by = join_by(Tool, circRNA_ID == circID_tool)
  )

# DE analysis per tool per tissue/sex ----

circ_DE_tisSex <- annotated_res |>
  filter(padj < 0.05) |>
  group_by(Tissue, Sex, Tool) |>
  summarize(
    UP = sum(log2FoldChange > 0, na.rm = TRUE),
    DOWN = sum(log2FoldChange < 0, na.rm = TRUE)
  ) |>
  pivot_longer(
    cols = c(UP, DOWN),
    names_to = "Change",
    values_to = "DE_genes"
  ) |>
  mutate(Change = factor(Change, levels = c("UP", "DOWN")))

ggplot(circ_DE_tisSex, aes(x = Tissue, y = DE_genes, fill = Tool)) +
  geom_col(position = position_dodge(width = 0.8, preserve = "single"), width = 0.8) +
  facet_grid(
    cols = vars(Change),
    rows = vars(Sex)
  ) +
  labs(y = "Differentially expressed circRNAs") +
  theme_classic()

ggsave("Out/LAKI/circRNA/DE_tissue_sex.pdf",
  width = 8,
  height = 5
)

# Consensus circRNAs ----

consensus_circDE <- annotated_res |>
  drop_na(padj) |>
  mutate(Change = factor(sign(log2FoldChange),
    levels = c(1, -1),
    labels = c("UP", "DOWN")
  )) |>
  filter(padj < 0.05) |>
  group_by(
    Tissue, Sex,
    circID, Cons_coords, type, gene_id, gene_type, gene_name,
    Change, circAtlas, circBase, circPedia
  ) |>
  summarize(
    NtoolsDE = n(),
    medianLFC = median(log2FoldChange),
    medianPadj = median(padj),
    medianStat = median(stat)
  ) |>
  filter(NtoolsDE > 2) # DE detected by more than 2 tools

consensus_circDE_tisSex <- consensus_circDE |>
  filter(NtoolsDE > 2) |> # DE detected by more than 2 tools
  group_by(Tissue, Sex, Change) |>
  summarize(
    DE_genes = n()
  )

ggplot(consensus_circDE_tisSex, aes(x = Tissue, y = DE_genes, fill = Change)) +
  geom_col(position = position_dodge(width = 0.8, preserve = "single"), width = 0.8) +
  facet_wrap(vars(Sex)) +
  scale_fill_manual(values = c("#ff3b30","#4cd964")) +
  theme_classic() +
  labs(y = "Differentially expressed circRNAs") +
  scale_y_continuous(label = scales::number_format(accuracy = 1))

ggsave("Out/LAKI/circRNA/DE_tissue_sex_consensus.pdf",
  width = 6,
  height = 3
)

consensus_circDE_tisSextype <- consensus_circDE |>
  filter(NtoolsDE > 2) |> # DE detected by more than 2 tools
  group_by(Tissue, Sex, Change, type) |>
  summarize(
    DE_genes = n()
  )

ggplot(consensus_circDE_tisSextype, aes(x = Tissue, y = DE_genes, fill = type)) +
  geom_col(position = position_stack(), width = 0.8) +
  facet_grid(
    cols = vars(Change),
    rows = vars(Sex)
  ) +
  labs(y = "Differentially expressed circRNAs") +
  theme_classic()

ggsave("Out/LAKI/circRNA/DE_tissue_sex_consensus_type.pdf",
  width = 8,
  height = 4
)

# Consistently-changed circRNAs ----

consistent_circDE <- consensus_circDE |>
  filter(NtoolsDE > 2) |> # DE detected by more than 2 tools
  group_by(
    Sex, circID, Cons_coords, type, gene_id, gene_type, gene_name,
    Change, circAtlas, circBase, circPedia
  ) |>
  summarize(
    NtissuesDE = n_distinct(Tissue),
    Tissues = paste(unique(Tissue),
      collapse = ", "
    ),
    medianLFC = median(medianLFC),
    medianPadj = median(medianPadj)
  )

# Detected circRNAs ----

detectedCircs <- circ_annotation |> 
  dplyr::select(circID,Tool,circID_tool,circAtlas,circBank,circBase,circPedia,type) |> 
  filter(!is.na(circID_tool)) |> 
  mutate(isAnnotated = !is.na(circAtlas) | !is.na(circBank) | !is.na(circBase) | !is.na(circPedia))

dbCirc_count <- detectedCircs |> 
  group_by(Tool, isAnnotated) |> 
  summarize(
    n_Circ = n_distinct(circID_tool, na.rm = TRUE)
  )

typeCirc_count <- detectedCircs |> 
  group_by(Tool, type) |> 
  summarize(
    n_Circ = n_distinct(circID_tool, na.rm = TRUE)
  ) |> 
  mutate(type = if_else(n_Circ < 50, "Other",type)) |> 
  group_by(Tool, type) |> 
  summarize(
    n_Circ = sum(n_Circ, na.rm = TRUE)
  )


ggplot(dbCirc_count, aes(x = Tool, y = n_Circ, fill = isAnnotated)) +
  geom_col() +
  theme_classic() +
  labs(y = "Number of detected circRNAs", fill = "Previously\nannotated")

ggsave("Out/LAKI/circRNA/annotatedCirc_db.pdf",
       width = 6,
       height = 3
)

ggplot(typeCirc_count, aes(x = Tool, y = n_Circ, fill = type)) +
  geom_col() +
  theme_classic() +
  labs(y = "Number of detected circRNAs", fill = "circRNA\ntype")

ggsave("Out/LAKI/circRNA/annotatedCirc_type.pdf",
       width = 6,
       height = 3
)

# Save data ----

write_tsv(annotated_res,
  file = "Out/LAKI/circRNA/circRNA_results_annotated.tsv"
)

write_tsv(consensus_circDE,
  file = "Out/LAKI/circRNA/circRNA_results_consensus.tsv"
)

write_tsv(consistent_circDE,
  file = "Out/LAKI/circRNA/circRNA_results_consistent.tsv"
)
