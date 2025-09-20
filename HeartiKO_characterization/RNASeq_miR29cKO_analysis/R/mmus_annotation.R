
# Load GENCODE data ---- 

RefSeq_gencode <- read_delim("/data/genomes/GRCm39_GENCODE/M32/gencode.vM32.metadata.RefSeq",col_names = c("transcript_id","RefSeq_gene","RefSeq_protein"))
GTF <- rtracklayer::import("/data/genomes/GRCm39_GENCODE/M32/gencode.vM32.primary_assembly.annotation.gtf")

# Import predictions ----

TargetScan_pred <- read_delim("/data/miR29/resources/Predicted_Targets_Info.default_predictions.txt", col_names = T)
miRDB_pred <- read_delim("/data/miR29/resources/miRDB_v6.0_prediction_result.txt", col_names = F)
microT_interaction <- read_delim("/data/miR29/resources/interactions_mouse.microT.mirbase.txt", col_names = T)

# Process predictions ----

## TargetScan ----

TargetScan_seeds <- TargetScan_pred |>
  filter(`Species ID` == "10090" & `miR Family` == "miR-29-3p") |> # Select only Mmusculus predictions and miR-29 family targets
  dplyr::count(`Gene Symbol`) |> # Now we count the number of seed sequences of each transcript
  dplyr::rename(TargetScan_nseed = n)

## miRDB ----

RefSeq_gencode <- RefSeq_gencode |>
  separate_wider_delim(RefSeq_gene,delim = ".",too_many = "drop",names = "RefSeq_gene")

miRDB_seeds <- miRDB_pred |>
  filter(X1 == "mmu-miR-29a-3p" | X1 == "mmu-miR-29b-3p" | X1 == "mmu-miR-29c-3p") |>
  left_join(RefSeq_gencode, by = join_by(X2 == RefSeq_gene), relationship = "many-to-many") |>
  dplyr::rename(RefSeq_gene = X2, miRNA = X1) |>
  arrange(desc(X3)) |>
  left_join(as.data.frame(GTF) |> 
              filter(type == "transcript")) |>
  dplyr::select(miRNA, gene_name, X3) |>
  arrange(desc(X3)) |>
  distinct(miRNA, gene_name, .keep_all = T) |>
  drop_na() |>
  dplyr::rename(miRDB_score = X3) |>
  group_by(gene_name,miRDB_score) |>
  reframe(miRDB_miRNA = paste(miRNA, collapse = ", ")) |>
  ungroup() |>
  arrange(desc(miRDB_score)) |>
  distinct(gene_name, .keep_all = T)

## microT ----

microT_seeds <- microT_interaction |>
  filter(mirna %in% c("mmu-miR-29a-3p","mmu-miR-29b-3p","mmu-miR-29c-3p")) |>
  filter(interaction_score > 0.7) |> #CutOff for valid targets
  arrange(desc(interaction_score)) |>
  group_by(ensembl_gene_id) |>
  summarise(microT_miRNA = paste(mirna,collapse = ","),
            microT_score = max(interaction_score)) |>
  ungroup() |>
  distinct(ensembl_gene_id, .keep_all = T) |>
  dplyr::rename(gene_id = ensembl_gene_id)

# Annotation 

annotation_data <- as.data.frame(GTF) |>
  filter(type=="gene") |>
  dplyr::select(gene_id, gene_name, gene_type, tag) |>
  mutate(gene_id_noVersion = str_remove(gene_id,"\\.\\d+$")) |>
  left_join(TargetScan_seeds, by = join_by(gene_name == 'Gene Symbol')) |>
  left_join(miRDB_seeds) |>
  left_join(microT_seeds, by = join_by(gene_id_noVersion == gene_id)) |>
  dplyr::select(-gene_id_noVersion) |>
  mutate(isTarget = (!is.na(miRDB_score) | !is.na(TargetScan_nseed) | !is.na(microT_miRNA)))

# Save annotation ----

saveRDS(annotation_data, file = "Data/mmus_annotation_data.rds")

#Remove environment
rm(list = ls())


