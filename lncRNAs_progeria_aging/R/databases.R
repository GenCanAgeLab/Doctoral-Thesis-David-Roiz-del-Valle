# Create a gene metadata table from gtf annotation ----

GTF <- import('/data/genomes/GRCm39_GENCODE/M32/gencode.vM32.primary_assembly.annotation.gtf')

geneMetadata <- as.data.frame(GTF) |>
  filter(type == "gene") |>
  dplyr::select(Chr = seqnames, start, end, width, strand, gene_id, gene_type, gene_name, mgi_id) |>
  mutate(gene_id2 = str_remove(gene_id,"\\.[0-9]+")) |>
  arrange(Chr, start)

write_tsv(geneMetadata, file = "Data/geneMetadata.tsv")

# lncRNA info ----

# Load EVLncRNA database ----

# I will not load this database as it seems useless

# library(readxl)
# 
# EVLnc_info <- read_xls("Data/resources/EVLncRNAs3_alldata/lncRNA_information.xls")
# EV_Lnc_function <- read_xls("Data/resources/EVLncRNAs3_alldata/function_information.xls")
# 
# EVLnc_all <- nest_join(EVLnc_info,
#                        EV_Lnc_function,
#                        by = join_by(ID, `LncRNA name`, Species))

# Load NPInteraction DB ----

NPInter <- read_delim("Data/resources/interaction_NPInterv5.txt.gz")

# Full mouse DB ----

lncRNA_mouseDB <- NPInter |>
  filter(organism == "Mus musculus") |>
  group_by(ncName) |>
  nest(.key = "NPInter_data")


# lncRNA_mouseDB <- full_join(
#   EVLnc_all |>
#     filter(Species == "Mus musculus"),
#   nestNPInter |>
#     filter(organism == "Mus musculus"),
#   by = join_by('LncRNA name' == ncName))

# Save database ----

saveRDS(lncRNA_mouseDB, file = "Data/resources/lncRNA_mouseDB.rds")

#Remove environment
rm(list = setdiff(ls(), lsf.str()))
