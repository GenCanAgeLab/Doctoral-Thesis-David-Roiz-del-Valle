#File paths and sample table  ----

directory <- "data/Salmon_quants"

sampleData <- read_delim("data/Sample_info.tsv")

liverSamples <- sampleData |>
  filter(Tissue == "Liver") |>
  mutate(files = file.path(directory, FileName, "quant.sf")) |>
  dplyr::select(files, names = SampName, Genotype)

heartSamples <- sampleData |>
  filter(Tissue == "Heart") |>
  mutate(files = file.path(directory, FileName, "quant.sf")) |>
  dplyr::select(files, names = SampName, Genotype)


# Import  -----
# Import with tximeta

## Liver

liver_se <- tximeta(liverSamples)
liver_gse <- summarizeToGene(liver_se)

## Heart

heart_se <- tximeta(heartSamples)
heart_gse <- summarizeToGene(heart_se)

# Saving  -----

# Save results for metaanalysis

saveRDS(liver_gse, file = "data/R_data/liver_gse.rds")
saveRDS(heart_gse, file = "data/R_data/heart_gse.rds")


#Remove environment
rm(list = ls())
