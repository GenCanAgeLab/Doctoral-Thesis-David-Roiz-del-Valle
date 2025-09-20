# Read RDS ----

res_Liver <- readRDS(file = "data/R_data/liver_results.rds")
res_miR29TG <- readRDS(file = "/home/droiva/miR29/Monahish_RNASeq_miR29TG/Data/results.rds")

# Gene sets ----

UP_DN_Liver <- res_Liver |>
    filter(padj < 0.05) |>
    mutate(Change = factor(log2FoldChange > 0, labels = c("DOWN_KO", "UP_KO"))) |>
    dplyr::select(gene_name, Change) |>
    nest(data = gene_name) |>
    pull(data, name = Change)

UP_DN_Liver <- lapply(UP_DN_Liver, function(x) {
    unlist(as.vector(x))
})

UP_DN_miR29TG <- res_miR29TG |>
    filter(padj < 0.05) |>
    mutate(Change = factor(log2FoldChange > 0, labels = c("DOWN_TG", "UP_TG"))) |>
    dplyr::select(gene_name, Change) |>
    nest(data = gene_name) |>
    pull(data, name = Change)

UP_DN_miR29TG <- lapply(UP_DN_miR29TG, function(x) {
    unlist(as.vector(x))
})

geneSets <- c(UP_DN_Liver, UP_DN_miR29TG)


plotVenn(geneSets, outFile = "Out/KO_TG.svg")
