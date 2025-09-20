source("R/init.R")

# Load gene metadata ----

geneMetadata <- read_tsv("Data/geneMetadata.tsv")

# Load ncRNA info ----

LAKI_ncRNA <- read_tsv(
    "Out/LAKI/ncRNA_LAKI.tsv"
)
TMS_ncRNA <- read_tsv("Out/TMS/ncRNA_TMS.tsv.gz")
TMS_corr <- read_tsv("Out/TMS/Correlation_results_TMS_all.tsv.gz") |>
    select(Tissue, Sex, gene_id, gene_name, gene_type, p_value, estimate) |>
    dplyr::rename(AgePval_TMS = p_value, AgeEstimate_TMS = estimate)

Zhou_ncRNA <- read_tsv("Out/Zhou/ncRNA_Zhou.tsv.gz") |>
    mutate(Sex = "Male")
Zhou_corr <- read_tsv("Out/Zhou/Correlation_results_Zhou_all.tsv.gz") |>
    select(Tissue, Sex, gene_id, gene_name, gene_type, p_value, estimate) |>
    dplyr::rename(AgePval_Zhou = p_value, AgeEstimate_Zhou = estimate)

# Adapt LAKI names

## Change Ileon to small intestine to compare to TMS

LAKI_ncRNA <- LAKI_ncRNA |>
    mutate(Tissue = if_else(Tissue == "Ileon", "Small_intestine", Tissue))

# Change Tissue names in Zhou to make it comparable to TMS

# Merge correlation data ----

TMS_ncRNA <- TMS_ncRNA |>
    full_join(TMS_corr, by = join_by(Sex, Tissue, gene_id, gene_type, gene_name)) |>
    mutate(gene_id = gene_id.y) |>
    dplyr::select(-gene_id.y)


Zhou_ncRNA <- Zhou_ncRNA |>
    full_join(Zhou_corr, by = join_by(Sex, Tissue, gene_id, gene_type, gene_name)) |>
    mutate(gene_id = gene_id.y) |>
    dplyr::select(-gene_id.y)

## Merge PhysioAging ----

PhysioAging <- full_join(
    TMS_ncRNA,
    Zhou_ncRNA,
    by = join_by(
        Sex,
        Tissue,
        gene_id,
        gene_name,
        gene_type
    ),
    suffix = c("", "_Zhou")
)

## Merge PhysioAging with LAKI ----

commonTissues <- Reduce(
    intersect,
    list(
        unique(TMS_ncRNA$Tissue),
        unique(Zhou_ncRNA$Tissue),
        unique(LAKI_ncRNA$Tissue)
    )
)

ncRNA_all <- PhysioAging |>
    full_join(
        LAKI_ncRNA,
        by = join_by(
            gene_id,
            gene_name,
            gene_type,
            Tissue,
            Sex
        ),
        suffix = c("_TMS", "_LAKI")
    )


getGeneLists <- function(tissue, sex) {
    # Progeria
    LAKI_UP <- ncRNA_all |>
        filter(Tissue == tissue) |>
        filter(Sex == sex) |>
        filter(padj_LAKI < 0.05) |>
        filter(log2FoldChange_LAKI > 0) |>
        pull(gene_name)

    LAKI_DOWN <- ncRNA_all |>
        filter(Tissue == tissue) |>
        filter(Sex == sex) |>
        filter(padj_LAKI < 0.05) |>
        filter(log2FoldChange_LAKI < 0) |>
        pull(gene_name)

    # PhysioAging
    TMS_UP <- ncRNA_all |>
        filter(Tissue == tissue) |>
        filter(Sex == sex) |>
        filter(padj_TMS < 0.05 | AgePval_TMS < 0.05) |>
        filter(log2FoldChange_TMS > 0) |>
        pull(gene_name)

    TMS_DOWN <- ncRNA_all |>
        filter(Tissue == tissue) |>
        filter(Sex == sex) |>
        filter(padj_TMS < 0.05 | AgePval_TMS < 0.05) |>
        filter(log2FoldChange_TMS < 0) |>
        pull(gene_name)

    Zhou_UP <- ncRNA_all |>
        filter(Tissue == tissue) |>
        filter(Sex == sex) |>
        filter(padj_Zhou < 0.05 | AgePval_Zhou < 0.05) |>
        filter(log2FoldChange_Zhou > 0) |>
        pull(gene_name)

    Zhou_DOWN <- ncRNA_all |>
        filter(Tissue == tissue) |>
        filter(Sex == sex) |>
        filter(padj_Zhou < 0.05 | AgePval_Zhou < 0.05) |>
        filter(log2FoldChange_Zhou < 0) |>
        pull(gene_name)

    list(
        "LAKI UP" = LAKI_UP,
        "LAKI DOWN" = LAKI_DOWN,
        "TMS UP" = TMS_UP,
        "TMS DOWN" = TMS_DOWN,
        "Zhou UP" = Zhou_UP,
        "Zhou DOWN" = Zhou_DOWN
    )
}

commonTissues <- c("Kidney", "Liver", "Heart", "Small_intestine")

tissueGeneLists <- mapply(
    getGeneLists,
    rep(commonTissues, times = 2),
    rep(c("Male", "Female"), each = length(commonTissues)),
    SIMPLIFY = FALSE
)

names(tissueGeneLists) <- paste(
    rep(commonTissues, times = 2),
    rep(c("Male", "Female"), each = length(commonTissues)),
    sep = "_"
)

# Dataframe with highlighted common LAKI TMS----

Sign_NcRNA_LAKI_TMS <- ncRNA_all |> 
  filter(padj_LAKI < 0.05) |>
  filter(padj_TMS < 0.05 | AgePval_TMS < 0.05) |> 
  filter(abs(log2FoldChange_TMS) > 1 &
           abs(log2FoldChange_LAKI) > 1) |>
  mutate(Change = case_when(
    log2FoldChange_LAKI > 0 & log2FoldChange_TMS > 0 ~ "UP UP",
    log2FoldChange_LAKI > 0 & log2FoldChange_TMS < 0 ~ "UP DN",
    log2FoldChange_LAKI < 0 & log2FoldChange_TMS > 0 ~ "DN UP",
    log2FoldChange_LAKI < 0 & log2FoldChange_TMS < 0 ~ "DN DN"
  )) |> 
  dplyr::select(Sex, Tissue, gene_name, Change, gene_type,
                log2FoldChange_LAKI, log2FoldChange_TMS, AgeEstimate_TMS)

write_tsv(Sign_NcRNA_LAKI_TMS, "Out/sign_LAKI_TMS.tsv")

# Out file for nVenn2 ----

genelists <- unlist(tissueGeneLists, recursive = FALSE)

# Function to pad lists
pad_list <- function(lst) {
    max_len <- max(sapply(lst, length))
    lapply(lst, function(x) {
        length(x) <- max_len
        x
    })
}

padded_gene_lists <- pad_list(genelists)
genelistdf <- as_tibble(padded_gene_lists) # This will create columns with NA for shorter lists

write_tsv(genelistdf, file = "Out/nVenn_TMS_LAKI_Zhou.tsv")

# Overlap analysis

calcOverlap <- function(Tissue, TisName) {
    return(
        tibble(
            Intersect = paste0(TisName, c(
                "_LAKI_TMS_UP",
                "_LAKI_TMS_DOWN",
                "_LAKI_TMS_DIV",
                "_LAKI_Zhou_UP",
                "_LAKI_Zhou_DOWN",
                "_LAKI_Zhou_DIV"
            )),
            Value = c(
                length(intersect(Tissue$"LAKI UP", Tissue$"TMS UP")) / length(Tissue$"TMS UP"),
                length(intersect(Tissue$"LAKI DOWN", Tissue$"TMS DOWN")) / length(Tissue$"TMS DOWN"),
                (length(intersect(Tissue$"LAKI UP", Tissue$"TMS DOWN")) +
                    length(intersect(Tissue$"LAKI DOWN", Tissue$"TMS UP"))) / (length(Tissue$"TMS UP") +
                    length(Tissue$"TMS DOWN")),
                length(intersect(Tissue$"LAKI UP", Tissue$"Zhou UP")) / length(Tissue$"Zhou UP"),
                length(intersect(Tissue$"LAKI DOWN", Tissue$"Zhou DOWN")) / length(Tissue$"Zhou DOWN"),
                (length(intersect(Tissue$"LAKI UP", Tissue$"Zhou DOWN")) +
                    length(intersect(Tissue$"LAKI DOWN", Tissue$"Zhou UP"))) / (length(Tissue$"Zhou UP") +
                    length(Tissue$"Zhou DOWN"))
            ),
            Common = c(
                length(Reduce(intersect, list(Tissue$"LAKI UP", Tissue$"TMS UP", Tissue$"Zhou UP"))) / length(Tissue$"TMS UP"),
                length(Reduce(intersect, list(Tissue$"LAKI DOWN", Tissue$"TMS DOWN", Tissue$"Zhou DOWN"))) / length(Tissue$"TMS DOWN"),
                (length(Reduce(intersect, list(Tissue$"LAKI UP", Tissue$"TMS DOWN", Tissue$"Zhou DOWN"))) +
                    length(Reduce(intersect, list(Tissue$"LAKI DOWN", Tissue$"TMS UP", Tissue$"Zhou UP")))) / (length(Tissue$"TMS UP") +
                    length(Tissue$"TMS DOWN")),
                length(Reduce(intersect, list(Tissue$"LAKI UP", Tissue$"TMS UP", Tissue$"Zhou UP"))) / length(Tissue$"Zhou UP"),
                length(Reduce(intersect, list(Tissue$"LAKI DOWN", Tissue$"TMS DOWN", Tissue$"Zhou DOWN"))) / length(Tissue$"Zhou DOWN"),
                (length(Reduce(intersect, list(Tissue$"LAKI UP", Tissue$"TMS DOWN", Tissue$"Zhou DOWN"))) +
                    length(Reduce(intersect, list(Tissue$"LAKI DOWN", Tissue$"TMS UP", Tissue$"Zhou UP")))) / (length(Tissue$"Zhou UP") +
                    length(Tissue$"Zhou DOWN"))
            )
        )
    )
}

overlaps <- mapply(calcOverlap,
    tissueGeneLists,
    names(tissueGeneLists),
    SIMPLIFY = FALSE
) |>
    bind_rows() |>
    mutate(Value = Value * 100) |>
    mutate(Common = Common * 100) |>
    mutate(Intersect = str_replace(Intersect, "Small_intestine", "Small intestine")) |>
    mutate(Intersect = str_replace(Intersect, "LAKI_TMS_Zhou", "LAKI_TMS Zhou")) |>
    separate_wider_delim(Intersect, delim = "_", names = c("Tissue", "Sex", "Set1", "Set2", "Change")) |>
    mutate(Change = factor(Change, levels = c("UP", "DOWN", "DIV"))) |>
    mutate(Tissue = factor(Tissue, levels = c("Liver", "Kidney", "Heart", "Small intestine"))) |>
    mutate(Intersect = paste(Set1, Set2)) |>
    mutate(Intersect = factor(Intersect, levels = c("LAKI TMS", "LAKI Zhou", "LAKI TMS Zhou"))) |>
    filter(Sex == "Male")


ggplot(overlaps, aes(x = Change, y = Value, fill = Intersect)) +
    geom_col(position = position_dodge(width = 0.8, preserve = "single"), width = 0.8) +
    scale_fill_manual(values = c("#CC79A7", "#0072B2")) +
    geom_col(aes(x = Change, y = Common, group = Intersect),
        position = position_dodge(width = 0.8, preserve = "single"), width = 0.8, inherit.aes = FALSE,
        fill = "#009E73"
    ) +
    facet_wrap(vars(Tissue), scales = "fixed", nrow = 1) +
    theme_classic() +
    labs(y = "Percentage of overlap (%)")

ggsave("Out/Overlap_TMS_LAKI_Zhou.pdf", width = 6, height = 3)
