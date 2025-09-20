
# Load gene metadata ----

geneMetadata <- read_tsv("Data/geneMetadata.tsv")

# Load and arrange sample info ----

## Some samples are from a different experiment, so they will be removed.

# Load Sample info
sampdf <- read.table("Data/Info/Sample_info_complete_genotiped.tsv",
                     header=T,sep='\t') |> # Cannot use readr as tibbles do not support rownames
  as_tibble() |> # We convert it to tibble as tximeta has problems with rownames
  
  # Add a column with file name
  mutate(
    files = paste0("/home/jandrulas/20231017_RNASeq_LAKI/STAR_aligns/", SampleID,"_ReadsPerGene.out.tab")) |>
  
  # Remove some animals from another experiment  
  filter(Notes != 'Ratones de Face-NR') |>
  
  # Fix wrong genotypes and remove HT animals (see Bash/Regenotyping) 
  dplyr::rename(Regenotype = Infered_Genotype)  |>
  filter(Regenotype != "HT")


# Create names column
sampdf <- sampdf |>
  dplyr::rename(names = SampleID) |>
  dplyr::select(names, files,
                MouseNo, Genotype = Regenotype, Gender, Tissue, DAB)

# Import

dds <- DESeqDataSetFromHTSeqCount(sampdf,
                                  directory = "/",
                                  design = ~ Genotype + Gender)


# We need to remove the first rows from the HTSeq-count files
keepHTSeq <- rownames(dds)[5:length(rownames(dds))]
dds <- dds[keepHTSeq,]


# DESeq2 ----

analysisDESeq <- function(dds,Tissue) {
  
  tmpDds <- dds[,colData(dds)[,"Tissue"] == Tissue] #First, select the samples from our target tissue
  
  # Keep it simple: 116 / 2 genotypes / 2 sexes / 6 tissues = almost 5 
  smallestGroupSize = 5
  
  keep <- rowSums(counts(tmpDds) >= 10) >= smallestGroupSize
  
  tmpDds <- tmpDds[keep,]
  
  # DESeq Wald test
  DESeq(tmpDds)
}

tissueList <- unique(colData(dds)$Tissue)

ddsList <- lapply(tissueList,analysisDESeq, dds = dds)
names(ddsList) <- tissueList

# PCA plots ----

vsdList <- lapply(ddsList, vst, blind = TRUE)

pcaList <- lapply(vsdList, plotPCA, intgroup=c("Genotype","Gender"))


# Full results ----

resultsList <- lapply(ddsList, results, 
                      contrast = c("Genotype", "KO", "WT"))

resultsList <- lapply(resultsList, function(x){
  x |> 
    as_tibble(rownames = "gene_id") |> 
    left_join(geneMetadata,
              by = join_by(gene_id))
})

names(resultsList) <- tissueList

# When we use name = scaleAge, the log2FoldChange is the increase/decrease
# in expression per unit of scaleAge.

# Save results ----

saveRDS(resultsList, file = "Out/LAKI/DESeq_results_STAR_LAKI.rds")

#Remove environment objects but no functions
rm(list = setdiff(ls(), lsf.str()))



