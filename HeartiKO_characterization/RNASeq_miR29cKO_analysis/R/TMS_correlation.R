# Create GOI set ----

## Load RNASeqdata

results <- read_delim('Out/results.tsv')

## Filter and GOI

Bulk_genes <- results |>
  filter(log2FoldChange_Bulk > 2 & padj_Bulk < 0.05) |>
  filter(isTarget_Bulk == T) |>
  pull(gene_name)

CMspec_genes <- results |>
  filter(log2FoldChange_CMspec > 2 & padj_CMspec < 0.05) |>
  filter(isTarget_Bulk == T) |>
  pull(gene_name)

GOIs <- union(Bulk_genes, CMspec_genes)

# Load TMS data ----

## Load miRNA expression (reads per million reads)
mirexp <- read_csv('/data/miR29/resources/GSE217458_TMS_rpmm_Complete_Filtered_global_GEO_18012023.csv.gz')

## Transform

mirexp <- mirexp |>
  filter(str_detect(ncRNA,"mmu-miR-29[abc]{1}-3p")) |>
  mutate(ncRNA = str_remove(ncRNA, "mmu-")) |>
  pivot_longer(-ncRNA, names_to = "Sample", values_to = "miRexp") |>
  pivot_wider(names_from = ncRNA, values_from = miRexp)

# Load TMS bulk RNASeq counts ----

# BiocManager::install("TabulaMurisSenisData")
# require(TabulaMurisSenisData)
# tms_bulk <- TabulaMurisSenisBulk()
# Create DESeq object for transformation
# dds <- DESeqDataSet(tms_bulk, ~organ)
# # Filter out low quality samples and very lowly expressed genes
# dds <- dds[,colSums(counts(dds))> 1E6]
# dds <- estimateSizeFactors(dds)
# dds <- dds[rowSums(counts(dds,normalized=T)) > 50,]
# vsd <- vst(dds,blind = T)
# saveRDS(vsd, '../resources/TMS_Bulk.rds')

vsd <- readRDS('/data/miR29/resources/TMS_Bulk.rds')

## Extract GOI and arrange ----

colData <- colData(vsd) |>
  as_tibble() |>
  dplyr::select(Sample = Sample.name,
                Name = source.name,
                age = characteristics..age,
                sex = characteristics..sex,
                organ)

genexp <- assay(vsd) |>
  as_tibble(rownames = "gene_name") |>
  filter(gene_name %in% GOIs) |>
  pivot_longer(-gene_name, names_to = "Sample", values_to = "genexp") |>
  left_join(colData) |>
  mutate(Sample = paste(Name, sex, age, sep = "-"))

## Correlation ----
  
correlation <- mirexp |>
  rowwise(Sample) |>
  summarise(miR29 = sum(
    c_across(everything())
    )) |>
  full_join(genexp) |>
  drop_na() |>
  group_by(gene_name, organ) |>
  summarize(cor = cor(log2(miR29),
                      genexp),
            pval = cor.test(log2(miR29),
                            genexp)[['p.value']]) |>
  mutate(adj.pval = p.adjust(pval, "bonferroni")) |>
  arrange(gene_name)

# Save correlation ----
write_delim(correlation, file = "Out/TMS_miR29_correlation.tsv",
            delim = "\t")


# Plot correlation ----

ggplot(correlation |>
         filter(organ == "Heart"),
       aes(y=gene_name,x=organ,color=cor,size=-log10(adj.pval)))+
  geom_point()+
  scale_color_gradient2('Pearson',low=muted("blue"),high=muted("red"))+
  coord_fixed()+
  theme_bw(10)+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position = 'top')

ggsave("Out/miR_29_gene_correlation.pdf", width = 6, height = 10)

# Gene expression over time ---- 

plotdt <- genexp |>
  group_by(gene_name, organ) |>
  mutate(z = scale(genexp),
         age = as.factor(
           as.numeric(age)
         )) |>
  filter(organ == "Heart")
  
ggplot(plotdt, aes(x=age,y=fct_rev(gene_name)))+
  geom_tile(aes(fill=z),color='gray')+
  coord_fixed()+
  scale_fill_gradient2('Scaled expression',low=muted("blue"),high=muted("red"))+
  theme_bw() +
  labs(x="Age",y="Gene",fill="Scaled expression")

ggsave("Out/TMS_gene_expression.pdf", width = 6, height = 10)


