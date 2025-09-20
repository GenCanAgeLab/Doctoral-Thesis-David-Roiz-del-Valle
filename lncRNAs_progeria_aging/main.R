# Run subscripts ----

source("R/init.R")
source("R/functions.R")
source("R/databases.R")

#GSEA
source("R/LAKI/GSEA_analysis.R")


## Sex divided analysis ----

source("R/LAKI/Sex_divided/DESeq_sexDivided.R")
# source("R/LAKI/Cis_reg.R")
source("R/LAKI/Sex_divided/ncRNA_analysis_sexDiv.R")
source("R/LAKI/Sex_divided/ncRNA_clustering_sexDiv.R")
source("R/LAKI/LAKI_metaanalysis.R")

## STAR ----
source("R/LAKI/STAR/DESeq_STAR.R")
source("R/LAKI/STAR/GSEA_analysis_STAR.R")

# TMS ----

source("R/TMS/TMS_DESeq_Corr.R")
source("R/TMS/TMS_ncRNA_analysis.R")
source("R/TMS/TMS_metaanalysis.R")
# source("R/TMS/TMS_gene_trajectories.R")

# Zhou ----

source("R/Zhou/Zhou_DESeq_Corr.R")
source("R/Zhou/Zhou_ncRNA_analysis.R")


# Metaanalysis ----

source("R/metaanalysis_physioAging.R")