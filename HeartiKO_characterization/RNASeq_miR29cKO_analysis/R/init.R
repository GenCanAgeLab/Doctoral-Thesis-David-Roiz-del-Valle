# Loading libraries ------------------
suppressMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(pheatmap)
  library(tximeta)
  library(DESeq2)
  library(fgsea)
  library(fgsea)
  library(nVennR)
  library(tibble)
  library(scales)
  library(purrr)
}
)

# Multicore parameters ---------

library(BiocParallel)
register(MulticoreParam(8))