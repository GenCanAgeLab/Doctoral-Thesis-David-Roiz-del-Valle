#!/usr/bin/env python3

# Library import

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path


# Path options
project_dir = Path(__file__).parent.parent
results_dir = project_dir / "Out" / "results.tsv"
out_path = project_dir / "Out" / "single_cell_targets.pdf"
out_path2 = project_dir / "Out" / "single_cell_targets_24m.pdf"

# Scanpy options
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

#better export text to PDF
plt.rcParams["pdf.use14corefonts"] = True

# Importing data

results = pd.read_table(results_dir, sep='\t',decimal='.')

Bulk_genes = results[(results["padj_Bulk"] < 0.05) &
                     (results["log2FoldChange_Bulk"] > 2) &
                     (results["isTarget_Bulk"] == True) &
                     (results["gene_type"] == "protein_coding") &
                     (results["tag"].isna())]['gene_name']
CM_genes = results[(results["padj_CMspec"] < 0.05) &
                   (results["log2FoldChange_CMspec"] > 1.5) &
                   (results["isTarget_Bulk"] == True) &
                   (results["gene_type"] == "protein_coding") &
                   (results["tag"].isna())]['gene_name']

GOI = set(CM_genes).union(set(Bulk_genes))

GOI = list(GOI)
GOI = sorted(GOI)

heart_sc = sc.read_h5ad('/data/miR29/resources/tabula-muris-senis-facs-processed-official-annotations-Heart.h5ad')

heart_sc_3m = heart_sc[heart_sc.obs['age'] == '3m']
heart_sc_24m = heart_sc[heart_sc.obs['age'] == '24m']

GOI_3m = [gene for gene in GOI if gene in heart_sc_3m.raw.var_names]

sc.pl.matrixplot(heart_sc_3m, GOI_3m, 'cell_ontology_class', dendrogram=False, standard_scale= 'var', cmap='Blues', colorbar_title='Gene scaled\nexpression',swap_axes=True)
plt.savefig(out_path, bbox_inches="tight")

GOI_24m = [gene for gene in GOI if gene in heart_sc_24m.raw.var_names]

sc.pl.matrixplot(heart_sc_24m, GOI_24m, 'cell_ontology_class', dendrogram=False, standard_scale= 'var', cmap='Blues', colorbar_title='Gene scaled\nexpression',swap_axes=True)
plt.savefig(out_path2, bbox_inches="tight")
