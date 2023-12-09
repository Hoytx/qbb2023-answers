#!/usr/bin/env python

import sys

import scanpy as sc
import numpy
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context


# Read the 10x dataset filtered down to just the highly-variable genes
adata = sc.read_h5ad("variable_data.h5")
adata.uns['log1p']['base'] = None # This is needed due to a bug in scanpy

#step 1.1
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

#step 1.2
sc.tl.leiden(adata)

#step 1.3
sc.tl.umap(adata, maxiter = 900)

sc.tl.tsne(adata)

fig, axes = plt.subplots(ncols=2)

sc.pl.umap(adata, color='leiden', ax = axes[0], title = "UMAP", show=False)
sc.pl.tsne(adata, color='leiden', ax = axes[1], title = "t-SNE", show=False)

fig.tight_layout()
fig.savefig("UMAP_and_tSNE_Plot.png", dpi=1000)

#step 2.1
wilcoxon_adata = sc.tl.rank_genes_groups(adata, groupby='leiden', method = "wilcoxon", use_raw = True, copy = True)

logreg_adata = sc.tl.rank_genes_groups(adata, groupby='leiden', method = "logreg", use_raw = True, copy = True)

#making two plots

sc.pl.rank_genes_groups(wilcoxon_adata, sharey=False, show=False, use_raw=True, title = "Wilcoxon", n_genes = 25, save = "_wilcoxon.png")
sc.pl.rank_genes_groups(logreg_adata, sharey=False, show=False, use_raw=True, title = "logreg", n_genes = 25, save = "_logreg.png")


leiden = adata.obs['leiden']
umap = adata.obsm['X_umap']
tsne = adata.obsm['X_tsne']
adata = sc.read_h5ad('filtered_data.h5')
adata.obs['leiden'] = leiden
adata.obsm['X_umap'] = umap
adata.obsm['X_tsne'] = tsne

with rc_context({'figure.figsize': (3, 3)}):
    sc.pl.umap(adata, color=['CD79A', 'MS4A1', 'GNLY', 'NKG7', 'CD3D', 'CST3', 'LYZ'], s=50, frameon=False, ncols=4, vmax='p99', save="Cell_type_determination.png")

adata.rename_categories('leiden', ['T-cells', 'Myeloid', 'B-cells', 3, 4, 'NK', 6, 7])

sc.pl.umap(adata, color='leiden', title = "Cell Types", save = "_Cell_Types.png")