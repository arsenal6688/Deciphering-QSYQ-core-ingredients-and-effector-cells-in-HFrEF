import scanpy as sc
import drug2cell as d2c
import blitzgsea as blitz
import pandas as pd
import omicverse as ov
import matplotlib.pyplot as plt
import pickle
import blitzgsea as blitz
import os
import numpy as np
import re

output_dir = '$Dir/result'

adata = ov.read('$Dir/result/04.celltype_manual.h5ad')
adata

with open("$Dir/herb-target_dictionary/QSYQ_herb.pkl", "rb") as f:
    chembl_dict_TCM = pickle.load(f)

filtered_dict = chembl_dict_TCM.copy()
filtered_dict.pop('QSYQ', None)  

sc.pp.normalize_per_cell(adata, counts_per_cell_after = 1e4)
sc.pp.log1p(adata)
adata.raw = adata

adata.obs['combined_categ'] = [
    f"{ct}_{g}"
    for ct, g in zip(adata.obs['major_celltype'], adata.obs['group'])
]

target_celltypes = ['Endothelial cell','Fibroblast','Cardiomyocyte','Macrophage','Perivascular cell']
adata_sub = adata[adata.obs['major_celltype'].isin(target_celltypes)].copy()

print(adata_sub.obs['major_celltype'].value_counts())
print(adata_sub.obs['group'].value_counts())

d2c.score(adata_sub, targets=filtered_dict, method="seurat", use_raw=True)

sc.tl.rank_genes_groups(adata_sub.uns['drug2cell'], method="wilcoxon", groupby="group",groups=['HF'], reference='CTRL')
deg_df = sc.get.rank_genes_groups_df(adata_sub.uns['drug2cell'],group='HF')
filename = os.path.join(output_dir, f'07-1.herb_target_celltypes.csv')
deg_df.to_csv(filename, index=False)

sc.pl.rank_genes_groups_dotplot(adata_sub.uns['drug2cell'],swap_axes=True, key="rank_genes_groups",dendrogram=False,groupby = 'major_celltype',cmap = 'Blues',standard_scale='var',categories_order=['Fibroblast','Cardiomyocyte','Endothelial cell','Perivascular cell','Macrophage'],show=False)
filename = os.path.join(output_dir, f'07-1.herb_target_celltypes_dotplot.svg')
plt.savefig(filename, dpi=300, bbox_inches='tight',format='svg')
plt.close()
