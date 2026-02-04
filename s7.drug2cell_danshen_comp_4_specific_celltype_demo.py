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
from scipy import sparse
import pandas as pd

output_dir = '$Dir/result'

adata = ov.read('$Dir/result/04.celltype_manual.h5ad')
adata

with open("$Dir/ingredient-target_dictionary/DanShen.pkl", "rb") as f:
    chembl_dict_TCM = pickle.load(f)

num_keys = len(chembl_dict_TCM)
num_keys

sc.pp.normalize_per_cell(adata, counts_per_cell_after = 1e4)
sc.pp.log1p(adata)
adata.raw = adata

target_celltypes = ['Fibroblast','Cardiomyocyte','Endothelial cell','Perivascular cell','Macrophage']

adata_sub = adata[adata.obs['major_celltype'].isin(target_celltypes)].copy()
print(adata_sub.obs['major_celltype'].value_counts())
print(adata_sub.obs['group'].value_counts())

d2c.score(adata_sub, targets=chembl_dict_TCM, method="seurat", use_raw=True)
sc.tl.rank_genes_groups(adata_sub.uns['drug2cell'], method="wilcoxon", groupby="group",groups=['HF'],reference='CTRL')
deg_df = sc.get.rank_genes_groups_df(adata_sub.uns['drug2cell'],group='HF')
filename = os.path.join(output_dir, f'07-2.danshen_comp_rank_genes.csv')
deg_df.to_csv(filename, index=False)

sc.pl.rank_genes_groups_dotplot(adata_sub.uns['drug2cell'],swap_axes=True, key="rank_genes_groups",dendrogram=False,groupby = 'major_celltype',cmap = 'Blues',standard_scale='var',categories_order=['Fibroblast','Cardiomyocyte','Endothelial cell','Perivascular cell','Macrophage'],show=False)
filename = os.path.join(output_dir, f'07-2.danshen_target_celltypes_components_dotplot.svg')
plt.savefig(filename, dpi=300, bbox_inches='tight',format='svg')
plt.close()

d2c_data = adata_sub.uns['drug2cell']
cell_types = d2c_data.obs['major_celltype'].values
groups = d2c_data.obs['group'].values
gene_names = d2c_data.var_names

hfd_mask = groups == 'HF'

results = {}
for cell_type in target_celltypes:
    cell_type_mask = (cell_types == cell_type) & hfd_mask
    if np.sum(cell_type_mask) > 0:
        cell_type_data = d2c_data.X[cell_type_mask, :]
        if sparse.issparse(cell_type_data):
            fracs = (cell_type_data > 0).mean(axis=0).A1
        else:
            fracs = (cell_type_data > 0).mean(axis=0)
        means = cell_type_data.mean(axis=0)
        if sparse.issparse(means):
            means = means.A1
        results[cell_type] = {
                'fracs': fracs,
                'means': means,
                'cell_count': np.sum(cell_type_mask)
        }

result_dfs = []
for cell_type in target_celltypes:
    if cell_type in results:
        df_temp = pd.DataFrame({
            'gene': gene_names,
            'cell_type': cell_type,
            'frac_hfd': results[cell_type]['fracs'],
            'mean_hfd': results[cell_type]['means'],
            'cell_count': results[cell_type]['cell_count']
        })
        result_dfs.append(df_temp)

final_df = pd.concat(result_dfs, ignore_index=True)
filename = os.path.join(output_dir, f'07-2.danshen_manual_frac_and_mean_count.csv')
final_df.to_csv(filename, index=False)

for cell_type in target_celltypes:
    if cell_type in results:
        print(f"{cell_type}: {results[cell_type]['cell_count']}")

sorted_results = final_df.sort_values(['cell_type', 'mean_hfd'], ascending=[True, False])
file_path = os.path.join(output_dir, '07-2.danshen_manual_frac_and_mean_count.csv')
df = pd.read_csv(file_path)

filtered_df = df[
    (df.groupby('gene')['frac_hfd'].transform('max') >= 0.25) &
    (df.groupby('gene')['mean_hfd'].transform('max') >= 0)
]

max_mean_df = filtered_df.groupby('gene').apply(
    lambda x: x.loc[x['mean_hfd'].idxmax()]
).reset_index(drop=True)
max_mean_df['cell_type'] = pd.Categorical(max_mean_df['cell_type'],
                                      categories=target_celltypes,
                                      ordered=True)
max_mean_df_sorted = max_mean_df.sort_values('cell_type')
max_mean_df_sorted.to_csv(os.path.join(output_dir, '07-2.danshen_manual_frac_and_mean_count_filterd_sorted.csv'), index=False)

top_genes = []
for cell_type in target_celltypes:
    cell_type_mask = max_mean_df_sorted['cell_type'] == cell_type
    top_genes.extend(max_mean_df_sorted[cell_type_mask]['gene'].tolist())
top_genes = list(dict.fromkeys(top_genes))

sc.pl.rank_genes_groups_dotplot(adata_sub.uns['drug2cell'],swap_axes=True, key="rank_genes_groups", dendrogram=False,groupby = 'major_celltype',cmap = 'Blues',standard_scale='var',categories_order=['Fibroblast','Cardiomyocyte','Endothelial cell','Perivascular cell','Macrophage'],groups=['HF'],var_names=top_genes,show=False)
filename = os.path.join(output_dir, f'07-2.danshen_top_genes_targeted_cells_dotplot.svg')
plt.savefig(filename, dpi=300, bbox_inches='tight',format='svg')
plt.close()
