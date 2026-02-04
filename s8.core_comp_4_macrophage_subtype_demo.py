import scanpy as sc
import drug2cell as d2c
import pandas as pd
import omicverse as ov
import pickle
import os
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt

output_dir = '$Dir/result'

adata = ov.read('$Dir/result/06-1.macrophage_subtypes_manual.h5ad')

pkl_files = [
    "$Dir/ingredient-target_dictionary/DanShen.pkl",
    "$Dir/ingredient-target_dictionary/HuangQi.pkl",
    "$Dir/ingredient-target_dictionary/JiangXiang.pkl",
    "$Dir/ingredient-target_dictionary/SanQi.pkl"
]

chembl_dict_TCM = {}
for file in pkl_files:
    with open(file, "rb") as f:
        herb_dict = pickle.load(f)
        chembl_dict_TCM.update(herb_dict)

file_paths = [
    '$Dir/result/07-2.danshen_manual_frac_and_mean_count_filterd_sorted.csv',
    '$Dir/result/07-2.huangqi_manual_frac_and_mean_count_filterd_sorted.csv',
    '$Dir/result/07-2.sanqi_manual_frac_and_mean_count_filterd_sorted.csv',
    '$Dir/result/07-2.jiangxiang_manual_frac_and_mean_count_filterd_sorted.csv'
]

gene_list = []
for file_path in file_paths:
    df = pd.read_csv(file_path)
    macrophage_genes = df[df['cell_type'] == 'Macrophage']['gene'].tolist()
    gene_list.extend(macrophage_genes)

target_celltypes = ['Mac-LYVE1','Mac-IL1β','Mac-undefined']

adata_sub = adata[adata.obs['major_celltype'].isin(target_celltypes)].copy()
cell_counts = adata_sub.obs['major_celltype'].value_counts().reindex(target_celltypes, fill_value=0)
for cell_type, count in cell_counts.items():
    print(f"{cell_type}: {count} 个细胞")

d2c.score(adata_sub, targets=chembl_dict_TCM, method="seurat", use_raw=True)

drug2cell_data = adata_sub.uns['drug2cell'].copy()

apoE_HFD_subset = drug2cell_data[drug2cell_data.obs['group'] == 'HF'].copy()

sc.tl.rank_genes_groups(apoE_HFD_subset, method="wilcoxon", groupby="group", groups=['HF'])

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
        gene_indices = [i for i, gene in enumerate(gene_names) if gene in gene_list]
        if len(gene_indices) > 0:
            cell_type_data_subset = cell_type_data[:, gene_indices]
            if sparse.issparse(cell_type_data_subset):
                fracs = (cell_type_data_subset > 0).mean(axis=0).A1
            else:
                fracs = (cell_type_data_subset > 0).mean(axis=0)
            means = cell_type_data_subset.mean(axis=0)
            if sparse.issparse(means):
                means = means.A1
            selected_genes = [gene_names[i] for i in gene_indices]
            results[cell_type] = {
                'genes': selected_genes,
                'fracs': fracs,
                'means': means,
                'cell_count': np.sum(cell_type_mask)
            }

result_dfs = []

for cell_type in target_celltypes:
    if cell_type in results:
        df_temp = pd.DataFrame({
            'gene': results[cell_type]['genes'],
            'cell_type': cell_type,
            'frac_hfd': results[cell_type]['fracs'],
            'mean_hfd': results[cell_type]['means'],
            'cell_count': results[cell_type]['cell_count']
        })
        result_dfs.append(df_temp)

final_df = pd.concat(result_dfs, ignore_index=True)

max_mean_df = final_df.groupby('gene').apply(
    lambda x: x.loc[x['mean_hfd'].idxmax()]
).reset_index(drop=True)
max_mean_df['cell_type'] = pd.Categorical(max_mean_df['cell_type'],
                                      categories=target_celltypes,
                                      ordered=True)
max_mean_df_sorted = max_mean_df.sort_values('cell_type')
max_mean_df_sorted.to_csv(os.path.join(output_dir, '08-1.core_comp_macrophage_sub_manual_frac_and_mean_count_filterd_sorted.csv'), index=False)

top_genes = []

for cell_type in target_celltypes:
    cell_type_mask = max_mean_df['cell_type'] == cell_type
    cell_type_data = max_mean_df[cell_type_mask]
    cell_type_sorted = cell_type_data.sort_values('mean_hfd', ascending=False)
    top_genes.extend(cell_type_sorted['gene'].tolist())

top_genes = list(dict.fromkeys(top_genes))

print("新生成的top_genes列表：")
for i, gene in enumerate(top_genes, 1):
    print(f"{i}. {gene}")

print(f"top_genes列表的个数为: {len(top_genes)}")

sc.pl.rank_genes_groups_dotplot(
    apoE_HFD_subset,
    swap_axes=True,
    key="rank_genes_groups",
    dendrogram=False,
    groupby='major_celltype',
    cmap='Blues',
    standard_scale='var',
    var_names=top_genes,
    categories_order=['Mac-LYVE1','Mac-IL1β','Mac-undefined'],
    show=False
)
filename = os.path.join(output_dir, f'08-1.core_comp_macrophage_sub_HFD_dotplot.svg')
plt.savefig(filename, dpi=300, bbox_inches='tight',format='svg')
plt.close()
