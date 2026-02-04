import scanpy as sc
import gseapy as gp
from gseapy.plot import barplot, dotplot
import matplotlib.pyplot as plt
import omicverse as ov
import pandas as pd
import numpy as np
import pickle

output_dir = '$Dir/result'

adata_merged = ov.read('$Dir/result/04.celltype_manual.h5ad')

adata = adata_merged[adata_merged.obs['major_celltype'] == 'Macrophage'].copy()
print(f"提取到巨噬细胞数量: {adata.n_obs}")

groups = ['CTRL', 'HF']

sc.tl.rank_genes_groups(
    adata,
    groupby='group',          
    groups=['HF'],           
    reference='CTRL',        
    method='wilcoxon',       
    corr_method='benjamini-hochberg',  
    tie_correct=True,         
    pts=True                  
)

deg_df = sc.get.rank_genes_groups_df(adata, group='HF')

filter_mask = (
    (deg_df['pvals_adj'] < 0.05) &
    (deg_df['logfoldchanges'].abs() > 1)

)
deg_genes = deg_df[filter_mask]['names'].tolist()

print(f"筛选获得显著差异基因数量: {len(deg_genes)}")

file_paths = [
    '$Dir/result/07-2.danshen_manual_frac_and_mean_count_filterd_sorted.csv',
    '$Dir/result/07-2.huangqi_manual_frac_and_mean_count_filterd_sorted.csv',
    '$Dir/result/07-2.sanqi_manual_frac_and_mean_count_filterd_sorted.csv',
    '$Dir/result/07-2.jiangxiang_manual_frac_and_mean_count_filterd_sorted.csv'
]

key_list = []
for file_path in file_paths:
    df = pd.read_csv(file_path)
    macrophage_genes = df[df['cell_type'] == 'Macrophage']['gene'].tolist()
    key_list.extend(macrophage_genes)
key_set = set(key_list)

pkl_files = [
    "$Dir/ingredient-target_dictionary/DanShen.pkl",
    "$Dir/ingredient-target_dictionary/HuangQi.pkl",
    "$Dir/ingredient-target_dictionary/JiangXiang.pkl",
    "$Dir/ingredient-target_dictionary/SanQi.pkl"
]

all_pkl_keys = set()
for pkl_path in pkl_files:
    with open(pkl_path, 'rb') as f:
        pkl_data = pickle.load(f)
    all_pkl_keys.update(pkl_data.keys())

global_match_keys = all_pkl_keys & key_set

all_genes = []
for pkl_path in pkl_files:
    with open(pkl_path, 'rb') as f:
        pkl_data = pickle.load(f)
    for key in global_match_keys:
        if key in pkl_data:
            val = pkl_data[key]
            all_genes.extend(val)

unique_genes = list(set(all_genes))
print(f"成分提取的基因数量: {len(unique_genes)}")

overlap_genes = list(set(unique_genes) & set(deg_genes))
print(f"overlap基因数量: {len(overlap_genes)}")

local_gene_sets = [
    "$Dir/database/KEGG_2021_Human.gmt",          
    "$Dir/database/Reactome_Pathways_2024.gmt",   
    "$Dir/database/WikiPathways_2024_Human.gmt"   
]

enr = gp.enrichr(
    gene_list=overlap_genes,
    gene_sets=local_gene_sets,
    organism='human',  
    verbose=True  
)

enr.results.to_csv(f'{output_dir}/09-1.enrichment_macrophage.txt', index=False)

sig_pathways = enr.results[
    (enr.results['Adjusted P-value'] < 0.05) &
    (enr.results['Odds Ratio'] > 1) &
    (enr.results['Overlap'].str.split('/').apply(lambda x: int(x[0]) > 5))
].sort_values('Adjusted P-value')

top_20_by_pvalue = sig_pathways.head(20)
sig_pathways = top_20_by_pvalue.sort_values('Odds Ratio', ascending=False)

sig_pathways.to_csv(f'{output_dir}/09-1.enrichment_macrophage_sig_top20_pathways.csv', index=False)
print(f"显著富集的top20: {len(sig_pathways)}")

df = sig_pathways.copy()

df['Gene_ratio'] = df['Overlap'].str.split('/').apply(
    lambda x: int(x[0])/int(x[1])
)
df['-log10(padj)'] = -np.log10(df['Adjusted P-value'])

sig_df = df[
    (df['Adjusted P-value'] < 0.05) &
    (df['Odds Ratio'] > 1)
].sort_values('Odds Ratio').head(20)

sig_df_sorted = sig_df.sort_values(by='Odds Ratio', ascending=False).head(20)
sig_df_sorted.to_csv(f'{output_dir}/09-1.enrichment_macrophage_sig_top20_pathways_filtered.txt', index=False)

plt.figure(figsize=(12,8))
sc = plt.scatter(
    x=sig_df['Odds Ratio'],
    y=sig_df['Term'].str[:50],  
    s=sig_df['Odds Ratio']*10,
    c=sig_df['-log10(padj)'],
    cmap='viridis_r',
    alpha=0.7,
    edgecolor='gray'
)

plt.colorbar(sc, label='-log10(adj.P)')

plt.xlabel('Odds ratio', fontsize=12)

plt.grid(True, linestyle=':', alpha=0.3)
plt.tight_layout()

plt.savefig(f'{output_dir}/09-1.enrichment_macrophage.svg', dpi=300, bbox_inches='tight',format='svg')
plt.show()
