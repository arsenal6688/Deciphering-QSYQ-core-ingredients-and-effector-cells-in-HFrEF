import omicverse as ov
print(f'omicverse version: {ov.__version__}')
import scanpy as sc
print(f'scanpy version: {sc.__version__}')
ov.ov_plot_set()
import matplotlib.pyplot as plt
import os
import pandas as pd
import pickle

adata = ov.read('$Dir/result/04.celltype_manual.h5ad')

macrophage_adata = adata[adata.obs['major_celltype'] == 'Macrophage'].copy()
print(f"提取到巨噬细胞数量: {macrophage_adata.n_obs}")

sc.pp.normalize_total(macrophage_adata)
sc.pp.log1p(macrophage_adata)
sc.pp.highly_variable_genes(macrophage_adata)

sc.pp.pca(macrophage_adata)

sc.pp.neighbors(macrophage_adata, n_neighbors=15, n_pcs=30, use_rep='X_harmony')
sc.tl.leiden(macrophage_adata, key_added="macrophage_subtypes", resolution=0.3)
macrophage_adata.obsm["X_mde"] = ov.utils.mde(macrophage_adata.obsm["X_harmony"])

ov.utils.embedding(macrophage_adata,
                basis='X_mde',
                color=["macrophage_subtypes"],
                title=['Resolution:1'],
                palette=ov.palette()[:],
                   frameon='small',
                   show=False
                   )

small_marker_dict={
    'Mac-LYVE1':['LYVE1','MRC1','MAF'],
    'Mac-IL1β':['PDE4B','C5AR1','CLEC7A','NR4A3','XBP1','NLRP3','IL1B','C3','CLEC4E','EREG','TLR2','HLA-DRA','HLA-DQA1','HLA-DPB1','HLA-A','CD74','HLA-DPA1','B2M','HLA-C','THBS1','HLA-B'],
    'Mac-undefined':['DCN','LUM','COL1A1','COL1A2','CD9','TREM2','LGALS3','LGALS1','ANXA2','S100A6','FTH1']
}

cluster2annotation = {
     '0': 'Mac-LYVE1',
     '1': 'Mac-LYVE1',
     '2': 'Mac-LYVE1',
     '3': 'Mac-IL1β',
     '4': 'Mac-undefined',
     '5': 'Mac-IL1β',
     '6': 'Mac-undefined'
     }
macrophage_adata.obs['major_celltype'] = macrophage_adata.obs['macrophage_subtypes'].map(cluster2annotation).astype('category')

high_contrast_palette = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
    '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5'
    ]

ov.utils.embedding(macrophage_adata,
                   basis='X_mde',
                   color=["macrophage_subtypes", "major_celltype"],
                   title=['Clusters', 'Major Cell types'],
                   palette=high_contrast_palette,
                   wspace=0.55,
                   frameon='small',
                   show=False
)

macrophage_adata.write("$Dir/result/06-1.macrophage_subtypes_manual.h5ad")

valid_markers = {}
for subtype, markers in small_marker_dict.items():
    existing_markers = [m for m in markers if m in macrophage_adata.var.index]
    if existing_markers:
        valid_markers[subtype] = existing_markers

print("有效标记基因字典:")
for subtype, markers in valid_markers.items():
    print(f"{subtype}: {markers}")

sc.pl.dotplot(macrophage_adata,
              groupby="macrophage_subtypes",
              var_names=valid_markers,
              standard_scale="var",
              dendrogram=True,
              show=False
              )

plt.savefig(
    "$Dir/result/06-1.macrophage_subtype_dotplot.svg",
    dpi=300, 
    bbox_inches="tight",
    format='svg'
)
plt.close()

sc.tl.rank_genes_groups(
    macrophage_adata, groupby="macrophage_subtypes", use_raw=False,
    corr_method='benjamini-hochberg',
    method="wilcoxon", key_added="macrophage_subtypes")

sc.tl.filter_rank_genes_groups(
    macrophage_adata,
    min_in_group_fraction=0.1,
    max_out_group_fraction=0.2,
    key="macrophage_subtypes",
    key_added="macrophage_subtypes_filtered",
)

sc.pl.rank_genes_groups_dotplot(
    macrophage_adata, groupby="macrophage_subtypes", dendrogram=True,
    standard_scale="var", n_genes=3, key="macrophage_subtypes_filtered",
    show=False
)
save_path = '$Dir/result/06-1.macrophage_subtype_dotplot_manual_with_marker_diff_genes.png'
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

clusters = sorted(macrophage_adata.obs['macrophage_subtypes'].cat.categories)
with open('$Dir/result/06-1.macrophage_subtype_all_clusters_top10_markers.txt', 'w') as f:
    for cluster in clusters:
        degs = sc.get.rank_genes_groups_df(
            macrophage_adata,
            group=str(cluster),
            key='macrophage_subtypes_filtered'
            ).dropna().head(20)
        f.write(f'=== Cluster {cluster} ===\n')
        f.write('\n'.join(degs['names'].astype(str)) + '\n\n')

adata_raw=macrophage_adata.raw.to_adata()
adata_raw

pathway_dict=ov.utils.geneset_prepare('$Dir/database/CellMarker_Augmented_2021.txt',organism='Human')
adata_aucs=ov.single.pathway_aucell_enrichment(adata_raw,
                                                pathways_dict=pathway_dict,
                                                num_workers=1)
adata_aucs
adata_aucs.obs=adata_raw[adata_aucs.obs.index].obs
adata_aucs.obsm=adata_raw[adata_aucs.obs.index].obsm
adata_aucs.obsp=adata_raw[adata_aucs.obs.index].obsp
adata_aucs

sc.tl.rank_genes_groups(
    adata_aucs, groupby="macrophage_subtypes", use_raw=False,
    method="wilcoxon", key_added="dea_leiden_aucs_res1"
)

sc.pl.rank_genes_groups_dotplot(adata_aucs,groupby='macrophage_subtypes',
                                cmap='RdBu_r',key='dea_leiden_aucs_res1',
                                standard_scale='var',n_genes=3,
                                show=False
                                )

save_path = '$Dir/result/06-1.macrophage_dotplot_manual_with_marker_database.png'
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

sc.tl.umap(macrophage_adata)

sc.pl.umap(macrophage_adata, color='major_celltype',
          title='UMAP Macrophage',
          palette=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
          frameon='small',
           show=False)

plt.savefig(
    os.path.join("$Dir/result", "06-1.macrophage_umap.svg"),
    dpi=300,
    bbox_inches="tight",
    format='svg'
)
plt.close()
