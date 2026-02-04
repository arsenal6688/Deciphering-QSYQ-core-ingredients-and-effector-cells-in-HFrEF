import omicverse as ov
print(f'omicverse version: {ov.__version__}')
import scanpy as sc
print(f'scanpy version: {sc.__version__}')
ov.ov_plot_set()
import matplotlib.pyplot as plt
import os
import pandas as pd
import pickle

adata = ov.read('$Dir/result/03.cluster.h5ad')
adata

adata_raw=adata.raw.to_adata()
adata_raw

small_marker_dict={
    'Cardiomyocytes':['RYR2','TTN','FGF12','RBM20','CTNNA3','SORBS2','MLIP','FHL2','CDH2','TECRL','LDB3','MYBPC3','NEBL','MYOM1','FHOD3','AKAP6','MYOM2','EFNA5','SLC8A1','TRDN-AS1'],
    'Fibroblasts':['ACSM3','NEGR1','ABCA9','DCN','CDH19','BICC1','COL6A3','ABCA6','PID1','FBN1','COL8A1','BNC2','ABCA8','ELN','SCN7A','ABCA10','RORA','GLIS3','LINC01088','COL15A1'],
    'Endothelium':['ST6GALNAC3','VWF','LDB2','ANO2','BTNL9','PTPRB','FLT1','ADGRF5','CCDC85A','F8','EMCN','SEC14L1','EGFL7','CYYR1','LNX1','RASGRF2','AQP1','ITGA6','PITPNC1','ST8SIA6'],
    'Perivascular cells':['GUCY1A2','EGFLAM','NR2F2-AS1','RGS5','KCNAB1','AC012409.2','NOTCH3','BMP5','PTH1R','TBX2','NR2F2','CD38','LINC02237','AGT','AGAP2','GJA4','LINC01197','FHL5','MRVI1','CPE'],
    'Macrophages':['F13A1','RBM47','MRC1','FMN1','CD163','STAB1','ATP8B4','MS4A6A','SMAP2','MSR1','TBXAS1','DOCK2','RAB31','KYNU','MS4A4E','MARCH1','DAPK1','MAF','FYB1','DOCK8'],
    'Lymphocytes':['SKAP1','IL7R','THEMIS','BCL11B','ITK','CD247','CAMK4','SAMD3','SLFN12L','ITGA4','LINC01934','RHOH','PRKCQ','CD96','IKZF3','PBX4','LINC00861','CD2','RIPOR2','SLAMF6'],
    'Neurons':['L1CAM','FXYD3','SLITRK5','PTPRZ1','PLP1','SHISA9','XKR4','INSC','NRXN1','CHL1','KCNK5','PRRG4','S100B','SLC5A7','WNT6','CADM3','U91319.1','LYPD6','GRIK3','SOX10']
}

smarker_genes_in_data = dict()
for ct, markers in small_marker_dict.items():
    markers_found = list()
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    smarker_genes_in_data[ct] = markers_found

del_markers = list()
for ct, markers in smarker_genes_in_data.items():
    if markers==[]:
        del_markers.append(ct)

for ct in del_markers:
    del smarker_genes_in_data[ct]

smarker_genes_in_data

with open('$Dir/result/04.smarker_genes_in_data.txt', 'w') as f:
    for cell_type, markers in smarker_genes_in_data.items():
        f.write(f"{cell_type}: {','.join(markers)}\n")

sc.pl.dotplot(
    adata,
    groupby="leiden_res1",
    var_names=smarker_genes_in_data,
    dendrogram=True,
    standard_scale="var",
    show=False
)
save_path = '$Dir/result/04.dotplot_manual_with_marker_dict.png'
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

adata.uns['log1p']['base']=None
sc.tl.rank_genes_groups(
    adata, groupby="leiden_res1", use_raw=False,
    method="wilcoxon", key_added="dea_leiden_res1"
)

sc.tl.filter_rank_genes_groups(
    adata,
    min_in_group_fraction=0.1,
    max_out_group_fraction=0.2,
    key="dea_leiden_res1",
    key_added="dea_leiden_res1_filtered",
)

sc.pl.rank_genes_groups_dotplot(
    adata, groupby="leiden_res1", dendrogram=True,
    standard_scale="var", n_genes=3, key="dea_leiden_res1_filtered",
    show=False
)
save_path = '$Dir/result/04.dotplot_manual_with_marker_diff_genes.png'
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

clusters = sorted(adata.obs['leiden_res1'].cat.categories)
with open('$Dir/result/04.all_clusters_top10_markers.txt', 'w') as f:
    for cluster in clusters:
        degs = sc.get.rank_genes_groups_df(
            adata,
            group=str(cluster),
            key='dea_leiden_res1_filtered'
        ).dropna().head(10)
        f.write(f'=== Cluster {cluster} ===\n')
        f.write('\n'.join(degs['names'].astype(str)) + '\n\n')

pathway_dict=ov.utils.geneset_prepare('$Dir/CellMarkers/CellMarker_Augmented_2021.txt',organism='Human')

adata_aucs=ov.single.pathway_aucell_enrichment(adata_raw,
                                                pathways_dict=pathway_dict,
                                                num_workers=1)
adata_aucs

adata_aucs.obs=adata_raw[adata_aucs.obs.index].obs
adata_aucs.obsm=adata_raw[adata_aucs.obs.index].obsm
adata_aucs.obsp=adata_raw[adata_aucs.obs.index].obsp
adata_aucs

sc.tl.rank_genes_groups(
    adata_aucs, groupby="leiden_res1", use_raw=False,
    method="wilcoxon", key_added="dea_leiden_aucs_res1"
)
sc.pl.rank_genes_groups_dotplot(adata_aucs,groupby='leiden_res1',
                                cmap='RdBu_r',key='dea_leiden_aucs_res1',
                                standard_scale='var',n_genes=3,
                                show=False)
save_path = '$Dir/result/04.dotplot_manual_with_marker_database.png'
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

output_path = '$Dir/result/04.marker_database_cluster_markers.txt'
with open(output_path, 'w') as f:
    for cluster in adata_aucs.uns['dendrogram_leiden_res1']['categories_ordered']:
        special_cluster = str(cluster)
        degs = sc.get.rank_genes_groups_df(
            adata_aucs,
            group=special_cluster,
            key='dea_leiden_aucs_res1'
        ).dropna()
        f.write(f'{special_cluster}: {"|".join(degs.names[:2].tolist())}\n')

cluster2annotation = {
     '0': 'Fibroblast',
     '1': 'Fibroblast',
     '2': 'Cardiomyocyte',
     '3': 'Macrophage',
     '4': 'Cardiomyocyte',
     '5': 'Endothelial cell',
     '6': 'Cardiomyocyte',
     '7': 'Fibroblast',
     '8': 'Perivascular cell',
     '9': 'Endothelial cell',
     '10': 'Lymphoid cell',
     '11': 'Cardiomyocyte',
     '12': 'Perivascular cell',
     '13': 'Endothelial cell',
     '14': 'Fibroblast',
     '15': 'Cardiomyocyte',
     '16': 'Macrophage',
     '17': 'Fibroblast',
     '18': 'Smooth muscle cell',
     '19': 'Endothelial cell',
     '20': 'Cardiomyocyte',
     '21': 'Endothelial cell',
     '22': 'Endothelial cell',
     '23': 'Unknown',
     '24': 'Neurons',
     '25': 'Fibroblast',
     '26': 'Unknown',   
     '27': 'Perivascular cell',
     '28': 'Fibroblast'
}
adata.obs['major_celltype'] = adata.obs['leiden_res1'].map(cluster2annotation).astype('category')

ov.utils.embedding(adata,
                basis='X_mde',
                color=["leiden_res1","major_celltype"],
                title=['Clusters','Major Cell types'],
                palette=ov.palette()[:],wspace=0.55,
                show=False,frameon='small',)

plt.savefig(
    os.path.join("$Dir/result", "04.manual_major_ct.mde.png"),
    dpi=300,
    bbox_inches="tight"
)
plt.close()

if 'dea_leiden_res1_filtered' in adata.uns:
    with open('$Dir/result/04.dea_results.pkl', 'wb') as f:
        pickle.dump(adata.uns['dea_leiden_res1_filtered'], f)
    dea_backup = adata.uns.pop('dea_leiden_res1_filtered')

adata.write("$Dir/result/04.celltype_manual.h5ad")

sc.tl.umap(adata)

sc.pl.umap(adata, color='major_celltype',
          title='UMAP Major Cell Types',
          palette=ov.palette()[:], 
          frameon='small',
           show=False
)

plt.savefig(
    os.path.join("$Dir/result", "04.umap_major_celltypes.svg"),
    dpi=300,
    bbox_inches="tight",
    format='svg'
)
plt.close()

sc.tl.tsne(adata)

sc.pl.tsne(
    adata,
    color='major_celltype', 
    title='t-SNE Major Cell Types',
    palette=ov.palette()[:],  
    show=False,
    frameon='small',
)

plt.savefig(
    os.path.join("$Dir/result", "04.tsne_major_celltypes.svg"),
    dpi=300,
    bbox_inches="tight",
    format='svg'
)
plt.close()
