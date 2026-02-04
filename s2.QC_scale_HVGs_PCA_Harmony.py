import omicverse as ov
import scanpy as sc
import os
import matplotlib.pyplot as plt
import numpy as np

adata_merged = ov.read('$Dir/result/01.input_merge_mouse2human.h5ad')

adata_merged.var['mt'] = adata_merged.var_names.str.startswith('MT-')

sc.pp.calculate_qc_metrics(adata_merged,
                           qc_vars=['mt'],
                           percent_top=None,
                           log1p=False,
                           inplace=True)

min_genes = 200
max_genes = 5000
mito_thresh = 0.2
nUMIs_thresh = 500

adata_merged = adata_merged[
    (adata_merged.obs['n_genes_by_counts'] >= min_genes) &
    (adata_merged.obs['n_genes_by_counts'] <= max_genes) &
    (adata_merged.obs['pct_counts_mt'] <= mito_thresh * 100) &
    (adata_merged.obs['total_counts'] >= nUMIs_thresh)
].copy()

ov.utils.store_layers(adata_merged,layers='counts')

adata_merged=ov.pp.preprocess(adata_merged,mode='shiftlog|pearson',
                       n_HVGs=3000,batch_key='batch')

print(f"保留细胞数: {adata_merged.n_obs}, 基因数: {adata_merged.n_vars}")

adata_merged.raw = adata_merged
adata_merged = adata_merged[:, adata_merged.var.highly_variable_features]

ov.pp.scale(adata_merged)
ov.pp.pca(adata_merged,layer='scaled',n_pcs=50)

adata_merged.obsm["X_mde_pca"] = ov.utils.mde(adata_merged.obsm["scaled|original|X_pca"])

ov.utils.embedding(adata_merged,
                basis='X_mde_pca',frameon='small',
                color=['batch'],show=False)
plt.savefig(
    os.path.join("$Dir/result", "02.before_Harmony_for_each_sample.png"),
    bbox_inches="tight",
    dpi=300
)
plt.close()

adata_merged_harmony=ov.single.batch_correction(adata_merged,batch_key='batch',
                                        methods='harmony',n_pcs=50)

adata_merged.obsm["X_mde_harmony"] = ov.utils.mde(adata_merged.obsm["X_harmony"])
ov.utils.embedding(adata_merged,
                basis='X_mde_harmony',frameon='small',
                color=['batch'],show=False)
plt.savefig(
    os.path.join("$Dir/result", "02.after_Harmony_for_each_sample.png"),
    bbox_inches="tight",
    dpi=300
)
plt.close()

adata_merged.write("$Dir/result/02.QC_scale_HVGs_PCA_Harmony.h5ad")
