import omicverse as ov
print(f'omicverse version: {ov.__version__}')
import scanpy as sc
print(f'scanpy version: {sc.__version__}')
ov.ov_plot_set()
import os
import matplotlib.pyplot as plt

adata = ov.read('$Dir/result/02.QC_scale_HVGs_PCA_Harmony.h5ad')

sc.pp.neighbors(adata, n_neighbors=15,
                n_pcs=30,use_rep='X_harmony')

sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)

adata.obsm["X_mde"] = ov.utils.mde(adata.obsm["X_harmony"])  

ov.utils.embedding(adata,
                basis='X_mde',
                color=["leiden_res1"],
                title=['Resolution:1'],
                palette=ov.palette()[:],
                show=False,frameon='small',)

plt.savefig(
    os.path.join("$Dir/result", "03.leiden_resolution.png"),
    dpi=300,
    bbox_inches="tight"
)
plt.close()

adata.write("$Dir/result/03.cluster.h5ad")
