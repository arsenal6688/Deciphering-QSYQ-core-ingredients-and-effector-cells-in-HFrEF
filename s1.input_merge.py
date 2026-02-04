import scanpy as sc
import pandas as pd
import numpy as np

adata_merged = sc.AnnData()

#-GSE302337
adata_GSM9102205 = sc.read_10x_h5(
    "$Dir/GSE302337/GSM9102205_c9_raw_feature_bc_matrix_FPR_0.05_filtered.h5",
)
adata_GSM9102205.var_names_make_unique()
adata_GSM9102205.obs['batch'] = 'GSM9102205'
adata_GSM9102205.obs['group'] = 'CTRL'

adata_GSM9102206 = sc.read_10x_h5(
    "$Dir/GSE302337/GSM9102206_c22_raw_feature_bc_matrix_FPR_0.05_filtered.h5",
)
adata_GSM9102206.var_names_make_unique()
adata_GSM9102206.obs['batch'] = 'GSM9102206'
adata_GSM9102206.obs['group'] = 'CTRL'

adata_GSM9102207 = sc.read_10x_h5(
    "$Dir/GSE302337/GSM9102207_v30_raw_feature_bc_matrix_FPR_0.05_filtered.h5",
)
adata_GSM9102207.var_names_make_unique()
adata_GSM9102207.obs['batch'] = 'GSM9102207'
adata_GSM9102207.obs['group'] = 'HF'

adata_GSM9102208 = sc.read_10x_h5(
    "$Dir/GSE302337/GSM9102208_v51_601_raw_feature_bc_matrix_FPR_0.05_filtered.h5",
)
adata_GSM9102208.var_names_make_unique()
adata_GSM9102208.obs['batch'] = 'GSM9102208'
adata_GSM9102208.obs['group'] = 'HF'

adata_GSM9102209 = sc.read_10x_h5(
    "$Dir/GSE302337/GSM9102209_v56_raw_feature_bc_matrix_FPR_0.05_filtered.h5",
)
adata_GSM9102209.var_names_make_unique()
adata_GSM9102209.obs['batch'] = 'GSM9102209'
adata_GSM9102209.obs['group'] = 'HF'

adata_GSM9102210 = sc.read_10x_h5(
    "$Dir/GSE302337/GSM9102210_v109_710_raw_feature_bc_matrix_FPR_0.05_filtered.h5",
)
adata_GSM9102210.var_names_make_unique()
adata_GSM9102210.obs['batch'] = 'GSM9102210'
adata_GSM9102210.obs['group'] = 'HF'


#-GSE161470
adata_GSM4907712 = sc.read_10x_mtx(
    "$Dir/GSE161470/GSM4907712",
    var_names="gene_symbols",
    cache=True
)
adata_GSM4907712.obs['batch']='GSM4907712'
adata_GSM4907712.obs['group'] ='CTRL'

adata_GSM4907713 = sc.read_10x_mtx(
    "$Dir/GSE161470/GSM4907713",
    var_names="gene_symbols",
    cache=True
)
adata_GSM4907713.obs['batch']='GSM4907713'
adata_GSM4907713.obs['group'] ='CTRL'

adata_GSM4907714 = sc.read_10x_mtx(
    "$Dir/GSE161470/GSM4907714",
    var_names="gene_symbols",
    cache=True
)
adata_GSM4907714.obs['batch']='GSM4907714'
adata_GSM4907714.obs['group'] ='CTRL'

adata_GSM4907715 = sc.read_10x_mtx(
    "$Dir/GSE161470/GSM4907715",
    var_names="gene_symbols",
    cache=True
)
adata_GSM4907715.obs['batch']='GSM4907715'
adata_GSM4907715.obs['group'] ='CTRL'

adata_GSM5292817 = sc.read_10x_mtx(
    "$Dir/GSE161470/GSM5292817",
    var_names="gene_symbols",
    cache=True
)
adata_GSM5292817.obs['batch']='GSM5292817'
adata_GSM5292817.obs['group'] ='HF'

adata_merged = sc.concat([adata_GSM9102205, adata_GSM9102206, adata_GSM9102207, adata_GSM9102208, adata_GSM9102209, adata_GSM9102210, adata_GSM4907712,adata_GSM4907713,adata_GSM4907714,adata_GSM4907715,adata_GSM5292817], axis=0, join='outer', merge='same',
              index_unique="-",keys=['GSM9102205','GSM9102206','GSM9102207','GSM9102208','GSM9102209','GSM9102210','GSM4907712','GSM4907713','GSM4907714','GSM4907715','GSM5292817'])
adata_merged.obs['batch'].unique()

total_cells = adata_merged.shape[0]
original_genes = adata_merged.var_names.tolist()
n_original = len(original_genes)
print("Original gene number:", n_original)

adata_merged.write("$Dir/result/01.input_merge_mouse2human.h5ad")


