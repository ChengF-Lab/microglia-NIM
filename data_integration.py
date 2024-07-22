
import anndata as ad
import numpy as np
import scvi
import scanpy as sc
import matplotlib.pyplot as plt


scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

seurat2h5ad = ad.read_h5ad('../data/NN_PNAS_148822_632025.h5ad')


# print(seurat2h5ad.X.shape)
#
# fig = plt.figure(figsize=(10, 6))
# sc.pl.violin(seurat2h5ad, ['nFeature_RNA'],groupby = "orig.ident",multi_panel=True, jitter = False, rotation = 90)
# plt.savefig("../output/integration_nFeatureRNA_sample.png", bbox_inches='tight', pad_inches=0)
#
# fig = plt.figure(figsize=(10, 6))
# sc.pl.violin(seurat2h5ad, ['nCount_RNA'],groupby = "orig.ident",multi_panel=True, jitter = False, rotation = 90)
# plt.savefig("../output/integration_nCountRNA_sample.png", bbox_inches='tight', pad_inches=0)
#
# fig = plt.figure(figsize=(10, 6))
# sc.pl.violin(seurat2h5ad, ['percent.mt'],groupby = "orig.ident",multi_panel=True, jitter = False, rotation = 90)
# plt.savefig("../output/integration_percent_mt_sample.png", bbox_inches='tight', pad_inches=0)
#
# fig = plt.figure(figsize=(10, 6))
# sc.pl.violin(seurat2h5ad, ['percent.ribo'],groupby = "orig.ident",multi_panel=True, jitter = False, rotation = 90)
# plt.savefig("../output/integration_percent_ribo_sample.png", bbox_inches='tight', pad_inches=0)
#
seurat2h5ad = seurat2h5ad[seurat2h5ad.obs.nFeature_RNA < 4000, :]
seurat2h5ad = seurat2h5ad[seurat2h5ad.obs.nCount_RNA < 8000, :]

nuclei = list(seurat2h5ad.obs.index)
f_out = open('../output/nuclei.txt', "w")
for n in nuclei:
    f_out.write(n + "\n")
f_out.close()

# print(seurat2h5ad.X.shape)
#
# fig = plt.figure(figsize=(10, 6))
# sc.pl.violin(seurat2h5ad, ['nFeature_RNA'],groupby = "orig.ident",multi_panel=True, jitter = False, rotation = 90)
# plt.savefig("../output/integration_nFeatureRNA_sample_qc_4000.png", bbox_inches='tight', pad_inches=0)
#
# fig = plt.figure(figsize=(10, 6))
# sc.pl.violin(seurat2h5ad, ['nCount_RNA'],groupby = "orig.ident",multi_panel=True, jitter = False, rotation = 90)
# plt.savefig("../output/integration_nCountRNA_sample_qc_8000.png", bbox_inches='tight', pad_inches=0)
#


print("seurat2h5ad.X.shape = ", seurat2h5ad.X.shape)
seurat2h5ad.raw = seurat2h5ad  # freeze the state in `.raw`
seurat2h5ad.layers["counts"] = seurat2h5ad.X.copy()  # preserve counts
sc.pp.normalize_total(seurat2h5ad, target_sum=1e4)
sc.pp.log1p(seurat2h5ad)


sc.pp.highly_variable_genes(
    seurat2h5ad,
    n_top_genes=3000,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="orig.ident",
    span = 1.0
)


scvi.model.SCVI.setup_anndata(seurat2h5ad, batch_key = 'orig.ident', layer = "counts",
                              continuous_covariate_keys=['nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo'])


latent_dim = 30
max_epochs = 100

model = scvi.model.SCVI(seurat2h5ad, n_latent = latent_dim)

model.train(max_epochs = max_epochs, plan_kwargs={'lr':5e-4})

latent = model.get_latent_representation()

seurat2h5ad.obsm["X_scVI"] = latent

file_name = "../output/scVI_hvg_{}_correction_all_{}_latent_{}_lr_{}_epoch_{}_latent.txt".format(3000, True, latent_dim, 5e-4, max_epochs)

np.savetxt(file_name, seurat2h5ad.obsm["X_scVI"], delimiter ='\t')
