#%%
# Core scverse libraries
import os
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

sc.settings.set_figure_params(dpi=50, facecolor="white")


import pooch

# adata = sc.read_10x_h5('/.mounts/labs/reimandlab/private/users/k2zhang/TIA18674.20211202/Ai_10xG_library/Ai_10xG_library.filtered_feature_bc_matrix.h5')

# Define the path to your data directory
data_dir = '/.mounts/labs/reimandlab/private/users/k2zhang/TIA18674.20211202/Ai_10xG_library/'

# Define the sample files
samples = {
    "sample1": "Ai_10xG_library.filtered_feature_bc_matrix.h5",
    # Add more samples as needed
}

# Initialize a dictionary to store AnnData objects
adatas = {}

# Load each sample into an AnnData object and store in the dictionary
for sample_id, filename in samples.items():
    path = os.path.join(data_dir, filename)
    sample_adata = sc.read_10x_h5(path)
    sample_adata.var_names_make_unique()
    adatas[sample_id] = sample_adata

# Concatenate all the AnnData objects
adata = ad.concat(adatas, label="sample")

# Ensure observation names are unique
adata.obs_names_make_unique()

# Print the count of observations per sample
print(adata.obs["sample"].value_counts())


#%%

# Basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter out cells with high mitochondrial content and low gene count
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# Regress out unwanted sources of variation
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

#%%
# Perform PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

# Plot PCA
sc.pl.pca(adata)


# Show the plot
plt.show()

#%%
# Compute the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
# do pca varaiance ratio plot and then decide cluster number




# Embed the neighborhood graph using UMAP
sc.tl.umap(adata)

# Cluster the cells
sc.tl.leiden(adata)

# Plot the UMAP
sc.pl.umap(adata, color=['leiden', 'pct_counts_mt', 'n_genes_by_counts'])

# Save the processed AnnData object
adata.write('processed_adata.h5ad')



# %%
# Plot the UMAP and save to file
sc.pl.umap(adata, color=['leiden', 'pct_counts_mt', 'n_genes_by_counts'], save='_umap_plot.png')

# %%
# Viewing the first few cells (observations)
print("First few cells (observations):")
print(adata.obs.head())

# Viewing the first few genes (features)
print("\nFirst few genes (features):")
print(adata.var.head())

# %%
