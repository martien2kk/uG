#%%
# Core libraries
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
adata

#%%
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

"""
Violin plots of some of the computed QC metrics:
- number of genes expressed in the count matrix
- total counts per cell
- percentage of counts in mitochondrial genes
"""

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

# scatter plot 
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)

# Doublet Score
sc.pp.scrublet(adata, batch_key="sample")

#%%
# Normalization

# Saving count data
adata.layers["counts"] = adata.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

# %%

# Feature Selection
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")

sc.pl.highly_variable_genes(adata)

# %%

# PCA / Dimentionality Reduction
sc.tl.pca(adata)

sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

# principal components
sc.pl.pca(
    adata,
    color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
)
# %%

# Nearest neighbor graph constuction and visualization
sc.pp.neighbors(adata)

sc.tl.umap(adata)

sc.pl.umap(
    adata,
    color="sample",
    # Setting a smaller point size to get prevent overlap
    size=2,
)
# %%

#Clustering

# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

sc.pl.umap(adata, color=["leiden"])

# Re-assess quality
sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3,
)

sc.pl.umap(
    adata,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
)

# %%

# Manual Cell Tyoe Annotation
# generate a set of clustering solutions which we can then use to annotate our cell types
for res in [0.02, 0.5, 2.0]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )

sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",
)

sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",
)

#%%

# Marker gene set
marker_genes = {
    "Microglia": ["AIF1", "C1QA", "CSF1R"],
    "Neuron": ["STMN2", "TUBB3"],
    "Astrocyte": ["GFAP", "AQP4", "ALDH1L1"],
    "oligodendrocyte lineage": ["OLIG1", "OLIG2", "PDGFRA"],
    "Myelinating Oligodendrocyte": ["PLP1", "MBP"],
}

sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.02", standard_scale="var")

adata.obs["cell_type_lvl1"] = adata.obs["leiden_res_0.02"].map(
    {
        "0": "Microglia",
        "1": "Neuron",
        "2": "Astrocyte",
        "3": "oligodendrocyte lineage",
        "4": "Myelinating Oligodendrocyte",
    }
)

sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.50", standard_scale="var")
# %%

# Differentially-expressed Genes as Markers

# Obtain cluster-specific differentially expressed genes
sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")

sc.pl.rank_genes_groups_dotplot(
    adata, groupby="leiden_res_0.50", standard_scale="var", n_genes=5
)

sc.get.rank_genes_groups_df(adata, group="7").head(5)

dc_cluster_genes = sc.get.rank_genes_groups_df(adata, group="7").head(5)["names"]
sc.pl.umap(
    adata,
    color=[*dc_cluster_genes, "leiden_res_0.50"],
    legend_loc="on data",
    frameon=False,
    ncols=3,
)

# %%
