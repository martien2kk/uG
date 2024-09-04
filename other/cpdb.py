#%%
import pandas as pd
import sys
import os

pd.set_option('display.max_columns', 100)
current_directory = os.getcwd()
print("Current Working Directory:", current_directory)

#%%
os.chdir('/.mounts/labs/reimandlab/private/users/k2zhang/cpdb')

cpdb_file_path = '/db/test/v5.0.0/cellphonedb.zip'
meta_file_path = '../Ai_10xG_library/Ai_10xG_library.metrics_summary.csv'
counts_file_path = 'T../Ai_10xG_library/Ai_10xG_library.filtered_feature_bc_matrix.h5'
# microenvs_file_path = 'data/microenvironment.tsv'
out_path = 'results/method1'

metadata = pd.read_csv(meta_file_path, sep = ',')
metadata.head(3)

from cellphonedb.src.core.methods import cpdb_analysis_method

cpdb_results = cpdb_analysis_method.call(
    cpdb_file_path = cpdb_file_path,           # mandatory: CellphoneDB database zip file.
    meta_file_path = meta_file_path,           # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,       # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    counts_data = 'hgnc_symbol',               # defines the gene annotation in counts matrix.
    # microenvs_file_path = microenvs_file_path, # optional (default: None): defines cells per microenvironment.
    score_interactions = True,                 # optional: whether to score interactions or not. 
    output_path = out_path,                    # Path to save results    microenvs_file_path = None,
    separator = '|',                           # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    threads = 5,                               # number of threads to use in the analysis.
    threshold = 0.1,                           # defines the min % of cells expressing a gene for this to be employed in the analysis.
    result_precision = 3,                      # Sets the rounding for the mean values in significan_means.
    debug = False,                             # Saves all intermediate tables emplyed during the analysis in pkl format.
    output_suffix = None                       # Replaces the timestamp in the output files by a user defined string in the  (default: None)
)


# %%
