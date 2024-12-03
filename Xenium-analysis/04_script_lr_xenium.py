import scanpy as sc
import squidpy as sq
import anndata as ad

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

from scipy.sparse import coo_matrix

from scipy.sparse import issparse
from scipy.stats import norm
import itertools

import sys
sys.path.insert(0, '/mnt/isilon/tan_lab/yuw1/R_work_dir/NB/Revision/xenium_venv/lib/python3.10/site-packages')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("mid", help="mouse id", type=str)
parser.add_argument("dist_thr", help="distance_thr", type=int)
args = parser.parse_args()
mid = args.mid
dist_thr = args.dist_thr

# ## Define functions

# In[ ]:


def get_connectivity_chart(adata, distance_threshold):
    """
    Function to generate connectivity chart between any two cells in spatial data.

    Parameters:
    adata : AnnData
      AnnData object containing single-cell data with spatial coordinates and cell type annotations.
      
    distance_threshold : float
      Distance threshold to define neighboring cells.

    Returns:
    connectivity_chart : pandas.DataFrame
      DataFrame containing pairs of cells of specified types that are within the distance threshold.
    """


    print("Building spatial graph.")
    # Build the spatial neighbors graph using the specified distance threshold
    sq.gr.spatial_neighbors(
      adata,
      radius=distance_threshold,
      coord_type="generic",
      spatial_key="spatial",  # Ensure your spatial coordinates are stored in adata.obsm['spatial']
    )

    print("Generating connectivity chart.")
    # Extract the adjacency matrix
    adjacency = adata.obsp["spatial_connectivities"].tocoo()
    edge_list = pd.DataFrame({"Cell_Number_A": adjacency.row, "Cell_Number_B": adjacency.col})

    # Get cell IDs
    all_cells = adata.obs_names.values
    edge_list["Cell_A"] = all_cells[edge_list["Cell_Number_A"]]
    edge_list["Cell_B"] = all_cells[edge_list["Cell_Number_B"]]

    # Get cell type identities
    edge_list["cell_state_A"] = adata.obs["cell_state"].values[edge_list["Cell_Number_A"]]
    edge_list["cell_state_B"] = adata.obs["cell_state"].values[edge_list["Cell_Number_B"]]

    return edge_list


def filter_connectivity_chart(edge_list, cell_state_A, cell_state_B):
    '''
    cell_state_A : str
      The identity of the first cell type.
    cell_state_B : str
      The identity of the second cell type.
    '''
    # Filter connections between cell_state_A and cell_state_B
    connectivity_chart = edge_list[
      (edge_list["cell_state_A"] == cell_state_A) &
      (edge_list["cell_state_B"] == cell_state_B)
    ].copy()

    return connectivity_chart


def calculate_ligand_receptor_score(adata, ligand_gene, receptor_gene, connectivity_chart):

    # Check if the ligand and receptor genes are in adata.var_names
    if ligand_gene not in adata.var_names:
        raise ValueError(f"Gene '{ligand_gene}' not found in the dataset.")
    if receptor_gene not in adata.var_names:
        raise ValueError(f"Gene '{receptor_gene}' not found in the dataset.")

    # Check if all cells in connectivity_chart are present in adata.obs_names
    missing_cells_A = set(connectivity_chart['Cell_A']) - set(adata.obs_names)
    if missing_cells_A:
        raise ValueError(f"Cells {missing_cells_A} in 'cell_state_A' not found in the dataset.")
    missing_cells_B = set(connectivity_chart['Cell_B']) - set(adata.obs_names)
    if missing_cells_B:
        raise ValueError(f"Cells {missing_cells_B} in 'cell_state_B' not found in the dataset.")

    # Helper function to extract gene expression values
    def get_gene_expression(adata, gene_name):
        data = adata[:, gene_name].X
        if issparse(data):
            data = data.toarray().flatten()
        else:
            data = np.array(data).flatten()
        return pd.Series(data, index=adata.obs_names)

    # Get expression values for the ligand and receptor genes
    ligand_expr = get_gene_expression(adata, ligand_gene)
    receptor_expr = get_gene_expression(adata, receptor_gene)

    # Extract expression values for the specified cells
    ligand_expr_A = ligand_expr.loc[connectivity_chart['Cell_A']].values
    receptor_expr_B = receptor_expr.loc[connectivity_chart['Cell_B']].values

    # Calculate the ligand-receptor score
    lr_score = np.sum(ligand_expr_A * receptor_expr_B)

    return lr_score



def get_non_connecting_chart(adata, cell_state_A, cell_state_B, connectivity_chart):

    # Get all cells of cell_state_A and cell_state_B
    all_cells_A = adata.obs[adata.obs['cell_state'] == cell_state_A].index.tolist()
    all_cells_B = adata.obs[adata.obs['cell_state'] == cell_state_B].index.tolist()

    print("Generating all possible non-connecting pairs.")

    # Ensure connectivity_chart has 'Cell_A' and 'Cell_B' columns
    if not {'Cell_A', 'Cell_B'}.issubset(connectivity_chart.columns):
        raise ValueError("connectivity_chart must have 'Cell_A' and 'Cell_B' columns.")

    # Get unique connecting pairs
    connecting_pairs = connectivity_chart[['Cell_A', 'Cell_B']].drop_duplicates()

    # Generate all possible pairs between all_cells_A and all_cells_B
    all_possible_pairs = pd.DataFrame(
      itertools.product(all_cells_A, all_cells_B), 
      columns=['Cell_A', 'Cell_B']
    )

    # Perform an anti-join to get non-connecting pairs
    merged = all_possible_pairs.merge(
      connecting_pairs, 
      on=['Cell_A', 'Cell_B'], 
      how='left', 
      indicator=True
    )
    non_connecting_pairs = merged[merged['_merge'] == 'left_only'].drop(columns=['_merge'])

    return non_connecting_pairs


def get_non_connecting_chart_large(adata, cell_state_A, cell_state_B, connectivity_chart):

    # Get all cells of cell_state_A and cell_state_B
    all_cells_A = adata.obs[adata.obs['cell_state'] == cell_state_A].index
    all_cells_B = adata.obs[adata.obs['cell_state'] == cell_state_B].index

    print("Generating non-connecting pairs for large datasets.")

    # Ensure connectivity_chart has 'Cell_A' and 'Cell_B' columns
    if not {'Cell_A', 'Cell_B'}.issubset(connectivity_chart.columns):
        raise ValueError("connectivity_chart must have 'Cell_A' and 'Cell_B' columns.")

    # Convert connectivity pairs to a set for faster lookup
    connecting_pairs_set = set(zip(connectivity_chart['Cell_A'], connectivity_chart['Cell_B']))

    # Initialize a list to collect non-connecting pairs
    non_connecting_pairs = []

    # Iterate over all possible pairs -- reduce if it's too large
    if len(all_cells_A) > 5000:
        all_cells_A =  np.random.choice(all_cells_A, size=5000, replace=False)

        
    if len(all_cells_B) > 5000:
        all_cells_B = np.random.choice(all_cells_B, size=5000, replace=False)
    
    
    for cell_A in all_cells_A:
        for cell_B in all_cells_B:
            if (cell_A, cell_B) not in connecting_pairs_set:
                non_connecting_pairs.append({'Cell_A': cell_A, 'Cell_B': cell_B})

    # Convert list to DataFrame
    non_connecting_pairs_df = pd.DataFrame(non_connecting_pairs)

    return non_connecting_pairs_df


## connectivity_chart columns: Cell_A, Cell_B 
def empirical_shuffling(
  adata,
  ligand_gene,
  receptor_gene,
  source,
  target,
  n_shuffles=1000,
  connectivity_chart=None,
  non_connecting_pairs=None
):
    from scipy.stats import norm
    from scipy.sparse import issparse

    print("Getting LR score for contacting pairs.")

    # Calculate actual LR score for connected pairs
    actual_lr_score = calculate_ligand_receptor_score(
      adata,
      ligand_gene,
      receptor_gene,
      connectivity_chart
    )
    
    if actual_lr_score == 0: ## no need to compute null distribution
        mean_shuffled_score = 0
        sd_shuffled_score = 0
        z_score = 0
        p_value = 0.5
        empr_p_value = 0.5
        fold_change = 0
        
    else:
        

        # Get normalized data
        # Assuming the expression data is stored in adata.X
        # and that adata has been properly normalized
        if issparse(adata.X):
            normalized_data = adata.X.toarray()
        else:
            normalized_data = adata.X

        normalized_data = pd.DataFrame(
          normalized_data, index=adata.obs_names, columns=adata.var_names
        )

        print("Calculating shuffled distribution.")

        shuffled_scores = np.zeros(n_shuffles)

        # For progress tracking
        print_interval = max(1, n_shuffles // 5)

        for i in range(n_shuffles):
            if (i + 1) % print_interval == 0:
                print(f"Shuffles completed: {i + 1}/{n_shuffles}")

            # Number of pairs in the actual connectivity chart
            num_pairs = len(connectivity_chart)

            # Sample pairs from non-connecting pairs
            sampled_pairs = non_connecting_pairs.sample(n=num_pairs, replace=False)

            # Get expression values for ligand and receptor genes
            ligand_expr_A_shuffled = normalized_data.loc[
              sampled_pairs['Cell_A'], ligand_gene
            ].values

            receptor_expr_B_shuffled = normalized_data.loc[
              sampled_pairs['Cell_B'], receptor_gene
            ].values

            # Calculate shuffled LR score
            shuffled_scores[i] = np.sum(ligand_expr_A_shuffled * receptor_expr_B_shuffled)

        print("Computing statistical significance.")

        # Calculate mean and standard deviation of shuffled scores
        mean_shuffled_score = np.mean(shuffled_scores)
        sd_shuffled_score = np.std(shuffled_scores, ddof=1)  # Sample standard deviation

        # Calculate z-score
        z_score = (actual_lr_score - mean_shuffled_score) / sd_shuffled_score

        # One-tailed p-value for enrichment
        p_value = norm.sf(z_score)  # Survival function (1 - CDF)

        empr_p_value = np.mean(shuffled_scores > actual_lr_score)
        # Calculate fold change
        fold_change = actual_lr_score / mean_shuffled_score if mean_shuffled_score != 0 else np.nan

    results = {
      'ligand': [ligand_gene],
      'receptor': [receptor_gene],
      'source': [source],
      'target': [target],  
      'actual_score': [actual_lr_score],
      'mean_shuffled': [mean_shuffled_score],
      'sd_shuffled': [sd_shuffled_score],
      'z_score': [z_score],
      'p_value': [p_value],
      'empr_p_value': [empr_p_value],
      'fold_change': [fold_change]
    }
    results = pd.DataFrame(results)

    return results




# ## Prepare Data

# In[ ]:


dir0 = 'data/intermediate/Xenium/'


# In[ ]:


adata = sc.read_h5ad('nbl_xenium_anndata_final.h5ad')

adata = adata[adata.obs['cell_type'].isin(['Neuroblast', 'Macrophage'])]

adata.obs.set_index('cell_id', inplace=True)


lr_res = pd.read_csv(dir0 + 'ligand_receptor_cytotalk.csv') ## read ligand receptor results


## read label transfer results generated by step3

nbl_state = pd.read_csv(dir0 + 'label_transfer_neuroblast_rpca.csv')
macro_state = pd.read_csv(dir0 + 'label_transfer_macrophage_rpca.csv')

nbl_state = nbl_state[nbl_state['prediction.score.max'] > 0.4]
macro_state = macro_state[macro_state['prediction.score.max'] > 0.4]

nbl_state['predicted.id'] = nbl_state['predicted.id'].astype('category')
nbl_state['predicted.id'].value_counts()
macro_state['predicted.id'] = macro_state['predicted.id'].astype('category')
macro_state['predicted.id'].value_counts()

nbl_state['predicted.id'] = nbl_state['predicted.id'].cat.reorder_categories(['ADRN-Calcium', 'ADRN-Baseline',
                                                                            'Interm-OxPhos', 'ADRN-Dopaminergic',
                                                                            'ADRN-Proliferating', 'MES'])

macro_state['predicted.id'] = macro_state['predicted.id'].cat.reorder_categories(['THY1+', 'HS3ST2+', 'F13A1+',
                                                   'CCL4+', 'IL18+', 'VCAN+', 'C1QC+SPP1+', 'Proliferating'])



adata_nbl = adata[adata.obs.index.isin(nbl_state['cell_id'])]
adata_macro = adata[adata.obs.index.isin(macro_state['cell_id'])]

# Ensure nbl_state is indexed by 'cell_id'
nbl_state.set_index('cell_id', inplace=True)
nbl_state.rename(columns = {'predicted.id':'cell_state'}, inplace = True)

# Join the 'cell_state' column from nbl_state to adata.obs
adata_nbl.obs = adata_nbl.obs.join(nbl_state['cell_state'])


# Ensure nbl_state is indexed by 'cell_id'
macro_state.set_index('cell_id', inplace=True)
macro_state.rename(columns = {'predicted.id':'cell_state'}, inplace = True)

# Join the 'cell_state' column from macro_state to adata.obs
adata_macro.obs = adata_macro.obs.join(macro_state['cell_state'])


# Perform the merge
adata = ad.concat([adata_nbl, adata_macro], axis=0)
del adata_nbl, adata_macro

# ### Distance = 50

# In[ ]:
adata = adata[adata.obs['mouse_id'] == mid]

edges = get_connectivity_chart(adata,  distance_threshold = dist_thr)


## run all L-R pairs
np.random.seed(42)
result_dfs = []
for index, row in lr_res.iterrows():
    print(f"Working on :", index, "th L-R pair: ")
    conn_chart = filter_connectivity_chart(edges, row['source'], row['target'])
    non_conn_chart = get_non_connecting_chart_large(adata, row['source'], row['target'], edges)
    print(f"# of paired cells btw", row['source'], "and",  row['target'], ':', len(conn_chart))
    print(f"# of non-paired cells btw", row['source'], "and",  row['target'], ':', len(non_conn_chart))
    result = empirical_shuffling(adata, row['ligand'], row['receptor'],row['source'], row['target'],
                                 n_shuffles = 100, connectivity_chart = conn_chart, 
                                 non_connecting_pairs = non_conn_chart)
    #result.to_csv(dir0 + 'ligand_receptor_dist' + str(dist_thr) +'_xenium_validation_' + mid + '.csv', 
    #              mode='a', index=False, header=False)
    result_dfs.append(result)


# Concatenate all the DataFrames into a single DataFrame
final_result_dfs = pd.concat(result_dfs, ignore_index=True)
final_result_dfs.to_csv(dir0 + 'LR_validation_dist' + str(dist_thr) + '_mouse' +  mid + '_rpca.csv', index=False)


