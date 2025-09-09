import os
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import time
import tracemalloc

import cell2location
import scvi
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

def run_cell2location(st_counts_path, sc_counts_path, sc_labels_path, output_dir="."):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    sc_counts = pd.read_csv(sc_counts_path, header=0, index_col=0).astype("float32")
    sc_labels = np.array(pd.read_csv(sc_labels_path, header=0, index_col=0)["celltype"])

    use_labels = sc_labels

    st_counts = pd.read_csv(st_counts_path, header=0, index_col=0).astype("float32")

    # sc
    sc_start_time = time.time()
    tracemalloc.start()

    sc_adata = sc.AnnData(sc_counts)
    sc_adata.obs["cell_types"] = use_labels
    sc_adata.var['SYMBOL'] = sc_adata.var_names

    selected = filter_genes(sc_adata, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
    sc_adata = sc_adata[:, selected].copy()

    cell2location.models.RegressionModel.setup_anndata(adata=sc_adata, labels_key="cell_types")
    mod = RegressionModel(sc_adata)
    mod.train(max_epochs=250, accelerator="cpu")

    sc_adata = mod.export_posterior(sc_adata, sample_kwargs={'num_samples': 1000, 'batch_size': 2500})

    if 'means_per_cluster_mu_fg' in sc_adata.varm.keys():
        inf_aver = sc_adata.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in sc_adata.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = sc_adata.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in sc_adata.uns['mod']['factor_names']]].copy()
    inf_aver.columns = sc_adata.uns['mod']['factor_names']

    sc_size, sc_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    sc_end_time = time.time()
    sc_elapsed_time = sc_end_time - sc_start_time

    inf_aver.T.to_csv(output_dir + "/cell2location_signature.csv")


    # st
    st_start_time = time.time()
    tracemalloc.start()

    st_adata = sc.AnnData(st_counts)
    st_adata.var['SYMBOL'] = st_adata.var_names

    st_adata.var['MT_gene'] = [gene.startswith('MT-') for gene in st_adata.var['SYMBOL']]
    st_adata.obsm['MT'] = st_adata[:, st_adata.var['MT_gene'].values].X.toarray()
    st_adata = st_adata[:, ~st_adata.var['MT_gene'].values]


    intersect = np.intersect1d(st_adata.var_names, inf_aver.index)
    st_adata = st_adata[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    cell2location.models.Cell2location.setup_anndata(adata=st_adata)

    mod = cell2location.models.Cell2location(st_adata, cell_state_df=inf_aver, N_cells_per_location=30, detection_alpha=20)
    mod.train(max_epochs=30000, batch_size=None, train_size=1)

    st_adata = mod.export_posterior(st_adata, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs})
    weights = st_adata.obsm['q05_cell_abundance_w_sf']
    # weights.columns = [i.split("_")[-1] for i in list(weights.columns)]
    weights.columns = [i[len("q05cell_abundance_w_sf_"):] for i in list(weights.columns)]
    weights = (weights.T / weights.sum(axis=1)).T
    
    st_size, st_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    st_end_time = time.time()
    st_elapsed_time = st_end_time - st_start_time

    weights.to_csv(output_dir + "/cell2location_prop.csv")

    pd.DataFrame({"sc_elapsed_time": [sc_elapsed_time], "st_elapsed_time": [st_elapsed_time], "sc_peak": [sc_peak / 1024 / 1024], "st_peak": [st_peak / 1024 / 1024]}).to_csv(output_dir + "/cell2location_time_memory.csv")






if __name__ == "__main__":

    #################### Run cell2location for paired scenario
    for i in os.listdir("D:/Rshuju/shiyan/SMART_data/paired_scenario/st_data"):
        st_counts_path = "D:/Rshuju/shiyan/SMART_data/paired_scenario/st_data/" + i + "/st_counts.csv"
        output_dir = "D:/Rshuju/shiyan/SMART_results/simulation/paired_scenario/comparison_methods/cell2location_results/" + i
        run_cell2location(st_counts_path=st_counts_path,
                          sc_counts_path="D:/Rshuju/shiyan/SMART_data/paired_scenario/sc_data/counts_validation.csv",
                          sc_labels_path="D:/Rshuju/shiyan/SMART_data/paired_scenario/sc_data/labels_validation.csv",
                          output_dir=output_dir)


    #################### Run cell2location for unpaired scenario
    for i in os.listdir("D:/Rshuju/shiyan/SMART_data/unpaired_scenario/st_data"):
        st_counts_path = "D:/Rshuju/shiyan/SMART_data/unpaired_scenario/st_data/" + i + "/st_counts.csv"
        output_dir = "D:/Rshuju/shiyan/SMART_results/simulation/unpaired_scenario/comparison_methods/cell2location_results/" + i
        run_cell2location(st_counts_path=st_counts_path,
                          sc_counts_path="D:/Rshuju/shiyan/SMART_data/unpaired_scenario/sc_data/pbmc1_10x_v2_counts.csv",
                          sc_labels_path="D:/Rshuju/shiyan/SMART_data/unpaired_scenario/sc_data/pbmc1_10x_v2_labels.csv",
                          output_dir=output_dir)

