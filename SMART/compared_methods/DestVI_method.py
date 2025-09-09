import os
import scanpy as sc
import numpy as np
import pandas as pd
from scvi.model import CondSCVI, DestVI

from functools import reduce
import time
import tracemalloc



def run_DestVI(st_counts_path, sc_counts_path, sc_labels_path, output_dir="."):

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
    sc.pp.filter_genes(sc_adata, min_counts=10)
    sc_adata.layers["counts"] = sc_adata.X.copy()
    sc.pp.highly_variable_genes(
        sc_adata,
        n_top_genes=2000,
        subset=True,
        layer="counts",
        flavor="seurat_v3"
    )
    sc.pp.normalize_total(sc_adata, target_sum=10e4)
    sc.pp.log1p(sc_adata)
    sc_adata.raw = sc_adata

    st_adata = sc.AnnData(st_counts)
    st_adata.layers["counts"] = st_adata.X.copy()
    sc.pp.normalize_total(st_adata, target_sum=10e4)
    sc.pp.log1p(st_adata)
    st_adata.raw = st_adata

    intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
    st_adata = st_adata[:, intersect].copy()
    sc_adata = sc_adata[:, intersect].copy()

    CondSCVI.setup_anndata(sc_adata, layer="counts", labels_key="cell_types")
    sc_model = CondSCVI(sc_adata, weight_obs=False)
    sc_model.train()

    sc_size, sc_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    sc_end_time = time.time()
    sc_elapsed_time = sc_end_time - sc_start_time


    # st
    st_start_time = time.time()
    tracemalloc.start()

    DestVI.setup_anndata(st_adata, layer="counts")
    st_model = DestVI.from_rna_model(st_adata, sc_model)
    st_model.train(max_epochs=2500)

    weights = st_model.get_proportions()

    st_size, st_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    st_end_time = time.time()
    st_elapsed_time = st_end_time - st_start_time
    
    weights.to_csv(output_dir + "/DestVI_prop.csv")

    # mu_list = []
    signature_list = []
    unique_labels = np.unique(use_labels).tolist()
    px_o = st_model.module.px_o.detach().cpu().numpy()
    beta = st_model.module.beta.detach().cpu().numpy()
    for ct in unique_labels:
        mu_ct = st_model.get_scale_for_ct(ct)
        mu_ct.to_csv(output_dir + "/DestVI_" + ct + "_mu.csv")
        # mu_list.append(mu_ct)
        # signature_ct = mu_ct / (np.exp(px_o) * np.exp(beta))
        signature_ct = mu_ct / (np.exp(beta))
        signature_ct.to_csv(output_dir + "/DestVI_" + ct + "_signature.csv")
        signature_list.append(signature_ct)

    pd.DataFrame({"sc_elapsed_time": [sc_elapsed_time], "st_elapsed_time": [st_elapsed_time], "sc_peak": [sc_peak / 1024 / 1024], "st_peak": [st_peak / 1024 / 1024]}).to_csv(output_dir + "/DestVI_time_memory.csv")




if __name__ == "__main__":

    #################### Run DestVI for paired scenario
    for i in os.listdir("D:/Rshuju/shiyan/SMART_data/paired_scenario/st_data"):
        st_counts_path = "D:/Rshuju/shiyan/SMART_data/paired_scenario/st_data/" + i + "/st_counts.csv"
        output_dir = "D:/Rshuju/shiyan/SMART_results/simulation/paired_scenario/comparison_methods/DestVI_results/" + i
        run_DestVI(st_counts_path=st_counts_path,
                   sc_counts_path="D:/Rshuju/shiyan/SMART_data/paired_scenario/sc_data/counts_validation.csv",
                   sc_labels_path="D:/Rshuju/shiyan/SMART_data/paired_scenario/sc_data/labels_validation.csv",
                   output_dir=output_dir)


    #################### Run DestVI for unpaired scenario
    for i in os.listdir("D:/Rshuju/shiyan/SMART_data/unpaired_scenario/st_data"):
        st_counts_path = "D:/Rshuju/shiyan/SMART_data/unpaired_scenario/st_data/" + i + "/st_counts.csv"
        output_dir = "D:/Rshuju/shiyan/SMART_results/simulation/unpaired_scenario/comparison_methods/DestVI_results/" + i
        run_DestVI(st_counts_path=st_counts_path,
                   sc_counts_path="D:/Rshuju/shiyan/SMART_data/unpaired_scenario/sc_data/pbmc1_10x_v2_counts.csv",
                   sc_labels_path="D:/Rshuju/shiyan/SMART_data/unpaired_scenario/sc_data/pbmc1_10x_v2_labels.csv",
                   output_dir=output_dir)


