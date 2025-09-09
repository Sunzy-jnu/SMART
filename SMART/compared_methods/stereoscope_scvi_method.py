import os
import numpy as np
import pandas as pd
import scanpy as sc
import torch
import time
import tracemalloc

from scvi.external import RNAStereoscope, SpatialStereoscope



def run_stereoscope_scvi(st_counts_path, sc_counts_path, sc_labels_path, output_dir="."):

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

    non_mito_genes_list = [
        name for name in sc_adata.var_names if not name.startswith("MT-")
    ]
    sc_adata = sc_adata[:, non_mito_genes_list]

    sc_adata.layers["counts"] = sc_adata.X.copy()
    sc.pp.normalize_total(sc_adata, target_sum=1e5)
    sc.pp.log1p(sc_adata)
    sc_adata.raw = sc_adata

    sc.pp.highly_variable_genes(
        sc_adata,
        n_top_genes=7000,
        subset=True,
        layer="counts",
        flavor="seurat_v3",
        span=1,
    )

    st_adata = sc.AnnData(st_counts)
    st_adata.var_names_make_unique()

    intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
    st_adata = st_adata[:, intersect].copy()
    sc_adata = sc_adata[:, intersect].copy()

    RNAStereoscope.setup_anndata(sc_adata, layer="counts", labels_key="cell_types")

    sc_model = RNAStereoscope(sc_adata)
    sc_model.train(max_epochs=100)
    # sc_model.history["elbo_train"][10:].plot()

    sc_size, sc_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    sc_end_time = time.time()
    sc_elapsed_time = sc_end_time - sc_start_time

    sc_model.save(output_dir + "/scmodel", overwrite=True)

    # st
    st_start_time = time.time()
    tracemalloc.start()

    st_adata.layers["counts"] = st_adata.X.copy()
    SpatialStereoscope.setup_anndata(st_adata, layer="counts")

    spatial_model = SpatialStereoscope.from_rna_model(st_adata, sc_model)
    spatial_model.train(max_epochs=10000)
    # spatial_model.history["elbo_train"][10:].plot()

    st_size, st_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    st_end_time = time.time()
    st_elapsed_time = st_end_time - st_start_time

    spatial_model.save(output_dir + "/stmodel", overwrite=True)

    weights = spatial_model.get_proportions()

    weights.to_csv(output_dir + "/stereoscope_scvi_prop.csv")

    # signature = torch.nn.functional.softplus(spatial_model.module.W).detach().cpu().numpy().T
    signature = (torch.exp(spatial_model.module.px_o).unsqueeze(1) * torch.nn.functional.softplus(spatial_model.module.W)).detach().cpu().numpy().T
    signature = pd.DataFrame(signature, index=spatial_model.cell_type_mapping, columns=spatial_model.adata.var_names)
    signature.to_csv(output_dir + "/stereoscope_scvi_signature.csv")

    pd.DataFrame({"sc_elapsed_time": [sc_elapsed_time], "st_elapsed_time": [st_elapsed_time], "sc_peak": [sc_peak / 1024 / 1024], "st_peak": [st_peak / 1024 / 1024]}).to_csv(output_dir + "/stereoscope_scvi_time_memory.csv")





if __name__ == "__main__":

    #################### Run stereoscope_scvi for paired scenario
    for i in os.listdir("D:/Rshuju/shiyan/SMART_data/paired_scenario/st_data"):
        st_counts_path = "D:/Rshuju/shiyan/SMART_data/paired_scenario/st_data/" + i + "/st_counts.csv"
        output_dir = "D:/Rshuju/shiyan/SMART_results/simulation/paired_scenario/comparison_methods/stereoscope_scvi_results/" + i
        run_stereoscope_scvi(st_counts_path=st_counts_path,
                             sc_counts_path="D:/Rshuju/shiyan/SMART_data/paired_scenario/sc_data/counts_validation.csv",
                             sc_labels_path="D:/Rshuju/shiyan/SMART_data/paired_scenario/sc_data/labels_validation.csv",
                             output_dir=output_dir)


    #################### Run stereoscope_scvi for unpaired scenario
    for i in os.listdir("D:/Rshuju/shiyan/SMART_data/unpaired_scenario/st_data"):
        st_counts_path = "D:/Rshuju/shiyan/SMART_data/unpaired_scenario/st_data/" + i + "/st_counts.csv"
        output_dir = "D:/Rshuju/shiyan/SMART_results/simulation/unpaired_scenario/comparison_methods/stereoscope_scvi_results/" + i
        run_stereoscope_scvi(st_counts_path=st_counts_path,
                             sc_counts_path="D:/Rshuju/shiyan/SMART_data/unpaired_scenario/sc_data/pbmc1_10x_v2_counts.csv",
                             sc_labels_path="D:/Rshuju/shiyan/SMART_data/unpaired_scenario/sc_data/pbmc1_10x_v2_labels.csv",
                             output_dir=output_dir)


