# 定义 PANDA 运行函数
run_PANDA <- function(st_counts_path, sc_counts_path, sc_labels_path, output_dir = ".",
                      n_archetypes_vec = 10, n_hvgs_sc = 2000, n_hvgs_st = 5000,
                      n_markers = 20, n_sample_cells = 500, sigma = 0.3,
                      do_parallel = TRUE, n_cores = 20) {

  # 加载必要的库
  library(PANDA)
  library(parallel)
  library(foreach)
  library(doParallel)
  library(peakRAM)
  library(Seurat)

  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # 读取数据
  sc_counts <- read.csv(sc_counts_path, row.names = 1)
  sc_labels <- read.csv(sc_labels_path, row.names = 1)$cell_type
  st_counts <- read.csv(st_counts_path, row.names = 1)

  # 记录开始时间
  start_time <- Sys.time()
  ##sc_results <- sc_train(sc_counts, sc_labels, n_archetypes_vec = 10, n_cores = 6, save_dir = "D:/Rshuju/shiyan/NPANDA_results/simulation/unpaired_scenario/comparison_methods/PANDA_results/sc_results")
  file_path <- "D:/Rshuju/shiyan/NPANDA_results/simulation/unpaired_scenario/comparison_methods/PANDA_results/sc_results/sc_results.rds"
  sc_results <- readRDS(file_path)


    st_results <- st_train(
      st_counts = st_counts,
      sc_results = sc_results,
      n_hvgs = n_hvgs_st,
      sigma = sigma,
      save_res = TRUE,
      save_mu_csv = TRUE,
      save_dir = output_dir
    )


  # 记录结束时间
  end_time <- Sys.time()

  # 计算运行时间
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))


  # 保存结果
  write.csv(as.matrix(st_results$proportion), file = paste(output_dir, "PANDA_prop.csv", sep = "/"))
  write.csv(as.matrix(st_results$alpha), file = paste(output_dir, "PANDA_alpha.csv", sep = "/"))
  write.csv(t(as.matrix(c("elapsed_time" = elapsed_time))), file = paste(output_dir, "PANDA_time_memory.csv", sep = "/"))
}


# #################### Run PANDA for paired scenario
for (i in list.files("D:/Rshuju/shiyan/NPANDA_data/paired_scenario/st_data")) {
  st_counts_path <- paste(paste("D:/Rshuju/shiyan/NPANDA_data/paired_scenario/st_data", i, sep = "/"), "st_counts.csv", sep = "/")
  output_dir <- paste("D:/Rshuju/shiyan/NPANDA_results/simulation/paired_scenario/comparison_methods/PANDA_results", i, sep = "/")
  run_PANDA(st_counts_path = st_counts_path,
            sc_counts_path = "D:/Rshuju/shiyan/NPANDA_data/paired_scenario/sc_data/counts_validation.csv",
            sc_labels_path = "D:/Rshuju/shiyan/NPANDA_data/paired_scenario/sc_data/labels_validation.csv",
            output_dir = output_dir,
            n_archetypes_vec = 10,
            n_hvgs_sc = 2000,
            n_hvgs_st = 5000,
            n_markers = 20,
            n_sample_cells = 500,
            sigma = 0.3,
            do_parallel = TRUE,
            n_cores = 20)
}


# #################### Run PANDA for unpaired scenario
for (i in list.files("D:/Rshuju/shiyan/NPANDA_data/unpaired_scenario/st_data")) {
  st_counts_path <- paste(paste("D:/Rshuju/shiyan/NPANDA_data/unpaired_scenario/st_data", i, sep = "/"), "st_counts.csv", sep = "/")
  output_dir <- paste("D:/Rshuju/shiyan/NPANDA_results/simulation/unpaired_scenario/comparison_methods/PANDA_results", i, sep = "/")
  run_PANDA(st_counts_path = st_counts_path,
            sc_counts_path = "D:/Rshuju/shiyan/NPANDA_data/unpaired_scenario/sc_data/pbmc2_inDrops_counts.csv",
            sc_labels_path = "D:/Rshuju/shiyan/NPANDA_data/unpaired_scenario/sc_data/pbmc2_inDrops_labels.csv",
            output_dir = output_dir,
            n_archetypes_vec = 10,
            n_hvgs_sc = 2000,
            n_hvgs_st = 5000,
            n_markers = 20,
            n_sample_cells = 500,
            sigma = 0.3,
            do_parallel = TRUE,
            n_cores = 20)
}


# #################### Run PANDA for merfish
st_counts_path <- "Data/processed_data/simulations/merfish/st_data/st_counts.csv"
output_dir <- "Results/simulations/merfish/comparison_methods/PANDA_results"
run_PANDA(st_counts_path = st_counts_path,
          sc_counts_path = "Data/processed_data/simulations/merfish/sc_data/sc_counts_subset.csv",
          sc_labels_path = "Data/processed_data/simulations/merfish/sc_data/sc_labels.csv",
          output_dir = output_dir,
          n_archetypes_vec = 10,
          n_hvgs_sc = 2000,
          n_hvgs_st = 5000,
          n_markers = 20,
          n_sample_cells = 500,
          sigma = 0.3,
          do_parallel = TRUE,
          n_cores = 20)


# #################### Run PANDA for time and memory
st_counts_path <- paste(paste("Data/processed_data/simulations/paired_scenario/st_data", "uniform_ST", sep = "/"), "st_counts.csv", sep = "/")
output_dir <- "Results/simulations/time_memory/comparison_methods/PANDA_results"
run_PANDA(st_counts_path = st_counts_path,
          sc_counts_path = "Data/processed_data/simulations/paired_scenario/sc_data/counts_validation.csv",
          sc_labels_path = "Data/processed_data/simulations/paired_scenario/sc_data/labels_validation.csv",
          output_dir = output_dir,
          n_archetypes_vec = 10,
          n_hvgs_sc = 2000,
          n_hvgs_st = 5000,
          n_markers = 20,
          n_sample_cells = 500,
          sigma = 0.3,
          do_parallel = TRUE,
          n_cores = 20)
