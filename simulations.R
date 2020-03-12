#################################################################################
#################################################################################

# Simulations 
# This script contains code to reproduce simulation results 
# in Figure 2a.b of the main manuscript.

#################################################################################
#################################################################################

## Hard scenario 
# Simlulate RNA matrix

set.seed(2019)
seed <- sample(10000, 10)
seed

load("SymSim/data/gene_len_pool.RData")

ari_list <- list()

for (i in 1:length(seed)) {
    
    
    cat("=========================", i, "===========================")
    
    ngenes <- 10000
    
    gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
    
    true_counts_res <- SimulateTrueCounts(ncells_total = 500, 
                                          min_popsize = 20, 
                                          i_minpop = 2, 
                                          ngenes = ngenes, 
                                          nevf = 30, 
                                          evf_type = "discrete", 
                                          n_de_evf=  18, 
                                          vary = "s", 
                                          Sigma = 0.9, 
                                          gene_effects_sd = 0.2,
                                          scale_s = 0.3,
                                          mean_hge = 6,
                                          phyla = phyla6, 
                                          randseed = seed[i])
    
    
    
    
    tsne_true_counts <- PlotTsne(meta=true_counts_res[[3]], data=log2(true_counts_res[[1]]+1), evf_type="continuous", n_pc=20, label='pop', saving = F, plotname="continuous populations (true counts)")
    
    ggplot(tsne_true_counts[[1]], aes(x = x, y = y, 
                                      color = as.factor(true_counts_res$cell_meta$pop))) +
        geom_point() +
        scale_color_tableau() +
        theme_bw() +
        theme(aspect.ratio = 1)
    
    
    
    
    observed_counts_RNA <- True2ObservedCounts(true_counts = true_counts_res[[1]], 
                                               meta_cell = true_counts_res[[3]], 
                                               protocol = "UMI", 
                                               alpha_mean = 0.0025, 
                                               alpha_sd = 0.0001, 
                                               gene_len = gene_len, 
                                               depth_mean = 45000, 
                                               depth_sd = 3000,
                                               nPCR1 = 14)
    
    hist(colSums(observed_counts_RNA$counts))
    
    tsne_RNA_counts <- PlotTsne(meta = observed_counts_RNA[[2]], 
                                data = log2(observed_counts_RNA[[1]]+1), 
                                evf_type = "discrete", 
                                n_pc = 20, label = 'pop', saving = F, 
                                plotname = "observed counts UMI")
    
    
    ggplot(tsne_RNA_counts[[1]], aes(x = x, y = y, color =  as.factor(true_counts_res$cell_meta$pop))) +
        geom_point() +
        scale_color_tableau() +
        theme_bw() +
        theme(aspect.ratio = 1)
    
    cat("Done RNA simulation", i)
    
# Simlulate ADT count matrix
    
    ngenes <- 100
    # load("SymSim/data/gene_len_pool.RData")
    gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
    
    true_counts_res_100 <- SimulateTrueCounts(ncells_total = 500, 
                                              min_popsize = 20, 
                                              i_minpop = 2, 
                                              ngenes = ngenes, 
                                              nevf = 30, 
                                              evf_type = "discrete", 
                                              n_de_evf=  20, 
                                              vary = "s", 
                                              Sigma = 0.5, 
                                              gene_effects_sd = 0.1,
                                              scale_s = 0.7,
                                              mean_hge = 6,
                                              phyla = phyla6, 
                                              randseed = seed[i])
    
    
    summary(colSums(true_counts_res_100[[1]]))
    
    tsne_true_counts_100 <- PlotTsne(meta=true_counts_res_100[[3]], data=log2(true_counts_res_100[[1]]+1), evf_type="continuous", n_pc=20, label='pop', saving = F, plotname="continuous populations (true counts)")
    
    
    ggplot(tsne_true_counts_100[[1]], aes(x = x, y = y, 
                                          color = as.factor(true_counts_res_100$cell_meta$pop))) +
        geom_point() +
        scale_color_tableau() +
        theme_bw() +
        theme(aspect.ratio = 1)
    
    gene_len <- sample(gene_len_pool[gene_len_pool > 10000], ngenes, replace = FALSE)
    
    
    observed_counts_ADT <- True2ObservedCounts(true_counts = true_counts_res_100[[1]], 
                                               meta_cell = true_counts_res_100[[3]], 
                                               protocol = "UMI", 
                                               alpha_mean = 0.1, 
                                               alpha_sd = 0.02, 
                                               gene_len = gene_len, 
                                               depth_mean = 10e+04, 
                                               depth_sd = 5000,
                                               nPCR1 = 10
    )
    
    tsne_ADT_counts <- PlotTsne(meta = observed_counts_ADT[[2]], 
                                data = log2(observed_counts_ADT[[1]]+1), 
                                evf_type = "discrete", 
                                n_pc = 20, label = 'pop', saving = F, 
                                plotname = "observed counts UMI")
    
    
    ggplot(tsne_ADT_counts[[1]], aes(x = x, y = y, color =  as.factor(true_counts_res_100$cell_meta$pop))) +
        geom_point() +
        scale_color_tableau() +
        theme_bw() +
        theme(aspect.ratio = 1)
    
    summary(colSums(observed_counts_ADT[[1]]))
    
    cat("Done ADT simulation", i)
    
# Perform different clustering algorithms for comparison
    
    # CiteFuse
    
    ## Create `SingleCellExperiment`
    library(CiteFuse)
    library(scater)
    library(scran)
    library(SingleCellExperiment)
    # sce_sim <- SingleCellExperiment(assay = list(counts = observed_counts_RNA$counts))
    # sce_sim <- scater::normalize(sce_sim)
    # rna_mat_sim_snf <- as.matrix(logcounts(sce_sim))
    # 
    # rho_adt_mat_sim <- propr(as.matrix(observed_counts_ADT$counts))
    # clr_adt_mat_sim <- rho_adt_mat_sim@logratio
    # 
    # dist_adt_sim <- 1 - as.matrix(rho_adt_mat_sim@matrix)
    # dist_rna_sim <- as.matrix(1 - cor(rna_mat_sim_snf))
    
    adt_mat <- observed_counts_ADT$counts
    rna_mat <- observed_counts_RNA$counts
    
    colnames(adt_mat) <- colnames(rna_mat) <- paste("Cell", 1:ncol(rna_mat), sep = "_")
    
    rownames(adt_mat) <- paste("ADT", 1:nrow(adt_mat), sep = "_")
    rownames(rna_mat) <- paste("RNA", 1:nrow(rna_mat), sep = "_")
    
    sce_sim <- preprocessing(list(RNA = rna_mat,
                                  ADT = adt_mat))
    
    sce_sim
    
    
    
    ## Normalisation
    
    
    sce_sim <- normaliseExprs(sce = sce_sim, 
                              altExp_name = "ADT", 
                              transform = "log")
    
    sce_sim <- scater::normalize(sce_sim)
    
    
    ## Run SNF
    
    system.time(sce_sim <- CiteFuse::CiteFuse(sce_sim,
                                              K_knn = 30))
    
    
    
    
    W_cluster_sim <- spectralClustering(metadata(sce_sim)[["SNF_W"]], K = 6)
    table(W_cluster_sim$labels, true_counts_res$cell_meta$pop)
    
    
    
    
    system.time(W1_cluster_sim <- spectralClustering(metadata(sce_sim)[["ADT_W"]], K = 6))
    system.time(W2_cluster_sim <- spectralClustering(metadata(sce_sim)[["RNA_W"]], K = 6))
    
    library(mclust)
    
    ari_W1 <- adjustedRandIndex(W1_cluster_sim$labels, true_counts_res$cell_meta$pop)
    ari_W1
    
    ari_W2 <- adjustedRandIndex(W2_cluster_sim$labels, true_counts_res$cell_meta$pop)
    ari_W2
    
    ari_W <- adjustedRandIndex(W_cluster_sim$labels, true_counts_res$cell_meta$pop)
    ari_W
    
    W_louvain_sim <- CiteFuse::louvain(sce_sim, metadata = "SNF_W")
    W1_louvain_sim <- CiteFuse::louvain(sce_sim, metadata = "ADT_W")
    W2_louvain_sim <- CiteFuse::louvain(sce_sim, metadata = "RNA_W")
    
    ari_W_louvain <- adjustedRandIndex(W_louvain_sim, true_counts_res$cell_meta$pop)
    ari_W_louvain
    
    ari_W1_louvain <- adjustedRandIndex(W1_louvain_sim, true_counts_res$cell_meta$pop)
    ari_W1_louvain
    
    ari_W2_louvain <- adjustedRandIndex(W2_louvain_sim, true_counts_res$cell_meta$pop)
    ari_W2_louvain
    
    
    
    cat("\n")
    cat("Done SNF, ", "ARI: ", ari_W)
    
    ### Visualisation
    
    
    
    
    ggplot(tsne_RNA_counts[[1]], aes(x = x, y = y, color =  as.factor(W_cluster_sim$labels))) +
        geom_point() +
        scale_color_tableau() +
        theme_bw() +
        theme(aspect.ratio = 1)
    
    ggplot(tsne_ADT_counts[[1]], aes(x = x, y = y, color =  as.factor(W_cluster_sim$labels))) +
        geom_point() +
        scale_color_tableau() +
        theme_bw() +
        theme(aspect.ratio = 1)
    
    
    
    
    #### Seurat
    
    
    library(Seurat)
    # colnames(observed_counts_RNA$counts) <- observed_counts_RNA$cell_meta$cellid
    # 
    # rownames(observed_counts_RNA$counts) <- paste("gene", 1:nrow(sce_sim), sep = "_")
    
    
    sim <- CreateSeuratObject(counts = rna_mat,
                              project = "sim",
                              min.cells = 0,
                              min.features = 0)
    sim <- NormalizeData(sim)
    
    sim <- FindVariableFeatures(sim, selection.method = "vst", nfeatures = 2000)
    
    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(sim), 10)
    
    # plot variable features with and without labels
    # plot1 <- VariableFeaturePlot(sim)
    # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    # CombinePlots(plots = list(plot1, plot2))
    
    all.genes <- rownames(sim)
    sim <- ScaleData(sim, features = all.genes)
    
    sim <- RunPCA(sim, features = VariableFeatures(object = sim))
    
    sim <- FindNeighbors(sim, dims = 1:10)
    sim <- FindClusters(sim, resolution = 1.8)
    
    table(sim$seurat_clusters, true_counts_res$cell_meta$pop)
    
    
    adjustedRandIndex(sim$seurat_clusters, true_counts_res$cell_meta$pop)
    
    DimPlot(sim, reduction = "pca")
    
    
    ari_seurat_rna <- adjustedRandIndex(sim$seurat_clusters, true_counts_res$cell_meta$pop)
    ari_seurat_rna
    
    
    
    
    
    library(Seurat)
    # colnames(observed_counts_RNA$counts) <- observed_counts_RNA$cell_meta$cellid
    # 
    # rownames(observed_counts_RNA$counts) <- paste("gene", 1:nrow(sce_sim), sep = "_")
    
    
    sim_adt <- CreateSeuratObject(counts = adt_mat,
                                  project = "sim_adt",
                                  min.cells = 0,
                                  min.features = 0)
    sim_adt <- NormalizeData(sim_adt)
    
    sim_adt <- FindVariableFeatures(sim_adt, selection.method = "vst", nfeatures = 2000)
    
    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(sim_adt), 10)
    
    # plot variable features with and without labels
    # plot1 <- VariableFeaturePlot(sim_adt)
    # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    # CombinePlots(plots = list(plot1, plot2))
    
    all.genes <- rownames(sim_adt)
    sim_adt <- ScaleData(sim_adt, features = all.genes)
    
    sim_adt <- RunPCA(sim_adt, features = VariableFeatures(object = sim_adt))
    
    sim_adt <- FindNeighbors(sim_adt, dims = 1:10)
    sim_adt <- FindClusters(sim_adt, resolution = 2.1)
    
    table(sim_adt$seurat_clusters, true_counts_res$cell_meta$pop)
    
    
    adjustedRandIndex(sim_adt$seurat_clusters, true_counts_res$cell_meta$pop)
    
    DimPlot(sim_adt, reduction = "pca")
    
    
    ari_seurat_adt <- adjustedRandIndex(sim_adt$seurat_clusters, true_counts_res$cell_meta$pop)
    ari_seurat_adt
    
    cat("Done Seurat, ", "ARI (ARI): ", ari_seurat_rna, "ARI (ADT): ", ari_seurat_adt)
    
    ### Benchmark with other
    
    
    #### SIMLR
    
    
    library(SIMLR)
    simlr_res <- SIMLR(X = as.matrix(logcounts(sce_sim)), c = 6, cores.ratio = 0)
    
    adjustedRandIndex(simlr_res$y$cluster, true_counts_res$cell_meta$pop)
    
    
    
    ari_simlr_rna <- adjustedRandIndex(simlr_res$y$cluster, true_counts_res$cell_meta$pop)
    
    
    simlr_res_adt <- SIMLR(X = as.matrix(assay(altExp(sce_sim, "ADT"), "logcounts")), c = 6, cores.ratio = 0)
    
    adjustedRandIndex(simlr_res_adt$y$cluster, true_counts_res$cell_meta$pop)
    
    
    
    ari_simlr_adt <- adjustedRandIndex(simlr_res_adt$y$cluster, true_counts_res$cell_meta$pop)
    
    
    
    #### pca + kmeans
    
    
    prcomp_res <- prcomp(t(logcounts(sce_sim)))
    km_res_rna <- kmeans(prcomp_res$x[,1:10], centers = 6)
    adjustedRandIndex(km_res_rna$cluster, true_counts_res$cell_meta$pop)
    
    ari_pca_rna <- adjustedRandIndex(km_res_rna$cluster, true_counts_res$cell_meta$pop)
    
    
    prcomp_res <- prcomp(t(as.matrix(assay(altExp(sce_sim, "ADT"), "logcounts"))))
    km_res_adt <- kmeans(prcomp_res$x[,1:10], centers = 6)
    adjustedRandIndex(km_res_adt$cluster, true_counts_res$cell_meta$pop)
    
    ari_pca_adt <- adjustedRandIndex(km_res_adt$cluster, true_counts_res$cell_meta$pop)
    
    
    
    
    
    df <- data.frame(ARI = c(ari_W, ari_W1, ari_W2, 
                             ari_W_louvain, ari_W1_louvain, ari_W2_louvain,
                             ari_simlr_rna, ari_simlr_adt, 
                             ari_pca_rna, ari_pca_adt,
                             ari_seurat_rna, ari_seurat_adt),
                     method = c("SNF (Fused)", "SNF (ADT)", "SNF (RNA)",
                                "SNF (Louvain)", "SNF (Louvain-ADT)", "SNF (Louvain-RNA)",
                                "SIMLR (RNA)", "SIMLR (ADT)",
                                "PCA (RNA)", "PCA (ADT)",
                                "Seurat (RNA)", "Seurat (ADT)"))
    
    df$method <- as.character(df$method)
    df$method <- factor(df$method, 
                        levels = c("SNF (Fused)", "SNF (ADT)", "SNF (RNA)",
                                   "SNF (Louvain)", "SNF (Louvain-ADT)", "SNF (Louvain-RNA)",
                                   "SIMLR (RNA)", "SIMLR (ADT)",
                                   "PCA (RNA)", "PCA (ADT)",
                                   "Seurat (RNA)", "Seurat (ADT)"))
    
    method_col <- tableau_color_pal("Tableau 20")(20)
    method_col <- method_col[c(11, 1:8)]
    
    library(moon)
    
    ggplot(df, aes(x = method, y = ARI, fill = method)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = method_col) +
        theme_bw() +
        theme(text = element_text(size = 12),
              axis.text.x = element_text(angle = 90))
    
    df$simulation <- i
    
    
    ari_list[[i]] <- df
    
}
ari_list <- lapply(ari_list, function(x) {
    x$method <- c("SNF (Fused)", "SNF (ADT)", "SNF (RNA)",
                  "SNF (Louvain)", "SNF (Louvain-ADT)", "SNF (Louvain-RNA)",
                  "SIMLR (RNA)", "SIMLR (ADT)",
                  "PCA (RNA)", "PCA (ADT)",
                  "Seurat (RNA)", "Seurat (ADT)")
    x$method <- factor(x$method, levels = x$method)
    x
})

ari_list_sim2 <- do.call(rbind, ari_list)

ggplot(ari_list_sim2, aes(x = method, y = ARI, 
                          color = method)) +
    # geom_violin() +
    geom_boxplot() +
    #scale_fill_manual(values = method_col) +
    #scale_color_manual(values = method_col) +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 90))

ggsaveWithDate("simulation_n500_c6_boxplot_sim2", width = 8, height = 6)



## Easy scenario 
# Simlulate RNA matrix

ari_list <- list()

for (i in 1:length(seed)) {
    
    
    cat("=========================", i, "===========================")
    
    ngenes <- 10000
    
    gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
    
    true_counts_res <- SimulateTrueCounts(ncells_total = 500, 
                                          min_popsize = 50, 
                                          i_minpop = 2, 
                                          ngenes = ngenes, 
                                          nevf = 30, 
                                          evf_type = "discrete", 
                                          n_de_evf=  18, 
                                          vary = "s", 
                                          Sigma = 0.8, 
                                          gene_effects_sd = 0.2,
                                          scale_s = 0.3,
                                          mean_hge = 6,
                                          phyla = phyla6, 
                                          randseed = seed[i])
    
    
    
    
    tsne_true_counts <- PlotTsne(meta=true_counts_res[[3]], data=log2(true_counts_res[[1]]+1), evf_type="continuous", n_pc=20, label='pop', saving = F, plotname="continuous populations (true counts)")
    
    ggplot(tsne_true_counts[[1]], aes(x = x, y = y, 
                                      color = as.factor(true_counts_res$cell_meta$pop))) +
        geom_point() +
        scale_color_tableau() +
        theme_bw() +
        theme(aspect.ratio = 1)
    
    
    cat("Done RNA")
    
    
    
    observed_counts_RNA <- True2ObservedCounts(true_counts = true_counts_res[[1]], 
                                               meta_cell = true_counts_res[[3]], 
                                               protocol = "UMI", 
                                               alpha_mean = 0.0025, 
                                               alpha_sd = 0.0001, 
                                               gene_len = gene_len, 
                                               depth_mean = 45000, 
                                               depth_sd = 3000,
                                               nPCR1 = 14)
    
    
    hist(colSums(observed_counts_RNA$counts))
    
    tsne_RNA_counts <- PlotTsne(meta = observed_counts_RNA[[2]], 
                                data = log2(observed_counts_RNA[[1]]+1), 
                                evf_type = "discrete", 
                                n_pc = 20, label = 'pop', saving = F, 
                                plotname = "observed counts UMI")
    
    
    ggplot(tsne_RNA_counts[[1]], aes(x = x, y = y, color =  as.factor(true_counts_res$cell_meta$pop))) +
        geom_point() +
        scale_color_tableau() +
        theme_bw() +
        theme(aspect.ratio = 1)
    
# Simlulate ADT count matrix
    
    ngenes <- 100
    # load("SymSim/data/gene_len_pool.RData")
    gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
    
    true_counts_res_100 <- SimulateTrueCounts(ncells_total = 500, 
                                              min_popsize = 50, 
                                              i_minpop = 2, 
                                              ngenes = ngenes, 
                                              nevf = 30, 
                                              evf_type = "discrete", 
                                              n_de_evf=  20, 
                                              vary = "s", 
                                              Sigma = 0.2, 
                                              gene_effects_sd = 0.1,
                                              scale_s = 0.7,
                                              mean_hge = 6,
                                              phyla = phyla6, 
                                              randseed = seed[i])
    
    
    summary(colSums(true_counts_res_100[[1]]))
    
    tsne_true_counts_100 <- PlotTsne(meta=true_counts_res_100[[3]], data=log2(true_counts_res_100[[1]]+1), evf_type="continuous", n_pc=20, label='pop', saving = F, plotname="continuous populations (true counts)")
    
    
    ggplot(tsne_true_counts_100[[1]], aes(x = x, y = y, 
                                          color = as.factor(true_counts_res_100$cell_meta$pop))) +
        geom_point() +
        scale_color_tableau() +
        theme_bw() +
        theme(aspect.ratio = 1)
    
    gene_len <- sample(gene_len_pool[gene_len_pool > 10000], ngenes, replace = FALSE)
    
    
    observed_counts_ADT <- True2ObservedCounts(true_counts = true_counts_res_100[[1]], 
                                               meta_cell = true_counts_res_100[[3]], 
                                               protocol = "UMI", 
                                               alpha_mean = 0.1, 
                                               alpha_sd = 0.02, 
                                               gene_len = gene_len, 
                                               depth_mean = 5e+04, 
                                               depth_sd = 5000,
                                               nPCR1 = 10
    )
    
    cat("Done ADT simulation", i)
    
# Perform different clustering algorithms for comparison
    
    # CiteFuse
    
    ## Create `SingleCellExperiment`
    
    library(scater)
    library(scran)
    library(SingleCellExperiment)
    # sce_sim <- SingleCellExperiment(assay = list(counts = observed_counts_RNA$counts))
    # sce_sim <- scater::normalize(sce_sim)
    # rna_mat_sim_snf <- as.matrix(logcounts(sce_sim))
    # 
    # rho_adt_mat_sim <- propr(as.matrix(observed_counts_ADT$counts))
    # clr_adt_mat_sim <- rho_adt_mat_sim@logratio
    # 
    # dist_adt_sim <- 1 - as.matrix(rho_adt_mat_sim@matrix)
    # dist_rna_sim <- as.matrix(1 - cor(rna_mat_sim_snf))
    
    adt_mat <- observed_counts_ADT$counts
    rna_mat <- observed_counts_RNA$counts
    
    colnames(adt_mat) <- colnames(rna_mat) <- paste("Cell", 1:ncol(rna_mat), sep = "_")
    
    rownames(adt_mat) <- paste("ADT", 1:nrow(adt_mat), sep = "_")
    rownames(rna_mat) <- paste("RNA", 1:nrow(rna_mat), sep = "_")
    
    sce_sim <- preprocessing(list(RNA = rna_mat,
                                  ADT = adt_mat))
    
    sce_sim
    
    
    
    ## Normalisation
    
    
    sce_sim <- normaliseExprs(sce = sce_sim, 
                              altExp_name = "ADT", 
                              transform = "log")
    
    sce_sim <- scater::normalize(sce_sim)
    
    
    ## Run SNF
    
    
    system.time(sce_sim <- CiteFuse::CiteFuse(sce_sim,
                                              K_knn = 30))
    
    
    
    
    W_cluster_sim <- spectralClustering(metadata(sce_sim)[["SNF_W"]], K = 6)
    table(W_cluster_sim$labels, true_counts_res$cell_meta$pop)
    
    
    
    
    system.time(W1_cluster_sim <- spectralClustering(metadata(sce_sim)[["ADT_W"]], K = 6))
    system.time(W2_cluster_sim <- spectralClustering(metadata(sce_sim)[["RNA_W"]], K = 6))
    
    library(mclust)
    
    ari_W1 <- adjustedRandIndex(W1_cluster_sim$labels, true_counts_res$cell_meta$pop)
    ari_W1
    
    ari_W2 <- adjustedRandIndex(W2_cluster_sim$labels, true_counts_res$cell_meta$pop)
    ari_W2
    
    ari_W <- adjustedRandIndex(W_cluster_sim$labels, true_counts_res$cell_meta$pop)
    ari_W
    
    W_louvain_sim <- CiteFuse::louvain(sce_sim, metadata = "SNF_W")
    W1_louvain_sim <- CiteFuse::louvain(sce_sim, metadata = "ADT_W")
    W2_louvain_sim <- CiteFuse::louvain(sce_sim, metadata = "RNA_W")
    
    ari_W_louvain <- adjustedRandIndex(W_louvain_sim, true_counts_res$cell_meta$pop)
    ari_W_louvain
    
    ari_W1_louvain <- adjustedRandIndex(W1_louvain_sim, true_counts_res$cell_meta$pop)
    ari_W1_louvain
    
    ari_W2_louvain <- adjustedRandIndex(W2_louvain_sim, true_counts_res$cell_meta$pop)
    ari_W2_louvain
    
    
    cat("\n")
    cat("Done SNF, ", "ARI: ", ari_W)
    
    ### Visualisation
    
    
    
    
    ggplot(tsne_RNA_counts[[1]], aes(x = x, y = y, color =  as.factor(W_cluster_sim$labels))) +
        geom_point() +
        scale_color_tableau() +
        theme_bw() +
        theme(aspect.ratio = 1)
    
    ggplot(tsne_ADT_counts[[1]], aes(x = x, y = y, color =  as.factor(W_cluster_sim$labels))) +
        geom_point() +
        scale_color_tableau() +
        theme_bw() +
        theme(aspect.ratio = 1)
    
    
    
    
    #### Seurat
    
    
    library(Seurat)
    # colnames(observed_counts_RNA$counts) <- observed_counts_RNA$cell_meta$cellid
    # 
    # rownames(observed_counts_RNA$counts) <- paste("gene", 1:nrow(sce_sim), sep = "_")
    
    
    sim <- CreateSeuratObject(counts = rna_mat,
                              project = "sim",
                              min.cells = 0,
                              min.features = 0)
    sim <- NormalizeData(sim)
    
    sim <- FindVariableFeatures(sim, selection.method = "vst", nfeatures = 2000)
    
    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(sim), 10)
    
    # plot variable features with and without labels
    # plot1 <- VariableFeaturePlot(sim)
    # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    # CombinePlots(plots = list(plot1, plot2))
    
    all.genes <- rownames(sim)
    sim <- ScaleData(sim, features = all.genes)
    
    sim <- RunPCA(sim, features = VariableFeatures(object = sim))
    
    sim <- FindNeighbors(sim, dims = 1:10)
    sim <- FindClusters(sim, resolution = 1.8)
    
    table(sim$seurat_clusters, true_counts_res$cell_meta$pop)
    
    
    adjustedRandIndex(sim$seurat_clusters, true_counts_res$cell_meta$pop)
    
    DimPlot(sim, reduction = "pca")
    
    
    ari_seurat_rna <- adjustedRandIndex(sim$seurat_clusters, true_counts_res$cell_meta$pop)
    ari_seurat_rna
    
    
    
    
    
    library(Seurat)
    # colnames(observed_counts_RNA$counts) <- observed_counts_RNA$cell_meta$cellid
    # 
    # rownames(observed_counts_RNA$counts) <- paste("gene", 1:nrow(sce_sim), sep = "_")
    
    
    sim_adt <- CreateSeuratObject(counts = adt_mat,
                                  project = "sim_adt",
                                  min.cells = 0,
                                  min.features = 0)
    sim_adt <- NormalizeData(sim_adt)
    
    sim_adt <- FindVariableFeatures(sim_adt, selection.method = "vst", nfeatures = 2000)
    
    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(sim_adt), 10)
    
    # plot variable features with and without labels
    # plot1 <- VariableFeaturePlot(sim_adt)
    # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    # CombinePlots(plots = list(plot1, plot2))
    
    all.genes <- rownames(sim_adt)
    sim_adt <- ScaleData(sim_adt, features = all.genes)
    
    sim_adt <- RunPCA(sim_adt, features = VariableFeatures(object = sim_adt))
    
    sim_adt <- FindNeighbors(sim_adt, dims = 1:10)
    sim_adt <- FindClusters(sim_adt, resolution = 2.1)
    
    table(sim_adt$seurat_clusters, true_counts_res$cell_meta$pop)
    
    
    adjustedRandIndex(sim_adt$seurat_clusters, true_counts_res$cell_meta$pop)
    
    DimPlot(sim_adt, reduction = "pca")
    
    
    ari_seurat_adt <- adjustedRandIndex(sim_adt$seurat_clusters, true_counts_res$cell_meta$pop)
    ari_seurat_adt
    
    cat("Done Seurat, ", "ARI (ARI): ", ari_seurat_rna, "ARI (ADT): ", ari_seurat_adt)
    
    ### Benchmark with other
    
    
    #### SIMLR
    
    
    library(SIMLR)
    simlr_res <- SIMLR(X = as.matrix(logcounts(sce_sim)), c = 6, cores.ratio = 0)
    
    adjustedRandIndex(simlr_res$y$cluster, true_counts_res$cell_meta$pop)
    
    
    
    ari_simlr_rna <- adjustedRandIndex(simlr_res$y$cluster, true_counts_res$cell_meta$pop)
    
    
    simlr_res_adt <- SIMLR(X = as.matrix(assay(altExp(sce_sim, "ADT"), "logcounts")), c = 6, cores.ratio = 0)
    
    adjustedRandIndex(simlr_res_adt$y$cluster, true_counts_res$cell_meta$pop)
    
    
    
    ari_simlr_adt <- adjustedRandIndex(simlr_res_adt$y$cluster, true_counts_res$cell_meta$pop)
    
    
    
    #### pca + kmeans
    
    
    prcomp_res <- prcomp(t(logcounts(sce_sim)))
    km_res_rna <- kmeans(prcomp_res$x[,1:10], centers = 6)
    adjustedRandIndex(km_res_rna$cluster, true_counts_res$cell_meta$pop)
    
    ari_pca_rna <- adjustedRandIndex(km_res_rna$cluster, true_counts_res$cell_meta$pop)
    
    
    prcomp_res <- prcomp(t(as.matrix(assay(altExp(sce_sim, "ADT"), "logcounts"))))
    km_res_adt <- kmeans(prcomp_res$x[,1:10], centers = 6)
    adjustedRandIndex(km_res_adt$cluster, true_counts_res$cell_meta$pop)
    
    ari_pca_adt <- adjustedRandIndex(km_res_adt$cluster, true_counts_res$cell_meta$pop)
    
    
    
    
    
    df <- data.frame(ARI = c(ari_W, ari_W1, ari_W2, 
                             ari_W_louvain, ari_W1_louvain, ari_W2_louvain,
                             ari_simlr_rna, ari_simlr_rna, 
                             ari_pca_rna, ari_pca_adt,
                             ari_seurat_rna, ari_seurat_adt),
                     method = c("SNF (Fused)", "SNF (ADT)", "SNF (RNA)",
                                "SNF (Louvain)", "SNF (Louvain-ADT)", "SNF (Louvain-RNA)",
                                "SIMLR (RNA)", "SIMLR (ADT)",
                                "PCA (RNA)", "PCA (ADT)",
                                "Seurat (RNA)", "Seurat (ADT)"))
    
    df$method <- as.character(df$method)
    df$method <- factor(df$method, 
                        levels = c("SNF (Fused)", "SNF (ADT)", "SNF (RNA)",
                                   "SNF (Louvain)", "SNF (Louvain-ADT)", "SNF (Louvain-RNA)",
                                   "SIMLR (RNA)", "SIMLR (ADT)",
                                   "PCA (RNA)", "PCA (ADT)",
                                   "Seurat (RNA)", "Seurat (ADT)"))
    
    method_col <- tableau_color_pal("Tableau 20")(20)
    method_col <- method_col[c(11, 1:8)]
    
    library(moon)
    
    ggplot(df, aes(x = method, y = ARI, fill = method)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = method_col) +
        theme_bw() +
        theme(text = element_text(size = 12),
              axis.text.x = element_text(angle = 90))
    
    df$simulation <- i
    
    
    ari_list[[i]] <- df
    
}


ari_list_sim1 <- do.call(rbind, ari_list)

ggplot(ari_list_sim1, aes(x = method, y = ARI, 
                          color = method)) +
    # geom_violin() +
    geom_boxplot() +
    # scale_fill_manual(values = method_col) +
    # scale_color_manual(values = method_col) +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 90))

ggsaveWithDate("simulation_n500_c6_boxplot_sim1", width = 8, height = 6)
