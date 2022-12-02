## ----libraries----------------------------------------------------------------
source("../bin/chunk-options.R")
knitr_fig_path("05-")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))

data_dir <- '../data'


## ----seed---------------------------------------------------------------------
# set a seed for reproducibility in case any randomness used below
set.seed(1418)


## ----load_data----------------------------------------------------------------
liver <- readRDS(file.path(data_dir, 'lesson04.rds'))


## ----normalization, message=FALSE---------------------------------------------
liver <- liver %>%
            NormalizeData(normalization.method = "LogNormalize")


## ----var_features, message=FALSE, warning=FALSE-------------------------------
liver <- liver %>% 
              FindVariableFeatures(nfeatures = 2000)

# Identify the 25 most highly variable genes
top25 <- head(VariableFeatures(liver), 25)

plot1 <- VariableFeaturePlot(liver)
LabelPoints(plot = plot1, points = top25, xnudge = 0, 
            ynudge = 0, repel = TRUE)


## ----cellcycle----------------------------------------------------------------
cc.genes <- readLines(file.path(data_dir,
  'regev_lab_cell_cycle_genes_mm.fixed.txt'))
s.genes   <- cc.genes[1:43]
g2m.genes <- cc.genes[44:98]

liver <- CellCycleScoring(liver, 
                          s.features   = s.genes, 
                          g2m.features = g2m.genes, 
                          set.ident    = FALSE)


## ----scaling, message=FALSE---------------------------------------------------
liver <- liver %>%
    ScaleData(vars.to.regress = c("percent.mt", "nFeature_RNA"))


## ----PCA, message=FALSE-------------------------------------------------------
liver <- liver %>%
              RunPCA(verbose = FALSE, npcs = 100)


## ----pcplot-------------------------------------------------------------------
DimPlot(liver, reduction = "pca")


## ----elbow--------------------------------------------------------------------
ElbowPlot(liver, ndims = 100)


## ----elbow2-------------------------------------------------------------------
ElbowPlot(liver, ndims = 50)


## ----elbow3-------------------------------------------------------------------
num_pc <- 24
ElbowPlot(liver, ndims = 40) + geom_vline(xintercept = num_pc)


## ----umap, warning=FALSE------------------------------------------------------
liver <- RunUMAP(liver, reduction = 'pca', dims = 1:num_pc, 
    verbose = FALSE)


## ----seurat3, message=FALSE---------------------------------------------------
liver <- FindNeighbors(liver, reduction = 'pca', 
                       dims = 1:num_pc, verbose = FALSE) %>%
           FindClusters(verbose = FALSE, resolution = 0.3)
UMAPPlot(liver, label = TRUE, label.size = 6)


## ----save_seurat--------------------------------------------------------------
saveRDS(liver, file = file.path(data_dir, 'lesson05.rds'))


## ----session_info-------------------------------------------------------------
sessionInfo()

