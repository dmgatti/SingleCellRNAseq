## ---- include=FALSE-----------------------------------------------------------
source("../bin/chunk-options.R")
knitr_fig_path("06-")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(enrichR))

data_dir <- '../data'


## ----seed---------------------------------------------------------------------
# set a seed for reproducibility in case any randomness used below
set.seed(1418)


## ----load_data----------------------------------------------------------------
liver <- readRDS(file.path(data_dir, 'lesson05.rds'))


## ----table_smple_clusters-----------------------------------------------------
table(liver$sample, liver$seurat_clusters)


## ----sample_effects, fig.width = 7, fig.height = 6----------------------------
UMAPPlot(liver, group.by = 'sample', pt.size = 0.1)


## ----find_markers1------------------------------------------------------------
markers13 <- FindMarkers(liver, '13', only.pos = TRUE, logfc.threshold = 1,
                         max.cells.per.ident = 500)
head(markers13, 6)


## ----cd79a_vln, fig.width=7, fig.height=4-------------------------------------
VlnPlot(liver, 'Cd79a')


## ----cd79a_fp, fig.width = 7, fig.height = 6----------------------------------
FeaturePlot(liver, "Cd79a", cols = c('lightgrey', 'red'), 
            label = TRUE, label.size = 6)


## ----c21----------------------------------------------------------------------
table(liver$sample[liver$seurat_clusters == '21'])


## ----harmony, message = FALSE, warning = FALSE--------------------------------
# Store old UMAP and old clusters
liver$before_harmony_clusters <- liver$seurat_clusters
liver@misc$noharmony_umap <- liver@reductions$umap

# Run harmony
liver <- RunHarmony(liver, 'sample', assay.use='RNA',
           theta=1, dims.use=1:40, max.iter.harmony=100)
ElbowPlot(liver, reduction = 'harmony', ndims = 40)


## ----finish_harmony, warning = FALSE------------------------------------------
liver <- FindNeighbors(liver, reduction = 'harmony', dims = 1:24) %>%
         FindClusters(verbose = FALSE, resolution = 0.3) %>%
         RunUMAP(dims = 1:24, reduction = 'harmony')
liver$after_harmony_clusters <- liver$seurat_clusters


## ----c1321--------------------------------------------------------------------
table(liver$before_harmony_clusters, 
      liver$after_harmony_clusters)[c('13', '21'), ]


## ----c8, fig.width = 7, fig.height = 6----------------------------------------
FeaturePlot(liver, 'Cd79a', cols = c('lightgrey', 'red'), label = T, 
            label.size = 6)


## ----c9, fig.width = 8, fig.height = 4----------------------------------------
VlnPlot(liver, 'Cd79a')


## ----renaming_clusters--------------------------------------------------------
genes <- c('Socs3', 'Gnmt', 'Timd4', 'Ms4a4b', 'S100a4',
           'Adgrg6', 'Cd79a', 'Dck', 'Siglech', 'Dcn', 
           'Wdfy4', 'Vwf', 'Spp1', 'Hdc')
Idents(liver) <- 'after_harmony_clusters'
a <- AverageExpression(liver, features = genes)[['RNA']]
highest_clu <- unname(colnames(a))[apply(a, 1, which.max)]

cluster_converter <- setNames(paste0('c', c(0:1, 3:14)), highest_clu)
remaining_clusters <- setdiff(as.character(unique(Idents(liver))),
                              highest_clu)
a <- AverageExpression(subset(liver, idents = remaining_clusters),
                       features = c("Lyve1", "Cd5l"))[['RNA']]
highest_clu <- unname(colnames(a))[apply(a, 1, which.max)]
cluster_converter[highest_clu] <- c('c2', 'c15')

liver$renamed_clusters <- cluster_converter[as.character(liver$after_harmony_clusters)]
Idents(liver) <- 'renamed_clusters'


## ----markers, message=FALSE---------------------------------------------------
liver_mini <- subset(liver, downsample = 300)
markers <- FindAllMarkers(liver_mini, only.pos = TRUE, 
                          logfc.threshold	= log2(1.25), min.pct = 0.2) 


## ----markers_massage----------------------------------------------------------
old_markers <- markers
markers <- as_tibble(markers) %>% 
              select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj)
head(markers, 6)


## ----top_markers--------------------------------------------------------------
group_by(markers, cluster) %>% 
  top_n(3, avg_log2FC) %>%
  mutate(rank = 1:n()) %>%
  pivot_wider(-c(avg_log2FC, pct.1, pct.2, p_val_adj), 
              names_from = 'rank', values_from = 'gene') %>%
  arrange(cluster)


## ----top_markers2, fig.width = 8, fig.height = 10-----------------------------
top_markers <- group_by(markers, cluster) %>% 
                 arrange(desc(avg_log2FC)) %>%
                 top_n(1, avg_log2FC) %>% pull(gene)
VlnPlot(liver, features = top_markers, stack = TRUE, flip = TRUE)


## ----expr---------------------------------------------------------------------
UMAPPlot(liver, label = TRUE, label.size = 6) + NoLegend()


## ----expr2--------------------------------------------------------------------
FeaturePlot(liver, features = "Kdr", cols = c('lightgrey', 'red'))


## ----expr3--------------------------------------------------------------------
FeaturePlot(liver, "Fabp1", cols = c('lightgrey', 'red'),
            label = TRUE, label.size = 6)


## ----expr4--------------------------------------------------------------------
FeaturePlot(liver, "Serpina1a", cols = c('lightgrey', 'red'),
            label = TRUE, label.size = 6)


## ----sp_markers, fig.width = 8, fig.height = 10-------------------------------
specific_markers <- group_by(markers, cluster) %>% 
  arrange(desc(avg_log2FC)) %>%
  filter(pct.2 < 0.2) %>%
  arrange(cluster) %>%
  top_n(1, avg_log2FC) %>% pull(gene)
VlnPlot(liver, features = specific_markers, stack = TRUE, flip = TRUE)


## ----vsig4--------------------------------------------------------------------
FeaturePlot(liver, "Vsig4", cols = c('lightgrey', 'red'), label = TRUE,
            label.size = 6)


## ----adgre1-------------------------------------------------------------------
FeaturePlot(liver, "Adgre1", cols = c('lightgrey', 'red'), label = TRUE,
            label.size = 6)


## ----doublets, fig.width = 8, fig.height = 4----------------------------------
VlnPlot(liver, c("Adgre1", "Fabp1"), idents = c('c3', 'c15', 'c1'), sort = T)


## ----labelling----------------------------------------------------------------
labels <- tibble(cluster_num = unique(liver$renamed_clusters)) %>%
  mutate(cluster_num = as.character(cluster_num)) %>%
  mutate(cluster_name = case_when(
         cluster_num %in% c('c0', 'c2', 'c6', 'c12') ~ 'ECs',   # endothelial cells
         cluster_num == 'c1' ~ 'hepatocytes',
         cluster_num %in% c('c3', 'c8') ~ 'Kupffer cells',
         cluster_num == 'c4' ~ 'T cells',
         cluster_num == 'c7' ~ 'B cells',
         cluster_num == 'c9' ~ 'pDCs',               # plasmacytoid dendritic cells
         cluster_num == 'c14' ~ 'neutrophils',
         cluster_num == 'c15' ~ 'KH doub.',          # Kupffer-hepatocyte doublets
         TRUE ~ cluster_num))

liver$labels <- deframe(labels)[as.character(liver$renamed_clusters)]
UMAPPlot(liver, label = TRUE, label.size = 6, group.by = 'labels') + NoLegend()


## ----fake, fig.height = 4.5, fig.width = 10-----------------------------------
libraries <- unique(liver$sample)
treatment_group <- setNames(c(rep('control', 5), rep('drug', 4)), libraries)
liver$trt <- treatment_group[liver$sample]

hepatocytes <- subset(liver, labels == "hepatocytes")
Idents(hepatocytes) <- "trt"
UMAPPlot(hepatocytes, split.by = 'trt', group.by = 'labels', label = T,
         label.size = 6)


## ----deg1---------------------------------------------------------------------
deg1 <- FindMarkers(hepatocytes, ident.1 = 'drug', ident.2 = 'control',
                    logfc.threshold = 0.2, only.pos = FALSE)


## ----deg_res------------------------------------------------------------------
head(deg1, 10)


## ----pseudobulk1--------------------------------------------------------------
# Make pseudobulks.
pseudobulk <- AggregateExpression(hepatocytes, slot = 'counts', 
                                  group.by = 'sample', assays = 'RNA')[['RNA']]
dim(pseudobulk)
head(pseudobulk, 6)

# Run DESeq2
pseudobulk_metadata <- hepatocytes[[c("sample", "trt")]] %>%
  as_tibble() %>% distinct() %>% as.data.frame() %>%
  column_to_rownames('sample') %>%
  mutate(trt = as.factor(trt))
pseudobulk_metadata <- pseudobulk_metadata[colnames(pseudobulk), , drop = F]
dds <- DESeqDataSetFromMatrix(pseudobulk, 
                              colData = pseudobulk_metadata, 
                              design = ~ trt)
trt <- DESeq(dds, test = "LRT", reduced = ~ 1)
res1 <- results(trt)
head(res1)
sum(!is.na(res1$padj) & res1$padj < 0.05)


## ----pway, fig.width = 9, fig.height = 5--------------------------------------
db_names <- c("KEGG"='KEGG_2019_Mouse',
              "GO"='GO_Biological_Process_2021',
              "MsigDB"='MSigDB_Hallmark_2020')
genes <- filter(markers, cluster == 'c14') %>%
  top_n(100, avg_log2FC) %>% pull(gene)
enrich_genes <- enrichr(genes, databases = db_names)
names(enrich_genes) <- names(db_names)
e <- bind_rows(enrich_genes, .id = 'database') %>%
  mutate(Term = paste0(database, ': ', Term))
plotEnrich(e, title = "Neutrophil pathway enrichment", 
           showTerms = 15, numChar = 50)


## ----session_info,collapse=TRUE-----------------------------------------------
sessionInfo()

