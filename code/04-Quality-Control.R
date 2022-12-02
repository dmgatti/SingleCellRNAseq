## ---- setup, include=FALSE----------------------------------------------------
source("../bin/chunk-options.R")
knitr_fig_path("04-")


## ----libs, warning=FALSE------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scds))
suppressPackageStartupMessages(library(Seurat))

data_dir <- '../data'


## ----load_data,include=FALSE--------------------------------------------------
# Only needed to build lessons.
# Don't use in course.
load(file.path(data_dir, 'lesson03.Rdata'))


## ----scds1--------------------------------------------------------------------
cell_ids <- filter(metadata, sample == 'CS52') %>% pull(cell)
sce <- SingleCellExperiment(list(counts = counts[, cell_ids]))
sce <- cxds_bcds_hybrid(sce)
doublet_preds <- colData(sce)


## ----gc, echo = FALSE, message = FALSE----------------------------------------
rm(sce)
gc()


## ----zero_gene_counts---------------------------------------------------------
gene_counts <- Matrix::rowSums(counts)
sum(gene_counts == 0)


## ----filter_gene_by_counts----------------------------------------------------
counts <- counts[gene_counts > 0,]


## ----gene_count_hist----------------------------------------------------------
gene_counts <- tibble(counts  = Matrix::rowSums(counts > 0))

gene_counts %>% 
  ggplot(aes(counts)) +
    geom_histogram(bins = 1000) +
    labs(title = 'Histogram of Number of Cells in which Gene was Detected',
         x     = 'Number of Genes',
         y     = 'Number of Cells in which Gene was Detected') +
  theme_bw(base_size = 14)


## ----gene_count_hist_2--------------------------------------------------------
gene_counts %>%
  ggplot(aes(counts)) +
    geom_histogram(bins = 1000) +
    labs(title = 'Histogram of Number of Cells in which Gene was Detected',
         x     = 'Number of Genes',
         y     = 'Number of Cells in which Gene was Detected') +
    lims(x = c(0, 50)) +
    theme_bw(base_size = 14) +
    annotate('text', x = 2, y = 1500, hjust = 0,
             label = str_c(sum(gene_counts == 1), ' genes were detected in only one cell')) +
    annotate('text', x = 3, y = 900, hjust = 0,
             label = str_c(sum(gene_counts == 2), ' genes were detected in two cells'))


## ----sum_cell_counts----------------------------------------------------------
tibble(counts  = Matrix::colSums(counts > 0)) %>%
  ggplot(aes(counts)) +
    geom_histogram(bins = 500) +
    labs(title = 'Histogram of Number of Genes per Cell',
         x     = 'Number of Genes with Counts > 0',
         y     = 'Number of Cells')


## ----seed---------------------------------------------------------------------
# set a seed for reproducibility in case any randomness used below
set.seed(1418)


## ----create_seurat_obj--------------------------------------------------------
metadata <- as.data.frame(metadata) %>%
              column_to_rownames('cell')
liver <- CreateSeuratObject(counts    = counts, 
                            project   = 'liver: scRNA-Seq',
                            meta.data = metadata,
                            min.cells = 5)


## ----add_doublets-------------------------------------------------------------
liver <- AddMetaData(liver, as.data.frame(doublet_preds))


## ----get_assays---------------------------------------------------------------
Seurat::Assays(liver)


## ----get_assay_data-----------------------------------------------------------
tmp = GetAssayData(object = liver, slot = 'data')
tmp[1:5,1:5]


## ----show_meta----------------------------------------------------------------
head(liver[[]])


## ----pct_mito-----------------------------------------------------------------
liver <- liver %>% 
              PercentageFeatureSet(pattern = "^mt-", col.name = "percent.mt")


## ----seurat_counts_plots------------------------------------------------------
VlnPlot(liver, features = "percent.mt", group.by = 'sample')


## ----seurat_counts_plots2,fig.height=4,fig.width=6----------------------------
VlnPlot(liver, features = "percent.mt", group.by = 'sample', pt.size = 0)


## ----mito_by_cell_type, include = FALSE, eval = FALSE-------------------------
## # DMG made the plot below to see if there are mitochondrial expression differences by annotated cell type. The students won't have this file at this stage of the analysis. But how do we discuss these differences? Since we're looking for high values, it may not be too important.
## liver[[c('annot', 'percent.mt')]] %>%
##     ggplot(aes(annot, percent.mt + 0.01)) +
##       geom_boxplot() +
##       scale_y_log10() +
##       coord_flip()


## ----subset_by mito-----------------------------------------------------------
#liver <- subset(liver, subset = percent.mt < 14)


## ----filter_gene_counts,fig.height=4,fig.width=6------------------------------
VlnPlot(liver, 'nFeature_RNA', group.by = 'sample', pt.size = 0)


## ----filter_gene_counts_5k,fig.height=4,fig.width=6---------------------------
VlnPlot(liver, 'nFeature_RNA', group.by = 'sample', pt.size = 0) +
  scale_y_log10() + 
  geom_hline(yintercept = 600) + 
  geom_hline(yintercept = 5000)
#liver <- subset(liver, nFeature_RNA > 600 & nFeature_RNA < 5000)


## ----genes_umi----------------------------------------------------------------
ggplot(liver@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point() +
  theme_bw(base_size = 16) +
  xlab("nUMI") + ylab("nGenes") +
  scale_x_log10() + scale_y_log10()


## ----filter_umi,fig.height=4,fig.width=6--------------------------------------
VlnPlot(liver, 'nCount_RNA', group.by = 'sample', pt.size = 0) +
  scale_y_log10() + 
  geom_hline(yintercept = 900) + 
  geom_hline(yintercept = 25000)
#liver <- subset(liver, nCount_RNA > 900 & nCount_RNA < 25000)


## ----liver_outliers, include = FALSE, eval = FALSE----------------------------
## # This is modeled after what the authors of the liver cell atlas did
## library(scuttle)
## genes <- scuttle::isOutlier(liver$nFeature_RNA, nmads = 3,
##                             batch = liver$sample, log=TRUE)
## umi <- scuttle::isOutlier(liver$nCount_RNA, nmads = 3,
##                           batch = liver$sample, log=TRUE)
## mt <- scuttle::isOutlier(liver$percent.mt, nmads = 3, log = FALSE,
##                          batch = liver$sample, type = 'lower')
## tapply(genes | umi | mt, liver$annot, mean, na.rm=T)
## attr(umi, 'thresholds')


## ----filtering----------------------------------------------------------------
liver$keep <- with(liver, percent.mt < 14 & nFeature_RNA > 600 &
  nFeature_RNA < 5000 & nCount_RNA > 900 & nCount_RNA < 25000)


## ----doublet_plot-------------------------------------------------------------
ggplot(mutate(liver[[]], class = ifelse(keep, 'QC singlet', 'QC doublet')),
  aes(x = class, y = hybrid_score)) + 
  geom_violin() + theme_bw(base_size = 18) +
  xlab("") + ylab("SCDS hybrid score")


## ----subsetting---------------------------------------------------------------
liver <- subset(liver, subset = percent.mt < 14 & nFeature_RNA > 600 &
  nFeature_RNA < 5000 & nCount_RNA > 900 & nCount_RNA < 25000)


## ----save_seurat,include=FALSE------------------------------------------------
# Only needed to build lessons.
# Don't use in course.
saveRDS(liver, file = file.path(data_dir, 'lesson04.rds'))


## ----session_info,collapse=TRUE-----------------------------------------------
sessionInfo()

