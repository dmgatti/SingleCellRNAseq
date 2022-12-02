## ----setup, include=FALSE-----------------------------------------------------
source("../bin/chunk-options.R")
knitr_fig_path("03-")


## ----libs,warning=FALSE-------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))

data_dir <- '../data'


## ----read_counts--------------------------------------------------------------
# uses the Seurat function Read10X()
counts <- Read10X(file.path(data_dir, 'mouseStSt_invivo'), gene.column = 1)


## ----dim_counts---------------------------------------------------------------
dim(counts)


## ----rownames_counts----------------------------------------------------------
head(rownames(counts), n = 10)


## ----duplicated_symbol--------------------------------------------------------
sum(duplicated(rownames(counts)))


## ----colnames_counts----------------------------------------------------------
head(colnames(counts), n = 10)


## ----barcodes_unique----------------------------------------------------------
sum(duplicated(colnames(counts)))


## ----view_counts--------------------------------------------------------------
counts[1:10, 1:20]


## ----counts_class-------------------------------------------------------------
str(counts)


## ----counts_image,fig.height=6------------------------------------------------
image(1:100, 400:600, t(as.matrix(counts[400:600,1:100]) > 0), 
      xlab = 'Cells', ylab = 'Genes')


## ----gene_sums----------------------------------------------------------------
gene_sums <- data.frame(gene_id = rownames(counts),
                        sums    = Matrix::rowSums(counts))
sum(gene_sums$sums == 0)


## ----cell_counts--------------------------------------------------------------
hist(Matrix::colSums(counts))

Matrix::colSums(counts) %>% enframe() %>%
  ggplot(aes(value)) + geom_histogram(bins = 30) + 
  theme_bw(base_size = 16) + scale_x_log10()


## ----read_metadata------------------------------------------------------------
metadata <- read_csv(file.path(data_dir, 'mouseStSt_invivo', 'annot_metadata_first.csv'))


## ----head_metadata------------------------------------------------------------
head(metadata)


## ----cell_classes-------------------------------------------------------------
dplyr::count(metadata, digest, typeSample)


## ----save_data,include=FALSE--------------------------------------------------
save(counts, metadata, file = file.path(data_dir, 'lesson03.Rdata'))


## ----session_info,collapse=TRUE-----------------------------------------------
sessionInfo()

