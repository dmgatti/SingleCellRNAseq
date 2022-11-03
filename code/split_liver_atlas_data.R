################################################################################
# Read in the liver atlas data and produce two subsets of data which maintain
# the same proportions of cell types and "bad" cells.
#
# 2022-07-26
# Daniel Gatti
# dan.gatti@jax.org
################################################################################

library(R.utils)
library(data.table)
library(tidyverse)
library(Matrix)
library(Seurat)

# Download the data from: https://www.livercellatlas.org/download.php
# Get "Liver Cell Atlas: Mouse StSt": All Liver Cells
# Gene-cell count matrix
# Cell annotation matrix
# Get a node with 64GB to read in file.

src_dir  = '/fastscratch/dgatti/scrnaseq/data'
dest_dir = '/projects/compsci/vmp/USERS/dgatti/projects/scrnaseq/data'

set.seed(295872)

# Read in cell metadata.
metadata = read_csv(file.path(src_dir, 'annot_mouseStStAll.csv'),
                    show_col_types = FALSE)
n_cells = nrow(metadata)

# Read in the counts.
counts = Read10X(file.path(src_dir, 'rawData_mouseStSt/countTable_mouseStSt'),
                 gene.column = 1)

# What are the types of cell digests?
count(metadata, digest)

# What are the proportions of cell types?
count(metadata, annot)

# What are the sample types?
count(metadata, typeSample)

# What proportion of the cells are un-annotated?
(ncol(counts) - nrow(metadata)) / ncol(counts)

# What are the cell type counts by digest?
count(metadata, digest, annot) %>%
  group_by(digest) %>%
  mutate(total = sum(n, na.rm = TRUE),
         prop = n / total) %>%
  select(-total, -n) %>%
  pivot_wider(names_from = digest, values_from = prop)

# Get the sample IDs.
sample_ids = metadata %>%
               count(sample, digest) %>%
               filter(digest != 'nuclei')

# Get barcodes for the cells without metadata.
bad_cells = colnames(counts)[!colnames(counts) %in% metadata$cell]
bad_cells = sample(bad_cells, size = floor(length(bad_cells)) / 4)
bad_cells = data.frame(UMAP_1  = 0,
                       UMAP_2  = 0,
                       cluster = 0, 
                       annot   = 'bad',
                       sample  = NA_character_,
                       cell    = bad_cells,
                       digest  = sample(c('inVivo', 'exVivo'), size = length(bad_cells), replace = TRUE),
                       typeSample = sample('scRnaSeq', 
                                           size = length(bad_cells), replace = TRUE))

# Add sample IDs to the inVivo and exVivo data.
wh_invivo = which(bad_cells$digest == 'inVivo')
bad_cells$sample[wh_invivo] = sample(sample_ids$sample[sample_ids$digest == 'inVivo'],
                                     size = length(wh_invivo), replace = TRUE)

wh_exvivo = which(bad_cells$digest == 'exVivo')
bad_cells$sample[wh_exvivo] = sample(sample_ids$sample[sample_ids$digest == 'exVivo'],
                                     size = length(wh_exvivo), replace = TRUE)

rm(wh_invivo, wh_exvivo)

# Add the bad cells to metadata and subset counts.
metadata = bind_rows(metadata, bad_cells)
counts   = counts[,metadata$cell]
rm(bad_cells)

# Remove citeSeq data.
metadata = subset(metadata, typeSample != 'citeSeq')
counts   = counts[,metadata$cell]

# Write single cell data.
write_sc = function(meta, counts, out_dir) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Write metadata.
  write_csv(meta, file = file.path(out_dir, 'annot_metadata.csv'))

  # Remove UMAP and cell type annotation from metadata file.
  meta = meta[,c('sample', 'cell', 'digest', 'typeSample')]
  write_csv(meta, file = file.path(out_dir, 'annot_metadata_first.csv'))
  
  # Note: 3 files: barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
  barcode_file = file.path(out_dir, 'barcodes.tsv')
  feature_file = file.path(out_dir, 'features.tsv')
  matrix_file  = file.path(out_dir, 'matrix.mtx')
  writeLines(text = colnames(counts), con = barcode_file, sep = '\n')
  writeLines(text = rownames(counts), con = feature_file, sep = '\n')
  Matrix::writeMM(obj = counts, file = matrix_file)
  gzip(barcode_file, overwrite = TRUE)
  gzip(feature_file, overwrite = TRUE)
  gzip(matrix_file,  overwrite = TRUE)

} # write_sc()

# Verify that we have 'bad' samples in the metadata.
stopifnotsum(metadata$annot == 'bad' > 0)
# Verifyg that we have only 'scRnaSeq' samples.
stopifnot(all(metadata$typeSample == 'scRnaSeq'))
# Verify that we only have inVivo and esVivo data.
stopifnot(all(metadata$digest %in% c('inVivo', 'exVivo')))

# Split invivo and exvivo data.
meta_iv   = subset(metadata, digest == 'inVivo')
counts_iv = counts[,meta_iv$cell]
write_sc(meta_iv, counts_iv, file.path(dest_dir, 'mouseStSt_invivo'))
rm(meta_iv, counts_iv)

meta_ev   = subset(metadata, digest == 'exVivo')
counts_ev = counts[,meta_ev$cell]
write_sc(meta_ev, counts_ev, file.path(dest_dir, 'mouseStSt_exvivo'))



################################################################################
# DMG: The code below was used for sub-sampling. We may not need to do that now.

# Subset the scRNAseq data, retaining the proportion of cells.
# Write results to the given directory.
# Assumes data read in using Seurat::Read10X().
# meta: data.frame containing cell metadata.
# counts: sparse matrix containing counts. Genes in rows, cell in columns.
# prop: float between 0 and 1 indicating the proportion of the cells to keep.
# out_dir: string containing the full path to the output directory. 
subsample_scrnaseq = function(meta, counts, prop, out_dir) {

  n_cells = nrow(meta)
  new_n_cells = round(n_cells * prop)

  cell_prop = meta %>%
                count(annot) %>%
                mutate(prop = n / n_cells)

  print(cell_prop)

  # Vector of cell type annotations, sampled according to abundances.
  keep = sample(cell_prop$annot, size = new_n_cells, replace = TRUE, 
                prob = cell_prop$prop)

  # Vector of cells to keep.
  cell_keep = NULL

  # Select cells to keep.
  for(i in 1:nrow(cell_prop)) {

    curr_type = cell_prop$annot[i]
    wh   = which(meta$annot == curr_type)
    rows = sample(wh, size = min(sum(keep == curr_type), length(wh)))
    stopifnot(unique(meta$annot[wh]) == curr_type)
    cell_keep = c(cell_keep, meta$cell[rows])
  
  } # for(i)

  stopifnot(length(cell_keep) == new_n_cells)
  stopifnot(sum(duplicated(cell_keep)) == 0)
  
  print(paste('Num cells:', length(cell_keep)))

   # Create output directory.
   dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
   
   # Subset and write metadata.
   meta_ss = meta[match(cell_keep, meta$cell),]
   stopifnot(all(cell_keep == meta_ss$cell))
   write_csv(meta_ss, file = file.path(out_dir, 'annot_metadata_paper.csv'))

   # Remove UMAP and annotation columns and write out file for
   # students to start with.
   meta_ss = meta_ss[,]
   write_csv(meta_ss, file = file.path(out_dir, 'annot_metadata.csv'))
   rm(meta_ss)
   
   # Subset counts and write out.
   # Note: 3 files: barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
   barcode_file = file.path(out_dir, 'barcodes.tsv')
   feature_file = file.path(out_dir, 'features.tsv')
   matrix_file  = file.path(out_dir, 'matrix.mtx')
   counts_ss    = counts[,cell_keep]
   writeLines(text = colnames(counts_ss), con = barcode_file, sep = '\n')
   writeLines(text = rownames(counts_ss), con = feature_file, sep = '\n')
   Matrix::writeMM(obj = counts_ss, file = matrix_file)
   rm(counts_ss)
   gzip(barcode_file, overwrite = TRUE)
   gzip(feature_file, overwrite = TRUE)
   gzip(matrix_file,  overwrite = TRUE)

} # subsample_scrnaseq()

# Filter to retain in vivo data and subset.
meta_ss   = filter(metadata, digest == 'inVivo')
counts_ss = counts[,meta_ss$cell]
subsample_scrnaseq(meta = meta_ss, counts = counts_ss, prop = 1.0, 
        out_dir = file.path(dest_dir, 'mouseStSt_invivo'))

# Filter to retain ex vivo data and subset.
meta_ss   = filter(metadata, digest == 'exVivo')
counts_ss = counts[,meta_ss$cell]
subsample_scrnaseq(meta = meta_ss, counts = counts_ss, prop = 1.0, 
        out_dir = file.path(dest_dir, 'mouseStSt_exvivo'))




