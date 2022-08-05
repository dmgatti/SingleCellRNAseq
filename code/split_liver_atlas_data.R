################################################################################
# Read in the liver atlas data and produce two subsets of data which maintain
# the same proportions of cell types and "bad" cells.
# 2022-07-26
# Daniel Gatti
# dan.gatti@jax.org
################################################################################

library(R.utils)
library(data.table)
library(tidyverse)
library(Matrix)
library(Seurat)

src_dir  = '/fastscratch/dgatti'
dest_dir = '/projects/compsci/USERS/dgatti/projects/scrnaseq/data'

# Read in cell metadata.
metadata = read_csv(file.path(src_dir, 'annot_mouseStStAll.csv'),
                    show_col_types = FALSE)
n_cells = nrow(metadata)

# Read in the counts.
counts = Read10X(file.path(src_dir, 'rawData_mouseStSt/countTable_mouseStSt'),
                 gene.column = 1)

# What are the types of cell digets?
count(metadata, digest)

# What are the proportions of cell types?
count(metadata, annot)

# What proportion of the cells are un-annotated?
(ncol(counts) - nrow(metadata)) / ncol(counts)

# What are the cell type counts by digest?
count(metadata, digest, annot) %>%
  group_by(digest) %>%
  mutate(total = sum(n, na.rm = TRUE),
         prop = n / total) %>%
  select(-total, -n) %>%
  pivot_wider(names_from = digest, values_from = prop)

# Get barcodes for the cells without metadata and keep  half of them.
bad_cells = colnames(counts)[!colnames(counts) %in% metadata$cell]
bad_cells = sample(bad_cells, size = floor(length(bad_cells) / 2))
bad_cells = data.frame(UMAP_1  = 0,
                       UMAP_2  = 0,
                       cluster = 0, 
                       annot   = 'bad',
                       sample  = NA_character_,
                       cell    = bad_cells,
                       digest  = NA_character_,
                       typeSample = sample(c('scRnaSeq', 'citeSeq'), 
                                           size = length(bad_cells), replace = TRUE))

# Add the bad cells to metadata and subset counts.
metadata = bind_rows(metadata, bad_cells)
counts   = counts[,metadata$cell]
rm(bad_cells)

# Subset the scRNAseq data, retaining the proportion of cells.
# Write results to the given directory.
# Assumes data read in using Seurat::Read10X().
# meta: data.frame containing cell metadata.
# counts: sparce matrix containing counts. Genes in rows, cell in columns.
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
    rows = sample(wh, size = sum(keep == curr_type))
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

#subsample_scrnaseq(meta = metadata, counts = counts, prop = 0.5, 
#        out_dir = file.path(dest_dir, 'mouseStSt_50'))

# Filter to retain scRNAseq data and subset.
meta_ss   = filter(metadata, typeSample %in% c('scRnaSeq', 'nucSeq'))
counts_ss = counts[,meta_ss$cell]
subsample_scrnaseq(meta = meta_ss, counts = counts_ss, prop = 0.75, 
        out_dir = file.path(dest_dir, 'mouseStSt_scrnaseq_75pct'))

# Filter to retain citeseq data and subset.
meta_ss   = filter(metadata, typeSample == 'citeSeq')
counts_ss = counts[,meta_ss$cell]
subsample_scrnaseq(meta = meta_ss, counts = counts_ss, prop = 0.75, 
        out_dir = file.path(dest_dir, 'mouseStSt_citeseq_75pct'))


#subsample_scrnaseq(meta = metadata, counts = counts, prop = 0.10, 
#        out_dir = file.path(dest_dir, 'mouseStSt_10'))



