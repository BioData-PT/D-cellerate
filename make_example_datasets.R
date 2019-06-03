library(Matrix)
library(Seurat)

# example datasets
example.datasets <- list(

  mca_lung_2000_50 = list(
    dataframe = Matrix(as.matrix(read.table("data/mca_lung1_2000_50.txt", header=TRUE)), sparse=TRUE),
    description = "Sample of lung cells from the Mouse Cell Atlas paper. Only cells with at least 2000 UMI and 50 genes expressed."
  ),
  mca_lung_1000_10 = list(
    dataframe = Matrix(as.matrix(read.table("data/mca_lung1_1000_10.txt", header=TRUE)), sparse=TRUE),
    description = "Sample of lung cells from the Mouse Cell Atlas paper. Only cells with at least 1000 UMI and 10 genes expressed."
  ),
  mca_lung_500 = list(
    dataframe = Matrix(as.matrix(read.table("data/mca_lung1_500.txt", header=TRUE)), sparse=TRUE),
    description = "Sample of lung cells from the Mouse Cell Atlas paper. Only cells with at least 500 UMI."
    ),
  mca_lung_full = list(
    dataframe = (function() {
      df <- read.table(gzfile("data/lung1_10k.dge.txt.gz"), header=TRUE)
      rownames(df) <- df$GENE
      df <- df[, -1]
      Matrix(as.matrix(df), sparse=TRUE)
    })(),
    description = "Sample of lung cells from the Mouse Cell Atlas paper. 10k barcodes."
  ),
  pbmc4k_filtered = list(
    dataframe = Read10X_h5("data/filtered_gene_bc_matrices_h5.h5"),
    description = "Filtered PBMC4k from 10x Genomics."
  ),
  pbmc4k_raw = list(
    dataframe = Read10X_h5("data/raw_gene_bc_matrices_h5.h5"),
    description = "Raw PBMC4k from 10x Genomics containing all barcodes."
  )
  
)

saveRDS(example.datasets, file = "data/example_datasets.rds")
