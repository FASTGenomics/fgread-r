canonize_id <- function(metadata, regex, to_name) {
  cols <- names(metadata)
  cell_id_names <- grep(regex, cols, value = T, ignore.case = T)
  if (length(cell_id_names) > 0) {
    cell_id_name <- sort(cell_id_names)[1]
    names(metadata)[cols == cell_id_name] <- to_name
  } else {
    stop(glue::glue("Looking for {to_name}, could not find a column matching {regex}."))
  }
  return(metadata)
}

#' a helper function that constructs a seurat object from an expression matrix, cell &
#' gene metadata.  Used to read loom and AnnData files
#' 
#' @param matrix The path to the expression matrix file.
#' 
#' @param cell_metadata The path to the cell metadata file.
#' 
#' @param gene_metadata The path to the gene metadata file.
#' 
#' @return A Seurat object
#' 
matrix_to_seurat <- function(matrix, cell_metadata, gene_metadata) {
  seurat <- Seurat::CreateSeuratObject(counts = matrix, min.cells = 0, min.features = 0)
  Seurat::AddMetaData(seurat, metadata = cell_metadata)
  seurat@assays$RNA@meta.features <- gene_metadata
  return(seurat)
}

#' Read a Seurat object.
#' 
#' @param ds_file The path to the dataset file.
#' 
#' @return A Seurat Object
#' 
read_seurat_to_seurat <- function(ds_file) {
  return(readRDS(ds_file))
}

#' Reading a loom dataset in Seurat v3 uses loomR and may not work as expected.
#' For your convenience we implemented experimental readers that you can use by setting "experimental_readers=TRUE"
#' in \code{\link{load_data}}.
#'
#' @param ds_file The path to the dataset file.
#' 
#' @return A seurat Object.
#' 
read_loom_to_seurat <- function(ds_file) {

  warning(glue::glue(
    '!!!\n',
    'Importing Loom is not always working as expected in Seurat v3.\n',
    'If you encounter any problems you can use our FASTGenomics experimental ',
    'reader by setting "experimental_readers=TRUE" in `load_data`.\n',
    'For more information please see {DOCSURL}.\n',
    '!!!\n'))
  loom <- loomR::connect(filename = ds_file, mode = "r")
  return(Seurat::as.Seurat(loom))
}

#' Read a loom dataset with the experimental FASTGenomics reader.
#' Importing Loom is not always working as expected in Seurat v3.
#' This beta loader provided by FASTGenomics only reads the count table
#' and the row/col attributes that can fit in a data frame structure (i.e. all higher
#' dimensional attributes are discarded). To keep the implementation simple we read the
#' whole object into memory, including the dense count matrix.  This could be a
#' potential bottleneck for larger datasets.
#' 
#' @param ds_file The path to the dataset file.
#' 
#' @return A seurat Object.
#' 
read_loom_to_seurat_exp <- function(ds_file) {

  warning(
      "!!!\n",
      "Importing Loom is not always working as expected in Seurat v3.\n",
      "For your convenience the FASTGenomics team provides this beta loading routine.\n",
      "In case of problems please consider using another format or implement your own loading routine.\n",
      "!!!\n",
      call. = TRUE, immediate. = TRUE)

  file <- rhdf5::H5Fopen(ds_file, flags = "H5F_ACC_RDONLY")
  contents <- rhdf5::h5dump(file)
  rhdf5::H5Fclose(file)

  cleanup <- function(cols) {
    Map(as.vector, Filter(function(x) length(dim(x)) == 1, cols))
  }

  gene_metadata <- data.frame(cleanup(contents$row_attrs))
  cell_metadata <- data.frame(cleanup(contents$col_attrs))

  ## this is necessary because the names for cell/gene ids changed over time for the
  ## loom format.  This way we make sure that the at least the cell/gene ids column
  ## has a consistent naming.  Both fields are also required to create a Seurat
  ## object.
  cell_metadata <- canonize_id(cell_metadata, regex = "^(cell[^a-z]*(id|name)?$|name|obs_names)$", to_name = "fg_cell_id")
  gene_metadata <- canonize_id(gene_metadata, regex = "^(gene[^a-z]*(id|name)?$|name|var_names)$", to_name = "fg_gene_id")

  rownames(gene_metadata) <- gene_metadata[["fg_gene_id"]]
  rownames(cell_metadata) <- cell_metadata[["fg_cell_id"]]

  matrix <- Matrix::Matrix(t(contents$matrix), sparse = T)
  dimnames(matrix) <- list(rownames(gene_metadata), rownames(cell_metadata))

  return(matrix_to_seurat(matrix, cell_metadata, gene_metadata))
}

#' Read AnnData to Seurat with Seurat's function \code{\link{ReadH5AD}}.
#' Please note that this might not work as expected.
#' For your convenience we implemented importing with SeuratDisk for h5ad files from Anndata>=0.7.3.
#' Alternatively, you can use our experimental reader by setting "experimental_readers=TRUE"
#' in \code{\link{load_data}}.
#' 
#' @param ds_file The path to the dataset file.
#' 
#' @return A seurat Object.
#' 
read_anndata_to_seurat <- function(ds_file) {
  warning(glue::glue(
    "!!!\n",
    'Importing AnnData is not always working as expected in Seurat v3.\n',
    'If you encounter any problems you can use our FASTGenomics experimental ',
    'reader by setting "experimental_readers=TRUE" in `load_data`.\n',
    'For more information please see {DOCSURL}.\n',
    "!!!\n"))

  intermediate <- file.path(getwd(), "intermediate_data.h5seurat")
  seurat <- tryCatch({
    return(Seurat::ReadH5AD(ds_file))
  }, error = function(e) {
    warning(
      "!!!\n",
      'Import with `Seurat::ReadH5AD` failed.\n',
      'Trying again with the `SeuratDisk` library.\n',
      "!!!\n")
    print(glue::glue('Converting h5ad to intermediate file `{intermediate}`.\n'))
    SeuratDisk::Convert(ds_file, dest = intermediate, overwrite = T)
    seurat <- SeuratDisk::LoadH5Seurat(intermediate)
    unlink(intermediate)
    return(seurat)
  })
  return(seurat)
}


#' Read AnnData to Seurat with experimental FASTGenomics reader.
#' Importing AnnData is not always working as expected in Seurat v3
#' For your convenience the FASTGenomics team provides this beta loading routine.
#' This custom reader, however, only reads the `.X`, `.obs` and `.var` components of the AnnData object.
#' 
#' @param ds_file The path to the dataset file.
#' 
#' @return A seurat Object.
#' 
read_anndata_to_seurat_exp <- function(ds_file) {

  warning(
    "!!!\n",
    "Importing AnnData is not always working as expected in Seurat v3.\n",
    "For your convenience the FASTGenomics team provides this beta loading routine.\n",
    "However, this reader only reads the `.X`, `.obs` and `.var` components of the AnnData object.\n",
    "!!!\n",
    call. = TRUE, immediate. = TRUE)

  file <- rhdf5::H5Fopen(ds_file, flags = "H5F_ACC_RDONLY")
  contents <- rhdf5::h5dump(file)
  rhdf5::H5Fclose(file)

  cell_metadata <- contents$obs
  gene_metadata <- contents$var
  if (class(cell_metadata) == "list") {
    stop("The data file seems to be from Anndata>0.7.2.\n",
            "Please use `experimental_reders = F`.")
  }

  if (class(contents$X) == "list") {
    matrix <- Matrix::sparseMatrix(
            i = contents$X$indices + 1,
            p = contents$X$indptr,
            x = as.vector(contents$X$data),
            dims = c(dim(gene_metadata)[1], dim(cell_metadata)[1]),
            dimnames = list(rownames(gene_metadata), rownames(cell_metadata))
        )
  } else if (class(contents$X) == "matrix") {
    matrix <- contents$X
    dimnames(matrix) <- list(rownames(gene_metadata), rownames(cell_metadata))
  } else {
    stop("Could not read matrix inside the anndata object.")
  }

  return(matrix_to_seurat(matrix, cell_metadata, gene_metadata))
}



#' Read 10x hdf5 dataset into seurat.
#' 
#' @param ds_file The path to the dataset file.
#' 
#' @return A seurat Object.
#' 
read_10xhdf5_to_seurat <- function(ds_file) {
  matrix <- Seurat::Read10X_h5(ds_file)
  seurat <- Seurat::CreateSeuratObject(counts = matrix, min.cells = 0, min.features = 0)
  return(seurat)
}


#' Read 10x mtx (mex) dataset to seurat.
#' 
#' @param ds_file The path to the dataset file.
#' 
#' @return A seurat Object.
#' 
read_10xmtx_to_seurat <- function(ds_file) {
  dir = dirname(ds_file)
  expression_matrix <- Seurat::Read10X(data.dir = dir)
  seurat <- Seurat::CreateSeuratObject(
  counts = expression_matrix, min.cells = 0, min.features = 0)

  return(seurat)
}




#' Read dense matrix in csv form
#' 
#' @param ds_file The path to the dataset file.
#' 
#' @return A seurat Object.
#' 
read_densecsv_to_seurat <- function(ds_file) {
  return(read_densemat_to_seurat(ds_file, ","))
}


#' Read dense matrix in tsv form
#' 
#' @param ds_file The path to the dataset file.
#' 
#' @return A seurat Object.
#' 
read_densetsv_to_seurat <- function(ds_file) {
  return(read_densemat_to_seurat(ds_file, "\t"))
}


#' Read dense matrix
#' here we need to unpack the dataset before reading it
#' 
#' @param ds_file The path to the dataset file.
#' 
#' @param sep The column seperator (e.g. tab or comma)
#' 
#' @return A seurat Object.
#' 
read_densemat_to_seurat <- function(ds_file, sep) {
  x <- data.table::fread(ds_file, sep = sep, header = F, skip = 1, na.strings = NULL)
  genes <- x[[1]]
  cells <- colnames(data.table::fread(ds_file, sep = sep, header = T, nrows = 0))
  cells <- utils::tail(cells, dim(x)[[2]] - 1) # ignore column name of genes, if present
  matrix <- as.matrix(x[, 2:dim(x)[2]])
  dimnames(matrix) <- list(genes, cells)
  spmatrix <- Matrix::Matrix(matrix, sparse = T)
  seurat <- Seurat::CreateSeuratObject(counts = spmatrix, min.cells = 0, min.features = 0)
  return(seurat)
}


EXPERIMENTAL_READERS <- list(
    "h5ad" = read_anndata_to_seurat_exp,
    "loom" = read_loom_to_seurat_exp
)

DEFAULT_READERS <- list(
    "loom" = read_loom_to_seurat,
    "rds" = read_seurat_to_seurat,
    "h5ad" = read_anndata_to_seurat,
    "hdf5" = read_10xhdf5_to_seurat,
    "h5" = read_10xhdf5_to_seurat,
    "tsv" = read_densetsv_to_seurat,
    "csv" = read_densecsv_to_seurat
)
