canonize_id <- function(metadata, regex, to_name){
    cols<- names(metadata)
    cell_id_names <- grep(regex, cols, value=T, ignore.case=T)
    if(length(cell_id_names) > 0){
        cell_id_name <- sort(cell_id_names)[1]
        names(metadata)[cols==cell_id_name] <- to_name
    } else {
        stop(glue::glue("Looking for {to_name}, could not find a column matching {regex}."))
    }
    return(metadata)
}

#' a helper function that constructs a seurat object from an expression matrix, cell &
#' gene metadata.  Used to read loom and AnnData files
matrix_to_seurat <- function(matrix, cell_metadata, gene_metadata){
    seurat <- Seurat::CreateSeuratObject(counts=matrix, min.cells=0, min.features=0)
    Seurat::AddMetaData(seurat, metadata=cell_metadata)
    seurat@assays$RNA@meta.features <- gene_metadata
    return(seurat)
}

#' Read a Seurat object.
read_seurat_to_seurat <- function(data_set){
    return(readRDS(data_set@file))
}

#' Reading a loom dataset currently not implemented in Seurat v3.
#' For your convenience we implemented experimental readers that you can use by setting "experimental_readers=TRUE"
#' in \code{\link{read_datasets}} or \code{\link{read_dataset}}.
read_loom_to_seurat <- function(data_set){
    stop(glue::glue(
        'Loading of loom files is currently not supported in Seurat v3. ',
        'You can use our FASTGenomics experimental reader by setting "experimental_readers=TRUE" in `read_datasets` or `read_dataset`. ',
        'For more information please see {BLOGURL}.'))
}

#' Importing loom files is currently unavailable in Seurat v3
#' This beta loader provided by our team only reads the count table
#' and the row/col attributes that can fit in a data frame structure (i.e. all higher
#' dimensional attributes are discarded). To keep the implementation simple we read the
#' whole object into memory, including the dense count matrix.  This could be a
#' potential bottleneck for larger datasets.
read_loom_to_seurat_exp <- function(data_set){
  
  warning("
          !! Importing loom files is currently unavailable in Seurat v3 !! 
          For your convenience the FASTGenomics team provides this beta loading routine. 
          In case of problems please consider using another format or implement your own loading routine.\n",
          call. = TRUE, immediate. = TRUE)
    
    file <- rhdf5::H5Fopen(data_set@file, flags="H5F_ACC_RDONLY")
    contents <- rhdf5::h5dump(file)
    rhdf5::H5Fclose(file)

    cleanup <- function(cols){
        Map(as.vector, Filter(function(x) length(dim(x))==1, cols))
    }

    gene_metadata <- data.frame(cleanup(contents$row_attrs))
    cell_metadata <- data.frame(cleanup(contents$col_attrs))

    ## this is necessary because the names for cell/gene ids changed over time for the
    ## loom format.  This way we make sure that the at least the cell/gene ids column
    ## has a consistent naming.  Both fields are also required to create a Seurat
    ## object.
    cell_metadata <- canonize_id(cell_metadata, regex="^(cell[^a-z]*(id|name)?$|name|obs_names)$", to_name="fg_cell_id")
    gene_metadata <- canonize_id(gene_metadata, regex="^(gene[^a-z]*(id|name)?$|name|var_names)$", to_name="fg_gene_id")

    rownames(gene_metadata) <- gene_metadata[["fg_gene_id"]]
    rownames(cell_metadata) <- cell_metadata[["fg_cell_id"]]

    matrix <- Matrix::Matrix(t(contents$matrix), sparse=T)
    dimnames(matrix) <- list(rownames(gene_metadata), rownames(cell_metadata))

    return(matrix_to_seurat(matrix, cell_metadata, gene_metadata))
}

#' Read AnnData to Seurat with Seurat's function \code{\link{ReadH5AD}}.
#' Pleas enote that this might not work as expected.
#' For your convenience we implemented experimental readers that you can use by setting "experimental_readers=TRUE"
#' in \code{\link{read_datasets}} or \code{\link{read_dataset}}.
read_anndata_to_seurat <- function(data_set){
    return(Seurat::ReadH5AD(data_set@file))
}


#' Read AnnData to Seurat.
#' Importing AnnData is not always working as expected in Seurat v3 
#' For your convenience the FASTGenomics team provides this beta loading routine.
#' Import of AnnData only works if there is a CSR matrix in the AnnData object.
read_anndata_to_seurat_exp <- function(data_set){
  
  warning("
          !! Importing AnnData is not generally available in Seurat v3 !! 
          Import of AnnData only works if there is a CSR matrix in the AnnData object.
          For your convenience the FASTGenomics team provides this beta loading routine.\n",
          call. = TRUE, immediate. = TRUE)
  
    file <- rhdf5::H5Fopen(data_set@file, flags="H5F_ACC_RDONLY")
    contents <- rhdf5::h5dump(file)
    rhdf5::H5Fclose(file)

    cell_metadata <- contents$obs
    gene_metadata <- contents$var

    matrix <- Matrix::sparseMatrix(
        i = contents$X$indices+1,
        p = contents$X$indptr,
        x = as.vector(contents$X$data),
        dims = c(dim(gene_metadata)[1], dim(cell_metadata)[1]),
        dimnames = list(rownames(gene_metadata), rownames(cell_metadata))
    )

    return(matrix_to_seurat(matrix, cell_metadata, gene_metadata))
}



#' Read 10x hdf5 dataset into seurat.
read_10xhdf5_to_seurat <- function(data_set){
    matrix <- Seurat::Read10X_h5(data_set@file)
    seurat <- Seurat::CreateSeuratObject(counts=matrix, min.cells=0, min.features=0)
    return(seurat)
}


#'read_10xhdf5_to_seuratad 10x mtx (mex) dataset to seurat.
read_10xmtx_to_seurat <- function(data_set){
  list_files <- list.files(data_set@path)
  suffix <- tail(list_files[[1]], 3)
  if(suffix == '.gz'){
    data <- Seurat::Read10X(data.dir = data_dir)
    seurat <- Seurat::CreateSeuratObject(
      counts = data$`Gene Expression`, min.cells=0, min.features=0)
  }
  else{
    expression_matrix <- Seurat::Read10X(data.dir = data_set@path)
    seurat <- Seurat::CreateSeuratObject(
      counts = expression_matrix, min.cells=0, min.features=0)  
  }
  return(seurat)
}


#' Read dense matrix in csv form
read_densecsv_to_seurat <- function(data_set){
    return(read_densemat_to_seurat(data_set, ","))
}


#' Read dense matrix in tsv form
read_densetsv_to_seurat <- function(data_set){
    return(read_densemat_to_seurat(data_set, "\t"))
}


#' Read dense matrix
#' here we need to unpack the dataset before reading it
read_densemat_to_seurat <- function(data_set, sep){
    file <- data_set@file
    x <- data.table::fread(file, sep=sep, header=F, skip=1, na.strings=NULL)
    genes <- x[[1]]
    cells <- colnames(data.table::fread(file, sep=sep, header=T, nrows=0))
    cells <- tail(cells, dim(x)[[2]]-1) # ignore column name of genes, if present
    matrix <- as.matrix(x[,2:dim(x)[2]])
    dimnames(matrix) <- list(genes, cells)
    spmatrix <- Matrix::Matrix(matrix, sparse=T)
    seurat <- Seurat::CreateSeuratObject(counts=spmatrix, min.cells=0, min.features=0)
    return(seurat)
}


DEFAULT_READERS <- list(
    "Loom"=read_loom_to_seurat,
    "Seurat Object"=read_seurat_to_seurat,
    "AnnData"=read_anndata_to_seurat,
    "10x (hdf5)"=read_10xhdf5_to_seurat,
    "10x (mtx)"=read_10xmtx_to_seurat,
    "tab-separated text"= read_densetsv_to_seurat,
    "comma-separated text"= read_densecsv_to_seurat
)

EXPERIMENTAL_READERS <- list(
    "AnnData"=read_anndata_to_seurat_exp,
    "Loom"=read_loom_to_seurat_exp
)
