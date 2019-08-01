data_file_name <- "data"

canonize_id <- function(metadata, regex, to_name){
    cols<- names(metadata)
    cell_id_names <- grep(regex, cols, value=T, ignore.case=T)
    if(length(cell_id_names) > 0){
        cell_id_name <- sort(cell_id_names)[1]
        names(metadata)[cols==cell_id_name] <- to_name
    } else {
        stop(stringr::str_interp("Looking for ${to_name}, could not find a column matching ${regex}."))
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


read_seurat <- function(directory, manifest){
    file_name <- file.path(directory, paste(data_file_name, "rds", sep="."))
    return(readRDS(file_name))
}


#' Unfortunaly there is no Seurat function that straight up reads the loom file into a
#' Seurat object so we have to implement our own.  This one only reads the count table
#' and the row/col attributes that can fit in a data frame structure (i.e.. all higher
#' dimensional attributes are discarded).  To keep the implementation simple we read the
#' whole object into memory, including the dense count matrix.  This could be a
#' potential bottleneck for larger data sets but can be optimized later.
read_loom <- function(directory, manifest){
    file_name <- file.path(directory, paste(data_file_name, "loom", sep="."))
    file <- rhdf5::H5Fopen(file_name, flags="H5F_ACC_RDONLY")
    contents <- rhdf5::h5dump(file)

    cleanup <- function(cols){
        Map(as.vector, Filter(function(x) length(dim(x))==1, cols))
    }

    gene_metadata <- data.frame(cleanup(contents$row_attrs))
    cell_metadata <- data.frame(cleanup(contents$col_attrs))

    ## this is necessary because the names for cell/gene ids changed over time for the
    ## loom format.  This way we make sure that the at least the cell/gene ids column
    ## has a consistent naming.  Both fields are also required to create a Seurat
    ## object.
    cell_metadata <- canonize_id(cell_metadata, regex="^cell[^a-z]*id$", to_name="CellID")
    gene_metadata <- canonize_id(gene_metadata, regex="^(gene|name|var_names)$", to_name="Name")

    rownames(gene_metadata) <- gene_metadata[["Name"]]
    rownames(cell_metadata) <- cell_metadata[["CellID"]]

    matrix <- Matrix::Matrix(t(contents$matrix), sparse=T)
    dimnames(matrix) <- list(rownames(gene_metadata), rownames(cell_metadata))

    return(matrix_to_seurat(matrix, cell_metadata, gene_metadata))
}


#' this only works if there's a CSR matrix in the AnnData object.  We would be happy to
#' use the Seurats ReadH5AD function but it's broken.
read_anndata <- function(directory, manifest){
    file_name <- file.path(directory, paste(data_file_name, "h5ad", sep="."))
    file <- rhdf5::H5Fopen(file_name, flags="H5F_ACC_RDONLY")
    contents <- rhdf5::h5dump(file)

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

read_10x_hdf5 <- function(directory, manifest){
    file_name <- file.path(directory, paste(data_file_name, "h5", sep="."))
    matrix <- Seurat::Read10X_h5(file_name)
    seurat <- Seurat::CreateSeuratObject(counts=matrix, min.cells=0, min.features=0)
    return(seurat)
}

#' here we need to unpack the data set before reading it
read_dropseq <- function(directory, manifest){
    file_name <- file.path(directory, paste(data_file_name, "tsv", sep="."))
    x <- data.table::fread(file_name, sep="\t", header=F, skip=1, na.strings=NULL)
    genes <- x[[1]]
    cells <- colnames(data.table::fread(file_name, sep="\t", header=T, nrows=0))
    matrix <- as.matrix(x[,2:dim(x)[2]])
    dimnames(matrix) <- list(genes, cells)
    spmatrix <- Matrix::Matrix(matrix, sparse=T)
    seurat <- Seurat::CreateSeuratObject(counts=spmatrix, min.cells=0, min.features=0)
    return(seurat)
}

DEFAULT_READERS <- list(
    "Loom"=read_loom,
    "Seurat Object"=read_seurat,
    "AnnData"=read_anndata,
    "10x h5"=read_10x_hdf5,
    "Drop-Seq"=read_dropseq
)
