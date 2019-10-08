ROOT_DIR <- Sys.getenv("FGROOT","/fastgenomics")
DATA_DIR <- file.path(ROOT_DIR, "data")
INFO_FILE_NAME <- "dataset_info.json"

#' This class stores the information about a dataset.
#'
#' @slot path Path to the dataset, typically of the form
#'     \code{"/fastgenomics/data/dataset_xxxx"}
#' @slot id Id of the dataset.  Defined as the \code{xxxx} part of the directory
#'     \code{".../dataset_xxxx"} under which the dataset is stored.)
#' @slot metadata Information extracted from \code{.../dataset_xxxx/dataset_info.json}
#' @slot file Absolute path to the file under which the dataset is stored, something
#'     like \code{".../dataset_xxxx/data.loom"}.
setClass("DataSet",
         slots = list(id="integer", metadata="list", path="character", file="character")
         )

#' constructor for the DataSet Object, takes dataset path and extracts all information
#' from it.
DataSet <- function(path){
    id <- as.integer(tail(strsplit(path, "_")[[1]], n=1))
    metadata <- jsonlite::read_json(file.path(path, INFO_FILE_NAME))
    return(new("DataSet",
               id=id,
               metadata=metadata,
               path=path,
               file=file.path(path, metadata$file)))
}

setMethod(
    "show", "DataSet",
    function(object){
        print(glue::glue(
                        "id: {object@id}\n",
                        "title: {object@metadata$title}\n",
                        "format: {object@metadata$format}\n",
                        "path: {object@path}"
                    ))
    }
)

#' Lists datasets provided by FASTGenomics.
#'
#' The optional argument available under \code{data_dir}.  The function then looks for
#' all directories matching the \code{".*/dataset_\\d{4}$"} pattern and constructs a
#' \code{DataSet} object from each of them.
#'
#' @param data_dir The directory containing sub-folders of the \code{"dataset_xxxx"}
#'     form.  Defaults to the directory provided by the FASTGenomics platform.
#'
#' @return List of datasets (as \code{DataSet} objects) indexed by dataset id.
#'
#' @examples
#' dsets_list <- get_datasets()
#' dsets_list[[1]]  # gives you the dataset with id = 1
#'
#' @export
get_datasets <- function(data_dir=DATA_DIR){
    dirs <- list.dirs(path=data_dir, full.names=T)
    dirs <- dirs[grepl(".*/dataset_\\d{4}$", dirs)]

    data_sets = list()
    for (dir in dirs){
        data_set <- DataSet(dir)
        data_sets[[data_set@id]] <- data_set
    }
    return(data_sets)
}

#' Prints the list of available datasets'
#'
#' @export
print_datasets <- function(data_dir=DATA_DIR){
    datasets <- get_datasets(data_dir=data_dir)

    print(datasets)
}

#' Adds some data from the metadata directly to the meta.data of the seurat object.
add_metadata <- function(seurat, data_set){
    seurat@project.name <- data_set@metadata$title
    seurat@meta.data$fg_dataset_id <- as.factor(data_set@id)
    seurat@meta.data$fg_dataset_title <- as.factor(data_set@metadata$title)
    seurat@misc$metacolumns <- c(seurat@misc$metacolumns, "fg_dataset_id", "fg_dataset_title")
    seurat@misc$fastgenomics = list(metadata=data_set@metadata, id=data_set@id)
    return(seurat)
}


#' Reads a single dataset.
#'
#' Takes a \code{\link{DataSet}} object as an argument.  This function should be used
#' with conjunction to \code{\link{get_datasets}}.
#'
#' @param data_set The dataset to load (passed as a DataSet object)
#'
#' @return dataset loaded as a Seurat Object
#'
#' @examples
#' dsets_list <- get_datasets()
#' read_dataset(dsets_list[[1]])  # returns the Seurat object constructed from the first dataset
#'
#' @export
read_dataset <- function(data_set, readers=DEFAULT_READERS){
    force(data_set)
    format <- data_set@metadata$format

    ## find a matching reader
    supported_readers_str <- paste(names(readers), collapse=", ")
    if(format == "Other"){
        stop(glue::glue('The format of the dataset "{data_set@metadata$title}" is "{format}".  datasets with the "{format}" format are unsupported by this module and have to be loaded manually.'))
    }
    else if(format == "Not set"){
        stop(glue::glue('The format of the dataset "{data_set@metadata$title}" is "{format}".  Please specify the data format in the Details of this dataset if you can modify the dataset or ask the dataset owner to do that.'))
    }
    else if(format %in% names(readers)){
        title <- data_set@metadata$title
        path <- data_set@path
        print(glue::glue('Loading dataset "{title}" in format "{format}" from directory "{path}"...'))
        seurat <- readers[[format]](data_set)

        ## Calling this function here provides compatibility between various readers,
        ## e.g. every seurat dataset will have @project.name coming from the manifest.
        ## On the downside, with custom readers this may lead to overwriting
        ## user-defined data in the seurat object.
        seurat <- add_metadata(seurat, data_set)
        n_genes <- dim(seurat)[[1]]
        n_cells <- dim(seurat)[[2]]
        print(glue::glue('Loaded dataset "{title}" with {n_genes} genes and {n_cells} cells'))
        return(seurat)
    } else {
        stop(glue::glue('Unsupported format: "{format}".'))
    }
}


#' Loads all datasets available for this analysis.
#'
#' Optionally, you can provide a list of datasets to load to load only selected ones.
#'
#' @param data_sets Optional argument with a list of the datasets to load.  Typically a
#'     list returned by \code{\link{get_datasets}} filtered to include only specific
#'     datasets.
#'
#' @return A list of datasets loaded into seurat objects, indexed by dataset id (same
#'     indices as in the \code{\link{get_datasets}}).
#'
#' @examples
#' # loads all datasets as seurat objects
#' seurat <- read_datasets()
#'
#' # loads only the first and second datasets
#' dsets_list <- get_datasets()
#' seurat <- read_datasets(dsets_list[c(1,2)])
#'
#' @export
read_datasets <- function(data_sets=get_datasets(DATA_DIR), readers=DEFAULT_READERS){
    loaded = list()
    for (dset in data_sets){
        print(glue::glue('Loading dataset "{dset@metadata$title}" in format {dset@metadata$format} from {dset@path}.'))
        loaded[[dset@id]] <- read_dataset(data_set=dset, readers=readers)
    }
    return(loaded)
}
