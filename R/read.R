ROOT_DIR <- Sys.getenv("FGROOT","/fastgenomics")
DATA_DIR <- file.path(ROOT_DIR, "data")
INFO_FILE_NAME <- "dataset_info.json"

#' This class stores the information about a data set.
#'
#' @slot path Path to the data set, typically of the form
#'     \code{"/fastgenomics/data/dataset_xxxx"}
#' @slot id Id of the data set.  Defined as the \code{xxxx} part of the directory
#'     \code{".../dataset_xxxx"} under which the data set is stored.)
#' @slot metadata Information extracted from \code{.../dataset_xxxx/dataset_info.json}
#' @slot file Absolute path to the file under which the data set is stored, something
#'     like \code{".../dataset_xxxx/data.loom"}.
setClass("DataSet",
         slots = list(id="integer", metadata="list", path="character", file="character")
         )

#' constructor for the DataSet Object, takes data set path and extracts all information
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

#' Lists data sets provided by FASTGenomics.
#'
#' The optional argument available under \code{data_dir}.  The function then looks for
#' all directories matching the \code{".*/dataset_\\d{4}$"} pattern and constructs a
#' \code{DataSet} object from each of them.
#'
#' @param data_dir The directory containing sub-folders of the \code{"dataset_xxxx"}
#'     form.  Defaults to the directory provided by the FASTGenomics platform.
#'
#' @return List of data sets (as \code{DataSet} objects) indexed by data set id.
#'
#' @examples
#' dsets_list <- list_datasets()
#' dsets_list[[1]]  # gives you the data set with id = 1
#'
#' @export
list_datasets <- function(data_dir=DATA_DIR){
    dirs <- list.dirs(path=data_dir, full.names=T)
    dirs <- dirs[grepl(".*/dataset_\\d{4}$", dirs)]

    data_sets = list()
    for (dir in dirs){
        data_set <- DataSet(dir)
        data_sets[[data_set@id]] <- data_set
    }
    return(data_sets)
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


#' Reads a single data set.
#'
#' Takes a \code{\link{DataSet}} object as an argument.  This function should be used
#' with conjunction to \code{\link{list_datasets}}.
#'
#' @param data_set The data set to load (passed as a DataSet object)
#'
#' @return Data set loaded as a Seurat Object
#'
#' @examples
#' dsets_list <- list_datasets()
#' read_dataset(dsets_list[[1]])  # returns the Seurat object constructed from the first data set
#'
#' @export
read_dataset <- function(data_set, readers=DEFAULT_READERS){
    force(data_set)
    format <- data_set@metadata$format

    ## find a matching reader
    supported_readers_str <- paste(names(readers), collapse=", ")
    if(format == "Other"){
        stop(glue::glue('The format of the data set "{data_set@metadata$title}" is "{format}".  Data sets with the "{format}" format are unsupported by this module and have to be loaded manually.'))
    }
    else if(format == "Not set"){
        stop(glue::glue('The format of the data set "{data_set@metadata$title}" is "{format}".  Please specify the data format in the Details of this data set if you can modify the data set or ask the data set owner to do that.'))
    }
    else if(format %in% names(readers)){
        seurat <- readers[[format]](data_set)

        ## Calling this function here provides compatibility between various readers,
        ## e.g. every seurat data set will have @project.name coming from the manifest.
        ## On the downside, with custom readers this may lead to overwriting
        ## user-defined data in the seurat object.
        seurat <- add_metadata(seurat, data_set)
        return(seurat)
    } else {
        stop(glue::glue('Unsupported format: "{format}".'))
    }
}


#' Loads all data sets available for this analysis.
#'
#' Optionally, you can provide a list of data sets to load to load only selected ones.
#'
#' @param data_sets Optional argument with a list of the datasets to load.  Typically a
#'     list returned by \code{\link{list_datasets}} filtered to include only specific
#'     data sets.
#'
#' @return A list of data sets loaded into seurat objects, indexed by data set id (same
#'     indices as in the \code{\link{list_datasets}}).
#'
#' @examples
#' # loads all data sets as seurat objects
#' seurat <- read_datasets()
#'
#' # loads only the first and second data sets
#' dsets_list <- list_datasets()
#' seurat <- read_datasets(dsets_list[c(1,2)])
#'
#' @export
read_datasets <- function(data_sets=list_datasets(DATA_DIR), readers=DEFAULT_READERS){
    loaded = list()
    for (dset in data_sets){
        print(glue::glue('Loading data set "{dset@metadata$title}" in format {dset@metadata$format} from {dset@path}.'))
        loaded[[dset@id]] <- read_dataset(data_set=dset, readers=readers)
    }
    return(loaded)
}
