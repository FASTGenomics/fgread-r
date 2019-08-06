library(glue)

ROOT_DIR <- Sys.getenv("FGROOT","/fastgenomics")
DATA_DIR <- file.path(ROOT_DIR, "data")

setClass("DataSet",
         slots = list(id="integer", manifest="list", path="character")
         )

## constructor for DataSet
DataSet <- function(path){
    id <- as.integer(tail(strsplit(path, "_")[[1]], n=1))
    manifest <- jsonlite::read_json(file.path(path, "manifest.json"))
    return(new("DataSet", id=id, manifest=manifest, path=path))
}

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

#' adds some data from the manifest directly to the meta.data of the seurat object.
add_metadata <- function(seurat, data_set){
    seurat@meta.data$fg_dataset_id <- as.factor(data_set@id)
    seurat@meta.data$fg_dataset_title <- as.factor(data_set@manifest$title)
    seurat@misc$metacolumns <- c(seurat@misc$metacolumns, "fg_dataset_id", "fg_dataset_title")
    return(seurat)
}


#' Reads a single data set
#'
#' @export
read_dataset <- function(data_set, readers=DEFAULT_READERS){
    force(data_set)
    format <- data_set@manifest$format

    ## find a matching reader
    if(format %in% names(readers)){
        seurat <- readers[[format]](data_set)
        seurat <- add_metadata(seurat, data_set)
        return(seurat)

    } else {
        stop(glue("Unsupported format: {format}. Use one of {names(readers)}"))
    }
}


#' Loads multiple data sets
#'
#'
#' @export
read_datasets <- function(data_sets=list_datasets(data_dir), readers=DEFAULT_READERS){
    Map(
        function(dset){
            print(glue("Loading data set \"${dset@manifest$title}\" in format ${dset@manifest$format} from ${dset@path}."))
            read_dataset(data_set=dset, readers=readers)
        },
        data_sets
    )
}
