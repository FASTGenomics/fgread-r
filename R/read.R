ROOT_DIR <- Sys.getenv("FGROOT","/fastgenomics")
DATA_DIR <- file.path(ROOT_DIR, "data")
INFO_FILE_NAME <- "dataset_info.json"

setClass("DataSet",
         slots = list(id="integer", metadata="list", path="character", file="character")
         )

## constructor for DataSet
DataSet <- function(path){
    id <- as.integer(tail(strsplit(path, "_")[[1]], n=1))
    metadata <- jsonlite::read_json(file.path(path, INFO_FILE_NAME))
    return(new("DataSet",
               id=id,
               metadata=metadata,
               path=path,
               file=file.path(path, metadata$file)))
}

#' Lists data sets
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

#' adds some data from the metadata directly to the meta.data of the seurat object.
add_metadata <- function(seurat, data_set){
    seurat@project.name <- data_set@metadata$title
    seurat@meta.data$fg_dataset_id <- as.factor(data_set@id)
    seurat@meta.data$fg_dataset_title <- as.factor(data_set@metadata$title)
    seurat@misc$metacolumns <- c(seurat@misc$metacolumns, "fg_dataset_id", "fg_dataset_title")
    seurat@misc$fg_metadata <- data_set@metadata
    return(seurat)
}


#' Reads a single data set
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
        seurat <- add_metadata(seurat, data_set)
        return(seurat)
    } else {
        stop(glue::glue('Unsupported format: "{format}".'))
    }
}


#' Loads multiple data sets
#'
#'
#' @export
read_datasets <- function(data_sets=list_datasets(data_dir), readers=DEFAULT_READERS){
    Map(
        function(dset){
            print(glue::glue('Loading data set "{dset@metadata$title}" in format {dset@metadata$format} from {dset@path}.'))
            read_dataset(data_set=dset, readers=readers)
        },
        data_sets
    )
}
