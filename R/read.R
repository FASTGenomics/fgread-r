ROOT_DIR <- Sys.getenv("FGROOT","/fastgenomics")
DATA_DIR <- file.path(ROOT_DIR, "data")

#' @export
list_datasets <- function(data_dir=DATA_DIR){
    data_sets <- list.dirs(path=data_dir, full.names=T)
    data_sets <- data_sets[grepl(".*/dataset_\\d{4}$", data_sets)]

    Map(
        function(path) {
            list(
                path=path,
                manifest=jsonlite::read_json(file.path(path, "manifest.json"))
            )
        },
        data_sets
    )
}


#' Reads a single data set from a given path
#'
#' @param path Location of the data set
#'
#' @export
read_dataset <- function(path, readers=DEFAULT_READERS){
    manifest_file <- file.path(path, "manifest.json")
    manifest <- jsonlite::read_json(manifest_file)
    format <- manifest$format

    if(format %in% names(readers)){
        return(readers[[format]](directory=path, manifest=manifest))
    } else {
        stop(stringr::str_interp("Unsupported format: ${format}. Use one of ${names(readers)}"))
    }
}


#'
#'
#'
#' @export
read_datasets <- function(data_dir=DATA_DIR, readers=DEFAULT_READERS){
    data_sets <- list_datasets(data_dir=data_dir)
    Map(
        function(dset){
            print(stringr::str_interp("Loading data set \"${dset$manifest$title}\" in format ${dset$manifest$format} from ${dset$path}."))
            read_dataset(path=dset$path, readers=readers)
        },
        data_sets
    )
}
