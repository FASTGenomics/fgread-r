ROOT_DIR <- Sys.getenv("FGROOT", "/fastgenomics")
DATA_DIR <- file.path(ROOT_DIR, "data")
INFO_FILE_NAME <- "dataset_info.json"

# Set blog url for more information
BLOGURL <- "https://www.fastgenomics.org/blog_posts/readers/"


#' Get information on all available datasets in this analysis.
#'
#' The optional argument available under \code{data_dir}.  The function then looks for
#' all directories matching the \code{".*/dataset_\\d{4}$"} pattern and constructs a
#' data frame with all data sets.
#'
#' @param ds A single dataset ID or dataset title. If set, only this dataset will be displayed.
#'
#' @param output Boolean whether to return a DataFrame or not, by default True
#'
#' @param data_dir The directory containing sub-folders of the \code{"dataset_xxxx"}
#'     form.  Defaults to the directory provided by the FASTGenomics platform.
#'
#' @return A data frame containing all, or a single dataset (depends on ``ds`` and ``output``)
#'
#' @examples
#' dsets <- ds_info()
#' dsets <- ds_info('Test loom data')
#'
#' @export
ds_info <- function(ds, output = TRUE, data_dir = DATA_DIR) {
  dirs <- list.dirs(path = data_dir, full.names = T)
  dirs <- dirs[grepl(".*/dataset_\\d{4}$", dirs)]

  ds_list = list()
  for (dir in dirs) {
    ds_path <- file.path(dir, INFO_FILE_NAME)
    ds_info <- jsonlite::read_json(ds_path)
    ds_info["path"] <- ds_path
    ds_info["schemaVersion"] <- NULL
    ds_list <- rbind(ds_list, ds_info)
  }
  ds_df <- data.frame(ds_list, row.names=seq_along(dirs))

  # sort colnames
  sort_order <- c("title", "id", "format", "organism", "tissue", "numberOfCells","numberOfGenes")
  col_names_sorted <- c(sort_order, sort(setdiff(colnames(ds_df), sort_order)))
  ds_df <- ds_df[col_names_sorted]

  # TODO: check display pretty in notebook
  # TODO: display data frame independent of output variable
  
  if (!missing(ds)){
    if (output){
      return(select_ds_id(ds, ds_df))
    }
  } else {
    if (output){
      return(ds_df)
    }
  }
}

#' Select a single dataset from a pandas DataFrame by its ID or title
#'
#' @param ds A single dataset ID or dataset title for selection
#'
#' @param df A data frame from which a single entry is selected, by default None
#'
#' @return A data frame with only the selected dataset
select_ds_id <- function(ds, df){
  single_df = df[df$title == ds | df$id == ds,]
  len_df <- dim(single_df)[1]

  if (len_df == 1){
    return(single_df)
  } else {
    stop(glue::glue("Your selection matches {len_df} datasets. Please make sure to select exactly one"))
  }
}


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
         slots = list(id = "integer", metadata = "list", path = "character", file = "character")
         )

#' constructor for the DataSet Object, takes dataset path and extracts all information
#' from it.
DataSet <- function(path) {
  id <- as.integer(tail(strsplit(path, "_")[[1]], n = 1))
  metadata <- jsonlite::read_json(file.path(path, INFO_FILE_NAME))
  return(new("DataSet",
               id = id,
               metadata = metadata,
               path = path,
               file = file.path(path, metadata$file)))
}

setMethod(
    "show", "DataSet",
    function(object) {
      print(glue::glue(
                        "id: {object@id}\n",
                        "title: {object@metadata$title}\n",
                        "format: {object@metadata$format}\n",
                        "path: {object@path}\n",
                        "file: {object@file}"
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
get_datasets <- function(data_dir = DATA_DIR) {
  dirs <- list.dirs(path = data_dir, full.names = T)
  dirs <- dirs[grepl(".*/dataset_\\d{4}$", dirs)]

  data_sets = list()
  for (dir in dirs) {
    data_set <- DataSet(dir)
    data_sets[[data_set@id]] <- data_set
  }
  return(data_sets)
}


#' Adds some data from the metadata directly to the meta.data of the seurat object.
add_metadata <- function(seurat, data_set) {
  seurat@project.name <- data_set@metadata$title
  seurat@meta.data$fg_dataset_id <- as.factor(data_set@id)
  seurat@meta.data$fg_dataset_title <- as.factor(data_set@metadata$title)
  seurat@misc$metacolumns <- c(seurat@misc$metacolumns, "fg_dataset_id", "fg_dataset_title")
  seurat@misc$fastgenomics = list(metadata = data_set@metadata, id = data_set@id)
  return(seurat)
}


#' Reads a single dataset.
#'
#' Takes a \code{\link{DataSet}} object as an argument. This function should be used
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
read_dataset <- function(data_set, additional_readers = list(), experimental_readers = F) {
  force(data_set)
  format <- data_set@metadata$format

  if (experimental_readers) {
    readers <- utils::modifyList(DEFAULT_READERS, EXPERIMENTAL_READERS)
    readers <- utils::modifyList(readers, additional_readers)
  }
  else {
    readers <- utils::modifyList(DEFAULT_READERS, additional_readers)
  }

  ## find a matching reader
  supported_readers_str <- paste(names(readers), collapse = ", ")
  if (format %in% names(readers)) {
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
    print(glue::glue('Loaded dataset "{title}" with {n_genes} genes and {n_cells} cells\n\n'))
    return(seurat)
  }
  else if (format == "Other") {
    stop(glue::glue(
            'The format of the dataset "{data_set@metadata$title}" is "{format}". ',
            'Datasets with the "{format}" format are unsupported by this module and have to be loaded manually. ',
            'For more information please see {BLOGURL}.'))
  }
  else if (format == "Not set") {
    stop(glue::glue(
            'The format of the dataset "{data_set@metadata$title}" is "{format}". ',
            'Please specify the data format in the details of this dataset if you can modify the dataset or ask the dataset owner to do that. ',
            'For more information please see {BLOGURL}.'))
  }
  else {
    stop(glue::glue(
            'Unsupported format: "{format}". ',
            'For more information please see {BLOGURL}.'))
  }
}


#' Loads all datasets available for this analysis.
#'
#' Optionally, you can provide a list of datasets or a single dataset to load to load only selected ones.
#'
#' @param data_sets Optional argument with a list of the datasets to load.  Typically a
#'     list returned by \code{\link{get_datasets}} filtered to include only specific
#'     datasets. You can also pass a single dataset.
#'
#' @return A list of datasets loaded into seurat objects, indexed by dataset id (same
#'     indices as in the \code{\link{get_datasets}}). When you pass a single dataset to the function
#'     a single Seurat object is returned.
#'
#' @examples
#' # loads all datasets as seurat objects
#' seurat <- read_datasets()
#'
#' # loads only the first and second datasets
#' dsets_list <- get_datasets()
#' seurat <- read_datasets(dsets_list[c(1,2)])
#' 
#' # loads only the firt dataset and returns a Seurat object
#' seurat <- read_datasets(dsets_list[[1]])
#'
#' @export
read_datasets <- function(data_sets = get_datasets(DATA_DIR), additional_readers = list(), experimental_readers = F) {
  if (typeof(data_sets) == "S4") {
    return(read_dataset(data_sets, additional_readers = additional_readers, experimental_readers = experimental_readers))
  }
  else if (typeof(data_sets) == "list") {
    loaded = list()
    for (dset in data_sets) {
      loaded[[dset@id]] <- read_dataset(data_set = dset, additional_readers = additional_readers, experimental_readers = experimental_readers)
    }
    return(loaded)
  }
  else {
    stop(glue::glue(
            'You have to pass a dataset list or a single dataset. ',
            'For more information please see {BLOGURL}.'))
  }

}
