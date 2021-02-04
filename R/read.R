ROOT_DIR <- Sys.getenv("FGROOT", "/fastgenomics")
DATA_DIR <- file.path(ROOT_DIR, "data")
INFO_FILE_NAME <- "dataset_info.json"

# Set blog url for more information
if (Sys.getenv("FG_URL") != "") {
  ENV_FGURL <- Sys.getenv("FG_URL")
  FGURL <- paste(strsplit(ENV_FGURL, ":")[[1]][1], strsplit(ENV_FGURL, ":")[[1]][2], sep = ":")
} else {
  FGURL <- "https://beta.fastgenomics.org"
}

DOCSURL <- paste(FGURL, "/docs/", sep = "")
DS_URL_PREFIX <- paste(FGURL, "/webclient/ui/#/datasets/detail-", sep = "")

#' Construct a dataframe of all availabe datasets
#' 
#' @param data_dir The directory containing sub-folders of the \code{"dataset_xxxx"}
#'     form.  Defaults to the directory provided by the FASTGenomics platform.
#' 
#' @param ignore_empty Whether to ignore if no datasets are attached or not.
#'      Defaults to TRUE.
#' 
#' @return A dataframe containing all dataset metadata
#' 
get_datasets_df <- function(data_dir = DATA_DIR, ignore_empty = TRUE) {
  # get all data set folders
  dirs <- get_ds_paths(data_dir = data_dir, ignore_empty = ignore_empty)

  # create data frame with all data set informations
  ds_list = list()
  for (dir in dirs) {
    ds_path <- file.path(dir, INFO_FILE_NAME)
    ds_info <- jsonlite::read_json(ds_path)
    ds_info["path"] <- dir
    ds_info["numberOfExpressionDataFiles"] <- length(ds_info["expressionDataFileInfos"][[1]])
    ds_info["numberOfMetaDataFiles"] <- length(ds_info["metaDataFileInfos"][[1]])
    ds_info["schemaVersion"] <- NULL
    ds_info["expressionDataFileNames"] <- ""
    ds_info["metaDataFileNames"] <- ""
    ds_list <- rbind(ds_list, ds_info)
  }
  ds_df <- data.frame(ds_list, row.names = seq_along(dirs))

  # sort colnames
  sort_order <- c(
        "title",
        "id",
        "organism",
        "tissue",
        "numberOfCells",
        "numberOfGenes",
        "path",
        "numberOfExpressionDataFiles",
        "expressionDataFileNames",
        "numberOfMetaDataFiles",
        "metaDataFileNames",
        "expressionDataFileInfos",
        "metaDataFileInfos"
  )
  col_names_sorted <- c(sort_order, sort(setdiff(colnames(ds_df), sort_order)))
  # construct empty dataframe if no datasets attached
  if (length(dirs) == 0) {
    ds_df <- stats::setNames(data.frame(matrix(ncol = length(col_names_sorted), nrow = 0)), col_names_sorted)
  } else {
    ds_df <- ds_df[col_names_sorted]
  }

  return(ds_df)
}


#' Get information on all available datasets in this analysis.
#'
#' The optional argument available under \code{data_dir}.  The function then looks for
#' all directories matching the \code{".*/dataset_\\d{4}$"} pattern and constructs a
#' data frame with all data sets.
#'
#' @param ds A single dataset ID or dataset title. If set, only this dataset will be displayed.
#'
#' @param pretty Boolean whether to display some nicely formatted output, by default True
#'
#' @param output Boolean whether to return a DataFrame or not, by default True
#'
#' @param data_dir The directory containing sub-folders of the \code{"dataset_xxxx"}
#'     form.  Defaults to the directory provided by the FASTGenomics platform.
#' 
#' @param ignore_empty Whether to ignore if no datasets are attached or not.
#'      Defaults to TRUE.
#'
#' @return A data frame containing all, or a single dataset (depends on ``ds`` and ``output``)
#'
#' @export
ds_info <- function(ds = NULL, pretty = NULL, output = NULL, data_dir = DATA_DIR, ignore_empty = TRUE) {

  if (is.null(pretty)) {
    pretty = !is.null(ds)
  }
  if (is.null(output)) {
    output = is.null(ds)
  }

  if (!pretty & !output) {
    warning('You have set "pretty" and "output" to false. Hence, this function will do/return nothing.', immediate = T)
  }

  # create output and display 
  if (!missing(ds)) {
    # if ds is specified
    ds_df <- get_datasets_df(data_dir = data_dir, ignore_empty = FALSE)
    single_ds_df = select_ds_id(ds, ds_df)

    single_ds_df["expressionDataFileNames"] <- expr_file_infos_to_list(single_ds_df, "expressionDataFileInfos", as_string = T)
    single_ds_df["metaDataFileNames"] <- expr_file_infos_to_list(single_ds_df, "metaDataFileInfos", as_string = T)

    if (pretty) {
      pretty_df <- single_ds_df
      pretty_df$title <- paste0("<a href='", DS_URL_PREFIX, pretty_df$id, "' target='_blank'>", pretty_df$title, "</a>")

      pretty_df["expressionDataFileNames"] <- gsub(", ", "<br>", single_ds_df["expressionDataFileNames"])
      pretty_df["metaDataFileNames"] <- gsub(", ", "<br>", single_ds_df["metaDataFileNames"])

      drop = c("expressionDataFileInfos", "metaDataFileInfos")
      pretty_df = pretty_df[, !(names(pretty_df) %in% drop)]
      pretty_df = pretty_df[!sapply(pretty_df, function(x) is.null(x[[1]]))]
      pretty_df = pretty_df[!sapply(pretty_df, function(x) x[[1]] == "")]

      dt <- htmlTable::htmlTable(t(data.frame(pretty_df, row.names = NULL)), align = "l")
      IRdisplay::display_html(dt)
    }

    if (output) {
      return(single_ds_df)
    }

  } else {
    ds_df <- get_datasets_df(data_dir = data_dir, ignore_empty = ignore_empty)
    drop = c("expressionDataFileNames",
             "metaDataFileNames"
    )
    ds_df = ds_df[, !(names(ds_df) %in% drop)]

    if (pretty) {
      drop = c("description",
               "license",
               "preprocessing",
               "citation",
               "webLink",
               "file",
               "expressionDataFileInfos",
               "metaDataFileInfos"
      )
      df = ds_df[, !(names(ds_df) %in% drop)]

      if (dim(df)[1] > 0) {
        df$title <- paste0("<a href='", DS_URL_PREFIX, df$id, "' target='_blank'>", df$title, "</a>")
      }
      dt <- htmlTable::htmlTable(df)
      IRdisplay::display_html(dt)
    }

    if (output) {
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
select_ds_id <- function(ds, df) {
  single_df = df[df$title == ds | df$id == ds,]
  len_df <- dim(single_df)[1]
  if (len_df == 1) {
    return(single_df)
  } else if (len_df == 0) {
    add_err <- ""
    if (!startsWith(ds, "dataset-")) {
      add_err <- " Please note that dataset titles can be changed by the owner. To be safe, you might want to consider dataset IDs instead."
    }
    stop(paste("Your selection matches no datasets.", add_err, sep = ""))
  } else {
    stop(glue::glue("Your selection matches {len_df} datasets. Please make sure to select exactly one"))
  }
}


#' Loads a single dataset to a Seurat object.
#'
#' If there are multiple datasets available you need to specify one by setting
#' \code{ds} to a dataset \code{id} or dataset \code{title}.
#' To get an overview of availabe dataset use the fucntion \code{ds_info}.
#'
#' @param ds A single dataset ID or dataset title.
#'
#' @param data_dir The directory containing sub-folders of the \code{"dataset_xxxx"}
#'     form. Defaults to the directory provided by the FASTGenomics platform. 
#'
#' @param additional_readers List to specify your own readers for the specific dataset 
#'      format. Still experimental, by default an empty list
#'
#' @param experimental_readers Boolean whether FASTGenomics in-house experimental
#'       readers should be used. By default set to False. 
#' @param expression_file The name of the expression file you want to load.
#' Only needed when there are multiple expression files in a dataset. 
#' 
#' @param as_format Specifies which reader should be uses for this dataset. Overwrites the
#' auto-detection of the format. Possible parameters are the file extensions of our
#' supported data formats: "h5ad", "h5", "hdf5", "loom", "rds", "csv", "tsv".
#'
#' @return dataset loaded as a Seurat Object
#'
#'
#' @export
load_data <- function(ds, data_dir = DATA_DIR, additional_readers = list(), experimental_readers = F, expression_file, as_format) {
  # update list of readers
  if (experimental_readers) {
    readers <- utils::modifyList(DEFAULT_READERS, EXPERIMENTAL_READERS)
    readers <- utils::modifyList(readers, additional_readers)
  } else {
    readers <- utils::modifyList(DEFAULT_READERS, additional_readers)
  }

  # get single dataset
  if (missing(ds)) {
    single_df <- ds_info(data_dir = data_dir, pretty = FALSE, ignore_empty = F)

    if (dim(single_df)[1] != 1) {
      stop("There is more than one dataset available. Please select one by its ID or title.")
    }
  } else {
    single_df <- select_ds_id(ds, df = ds_info(data_dir = data_dir, pretty = FALSE, ignore_empty = F))
  }

  exp_count <- single_df$numberOfExpressionDataFiles[[1]]
  meta_count <- single_df$numberOfMetaDataFiles[[1]]

  if (exp_count == 0) {
    stop("There is no expression data available in this data set.")
  } else if (exp_count > 1) {
    expr_names <- expr_file_infos_to_list(single_df, "expressionDataFileInfos")
    str_expr_names <- paste(expr_names, collapse = ", ")

    if (missing(expression_file)) {
      stop(glue::glue("There is more than one expression data available in this data set.\n",
        'Please specifiy which you want to load by setting "expression_file".\n',
        "Available files are: {str_expr_names}."))
    } else {
      if (is.element(expression_file, expr_names)) {
        file <- expression_file
      } else {
        stop(glue::glue('File "{expression_file}" not found in dataset expresison files ({str_expr_names}).'))
      }

    }
  } else {
    file <- single_df["expressionDataFileInfos"][[1]][[1]][[1]]$name

    if (!missing(expression_file)) {
      if (expression_file != file) {
        stop(glue::glue('Expression file "{expression_file}", not found in this dataset.\n',
        "There is only one expression file ({file}) in this dataset.\n",
        'You can skip "expression_file".'))
      }
    }
  }


  title <- single_df$title[[1]]
  path <- single_df$path[[1]]

  if (missing(as_format)) {
    tryCatch({ format <- tolower(tools::file_ext(file)) },
    error = function(e) stop(glue::glue('The expression file "{file}" has no suffix.'))
    )
  } else {
    format <- tolower(as_format)
  }

  ## find a matching reader
  supported_readers_str <- paste(names(readers), collapse = ", ")
  if (format %in% names(readers)) {
    if (meta_count >= 1) {
      if (meta_count == 1) {
        print(glue::glue("There is {meta_count} metadata file in this dataset.\n"))
      } else {
        print(glue::glue("There are {meta_count} metadata files in this dataset.\n"))
      }
      print(glue::glue("This metadata will not be integrated into the anndata object."))
    }
    print(glue::glue('Loading file "{file}" from dataset "{title}" in format "{format}" from directory "{path}"...'))
    file_path = file.path(path, file)
    seurat <- readers[[format]](file_path)

    ## Calling this function here provides compatibility between various readers,
    ## e.g. every seurat dataset will have @project.name coming from the manifest.
    ## On the downside, with custom readers this may lead to overwriting
    ## user-defined data in the seurat object.
    seurat <- add_metadata(seurat, single_df)
    n_genes <- dim(seurat)[[1]]
    n_cells <- dim(seurat)[[2]]
    print(glue::glue('Loaded dataset "{title}" with {n_genes} genes and {n_cells} cells\n\n'))
    return(seurat)
  } else {
    stop(glue::glue(
            'Unsupported file format "{format}", use one of {supported_readers_str} or implement your ',
            "own reading function. See {DOCSURL} for more information."))
  }
}

#' Adds some data from the metadata directly to the meta.data of the seurat object.
#' 
#' @param seurat A seurat object.
#' 
#' @param ds_df The dataframe containing the metadata for the Seurat object.
#' 
#' @return A Seurat object with added metadata.
#' 
add_metadata <- function(seurat, ds_df) {

  title = ds_df$title[[1]]
  ds_id = ds_df$id[[1]]
  format = ds_df$format[[1]]
  path = ds_df$path[[1]]
  file = ds_df$file[[1]]
  metadata <- jsonlite::read_json(file.path(path, INFO_FILE_NAME))

  seurat@project.name <- ds_df$title[[1]]
  seurat@meta.data$fg_dataset_id <- as.factor(ds_df$id[[1]])
  seurat@misc$fastgenomics = list(metadata = metadata, id = ds_df$id[[1]])

  return(seurat)
}

#' Returns DataFileInfos as a vecotor or as a string.
#' 
#' @param single_ds_df A dataframe with a single dataset.
#' 
#' @param key The column to convert (e.g. "expressionDataFileInfos").
#' 
#' @param as_string Whether to return a formatted string.
#' 
#' @return A vecotr or a formatted string.
#' 
expr_file_infos_to_list <- function(single_ds_df, key, as_string = F) {
  z <- NULL
  for (expr in single_ds_df[key][[1]][[1]]) {
    z <- c(z, expr$name)
  }
  if (as_string) {
    s <- paste(z, collapse = ", ")
    return(s)
  } else {
    return(z)
  }
}


#' Get paths of all available datasets
#' 
#' @param data_dir The directory containing sub-folders of the \code{"dataset_xxxx"}
#'     form. Defaults to the directory provided by the FASTGenomics platform.
#' 
#' @param ignore_empty Whether to ignore if no datasets are attached or not.
#'     Defaults to TRUE.
#' 
#' @return A list of dataset directories
#' 
get_ds_paths <- function(data_dir, ignore_empty = T) {
  dirs <- list.dirs(path = data_dir, full.names = T, recursive = F)
  dirs <- dirs[grepl(".*/dataset_\\d{4}$", dirs)]

  if (length(dirs) == 0) {
    if (ignore_empty) {
      warning("There are no datasets attached to this analysis.")
    } else {
      stop("There are no datasets attached to this analysis.")
    }
  }
  return(dirs)
}
