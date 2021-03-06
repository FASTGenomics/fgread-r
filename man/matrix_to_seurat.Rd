% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readers.R
\name{matrix_to_seurat}
\alias{matrix_to_seurat}
\title{a helper function that constructs a seurat object from an expression matrix, cell &
gene metadata.  Used to read loom and AnnData files}
\usage{
matrix_to_seurat(matrix, cell_metadata, gene_metadata)
}
\arguments{
\item{matrix}{The path to the expression matrix file.}

\item{cell_metadata}{The path to the cell metadata file.}

\item{gene_metadata}{The path to the gene metadata file.}
}
\value{
A Seurat object
}
\description{
a helper function that constructs a seurat object from an expression matrix, cell &
gene metadata.  Used to read loom and AnnData files
}
