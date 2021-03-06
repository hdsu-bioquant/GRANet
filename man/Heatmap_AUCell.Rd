% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{Heatmap_AUCell}
\alias{Heatmap_AUCell}
\title{Heatmap of cssRegulons activities}
\usage{
Heatmap_AUCell(
  GRANetObject,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_column_dend = TRUE,
  show_row_dend = TRUE,
  cluster_rows = TRUE
)
}
\arguments{
\item{GRANetObject}{GRANet object with calculated cssRegulon activities.
See function "cssRegulons_activity".}

\item{show_column_names}{Display cell IDs.}

\item{show_row_names}{Display transcription factor symbol.}

\item{show_column_dend}{Display column dendrogram.}

\item{show_row_dend}{Display row dendrogram}

\item{cluster_rows}{Cluster rows (cssRegulons).}
}
\value{

}
\description{
Heatmap of cssRegulons activities
}
