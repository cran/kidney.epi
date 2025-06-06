% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix.R
\name{matrix.save_named_matrix_to_excel}
\alias{matrix.save_named_matrix_to_excel}
\title{Save a named matrix as an Excel file.}
\usage{
matrix.save_named_matrix_to_excel(
  matrix_data,
  file_path,
  sheet_name = "Sheet1",
  save_type = "new"
)
}
\arguments{
\item{matrix_data}{matrix. Matrix for saving.}

\item{file_path}{character. Path to the Excel file.}

\item{sheet_name}{character. Name of the Excel sheet.}

\item{save_type}{character. Defines whether the Excel file should be created or overwritten (save_type = "new"), or new sheet should be added to the existing Excel file (save_type = "add").}
}
\value{
Excel file saved to a specified folder.
}
\description{
Save a named matrix as an Excel file.
}
\details{
Save a named matrix as an Excel file.
}
