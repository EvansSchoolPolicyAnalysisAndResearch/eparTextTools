% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/textFunctions.R
\name{allDocs}
\alias{allDocs}
\title{Calls getTextR on all files in a directory joining into corpus.}
\usage{
allDocs(directory, SkiponError = FALSE, gen_pdf_tools = TRUE)
}
\arguments{
\item{directory}{Folder to read files from'}

\item{gen_pdf_tools}{if using a pc set to true. if on a unix server, setting to false will preserve pdf metadata.}

\item{onError}{if skip, skip documents that read wrong.}
}
\value{
Corpus of text documents.
}
\description{
This function will read from a folder documents of the class pdf, docx, doc or txt.
}
\examples{
allDocs("folder")
}
\seealso{
\code{\link{corpus}}
}
