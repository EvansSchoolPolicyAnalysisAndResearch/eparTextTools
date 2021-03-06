% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TopicFunctions.R
\name{PreTopicFrame}
\alias{PreTopicFrame}
\title{Process a corpus into a topic model ready object}
\usage{
PreTopicFrame(CORPUS_A, howmanyentities = 10)
}
\arguments{
\item{CORPUS_A}{Document corpus}

\item{howmanyentities}{Count of entities you would like to add as potential covariates for topic model}
}
\value{
SentFrame data frame with one row for each paragraph chunk

Annotations openNLP annotation object

processed stm processed documents

out stm out documents for use in topic model
}
\description{
This function takes a corpus and creates a processed version of that corpus with entities removed for use in a topic model. Additionally it allows you to specify common entities to count across documents for use as a covariate in the topic model. The object it returns includes a frame of the text, an annotation object, a processed version of the corpus with stems and stopwords removed, and an out object which is the input object for fitting a topic model within the stm package.
}
\examples{
BASE_INPUT<-PreTopicFrame(corpus1,1)
}
\seealso{
\code{\link{stm}}
}
