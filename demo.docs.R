#' This is data to be included in my package
#'
#' @name demo.docs
#' @docType data
#' @author Ryan Scott \email{ryscott@uw.edu}
#' @keywords data
NULL

demo.docs<-read.csv("demodoc_links.csv",stringsAsFactors=FALSE)
demo.docs<-demo.docs$links
if(file.exists("demo.docs.folder")==FALSE){dir.create("demo.docs.folder")}
demo.docs<-demo.docs[which(file.exists(paste("demo.docs.folder/",sapply(demo.docs, function(X){tail(unlist(strsplit(X,split="/",fixed=TRUE)),n=1)}),sep=""))==FALSE)]
lapply(demo.docs, function(X){download(X, destfile=file.path("demo.docs.folder",tail(unlist(strsplit(X,split="/",fixed=TRUE)),n=1)))})


