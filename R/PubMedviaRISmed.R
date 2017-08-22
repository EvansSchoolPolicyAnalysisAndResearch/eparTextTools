

#' Return structured dataframe based on PMID info.
#'
#' @return a dataframe of info on requested publications
#' @param pmids character vector of pubmed ids
#' @export
#' @description  This function extracts pubmed data and returns a dataframe based on a pubmed id. More variables could be added by adding additional names and commands from the RISmed package
#' @examples
#' example_documents()
searchPMIDs<-function(pmids){
  pm1<-RISmed::EUtilsGet(pmids,type="efetch",db="pubmed")
  toframe<-data.frame("Abstact"=RISmed::AbstractText(pm1),"Title"=RISmed::Title(pm1),"Affiliation"=RISmed::Affiliation(pm1),"YearAccepted"=RISmed::YearAccepted(pm1),"Agency"=RISmed::Agency(pm1),"id"=RISmed::ArticleId(pm1),"Author"=sapply(RISmed::Author(pm1),function(X) paste(X$LastName[X$order],X$Initials[X$order],sep=" ",collapse=";")),"ISSN"=RISmed::ISSN(pm1),"Acronym"=RISmed::Acronym(pm1),stringsAsFactors=FALSE)
}


#' Return structured dataframe based on search terms.
#'
#' @return a dataframe of info on requested publications
#' @export
#' @param searchterms a query to search
#' @email users email address
#' @retmax how many arguments to return
#' @... all other arguments passed to EUtilsQuery
#' @description  This function parses a search and returns a dataframe based on a pubmed id. More variables could be added by adding additional names and commands from the RISmed package
#' @examples
#'test1<-searchWordsPubMED("'gender equity'+environment","ryscott5@uw.edu",retmax=10)
#'head(test1)
searchWordsPubMED<-function(searchterms,email,retmax=25){
  qout<-EUtilsQuery(searchterms,type="esearch",db="pubmed",retmax=retmax) %>% gsub("s.a.kovalchik@gmail.com",email,.,fixed=T,...)
  temp<-httr::GET(qout) %>% httr::content(.,"text") %>% xml2::as_xml_document(.) %>% xml2::as_list()
  t2<-unlist(temp$IdList)
  pm1<-RISmed::EUtilsGet(t2,type="efetch",db="pubmed")
  toframe<-data.frame("Abstract"=RISmed::AbstractText(pm1),"Title"=RISmed::Title(pm1),"Affiliation"=RISmed::Affiliation(pm1),"YearAccepted"=RISmed::YearAccepted(pm1),"Agency"=RISmed::Agency(pm1),"id"=RISmed::ArticleId(pm1),"Author"=sapply(RISmed::Author(pm1),function(X) paste(X$LastName[X$order],X$Initials[X$order],sep=" ",collapse=";")),"ISSN"=RISmed::ISSN(pm1),"Acronym"=RISmed::Acronym(pm1),stringsAsFactors=FALSE)
  toframe
}


