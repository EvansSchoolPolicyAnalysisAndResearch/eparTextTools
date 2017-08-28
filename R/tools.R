#==================
#FindHiddenHTML
#==================
#We had a large set of HTML files that were all txt files.  This converts them to proper txt files. This parses the files and removes all HTML properly. Doing so by regex is much faster but prone to missing html bits and creating bugs.
#' @param path A directory containing .txt files with some HTML in them we want removed
#' @param moveOldFiles moves the source files to a newsource "HTML_Sources" directory *alongside* the folder designated by the path. This folder holds all the text files.  If this moveOldFiles is set to FALSE it will not move the files and they will be left where they are.
#' @return A folder location containing the finished files
#' @seealso \code{\link{tm.plugin.webmining::extractHTMLStrip}} 
#' @export
#' @examples
#' outputFileLocation <- findHiddenHTML("C:/Users/Graham/Documents/EPAR/HTML_TXT_TestFiles")

findHiddenHTML<-function(path, moveOldFiles=TRUE){

#library 
library(tm.plugin.webmining)
  
#save the working directory
oldWD <- getwd()
#change working directory
setwd(path)
#build list of all PDFs at "path" and put that list in "listOfPDFs"
listOfFiles <- list.files(path,"\\.txt", recursive = TRUE, include.dirs = FALSE)
#create the output folder

outputPath <- toString(path)
setwd("..")
storeOriginalsPath <- toString(paste(getwd(),"/HTML_Sources/", sep=""))
setwd(path) #put the working directory back after creating the sources path.

dir.create(storeOriginalsPath)
outputCounter <- 1
for(eachFilePath in listOfFiles){
  fileName<-eachFilePath
  if(.Platform$OS.type == "windows"){
    fileName<-strsplit(eachFilePath, "/") #if we have a full path, cut it into chunks based on the /
    fileName<-fileName[length(fileName)] #set the filename to the last chunk of the path
    fileName <- paste0(fileName,collapse="_")
  }
  else{print("Warning: OS is not windows, so remember that the findHiddenHTML function must be called with filenames and not paths on other OS's.")}
  
  print(paste("Beginning with file ", toString(outputCounter), " out of ", length(listOfFiles), ", named ", fileName, sep=""))
  outputCounter <- outputCounter + 1
  
  fileName <- sapply(fileName,tolower)
  outputFileName <- gsub(".txt", "", fileName)  
  #prep a file connection to output the raw text
  outputFile <- file(description=paste("noHTML_",toString(outputFileName),".txt",sep=""), open="a", blocking=FALSE, encoding="", raw=FALSE, method="internal")
  
  #open the file to read
  fileConnection = file(eachFilePath, "r")
  
  #readline
#  while(TRUE){
#    line = readLines(fileConnection, n = 1)
#    if ( length(line) == 0 ) {
#      break
#    }
#  
#    #strip HTML tags
#    #outputLine <-gsub ("<.*?>", "", line)
#    outputLine <- tm.plugin.webmining::extractHTMLStrip(line)
#    #writeline
#    write(outputLine,outputFile)
#  }
  outputLine <- tm.plugin.webmining::extractHTMLStrip(eachFilePath,asText = FALSE)
  write(outputLine,outputFile)
        
  close(fileConnection)
  closeAllConnections() #ended up with files not closed, so lets really close them
  
  if(moveOldFiles){
  #and move the source file out of the main directory to avoid showing two similar files to the next bit of code
      file.copy(eachFilePath, paste(storeOriginalsPath, fileName, sep = ""), copy.date = TRUE, copy.mode = TRUE)      
  }
  file.remove(eachFilePath)
}

#return the working directory
setwd(oldWD)
closeAllConnections()
return(outputPath)

}

#===================================
#Function getCitations - First version 6/8/2017 by D. Graham Andrews
#===================================
#'
#' @return a dataframe with details on the documents which cite all results of the search, col1 = PubMedID, col2=Title, col3=date published, col4=PubMedID of doc cited
#' @export
#' @param searchterms searchResultsDataFrame: a dataframe as output by searchWordsPubMed()
#' @retmax how many arguments to return
#' @description Takes output directly from epartexttools:searchWordsPubMed() and returns a matrix of the pubmed IDs titles, and dates of all citations
#' @examples
#'output_text <- getCitations(searchWordsPubMed("influenza"),retmax=10))
#
#Uses searches based on results[4] (year) and results[6] (PMID)
#
getCitations <- function(searchResultsDataFrame, retmax=25){
  pubmedSearch <- searchResultsDataFrame
  retmax <- retmax #the maximum number of results returned
  #  if(!inherits(class(data.frame()))){
  #    warning("input class to function getCitations is not a dataframe. It's supposed to be a data frame.")
  #  }
  x <-0
  publicationIds <- pubmedSearch[6]
  for(eachID in publicationIds){
    queryText <- paste(toString(eachID),"[uid]",sep="")
    queryTemp <- RISmed::EUtilsQuery(queryText, type ="esearch", db="pubmed", retmax=retmax)
    temp<-httr::GET(queryTemp) %>% httr::content(.,"text") %>% xml2::as_xml_document(.) %>% xml2::as_list()
    t2 <- unlist(temp$IdList)
    pm1<-RISmed::EUtilsGet(t2,type="efetch", db="pubmed")
    
    result <- data.frame("Citations"=RISmed::Cited(pm1),"Title"=RISmed::ArticleTitle(pm1),"Grant Number"=RISmed::GrantID(pm1),"Affiliation"=RISmed::Affiliation(pm1),"YearAccepted"=RISmed::YearAccepted(pm1),"Agency"=RISmed::Agency(pm1),"id"=RISmed::ArticleId(pm1),"Author"=sapply(RISmed::Author(pm1),function(X) paste(X$LastName[X$order],X$Initials[X$order],sep=" ",collapse=";")),"ISSN"=RISmed::ISSN(pm1),"Acronym"=RISmed::Acronym(pm1),stringsAsFactors=FALSE)
    
    #pull out the citations list
    
  }
  for(i in 1:nrow(result[2])){ #result[2] will list the titles.  For each title we have, list print the data.
    print(paste("This paper, titled :", result[i,2], " has ", result[i,1], " citations.",sep=""))
  }
  return(result)
}



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


