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

findHiddenHTML("C:/Users/Graham/Documents/EPAR/HTML_TXT_TestFiles", moveOldFiles = FALSE)