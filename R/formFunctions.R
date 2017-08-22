
#' This function allows building a form extactor using the simple approach of selecting the text string that occurs before and after a desired string.
#'
#' @param docxml The xml table for a word document.
#' @param string1 String to start capture on.
#' @param string2 String to end capture on.
#' @return a character vector of length one
#' @export
#' @description Read the primary table from a Gates proposal word document, keeping column names. 
#' @examples
#' gates3<-bind_rows(lapply(tfiles,outcome_extractor))
cell_extractor<-function(docxml,string1,string2){paste(docxml[c(which(stringr::str_detect(docxml, string1))+1):c(which(str_detect(docxml, string2))-1)],collapse="\n",sep="\n")}


#' Read in a Gates proposal template.
#'
#' @param tfile File to read in for scraping table
#' @return a data frame, if only file exists, there was an error.
#' @export
#' @description Read the primary table from a Gates proposal word document, keeping column names. 
#' @examples
#' gates<-bind_rows(lapply(tfiles,readfiletab_gates))
readfiletab_gates<-function(tfile){
print(tfile)
propnar1<-read_docx(tfile)
docx_extract_tbl(propnar1,1)
propnartab<-docx_extract_tbl(propnar1,2)
if(length(propnartab)>0){  
  keepr<-which(nchar(propnartab[,1])>0)
  outro<-t(propnartab[keepr,2])
  colnames(outro)<-propnartab[keepr,1]
  keepr<-which(nchar(propnartab[,3])>0)
  outro2<-t(propnartab[keepr,4])
  colnames(outro2)<-propnartab[keepr,3]
  outro<-data.frame(cbind(outro,outro2))
  outro$file<-tfile
  outro} else {data.frame("file"=tfile)}}



#' Extract Outcomes from Gates proposal documents
#'
#' @param file File to read in for scraping.
#' @return a data frame, if only file exists, there was an error.
#' @export
#' @description Extracts the outcome table from a Gates proposal document.
#' @examples
#' outcome_extractor()
outcome_extractor<-function(file){
  #file<-"../epar/PracticeDocs2/files/OPP1140962_2015_Ombech_Proposal_Narrative.docx"
  # file<-tfiles[14]
  print(file)
  temptab<-docx_extract_all_tbls(read_docx(file)) %>% .[str_detect(.,"Primary Outcomes")]
  if(length(temptab)>=1){
    temptab<-temptab[[1]]
    temptab$Unstandard<-NA
  }
  if(length(temptab)<1){temptab=data.frame("Primary Outcomes"=NA,"Intermediate Outcomes"=NA,"Outputs"=NA,"Unstandard"=NA)}
  if(ncol(temptab)<4){temptab=data.frame("Primary Outcomes"=NA,"Intermediate Outcomes"=NA,"Outputs"=NA,"Unstandard"=paste(temptab[,1:ncol(temptab)],collapse="\n",sep="\n"))}
  colnames(temptab)<-c("Primary Outcomes","Intermediate Outcomes","Outputs","Unstandard")
  for(i in 1:ncol(temptab)){temptab[,i]<-as.character(temptab[,i])}
  temptab$file<-file
  temptab
}



#' See xml table contents of a docx file in a viewer window.
#'
#' @param file file to view xml table contents of
#' @param export_frame should the data frame be saved to output?
#' @return view of file
#' @export
#' @description a method for viewing the content of a word document via parsing the xml file to a data frame.
docx_table_view<-function(file,export_frame=FALSE,showView=T){
  #file<-tfiles[1]
  tmp <- tempfile()
  if (!dir.create(tmp)) 
    stop("Temporary directory could not be established.")
  utils::unzip(file, exdir = tmp)
  xmlfile <- file.path(tmp, "word", "document.xml")
  doc <- XML::xmlTreeParse(xmlfile, useInternalNodes = TRUE)
  nodeSet <- XML::getNodeSet(doc, "//w:p")
  unlink(tmp, recursive = TRUE)
  docxml <- sapply(nodeSet, XML::xmlValue)
  wordcontent<-dplyr::filter(data.frame("content"=na.omit(stringr::str_trim(docxml))),content!="")
  if(showView==T){View(wordcontent)}
  if(export_frame==TRUE){wordcontent}
}


#' Try to auto-extract form information.
#'
#' @param file_list list of word document files
#' @param multiple allow selecting multiple words
#' @param graphics graphical or boring interface.
#' @return 
#' @export
#' @description a method for viewing the content of a word
datapicker<-function(file_list,multiple=TRUE,graphics=FALSE){
temp<-lapply(file_list,function(X) docx_table_view(X,export_frame=T,showView=F))
#temp[[1]]$content %in% temp[[2]]$content %in% temp[[3]]$content
commonlist1<-lapply(temp,function(X) as.character(temp[[1]]$content[which(temp[[1]]$content %in%X$content)]))
temp<-matrix(commonlist1[[which.min(sapply(commonlist1,length))]])
whichitem<-select.list(temp,graphics=graphics,multiple=TRUE)
do.call(rbind,lapply(file_list,function(file1){ 
  temp1<-data.frame(matrix(ncol=length(whichitem)))
  colnames(temp1)<-whichitem
  temp1[,1:ncol(temp1)]<-as.character(temp1[,1:ncol(temp1)])
  temp1[1,]<-sapply(which(temp%in%whichitem),function(i){
    #cell_extractor_fixed_raw(file1,temp[i],temp[i+1])
    tmp <- tempfile()
    if (!dir.create(tmp)) 
      stop("Temporary directory could not be established.")
    utils::unzip(file1, exdir = tmp)
    xmlfile <- file.path(tmp, "word", "document.xml")
    doc <- XML::xmlTreeParse(xmlfile, useInternalNodes = TRUE)
    nodeSet <- XML::getNodeSet(doc, "//w:p")
    unlink(tmp, recursive = TRUE)
    docxml <- sapply(nodeSet, XML::xmlValue)
    startint<-which(str_detect(docxml, fixed(temp[i])))+1
    stopint<-which(str_detect(docxml,fixed(temp[i+1])))-1
    stopint<-if(length(stopint>1)){subset(stopint, stopint-startint>=0) %>% .[which.min(.-startint)]}
    if(length(stopint)<1){NA} else 
    {paste(docxml[c(startint:stopint)],collapse="\n",sep="\n")}
    })
  temp1$doc<-file1
  temp1}))}

cell_extractor<-function(docxml,string1,string2){
startint<-which(str_detect(docxml, string1))+1
stopint<-which(str_detect(docxml,string2))-1
stopint<-if(length(stopint>1)){subset(stopint, stopint-startint>=0) %>% .[which.min(.-startint)]}
paste(docxml[c(startint:stopint)],collapse="\n",sep="\n")}

docx_man_ex<-function(file,string1,string2,fixedString=TRUE){
  tmp <- tempfile()
  if (!dir.create(tmp)) 
    stop("Temporary directory could not be established.")
  utils::unzip(file, exdir = tmp)
  xmlfile <- file.path(tmp, "word", "document.xml")
  doc <- XML::xmlTreeParse(xmlfile, useInternalNodes = TRUE)
  nodeSet <- XML::getNodeSet(doc, "//w:p")
  unlink(tmp, recursive = TRUE)
  docxml <- sapply(nodeSet, XML::xmlValue)
  if(fixedString==TRUE){
  cell_extractor_fixed(docxml,string1,string2)} else {cell_extractor(docxml,string1,string2)}}

cell_extractor_fixed<-function(docxml,string1,string2){
  startint<-which(str_detect(docxml, fixed(string1)))+1
  stopint<-which(str_detect(docxml,fixed(string2)))-1
  stopint<-if(length(stopint>1)){subset(stopint, stopint-startint>=0) %>% .[which.min(.-startint)]}
  if(length(stopint)<1){NA} else 
  {paste(docxml[c(startint:stopint)],collapse="\n",sep="\n")}
  }

cell_extractor_fixed_raw<-function(file,string1,string2){
  tmp <- tempfile()
  if (!dir.create(tmp)) 
    stop("Temporary directory could not be established.")
  utils::unzip(file, exdir = tmp)
  xmlfile <- file.path(tmp, "word", "document.xml")
  doc <- XML::xmlTreeParse(xmlfile, useInternalNodes = TRUE)
  nodeSet <- XML::getNodeSet(doc, "//w:p")
  unlink(tmp, recursive = TRUE)
  docxml <- sapply(nodeSet, XML::xmlValue)
  startint<-which(str_detect(docxml, fixed(string1)))+1
  stopint<-which(str_detect(docxml,fixed(string2)))-1
  stopint<-if(length(stopint>1)){subset(stopint, stopint-startint>=0) %>% .[which.min(.-startint)]}
  if(length(stopint)<1){NA} else 
  {paste(docxml[c(startint:stopint)],collapse="\n",sep="\n")}
}

#' Cluster forms based on xml content.
#'
#' @param file_list list of word document files
#' @export
#' @description a method for viewing the content of a word
formcluster<-function(file_list){
alltables<-lapply(file_list, docx_table_view, export_frame=T,showView=F)
commatch<-sort(table(do.call(rbind,alltables)$content),decreasing=T) %>% .[nchar(names(.))>5] %>% names(.)
docmat<-matrix(do.call(c,lapply(alltables,function(AT) sapply(commatch,function(X) X%in%AT$content,USE.NAMES=F)/length(AT))),nrow=length(alltables))
hc<-hclust(dist(docmat))
plot(hc)}


