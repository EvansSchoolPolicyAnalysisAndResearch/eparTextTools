
expandJava<-function(){options(java.parameters = "-Xmx6000m")}

#list.of.packages <- c("ggplot2", "Rcpp","tm","ggthemes","SnowballC","rvest","downloader","DT","wordcloud","d3heatmap","plyr","reshape2","dplyr","qdapTools","stringr","openNLP","NLP","stm","LDAvis","servr","Rtsne","geometry","downloader","corrplot","pryr","openNLPmodels.en","lubridate","pbapply","devtools","tm.plugin.mail","plotly","data.table","igraph","networkD3")

#if("StanfordCoreNLP"%in%c(installed.packages()[,"Package"])==FALSE){install.packages('StanfordCoreNLP',repos="http://datacube.wu.ac.at/",type="source")}
#if("openNLPmodels.en"%in%c(installed.packages()[,"Package"])==FALSE){install.packages('openNLPmodels.en',repos="http://datacube.wu.ac.at/",type="source")}

#packages.Req <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

#if(length(packages.Req)) install.packages(packages.Req)

#lapply(list.of.packages, function(X) library(X,character=TRUE))

#package depends on xpdf... on mac, install via brew install xpdf and follow instructions.

demandwordvec<-function(){
  try({if (!require(wordVectors)) {
  if (!(require(devtools))) {
    install.packages("devtools")
  }
  devtools::install_github("bmschmidt/wordVectors")
}})
}


#' Build Example Data.
#'
#' @return a folder of documents
#' @export
#' @description  This function creates a folder if documents in the working directory.
#' @examples
#' example_documents()
example_documents<-function(){
 data(demo.docs)
if(file.exists("demo.docs.folder")==FALSE){dir.create("demo.docs.folder")}
demo.docs2<-demo.docs[which(file.exists(file.path("demo.docs.folder",sapply(demo.docs$links, function(X){tail(unlist(strsplit(X,split="/",fixed=TRUE)),n=1)})))==FALSE),]
lapply(demo.docs2$links, function(X){try(download.file(X, destfile=file.path("demo.docs.folder",tail(unlist(strsplit(X,split="/",fixed=TRUE)),n=1)),mode="wb"))})
}

#' Check to see if Tika is available.
#'
#' @param directory directory to look for tika in. 
#' @return If false will download tika to working directory.
#' @seealso \code{\link{http://apache.claz.org/tika}} 
#' @export
#' @description The check for tika command will install tika to the current working directory.
#' @examples
#' checkForTika()
checkForTika<-function(directory=getwd()){if("tika-app-1.13.jar"%in%list.files(path=directory)) {cat("success")} else {download.file("http://apache.claz.org/tika/tika-app-1.13.jar",file.path(directory,"tika-app-1.13.jar"))}}
#' Loads and processes doc,docx,pdf,and txt files into tm corpus.
#'
#'This command loads files into R from a directory into a corpus. Currently it reads doc, docx, pdf, and txt files.
#' @param fname Directory name
#' @param tika Should tika be used? 
#' @param gen_pdf_tools if using a pc, set to true. if on linux, setting to false and installing the pdf2text library will be more useful. See readPDF in the tm package for more info.
#' @param tikapath Path to the Tika application
#' @return TM Text Corpus.
#' @seealso \code{\link{tm::corpus}} 
#' @export
#' @examples
#' corpus1<-getTextR(file.path())
getTextR<-function(fname,tika=FALSE,gen_pdf_tools=T,tikapath="tika-app-1.13.jar"){
  fname <- sapply(fname,tolower)
  if(tika==TRUE){
    pdoc<-system(command=paste("java -jar",tikapath,"-t",gsub(" ","\\ ",fname,fixed=TRUE)),intern=TRUE,wait=TRUE)
  } else {
    pdoc<-if(stringr::str_detect(fname,".docx+$")==TRUE){read_docxtm(fname)} #if it's a docx, use read_docx
    else {
      if(stringr::str_detect(fname,".doc+$")==TRUE){tm::readDOC()(language="en",elem=list(uri=fname))} #if its a .doc, use readDOC
      else {pdoc<-if(stringr::str_detect(fname, stringr::fixed(".pdf"))==TRUE){ #if it's a PDF, we choose from two
        if(gen_pdf_tools==F){readPDF2(engine="xpdf")(elem=list(uri=fname), language="en")} #read pdf with the xpdf engine
        else {readPDF2(engine="pdftools")(elem=list(uri=fname), language="en")}} #or read pdf with the pdftools engine
      else {if(stringr::str_detect(fname,stringr::fixed(".txt"))==TRUE){ #if it's a raw text file
        tm::readPlain(elem=list(uri=fname,content=iconv(enc2utf8(readLines(fname)), sub = "byte")),language="en")}
        else {"getTextR: FILETYPE UNKNOWN"}}}}}
  pdoc
}

#' Read Emails into Corpus.
#'
#' @param folder_in The path to the folder where the emails are located
#' @param newmailsdirectory The path to the folder where the processed emails will be stored.
#' @return A tm Corpus of processed emails.
#' @seealso \code{\link{http://www.matijs.net/software/msgconv/}} 
#' @export
#' @description  A tm Corpus of processed emails which can be combined with a corpus created by tm using the c() command. The command calls to msgconvert, and will only work if you have msgconvert installed.
#' @examples
#' readMails("../doc","../procdocs")
readMails<-function(folder_in,newmailsdirectory){
  dir.create(newmailsdirectory)
  messages<-list.files(folder_in,full.names=TRUE)[str_detect(tools::file_ext(list.files(folder_in)),"msg")]
  plyr::l_ply(messages,function(X) system(paste(paste("msgconvert --mbox",file.path(folder_in,"msgs.mbox",sep=""),X))))
  file.remove(messages)
  mbases<-list.files(folder_in)[tools::file_ext(list.files(folder_in))%in%"mbox"==FALSE]
  file.copy(file.path(folder_in,"msgs.mbox"),file.path(newmailsdirectory,"msgs.mbox"))
  mails<-MBoxSource(file.path(newmailsdirectory,"msgs.mbox"))
  file.remove(file.path(folder_in,"msgs.mbox"))
  mailsc<-Corpus(mails)
  mailsc}

#' Tests if save directory exists and creates new one if not
#'
#' @param dir directory
#' @export
#' @description  This function does the quite useful task of making a workingfolder object in the environment.
makeworking<-function(dir){if(dir.exists(dir)==FALSE) {dir.create(dir)} else {cat("Good News, Directory Already Exists Containing:",paste(list.files(dir),collapse="\n"))}
  workingfolder<<-dir}


#' Creates directory and copies from old to new.
#'
#' @param from Old directory
#' @param to New directory
#' @return NA, but he file will be relocated
#' @seealso \code{\link{file.copy}} 
#' @export
#' @description  This function is just a basic wrapper to copy and move folders recursively in linux. Should probably just use file.copy instead.
#' @examples
#' copyDir("../doc","../newdocs")
copyDir<-function(from,to){
  dir.create(to)
  system(paste("cp -r", from, to))}


#' Reads docx files into tm corpus.
#'
#' @param file File to read. Should end in docx
#' @param skip Value is settable for debugging. Skips lines. Set at 0 for most applications.
#' @return a plain text document corpus
#' @seealso \code{\link{PlainTextDocument}} 
#' @export
#' @description  This function reads a docx file into a tm corpus.
#' @examples
#' read_docxtm("demodoc.docx",skip=0)
read_docxtm<-function (file, skip = 0) {
  tmp <- tempfile()
  if (!dir.create(tmp)) 
    stop("Temporary directory could not be established.")
  utils::unzip(file, exdir = tmp)
  xmlfile <- file.path(tmp, "word", "document.xml")
  doc <- XML::xmlTreeParse(xmlfile, useInternalNodes = TRUE)
  xml_metadata<-read_xml(file.path(tmp, "docProps", "core.xml")) 
  unlink(tmp, recursive = TRUE)
  nodeSet <- XML::getNodeSet(doc, "//w:p")
  pvalues <- sapply(nodeSet, XML::xmlValue)
  pvalues <- pvalues[pvalues != ""]
  if (skip > 0) 
    pvalues <- pvalues[-seq(skip)]
  tempdes<-xml_text(xml_metadata)
  xml_metadata<-xml_metadata %>% as_list() 
  PlainTextDocument(x=pvalues,author=xml_metadata$creator[[1]],description=tempdes,datetimestamp=xml_metadata$modified[[1]], id=file)}



#' Worker function for processing URI info. From tm package.
#'
#' @param uri the path to the uri to be processed.
#' @return uri
#' @seealso \code{\link{tm}} 
#' @export
#' @description  Worker function for processing URI info. From tm package.
processURI2 <-function(uri) {
  uri <- as.character(uri)
  if (identical(substr(uri, 1, 7), "file://"))
    uri <- substr(uri, 8, nchar(uri))
  uri
}


#' Function call to read pdf files into tm corpus.
#'
#' @param engine which pdf library to use. see tm readPDF for more details
#' @return engine for reading pdf
#' @seealso \code{\link{readPDF}} 
#' @export
#' @description  This function is a port of readpdf.
#' @examples
#' readPDF2(engine="xpdf")(elem=list(uri=fname), language="en")
readPDF2<-function (engine = c("xpdf", "Rpoppler", "ghostscript", "Rcampdf","custom","pdftools"), control = list(info = NULL, text = NULL)) {
  stopifnot(is.character(engine), is.list(control))
  engine <- match.arg(engine)
  pdf_info <- switch(engine, xpdf = function(x) tm:::pdf_info_via_xpdf(x,control$info), Rpoppler = Rpoppler::PDF_info, pdftools=pdftools::pdf_info, ghostscript = tm:::pdf_info_via_gs,Rcampdf = Rcampdf::pdf_info, custom = control$info)
  pdf_text <- switch(engine, xpdf = function(x) system2("pdftotext", c(control$text, shQuote(x), "-"), stdout = TRUE), Rpoppler = Rpoppler::PDF_text, pdftools=pdftools::pdf_text, ghostscript = pdf_text_via_gs, Rcampdf = Rcampdf::pdf_text, custom = control$text)
  if (!is.function(pdf_info) || !is.function(pdf_text)) 
    stop("invalid function for PDF extraction")
  function(elem, language, id) {
    uri <- processURI2(elem$uri)
    meta <- pdf_info(uri)
    content <- pdf_text(uri)
    content<-iconv(enc2utf8(content), sub = "byte")
    tm::PlainTextDocument(content, meta$Author, meta$CreationDate, 
                          meta$Subject, meta$Title, basename(elem$uri), language, 
                          meta$Creator)
  }
}


#' Look up how to open a file. 
#'
#' @param ending the ending of the file in question.
#' @export
how_do_i_open<-function(ending){
  browseURL(paste("http://lmgtfy.com/?q=how+do+i+open+a+",ending,"+file+in+R",sep=""))
}

#' Look up how to do something in R
#'
#' @param theproblem the problem.
#' @description  A joke command to look up a method in R.
#' @export
how_do_i<-function(theproblem){
  browseURL(paste("http://lmgtfy.com/?q=how+do+i+",theproblem,"+in+R",sep=""))
}


#' Calls getTextR on all files in a directory joining into corpus. 
#'
#' @param directory Folder to read files from'
#' @param onError if skip, skip documents that read wrong.
#' @param gen_pdf_tools if using a pc set to true. if on a unix server, setting to false will preserve pdf metadata.
#' @return Corpus of text documents.
#' @seealso \code{\link{corpus}} 
#' @export
#' @description  This function will read from a folder documents of the class pdf, docx, doc or txt.
#' @examples
#' allDocs("folder")
allDocs<-function (directory, SkiponError = FALSE,gen_pdf_tools=TRUE) {
  if(SkiponError==TRUE){
    temp <- lapply(file.path(directory, list.files(directory, recursive = TRUE, include.dirs = FALSE)), 
                   function(FILENAME) {
                     try(getTextR(FILENAME, gen_pdf_tools=gen_pdf_tools))
                   })
    temp <- temp[lapply(temp, class) != "character"]
    temp <- temp[lapply(temp, class) != "error"]
    temp <- temp[lapply(temp, is.null) != TRUE]
    temp <- temp[lapply(temp,class)!="try-error"]
    do.call(c, temp)
  } else {
    do.call(c, lapply(file.path(directory, list.files(directory, recursive = TRUE, include.dirs = FALSE)),getTextR))
  }
}

#' Cleans documents performing many common tasks .
#'
#' @param corpus corpus to clean
#' @return a cleaned tm corpus
#' @seealso \code{\link{tm_map}} 
#' @export
#' @description  This function strips whitespace, removes stop words,removes punctuation,and stems documents.
#' @examples
#' doc_clean_process(corpus1)
doc_clean_process<-function(corpusname){
  stopWords <- function(x) tm::removeWords(x, tm::stopwords("en"))
  funs <- list(tm::stripWhitespace,
               stopWords,
               tm::removeNumbers,
               tm::removePunctuation,
               tm::stemDocument,
               tm::content_transformer(tolower))
  corpus2<-tm::tm_map(corpusname, FUN = tm::tm_reduce, tmFuns = funs, mc.cores=1)
  corpus2}

#' Makes a pretty word association table .
#'
#' @param assoctable rough word association table
#' @param corpus word corpus for table
#' @param ngram
#' @return a cleaned tm corpus
#' @seealso \code{\link{tm_map}} 
#' @export
#' @description  This function takes a findAssoc object and returns a datatable. assocPrettyOneStep is wraps this to allow everything in one step.
#' @examples
#' assocPTable(findAssoc(corpus1,wordlist,corrVal),corpus)
assocPTable<-function(assoctable,corpus,ngram=FALSE){
  #assoctable<-assoctable[sapply(assoctable,length)>0]
  dft<-do.call(rbind,lapply(1:length(assoctable),function(i){tryCatch({data.frame("Word"=names(assoctable)[i],"Match"=names(assoctable[[i]]),"Association"=c(assoctable[[i]]))},error=function(e){data.frame("Word"=names(assoctable)[i],"Match"="too few words","Association"=c(0))})}))
  dft$Word<-as.character(dft$Word)
  dft$Match<-as.character(dft$Match)
  #dft$Word[dft$Word%in%c(names(assoctable)[sapply(assoctable,length)>0])]<-stemCompletion(dft$Word[dft$Word%in%c(names(assoctable)[sapply(assoctable,length)>0])],dictionary=corpus,type="prevalent")
  dft$Match[dft$Word%in%c(names(assoctable)[sapply(assoctable,length)>0])]<-stemCompletion(dft$Match[dft$Word%in%c(names(assoctable)[sapply(assoctable,length)>0])],dictionary=corpus,type="prevalent")
  DT::datatable(data=dft,rownames=FALSE,filter="top")}

#' Makes a pretty word association table all in one step.
#'
#' @param wordlist list of words for association
#' @param termDocumentMatrix tdm of corpus
#' @param corpus the corpus the words are from
#' @param corrVal the correlation value to report at or above
#' @return a word association table
#' @seealso \code{\link{findAssoc}} 
#' @export
#' @description  This function takes a wordlist, tdm, coprus, and settable correlation value and returns a datatable. 
#' @examples
#' assocPrettyOneStep(c("research","find"),termDocumentMatrix(corpus),corpus1,.8)
assocPrettyOneStep<-function(wordlist,termDocumentMatrix,corpus,corrVal=.8){
  assocPTable(findAssocs(termDocumentMatrix,wordlist,corrVal),corpus)
}


#' Creates plots of word frequencies
#'
#' @param termDocumentMatrix tm TermDocumentMatrix object
#' @param typeplot 1=word on x, 2=variable on x
#' @param wordcount How many words to graph 
#' @param minfreq minimum number of times a word should occur to be reported
#' @param shortendoc Should document names be shortened?
#' @param allCorpus For plot type 1, would you like the plot for all documents combined?
#' @return a ggplot2 object.
#' @seealso \code{\link{ggplot}} 
#' @export
#' @description  This function creates bar graphs in ggplot. It can be customized by adding items to the returned object.
#' @examples
#' wfplots(TermDocumentMatrix(corpus1),typePlot=1,5,minfreq=5,shortendoc=T)
wfplots<-function(termDocumentMatrix,typePlot=1,wordcount,minfreq=5,shortendoc=FALSE, allCorpus=FALSE){
  mfcomframe<-data.frame(as.matrix(termDocumentMatrix[findFreqTerms(termDocumentMatrix, lowfreq=minfreq),]))
  mfcomframe<-mfcomframe[sort(rowSums(mfcomframe),index.return=TRUE,decreasing=TRUE)$ix[1:wordcount],]
  mfcomframe$word<-row.names(mfcomframe)
  mfcomframe$word<-factor(mfcomframe$word, levels = mfcomframe$word)
  mfcomframe<-reshape2::melt(mfcomframe,id=c("word"))
  if(shortendoc==TRUE){
    mfcomframe$variable<-as.character(mfcomframe$variable)
    mfcomframe$variable<-stringr::str_trunc(mfcomframe$variable,side="left",width=20,ellipsis="...")
  }
  if(typePlot==1){
    if(allCorpus==TRUE){
      plotout<-ggplot(mfcomframe)+geom_bar(aes(x=word,y=value),position="stack",stat="identity")+coord_flip()+ggthemes::theme_pander()+ylab("Frequency")+xlab("Word")
      } else {
    plotout<-ggplot(mfcomframe)+geom_bar(aes(x=word,y=value,fill=variable),position="stack",stat="identity")+coord_flip()+ggthemes::theme_pander()+ggthemes::scale_fill_tableau(name="Document")+ylab("Frequency")+xlab("Word")
    }} else {
    plotout<-ggplot(mfcomframe)+geom_bar(aes(x=variable,y=value,fill=word),position="stack",stat="identity")+coord_flip()+ggthemes::theme_pander()+ggthemes::scale_fill_tableau(name="Word",palette="tableau10")+ylab("Frequency")+xlab("Document")
    }
  plotout
}



#' TermDocumentMatrix to interactive heatmap.
#'
#' @param termDocumentMatrix tm TermDocumentMatrix object
#' @param wordcount How many words to graph 
#' @param minfreq minimum number of times a word should occur to be reported
#' @param pickwords list of words to be included in heatmap. Default is all words sorted by frequency.
#' @param col_labels document labels in order of tdm columns
#' @param dendrows draw a dendogram for rows 
#' @param dendcolumns draw a dendogram fo columns
#' @return a d3heatmap object.
#' @seealso \code{\link{d3heatmap}} 
#' @export
#' @description  This function creates a d3 interactive heatmap
#' @examples
#' word_heatmap(TermDocumentMatrix(corpus1),10)
word_heatmap<-function(termDocumentMatrix,wordcount,minfreq=2,col_labels=paste(1:ncol(termDocumentMatrix),stringr::str_sub(colnames(termDocumentMatrix),-10,-1)),dendrows=TRUE,dendcolumns=TRUE,pickwords=c()){
  if(length(pickwords)<=0){
  mfcomframe<-data.frame(as.matrix(termDocumentMatrix[findFreqTerms(termDocumentMatrix, lowfreq=minfreq),]))
  mfcomframe<-mfcomframe[sort(rowSums(mfcomframe),index.return=TRUE,decreasing=TRUE)$ix[1:wordcount],]
  } else { 
    mfcomframe<-data.frame(as.matrix(termDocumentMatrix[pickwords,]))
  }                           
  mfcomframe$word<-row.names(mfcomframe)
  mtcells<-as.matrix(termDocumentMatrix[mfcomframe$word,])
  d3heatmap::d3heatmap(termDocumentMatrix[mfcomframe$word,], scale="column",colors="Purples",Rowv=dendrows,Colv=dendcolumns, labCol=col_labels)}

#' Creates ggplot of wordcounts!
#'
#' @param termDocumentMatrix tm TermDocumentMatrix object
#' @param wordlist list of words
#' @param by.var can specify variable to group documents by
#' @param name plottable name of by.var 
#' @return a ggplot2 object
#' @seealso \code{\link{ggplot2}} 
#' @export
#' @description  This function ggplot of wordcounts
#' @examples
#' interest_plot("research",TermDocumentMatrix(corpus1),by.var=variable,"DocID")
interest_plot<-function(wordlist,termDocumentMatrix,by.var=NULL,byvarname=""){
  tempframe<-data.frame(as.matrix(termDocumentMatrix[wordlist,]))
  tempframe<-data.frame("Count"=rowSums(tempframe),"word"=row.names(tempframe))
  if(is.null(by.var)){
    ggplot(tempframe, aes(word, Count)) + geom_bar(fill="#8ebfad", position = "dodge", stat="identity") + theme(axis.text.x=element_text(color="#000000",angle=50, hjust=1, size=12),panel.background=element_blank())+xlab("")+ylab("Frequency")}
  else {
    ggplot(tempframe, aes(word, Count, fill=by.var)) + geom_bar(position = "dodge", stat="identity") + theme(axis.text.x=element_text(color="#000000",angle=50, hjust=1, size=12),panel.background=element_blank())+xlab("")+ylab("Frequency")+ggthemes::scale_fill_pander(name=as.character(byvarname))}
}


#' Creates ggplot of wordcounts by document
#'
#' @param termDocumentMatrix tm TermDocumentMatrix object
#' @param wordlist list of words
#' @return a ggplot2 object
#' @seealso \code{\link{ggplot2}} 
#' @export
#' @description  This function ggplot of wordcounts by document id.
#' @examples
#' interest_plot_bydoc("research",TermDocumentMatrix(corpus1))
interest_plot_bydoc<-function(wordlist,termDocumentMatrix){
  tempframe<-data.frame(as.matrix(termDocumentMatrix[wordlist,]))
  tempframe$word<-row.names(tempframe)
  tempframe<-reshape2::melt(tempframe,id=c("word"))
  tempframe$variable<-as.character(tempframe$variable) %>% stringr::str_trunc(.,20,side="left")
  if(length(wordlist)>1){
    ggplot(tempframe, aes(variable, value, fill=word)) + geom_bar(position = "stack", stat="identity") + theme(axis.text.x=element_text(color="#000000",angle=50, hjust=1, size=12),panel.background=element_blank())+xlab("")+ylab("Frequency")+ggthemes::scale_fill_pander()
  } else {
    ggplot(tempframe, aes(variable, value)) + geom_bar(fill="#8ebfad", position = "stack", stat="identity") + theme(axis.text.x=element_text(color="#000000",angle=50, hjust=1, size=12), panel.background=element_blank())+xlab("")+ylab("Frequency")+ggthemes::scale_fill_pander()
  }}

#' Generates tornado plot
#'
#' @param termmatrix tm TermDocumentMatrix object
#' @param pickword a vector of words to select based on
#' @param frequency how often should the words occur in the document
#' @param nwords how many frequent words should be graphed?
#' @return a ggplot2 object
#' @seealso \code{\link{ggplot2}} 
#' @export
#' @description  This function makes a tornado plot comparing documents based on selected word frequencies. For example, can compare documents containing the word "gender" to all other documents
#' @examples
#' interest_plot_bydoc("research",TermDocumentMatrix(corpus1))
tornadoCompare<-function(termmatrix,pickword,frequency,nwords) {
  words_gender<-as.matrix(termmatrix[,as.vector(termmatrix[pickword,])>frequency]) %>% rowSums()
  words_gender<-words_gender/sum(words_gender)
  words_notgender<-as.matrix(tdm[,as.vector(termmatrix[pickword,])<=frequency]) %>% rowSums()
  words_notgender<-words_notgender/sum(words_notgender)
  bothgen<-cbind2(words_gender,words_notgender)
  bothgen<-as.TermDocumentMatrix(bothgen,weighting=weightTf) 
  bothgen$dimnames$Docs<-c("Chosen","Inverse")
  top10<-unique(c(sapply(1:2,function(X) sort(bothgen[,X]$v,decreasing=TRUE,index.return=T)$ix[1:nwords])))
  bothgen<-reshape2::melt(as.matrix(bothgen[top10,]))
  bothgen$value[bothgen$Docs=="Inverse"]<-bothgen$value[bothgen$Docs=="Inverse"]*-1
  ggplot(bothgen)+geom_bar(aes(x=Terms,y=value,fill=as.factor(Docs)),stat="identity",position=position_dodge(width=0.0))+coord_flip()+theme_minimal()+scale_fill_brewer("Document",palette="Dark2")+ylab("relative frequency")
}


#' Creates ggplot of wordcounts by document, allowing external characteristics.
#'
#' @param termDocumentMatrix tm TermDocumentMatrix object
#' @param wordlist list of words
#' @param doccharacteristic a vector the same length as tm specifying the trait you'd like to plot by
#' @return a ggplot2 object
#' @seealso \code{\link{ggplot2}} 
#' @export
#' @description  This function ggplot of wordcounts by document characteristic from an external list.
#' @examples
#' interest_plot_bydoc("research",TermDocumentMatrix(corpus1),rep("bacon",length(corpus1)))
interest_plot_bydoc_char<-function(wordlist,termDocumentMatrix,doccharacteristic){
  termDocumentMatrix$dimnames$Docs<-doccharacteristic
  tempframe<-data.frame(as.matrix(termDocumentMatrix[wordlist,]))
  tempframe$word<-row.names(tempframe)
  tempframe<-reshape2::melt(tempframe,id=c("word"))
  if(length(wordlist)>1){
    ggplot(tempframe, aes(variable, value, fill=word)) + geom_bar(position = "stack", stat="identity") + theme(axis.text.x=element_text(color="#000000",angle=50, hjust=1, size=12),panel.background=element_blank())+xlab("")+ylab("Frequency")+ggthemes::scale_fill_pander()
  } else {
    ggplot(tempframe, aes(variable, value)) + geom_bar(fill="#8ebfad", position = "stack", stat="identity") + theme(axis.text.x=element_text(color="#000000",angle=50, hjust=1, size=12), panel.background=element_blank())+xlab("")+ylab("Frequency")+ggthemes::scale_fill_pander()
  }}


#' Suggests keywords based on the conceptnet API
#'
#' @param keyword word to search conceptnet for
#' @param path_root if true, can put in word path returned from previous calls otherwise defaults to keyword.
#' @export
#' @example suggestkeywords("polio") %>% suggestkeywords(.$word[3],path_root=T)
suggestkeywords<-function(keyword,path_root=F){
  cnetpath<-if(path_root==T){"http://api.conceptnet.io"} else {"http://api.conceptnet.io/c/en"}
  temp<-jsonlite::fromJSON(file.path(cnetpath,keyword))
  temp<-lapply(temp$edges$start$term, function(keyword1) {
    temp2<-jsonlite::fromJSON(file.path("http://api.conceptnet.io",keyword1))
    data.frame("label"=temp2$edges$start$label,"term"=temp2$edges$start$term)})
  list("word"=c(c(sapply(temp, function(X) as.character(X$term)) %>% unlist(.) %>% unique()) %>% unique),"term"=c(c(sapply(temp, function(X) as.character(X$label)) %>% unlist(.) %>% unique())%>% unique))
}




#' Creates a table of word counts within a set of documents.
#'
#' @param termDocumentMatrix tm TermDocumentMatrix object
#' @param wordlist list of words
#' @param doccorpus original document corpus
#' @param trunc should filenames be truncated in output table
#' @param weighting passed to tm::TermDocumentMatrix (by default is weightTf)
#' @param raw if true then outputs a data.table if false then produces an htmltable allowing csv/excel export. the resulting table can be saved using the saveWidgets command in htmlwidgets.
#' @param onlywordlist if true, bins words, otherwise returns all unique matches
#' @return a data.table or DT datatable depending on if raw is set to true or false.
#' @description  This function is useful for looking at occurences of words within documents
#' @export
#' @examples
#' wordcount_table(c("gender","nutrition"),TermDocumentMatrix(clcorp),corpus1)
wordcount_table<-function(wordlist,termDocumentMatrix,doccorpus,trunc=FALSE,raw=FALSE,onlywordlist=F){
  truelist<-unique(unlist(sapply(wordlist, function(X) which(str_detect(termDocumentMatrix$dimnames$Terms,X)),USE.NAMES = FALSE)))
  tempframe<-data.frame(as.matrix(termDocumentMatrix[truelist,]))
  tempframe$word<-row.names(tempframe)
  tempframe<-reshape2::melt(tempframe,id=c("word"))
  tempframe$variable<-as.character(tempframe$variable)
  tempframe<-dplyr::filter(tempframe, value>0)
  tempframe<-tempframe[order(nchar(tempframe$word),tempframe$value,decreasing=c(F,T)),]
  sents<-adply(tempframe,.margins=1,.fun=function(X){
    sframe<-unlist(tokenizers::tokenize_sentences(paste(textreg::convert.tm.to.character(doccorpus[which(gsub("%",".",names(doccorpus),fixed=TRUE)==X[,2])]),collapse=" ")))
    sframe<-sframe[unique(which(str_detect(sframe, X[,1])))] %>% paste(str_trim(.),collapse="...   ")
    sframe
  })
  if(trunc==TRUE) {sents$variable<-sents$variable %>% stringr::str_trunc(.,20,side="left")}
  colnames(sents)<-c("Word","Document","Count","String")
  sents<-filter(sents,nchar(String)>0)
  if(onlywordlist==TRUE){
    sents<-lapply(wordlist,function(W){
      temp1<-filter(sents, str_detect(Word,W))
      temp1<-ddply(temp1,.(Document), summarize,Count=sum(Count),String=paste(String,collapse="...   "))
    temp1
  })
    names(sents)<-wordlist
    sents<-bind_rows(sents,.id="Word")
}
if(raw==TRUE){data.table::data.table(sents)} else {
  DT::datatable(sents, options = list(columnDefs = list(list(
    targets = 4,
    render = DT::JS(
      "function(data, type, row, meta) {",
      "return type === 'display' && data.length > 50 ?",
      "'<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>' : data;",
      "}")
  ),list(
    targets = 2,
    render = DT::JS(
      "function(data, type, row, meta) {",
      "return type === 'display' && data.length > 20 ?",
      "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
      "}")
  )),
  dom = 'Bfrtip',
  buttons = c('csv', 'excel')), extensions = 'Buttons', callback = DT::JS('table.page(3).draw(false);'),filter="top")
  }}

#==============================
#The OCR_DOCS() function.  Original version By D. Graham Andrews and Ryan Scott, first finished 5/16/2017
#==============================
#' @Title \code{OCR_DOCS} will read a given folder for all PDFs in it
#' @description This function uses Google's "parsey mcparseface" neural network to read PDFs and output machine readable text files of the content.  It moves the source files up one directory into a director called OCR_Sources, and replaces them in the working directory with text versions.
#' @param path: A path to a directory containing mixed pdf/other documents.  Runs OCR on all the PDF documents. (Does not detect which PDFs have existing text.)
#' @returns output path: creates text files in the local working directory with the same names as the PDFs to be OCR'd.
#' @examples
#' text_files_location <- OCR_DOCS(path_to_pdf_files)
#' @export
OCR_DOCS<-function(path){
  #save the working directory
  oldWD <- getwd()
  #change working directory
  setwd(path)
  #build list of all PDFs at "path" and put that list in "listOfPDFs"
  listOfPDFs <- list.files(path,"\\.pdf")
  
  #create the output folder and the source folder
  #outputPath <- toString(paste(path,"OCR_Results/",sep="")) This worked, but instead of putting the output in a new folder we're moving the sources away and keeping the output where the sources were
  outputpath <- toString(path)
  setwd('..')
  sourcePath <- toString(paste(getwd(),"OCR_Sources/", sep=""))
  setwd(path) #put the working directory back after creating the sources path.
  
  dir.create(outputPath)
  
  for(eachPDFpath in listOfPDFs){
    pdfInfo<-pdftools::pdf_info(eachPDFpath)
    numPages<-pdfInfo$pages
    fileName<-eachPDFpath
    if(.Platform$OS.type == "windows"){
      fileName<-strsplit(eachPDFpath, "/") #if we have a full path, cut it into chunks based on the /
      fileName<-fileName[length(fileName)] #set the filename to the last chunk of the path
    }
    else{print("Warning: OS is not windows, so remember that OCR_DOCs() function must be called with filenames and not paths on other OS's.")}
    
    #prep the output filename (remove the ".pdf" so as not to confuse later parts of texttools
    fileName <- sapply(fileName,tolower)
    outputFileName <- gsub(".pdf", "", fileName)
    #prep a file connection to output the raw text
    outputFile <- file(description=paste("rawText_",toString(outputFileName),".txt",sep=""), open="a", blocking=FALSE, encoding="", raw=FALSE, method="internal")
    
    tout<-lapply(1:numPages,function(PageNumber){
      
      #make a new temporary png image file
      tmpf<-tempfile("page", fileext = c(".png"))
      
      #Render a page and save it to a new png file
      print(paste("Now rendering page ", PageNumber, " of ", numPages, " in  document ", fileName, ".", sep=""))
      pdftools::pdf_render_page(toString(fileName),page=PageNumber,dpi=300, numeric=TRUE) %>% png::writePNG(tmpf)
      
      #run the OCR and then discard the temp file
      raw_text<-tesseract::ocr(tmpf,engine=tesseract("eng"))
      unlink(tmpf)
      
      #save the raw text
      
      write(raw_text,outputFile)
    }) #this is the end of the lapply function definition, so a ) looks odd here but belongs
    close(outputFile)
    #and move the source file out of the main directory to avoid showing two similar files to the next bit of code
    file.copy(eachPDFpath, paste(sourcePath, fileName, sep = ""), copy.date = TRUE, copy.mode = TRUE)
    file.remove(eachPDFpath)
  }
  #Move all files except the PDFs to a new folder
  #listOfFiles <- setdiff(list.files(path),listOfPDFs)
  #file.copy(file.path(path, listOfFiles),outputPath)
  
  #return the working directory
  setwd(oldWD)
  
  return(outputPath)
}

#============================
#readDocInitWindows()
#============================
#by D. Graham Andrews.  Original ver 5.30.2017
#Inputs: None
#outputs: No return value. Writes file to disk at c:/antiword/8859-1.txt if it doesn't exist.
#The purpose of this function is to cover for incomplete fresh installations of antiword, which don't populate environment variables with suitable mapping files to read the .doc files and end up failing.
#When the environment variables for antiword are missing, it defaults to checking in c:/antiword/ so by populating our file there it will succeed.
#' @Title \code{readDocInitWindows} will add the required parsing files to call readDoc() successfully on a machine that has just installed antiword but not placed the needed parser files in the needed directories.
#' @description Only needed for windows installations: Places a copy of the file 8859-1-2.txt into the c:/antiword folder.  Antiword is hard coded to look for this file when parsing .doc files, but does not place the file there itself. Antiword fails silently when this file is not present. This function is for our package to detect this case and place the file before calling antiword so that it won't fail.
#' @param none
#' @returns none
#' @examples
#' if(os_is_windows & parser_file_missing){readDocInitWindows()}
#' @export
readDocInitWindows <- function(){
  if(.Platform$OS.type == "windows" & file.exists("c:/antiword")==FALSE){
    warning("Antiword tool for reading .doc files was not fully installed.  Creating default mapping file (ENGLISH UNICODE) for antiword at c:/antiword where antiword searches for them by default")
    dir.create("C:/antiword")
    outputFile <- file(description="c:/antiword/8859-1-2.txt", open="a", blocking=FALSE, encoding="", raw=FALSE, method="internal")
    fileData = "#
    #	Name:             ISO/IEC 8859-1:1998 to Unicode
    #	Unicode version:  3.0
    #	Table version:    1.0
    #	Table format:     Format A
    #	Date:             1999 July 27
    #	Authors:          Ken Whistler <kenw@sybase.com>
    #
    #	Copyright (c) 1991-1999 Unicode, Inc.  All Rights reserved.
    #
    #	This file is provided as-is by Unicode, Inc. (The Unicode Consortium).
    #	No claims are made as to fitness for any particular purpose.  No
    #	warranties of any kind are expressed or implied.  The recipient
    #	agrees to determine applicability of information provided.  If this
    #	file has been provided on optical media by Unicode, Inc., the sole
    #	remedy for any claim will be exchange of defective media within 90
    #	days of receipt.
    #
    #	Unicode, Inc. hereby grants the right to freely use the information
    #	supplied in this file in the creation of products supporting the
    #	Unicode Standard, and to make copies of this file in any form for
    #	internal or external distribution as long as this notice remains
    #	attached.
    #
    #	General notes:
    #
    #	This table contains the data the Unicode Consortium has on how
    #       ISO/IEC 8859-1:1998 characters map into Unicode.
    #
    #	Format:  Three tab-separated columns
    #		 Column #1 is the ISO/IEC 8859-1 code (in hex as 0xXX)
    #		 Column #2 is the Unicode (in hex as 0xXXXX)
    #		 Column #3 the Unicode name (follows a comment sign, '#')
    #
    #	The entries are in ISO/IEC 8859-1 order.
    #
    #	Version history
    #	1.0 version updates 0.1 version by adding mappings for all
    #	control characters.
    #
    #	Updated versions of this file may be found in:
    #		<ftp://ftp.unicode.org/Public/MAPPINGS/>
    #
    #	Any comments or problems, contact <errata@unicode.org>
    #	Please note that <errata@unicode.org> is an archival address;
    #	notices will be checked, but do not expect an immediate response.
    #
    0x00	0x0000	#	NULL
    0x01	0x0001	#	START OF HEADING
    0x02	0x0002	#	START OF TEXT
    0x03	0x0003	#	END OF TEXT
    0x04	0x0004	#	END OF TRANSMISSION
    0x05	0x0005	#	ENQUIRY
    0x06	0x0006	#	ACKNOWLEDGE
    0x07	0x0007	#	BELL
    0x08	0x0008	#	BACKSPACE
    0x09	0x0009	#	HORIZONTAL TABULATION
    0x0A	0x000A	#	LINE FEED
    0x0B	0x000B	#	VERTICAL TABULATION
    0x0C	0x000C	#	FORM FEED
    0x0D	0x000D	#	CARRIAGE RETURN
    0x0E	0x000E	#	SHIFT OUT
    0x0F	0x000F	#	SHIFT IN
    0x10	0x0010	#	DATA LINK ESCAPE
    0x11	0x0011	#	DEVICE CONTROL ONE
    0x12	0x0012	#	DEVICE CONTROL TWO
    0x13	0x0013	#	DEVICE CONTROL THREE
    0x14	0x0014	#	DEVICE CONTROL FOUR
    0x15	0x0015	#	NEGATIVE ACKNOWLEDGE
    0x16	0x0016	#	SYNCHRONOUS IDLE
    0x17	0x0017	#	END OF TRANSMISSION BLOCK
    0x18	0x0018	#	CANCEL
    0x19	0x0019	#	END OF MEDIUM
    0x1A	0x001A	#	SUBSTITUTE
    0x1B	0x001B	#	ESCAPE
    0x1C	0x001C	#	FILE SEPARATOR
    0x1D	0x001D	#	GROUP SEPARATOR
    0x1E	0x001E	#	RECORD SEPARATOR
    0x1F	0x001F	#	UNIT SEPARATOR
    0x20	0x0020	#	SPACE
    0x21	0x0021	#	EXCLAMATION MARK
    0x22	0x0022	#	QUOTATION MARK
    0x23	0x0023	#	NUMBER SIGN
    0x24	0x0024	#	DOLLAR SIGN
    0x25	0x0025	#	PERCENT SIGN
    0x26	0x0026	#	AMPERSAND
    0x27	0x0027	#	APOSTROPHE
    0x28	0x0028	#	LEFT PARENTHESIS
    0x29	0x0029	#	RIGHT PARENTHESIS
    0x2A	0x002A	#	ASTERISK
    0x2B	0x002B	#	PLUS SIGN
    0x2C	0x002C	#	COMMA
    0x2D	0x002D	#	HYPHEN-MINUS
    0x2E	0x002E	#	FULL STOP
    0x2F	0x002F	#	SOLIDUS
    0x30	0x0030	#	DIGIT ZERO
    0x31	0x0031	#	DIGIT ONE
    0x32	0x0032	#	DIGIT TWO
    0x33	0x0033	#	DIGIT THREE
    0x34	0x0034	#	DIGIT FOUR
    0x35	0x0035	#	DIGIT FIVE
    0x36	0x0036	#	DIGIT SIX
    0x37	0x0037	#	DIGIT SEVEN
    0x38	0x0038	#	DIGIT EIGHT
    0x39	0x0039	#	DIGIT NINE
    0x3A	0x003A	#	COLON
    0x3B	0x003B	#	SEMICOLON
    0x3C	0x003C	#	LESS-THAN SIGN
    0x3D	0x003D	#	EQUALS SIGN
    0x3E	0x003E	#	GREATER-THAN SIGN
    0x3F	0x003F	#	QUESTION MARK
    0x40	0x0040	#	COMMERCIAL AT
    0x41	0x0041	#	LATIN CAPITAL LETTER A
    0x42	0x0042	#	LATIN CAPITAL LETTER B
    0x43	0x0043	#	LATIN CAPITAL LETTER C
    0x44	0x0044	#	LATIN CAPITAL LETTER D
    0x45	0x0045	#	LATIN CAPITAL LETTER E
    0x46	0x0046	#	LATIN CAPITAL LETTER F
    0x47	0x0047	#	LATIN CAPITAL LETTER G
    0x48	0x0048	#	LATIN CAPITAL LETTER H
    0x49	0x0049	#	LATIN CAPITAL LETTER I
    0x4A	0x004A	#	LATIN CAPITAL LETTER J
    0x4B	0x004B	#	LATIN CAPITAL LETTER K
    0x4C	0x004C	#	LATIN CAPITAL LETTER L
    0x4D	0x004D	#	LATIN CAPITAL LETTER M
    0x4E	0x004E	#	LATIN CAPITAL LETTER N
    0x4F	0x004F	#	LATIN CAPITAL LETTER O
    0x50	0x0050	#	LATIN CAPITAL LETTER P
    0x51	0x0051	#	LATIN CAPITAL LETTER Q
    0x52	0x0052	#	LATIN CAPITAL LETTER R
    0x53	0x0053	#	LATIN CAPITAL LETTER S
    0x54	0x0054	#	LATIN CAPITAL LETTER T
    0x55	0x0055	#	LATIN CAPITAL LETTER U
    0x56	0x0056	#	LATIN CAPITAL LETTER V
    0x57	0x0057	#	LATIN CAPITAL LETTER W
    0x58	0x0058	#	LATIN CAPITAL LETTER X
    0x59	0x0059	#	LATIN CAPITAL LETTER Y
    0x5A	0x005A	#	LATIN CAPITAL LETTER Z
    0x5B	0x005B	#	LEFT SQUARE BRACKET
    0x5C	0x005C	#	REVERSE SOLIDUS
    0x5D	0x005D	#	RIGHT SQUARE BRACKET
    0x5E	0x005E	#	CIRCUMFLEX ACCENT
    0x5F	0x005F	#	LOW LINE
    0x60	0x0060	#	GRAVE ACCENT
    0x61	0x0061	#	LATIN SMALL LETTER A
    0x62	0x0062	#	LATIN SMALL LETTER B
    0x63	0x0063	#	LATIN SMALL LETTER C
    0x64	0x0064	#	LATIN SMALL LETTER D
    0x65	0x0065	#	LATIN SMALL LETTER E
    0x66	0x0066	#	LATIN SMALL LETTER F
    0x67	0x0067	#	LATIN SMALL LETTER G
    0x68	0x0068	#	LATIN SMALL LETTER H
    0x69	0x0069	#	LATIN SMALL LETTER I
    0x6A	0x006A	#	LATIN SMALL LETTER J
    0x6B	0x006B	#	LATIN SMALL LETTER K
    0x6C	0x006C	#	LATIN SMALL LETTER L
    0x6D	0x006D	#	LATIN SMALL LETTER M
    0x6E	0x006E	#	LATIN SMALL LETTER N
    0x6F	0x006F	#	LATIN SMALL LETTER O
    0x70	0x0070	#	LATIN SMALL LETTER P
    0x71	0x0071	#	LATIN SMALL LETTER Q
    0x72	0x0072	#	LATIN SMALL LETTER R
    0x73	0x0073	#	LATIN SMALL LETTER S
    0x74	0x0074	#	LATIN SMALL LETTER T
    0x75	0x0075	#	LATIN SMALL LETTER U
    0x76	0x0076	#	LATIN SMALL LETTER V
    0x77	0x0077	#	LATIN SMALL LETTER W
    0x78	0x0078	#	LATIN SMALL LETTER X
    0x79	0x0079	#	LATIN SMALL LETTER Y
    0x7A	0x007A	#	LATIN SMALL LETTER Z
    0x7B	0x007B	#	LEFT CURLY BRACKET
    0x7C	0x007C	#	VERTICAL LINE
    0x7D	0x007D	#	RIGHT CURLY BRACKET
    0x7E	0x007E	#	TILDE
    0x7F	0x007F	#	DELETE
    0x80	0x0080	#	<control>
    0x81	0x0081	#	<control>
    0x82	0x0082	#	<control>
    0x83	0x0083	#	<control>
    0x84	0x0084	#	<control>
    0x85	0x0085	#	<control>
    0x86	0x0086	#	<control>
    0x87	0x0087	#	<control>
    0x88	0x0088	#	<control>
    0x89	0x0089	#	<control>
    0x8A	0x008A	#	<control>
    0x8B	0x008B	#	<control>
    0x8C	0x008C	#	<control>
    0x8D	0x008D	#	<control>
    0x8E	0x008E	#	<control>
    0x8F	0x008F	#	<control>
    0x90	0x0090	#	<control>
    0x91	0x0091	#	<control>
    0x92	0x0092	#	<control>
    0x93	0x0093	#	<control>
    0x94	0x0094	#	<control>
    0x95	0x0095	#	<control>
    0x96	0x0096	#	<control>
    0x97	0x0097	#	<control>
    0x98	0x0098	#	<control>
    0x99	0x0099	#	<control>
    0x9A	0x009A	#	<control>
    0x9B	0x009B	#	<control>
    0x9C	0x009C	#	<control>
    0x9D	0x009D	#	<control>
    0x9E	0x009E	#	<control>
    0x9F	0x009F	#	<control>
    0xA0	0x00A0	#	NO-BREAK SPACE
    0xA1	0x00A1	#	INVERTED EXCLAMATION MARK
    0xA2	0x00A2	#	CENT SIGN
    0xA3	0x00A3	#	POUND SIGN
    0xA4	0x00A4	#	CURRENCY SIGN
    0xA5	0x00A5	#	YEN SIGN
    0xA6	0x00A6	#	BROKEN BAR
    0xA7	0x00A7	#	SECTION SIGN
    0xA8	0x00A8	#	DIAERESIS
    0xA9	0x00A9	#	COPYRIGHT SIGN
    0xAA	0x00AA	#	FEMININE ORDINAL INDICATOR
    0xAB	0x00AB	#	LEFT-POINTING DOUBLE ANGLE QUOTATION MARK
    0xAC	0x00AC	#	NOT SIGN
    0xAD	0x00AD	#	SOFT HYPHEN
    0xAE	0x00AE	#	REGISTERED SIGN
    0xAF	0x00AF	#	MACRON
    0xB0	0x00B0	#	DEGREE SIGN
    0xB1	0x00B1	#	PLUS-MINUS SIGN
    0xB2	0x00B2	#	SUPERSCRIPT TWO
    0xB3	0x00B3	#	SUPERSCRIPT THREE
    0xB4	0x00B4	#	ACUTE ACCENT
    0xB5	0x00B5	#	MICRO SIGN
    0xB6	0x00B6	#	PILCROW SIGN
    0xB7	0x00B7	#	MIDDLE DOT
    0xB8	0x00B8	#	CEDILLA
    0xB9	0x00B9	#	SUPERSCRIPT ONE
    0xBA	0x00BA	#	MASCULINE ORDINAL INDICATOR
    0xBB	0x00BB	#	RIGHT-POINTING DOUBLE ANGLE QUOTATION MARK
    0xBC	0x00BC	#	VULGAR FRACTION ONE QUARTER
    0xBD	0x00BD	#	VULGAR FRACTION ONE HALF
    0xBE	0x00BE	#	VULGAR FRACTION THREE QUARTERS
    0xBF	0x00BF	#	INVERTED QUESTION MARK
    0xC0	0x00C0	#	LATIN CAPITAL LETTER A WITH GRAVE
    0xC1	0x00C1	#	LATIN CAPITAL LETTER A WITH ACUTE
    0xC2	0x00C2	#	LATIN CAPITAL LETTER A WITH CIRCUMFLEX
    0xC3	0x00C3	#	LATIN CAPITAL LETTER A WITH TILDE
    0xC4	0x00C4	#	LATIN CAPITAL LETTER A WITH DIAERESIS
    0xC5	0x00C5	#	LATIN CAPITAL LETTER A WITH RING ABOVE
    0xC6	0x00C6	#	LATIN CAPITAL LETTER AE
    0xC7	0x00C7	#	LATIN CAPITAL LETTER C WITH CEDILLA
    0xC8	0x00C8	#	LATIN CAPITAL LETTER E WITH GRAVE
    0xC9	0x00C9	#	LATIN CAPITAL LETTER E WITH ACUTE
    0xCA	0x00CA	#	LATIN CAPITAL LETTER E WITH CIRCUMFLEX
    0xCB	0x00CB	#	LATIN CAPITAL LETTER E WITH DIAERESIS
    0xCC	0x00CC	#	LATIN CAPITAL LETTER I WITH GRAVE
    0xCD	0x00CD	#	LATIN CAPITAL LETTER I WITH ACUTE
    0xCE	0x00CE	#	LATIN CAPITAL LETTER I WITH CIRCUMFLEX
    0xCF	0x00CF	#	LATIN CAPITAL LETTER I WITH DIAERESIS
    0xD0	0x00D0	#	LATIN CAPITAL LETTER ETH (Icelandic)
    0xD1	0x00D1	#	LATIN CAPITAL LETTER N WITH TILDE
    0xD2	0x00D2	#	LATIN CAPITAL LETTER O WITH GRAVE
    0xD3	0x00D3	#	LATIN CAPITAL LETTER O WITH ACUTE
    0xD4	0x00D4	#	LATIN CAPITAL LETTER O WITH CIRCUMFLEX
    0xD5	0x00D5	#	LATIN CAPITAL LETTER O WITH TILDE
    0xD6	0x00D6	#	LATIN CAPITAL LETTER O WITH DIAERESIS
    0xD7	0x00D7	#	MULTIPLICATION SIGN
    0xD8	0x00D8	#	LATIN CAPITAL LETTER O WITH STROKE
    0xD9	0x00D9	#	LATIN CAPITAL LETTER U WITH GRAVE
    0xDA	0x00DA	#	LATIN CAPITAL LETTER U WITH ACUTE
    0xDB	0x00DB	#	LATIN CAPITAL LETTER U WITH CIRCUMFLEX
    0xDC	0x00DC	#	LATIN CAPITAL LETTER U WITH DIAERESIS
    0xDD	0x00DD	#	LATIN CAPITAL LETTER Y WITH ACUTE
    0xDE	0x00DE	#	LATIN CAPITAL LETTER THORN (Icelandic)
    0xDF	0x00DF	#	LATIN SMALL LETTER SHARP S (German)
    0xE0	0x00E0	#	LATIN SMALL LETTER A WITH GRAVE
    0xE1	0x00E1	#	LATIN SMALL LETTER A WITH ACUTE
    0xE2	0x00E2	#	LATIN SMALL LETTER A WITH CIRCUMFLEX
    0xE3	0x00E3	#	LATIN SMALL LETTER A WITH TILDE
    0xE4	0x00E4	#	LATIN SMALL LETTER A WITH DIAERESIS
    0xE5	0x00E5	#	LATIN SMALL LETTER A WITH RING ABOVE
    0xE6	0x00E6	#	LATIN SMALL LETTER AE
    0xE7	0x00E7	#	LATIN SMALL LETTER C WITH CEDILLA
    0xE8	0x00E8	#	LATIN SMALL LETTER E WITH GRAVE
    0xE9	0x00E9	#	LATIN SMALL LETTER E WITH ACUTE
    0xEA	0x00EA	#	LATIN SMALL LETTER E WITH CIRCUMFLEX
    0xEB	0x00EB	#	LATIN SMALL LETTER E WITH DIAERESIS
    0xEC	0x00EC	#	LATIN SMALL LETTER I WITH GRAVE
    0xED	0x00ED	#	LATIN SMALL LETTER I WITH ACUTE
    0xEE	0x00EE	#	LATIN SMALL LETTER I WITH CIRCUMFLEX
    0xEF	0x00EF	#	LATIN SMALL LETTER I WITH DIAERESIS
    0xF0	0x00F0	#	LATIN SMALL LETTER ETH (Icelandic)
    0xF1	0x00F1	#	LATIN SMALL LETTER N WITH TILDE
    0xF2	0x00F2	#	LATIN SMALL LETTER O WITH GRAVE
    0xF3	0x00F3	#	LATIN SMALL LETTER O WITH ACUTE
    0xF4	0x00F4	#	LATIN SMALL LETTER O WITH CIRCUMFLEX
    0xF5	0x00F5	#	LATIN SMALL LETTER O WITH TILDE
    0xF6	0x00F6	#	LATIN SMALL LETTER O WITH DIAERESIS
    0xF7	0x00F7	#	DIVISION SIGN
    0xF8	0x00F8	#	LATIN SMALL LETTER O WITH STROKE
    0xF9	0x00F9	#	LATIN SMALL LETTER U WITH GRAVE
    0xFA	0x00FA	#	LATIN SMALL LETTER U WITH ACUTE
    0xFB	0x00FB	#	LATIN SMALL LETTER U WITH CIRCUMFLEX
    0xFC	0x00FC	#	LATIN SMALL LETTER U WITH DIAERESIS
    0xFD	0x00FD	#	LATIN SMALL LETTER Y WITH ACUTE
    0xFE	0x00FE	#	LATIN SMALL LETTER THORN (Icelandic)
    0xFF	0x00FF	#	LATIN SMALL LETTER Y WITH DIAERESIS
    "
    write(fileData, outputFile)
    close(outputFile)
  }
}
