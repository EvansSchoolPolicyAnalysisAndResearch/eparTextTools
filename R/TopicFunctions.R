#TopicFunctions

#' Clean the java environment.
#'
#' Function cleans the java environment. This is a great command. Run it whenever you feel like it.
#' @seealso \code{\link{rJava}} 
#' @examples
#' jgc()
jgc <- function(){
  gc()
  rJava::.jcall("java/lang/System", method = "gc")
} 

#' Loads nlp parsers.
#'
#'This function loads the parsers for openNLP
#' @examples
#' loadparsers()
loadparsers<-function(){
sent_token_annotator <<- openNLP::Maxent_Sent_Token_Annotator()
word_token_annotator <<- openNLP::Maxent_Word_Token_Annotator()
pos_tag_annotator<<- openNLP::Maxent_POS_Tag_Annotator()
org.annotate<<-openNLP::Maxent_Entity_Annotator(language = "en", kind="organization", probs = FALSE,model = NULL)
pers.annotate<<-openNLP::Maxent_Entity_Annotator(language = "en", kind="person", probs = FALSE,model = NULL)
location.annotate<<-openNLP::Maxent_Entity_Annotator(language = "en", kind="location", probs = FALSE,model = NULL)
#money.annotate<<-openNLP::Maxent_Entity_Annotator(language = "en", kind="money", probs = FALSE,model = NULL)
parse_annotator <<- openNLP::Parse_Annotator()
t1<-gc()
jgc()
gc()-t1}

#'Process a corpus into a topic model ready object
#'
#'  This function takes a corpus and creates a processed version of that corpus with entities removed for use in a topic model. Additionally it allows you to specify common entities to count across documents for use as a covariate in the topic model. The object it returns includes a frame of the text, an annotation object, a processed version of the corpus with stems and stopwords removed, and an out object which is the input object for fitting a topic model within the stm package.
#' @param CORPUS_A Document corpus
#' @param howmanyentities Count of entities you would like to add as potential covariates for topic model
#' @return SentFrame data frame with one row for each paragraph chunk
#' @return Annotations openNLP annotation object
#' @return processed stm processed documents
#' @return out stm out documents for use in topic model
#' @seealso \code{\link{stm}} 
#' @export
#' @examples
#' BASE_INPUT<-PreTopicFrame(corpus1,1)
PreTopicFrame<-function(CORPUS_A,howmanyentities=10){
  loadparsers()
  corpEntN<-CORPUS_A[sapply(CORPUS_A,function(x) length(content(x)))>0]
  corpus1b<-CORPUS_A[sapply(CORPUS_A,function(x) length(content(x)))>0]
  corpEntN<-lapply(names(corpEntN[[1]]$meta),function(K){meta(corpEntN,K)})
  corpEntN[[1]]<-melt(sapply(1:length(corpEntN[[1]]),function(x) paste(corpEntN[[1]][[x]],collapse="")))[,1]
  corpEntN[[2]]<-melt(sapply(1:length(corpEntN[[2]]), function(x) lubridate::ymd_hms(corpEntN[[2]][[x]]),simplify=T))
  corpEntN[[2]]<-join(data.frame("L1"=1:length(corpEntN[[1]])),corpEntN[[2]])$value
  corpEntN[[3]]<-melt(sapply(1:length(corpEntN[[3]]),function(x) paste(corpEntN[[3]][[x]],collapse="")))[,1]
  corpEntN[[4]]<-melt(sapply(1:length(corpEntN[[4]]),function(x) paste(corpEntN[[4]][[x]],collapse="")))[,1]
  corpEntN[[5]]<-melt(sapply(1:length(corpEntN[[5]]),function(x) paste(corpEntN[[5]][[x]],collapse="")))[,1]
  corpEntN[[6]]<-melt(sapply(1:length(corpEntN[[6]]),function(x) paste(corpEntN[[6]][[x]],collapse="")))[,1]
  corpEntN[[7]]<-melt(sapply(1:length(corpEntN[[6]]),function(x) paste(corpEntN[[7]][[x]],collapse="")))[,1]
  corpEntN<-data.frame(corpEntN)
  colnames(corpEntN)<-names(meta(CORPUS_A[[1]]))
  para_token_annotator <-Annotator(function(s, a = Annotation()) {
    spans <- blankline_tokenizer(s)
    n <- length(spans)
    ## Need n consecutive ids, starting with the next "free"
    ## one:
    from <- next_id(a$id)
    Annotation(seq(from = from, length.out = n),
               rep.int("paragraph", n),
               spans$start,
               spans$end)}, list(description ="A paragraph token annotator based on blankline_tokenizer()."))
  par1<-lapply(corpus1b,function(x) NLP::annotate(as.String(x$content),para_token_annotator))
  jgc()
  corpEntN$Orig<-1:nrow(corpEntN)
  par2<-lapply(1:length(par1),function(i) {
    dfout<-join(data.frame("S1"=sapply(par1[i],function(x) as.character(as.String(content(corpus1b[[i]]))[x]),USE.NAMES=FALSE),"Orig"=i),corpEntN[i,],by="Orig")
    colnames(dfout)[1]<-"S"
    dfout})
  jgc()
  par2<-do.call(rbind,par2)
  par2$SC<-as.character(par2$S)
  par2<-par2[nchar(par2$SC)>=100,]
  jgc()
  allans<-pbapply::pblapply(par2$S,function(X) NLP::annotate(X,sent_token_annotator))
  par2<-par2[which(sapply(allans,function(X) length(X)>0)==TRUE),]
  jgc()
  allans<-allans[sapply(allans,function(X) length(X)>0)]
  allans<-pbapply::pblapply(1:nrow(par2),function(i){NLP::annotate(par2$S[i],word_token_annotator,allans[[i]])})
  jgc()
  allans<-pbapply::pblapply(1:nrow(par2),function(i){NLP::annotate(par2$S[i],list(org.annotate,pers.annotate,location.annotate),allans[[i]])})
  jgc()
  par2$SnE<-NA
  par2$ents<-NA
  for(i in 1:nrow(par2)){
    entit1<-c(as.String(par2$S[i])[subset(allans[[i]], type=="entity")])
    par2$ents[i]<-paste(entit1,sep="",collapse=",")
    sn<-as.character(par2$S[i])
    if(length(entit1)>0){
      for(k in 1:length(entit1)){
        sn<-gsub(entit1[k],"",sn,fixed=TRUE)
      }
      par2$SnE[i]<-sn}}
  allents<-sapply(par2$ents,function(x) c(unlist(strsplit(x,","))),USE.NAMES=FALSE)
  allents<-unlist(allents)
  top.ents<-names(sort(table(allents),decreasing=TRUE)[1:howmanyentities])
  entp<-lapply(top.ents, function(X) str_detect(par2$SC,fixed(X)))
  names(entp)<-top.ents
  entp<-do.call(cbind,entp)
  processed <-stm::textProcessor(par2$SnE,metadata=cbind(par2[,2:length(par2)],entp),sparselevel=1)
  out <- stm::prepDocuments(processed$documents,processed$vocab,processed$meta,lower.thresh=4)
  colnames(out$meta)<-gsub("\\W",".",colnames(out$meta))
  list("SentFrame"=par2,"Annotations"=allans,"processed"=processed,"out"=out)
}



#'Process a corpus into a topic model ready object, but do so while using the syntaxnet or spacy Parser rather than the openNLP parser.
#'
#'  This function takes a corpus and creates a processed version of that corpus with entities removed for use in a topic model. Additionally it allows you to specify common entities to count across documents for use as a covariate in the topic model. The object it returns includes a frame of the text, an annotation object, a processed version of the corpus with stems and stopwords removed, and an out object which is the input object for fitting a topic model within the stm package.
#' @param CORPUS_A Document corpus
#' @param sample_num If sampling from documents, how many chunks to sample from each document? If 0, all chunks are kept.
#' @param removeentities TRUE removes all proper nouns. FALSE keeps all proper nouns. "ONLY" will keep only proper nouns. This means you could fit a topic model to only proper nouns contained in documents
#' @param spcy Defaults to true. Uses spacy to create a database.
#' @param syntaxnet Defaults to true but if it is false, this skips annotating for NLP and goes directly to text processing. This might be useful if you dont care about limiting verb-phrases in later steps and want to fit a topic model to all documents.
#' @param workingfolder folder where files are to be saved.
#' @return SentFrame data frame with one row for each paragraph chunk
#' @return Annotations Rather than being stored in memory, annotations are stored in a sqllite database in the workingfolder in a file named textDB. This saves on memory. For each document, annotations are stored in a table named according to the document id. There are three additional tables, one of just propeer nouns called propers, one of nonproper nouns called nonprops, and one with everything called fulltable.
#' @return processed stm processed documents
#' @return out stm out documents for use in topic model
#' @seealso \code{\link{stm}} 
#' @export
#' @examples
#' BASE_INPUT<-PreTopicFrame2(CORPUS_A,0,syntaxnet=F,workingfolder)
PreTopicFrame2<-function(CORPUS_A,sample_num=0,workingfolder,removeentities=T,dbexists=FALSE,spcy=T,syntaxnet=F){
  cstrings<-lapply(CORPUS_A,function(X) paste(NLP::content(X),collapse="\n\n"))
  pstrings<-lapply(cstrings,tokenizers::tokenize_sentences,simplify=T)
  pstrings<-dplyr::bind_rows(lapply(pstrings,function(X) data.frame("string"=X)),.id="Orig")
  dca<-as.data.frame(do.call(cbind,lapply(names(CORPUS_A[[1]]$meta),function(K) sapply(CORPUS_A, function(X) paste(NLP::meta(X,K),collapse=";"),USE.NAMES=FALSE))))
  colnames(dca)<-names(NLP::meta(CORPUS_A[[1]]))
  dca$Orig<-dca$id
  row.names(dca)<-1:nrow(dca)
  fullorig<-merge(pstrings,dca,by.x="Orig")
  if(syntaxnet==T){
    if(dbexists==FALSE){SQLtabSyntaxNET(fullorig,sample_num)}
    my_db<-src_sqlite(file.path(workingfolder,"textDB"),create=F)
    if(removeentities==T){
      temptab<-ddply(collect(tbl(my_db, sql("SELECT Orig,Sent,V2 FROM nonpropers")),n=Inf),.(Orig,Sent),text=paste(V2,collapse=" "),summarize)
    } else {if(removeentities=="ONLY"){temptab<-ddply(collect(tbl(my_db, sql("SELECT Orig,Sent,V2 FROM propers")),n=Inf),.(Orig,Sent),text=paste(V2,collapse=" "),summarize)}
      else {
        temptab<-ddply(collect(tbl(my_db, sql("SELECT Orig,Sent,V2 FROM fulltable")),n=Inf),.(Orig,Sent),text=paste(V2,collapse=" "),summarize)
      }}
    processed <-stm::textProcessor(temptab$text,metadata=select(merge(temptab,dca,by.x="Orig"),-text),sparselevel=1)                                         
    out <- stm::prepDocuments(processed$documents,processed$vocab,processed$meta,lower.thresh=4)
  } 
  if(spcy==T){
    dtmidf <- tm::DocumentTermMatrix(CORPUS_A, control = list(weighting = tm::weightTfIdf))
    spcydb<-if(dbexists==F){
      dplyr::src_sqlite(file.path(workingfolder,'spacyframe.db'),create=T)} else {dplyr::src_sqlite(file.path(workingfolder,'spacyframe.db'),create=F)}
    spin<-spacyr::spacy_initialize()
    spinout<-spacyr::spacy_parse(pstrings$string,dependency=T,named_entity=T,full_parse = TRUE)
    spacyr::spacy_finalize()
    spinout$Orig<-fullorig$Orig[as.numeric(gsub("text","",spinout$doc_id))]
    dplyr::copy_to(spcydb,spinout,'parses_uw',temporary=F)
  }
  processed <-stm::textProcessor(fullorig$string,metadata=select(fullorig,-string),sparselevel=1)                                     
  out <- stm::prepDocuments(processed$documents,processed$vocab,processed$meta,lower.thresh=4)
  colnames(out$meta)<-gsub("\\W",".",colnames(out$meta))
  list("SentFrame"=fullorig,"Annotations"="please call textDB or spacyframe.db to access","processed"=processed,"out"=out)
}



#'Keep only those sentences with desired verbs.
#'
#'  This function parses sentences only keeping those sentences that include a verb that is in the list which is specified by the user. More verbs may be kept if they occur in the sentence alongside the desired verb.
#' @param WT a ilst of verbs
#' @param PROCESSED an stm processed object.
#' @return OUT an stm prepped object.
#' @return ANNOTATELIST openNLP annotation object
#' @return SENTENCEFRAME SentFrame object returned by PreTopicFrame
#' @return toptopics TopTopics extracted from an stm topic model.
#' @seealso \code{\link{stm}} 
#' @export
#' @examples
#' AnnotatesLarge<-AnnotateVerbsTopicJoin(allwords,BASE_INPUT$processed,BASE_INPUT$out,BASE_INPUT$Annotations,BASE_INPUT$SentFrame,BASE_INPUT$top.topics)
AnnotateVerbsTopicJoin<-function(WT,PROCESSED,OUT,ANNOTATELIST,SENTENCEFRAME,toptopics){
  loadparsers()
  pos_tag_annotator<- Maxent_POS_Tag_Annotator()
  #OUT<-BASE_INPUT$out
  #ANNOTATELIST<-BASE_INPUT$Annotations
  #PROCESSED<-BASE_INPUT$processed
  #WT<-allwords
    toChar<-SENTENCEFRAME$S[-PROCESSED$docs.removed] %>% .[-OUT$docs.removed]
    toCharANS<-ANNOTATELIST[-PROCESSED$docs.removed] %>% .[-OUT$docs.removed]
    KeepA<-unique(unlist(sapply(WT, function(W) which(str_detect(toChar,tolower(W))))))
      anns<-lapply(KeepA,function(i){
      ANT<-annotate(as.String(toChar[i]),pos_tag_annotator,toCharANS[i][[1]])
      wans<-subset(ANT, type=="word")
      vans<-wans[which(unlist(wans$features)%in%c("VB","VBN","VBZ","VBG","VBN","VBP"))]
      vanstf<-sapply(1:length(vans), function(i2){TRUE%in%str_detect(as.character(as.String(toChar[i])[vans[i2]]),WT)})
      vans2<-vans[vanstf==TRUE]
      vanstf<-do.call(rbind,lapply(1:length(vans2), function(i3){
        vsent<-subset(ANT,start<=vans2[i3]$start) %>% subset(.,end>=vans2[i3]$end) %>% subset(type=="sentence")
        #vsent<-subset(anns[[1]], start>=vsent$start) %>% subset(.,end<=vsent$end)
        #psent<-parse_annotator(as.String(toChar),vsent)
        #list(vsent,psent)
        as.String(toChar[i])[vsent]
        data.frame("Sent"=ifelse(length(vsent)>0,as.character(as.String(toChar[i])[vsent]),NA),"rowid"=i)
      }))
      vanstf})
    anf<-do.call(rbind,anns)  
    anf<-cbind(anf,OUT$meta[anf$rowid,])
    anf$top.topics<-toptopics[anf$rowid]
    anf<-subset(anf, is.na(anf$Sent)==FALSE)
    anf<-unique(anf)
    anf}

#'Keep only those annotations with tagged verbs
#'
#' Combines sentences into blocks reducing the size of the dataframe and only keeping those that contain desired verbs. Think of it like a refilter from AnnotateVerbsTopicJoin
#' @param PREPFRAME the frame object regurned by AnnotateVerbsTopicJoin
CombinationFrame<-function(PREPFRAME){
  ddply(.data=PREPFRAME, .(Orig,rowid),summarize,"Sent"=paste(as.character(Sent),collapse="\n",sep="\n"))
}


#'Write the formula for the call to the topic model.
#'
#' Write the formula for the prevelence covariates for the topic model.
#' @param BASE_INPUT the baseinput object returned by PreTopicFrame
#' @param workingfolder the workingfolder where you are saving objects.
#' @return writes the formula to a text file in the working formula named formula1.
#' @seealso \code{\link{stm}} 
#' @export
#' @examples
#' writeFormulaforSTM(BASE_INPUT,"../Research.Grants")
writeFormulaforSTM<-function(BASE_INPUT,workingfolder){
  form1<-paste("~as.factor(Orig)",paste(select.list(colnames(BASE_INPUT$out$meta),multiple=TRUE),sep="",collapse="+"),sep="+")
  writeLines(form1,file.path(workingfolder,"formula1.txt"))
}

#'Run STM in a seperate R session.
#'
#' this function will create .R file in the working folder which it then will call, resulting in fitting an stm model within a seperate R session. The output should be viewable in an Rout file.
#' @param workingfolder the workingfolder where you are saving objects.
#' @param filename name of the preprocessed rds file from pretopicprocess.
#' @return writes file and starts new r process in background
#' @seealso \code{\link{stm}} 
#' @export
#' @examples
#' runSTM("../Research.Grants")
runSTM<-function(workingfolder,filename){
writeLines(text=paste('
library(plyr)
library(dplyr)
library(stm)
args = commandArgs(TRUE)
workingfolder<-args[1]
baseinput<-readRDS(file.path(workingfolder,','"',filename,'"','))
  st1<-stm(baseinput$out$documents,baseinput$out$vocab,data=baseinput$out$meta,prevalence=eval(parse(text=readLines(file.path(workingfolder,"formula1.txt")))),K=0, init.type="Spectral",max.em.its=1000)
saveRDS(st1, file.path(workingfolder,"topicmodel.rds"))',sep=""),con=file.path(workingfolder,"callstm.R"))
  system(paste("R CMD BATCH --no-restore ","'--args ",workingfolder,"' ",file.path(workingfolder,"callstm.R"),sep=""),wait=FALSE)}

AnnotateVerbsByTopic<-function(MAXTOPS,WT,PROCESSED,OUT,ANNOTATELIST,SENTENCEFRAME){
  loadparsers()
  pos_tag_annotator<- Maxent_POS_Tag_Annotator()
  tcs<-pblapply(1:length(unique(MAXTOPS)), function(k){
    toChar<-SENTENCEFRAME$SC[-PROCESSED$docs.removed][-OUT$docs.removed][which(MAXTOPS==k)]
    toCharANS<-ANNOTATELIST[-PROCESSED$docs.removed] %>% .[-OUT$docs.removed] %>% .[which(MAXTOPS==k)]
    KeepA<-unique(unlist(sapply(WT, function(W) which(str_detect(toChar,tolower(W))))))
    anns<-lapply(KeepA,function(i){
      ANT<-annotate(as.String(toChar[i]),pos_tag_annotator,toCharANS[i][[1]])
      wans<-subset(ANT, type=="word")
      vans<-wans[which(unlist(wans$features)%in%c("VB","VBN","VBZ","VBG","VBN","VBP"))]
      vanstf<-sapply(1:length(vans), function(i2) TRUE%in%str_detect(as.character(as.String(toChar[i])[vans[i2]]),WT))
      droppedLs<-which(vanstf==FALSE)
      vans2<-vans[vanstf==TRUE]
      vanstf<-do.call(rbind,lapply(1:length(vans2), function(i3){
        vsent<-subset(ANT,start<=vans2[i3]$start) %>% subset(.,end>=vans2[i3]$end) %>% subset(type=="sentence")
        #vsent<-subset(anns[[1]], start>=vsent$start) %>% subset(.,end<=vsent$end)
        #psent<-parse_annotator(as.String(toChar),vsent)
        #list(vsent,psent)
        data.frame("Sent"=ifelse(length(vsent)>0,as.character(as.String(toChar[i])[vsent]),NA),"rowid"=which(SENTENCEFRAME[-PROCESSED$docs.removed,][-OUT$docs.removed,][which(MAXTOPS==k),][i,]$SC==SENTENCEFRAME[-PROCESSED$docs.removed,][-OUT$docs.removed,]$SC))
      }))
      vanstf})})
  names(tcs)<-c(1:MAXTOPS)[-droppedLs]
  droppedLs2<-which(sapply(tcs,length)<1)
  tcs2<-tcs[sapply(tcs, length)>0]
  allcs<-na.omit(melt(tcs2))
  allcs}

#'Clean sentences for sending to alchemy API.
#'
#' This function cleans the sentences for sending to alchemy API
#' @param PREPFRAME object returned by combinationframe or annotateframebytopicjoin
#' @return a data frame
#' @export
ProcessforAPI<-function(PREPFRAME){
PREPFRAME$Sent<-str_trim(PREPFRAME$Sent)
PREPFRAME$Sent<-gsub("\n"," ",PREPFRAME$Sent)
badcs<-paste(gsub("\\w|\\s|\\.|\\,","",PREPFRAME$Sent),collapse="")
badcs<-unique(unlist(strsplit(badcs,split="")))
for(i in 1:length(badcs)){
  PREPFRAME$Sent<-gsub(badcs[i]," ",PREPFRAME$Sent,fixed=T)
}
PREPFRAME$Sent<-gsub("  "," ",PREPFRAME$Sent)
PREPFRAME}



#'Identify topics with highest term probability
#'
#' taking a topic model, this function returns the top topics for that word.
#' @param STMOBJ stm object
#' @param TERM  word in character format
#' @return N number of topics to return
#' @export
#' @examples idtopics(st1,"research",10)
idtopics<-function(STMOBJ,TERM,N){
  #TERM<-"cassava"
  #STMOBJ<-st1
  ptops<-findTopic(STMOBJ,list(TERM),n=N,type="prob")
  if(length(ptops)>1){
  ptopsl<-sapply(1:length(ptops), function(i) paste(sageLabels(STMOBJ,n=N)$marginal$prob[ptops,][i,],collapse=" "))
  mselect<-ptops[which(ptopsl%in%select.list(ptopsl,multiple=T))]
  } else {mselect<-ptops}
  cat(paste("\nmodel",mselect,"selected",sep=" ",collapse=","))
  return(mselect)
}

#'Pick topics with suggested words
#'
#' Pulls similar words, expanding topic selection ability. Interactive.
#' @param STMOBJ stm object
#' @param TERM  word in character format
#' @return N number of topics to return
#' @export
#' @examples idtopics(st1,"research",10)
picktopics<-function(STMOBJ,TERM,N){
  test<-idtopics(STMOBJ,TERM,N)
  termlist<-TERM
  if(menu(c("yes","no"),title="\nCheck For Similar Words?")==1){
    terms <- getIndexTerms("NOUN", 1, getTermFilter("ExactMatchFilter", TERM, TRUE))
    newterms<-getSynonyms(terms[[1]]) %>% select.list(.,multiple=TRUE)
    termlist<-c(termlist,newterms)
    test<-append(test,sapply(newterms,function(X) tryCatch({idtopics(STMOBJ,X,N)},error=function(e){NA})))
    test<-na.omit(test)
  }
  cat("\nYou are using topics\n",test,"\nbased on selection of",termlist,sep=" ",collapse=" ")
  list("searchterms"=termlist,"topics"=test)
  }
  
#'Make a map based on the cliff-claven tagging
#'
#' Using the countrypredictions object, retuns a map for an opportunity number.
#' @param CountryPredictions object returned by PredictCountryByDoc
#' @param OPPORTUNITY OpID to search for
#' @return data frame ready for mapping in plotly
data_mapper<-function(CountryPredictions,OPPORTUNITY){
  gchars<-ddply(CountryPredictions,.(OpID),summarise,"charsum"=sum(nchars))
  CountryPredictions<-plyr::join(CountryPredictions,gchars)
  CountryPredictions$weight<-CountryPredictions$nchars/CountryPredictions$charsum
  CountryPredictions<-cbind(CountryPredictions[,1:2],CountryPredictions[,3:c(ncol(CountryPredictions)-3)]*CountryPredictions$weight)
  ftemp<-dplyr::filter(CountryPredictions, OpID%in%OPPORTUNITY)
  tframe<-reshape2::melt(colSums(ftemp[,3:ncol(ftemp)],na.rm=TRUE)/sum(ftemp[,3:ncol(ftemp)],na.rm=TRUE))
  tframe$countryids<-row.names(tframe)
  tframe$ccode2<-plyr::join(data.frame("ISO3166.1.Alpha.2"=tframe$countryids),ccodes)$ISO3166.1.Alpha.3
  tframe$nameC<-plyr::join(data.frame("ISO3166.1.Alpha.2"=tframe$countryids),ccodes)$official_name_en
  tframe$hover<-paste(tframe$nameC,": ",round(tframe$value*100)/100,sep="")
  tframe}


FillFolderMan<-function(PREPFRAME,FOLDERNAME){
  library(httr)
  #if(dir.exists("getAlchemy")==FALSE) {dir.create("getAlchemy")}
  #if(dir.exists(file.path(FOLDERNAME,"ALCHEMY"))==FALSE) {dir.create(file.path(FOLDERNAME,"ALCHEMY"))}
  for(i in 1:nrow(PREPFRAME)){
    X<-PREPFRAME$Sent[i]
    req <- POST("http://access.alchemyapi.com/calls/text/TextGetRelations", 
                body = list(apikey="6837e8ae18678cadd3c42fc55ed938b3818ce470",text=X,keywords=1,outputMode="json",disambiguate=0), encode = "form",write_disk(file.path("getAlchemy",FOLDERNAME,paste("al",i,".json",sep=""))))
    Sys.sleep(.5)
  }}

findBroken<-function(FOLDERNAME){
flist<-list.files(file.path(FOLDERNAME,"ALCHEMY"))
bad<-sapply(flist,function(X){
tfile<-readLines(file.path("getAlchemy",FOLDERNAME,X))
p1<-str_detect(paste(tfile,collapse=" ",sep=" "),"If you are seeing this message, you are likely making an excessive number of concurrent HTTP connections to this service.  Please check the concurrency limits for your assigned service tier.")
a<-tryCatch(jsonlite::fromJSON(tfile), error=function(e){"error"})
p2<-a=="error"
TRUE%in%c(p1,p2)
})        
flist[bad]
}

FillFolder<-function(PREPFRAME,FOLDERNAME){
  library(httr)
  if(dir.exists(file.path(FOLDERNAME,"ALCHEMY"))==FALSE) {dir.create(file.path(FOLDERNAME,"ALCHEMY"))}
  for(i in args[1]:nrow(PREPFRAME)){
    X<-PREPFRAME$Sent[i]
    req <- POST("http://access.alchemyapi.com/calls/text/TextGetRelations", 
                body = list(apikey=ALKEY,text=X,keywords=1,outputMode="json",disambiguate=0), encode = "form",write_disk(file.path("getAlchemy",FOLDERNAME,paste("al",i,".txt",sep=""))))
    Sys.sleep(1)
  }}

FillRetry<-function(PREPFRAME,FOLDERNAME,FILENAMES){
  library(httr)
  for(X in FILENAMES){
    Num1<-gsub(".json","",X) %>% gsub("al","",.) %>% as.numeric()
    Sent<-PREPFRAME$Sent[Num1]
    file.remove(file.path("getAlchemy",FOLDERNAME,X))
    req <- POST("http://access.alchemyapi.com/calls/text/TextGetRelations", 
                body = list(apikey="6837e8ae18678cadd3c42fc55ed938b3818ce470",text=Sent,keywords=1,outputMode="json",disambiguate=0), encode = "form",write_disk(file.path(FOLDERNAME,"ALCHEMY",X)))
    Sys.sleep(.5)
    
  }}


#'Process Documents in the Alchemy Folder
#'
#' This function parses the entries in the alchemy api folder and creates a dataframe
#' @param workingfolder Folder for working in
#' @param PREPFRAME prepframe object
#' @return data frame of processed sentences
processFolder<-function(workingfolder, PREPFRAME){
  tempfname<-list.files(file.path(workingfolder,"ALCHEMY"),full.names=T)
  nosize<-sapply(tempfname, function(X) file.info(X)$size==0)
  templist<-pblapply(tempfname[which(nosize==FALSE)], function(X){fromJSON(X,flatten=T)$relations %>% as.data.frame() %>% mutate(.,"filename"=as.character(X))})
  tempframe<-rbind.pages(templist)
  tempframe$comboID=str_extract(tempframe$filename,"[0-9]+")
  PREPFRAME$comboID<-1:nrow(PREPFRAME)
  recombineOut<-plyr::join(tempframe,PREPFRAME, by="comboID",type="left",match="all")
  return(recombineOut)
}


ParseFolderToFrame<-function(FOLDERNAME,PREPFRAME,WT){
  rels<-lapply(1:length(list.files(file.path(FOLDERNAME,"ALCHEMY"))),function(i) fromJSON(readLines(file.path(FOLDERNAME,"ALCHEMY",list.files(file.path(FOLDERNAME,"ALCHEMY"))[i])))$relations)
  rels2<-rels[which(sapply(rels,function(x) length(x)>0))]
  joinkey<-lapply(1:length(rels),function(i)  {
    mrow<-max(c(nrow(rels[[i]]$subject$keywords[[1]]),nrow(rels[[i]]$object$keywords[[1]]),length(rels[[i]]$action$verb$text)))
    emp<-data.frame(matrix(rep(NA,mrow*4),ncol=4))
    #lapply(1:4,function(q) emp[,q]<-as.character(emp[,q]))
    colnames(emp)<-c("subkey","obkey","verbkey","whichi")
    emp$subkey<-as.character(emp$subkey)
    emp$obkey<-as.character(emp$obkey)
    emp$verbkey<-as.character(emp$verbkey)
    emp$whichi<-as.character(emp$whichi)
    if(length(rels[[i]])>0){
      if(is.null(rels[[i]]$subject$keywords[[1]])==FALSE){
        emp$subkey[1:nrow(rels[[i]]$subject$keywords[[1]])]<-unlist(rels[[i]]$subject$keywords[[1]])
      }
      if(is.null(rels[[i]]$object$keywords[[1]])==FALSE){
        emp$obkey[1:nrow(rels[[i]]$object$keywords[[1]])]<-unlist(rels[[i]]$object$keywords[[1]])
      }
      if(is.null(rels[[i]]$action$verb$text)==FALSE){
        emp$verbkey[1:length(rels[[i]]$action$lemmatized)]<-rels[[i]]$action$verb$text
      }
      emp$whichi<-i
    }
    emp
  })
  joinkey<-do.call(rbind.fill,joinkey)
  joinkey$subkey<-tolower(joinkey$subkey)
  joinkey$obkey<-tolower(joinkey$obkey)
  joinkey$verbkey<-tolower(joinkey$verbkey)
  joinkey_sub<-subset(joinkey, verbkey%in%WT)
  library(dplyr)
  subcs<-cbind(joinkey_sub,PREPFRAME[as.numeric(list.files(file.path(FOLDERNAME,"ALCHEMY"))[as.numeric(joinkey_sub$whichi)] %>% gsub("al","",.) %>% gsub('.json',"",.,fixed=T)),])
  subcs_full<-cbind(joinkey,PREPFRAME[as.numeric(list.files(file.path(FOLDERNAME,"ALCHEMY"))[as.numeric(joinkey$whichi)] %>% gsub("al","",.) %>% gsub('.json',"",.,fixed=T)),])
  list("SmallVerb"=subcs,"FullVerb"=subcs_full)}


frametable1<-function(PARSEFRAME,BASEINPUT,FOLDERNAME,PREPFRAME){
  OV<-PREPFRAME[as.numeric(list.files(file.path(FOLDERNAME,"ALCHEMY")) %>% gsub(".json","",.) %>% gsub("al","",.) %>% .[as.numeric(PARSEFRAME$whichi)]),]
  OVF<-BASEINPUT$SentFrame[-BASEINPUT$processed$docs.removed,][-BASEINPUT$out$docs.removed,][OV$value,] 
  PARSEFRAME<-cbind(PARSEFRAME,OVF[,c("author","datetimestamp","id","ents")])
  PARSEFRAME<-PARSEFRAME[,c(1,2,3,5,7,8,9,10,11)]
  colnames(PARSEFRAME)<-c('Subject','Object',"Verb","Sentence","Topic","Author","Date","Document","Entities")
  PARSEFRAME$Topic<-sapply(1:10,function(i) paste(labelTopics(st1,n=5)[[1]][i,],collapse=","))[PARSEFRAME$Topic]
  PARSEFRAME$Topic<-as.factor(PARSEFRAME$Topic)
  PARSEFRAME$Verb<-as.factor(PARSEFRAME$Verb)
  PARSEFRAME}

frametable2<-function(PARSEFRAME,BASEINPUT,FOLDERNAME,PREPFRAME,TOPICMOD){
  # PARSEFRAME<-Frame1[[1]]
  # rm(PARSEFRAME)
  # BASEINPUT<-BASE_INPUT
  # rm(BASEINPUT)
  # FOLDERNAME<-"test2"
  # rm(FOLDERNAME)
  # PREPFRAME<-Pf1
  # rm(PREPFRAME)
  # TOPICMOD<-st1
  # rm(TOPICMOD)
  OV<-PREPFRAME[as.numeric(list.files(file.path(FOLDERNAME,"ALCHEMY")) %>% gsub(".json","",.) %>% gsub("al","",.) %>% .[as.numeric(PARSEFRAME$whichi)]),]
  OVF<-BASEINPUT$SentFrame[-BASEINPUT$processed$docs.removed,][-BASEINPUT$out$docs.removed,][OV$value,] 
  PARSEFRAME<-cbind(PARSEFRAME,OVF[,c("author","datetimestamp","id","ents")])
  colnames(PARSEFRAME)
  PARSEFRAME<-PARSEFRAME[,c("subkey","obkey","verbkey","Sent","L1","author","datetimestamp",'id',"ents")]
  colnames(PARSEFRAME)<-c('Subject','Object',"Verb","Sentence","Topic","Author","Date","Document","Entities")
  PARSEFRAME$Topic<-sapply(1:TOPICMOD$settings$dim$K,function(i) paste(labelTopics(TOPICMOD,n=5)[[1]][i,],collapse=","))[PARSEFRAME$Topic]
  PARSEFRAME$Topic<-as.factor(PARSEFRAME$Topic)
  PARSEFRAME$Verb<-as.factor(PARSEFRAME$Verb)
  PARSEFRAME}


frametable.html<-function(PARSEFRAME,BASEINPUT,FOLDERNAME,PREPFRAME,TOPICMOD){
  OV<-PREPFRAME[as.numeric(list.files(file.path(FOLDERNAME,"ALCHEMY")) %>% gsub(".json","",.) %>% gsub("al","",.) %>% .[as.numeric(PARSEFRAME$whichi)]),]
  OVF<-BASEINPUT$SentFrame[-BASEINPUT$processed$docs.removed,][-BASEINPUT$out$docs.removed,][OV$value,] 
  PARSEFRAME<-cbind(PARSEFRAME,OVF[,c("author","datetimestamp","id","ents")])
  PARSEFRAME<-PARSEFRAME[,c(1,2,3,5,7,8,9,10,11)]
  colnames(PARSEFRAME)<-c('Subject Keyword','Object Keyword',"Verb","Sentence","Topic","Author","Date and Time","Document","Entities")
  PARSEFRAME$Topic<-sapply(1:TOPICMOD$settings$dim$K,function(i) paste(labelTopics(TOPICMOD,n=5)[[1]][i,],collapse=","))[PARSEFRAME$Topic]
  PARSEFRAME$Topic<-as.factor(PARSEFRAME$Topic)
  PARSEFRAME$Verb<-as.factor(PARSEFRAME$Verb)
  datatable(data=PARSEFRAME,rownames=FALSE,filter="top")}


TopicCoreFrame<-function(MAXTOPS,PARSEFRAME){
  library(stringdist)
  cormatsbysub<-lapply(unique(MAXTOPS),function(X) {
    topiccs<-filter(PARSEFRAME,L1==X)
    if(nrow(topiccs)>0){
    dmh<-hclust(stringdistmatrix(unique(na.omit(topiccs$subkey)),method="lcs"))
    dmh$labels<-unique(na.omit(topiccs$subkey))
    ids1<-rect.hclust(dmh,h=mean(dmh$height))
    dmo<-hclust(stringdistmatrix(unique(na.omit(topiccs$obkey)),method="lcs"))
    dmo$labels<-unique(na.omit(topiccs$obkey))
    ids2<-rect.hclust(dmo,h=min(dmo$height))
    dmv<-hclust(stringdistmatrix(unique(na.omit(topiccs$verbkey)),method="lcs"))
    dmv$labels<-unique(na.omit(topiccs$verbkey))
    
    bysubject<-lapply(ids1,function(Y) topiccs[which(topiccs$subkey%in%names(Y)),])
    names(bysubject)<-sapply(ids1,function(Y) names(Y[1]))
    head(bysubject)
    byobject<-lapply(ids2,function(Y) topiccs[which(topiccs$obkey%in%names(Y)),])
    names(byobject)<-sapply(ids2,function(Y) names(Y[1]))
    
    byverb<-lapply(unique(na.omit(topiccs$verbkey)),function(Y) PARSEFRAME[which(topiccs$verbkey==Y),])
    names(byverb)<-unique(na.omit(topiccs$verbkey))
    
    combine.ds<-function(bysubject,byobject,byverb){
      byobject2<-unique(melt(byobject)[,-7])
      byverb2<-unique(melt(byverb))
      bysubject2<-unique(melt(bysubject)[,-7])
      head(bysubject2)
      colnames(byverb2)[5]<-"VerbCat"
      colnames(byobject2)[7]<-"ObjCat"
      colnames(bysubject2)[7]<-"SubCat"
      merge(bysubject2[,c(4,7)],byverb2[,c(4,5)], all=T) %>% merge(byobject2[,c(4,7)],all=T)
    }
    obcats<-combine.ds(bysubject,byobject,byverb)
    head(obcats)
    tst<-obcats[,c("SubCat","ObjCat")]
    ts.e<-matrix(rep(NA,length(unique(tst$ObjCat))*length(unique(tst$SubCat))),nrow=length(unique(tst$SubCat)))
    rownames(ts.e)<-unique(tst$SubCat)
    colnames(ts.e)<-unique(tst$ObjCat)
    for(i in 1:nrow(ts.e)){
      ts.e[i,]<-sapply(names(ts.e[i,]),function(K){length(which(subset(tst,tst$SubCat==rownames(ts.e)[i])$ObjCat==K))})
    }
    ts.e} else {matrix(NA)}
  })
  cormatsbysub
}


SimpFrameList<-function(TCOREFRAME){
  lapply(TCOREFRAME,function(X){t1<-melt(as.matrix(X)) %>% filter(value>0)
  colnames(t1)<-c("Subject","Object","Count")
  t1
  })}



clusterWords<-function(COLUMN,NUM,TOPICMODEL){
  subprobs<-do.call(rbind.fill,lapply(unique(COLUMN), function(X){
    k<-data.frame(t(exp(TOPICMODEL$beta$logbeta[[1]][,which(TOPICMODEL$vocab==tm::stemDocument(tolower(X)))])))
    if(nrow(k)>0){
      k$word<-X
    }
    k
  }))
  row.names(subprobs)<-subprobs$word
  subprobs<-subprobs[,1:10]
  d1<-kmeans(subprobs,NUM)
  optName<-sapply(1:NUM, function(k) {names(d1$cluster[d1$cluster==k][which.min(sapply(1:length(which(d1$cluster==k)),function(i) sum(abs(subprobs[which(d1$cluster==k),]-fitted(d1)[which(d1$cluster==k),])[i,])))])})
  data.frame("Word"=row.names(subprobs),"optName"=optName[d1$cluster])
}

nearest_to2<-function(topicm,wordvec,n=10,fixword=FALSE,limitwords=NULL){ 
  mt<-t(exp(topicm$beta$logbeta[[1]]))
  row.names(mt)<-topicm$vocab
  if(fixword==FALSE){
    whichwords<-unique(unlist(sapply(wordvec,function(X){which(stringr::str_detect(row.names(mt),X))})))} else
    {whichwords<-which(row.names(mt)%in%wordvec)}
  if(length(limitwords)>0){
    whichwords<-whichwords[limitwords]
  }
  sims = cosineSimilarity(mt, matrix(as.vector(mt[whichwords,]), ncol = ncol(mt)))
  ords = order(-sims[, 1])
  list("searchwords"=row.names(mt)[whichwords],"simwords"=structure(1 - sims[ords[1:n]], names = rownames(sims)[ords[1:n]]))
}

idSimilar<-function(SEARCH,COLUMN,NUM,TOPICMODEL){
  subprobs<-do.call(rbind.fill,lapply(unique(COLUMN), function(X){
    k<-data.frame(t(exp(TOPICMODEL$beta$logbeta[[1]][,which(TOPICMODEL$vocab==tm::stemDocument(tolower(X)))])))
    if(nrow(k)>0){
      k$word<-X
    }
    k
  }))
  sid<-which(subprobs$word==SEARCH)
  row.names(subprobs)<-subprobs$word
  subprobs<-subprobs[,1:10]
  d1<-kmeans(subprobs,NUM)
  names(d1$cluster[d1$cluster==d1$cluster[sid]])
}


#nex<-read.csv("../nonExcel.csv",stringsAsFactors=FALSE)
#nexjoin<-plyr::join(data.frame("name"=BASE_INPUT$SentFrame$id),data.frame("name"=basename(as.character(nex$path)),"OpID"=as.character(nex$Opportunity.ID)),type="left",match="first")
#BASE_INPUT$SentFrame$OpID<-nexjoin$OpID



#'Function to write code to run alchemy
#'
#' This function creates a .R file to call to alchemy api.
#' @param alchemykey key to alchemy api as character
#' @param workingfolder workingfolder for project
#' @return a .R file saved in the working folder
#' @export
writerun_alch<-function(alchemykey,workingfolder){
writeLines(text=c(paste('
rm(list=ls())
library(httr)
library(plyr)
library(dplyr)
library(stringr)
library(jsonlite)
ALKEY<-',alchemykey,sep=""),paste('workingfolder<-','"',file.path(workingfolder),'"',sep=""),'args = commandArgs(trailingOnly=TRUE)
FillFolder<-function(PREPFRAME,FOLDERNAME){
library(httr)
if(dir.exists(file.path(FOLDERNAME,"ALCHEMY"))==FALSE) {dir.create(file.path(FOLDERNAME,"ALCHEMY"))}
for(i in args[1]:nrow(PREPFRAME)){
X<-PREPFRAME$Sent[i]
req <- POST("http://access.alchemyapi.com/calls/text/TextGetRelations", 
body = list(apikey=ALKEY,text=X,keywords=1,outputMode="json",disambiguate=0), encode = "form",write_disk(file.path(FOLDERNAME,"ALCHEMY",paste("al",i,".txt",sep=""))))
Sys.sleep(1)
}}
recombine<-readRDS(file.path(workingfolder,"ProcessedFrame.rds"))
FillFolder(recombine,workingfolder)'),con=file.path(workingfolder,"run_Alchemy.R"))
}

#'Function to run .R file returned by writerun_alch
#'
#' This function runs the .R file returned by writerun_alch. It takes as an argument a number which corresponds to a row in the BASE_INPUT dataframe
#' @param num row number to start function at (this is to allow for restarts in case an error occurs?)
#' @return NA fills folder with alchemy api calls
#' @export
RunAlchy<-function(num){
  system(paste("R CMD BATCH --args",file.path(workingfolder,"run_Alchemy.R"),num),wait=FALSE)
}



#'Run this function to install a docker server from within R
#'
#' This function calls out using the system command. make sure to clear your workspace and history if you use this. However, what a cool function, right?
#' @param password you should type in your password 3 times.
#' @return NA installs a docker server
#' @export
makeadockercliff<-function(){
  cwd<-getwd()
  setwd("~")
  system('sudo apt-get -y install docker.io',input=readline("Enter your password: "))
  system('sudo git clone https://github.com/johnb30/cliff-docker
',input=readline("Enter your password: "))
  setwd("cliff-docker")
  system('sudo docker build -t cliff:2.1.1 .',input=readline("Enter your password: "))
  setwd("~")
  setwd(cwd)}


#'These commands start and run a cliff-docker server. 
#'
#'These commands collectively build, ensure a server is working, and stop a running docker server for gecoding. Once you have run buildcliff and startcliff, you can run the PredictCountryByDoc command to predict which countries documents correspond to. For these commands, we need to have the opportunity id labels as OpID in the SentFrame part of BASE_INPUT
#'@export
buildcliff<-function() {system('sudo docker run -p "8080:8080" -d --name cliff cliff:2.1.1')}
startcliff<-function() {system('sudo docker start cliff')}
checkcliff<-function(){system('sudo docker ps')}
stopcliff<-function(){system('sudo docker stop cliff')}
PredictCountryByDoc<-function(BASE_INPUT){
  fullc<-ddply(BASE_INPUT$SentFrame, .(Orig), summarise, "fullc"=paste(SC, collapse=" ",sep=" "))
  countries<-vector("list",nrow(fullc))
  for(i in 1:nrow(fullc)){
    try({
      TEXTI<-fullc$fullc[i]
      ncons<-ceiling(nchar(TEXTI)/ceiling(nchar(TEXTI)/4000))
      noutstop<-c(1,c(1:ceiling(nchar(TEXTI)/4000))*ncons)
      stsp<-lapply(1:c(length(noutstop)-1),function(X) c(noutstop[X],noutstop[X+1]))
      subers<-lapply(stsp,function(X) substr(TEXTI,X[1],X[2]))
      mres<-lapply(subers,function(X) httr::GET(url="http://localhost:8080/CLIFF-2.1.1/parse/text",query = list(replaceAllDemonyms="true",q=X)))
      mres<-mres[sapply(mres,function(X) X$status_code==200)]
      res<-do.call(rbind.fill, lapply(mres,function(X) jsonlite::fromJSON(X[[1]][[1]])$results$places$focus$countries))[,c("countryCode","score")]
      res<-unlist(sapply(1:nrow(res),function(i) rep(res$countryCode[i],res$score[i])))
      countries[[i]]<-data.frame(t(as.matrix(table(res)/length(res))))
      cat(i)})}
  for(i in which(sapply(countries,length)==0)){
    countries[[i]]=data.frame("none"=0)}
  countries<-do.call(rbind.fill,countries)
  countries$Orig<-fullc$Orig
  countries$nchars<-nchar(fullc$fullc)
  f1<-join(BASE_INPUT$SentFrame[,c("OpID","Orig")],countries)
  f1<-unique(f1)
  f1
}
reflectCountryCol<-function(MATCHFRAME,pred2,howmany,binomial=FALSE){
  pred2[,3:c(ncol(pred2)-1)][is.na(pred2[,3:c(ncol(pred2)-1)])]<-0
  if(binomial==TRUE){
    pred2[,3:c(ncol(pred2)-1)][pred2[,3:c(ncol(pred2)-1)]>0]<-1
  }
  maxCs<-colnames(pred2[,3:c(ncol(pred2)-1)])[sort(colSums(pred2[,3:c(ncol(pred2)-1)]),index=TRUE,decreasing=TRUE)$ix][1:howmany]
  cntjoin<-join(MATCHFRAME,pred2[,c("OpID","Orig",maxCs)], by=c("OpID","Orig"))
  return(cntjoin)}
runMap<-function(FILENAMEQUOTED,path.file=FALSE,titleC){
  shiny::shinyApp(
    ui = fluidPage(sidebarLayout(
      sidebarPanel(selectizeInput("OpporID", label = "Opportunity Number", choices = unique(pred2$OpID),multiple=TRUE, selected=pred2$OpID[1]),
                   selectInput("projection",label="Map Projection",choices=c('equirectangular','mercator', 'orthographic','natural earth','kavrayskiy7','miller', 'robinson','eckert4','azimuthal equal area','azimuthal equidistant','conic equal area','conic conformal','conic equidistant','gnomonic','stereographic','mollweide','hammer','transverse mercator'),selected='mollweide')
      ),
      mainPanel(
        plotlyOutput("plotly"),
        tableOutput("table")
      ))),
    server = function(input, output) {
      ccodes<<-read.csv(textConnection(RCurl::getURL("http://data.okfn.org/data/core/country-codes/r/country-codes.csv")),stringsAsFactors=FALSE)
      output$plotly<-renderPlotly({
        g <- list(showframe = FALSE,showcoastlines = TRUE,projection = list(type = input$projection))
        l <- list(color = toRGB("grey"), width = 0.5)
        d1<-data_mapper(pred2,input$OpporID)
        plot_ly(d1,z=value,hoverinfo="text", locations = ccode2, type = 'choropleth', marker = list(line = l),color = value, colors = 'Blues',zmax=max(value),zmin=min(value),text=hover,colorbar = list(title = 'Country Relevance'),source="select") %>% layout(title =paste(titleC,"<br>Source:<a href='https://github.com/c4fcm/CLIFF'>CLIFF</a>"), geo = g)
      })
      output$table<-renderTable({
        eventdata<-event_data("plotly_click", source = "select")$pointNumber+1
        d1<-data_mapper(pred2,input$OpporID)
        rframe<-data.frame("col1"=pred2$OpID[sort.int(pred2[,d1$countryids[eventdata]],decreasing=TRUE,index.return=TRUE)$ix[1:5]])
        colnames(rframe)[1]<-paste("Top 5 Opportunities for", d1$nameC[eventdata])
        rframe})
    },
    onStart=function(){
      suppressWarnings(library(plyr,warn.conflicts=FALSE,quietly=TRUE))
      suppressWarnings(library(dplyr,warn.conflicts=FALSE,quietly=TRUE))
      suppressWarnings(library(plotly,warn.conflicts=FALSE,quietly=TRUE))
      suppressWarnings(library(RCurl,warn.conflicts=FALSE,quietly=TRUE))
      suppressWarnings(library(DT,warn.conflicts=FALSE,quietly=TRUE))
      pred2<<-if(path.file==TRUE){read.csv(FILENAMEQUOTED,stringsAsFactors=FALSE) %>% .[,2:ncol(.)]} else {FILENAMEQUOTED}
      data_mapper<<-function(CountryPredictions,OPPORTUNITY){
        gchars<-ddply(CountryPredictions,.(OpID),summarise,"charsum"=sum(nchars))
        CountryPredictions<-plyr::join(CountryPredictions,gchars)
        CountryPredictions$weight<-CountryPredictions$nchars/CountryPredictions$charsum
        CountryPredictions<-cbind(CountryPredictions[,1:2],CountryPredictions[,3:c(ncol(CountryPredictions)-3)]*CountryPredictions$weight)
        ftemp<-dplyr::filter(CountryPredictions, OpID%in%OPPORTUNITY)
        tframe<-reshape2::melt(colSums(ftemp[,3:ncol(ftemp)],na.rm=TRUE)/sum(ftemp[,3:ncol(ftemp)],na.rm=TRUE))
        tframe$countryids<-row.names(tframe)
        tframe$ccode2<-plyr::join(data.frame("ISO3166.1.Alpha.2"=tframe$countryids),ccodes)$ISO3166.1.Alpha.3
        tframe$nameC<-plyr::join(data.frame("ISO3166.1.Alpha.2"=tframe$countryids),ccodes)$official_name_en
        tframe$hover<-paste(tframe$nameC,": ",round(tframe$value*100)/100,sep="")
        tframe}
    })}


nearest_to2<-function(topicm,wordvec,n=10,fixword=FALSE,limitwords=NULL){ 
  mt<-t(exp(topicm$beta$logbeta[[1]]))
  row.names(mt)<-topicm$vocab
  if(fixword==FALSE){
    whichwords<-unique(unlist(sapply(wordvec,function(X){which(stringr::str_detect(row.names(mt),X))})))} else
    {whichwords<-which(row.names(mt)%in%wordvec)}
  if(length(limitwords)>0){
    whichwords<-whichwords[limitwords]
  }
  sims = cosineSimilarity(mt, matrix(as.vector(mt[whichwords,]), ncol = ncol(mt)))
  ords = order(-sims[, 1])
  list("searchwords"=row.names(mt)[whichwords],"simwords"=structure(1 - sims[ords[1:n]], names = rownames(sims)[ords[1:n]]))
}



tableapp<-function(basic_table,TOPICMODEL){
  wordpossibles<-colnames(basic_table)[unlist(sapply(c("subject","action","object"),function(X) which(stringr::str_detect(colnames(basic_table),X))))]
  basic_table<-basic_table %>% select(-one_of("Sent","Orig"))
  for(i in which(sapply(1:ncol(basic_table), function(X) class(basic_table[,X]))=="list")){
    basic_table[,i]<-sapply(basic_table[,i],function(X){paste(unlist(X),collapse=";")})}
  basic_table$sentence<-stringr::str_trim(basic_table$sentence)
  shinyApp(ui = fluidPage(sidebarLayout(
    sidebarPanel(selectInput("NOV", label = "Select Word Type",choices = wordpossibles, selected = wordpossibles[1]),
                 textInput("word",label= "Word", value=""),
                 checkboxInput("sent",label="Show Sentence",value=FALSE),
                 sliderInput("Kclusters", label = "Number of Words",min = 1, max = 25, value = 5, step = 1),
                 checkboxInput('Shorten',label='Shorten Strings?', value=TRUE)),
    mainPanel(dataTableOutput("table1")))), 
    server = function(input, output) {
      output$table1<-renderDataTable({
        if(input$sent==TRUE){t1<-basic_table} else {t1<-basic_table %>% select(-one_of("sentence"))}
        if(nchar(input$word)>1){
          matchwords<-names(nearest_to2(TOPICMODEL,input$word,n=input$Kclusters,fixword=FALSE,limitwords=NULL)[[2]])
          t1<-t1[unlist(sapply(matchwords,function(X) {which(stringr::str_detect(t1[,input$NOV],X))})),]}
        else {t1<-t1} 
        if(input$Shorten==TRUE){
          datatable(t1,rownames=FALSE,filter ='bottom',extensions = 'Buttons')} else {datatable(t1,rownames=FALSE,filter ='bottom',extensions = 'Buttons', options = list(dom = 'Bfrtip',scrollX=TRUE,buttons = c('copy', 'csv', 'excel', 'pdf', 'print','colvis')), callback = JS('table.page(3).draw(false);'))}
      })})}


igraphob_object_force<-function(WORD,mtable,W,inputWord=TRUE,sankey=FALSE,verbfilter=c(),strictlimit=FALSE){
  # WORD='production'
  #mtable<-matchtable
  #W<-10
  if(inputWord==TRUE){
    ntr<-dplyr::filter(mtable,str_detect(tolower(subject.keywords),tolower(WORD)))
  }
  else {
    ntr<-dplyr::filter(mtable,str_detect(tolower(object.keywords),tolower(WORD)))
  }
  dt <- data.table(ntr)
  dt$object.keywords<-sapply(dt$object.keywords,function(X) paste(unlist(X),collapse=";"))
  dt$subject.keywords<-sapply(dt$subject.keywords,function(X) paste(unlist(X),collapse=";"))
  dt$action.lemmatized<-sapply(dt$action.lemmatized,function(X) paste(unlist(X),collapse=";"))
  dt$action.verb.text<-sapply(dt$action.verb.text,function(X) paste(unlist(X),collapse=";"))
  dt<-dt[,list(action.verb.text = unlist(strsplit(action.verb.text,";")),object.keywords=unlist(strsplit(object.keywords,";")),subject.keywords=unlist(strsplit(subject.keywords,";"))),by = sentence]
  dt<-na.omit(dt)
  if(length(verbfilter)>0){
    dt<-dplyr::filter(dt,action.verb.text%in%verbfilter)
  }
  if(strictlimit==TRUE){
    if(inputWord==TRUE){
      dt<-dplyr::filter(dt,str_detect(tolower(subject.keywords),tolower(WORD)))} else {dt<-dplyr::filter(dt,str_detect(tolower(object.keywords),tolower(WORD)))}
  }
  net1<-graph_from_data_frame(rbind(data.frame("In"=dt$subject.keywords,"Out"=dt$action.verb.text),data.frame("In"=dt$action.verb.text,"Out"=dt$object.keywords)))
  E(net1)$weight <- 1
  netsimp<-simplify(net1,edge.attr.comb=list(weight="sum","ignore"))
  E(netsimp)$width <- E(netsimp)$weight
  netsimp <- delete_edges(netsimp, E(netsimp)[weight<=W])
  netsimp<-simplify(netsimp,edge.attr.comb=list(weight="sum","ignore"))
  bad.vs<-V(netsimp)[degree(netsimp) == 0]
  netsimp <-delete.vertices(netsimp, bad.vs)
  d3ob<-netsimp %>% igraph_to_networkD3() 
  d3ob$nodes$group=sapply(d3ob$nodes$name,function(X) {if(X%in%dt$action.verb.text) {"verb"} else {if(X%in%dt$subject.keyword){"subject"} else {if(X%in%dt$object.keywords){"object"}}}})
  if(sankey==FALSE){
    forceNetwork(d3ob$links,d3ob$nodes,"source","target","value",NodeID="name",Group="group",fontSize=24,fontFamily="Arial", opacity = 0.9, bounded = TRUE, opacityNoHover = TRUE,colourScale=JS("d3.scale.ordinal().domain(['verb','subject','object']).range(['#d9d9d9','#b7a57a','#4b2e83'])"))} else {
      sankeyNetwork(d3ob$links,d3ob$nodes,"source","target","value",NodeID="name",NodeGroup="group",colourScale=JS("d3.scale.ordinal().domain(['verb','subject','object']).range(['#d9d9d9','#b7a57a','#4b2e83'])"),fontSize=24,fontFamily="Arial")}}


frametable<-function(PARSEFRAME,BASEINPUT,origent){
  basejoin<-BASEINPUT$out$meta
  basejoin$TopTopics<-BASEINPUT$top.topics
  joined<-cbind(PARSEFRAME,basejoin[PARSEFRAME$rowid,])
  colnames(joined)
  joined<-joined[,c(c(1:ncol(PARSEFRAME)),18,19,22,27,c(28+origent),ncol(joined))]
  joined<-joined[,c(1:ncol(joined))[-which(colnames(joined)%in%c("filename","comboID","rowid"))]]
  joined$ents<-gsub(",",";",joined$ents)
  joined
}

#'Function to run syntaxnet
#'
#' This function calls syntaxnet.
#' @export
SQLtabSyntaxNET<-function(TEXTIN,sample_num){
  TEXTIN$string<-stringi::stri_trans_general(TEXTIN$string, "latin-ascii")
  TEXTIN$string<-iconv(TEXTIN$string, "latin1", "ASCII",sub='') 
  TEXTIN$string<-gsub('[0-9]+', '', TEXTIN$string)
  TEXTIN$string<-gsub("\\s+"," ",TEXTIN$string)
  TEXTIN<-dplyr::filter(TEXTIN,nchar(TEXTIN$string)>25)
my_db<-dplyr::src_sqlite(file.path(workingfolder,"textDB"),create=T)
l_ply(unique(fullorig$Orig)[unique(fullorig$Orig)%in%src_tbls(my_db)==FALSE]
,function(X) {
  dft<-filter(fullorig,Orig==X)
  dft<-if(sample_num>0){sample_n(dft,sample_num)}
  copy_to(my_db,SNetTXTCALL(dft$string),X, temporary = FALSE)
  cat(paste("...", X))
  })
fulltable<-lapply(src_tbls(my_db)[which(src_tbls(my_db)%in%unique(fullorig$Orig))],function(X) collect(dplyr::tbl(my_db,X))) %>% bind_rows(.id="Orig")
head(fulltable)
fulltable$Orig<-src_tbls(my_db)[which(src_tbls(my_db)%in%unique(fullorig$Orig))][as.numeric(fulltable$Orig)]
propers<-filter(fulltable,V5%in%c("NNP","NNPS"))
nonpropers<-filter(fulltable,V5%in%c("NNP","NNPS")==FALSE,V4%in%c(".")==FALSE,tolower(V2)%in%stopwords("en")==FALSE)
copy_to(my_db,fulltable,"fulltable", temporary = FALSE)
copy_to(my_db,propers,"propers", temporary = FALSE)
copy_to(my_db,nonpropers,"nonpropers", temporary = FALSE)
cat(paste("tables saved to",file.path(workingfolder,"textDB")))
}

#'Function to write code to run syntaxnet
#'
#' This function calls syntaxnet.
#' @export
SNetTXTCALL<-function(doc,dontprintoutput=T){
  doc<-gsub("'","’",doc)
  writeLines(doc,"~/input.txt")
  file.create("~/output.txt")
  system(paste("cd; cd models/syntaxnet; syntaxnet/demo2.sh",sep=" "),intern=F,ignore.stderr = dontprintoutput)
  view_out<-lapply(unlist(strsplit(readr::read_file("~/output.txt"),split="\n\n")),function(X){
    temp<-read.delim(textConnection(X),header=F)
    closeAllConnections()
    temp}) %>% bind_rows(.id="Sent")
  view_out
}
