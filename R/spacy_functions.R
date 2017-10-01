
#' Parse sentences strings using SPACY as the parsing engine, and prep sentences for running through topic model
#'
#' @param workingfolder folder to save file in.
#' @param CORPUS_A corpus to analyze
#' @param sample_num number of rows to sample default is 0 which keeps all rows.
#' @param removeentitites should named entities be dropped from topic model
#' @param dbexits whether the spacy or sqlnet database exists
#' @param spcy use spacy to parse (set to True usually)
#' @param syntaxnet use syntaxnet to parse (set to False usually)
#' @return will create a database named spacyframe.db
#' @seealso the entire dplyr package 
#' @export
#' @description This command creates a sqlite database which parse information is saved in, while also creating the proper frames for running a topic model.

SPCY_PreTopicFrame<-function(CORPUS_A,sample_num=0,workingfolder,removeentities=T,dbexists=FALSE,spcy=T,syntaxnet=F){
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
    spcydb<-if(dbexists==F){
      dplyr::src_sqlite(file.path(workingfolder,'spacyframe.db'),create=T)} else {dplyr::src_sqlite(file.path(workingfolder,'spacyframe.db'),create=F)}
    spin<-spacyr::spacy_initialize()
    spinout<-spacyr::spacy_parse(pstrings$string,dependency=T,named_entity=T,full_parse = TRUE)
    spacyr::spacy_finalize()
    spinout$Orig<-fullorig$Orig[as.numeric(gsub("text","",spinout$docname))]
    dplyr::copy_to(spcydb,spinout,'parses_uw',temporary=F)
  }
  processed <-stm::textProcessor(fullorig$string,metadata=select(fullorig,-string),sparselevel=1)                                     
  out <- stm::prepDocuments(processed$documents,processed$vocab,processed$meta,lower.thresh=4)
  colnames(out$meta)<-gsub("\\W",".",colnames(out$meta))
  saveRDS(list("SentFrame"=fullorig,"Annotations"="please call textDB or spacyframe.db to access","processed"=processed,"out"=out),"PTF2.rds")
}



#' Access database for spacy parses.
#'
#' @param workingfolder folder to save file in.
#' @return will load a database names spacyframe.db
#' @seealso the entire dplyr package 
#' @export
#' @description This command creates a sqlite database which parse information is saved in.
#' @examples
#' accessSPCDB()
#' Create an access database for spacy parses.
#'
#' @param workingfolder folder to save file in.
#' @return will create a database names spacyframe.db
#' @seealso the entire dplyr package 
#' @export
#' @description This command creates a sqlite database which parse information is saved in.
#' @examples
#' accessSPCDB()
accessSPCDB<-function(workingfolder){
  spcydb<-dplyr::src_sqlite(file.path(workingfolder,"spacyframe.db"))
}

#' Create a tmidf matrix for keyword extraction
#'
#' @param CORPUS_A corpus to analyze
#' @param db database to save tfidf to
#' @return will create a database names spacyframe.db
#' @seealso the entire dplyr package 
#' @export
#' @description This command creates a sqlite database which parse information is saved in.
#' @examples
#' accessSPCDB()
maketermkeywords<-function(CORPUS_A,db){
  dtmidf <- tm::DocumentTermMatrix(CORPUS_A, control = list(weighting = weightTfIdf))
  dtmidf2<-as.matrix(dtmidf)
  colnames(dtmidf2)<-Terms(dtmidf)
  dtmidf2<-reshape2::melt(dtmidf2)
  colnames(dtmidf2)<-c("Orig","tokens","Weight")
  copy_to(db, dtmidf2,"tfidf",temporary=F)
}

#' Load a wordnet sqldatabase
#'
#' @export
#' @description This command loads a wordnet database named sqlite2.db
#' @examples
#' accessSPCDB()
loadwordnet<-function(){
  tempdir<-list.files(path="~",pattern="sqlite-31.db",full.names=T)
  if(length(tempdir)==0){if(menu(c('yes','no'),title="Okay to download and install wordnet db 650mb?")==1){
    download.file('https://sourceforge.net/projects/sqlunet/files/4.0.0/sqlite/sqlite-4.0.0-31-all.zip/download',"~/wndbf.zip")
    unzip("~/wndbf.zip",exdir="~")
    file.remove("~/wndbf.zip")
    tempdir<-list.files(path="~",pattern="sqlite2.db",full.names=T)}}
  dplyr::src_sqlite(tempdir,create=F)
}

#' Select relevant frames to keep
#'
#' @export
#' @description This command uses framenet to identify verb frames
#' @param VERBWORD regular expression for verb you would like to generally keep
#' @param parsecnnl table resulting from spacy parse 
#' @examples
#' IDframes("[cause|change],collect(tbl(scydb,"parses_uw"),n=Inf)")
IDframes<-function(VERBWORD,parsecnnl){
  wn<-loadwordnet()
  filter(tbl(wn, "fnwords"), word%in%filter(parsecnnl,pos=="VERB")$tokens) %>% left_join(tbl(wn,"fnlexemes")) %>% left_join(tbl(wn, "fnlexunits")) %>% left_join(tbl(wn, "fnframes")) %>% select(.,word,frame,framedefinition) %>% collect %>% filter(., stringr::str_detect(tolower(frame), VERBWORD))
}

#' Keep only sentences from relevant frames
#'
#' @export
#' @param parsecnnl table resulting from spacy parse 
#' @description This command uses the result of IDframes and keeps only sentences with relevant verb frames

sentkeeper<-function(idframes,parsecnnl,database=T){
  if(database==T){
    parsecnnl<-mutate(parsecnnl,"word_out"=tolower(token))
    keepers<-dplyr::filter(parsecnnl, word_out%in%idframes$word,pos=="VERB") %>% select(doc_id)
    dplyr::filter(parsecnnl, doc_id%in%collect(keepers)$doc_id) %>% collect(n=Inf)
    
  }  else {
    keepers<-dplyr::filter(parsecnnl, tolower(token)%in%idframes$word,pos=="VERB")$doc_id
    dplyr::filter(parsecnnl, doc_id%in%keepers)
  }
}


#' Keep only sentences from relevant frames
#'
#' @export
#' @param parseconnll spacy parsed frame
#' @param tmidf tfidf frame
#' @param scydb db with spacy in it
#' @description This command joins a parseconll and term matrix for a database.
join_parse_tmidf<-function(parseconnll, scydb){
  left_join(parseconnll,tbl(scydb,"tfidf"),copy=T)
  copy_to(scydb,parseconnll,"merged_wcuts",temporary=T)
  tbl(scydb,"merged_wcuts")
}

#' It makes a network graph of words. First part of list is a data frame showing connections, Second part I believe is an actual igraph object.
#'
#' @export
#' @param parseconnll spacy parsed frame
#' @param tmidf tfidf frame
#' @param scydb db with spacy in it
#' @description This command makes a network graph of words from a spacy output frame.
conligraph<-function(parseconll){
  spl<-lapply(unique(parseconll$doc_id),function(X) subset(parseconll, doc_id==X))
  library(igraph)
  splout<-lapply(spl,function(splt) {
    g1a<-data.frame("from"=splt$token, "to"=splt$token[sapply(splt$head_token_id,function(K) which(splt$token_id==K))])
    g1<-g1a %>% igraph:: graph_from_data_frame(directed=TRUE) 
    commst<-igraph::walktrap.community(g1)
    tg<-igraph::contract(g1,mapping=igraph::membership(igraph::walktrap.community(g1)),vertex.attr.comb=toString)
    # tokenizers::tokenize_word_stems(stopwords=tokenizers::stopwords("en")) 
    V(tg)$name <-V(tg)$name %>% tolower() %>% lapply(.,function(X) paste(unique(X),collapse=","))
    commst<-commst %>% membership() %>% as.matrix()
    splt2<-left_join(splt,data.frame("token"=names(commst[,1]),"community"=commst[,1]))
    splt2$words_out<-splt2$token
    list(tg,splt2)})
  list("frame"=lapply(splout,function(X) X[[2]]) %>% bind_rows(),"nets"=lapply(splout,function(X) X[[1]]))
}

