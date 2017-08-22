PreTopicFrame2<-function(CORPUS_A,sample_num=0,syntaxnet=T,workingfolder,removeentities=T,dbexists=FALSE){
  #CORPUS_A<-corpus1
  #sample_num=0
  #syntaxnet=T
  #removeentities=T
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
  } else {
    cat("WARNING: Skipping annotation annotation file will be empty")
    processed <-stm::textProcessor(fullorig$string,metadata=select(fullorig,-string),sparselevel=1)                                         
    out <- stm::prepDocuments(processed$documents,processed$vocab,processed$meta,lower.thresh=4)
  }
  colnames(out$meta)<-gsub("\\W",".",colnames(out$meta))
  list("SentFrame"=fullorig,"Annotations"="please call textDB to access","processed"=processed,"out"=out)
}


BASE_INPUT<-PreTopicFrame2(corpus1,sample_num=0,dbexists=FALSE,workingfolder=workingfolder)

dplyr::src_tbls(my_db)
tbl(my_db, sql("SELECT * FROM fulltable"),n=10)
