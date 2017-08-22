#where corpus1 is entering corpus
workingfolder="Research.Grants"

corpus1<-readRDS(file.path(workingfolder,"corpus.rds"))
PreTopicFrame2<-function(CORPUS_A,sample_num=0,workingfolder,removeentities=T,dbexists=FALSE,spcy=T,syntaxnet=F){
  CORPUS_A<-corpus1
  #sample_num=0
  syntaxnet=F
  spcy=T
  dbexists=F
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
  } 
  if(spcy==T){
    dtmidf <- tm::DocumentTermMatrix(CORPUS_A, control = list(weighting = weightTfIdf))
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
  list("SentFrame"=fullorig,"Annotations"="please call textDB or spacyframe.db to access","processed"=processed,"out"=out)
}

rm(spinout)

trial1<-PreTopicFrame2(corpus1,sample_num=0,workingfolder,dbexists=FALSE,spcy=T)


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

scydb<-accessSPCDB(workingfolder)

maketermkeywords<-function(CORPUS_A,db){
dtmidf <- tm::DocumentTermMatrix(CORPUS_A, control = list(weighting = weightTfIdf))
dtmidf2<-as.matrix(dtmidf)
colnames(dtmidf2)<-Terms(dtmidf)
dtmidf2<-reshape2::melt(dtmidf2)
colnames(dtmidf2)<-c("Orig","tokens","Weight")
copy_to(db, dtmidf2,"tfidf",temporary=F)
}
maketermkeywords(corpus1,scydb)

loadwordnet<-function(){
  tempdir<-list.files(path="~",pattern="sqlite2.db",full.names=T)
  if(length(tempdir)==0){if(menu(c('yes','no'),title="Okay to download and install wordnet db 650mb?")==1){
    download.file('https://sourceforge.net/projects/sqlunet/files/4.0.0/sqlite/sqlite-4.0.0-31-all.zip/download',"~/wndbf.db")
    tempdir<-list.files(path="~",pattern="sqlite2.db",full.names=T)}}
  dplyr::src_sqlite(tempdir)
}

wn<-loadwordnet()

IDframes<-function(VERBWORD,parsecnnl){
  wn<-loadwordnet()
  filter(tbl(wn, "fnwords"), word%in%filter(parsecnnl,pos=="VERB")$token) %>% left_join(tbl(wn,"fnlexemes")) %>% left_join(tbl(wn, "fnlexunits")) %>% left_join(tbl(wn, "fnframes")) %>% select(.,word,frame,framedefinition) %>% collect %>% filter(., stringr::str_detect(tolower(frame), VERBWORD))
}

idframes<-IDframes("[change|cause]", collect(tbl(scydb,"parses_uw"),n=Inf))

sentkeeper<-function(idframes,parsecnnl,database=T){
  if(database==T){
    parsecnnl<-mutate(parsecnnl,"word_out"=tolower(tokens))
    keepers<-dplyr::filter(parsecnnl, word_out%in%idframes$word,pos=="VERB") %>% select(docname)
    dplyr::filter(parsecnnl, docname%in%collect(keepers)$docname) %>% collect(n=Inf)
    
  }  else {
    keepers<-dplyr::filter(parsecnnl, tolower(tokens)%in%idframes$word,pos=="VERB")$docname
    dplyr::filter(parsecnnl, docname%in%keepers)
  }
}

spinsub<-sentkeeper(idframes,tbl(scydb,"parses_uw"),database=T)

spinsub<-left_join(spinsub,tbl(scydb,"tfidf"),copy=T)

join_parse_tmidf<-function(parseconnll, scydb){
  left_join(parseconnll,tbl(scydb,"tfidf"),copy=T)
  copy_to(scydb,parseconnll,"merged_wcuts",temporary=T)
  tbl(scydb,"merged_wcuts")
}

copy_to(scydb,spinsub,"merged_wcuts",temporary=T)

spinsub<-tbl(scydb,"merged_wcuts")

docs<-select(tbl(scydb,"merged_wcuts"),Orig) %>% distinct() %>% collect()

conligraph<-function(parseconll){
  spl<-lapply(unique(parseconll$docname),function(X) subset(parseconll, docname==X))
  library(igraph)
  splout<-lapply(spl,function(splt) {
    g1a<-data.frame("from"=splt$tokens, "to"=splt$tokens[sapply(splt$head_id,function(K) which(splt$id==K))])
    g1<-g1a %>% igraph:: graph_from_data_frame(directed=TRUE) 
    commst<-igraph::walktrap.community(g1)
    tg<-igraph::contract(g1,mapping=igraph::membership(igraph::walktrap.community(g1)),vertex.attr.comb=toString)
    # tokenizers::tokenize_word_stems(stopwords=tokenizers::stopwords("en")) 
    V(tg)$name <-V(tg)$name %>% tolower() %>% lapply(.,function(X) paste(unique(X),collapse=","))
    commst<-commst %>% membership() %>% as.matrix()
    splt2<-left_join(splt,data.frame("tokens"=names(commst[,1]),"community"=commst[,1]))
    splt2$words_out<-splt2$tokens
    list(tg,splt2)})
  list("frame"=lapply(splout,function(X) X[[2]]) %>% bind_rows(),"nets"=lapply(splout,function(X) X[[1]]))
}


library(igraph)
chg<-conligraph(collect(dplyr::filter(spinsub, Orig==docs$Orig[2])))

library(plyr)
library(dplyr)

nets.to.ig.d3<-function(chg,cutprob){
  #chg$frame<-subset(chg$frame, chg$frame$weight>mean(chg$frame$weight,na.rm=T))
  netout<-lapply(chg$nets, function(X) igraph::as_edgelist(X) %>% unique(.) %>% data.frame()) %>% dplyr::bind_rows(.id="Doc") %>% unique()
  simpidf<-function(col1,cutprob){
    temp<-stringr::str_split(col1,pattern=",") 
    temp<-sapply(temp, function(kt1){
      stringr::str_trim(kt1)[which.max(chg$frame$Weight[sapply(stringr::str_trim(kt1),function(K) which(K==chg$frame$words_out)[1])])]
    })
    temp[which(sapply(temp, length)==0)]<-NA
    temp<-do.call(c,temp)
    temp}
  netout$X1<-simpidf(netout$X1)
  netout$X2<-simpidf(netout$X2)
  netout$inW<-sapply(netout$X1,function(X) chg$frame$Weight[which(chg$frame$words_out==X)[1]])
  netout$outW<-sapply(netout$X2,function(X) chg$frame$Weight[which(chg$frame$words_out==X)[1]])
  netout2<-dplyr::filter(netout,inW>quantile(inW,probs=cutprob,na.rm=T),outW>quantile(outW,probs=cutprob,na.rm=T))
  net1<-igraph::graph_from_data_frame(data.frame("In"=netout2$X1,"Out"=netout2$X2),directed=T)
  igraph::E(net1)$weight <- 1
  netsimp<-igraph::simplify(net1,edge.attr.comb=list(weight="sum","ignore"))
  netsimp1<-igraph::set_vertex_attr(netsimp, "PageRank", value = igraph::page_rank(netsimp)$vector)
  netsimp1<-igraph::set_edge_attr(netsimp1,"Orig", value = levels(chg$frame$Orig))
  d3ob<-netsimp %>% networkD3::igraph_to_networkD3() 
  d3ob$nodes$PageRank<-igraph::page_rank(netsimp)$vector*12
  d3ob$nodes$Membership<-igraph::cluster_walktrap(netsimp) %>% igraph::membership() %>% as.character()
  d3ob$nodes$Orig<-levels(chg$frame$Orig)[1]
  list(netsimp1,d3ob)}


netalldocs<-function(workingfolder,prob,connllsub,docs,starton=1,skips=69){
if(dir.exists(file.path(workingfolder,"networks"))==FALSE){
  dir.create(file.path(workingfolder,"networks"))
}
for(DOC1 in docs$Orig[starton:length(docs$Orig)] %>% .[.!=docs$Orig[skips]]){
  outnet<-nets.to.ig.d3(conligraph(collect(dplyr::filter(connllsub, Orig==DOC1,prob))))
  saveRDS(outnet,file.path(workingfolder,"networks",paste("net",which(docs$Orig==DOC1),".rds",sep="")))
  cat(DOC1)
}
}


library(plyr)
library(dplyr)

setwd("vignettes")
#doc 69,101 didnt work
netalldocs(workingfolder,.5,spinsub,docs,starton=102,skips=101)
as.data.frame.igraph= function(g) {
  # prepare data frame
  res= cbind(as.data.frame(igraph::get.edgelist(g)),
             intergraph::asDF(g)$edges)[ , c(-3, -4)]
  # unfactorize
  res$V1= as.character(res$V1)
  res$V2= as.character(res$V2)
  # return df
  res
}

bignet<-lapply(list.files(file.path(workingfolder,"networks"),full.names=F),function(X){readRDS(file.path(workingfolder,"networks",X))[[1]]})
names(bignet)<-docs$Orig[-c(69,101)]
bignet2<-lapply(bignet,function(X){as.data.frame.igraph(X)})
names(bignet2)<-names(bignet)
#bignet2<-dplyr::bind_rows(bignet2,.id="Origindex") %>%  group_by(V1, V2) %>% summarise(weight= sum(weight))
bignet2<-dplyr::bind_rows(bignet2,.id="Origindex")
head(bignet2)
bignet2<-bignet2[,c(2,3,4,1)]
head(bignet2)

bignet2<-filter(bignet2,V1%in%tm::stopwords()==FALSE,V2%in%tm::stopwords()==FALSE,stringr::str_detect(V1,"[0-9]+")==FALSE,stringr::str_detect(V2,"[0-9]+")==FALSE,weight>5)

bignet<-igraph::graph.data.frame(bignet2, directed = T)
table(bignet2$Origindex)

bigcoms<-igraph::cluster_walktrap(bignet) %>% igraph::membership()

bignet<-networkD3::igraph_to_networkD3(bignet)



bignet[[1]]$Origindex<-bignet2$Origindex
bignet[[1]]$Weight<-bignet2$weight
bignet[[2]]$Community<-bigcoms

bignet[[1]]

networkD3::sankeyNetwork(bignet$links,bignet$nodes,"source","target","value",NodeID="name",NodeGroup="Community")

head(bignet$links)

networkD3::saveNetwork(bigsankey,"sankeyout.html",selfcontained=TRUE)


bigframe<-filter(bigframe,V1%in%tm::stopwords()==FALSE,V2%in%tm::stopwords()==FALSE,stringr::str_detect(V1,"[0-9]+")==FALSE,stringr::str_detect(V2,"[0-9]+")==FALSE)



bignet=igraph::simplify(igraph::graph.data.frame(bigframe, directed = T),remove.multiple=T,remove.loops=FALSE)
igraph::E(bignet)$weight=bigframe$weight
igraph::V(bignet)$degree=igraph::degree(bignet)
igraph::V(bignet)$pagerank=igraph::page_rank(bignet)$vector
biglabs<-igraph::cluster_walktrap(bignet)
memvec<-igraph::membership(biglabs)
igraph::V(bignet)$Walk.membership=memvec

sub1<-igraph::induced_subgraph(bignet,igraph::V(bignet)[igraph::V(bignet)$Walk.membership==3]) %>% list(.,networkD3::igraph_to_networkD3(.)) 

sub1[[2]]$nodes$PageRank<-igraph::page_rank(sub1[[1]])$vector
sub1[[2]]$nodes$Membership<-igraph::cluster_walktrap(sub1[[1]])%>% igraph::membership()

sankeysingle(sub1)


colnames(bd)[1]<-"degree"
bd$vertex<-row.names(bd)
bd<-bd[order(bd$degree,decreasing=T),]
bd<-dplyr::filter(bd, vertex%in%tm::stopwords()==FALSE)

bd[1:15,]




bignetd3<-lapply(list.files(file.path(workingfolder,"networks"),full.names=T),function(X) readRDS(X)[[2]])


  sankeysingle<-function(net.to.object){
    networkD3::sankeyNetwork(net.to.object$links,net.to.object$nodes,"source","target","value","name",NodeGroup="Membership",fontSize="PageRank")
  }
  
  sankeysingle(bignetd3[[90]])

  
  
  


networkD3::bignet

join(as.data.frame.igraph(bignet[[1]]),as.data.frame.igraph(bignet[[2]])
?join
union(bignet[[1]],bignet[[2]],byname=T)
?igraph::union

temp<-nets.to.ig.d3(chg,.5)
temp[[2]][[1]]
temp[[2]][[2]]

igraph::cluster_


saveRDS(do.call(igraph::union, lapply(list.files(file.path(workingfolder,network),full.names=T),function(X) readRDS(X)[[1]])),file.path(workingfolder,"bignetwork.rds"))

spinout$Terms<-tolower(spinout$tokens)

spinout<-left_join(spinout,dtmidf2)
colnames(dtmidf2)[3]<-"weight"
spinout<-plyr::join(spinout, dtmidf2,by="Terms",type="left",match="first")

conligraph<-function(parseconll){
  spl<-lapply(unique(parseconll$docname),function(X) subset(parseconll, docname==X))
  library(igraph)
  splout<-lapply(spl,function(splt) {
    g1a<-data.frame("from"=splt$tokens, "to"=splt$tokens[sapply(splt$head_id,function(K) which(splt$id==K))])
    g1<-g1a %>% igraph:: graph_from_data_frame(directed=TRUE) 
    commst<-igraph::walktrap.community(g1)
    tg<-igraph::contract(g1,mapping=igraph::membership(igraph::walktrap.community(g1)),vertex.attr.comb=toString)
    # tokenizers::tokenize_word_stems(stopwords=tokenizers::stopwords("en")) 
    V(tg)$name <-V(tg)$name %>% tolower() %>% lapply(.,function(X) paste(unique(X),collapse=","))
    commst<-commst %>% membership() %>% as.matrix()
    splt2<-left_join(splt,data.frame("tokens"=names(commst[,1]),"community"=commst[,1]))
    splt2$words_out<-splt2$tokens
    list(tg,splt2)})
  list("frame"=lapply(splout,function(X) X[[2]]) %>% bind_rows(),"nets"=lapply(splout,function(X) X[[1]]))
}

chg<-conligraph(filter(spinsub,Doc==)
                
nets.to.ig.d3.old<-function(chg){
#chg$frame<-subset(chg$frame, chg$frame$weight>mean(chg$frame$weight,na.rm=T))
netout<-lapply(chg$nets, function(X) igraph::as_edgelist(X) %>% unique(.) %>% data.frame()) %>% bind_rows(.id="Doc") %>% unique()
simpidf<-function(col1){
temp<-stringr::str_split(col1,pattern=",") 
temp<-sapply(temp, function(kt1){
str_trim(kt1)[which.max(chg$frame$weight[sapply(str_trim(kt1),function(K) which(K==chg$frame$Terms)[1])])]
})
temp[which(sapply(temp, length)==0)]<-NA
 temp<-do.call(c,temp)
 temp}
 netout$X1<-simpidf(netout$X1)
 netout$X2<-simpidf(netout$X2)
 netout$inW<-sapply(netout$X1,function(X) chg$frame$weight[which(chg$frame$words_out==X)[1]])
 netout$outW<-sapply(netout$X2,function(X) chg$frame$weight[which(chg$frame$words_out==X)[1]])
 netout2<-dplyr::filter(netout,inW>quantile(inW,,probs=.5,na.rm=T),outW>quantile(outW,probs=.5,na.rm=T))
 net1<-igraph::graph_from_data_frame(data.frame("In"=netout2$X1,"Out"=netout2$X2),directed=T)
 E(net1)$weight <- 1
 netsimp<-simplify(net1,edge.attr.comb=list(weight="sum","ignore"))
 netsimp1<-set_vertex_attr(netsimp, "PageRank", value = page_rank(netsimp)$vector)
 netsimp1<-set_edge_attr(netsimp1,"Orig", value = levels(chg$frame$Docs))
 d3ob<-netsimp %>% networkD3::igraph_to_networkD3() 
 d3ob$nodes$PageRank<-page_rank(netsimp)$vector*12
 d3ob$nodes$Membership<-cluster_edge_betweenness(netsimp) %>% igraph::membership() %>% as.character()
 list(netsimp1,d3ob)}


cutprob<-.8
nets.to.ig.d3<-function(chg,cutprob){
  #chg$frame<-subset(chg$frame, chg$frame$weight>mean(chg$frame$weight,na.rm=T))
  netout<-lapply(chg$nets, function(X) igraph::as_edgelist(X) %>% unique(.) %>% data.frame()) %>% dplyr::bind_rows(.id="Doc") %>% unique()
  simpidf<-function(col1,cutprob){
    temp<-stringr::str_split(col1,pattern=",") 
    temp<-sapply(temp, function(kt1){
      stringr::str_trim(kt1)[which.max(chg$frame$Weight[sapply(stringr::str_trim(kt1),function(K) which(K==chg$frame$words_out)[1])])]
    })
    temp[which(sapply(temp, length)==0)]<-NA
    temp<-do.call(c,temp)
    temp}
  netout$X1<-simpidf(netout$X1)
  netout$X2<-simpidf(netout$X2)
  netout$inW<-sapply(netout$X1,function(X) chg$frame$Weight[which(chg$frame$words_out==X)[1]])
  netout$outW<-sapply(netout$X2,function(X) chg$frame$Weight[which(chg$frame$words_out==X)[1]])
  netout2<-dplyr::filter(netout,inW>quantile(inW,probs=cutprob,na.rm=T),outW>quantile(outW,probs=cutprob,na.rm=T))
  net1<-igraph::graph_from_data_frame(data.frame("In"=netout2$X1,"Out"=netout2$X2),directed=T)
  igraph::E(net1)$weight <- 1
  netsimp<-igraph::simplify(net1,edge.attr.comb=list(weight="sum","ignore"))
  netsimp1<-igraph::set_vertex_attr(netsimp, "PageRank", value = page_rank(netsimp)$vector)
  netsimp1<-igraph::set_edge_attr(netsimp1,"Orig", value = levels(chg$frame$Orig))
  d3ob<-netsimp %>% networkD3::igraph_to_networkD3() 
  d3ob$nodes$PageRank<-igraph::page_rank(netsimp)$vector*12
  d3ob$nodes$Membership<-igraph::cluster_edge_betweenness(netsimp) %>% igraph::membership() %>% as.character()
  list(netsimp1,d3ob)}

temp<-nets.to.ig.d3(chg)
                


                
                
                
library(spacyr)

install_spacy<-function(){if("spacy"%in%system("pip list")==FALSE){system("cd; sudo pip install -U spacy; sudo python -m spacy.en.download")}}
install_spacy()
devtools::install_github("kbenoit/spacyr")

spacy_test<-function() {system("python -c 'import spacy; spacy.load('en'); print('OK')'")}

spin<-spacy_initialize()

dtmidf <- tm::DocumentTermMatrix(corpus1, control = list(weighting = weightTfIdf))

parselist<-function(pstrings,workingfolder){
  spcydb<-dplyr::src_sqlite(workingfolder,'spacyframe.db',create=T)
  l_ply(pstrings$Orig,function(X){
    ptemp<-dplyr::filter(pstrings,Orig==X)
    spinout<-spacy_parse(ptemp$string,dependency=T,named_entity=T,full_parse = TRUE)
    keep2s<-which(Terms(dtmidf)%in%tolower(spinout$tokens))
    dtmidf2<-as.matrix(dtmidf[X,keep2s])
    colnames(dtmidf2)<-Terms(dtmidf[X,keep2s])
    dtmidf2<-reshape2::melt(dtmidf2)
    head(dtmidf2)
    spinout$Terms<-tolower(spinout$tokens)
    head(spinout)
    spinout<-left_join(spinout,dtmidf2)
    colnames(dtmidf2)[3]<-"weight"
    spinout<-plyr::join(spinout, dtmidf2,by="Terms",type="left",match="first")
    dplyr::copy_to(spcydb,spinout,X)
  })}

loadwordnet<-function(){
  tempdir<-list.files(path="~",pattern="sqlite2.db",full.names=T)
  if(length(tempdir)==0){if(menu(c('yes','no'),title="Okay to download and install wordnet db 650mb?")==1){
    download.file('https://sourceforge.net/projects/sqlunet/files/4.0.0/sqlite/sqlite-4.0.0-31-all.zip/download',"~/wndbf.db")
    tempdir<-list.files(path="~",pattern="sqlite2.db",full.names=T)}}
  dplyr::src_sqlite(tempdir)
}

wn<-loadwordnet()

IDframes<-function(VERBWORD,parsecnnl){
  wn<-loadwordnet()
  filter(tbl(wn, "fnwords"), word%in%filter(parsecnnl,penn=="VB")$tokens) %>% left_join(tbl(wn,"fnlexemes")) %>% left_join(tbl(wn, "fnlexunits")) %>% left_join(tbl(wn, "fnframes")) %>% select(.,word,frame,framedefinition) %>% collect %>% filter(., stringr::str_detect(tolower(frame), VERBWORD))
}

idframes<-IDframes("[change|cause]",spinout)

sentkeeper<-function(idframes,parsecnnl){
  keepers<-dplyr::filter(parsecnnl, tolower(tokens)%in%idframes$word,penn=="VB")$docname
  dplyr::filter(parsecnnl, docname%in%keepers)
}

spinsub<-sentkeeper(idframes,spinout)

frameselect<-function(word_input){
  temp<-filter(wordtable,word==word_input)
  select.list(temp$frame,multiple=T)
}

conligraph<-function(parseconll){
  spl<-lapply(unique(parseconll$docname),function(X) subset(parseconll, docname==X))
  library(igraph)
  splout<-lapply(spl,function(splt) {
    g1a<-data.frame("from"=splt$tokens, "to"=splt$tokens[sapply(splt$head_id,function(K) which(splt$id==K))])
    g1<-g1a %>% igraph:: graph_from_data_frame(directed=TRUE) 
    commst<-igraph::walktrap.community(g1)
    tg<-igraph::contract(g1,mapping=igraph::membership(igraph::walktrap.community(g1)),vertex.attr.comb=toString)
    # tokenizers::tokenize_word_stems(stopwords=tokenizers::stopwords("en")) 
    V(tg)$name <-V(tg)$name %>% tolower() %>% lapply(.,function(X) paste(unique(X),collapse=","))
    commst<-commst %>% membership() %>% as.matrix()
    splt2<-left_join(splt,data.frame("tokens"=names(commst[,1]),"community"=commst[,1]))
    splt2$words_out<-splt2$tokens
    list(tg,splt2)})
  list("frame"=lapply(splout,function(X) X[[2]]) %>% bind_rows(),"nets"=lapply(splout,function(X) X[[1]]))
}

chg<-conligraph(spinsub)

library(igraph)
?str_trim
nets.to.ig.d3<-function(chg){
#chg$frame<-subset(chg$frame, chg$frame$weight>mean(chg$frame$weight,na.rm=T))
netout<-lapply(chg$nets, function(X) igraph::as_edgelist(X) %>% unique(.) %>% data.frame()) %>% dplyr::bind_rows(.id="Doc") %>% unique()
simpidf<-function(col1){
  temp<-stringr::str_split(col1,pattern=",") 
  temp<-sapply(temp, function(kt1){
    stringr::str_trim(kt1)[which.max(chg$frame$weight[sapply(stringr::str_trim(kt1),function(K) which(K==chg$frame$Terms)[1])])]
  })
  temp[which(sapply(temp, length)==0)]<-NA
  temp<-do.call(c,temp)
  temp}
netout$X1<-simpidf(netout$X1)
netout$X2<-simpidf(netout$X2)
netout$inW<-sapply(netout$X1,function(X) chg$frame$weight[which(chg$frame$words_out==X)[1]])
netout$outW<-sapply(netout$X2,function(X) chg$frame$weight[which(chg$frame$words_out==X)[1]])
netout2<-dplyr::filter(netout,inW>quantile(inW,probs=.5,na.rm=T),outW>quantile(outW,probs=.5,na.rm=T))
net1<-igraph::graph_from_data_frame(data.frame("In"=netout2$X1,"Out"=netout2$X2),directed=T)
E(net1)$weight <- 1
netsimp<-simplify(net1,edge.attr.comb=list(weight="sum","ignore"))
netsimp1<-set_vertex_attr(netsimp, "PageRank", value = page_rank(netsimp)$vector)
netsimp1<-set_edge_attr(netsimp1,"Orig", value = levels(chg$frame$Docs))
d3ob<-netsimp %>% networkD3::igraph_to_networkD3() 
d3ob$nodes$PageRank<-page_rank(netsimp)$vector*12
d3ob$nodes$Membership<-cluster_edge_betweenness(netsimp) %>% igraph::membership() %>% as.character()
list(netsimp1,d3ob)}

temp<-nets.to.ig.d3(chg)

sankeysingle<-function(net.to.object){
sankeyNetwork(net.to.object$links,net.to.object$nodes,"source","target","value","name",NodeGroup="Membership",fontSize="PageRank")
}


sankeysingle(temp[[2]])



parselist<-function(pstrings,workingfolder1,dbexists1){
  spcydb<-if(dbexists1==F){dplyr::src_sqlite(file.path(workingfolder,'spacyframe.db'),create=T)} else {dplyr::src_sqlite(file.path(workingfolder1,'spacyframe.db'),create=F)}
  needlist<-if(dbexists1==T){unique(pstrings$Orig)[unique(pstrings$Orig)%in%src_tbls(spcydb)==FALSE]} else {unique(pstrings$Orig)}
  l_ply(needlist,function(X){
    ptemp<-dplyr::filter(pstrings,Orig==X)
    spin<-spacyr::spacy_initialize()
    spinout<-spacy_parse(ptemp$string,dependency=T,named_entity=T,full_parse = TRUE)
    spacyr::spacy_finalize()
    keep2s<-which(Terms(dtmidf)%in%tolower(spinout$tokens))
    dtmidf2<-as.matrix(dtmidf[X,keep2s])
    colnames(dtmidf2)<-Terms(dtmidf[X,keep2s])
    dtmidf2<-reshape2::melt(dtmidf2)
    spinout$Terms<-tolower(spinout$tokens)
    spinout<-left_join(spinout,dtmidf2)
    colnames(dtmidf2)[3]<-"weight"
    spinout<-plyr::join(spinout, dtmidf2,by="Terms",type="left",match="first")
    dplyr::copy_to(spcydb,spinout,X)
  })}
parselist(pstrings,workingfolder,dbexists)
