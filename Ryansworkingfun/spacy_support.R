#spacy_support



install_spacy<-function() {if("spacy"%in%system("pip list")==FALSE){system("cd; sudo pip install -U spacy; sudo python -m spacy.en.download")}}

install_spacy()

spacy_test<-function() {system("python -c 'import spacy; spacy.load('en'); print('OK')'")}

spacy_test()
install.packages("withr")
devtools::install_github("kbenoit/spacyr")

library(spacyr)
spin<-spacy_initialize()


dtmidf <- tm::DocumentTermMatrix(corpus1, control = list(weighting = weightTfIdf))

parselist<-function(pstrings){
lapply(pstrings$Orig,function(X){
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
spinout
})
}


dtmidf[which(dtmidf$dimnames$Docs==X),spinout$tokens[1]]
keeps<-which(tolower(spinout$tokens)%in%dtmidf$dimnames$Terms==FALSE)
keep2s<-which(dtmidf$dimnames$Terms%in%tolower(spinout$tokens))


tokweight<-matrix(dtmidf)[1,tolower(spinout$tokens[-keeps])]

tv<-sapply(spinout$tokens, function(K) try(matrix(dtmidf[1,tolower(K)])[1],silent=T))
?try
?tm::inspect
hist(tweights)



length(dtmidf$dimnames$Terms)


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

spinout
sentkeeper<-function(idframes,parsecnnl){
  keepers<-dplyr::filter(parsecnnl, tolower(tokens)%in%idframes$word,penn=="VB")$docname
  dplyr::filter(parsecnnl, docname%in%keepers)
}

spinsub<-sentkeeper(idframes,spinout)

review_dtm_tfidf <- tm::DocumentTermMatrix(review_corpus, control = list(weighting = weightTfIdf))



#how to create wordtable data
#wn<-dplyr::src_sqlite("~/sqlite2.db")
#wordtable<-tbl(wn, "fnwords") %>% left_join(tbl(wn,"fnlexemes")) %>% left_join(tbl(wn, "fnlexunits")) %>% left_join(tbl(wn, "fnframes")) %>% select(.,word,frame,wordid,noccurs,lexunit,ludefinition,framedefinition) %>% collect(.) 

wn<-dplyr::src_sqlite("~/sqlite-31.db")
wordtable<-tbl(wn, "fnwords") %>% left_join(tbl(wn,"fnlexemes")) %>% left_join(tbl(wn, "fnlexunits")) %>% left_join(tbl(wn, "fnframes")) %>% select(.,word,frame,wordid,noccurs,lexunit,ludefinition,framedefinition) %>% collect(.) 
saveRDS(wordtable,"data/wordtable.RData")
save(wordtable,file="data/wordtable.RData")

frameselect<-function(word_input){
  temp<-filter(wordtable,word==word_input)
  select.list(temp$frame,multiple=T)
}



  
spinsub<-filter(spinout, penn=="VB",tolower(tokens)%in%na.omit(wordtable[str_detect(tolower(wordtable$frame),"produce"),"word"]))

wordtable

suggestkeywords<-function(keyword,path_root=F){
  cnetpath<-if(path_root==T){"http://api.conceptnet.io"} else {"http://api.conceptnet.io/c/en"}
  temp<-jsonlite::fromJSON(file.path(cnetpath,keyword))
  temp<-lapply(temp$edges$start$term, function(keyword1) {
    temp2<-jsonlite::fromJSON(file.path("http://api.conceptnet.io",keyword1))
    data.frame("label"=temp2$edges$start$label,"term"=temp2$edges$start$term)})
  list("word"=c(c(sapply(temp, function(X) as.character(X$term)) %>% unlist(.) %>% unique()) %>% unique),"term"=c(c(sapply(temp, function(X) as.character(X$label)) %>% unlist(.) %>% unique())%>% unique))
}

trybase<-suggestkeywords("coffee")
suggestkeywords(trybase$word[3],path_root=T)

  




"fnlexunits"               

vframes
src_tbls(wn)

download.file('https://sourceforge.net/projects/sqlunet/files/4.0.0/sqlite/sqlite-4.0.0-31-all.zip/download',"~/wndbf.db")
spinout
plyr::
tgf<-spinout
https://github.com/Noahs-ARK/semafor
http://www.nltk.org/howto/framenet.html

plyr::l_ply(complexregex,function(RGX){
tempp<-phrasemachine::extract_ngram_filter(spinout$penn,regex=RGX,maximum_ngram_length = 5,minimum_ngram_length = 2)
spinout$phrase<-NA
if(length(tempp)>0){
for(i in 1:nrow(tempp)){
  tgf$phrase[tempp[i,1]:tempp[i,2]]<-tempp[i,1]
}}})


parseconll<-spinsub

#Function to make community structure

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

chg$frame<-subset(chg$frame, chg$frame$weight>mean(chg$frame$weight,na.rm=T))

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

netout2<-filter(netout,inW>quantile(inW,,probs=.5,na.rm=T),outW>quantile(outW,probs=.5,na.rm=T))
?quantile
head(netout2)
net1<-igraph::graph_from_data_frame(data.frame("In"=netout2$X1,"Out"=netout2$X2),directed=T)
E(net1)$weight <- 1
netsimp<-simplify(net1,edge.attr.comb=list(weight="sum","ignore"))
PR<-page_rank(netsimp)$vector
V(netsimp)$VPR<-PR/sd(PR)
d3ob<-netsimp %>% networkD3::igraph_to_networkD3() 
d3ob$nodes$PageRank<-page_rank(netsimp)$vector
d3ob$nodes$PageRank*12
d3ob$nodes$Membership<-cluster_edge_betweenness(netsimp) %>% membership %>% as.character()

sankeyNetwork(d3ob$links,d3ob$nodes,"source","target","value","name",NodeGroup="Membership",fontSize="PageRank",fontColor="Membership")


?sankeyNetwork

plot(netsimp)
reshape2::melt(degree(netsimp))
netsimp
plot.igraph(netsimp,vertex.size=PageRank)
reshape2::melt(betweenness(netsimp))
page_rank(netsimp)$vector %>% t() %>% reshape2::melt() %>% .[order(-.$value),] 
?page_rank

clap<-
igraph::cluster_label_prop(netsimp) %>% communities()
?cluster_label_prop

plot(netsimp)

d3ob<-netsimp %>% networkD3::igraph_to_networkD3() 
d3ob1<-as.data.frame(igraph::get.edgelist(netsimp))

tbl(wn,"words") %>% filter(wordid==1003) %>% select(.,wordid)
d3ob1$V1[100]

tbl(wn,"wordsXsensesXsynsets") %>% filter(lemma=="collective") %>% collect(.) %>% .[which.max(.$tagcount),"definition"]

tbl(wn,"sumoterms")
tbl(wn,"sumorelations")
src_tbls(wn)
s(d3ob1)
?sankeyNetwork
sankeyNetwork(d3ob$links,d3ob$nodes,"source","target","value","name",NodeGroup=1)

length(unique(c(d3ob$links$source,d3ob$links$target)))
simpleNetwork(d3ob1)
length(d3ob$links$target)
              ?simpleNetwork
              d3ob$links$value))


networkD3::forceNetwork(d3ob)
d3ob$nodes$name
d3ob$links

networkD3::forceNetwork(d3ob$links,d3ob$nodes,"source","target","value",NodeID="name",fontSize=24,fontFamily="Arial", opacity = 0.9, bounded = TRUE, opacityNoHover = TRUE,colourScale=JS("d3.scale.ordinal().domain(['verb','subject','object']).range(['#d9d9d9','#b7a57a','#4b2e83'])"))

#graph_from_data_frame(rbind(data.frame("In"=dt$subject.keywords,"Out"=dt$action.verb.text),data.frame("In"=dt$action.verb.text,"Out"=dt$object.keywords)))
E(net1)$weight <- 1

netsimp<-simplify(net1,edge.attr.comb=list(weight="sum","ignore"))

E(netsimp)$width <- E(netsimp)$weight

netsimp <- delete_edges(netsimp, E(netsimp)[weight<=W])
netsimp<-simplify(netsimp,edge.attr.comb=list(weight="sum","ignore"))
bad.vs<-V(netsimp)[degree(netsimp) == 0]
netsimp <-delete.vertices(netsimp, bad.vs)
d3ob<-netsimp %>% igraph_to_networkD3() 
plot()
library(networkD3)
networkD3::forceNetwork(d3ob$links,d3ob$nodes,"source","target","value",NodeID="name",Group="group",fontSize=24,fontFamily="Arial", opacity = 0.9, bounded = TRUE, opacityNoHover = TRUE,colourScale=JS("d3.scale.ordinal().domain(['verb','subject','object']).range(['#d9d9d9','#b7a57a','#4b2e83'])"))








chg$frame$Unigroup<-interaction(chg$frame$docname,chg$frame$community)

install.packages("RKEA")



temp1<-reshape2::melt(xtabs(~tolower(chg$frame$words_out)+chg$frame$docname))
colnames(temp1)<-c('word','doc','count')
ddply(temp1, .(doc),)
table(temp1$word)

netout<-lapply(chg$nets, function(X) igraph::as_edgelist(X) %>% unique(.) %>% data.frame()) %>% bind_rows(.id="Doc") %>% unique()

head(netout)

tfidf<-tm::DocumentTermMatrix(corpus1, control = list(weighting = weightTfIdf))



ddply(chg$frame,.(docname, community), summarize, words=paste(words_out,collapse=" "))


chg$frame$wo
V(splout[[5]])$name<-V(splout[[5]])$name %>% tolower() %>% tokenizers::tokenize_word_stems(stopwords=tokenizers::stopwords("en")) %>% lapply(.,unique)
}



vertex_attr(g1,"commune")<-membership(commst)

mgraph<-lapply(1:length(commst),function(i) igraph::subgraph(g1,vertex_attr(g1,"commune")==i))
page_rank(g1)
igraph::page_rank(mgraph[[1]])
library(igraph)
V(g1)[g1$commune]
igraph::contract(g1,mapping=igraph::membership(igraph::walktrap.community(g1)),vertex.attr.comb=toString) %>% plot()

igraph::contract(g1,mapping=igraph::membership(igraph::cluster_louvain(g1)),vertex.attr.comb=toString) %>% plot()

g1a$membership<-igraph::walktrap.community(g1) %>% igraph::membership()
spinout$membership<-g1$membership


g1$community<-igraph::membership()
?igraph::subgraph
igraph::subgraph(g1,"community"==1)
?spacy_parse
install.packages("phrasemachine")
phrasemachine::extract_phrases(spinout)
?extract_phrases
tgf<-spinout
spacy_finalize()
restructure_spacy<-function(tgf, tag="penn"){
  lapply(unique(tgf$docname),function(X){
  tgf1<-dplyr::filter(tgf,docname==X)
  POSout<-list("token"=tgf1[,3],"tags"=tgf1[,tag])
  POSout
  })
}


