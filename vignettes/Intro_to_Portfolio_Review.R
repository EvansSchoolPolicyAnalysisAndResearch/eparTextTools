## ----eval=FALSE----------------------------------------------------------
#  devtools::install_github("ryscott5/epartexttools")

## ------------------------------------------------------------------------
options(java.parameters = "-Xmx4g")
library(epartexttools)

## ---- echo=TRUE, results='hide', eval=F----------------------------------
#  example_documents()
#  corpus1<-allDocs(directory="demo.docs.folder",SkiponError=FALSE)

## ------------------------------------------------------------------------
#Make directory for saving and save a copy of the corpus
epartexttools::makeworking("Research.Grants")

## ---- eval=F-------------------------------------------------------------
#  saveRDS(corpus1,file.path(workingfolder, "corpus.rds"))

## ---- echo=T, eval=T-----------------------------------------------------
corpus1<-readRDS(file.path(workingfolder,"corpus.rds"))

## ---- echo=TRUE----------------------------------------------------------
lapply(corpus1,function(X){X$meta})[[1]]

## ---- echo=TRUE, results='asis', eval=FALSE------------------------------
#  #cleans elements of corpus
#  corpus2<-doc_clean_process(corpus1)
#  saveRDS(corpus2,file.path(workingfolder,"corpus_cleaned.rds"))

## ----echo=FALSE----------------------------------------------------------
corpus2<-readRDS(file.path(workingfolder,"corpus_cleaned.rds"))

## ------------------------------------------------------------------------
#converts corpus into term document matrix, removing terms that occur infrequently (can be adjusted by manipulating that .6)
tdm<-TermDocumentMatrix(corpus2) %>% removeSparseTerms(.,.8)

## ----echo=TRUE, results='asis'-------------------------------------------
#Searches corpus 1 for the words gender and access based and returns counts based on the term document matrix we built above.
wordcount_table(c("gender","access"),tdm,corpus1)

## ----echo=TRUE, results='asis'-------------------------------------------
tout<-wordcount_table(c("gender","access"),tm::TermDocumentMatrix(corpus2,control=list(weighting=function(X) tm::weightTfIdf(X, normalize=FALSE))),corpus1,raw=T)
head(tout[,1:3])

## ----echo=TRUE, results='asis'-------------------------------------------
tdm[,as.vector(tdm["gender",])>1] %>% .[,as.vector(.["women",])>1] %>% word_heatmap(.,20)

## ------------------------------------------------------------------------
wfplots(tdm[,c(2,10,12)],typePlot=2, 5,shortendoc=TRUE)

## ------------------------------------------------------------------------
tdm[,as.vector(tdm["gender",])>20] %>% .[,as.vector(.["women",])>10] %>% wfplots(.,typePlot=2,10,shortendoc=TRUE)

## ------------------------------------------------------------------------
tdm[,as.vector(tdm["gender",])>20] %>% .[,as.vector(.["women",])>10] %>% wfplots(.,typePlot=1,10,shortendoc=TRUE)+ggtitle("Change X versus Y")

## ------------------------------------------------------------------------
interest_plot_bydoc(c("women","farmer","school"),tdm[,1:5])+coord_flip() 

## ------------------------------------------------------------------------
interest_plot_bydoc(c("women","farmer","school"),tdm[,1:5]) %>% plotly::ggplotly() 

## ------------------------------------------------------------------------
TermDocumentMatrix(corpus2[1:10],control=list(weighting=function(x) weightSMART(x))) %>% interest_plot_bydoc(c("women","farmer","school"),.) %>% plotly::ggplotly() 

## ----echo=TRUE, results='asis'-------------------------------------------
word_heatmap(tdm,6)

word_heatmap(tdm,pickwords=c("women","gender","access","land","right","work","labor","yield","security"))

## ------------------------------------------------------------------------
assocPrettyOneStep("gender",tdm, corpus2,.5)

## ------------------------------------------------------------------------
tornadoCompare(tdm,c("gender","equal","femal"),3,10)

## ----eval=FALSE----------------------------------------------------------
#  BASE_INPUT<-PreTopicFrame2(corpus1,workingfolder=workingfolder,removeentities=F)
#  BASE_INPUT$out$meta$OpID<-BASE_INPUT$out$meta$Orig
#  #saves files so you can reload
#  saveRDS(BASE_INPUT,file.path(workingfolder,"base_input1.rds"))

## ---- eval=FALSE---------------------------------------------------------
#  buildcliff()
#  startcliff()
#  library(RCurl)
#  library(httr)
#  BASE_INPUT$SentFrame$OpID<-BASE_INPUT$SentFrame$Orig
#  pred1<-PredictCountryByDoc(BASE_INPUT)
#  stopcliff()
#  BASE_INPUT$out$meta<-reflectCountryCol(BASE_INPUT$out$meta,pred1,10,FALSE)
#  getwd()
#  saveRDS(BASE_INPUT,file.path(workingfolder,"base_input1.rds"))
#  write.csv(pred1,file.path(workingfolder,"countrypredictions1.csv"))

## ----eval=FALSE----------------------------------------------------------
#  library(plotly)
#  runMap(file.path(workingfolder,"countrypredictions1.csv"),path.file=T,"countries")

## ----eval=FALSE----------------------------------------------------------
#  writeFormulaforSTM(BASE_INPUT,workingfolder)

## ----eval=FALSE----------------------------------------------------------
#  runSTM(workingfolder)

