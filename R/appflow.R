#'Write a shiny App
#'
#' Write a shiny app that draws network connections between words.
#' @param workingfolder the working folder
#' @param app.name name of the application (will be a folder name)
#' @param run.app should the shiny app be run?
#' @return a shiny app in the app.name folder within the working folder
#' @example
#' writeShinyFlow(workingfolder,"testapp")
writeShinyFlow<-function(workingfolder,app.name,run.app=TRUE){
dir.create(file.path(workingfolder,app.name))
writeLines(text=c('library(data.table)
library(devtools)
library(shiny)
library(igraph)
library(shinyapps)
library(DT)
library(stringr)
library(networkD3)
igraphob_object_force<-function(WORD,mtable,W,inputWord=TRUE,sankey=FALSE,verbfilter=c(),strictlimit=FALSE){
  if(inputWord==TRUE){
    ntr<-dplyr::filter(mtable,str_detect(tolower(subject.keywords),tolower(WORD)))
  } else {
    ntr<-dplyr::filter(mtable,str_detect(tolower(object.keywords),tolower(WORD)))
  }
  dt <- data.table(ntr)
  dt$object.keywords<-sapply(dt$object.keywords,function(X) paste(unlist(X),collapse=";"))
  dt$subject.keywords<-sapply(dt$subject.keywords,function(X) paste(unlist(X),collapse=";"))
  dt$action.lemmatized<-sapply(dt$action.lemmatized,function(X) paste(unlist(X),collapse=";"))
  dt$action.verb.text<-sapply(dt$action.verb.text,function(X) paste(unlist(X),collapse=";"))
  suppressWarnings(dt<-dt[,list(action.verb.text = unlist(strsplit(action.verb.text,";")),object.keywords=unlist(strsplit(object.keywords,";")),subject.keywords=unlist(strsplit(subject.keywords,";"))),by = sentence])
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
    forceNetwork(d3ob$links,d3ob$nodes,"source","target","value",NodeID="name",Group="group",fontSize=24,fontFamily="Arial", opacity = 0.9, bounded = TRUE, opacityNoHover = TRUE,colourScale=JS("d3.scale.ordinal().domain(\'verb\',\'subject\',\'object\']).range([\'#d9d9d9\',\'#b7a57a\',\'#4b2e83\'])"))} else {
      sankeyNetwork(d3ob$links,d3ob$nodes,"source","target","value",NodeID="name",NodeGroup="group",colourScale=JS("d3.scale.ordinal().domain([\'verb\',\'subject\',\'object\']).range([\'#d9d9d9\',\'#b7a57a\',\'#4b2e83\'])"),fontSize=24,fontFamily="Arial")}}
matchtable<-readRDS("matchtable.rds")',
'ui <- shinyUI(fluidPage(
   sidebarLayout(
      sidebarPanel(
        selectizeInput("topics","Topics",unique(matchtable$TopTopics),multiple=TRUE),
         textInput("word","Word:","Potato"),
         sliderInput("fstrength","Filter Strength",0,20,0,step=.5),
         checkboxInput("inword","Subject",TRUE),
         selectizeInput("verbs","Links",unique(matchtable$action.verb.text),multiple=TRUE),
         checkboxInput("strict","Strict Word?",TRUE),
        downloadButton("downloadData", "Download")),
      mainPanel(
         sankeyNetworkOutput("sankplot"),
         DT::dataTableOutput("senttab")))))
server <- shinyServer(function(input, output) {
output$sankplot <- renderSankeyNetwork({
     req(input$word)
     mframe<-if(length(input$topics)>=1){dplyr::filter(matchtable,TopTopics%in%input$topics)} else {matchtable}
     igraphob_object_force(input$word,mframe,input$fstrength,inputWord=input$inword,sankey=TRUE,verbfilter=input$verbs,strictlimit=input$strict)})
   
output$senttab<-renderDataTable({
     req(input$word)
      mframe<-if(length(input$topics)>=1){dplyr::filter(matchtable,TopTopics%in%input$topics)} else {matchtable}
  
     if(input$inword==TRUE){
       ntr<-dplyr::filter(mframe,str_detect(tolower(subject.keywords),tolower(input$word)))
     }
     else {
       ntr<-dplyr::filter(mframe,str_detect(tolower(object.keywords),tolower(input$word)))
     }
     dt <- data.table(ntr)
     dt$object.keywords<-sapply(dt$object.keywords,function(X) paste(unlist(X),collapse=";"))
     dt$subject.keywords<-sapply(dt$subject.keywords,function(X) paste(unlist(X),collapse=";"))
     dt$action.lemmatized<-sapply(dt$action.lemmatized,function(X) paste(unlist(X),collapse=";"))
     dt$action.verb.text<-sapply(dt$action.verb.text,function(X) paste(unlist(X),collapse=";"))
     suppressWarnings(dt<-dt[,list(action.verb.text = unlist(strsplit(action.verb.text,";")),object.keywords=unlist(strsplit(object.keywords,";")),subject.keywords=unlist(strsplit(subject.keywords,";"))),by = sentence])
     dt<-na.omit(dt)
     if(length(input$verbs)>0){
       dt<-dplyr::filter(dt,action.verb.text%in%input$verbs)
     }
     if(input$strict==TRUE){
       if(input$inword==TRUE){
         dt<-dplyr::filter(dt,str_detect(tolower(subject.keywords),tolower(input$word)))} else {dt<-dplyr::filter(dt,str_detect(tolower(object.keywords),tolower(input$word)))}
     }
     datatable(unique(data.frame(dt)[,c("subject.keywords","sentence")]))
   })
output$downloadData <- downloadHandler(
  filename = function() {
    paste(input$word, "html", sep = ".")
  },
  content = function(file) {
    req(input$word)
    mframe<-if(length(input$topics)>=1){dplyr::filter(matchtable,TopTopics%in%input$topics)} else {matchtable}
    saveNetwork(igraphob_object_force(input$word,mframe,input$fstrength,inputWord=input$inword,sankey=TRUE,verbfilter=input$verbs,strictlimit=input$strict),file=file)
  })})','shinyApp(ui = ui, server = server)'),con=file.path(workingfolder,app.name,'app.R'))
file.copy(file.path(workingfolder,"matchtable.rds"),file.path(workingfolder,app.name,"matchtable.rds"))
if(run.app==TRUE){shiny::runApp(file.path(workingfolder,app.name))}}

