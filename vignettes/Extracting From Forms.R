## ------------------------------------------------------------------------
tfiles<-list.files("~/CRO_form_demo",full.names=T)
tfiles<-tfiles[str_detect(tolower(tfiles),"docx")]

## ------------------------------------------------------------------------
head(docx_table_view(tfiles[1], export_frame=T,showView=F))

## ------------------------------------------------------------------------
cell_extractor(docx_table_view(tfiles[1],export_frame=T,showView=F)$content,"CO Resiliency","SUMMARY")

## ------------------------------------------------------------------------
sapply(tfiles, function(X) tryCatch({cell_extractor(docx_table_view(X,export_frame=T,showView=F)$content,"CO Resiliency","SUMMARY")},error=function(e) {NA}))

## ------------------------------------------------------------------------
formcluster(tfiles)

## ---- eval=F-------------------------------------------------------------
#  docx_table_view(tfiles[8])
#  docx_table_view(tfiles[1])

