#dependencies: oligo

processCeleraFiles = function(fileList,cDNAFile,wdir,pipedir,sampleNames=NULL,nCores=4,nMethod="quantile",bgMethod="rma",pipeComponents="ande") {
  #browser()
  rawData = oligo::read.celfiles(fileList,sampleNames=sampleNames)
  qFile = file.path(wdir,"probeSet.fasta")
  print(head(pmSequence(rawData)))

  writeTagFasta(pmSequence(rawData),qFile)

  cDNAFileRaw = unlist(cDNAFile)
  iDir =file.path(wdir,"btIndex")
  idxLocation = file.path(iDir,basename(cDNAFileRaw))
  #if(grepl("a",pipeComponents)){
  iFile= build_bt_index(cDNAFile,wdir,pipedir,iDir,prefix=basename(cDNAFileRaw))
  #}

  alFile = file.path(wdir,"alignment")
  bgData = oligo::backgroundCorrect(rawData,method=bgMethod)
  normData = oligo::normalize(bgData,method=nMethod)
  #normData = oligo::rma(rawData)
  print(head(normData))
  print(typeof(normData))
  print(class(normData))
  normEset = assayData(normData)$exprs
  #transform to logspace
  normEset = log2(normEset)

  print("rownames")
  print(head(rownames(normEset)))
  print(head(normEset))
  Rbowtie::bowtie(qFile,iFile,type="single",outfile=alFile,f=TRUE,force=TRUE,a=TRUE,best=TRUE,strata=TRUE,m=10,v=3,p=nCores)
  return(list(normEset=normEset,alFile=alFile,normData=normData))
}

unifyAnnotation = function(annotation) {
  acols = colnames(annotation)
  acols = gsub("group","Group",acols)
  acols = gsub("filename","File",acols)
  acols = gsub("ID","Name",acols)
  colnames(annotation) = acols
  return(annotation)
}


readAnnotation = function(pipedir,x) {
  annoFile = file.path(pipedir,"input",x['eid'],"annotation.txt")
  annotation = read.delim(annoFile,sep="\t",stringsAsFactors = F)
  annotation = unifyAnnotation(annotation)
}



annotationDualToSingleChannel = function(annotation,controlColor) {
  annoBG = annotation %>% mutate(Group = paste0(Group,"_control"))

  if(controlColor=="green") {
    annotation =  rbind(annotation,annoBG)
  } else {
    annotation = rbind(annoBG,annotation)
  }
  return(annotation)

}
