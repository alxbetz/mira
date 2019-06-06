#dependencies: Rbowtie
# x: named vector containing parameters
# cDNAFile
# platformFile
# probeIdColName
# probeSeqColName
# x = c("cDNAFile" = "/path/to/cDNAFile","platformFile" = "/path/to/platformFile","probeIdColName" = "ID","probeIdSeqName"="SEQUENCE")
#wdir: where stuff will be written
#pipedir: where this script is located
#usage: mapProbes(x,wdir,pipedir)

#' Write probe sequences into a fasta file for celera genomics microarray files
#'
#' @param seqs
#' @param outFile
#'
#' @return
#' @export
#'
#' @examples
writeTagFasta = function(seqs,outFile) {
  out = file(outFile)
  nseq = length(seqs)
  print(nseq)
  idList = paste(">",as.character(seq(nseq)),sep="")
  seqList = sapply(seqs,as.character)
  print(head(seqList))
  toWrite = c(rbind(idList,seqList))
  writeLines(toWrite,con=out)
  close(out)
}




#' Build the bowtie index
#'
#' @param cDNAFile file with cDNA sequences in fasta format
#' @param wdir working directory
#' @param pipedir pipeline directory
#' @param iDir directory for writing the bowtie index
#' @param prefix prefix string for index file
#'
#' @return path to bowtie index
#' @export
#'
#' @examples
build_bt_index = function(cDNAFile,wdir,pipedir,iDir,prefix="index",verbose=FALSE){


  #  species = strsplit(basename(cDNAFile),".")[[1]][1]
  #  iDir=paste(".",species,sep="/")
  idxLocation = file.path(iDir,prefix)
  if(verbose) {
    print(idxLocation)
  }

  if(length(list.files(idxLocation,pattern="*.ebwt") > 0)) {
    return(idxLocation)
  }
  #Bowtie does not work with relative paths
  #make sure reference file location is an absolute path
  if(substr(cDNAFile,1,2) == "./") {
    references = file.path(pipedir,substr(cDNAFile,3,nchar(cDNAFile)))
  } else {
    references = cDNAFile
  }
  print(references)
  Rbowtie::bowtie_build(references, iDir, prefix = prefix, force = T, strict = T)

  return(idxLocation)
}


#' Write Probe sequences into fasta file for all microarray definition files using tabular text format
#'
#' @param x
#' @param wdir
#'
#' @return
#' @export
#'
#' @examples
writeProbesAsFasta = function(x,wdir) {


  pInfo = readLines(x[["platformFile"]])
  #nSkip = grep("!platform_table_begin",pInfo)
  #for agilent array file format and ArrayExpress array file format
  nSkip = grep("!platform_table_begin|\\[main\\]",pInfo)
  if (length(nSkip) == 0){
    nSkip=0
  }
  probeTab = read.delim(x[["platformFile"]],skip=nSkip,stringsAsFactors = F)
  if(length(nSkip)>0) {
    probeTab = probeTab[1:dim(probeTab)[1]-1,]
  }

  print(colnames(probeTab))
  idCol = x[["probeIdColName"]]
  seqCol = x[["probeSeqColName"]]
  idCol = gsub(" ",".",idCol)
  seqCol =gsub(" ",".",seqCol)
  print(idCol)
  print(seqCol)
  info = probeTab[,c(idCol,seqCol)]
  info[,idCol] = paste(">",info[,idCol],sep="")
  info =c(rbind(info[,1],info[,2]))
  #print("test")
  #print(typeof(x[["eid"]]))
  qFile=file.path(wdir,"probeSet.fasta")
  print(head(info))
  #write(paste(info,collapse="\n"),file=file.path(wdir,paste(x[["eid"]],".fasta",sep="")))
  write(paste(info,collapse="\n"),file=qFile)
  return(qFile)

}

#' map probe sequences to reference genome
#'
#' @param x named vector containing parameters
#' cDNAFile
#' platformFile
#' probeIdColName
#' probeSeqColName
#' @param wdir STRING working directory
#' @param pipedir STRING pipeline directory
#' @param prefix STRING
#' @param nCores INT number of cores bowtie should use
#'
#' @return
#' @export
#'
#' @examples
mapProbes = function(x,wdir,pipedir,prefix='index',nCores=4) {
  iDir =file.path(wdir,"btIndex")
  cDNAFile = x[['cDNAFile']]
  #function(cDNAFileDummy,wdir,pipedir,iDir,prefix="index")
  iFile = build_bt_index(cDNAFile,wdir,pipedir,iDir,prefix=prefix)
  if(!file.exists(x[['platformFile']])) {
    x[['platformFile']] = downloadPlatformFile(x[['eid']],wdir)
  }
  qFile = writeProbesAsFasta(x,wdir)
  alFile = file.path(wdir,"alignment")
  Rbowtie::bowtie(qFile,iFile,type="single",outfile=alFile,f=TRUE,force=TRUE,a=TRUE,best=TRUE,strata=TRUE,m=10,v=3,p=nCores)
  return(alFile)
}
