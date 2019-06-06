#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
getMapping = function(x) {
  spl = strsplit(x,split="\t")
  as.list(spl[[1]][c(3,1)])
}

#' Summarise intensity values for one gene locus
#'
#' @param eid "experiment identifier"
#' @param alFile "alignment file"
#' @param normEset "normalised Expressionset"
#' @param wdir "working directory"
#' @param sumFUN "summarsation function, usually median or mean"
#'
#' @return
#' @export "table of normalised expression values with remapped loci"
#'
#' @examples
summarizeLoci = function(eid,alFile,normEset,wdir,sumFUN) {
  normEset = as.data.frame(normEset)
  PNtoLOCUS = as.data.frame(matrix(unlist(lapply(readLines(alFile),getMapping)),byrow=T,ncol=2),stringsAsFactors=F)
  colnames(PNtoLOCUS) = c("Locus","PN")

  expIDs = colnames(normEset)
  normEsetRAC = rownames_to_column(normEset,var="PN")

  locusData = dplyr::right_join(PNtoLOCUS,normEsetRAC) %>%  dplyr::select(-PN) %>% group_by(Locus) %>%  dplyr::summarise_all(sumFUN) %>% dplyr::filter(!is.na(Locus))
  locusData = tibble::column_to_rownames(locusData,var="Locus")


  return(locusData)
}


#' summarise all transcript intensities for one gene
#'
#' @param eset
#' @param sFun
#'
#' @return
#' @export
#'
#' @examples
summarizeTranscripts =function(eset,sFun=mean){
  esetNames = tibble::rownames_to_column(eset,var='transcriptID')
  gEset = esetNames %>% rowwise() %>% dplyr::mutate(geneID = splitTranscriptID(transcriptID)) %>%
    dplyr::select(-transcriptID) %>% dplyr::group_by(geneID) %>%
    dplyr::summarise_all(sFun)
  rEset = tibble::column_to_rownames(gEset,var="geneID")
  return(as.data.frame(rEset))
}

#' Split transript id into gene id and transcript number
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
splitTranscriptID = function(x) {
  s = strsplit(x,'\\.')
  l=length(s[[1]])
  paste0(s[[1]][1:l-1],collapse=".")
}
