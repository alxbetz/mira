convertDualChannelToSingleChannelDataStructure = function(rawData,fType) {
  Rdat = rawData$R
  Rbdat = rawData$Rb
  Gdat = rawData$G
  Gbdat = rawData$Gb



  nFiles = dim(Rdat)[2]
  nProbes = dim(Rdat)[1]

  E = matrix(0,ncol = nFiles *2,nrow =  nProbes)
  Eb = matrix(0,ncol = nFiles *2,nrow = nProbes )
  print(nFiles)
  print(dim(E))
  print((nFiles+1):(nFiles*2))

  E[,1:nFiles] = Rdat
  E[,(nFiles+1):(nFiles*2)] = Gdat
  Eb[,1:nFiles] = Rbdat
  Eb[,(nFiles+1):(nFiles*2)] = Gbdat


  rawDataSingle <- list(E = E, Eb = Eb, targets = rbind(paste0(rawData$targets,"_red"), paste0(rawData$targets,"_green")), source = fType, genes = rawData$genes)

  class(rawDataSingle) <- "EListRaw"

  return(rawDataSingle)

}


#' Title
#'
#' @param filenames
#' @param fType
#' @param dualChannel
#' @param bgmethod
#' @param treatDualChannelAsSingleChannel
#'
#' @return
#' @export
#'
#' @examples
normalizeLimma = function(filenames,fType,dualChannel,bgmethod="normexp",treatDualChannelAsSingleChannel=F,verbose=F) {
  if(verbose) {
    print("STARTING PREPROCESSING")
  }

  if(fType == "agilent") {
    idCol =  "ProbeName"
    nmethod="quantile"
  } else if (fType =="genepix") {
    idCol = "ID"
    nmethod="Aquantile"
  } else {
    nmethod="quantile"
  }

  if(treatDualChannelAsSingleChannel) {
    nmethod="quantile"
  }
  if(verbose) {
    print(nmethod)
    print(fType)
    print(c(idCol))
    print(filenames)

  }

  rawData= limma::read.maimages(filenames, source = fType, annotation = c(idCol),green.only=!dualChannel)

  if(treatDualChannelAsSingleChannel) {
    rawData = convertDualChannelToSingleChannelDataStructure(rawData,fType)
  }


  bgData <- limma::backgroundCorrect(rawData,method=bgmethod)
  normData = limma::normalizeBetweenArrays(bgData, method = nmethod)
  if(dualChannel && !treatDualChannelAsSingleChannel) {
    normEset =  normData$M
    print("Treating as Dual Channel")
  } else {

    print("Treating as Single Channel")
    normEset = normData$E
  }

  rownames(normEset) <- normData$genes[,idCol]
  return(normEset)
}
