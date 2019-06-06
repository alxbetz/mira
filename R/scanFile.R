
#' Reads the magic number from a binary file and identifies the microarray file type
#'
#' @param filename
#'
#'
#'
#'
#' @return
#'
#' @export
#'
#' @examples
scanFileType = function(filename) {
  print(filename)
  ## Handling binary Files
  #magic numbers
  #CELERA
  # 0000013b        00000100
  #gzipped
  # 08088b1f        4b8d
  #gzMagic = c("1f","8b","08","08","f8","38","8d","4b")
  gzMagic = c("1f","8b")
  #celMagic = c("3b","01","00","00","00","01","00","00")
  celMagic = c("3b","01")
  celMagic2 = c("5b","43")
  celMagic3 = c("40","00")

  magicN = as.character(readBin(filename,"raw",n=2,endian="little"))
  print("MagicNumber")
  print(magicN)
  if (all(magicN == gzMagic)) {
    tmpfile = tempfile()
    GEOquery::gunzip(filename,destname = tmpfile,remove=F)
    return(scanFileType(tmpfile))
  } else if(all(magicN==celMagic) || all(magicN==celMagic2) || all(magicN==celMagic3)) {
    return("celera")
  }

  ##handling text Files

  idStrings = c(
    genepix = "Type=GenePix Results 3",
    agilent = "Agilent Technologies",
    agilent2 = "Begin Measurement parameters"
  )
  headLines = readLines(filename,n=20)
  binIDX = sapply(idStrings,function(pat,lib) any(grepl(pat,lib)),lib=headLines)
  print(binIDX)
  print(any(binIDX))
  if(any(binIDX)) {
    if(names(idStrings[binIDX]) == "agilent2")
    {
      return("agilent")
    } else {
      return(names(idStrings[binIDX]))
    }

  } else {
    return("unknown")
  }

}
