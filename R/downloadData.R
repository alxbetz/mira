
#' Download a microarray dataset from ArrayExpress or GEO
#'
#' @param x
#' @param wdir
#'
#' @return
#' @export
#'
#' @examples
downloadDataSet = function(x,wdir) {
  if(x[['database']] =="AE") {
    downloadAESet(x,wdir)
  } else {
    downloadGEOSet(x,wdir)
  }
}

#' Download an ArrayExpress Dataset
#'
#' @param x
#' @param wdir
#'
#' @import ArrayExpress
#' @return
#'
#' @examples
downloadAESet = function(x,wdir) {
  try(dat <- getAE(x[['eid']],path=wdir))
  celFileList = list.files(wdir,pattern=".CEL",full.names=TRUE)
  if(length(celFileList) > 0) {
    print("Renaming Files")
    #browser()
    newCelFileNames = sapply(celFileList,function(x) { paste0(strsplit(x,"\\.CEL")[[1]][1],".cel")  })
    file.rename(celFileList,newCelFileNames)
    print(celFileList)
    print(newCelFileNames)
    print(list.files(wdir,pattern=".cel"))
    print("Done Renaming Files")
  }

}

#' Download a GEO dataset
#'
#' @param x
#' @param wdir
#'
#' @import utils
#' @return
#'
#' @examples
downloadGEOSet = function(x,wdir) {
  #getGEOSuppFiles(x['eid'],baseDir = wdir)
  if(nchar(x[['eid']]) == 7) {
    ftpLink = paste("ftp://ftp.ncbi.nlm.nih.gov/geo/series/",substr(x['eid'],1,4),"nnn/",x['eid'],"/suppl/",x['eid'],"_RAW.tar",sep="")
  } else {
    ftpLink = paste("ftp://ftp.ncbi.nlm.nih.gov/geo/series/",substr(x['eid'],1,5),"nnn/",x['eid'],"/suppl/",x['eid'],"_RAW.tar",sep="")
  }

  downloadFile=file.path(wdir,paste(x['eid'],"RAW.tar",sep="_"))
  if(!file.exists(downloadFile) || file.info(downloadFile)$size == 0) {
    download.file(ftpLink,downloadFile)
  }
  print("File path for download")
  print(downloadFile)
  untar(file.path(wdir,paste(x[['eid']],"RAW.tar",sep="_")),exdir=wdir)

  print("fileList")
  print(list.files(wdir))
}



#' Download a microarray platform file
#'
#' @param eid
#' @param wdir
#'
#' @return
#' full path filename of the platform file
#'
#' @examples
downloadPlatformFile = function(eid,wdir) {
  if(nchar(eid) == 7) {
    famURL = paste('ftp://ftp.ncbi.nlm.nih.gov/geo/series/',substr(eid,1,4),"nnn/",eid,"/soft/",eid,"_family.soft.gz",sep="")

  } else {
    famURL = paste('ftp://ftp.ncbi.nlm.nih.gov/geo/series/',substr(eid,1,5),"nnn/",eid,"/soft/",eid,"_family.soft.gz",sep="")
  }


  famFile = file.path(wdir,paste(eid,"_family.soft.gz",sep=""))
  #print(grepl(pfFile,"\.gz"))
  if(!file.exists(famFile)) {
    download.file(famURL,destfile =  famFile )
  }
  famObj = gzfile(famFile)
  famTxt = readLines(famObj)


  r = grep("Series_platform_id",famTxt,value=T)
  pid = strsplit(r," ")[[1]][3]

  pURL = paste("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",pid,"&targ=self&view=data&form=text",sep="")
  #echo 'wget "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${pid}&targ=self&view=data&form=text" -O "${platformDir}/${pid}/${pid}.txt"'
  pFile=file.path(wdir,paste(pid,"_platformFile.txt",sep=''))
  download.file(pURL,destfile = pFile)
  return(pFile)
}

#' Convert creinhardtii ensembl geneID to phytozome gene id
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
creEnsembl2phytozome = function(x) {
  paste0(
    'Cre',
    substr(x,7,8),
    '.',
    substr(x,10,15))
}


#' Loads GO annotation from Biomart
#'
#' @param specimen
#' @param annoFile
#'
#' @return
#' list of gene id to go mappings:
#' listnames are genes
#' list values are vectors of GO IDs
#' @export
#'
#'
#' @examples
loadAnnotation = function(specimen,annoFile=NULL) {

  if(specimen=="creinhardtii") {
    mart <- biomaRt::useMart(biomart = "plants_mart", dataset = "creinhardtii_eg_gene", host = 'plants.ensembl.org')
    TTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id", "go_id"), mart = mart)
    TTOGO$phytozome_gene_id = sapply(TTOGO$ensembl_gene_id,creEnsembl2phytozome)
    transcriptID2GO <- by(TTOGO$go_id, TTOGO$phytozome_gene_id, function(x) as.character(x))
    return(transcriptID2GO)
  } else {

    if(specimen=="athaliana") {
      mart <- biomaRt::useMart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host = 'plants.ensembl.org')
    } else {
      dataset=paste(specimen,"gene","ensembl",sep="_")
      print(dataset)
      mart <- biomaRt::useMart(biomart = "ensembl", dataset = dataset)
    }

    TTOGO <- biomaRt::getBM(attributes = c( "ensembl_transcript_id", "go_id"), mart = mart)
    transcriptID2GO <- by(TTOGO$go_id, TTOGO$ensembl_transcript_id, function(x) as.character(x))
    #TTOGO = TTOGO[TTOGO$go_id != '',]

    return(transcriptID2GO)
  }
}


