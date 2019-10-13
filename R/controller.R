
#' MIcroarray Re-Analysis
#'
#' @param x
#' @param pipedir
#' @param pipeComponents
#' @param summarizeFunction
#' @param treatDualChannelAsSingleChannel
#' @param controlColor
#' @param verbose
#' @param forceAnnotation
#' @param summarizeToGeneLevel
#' @param deMethod
#'
#' @return
#' @export
#'
#' @examples
mira = function(x,
                pipedir,
                pipeComponents='ande',
                summarizeFunction=mean,
                treatDualChannelAsSingleChannel=F,
                controlColor="green",
                verbose=FALSE,
                forceAnnotation=FALSE,
                summarizeToGeneLevel=TRUE,
                deMethod = 'limma'
                ) {

  #create working directory

  #remove trailing slashes in pathname
  pipedir = gsub('/+$','',pipedir)
  wdir = file.path(pipedir,"output",x[['eid']])
  if(verbose) {
    print("data passed to function")
    print(x)
    print(paste("creating",wdir))
  }
  dir.create(wdir)

  #load annotation file
  annotation = readAnnotation(pipedir,x)

  #Download array image reader files
  # if(x[['database']] != 'AE') {
  #   if(all(file.exists(file.path(pipedir,"input",x['eid'],annotation$File)))) {
  #     file.copy(file.path(pipedir,"input",x['eid'],annotation$File),file.path(pipedir,"output",x['eid'],annotation$File))
  #   } else {
  #     downloadDataSet(x,wdir)
  #   }
  # }

  downloadDataSet(x,wdir)



  if(verbose) {
    print(annotation)
  }

  cFile = file.path(pipedir,"input",x['eid'],"contrast.txt")
  if(file.exists(cFile)) {
    contrast = trimws(readLines(cFile,n=1))

    x[['contrastString']] = contrast
    if(verbose) {
      print("Contrast is")
      print(contrast)
    }

  }


  #Process and normalize
  if(x[['database']] != 'AE') {
    filenames = unique(file.path(wdir,annotation$File))

    #if files are gzipped, add gz suffix to filenames
    if(all(file.exists(
      paste(filenames,".gz",sep="")
    ))) {
      filenames = paste(filenames,".gz",sep="")
    }
    print("CHECKING FILETYPE")
    fType= scanFileType(filenames[1])

    print("FILETYPE")
    print(fType)
  }

  if(treatDualChannelAsSingleChannel) {
    groupStrings = data.frame(trt=unique(annotation$Group),ctrl=paste0(unique(annotation$Group),"_control"))
    contrastList = apply(groupStrings,1,function(x) paste(x['trt'],x['ctrl'],sep="-"))
    x[['contrastString']] = Reduce(function(x,y) { paste(x,y,sep=",")},contrastList)

    annotation = annotationDualToSingleChannel(annotation,controlColor=controlColor)

  } #else {
  #fType = 'loadedFromArrayExpress'
  #}



  if(fType=="celera") {

    sampleNameList = annotation$Name


    celObj = processCeleraFiles(filenames,x[['cDNAFile']],wdir,pipedir,sampleNameList)
    normEset = celObj[['normEset']]
    alFile = celObj[['alFile']]
    normData = celObj[['normData']]
  } else {

    if(grepl("a",pipeComponents)){
      print("Aligning probe sequences to Reference")
      #function(cDNAFile,wdir,pipedir,iDir,prefix)
      alFile = mapProbes(x,wdir,pipedir)
    } else {
      alFile = file.path(wdir,"alignment")
    }
    if(grepl("n",pipeComponents)){
      print("Normalizing Counts")
      if(x[['database']] == 'AE') {
        normEset = processArrayExpress(x[['eid']])
      } else {
        normEset = normalizeLimma(filenames,fType,x[['channelCount']]==2,treatDualChannelAsSingleChannel=treatDualChannelAsSingleChannel,verbose = verbose)
      }

      normSumEset = summarizeLoci(x[['eid']],alFile,normEset,wdir,summarizeFunction)
      write.table(normSumEset, file.path(wdir, paste(x[['eid']],"meanNormWithin.txt",sep="_")), sep = '\t', quote = FALSE)

    }



    if(grepl("d",pipeComponents)) {

      if(!exists('normSumEset')) {
        normSumEset = read.table(file.path(wdir, paste(x[['eid']],"meanNormWithin.txt",sep="_")),quote = "")
      }

      if(summarizeToGeneLevel==TRUE)
      {
        normSumEset=summarizeTranscripts(normSumEset)
      }
      if(deMethod == 'limma') {
        dtTable = diffExpression(x,normSumEset,annotation,wdir)
        if(verbose) {
          print(head(dtTable))
        }
      } else if(deMethod == 'fcros') {
        fcResList = fcrosWrapper(x,normSumEset,annotation,controlColor,TRUE)
      } else {
        stop("invalid deMethod parameter, use one of ('limma,'fcros)")
      }


      #write.csv(summary.diffE,file=file.path(wdir, paste("diffE_",experimentID,".csv",sep="")))
    }

    if(grepl("e",pipeComponents)) {
      cFile = file.path(pipedir,"input",x['eid'],"contrast.txt")
      if(file.exists(cFile)) {
        contrast = trimws(readLines(cFile,n=1))

        x[['contrastString']] = contrast
        print("Contrast is")
        print(contrast)
      }
      browser()
      annoObjFname = paste(x[['eid']],"t2go.rds",sep="_")
      if(file.exists(annoObjFname) && forceAnnotation == FALSE) {
        transcript2GO = readRDS(annoObjFname)
      } else {
        transcript2GO= loadAnnotation((x[['species']]))
        saveRDS(transcript2GO,file=annoObjFname)
      }
      if(deMethod =='limma') {
        ha = rownames(dtTable) %in% names(transcript2GO)
      } else if(deMethod =='fcros') {
        ha = fcResList[[1]]$idnames %in% names(transcript2GO)
      }

      fa = sum(ha) / length(ha)

      if(verbose) {
        print("Fraction of IDs without GO annotation")
        print(fa)
      }

      if(fa < 0.1) {
        print(head(rownames(dtTable)))
        print(head(names(transcript2GO)))
        stop("IDspace of differential Expression Table and annotation not compatible")
      }
      if(deMethod == "limma") {
        resList = goEnrichment(dtTable,transcript2GO,nTopNodes=20,verbose=verbose)
      } else if(deMethod == 'fcros') {
        resList = fcrosGoEnrichmentWrapper(x,fcResList,transcript2GO,nTopNodes=100,verbose=F)
      }

      if(is.null(resList)) {
        stop(paste0("No DE genes found for ",x[['eid']],". Stopping Analysis"))
      }

      if(deMethod == "limma") {
        lapply(seq_along(resList),writeEnrichment,resList=resList,wdir=wdir,deMethod=deMethod)
      } else if(deMethod == 'fcros') {
        saveRDS(resList,file = file.path(wdir,paste0(deMethod,"_goe.rds")))
      }

    }

  }


  list(eMatrix = normSumEset,enrichment=resList)

}
#Summarize
#return(summarizeLoci(x[['eid']],alFile,normEset,wdir,summarizeFunction))
#return(list(eid=unlist(x['eid']),alFile=alFile,normEset=normEset,wdir=wdir,normData=normData))
#}

# 1. read annotation file
# 2. read microarray image data
# 3. realign probes
# 4. preprocess arrays


# pipeline management:
# check if a step was successfull
# optional: save intermediate results of step
# optional: load intermediate results


#' setup_folders
#' set up pipeline folders
#'
#' @param pipepath
#' path to pipeline folder
#'
#' @return
#' @export
#'
#' @examples
setup_folders = function(pipepath) {
  dir.create(file.path(pipepath,'input'))
  dir.create(file.path(pipepath,'output'))
  dir.create(file.path(pipepath,'genomicDB'))
  dir.create(file.path(pipepath,'platforms'))
}

