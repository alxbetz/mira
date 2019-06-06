
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
      if(summarizeToGeneLevel==TRUE && grepl("\\.",colnames(normSumEset)[1]))
      {
        normSumEset=summarizeTranscripts(normSumEset)
      }

      dtTable = diffExpression(x,normSumEset,annotation,wdir)
      if(verbose) {
        print(head(dtTable))
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
      annoObjFname = paste(x[['eid']],"t2go.Robj",sep="_")
      if(file.exists(annoObjFname) && forceAnnotation == FALSE) {
        load(annoObjFname)
      } else {
        transcript2GO= loadAnnotation((x[['species']]))
        save('transcript2GO',file=annoObjFname)
      }
      ha = rownames(dtTable) %in% names(transcript2GO)
      fa = sum(ha) / length(ha)

      if(verbose) {
        print("Fraction of IDs without GO annotation")
        print(fa)
      }

      if(fa < 0.6) {
        print(head(rownames(dtTable)))
        print(head(names(transcript2GO)))
        stop("IDspace of diffExp Table and annotation not compatible")
      }

      resList = goEnrichment(dtTable,transcript2GO,nTopNodes=20,verbose=verbose)
      if(is.null(resList)) {
        print(paste("No DE genes found for ",x[['eid']],". Stopping Analysis",sep=""))
      }
      lapply(seq_along(resList),writeEnrichment,resList=resList,wdir=wdir)
    }

  }


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


