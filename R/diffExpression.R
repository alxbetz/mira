
#' Differential Expression
#'
#' @param x
#' @param eset
#' @param annotation
#' @param wdir
#' @param volcanoplots
#'
#' @return
#' @export
#'
#' @examples
diffExpression = function(x,eset,annotation,wdir,volcanoplots=F) {
  compString=x[['contrastString']]
  experimentID = x[["eid"]]
  nContrasts = length(unique(annotation$Group))

  if (compString == "simple" && nContrasts==1) {
    print('simple only')
    fit = limma::lmFit(eset)
    fit = limma::eBayes(fit)
    fit2 = fit


  } else {
    fac = factor(annotation$Group, levels = unique(annotation$Group))
    design = model.matrix(~0 + fac)
    colnames(design) = levels(fac)

    if(compString == "simple") {
      rmax = dim(design)[1]
      rownames(design) <- colnames(eset)[1:rmax]
    } else {
      rownames(design) <- colnames(eset)
    }


    fit <- limma::lmFit(eset, design)
    if(compString=="simple") {
      print("simple compstring processing")
      cstring = paste(unique(annotation$Group),"VScontrol=",unique(annotation$Group),sep="")
      contrast.matrix <- limma::makeContrasts(contrasts=cstring,levels=design)
    } else {
      print("standard compstring:")
      print(compString)
      print(unique(annotation$Group))
      contrast.matrix <- limma::makeContrasts(contrasts=strsplit(compString,",")[[1]],levels=design)
    }

    print(contrast.matrix)
    fit2 <- limma::contrasts.fit(fit, contrast.matrix)
    fit2 <- limma::eBayes(fit2)

  }



  dt = limma::decideTests(fit2,adjust.method="fdr")
  #go = goana(fit2,coef=2,species="At")



  diffE <- list()
  lfct = as.double(x['lfcDE'])
  pvt = as.double(x['pvalDE'])


  if (compString == "simple" && nContrasts ==1) {
    #print(topTable(fit2,number=100,adjust="fdr"))
    diffE[[1]] <- limma::topTable(fit2,number=1000000,adjust="fdr")
    diffE[[1]]$locusID <- rownames(diffE[[1]])
    names(diffE)[1] <- c("trVSctrl")
    write.csv(diffE[[1]],file =file.path(wdir, "tt_simple.csv"))

    if(volcanoplots) {
      vplot = customVolcanoPlot(diffE[[1]],"simple contrast",pvt,lfct)
      ggplot2::ggsave(file.path(wdir, "tt_volcano.pdf"),plot=vplot)
    }

    #   plot(diffE[[1]][["logFC"]],-unlist(diffE[[1]]["adj.P.Val"]),xlim = c(-2,2),ylim=c(-1.1,0.1))
    #   abline(h = -pvt)
    #   abline(v = -lfct)
    #   abline(v = lfct)



  } else {
    print("N contrasts")
    print(ncol(fit2$contrasts))
    for (d in 1:ncol(fit2$contrasts)){
      diffE[[d]] <- topTable(fit2,coef=colnames(fit2$contrasts)[d],number=1000000,adjust="fdr")
      diffE[[d]]$locusID <- rownames(diffE[[d]])
      names(diffE)[d] <- colnames(fit2$contrasts)[d]
      #remove brackets from contrasts to avoid invalid file names
      ttFname = gsub("[()\\/]","",colnames(fit2$contrasts)[d])
      write.csv(diffE[[d]],file=file.path(wdir, paste("tt_",ttFname,".csv",sep="")))
      vplot = customVolcanoPlot(diffE[[d]],ttFname,pvt,lfct)
      ggsave(file.path(wdir, paste("tt_volcano_",ttFname,".pdf",sep="")),plot=vplot)
      #pdf(file.path(wdir, paste("tt_volcano_",colnames(fit2$contrasts)[d],".pdf",sep="")))
      #customVolcanoPlot(diffE[[d]])
      #plot(unlist(diffE[[d]]["logFC"]),-unlist(diffE[[d]]["adj.P.Val"]),xlim = c(-2,2),ylim=c(-1.1,0.1))
      #abline(h = -pvt)
      #abline(v = -lfct)
      #abline(v = lfct)
      #dev.off()
    }

  }



  # print(head(diffE))
  # print(length(diffE))
  # c1<- list()
  # for (i in 1:length(names(diffE))){
  #   c1[[i]] <- diffE[[i]][which(diffE[[i]][,"adj.P.Val"] <= pvt),"logFC"]
  # }
  # names(c1) <- names(diffE)
  #

  # make list with the locusID,FC and pval of diff. expressed genes with pVal <= 0.05 and 1.5 log2 FC threshold

  diffEsig<- list()
  print("diffE")
  print(length(diffE))
  #print(diffE)
  #print(dim(diffE[[1]]))
  for (i in 1:length(names(diffE))){
    diffEsig[[i]] <- diffE[[i]][which(diffE[[i]][,"adj.P.Val"] <= pvt & abs(diffE[[i]][,"logFC"]) >= lfct ),c("locusID","logFC")]
  }
  names(diffEsig) <- names(diffE)



  for (n in names(diffEsig)){
    colnames(diffEsig[[n]])[2]<-paste(n,'logFC',sep='_')
  }

  summary.diffE <- Reduce(function(x,y) merge(x,y, all=T,by.x='locusID',by.y='locusID'),diffEsig, accumulate=F)

  summary.diffE[is.na(summary.diffE)] <- 0
  rownames(summary.diffE) <- summary.diffE$locusID
  summary.diffE$locusID <- NULL

  colnames(summary.diffE) <- gsub("_logFC","",colnames(summary.diffE))
  write.csv(summary.diffE,file=file.path(wdir, paste("diffE_",experimentID,".csv",sep="")))
  return(dt)
}



fcrosHelper = function(meid,eset,annotation,contrast,verbose){
  if(verbose) {
    print("Processing id:")
    print(meid)
  }
  esetNames = tibble::rownames_to_column(eset,"idname")


  if(x[['channelCount']] == 2) {
    if(verbose){
      print("Handling Dual channel")
      print("Contrast is")
      print(contrast)
    }

    groupStrings = data.frame(trt=unique(annotation$Group),ctrl=paste0(unique(annotation$Group),"_control"))
    contrastList = apply(groupStrings,1,function(x) paste(x['trt'],x['ctrl'],sep="-"))
    contrast = Reduce(function(x,y) { paste(x,y,sep=",")},contrastList)

    if(verbose){
      print("DualToSingleChannelProcessing")
      print("New ContrastString is")
      print(contrast)
    }


    annotation = annotationDualToSingleChannel(annotation,controlColor)

  }

  sp = strsplit(contrast,",")[[1]]


  #parse control and treatment groups from contrast strings
  spf = lapply(spNP,function(x) { strsplit(x,"-")[[1]]})
  cdf =  do.call(rbind, spf)
  tr = unique(cdf[,1])
  co = unique(cdf[,2])
  #browser()

  #trNames = annotation %>%  filter(group==tr) %>% dplyr::pull(ID)
  #coNames = annotation %>%  filter(group==co) %>% dplyr::pull(ID)
  trNames = colnames(eset)[annotation$Group %in% tr]
  coNames = colnames(eset)[annotation$Group %in% co]


  fcrosResult = fcros(esetNames,cont=unlist(coNames),test=unlist(trNames))
  fcrosResult$idnames = gsub("\"","",fcrosResult$idnames)


  #get AVG expression
  allNames = c(trNames,coNames)
  #print(allNames)
  #print(head(eset))
  fcrosResult$meanExpression = rowMeans(eset[,allNames])
  fcrosResult
}
