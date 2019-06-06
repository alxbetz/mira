
#' Calculate GO enrichment for one contrast
#'
#' @param dtCol
#' @param transcript2GO
#' @param nTopNodes
#'
#' @return
#' data frame of enrichment results for one contrast
#' @export
#'
#' @examples
goEnrichmentHelper = function(dtCol,transcript2GO,nTopNodes=20) {
  #print(dtCol)
  geneList = as.vector(dtCol)
  #names(geneList) = names(dtCol)
  #print(geneList)
  geneListFactor = factor(as.integer(geneList != 0))
  names(geneListFactor) = names(dtCol)

  print(head(geneListFactor))
  nDEgenes = sum(abs(geneList))
  print(paste("Using",as.character(nDEgenes),"DE genes" ))
  print(head(transcript2GO))
  #print("browser")
  #browser()
  go.obj <- new("topGOdata", ontology='BP'
                , allGenes = geneListFactor
                , annot = annFUN.gene2GO
                , gene2GO = transcript2GO
  )
  resultFisher <- topGO::runTest(go.obj, algorithm = "classic", statistic = "fisher")
  resultKS <- topGO::runTest(go.obj, algorithm = "classic", statistic = "ks")
  resultKS.elim <- topGO::runTest(go.obj, algorithm = "elim", statistic = "ks")

  allRes <- topGO::GenTable(go.obj, classicFisher = resultFisher,
                            classicKS = resultKS, elimKS = resultKS.elim,
                            orderBy = "elimKS", ranksOf = "classicFisher", topNodes = nTopNodes)


  return(allRes)
}

#' Run GO enrichment
#'
#' @param dt
#' @param transcript2GO
#' @param nTopNodes
#' @param verbose
#'
#' @return
#' list of GO enrichments, one element per contrast
#' @export
#'
#' @examples
goEnrichment = function(dt,transcript2GO,nTopNodes,verbose) {

  hasDEGenes = apply(dt,2,function(x) sum(abs(x)) > 0)
  #browser()
  if(sum(hasDEGenes) == 0){
    return(NULL)
  } else {
    if(verbose) {
      print("differentiallyExpressedGenes")
      print(hasDEGenes)
    }


    if(sum(hasDEGenes) == 1){
      return(list(goEnrichmentHelper(dt[,hasDEGenes],transcript2GO=transcript2GO,nTopNodes=nTopNodes)))
    } else {
      if(verbose) {
        print(head(dt))
        print(head(dt[,hasDEGenes]))
      }

      return(apply(dt[,hasDEGenes],2,goEnrichmentHelper,transcript2GO=transcript2GO,nTopNodes=nTopNodes))
    }
  }


}

#' Write a single enrichment table to a file
#'
#' @param idx
#' @param resList
#' @param wdir
#'
#' @examples
writeEnrichment = function(idx,resList,wdir) {
  fname = paste(names(resList)[idx],"_goe.tsv",sep="")
  write.table(resList[[idx]],file=file.path(wdir,fname))
}

#rlang + dplyr imports
fcrosGoEnrichmentHelper = function(x,fcRes,nTopNodes=100) {

  species= x[['species']]
  goAnno = annoList[[species]]
  print(head(fcRes))
  xranks = fcRes %>% dplyr::filter(key=="ri")

  require(stringr)
  require(tidyr)
  require(rlang)

  if(species == 'drerio') {
    xranks = xranks %>% tidyr::separate(idnames,c("idnames","transcriptID"),sep="\\.") %>% dplyr::group_by(idnames) %>% dplyr::summarise(median = median(!!rlang::sym(colnames(x)[3])))

  }

  maxIDX = dim(xranks)[2]

  geneList = unlist(xranks %>% dplyr::select(maxIDX))
  gln = sapply(as.character(xranks$idnames),function(x) gsub("\"","",x) )
  if(verbose) {
    names(gln)

    print(length(gln))
    print(length(geneList))
    print(head(geneList))
  }

  names(geneList) = gln
  if(verbose) {print(head(geneList))}

  nx = topDiffGenes(geneList)
 if(verbose) { print(sum(nx)) }


  go.obj <- new("topGOdata", ontology='BP'
                , allGenes = geneList
                , geneSel = topDiffGenes
                , annot = annFUN.gene2GO
                , gene2GO = goAnno
                , nodeSize = 15
  )




  #resultFisher <- topGO::runTest(go.obj, algorithm = "classic", statistic = "fisher")
  #resultKS <- topGO::runTest(go.obj, algorithm = "classic", statistic = "ks")
  resultKS.elim <- topGO::runTest(go.obj, algorithm = "elim", statistic = "ks")

  allRes <- topGO::GenTable(go.obj, elimKS = resultKS.elim,
                            orderBy = "elimKS", topNodes = nTopNodes)


  allGO = topGO::genesInTerm(go.obj)
  if(species=="drerio") {

    expTable = fcRes %>% dplyr::filter(key=="FC")  %>% tidyr::separate(idnames,c("idnames","transcriptID"),sep="\\.") %>% dplyr::group_by(idnames) %>% dplyr::summarise(median = median(!!rlang::sym(colnames(x)[3])))

    FClist = unlist(expTable[,2])
    names(FClist) = sapply(as.character(unlist(expTable[,1])),function(x) gsub("\"","",x) )

  } else {
    expTable = fcRes  %>% dplyr::filter(key=="FC") %>% dplyr::select(-key )
    FClist = unlist(expTable[,2])
    names(FClist) = sapply(as.character(unlist(expTable[,1])),function(x) gsub("\"","",x) )
  }








  goGeneList = allGO[allRes$GO.ID]



  FCperGO = lapply(goGeneList,function(genes,FClist) { FClist[genes]  },FClist=FClist)

  medianFC = lapply(FCperGO,median)
  meanFC = lapply(FCperGO,mean)
  maxFC = lapply(FCperGO,max)
  minFC = lapply(FCperGO,min)
  allResFinal = allRes %>% mutate(medianFC= medianFC ,meanFC = meanFC, minFC =minFC,maxFC=maxFC)
  list(enrichment=allResFinal,genesPerGO=goGeneList,goData = go.obj)

}



