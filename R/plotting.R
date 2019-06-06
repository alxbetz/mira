
customVolcanoPlot = function(x,title,pval,lfc) {
  xlimv = max(abs(min(x[['logFC']])),max(x[['logFC']]))
  xdf = as.data.frame(x)
  print(pval)
  print(lfc)
  if(pval && lfc) {
    xdf = xdf %>% mutate(DE=adj.P.Val < get("pval") & abs(logFC) > lfc)
    print(table(xdf$DE))
    print(head(xdf))
    hiLine = ggplot2::geom_vline(xintercept=lfc)
    loLine =  ggplot2::geom_vline(xintercept=-lfc)
    pvLine = ggplot2::geom_hline(aes(yintercept=-log10(pval)))
  } else if(pval) {
    xdf = xdf %>% dplyr::mutate(DE=adj.P.Val < get("pval"))

    pvLine = ggplot2::geom_hline(aes(yintercept=-log10(pval)), show.legend=TRUE)
  }


  p = ggplot(xdf) + geom_point(aes(x=logFC,y=-log10(adj.P.Val),color=DE)) +
    scale_x_continuous(limits = c(-xlimv, xlimv)) + theme_bw() + ggtitle(title)
  if(pval) {
    p = p + pvLine
  }
  if(lfc) {
    p = p + loLine + hiLine
  }
  p
}
