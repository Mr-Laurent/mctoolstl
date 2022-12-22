#' Get metacells names from the heatmap of the metacells object
#'
#' 
#'
#' @param data Metacell z matrix
#' @param metadata metadata with cells corresponding to metacells
#' @param HMkc Number of k groups of columns
#' @param HMkr Number of k groups of rows
#' @param rmclu list of clusters to remove, ex : "1,2,4"
#' @param name Name of the file
#' @return A matrix of the infile
#' @export
#' 
#' 



cells_from_mc_hm<-function(data=NULL,metadata=NULL,HMkc=40,HMkr=0,rmclu=NULL,name=NULL){

  # selectrmclu="2,3,4,5,6,37"

  hm_zmod_rm<-data
  print(paste0("computing heatmap of ",HMkc,"cols and ",HMkr,"rows"))
  set.seed(42)
  ht = draw(Heatmap(t(hm_zmod_rm), column_title = "Module clustering based on corrected z score", name = "mat",
                    row_names_gp =gpar(fontsize=rel(400/length(colnames(hm_zmod_rm)))),show_row_dend = TRUE,
                    row_dend_reorder = F, cluster_rows = T, row_order=NULL,
                    column_dend_reorder = F, cluster_columns = T,
                    column_names_gp = gpar(fontsize=rel(400/length(rownames(hm_zmod_rm)))),column_km = HMkc,row_km = HMkr, use_raster = T) )

  metacells_kclust<- column_order(ht) 
  metacells_kclustrow<- row_order(ht) 
  
  htmxord<-as.matrix(t(hm_zmod_rm))[ as.vector(unlist(metacells_kclustrow)) , as.vector(unlist(metacells_kclust)) ]
  htpt1<-heatmaply(htmxord, dendrogram = "none",limits=range(-2,2),colors = "RdBu",Rowv=F,Colv=F,legendgroup="1st",showlegend = T,coloraxis = 'coloraxis',
                   scale_fill_gradient_fun =scale_fill_gradient2(low = "blue",mid = "white", high = "red", limits=c(-2, 2), oob=squish) )

  # Problem if only one metacell is in the cluster : it doesn't read as a matrix !! 
  # Use rownames(hm_zmod_rm)[metacells_kclust[[i]]]
  # instead of 
  #     rownames(hm_zmod_rm[metacells_kclust[[i]],])
  for (i in 1:length(metacells_kclust)){
    if (i == 1) {
      clu <- t(t(rownames(hm_zmod_rm)[metacells_kclust[[i]]]))
      out <- cbind(clu, paste("cluster", i, sep=""))
      colnames(out) <- c("Metacell", "Cluster")
    } else {
      clu <- t(t(rownames(hm_zmod_rm)[metacells_kclust[[i]]]))
      clu <- cbind(clu, paste("cluster", i, sep=""))
      out <- rbind(out, clu)
    }
  }
  rm(clu)
  out<-as.data.frame(out)
  
  rm_clu<-paste0("cluster",unlist(strsplit(rmclu,split=",")))
  Mc_to_subset<- out$Metacell[which(out$Cluster%in%rm_clu)]
  print("metacells identified")
  metadata$mc2<-paste0("mc",metadata$mc)
  write.csv(x = metadata$names[which(metadata$mc2%in%Mc_to_subset)], file=paste0("cells_to_rm_",name,".csv") ,quote = F,row.names = T)
  
  print("csv written")
  #save(Mc_to_exclude,file=paste0("mc_sub_",name,".rds"))
  #print(paste0("saved cell names in mc_sub_",name,".rds"))
}






#' Get the gene to use in module analysis and plot them
#' They are the most variable+expressed genes, following loess curve
#' 
#'
#' @param varmean_df dataframe of mean and variance of the genes
#' @param vm_MeanThres threshold for the log(mean) of expression
#' @param vm_VarmeanThres threshold for the curve start
#' @param x1min min value for the plot
#' @param x1max max value for the plot
#' @return TRUE/FALSE genes kept
#' @export
#' 
#' 



vm_genes_in_mod<-function(varmean_df=NULL,
                          vm_MeanThres=-1,
                          vm_VarmeanThres=1,
                          x1min=-1,
                          x1max=3.5){
  #Compute loess curve to get the variable genes
  x=log10(varmean_df$m)
  breaks=seq(min(x),max(x),.2)
  lv=log2(varmean_df$v/varmean_df$m)
  # erreur en metacell : des levels eleves n'avaient pas de valeurs dedans, 
  # (a cause d'un outlier: MALAT1, 7k en moy vs 1k pour les autres) donc ajouter "drop=T" a split
  z=sapply(split(lv,cut(x,breaks)),min,na.rm=T)
  maskinf=is.infinite(z)
  z=z[!maskinf]
  b=breaks[-length(breaks)]
  b=b[!maskinf]
  # En fait si on drop les levels, on ne drop pas dans breaks donc toujours probleme.
  # SOLUTION : VIRER MALAT1
  lo=loess(z~b)
  #Plot of the loess curve for variable gene selection 
  plot(log10(varmean_df$m),log2(varmean_df$v/varmean_df$m),xlab="Log10(mean)",ylab="log2(var/mean)",panel.first=grid())
  x1=seq(x1min,x1max,l=100)
  lline=predict(lo,newdata =x1)
  lines(x1,lline+as.numeric(vm_VarmeanThres),col=2)
  abline(v=vm_MeanThres,col=2)
  lline2=predict(lo,newdata =log10(varmean_df$m))
  
  
  ## Selection of the genes 
  
  # Use the same threshold but this time to subset the genes to keep from the matrix
  geneModuleMask<-log10(varmean_df$m)>as.numeric(vm_MeanThres)&log2(varmean_df$v/varmean_df$m)>lline2+as.numeric(vm_VarmeanThres)
  return(geneModuleMask)

}