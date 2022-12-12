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
  rm_clu<-paste0("cluster",unlist(strsplit(rmclu,split=",")))
  hm_zmod_rm<-data
  
  set.seed(42)
  print(paste0("computing heatmap of ",HMkc,"cols and ",HMkr,"rows"))
  ht = Heatmap(t(hm_zmod_rm), column_title = "Module clustering based on corrected z score", name = "mat",
                    row_names_gp =gpar(fontsize=rel(400/length(colnames(hm_zmod_rm)))),show_row_dend = TRUE,
                    row_dend_reorder = F, cluster_rows = T, row_order=NULL,
                    column_dend_reorder = F, cluster_columns = T,
                    column_names_gp = gpar(fontsize=rel(400/length(rownames(hm_zmod_rm)))),column_km = HMkc,row_km = HMkr, use_raster = T) 

  htmxord<-as.matrix(t(hm_zmod_rm))[ as.vector(unlist(row_order(ht))) , as.vector(unlist(column_order(ht))) ]
  htpt1<-heatmaply(htmxord, dendrogram = "none",limits=range(-2,2),colors = "RdBu",Rowv=F,Colv=F,legendgroup="1st",showlegend = T,coloraxis = 'coloraxis',
                   scale_fill_gradient_fun =scale_fill_gradient2(low = "blue",mid = "white", high = "red", limits=c(-2, 2), oob=squish) )
  metacells_kclust<- column_order(ht) 
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
  
  Mc_to_exclude<- out$Metacell[which(out$Cluster%in%rm_clu)]
  print("metacells identified")
  df$mc<-paste0("mc",df$mc)

  write.csv(x = df$names[which(df$mc%in%Mc_to_exclude)], file=paste0("cells_to_rm_",name,".csv") ,quote = F,row.names = T)
  print("csv written")
  #save(Mc_to_exclude,file=paste0("mc_sub_",name,".rds"))
  #print(paste0("saved cell names in mc_sub_",name,".rds"))
}




