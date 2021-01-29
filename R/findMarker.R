#' Identify marker genes in the given expression matrix.
#'
#' @title findMarker
#'
#' @param mat A (m x n) data matrix of gene expression measurements of individual cells with rows representing genes and columns representing cells.
#' @param cluster A cluster output (vector) from clustering 
#' @param cluster_id Id (numeric) cluster group in `cluster`
#'
#' @return A data frame with gene's ordered by their significance.
#'
#' 
#' @importFrom Hmisc rcorr
#' @import MAST
#'
#' @export

findMarker <- function(mat, cluster, cluster_id){
  
  #group 1 is the cluster of interest
  group1 = which(cluster == cluster_id)  
  group1 = colnames(mat)[group1]  #find out which cell belongs to group 1
  
  group_label <- ifelse(colnames(mat) %in% group1, "group1", "group2")  #if the cell is in group 1 ,create a label called group 1, else create a label called group 2
  
  #create a dataframe , with cell name and the group that this cell belongs to 
  group_label <- data.frame(group_label, colnames(mat))
  rownames(group_label) <- group_label[,2]
  colnames(group_label) <- c("group","wellKey")
  #MAST requires a wellkey (which essentially is the name of the cell)
  
  
  #create a dataframe, with gene name 
  gene_label <- data.frame(rownames( mat))
  colnames(gene_label) <- "primerid"
  rownames(gene_label) <- gene_label[, 1]
  #MAST requires a primerid (which essentially is the gene name)
  
  #create a MAST object, with the group id of each cell and the gene name 
  sca <- MAST::FromMatrix(
    exprsArray = mat,
    cData = group_label,
    fData = gene_label
  )
  
  #this is the differential expression model (called a hurdle model)
  #(see https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html#4_differential_expression_using_a_hurdle_model)
  zlmCond <- MAST::zlm(~group, sca) 
  summaryCond <- summary(object = zlmCond, doLRT = 'groupgroup2') 
  
  #now extract the information 
  summaryDT <- data.frame(summaryCond$datatable)
  return_val <-data.frame(summaryDT[summaryDT[, "component"] == "H", 1], summaryDT[summaryDT[, "component"] == "H", 4])
  #We select the "H" rows, the "H" means we want hurdle P values
  #extract both the gene name and the P value associated with this gene
  colnames(return_val) <- c("gene", "P_value")
  
  #order by the p value , from the most significant 
  return_val <- return_val[order(return_val$P_value), ]
  
  return (return_val)
  
}