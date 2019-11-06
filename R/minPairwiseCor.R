#' Calculate min Pairwise correlation for subregions probes
#'
#' @param betaCluster_mat matrix of beta values, with rownames = sample ids,
#'    column names = CpG ids. Note that the CpGs need to be ordered by their genomic positions,
#'    this can be accomplished by the \code{OrderCpGbyLocation} function.
#'
#' @param betaToM indicates if converting to mvalues before computing correlations
#'
#' @param minPairwiseCorr thershold for min pairwise correlation between any cpg the
#'    rest of the CpGs
#'
#' @param method correlation method, can be pearson or spearman
#'
#' @param probes.list list of regions with cpgs
#' 
#' @param subregions.annot A dataframe with probes, subregion that will be returned
#' with the minPairwise correlation and if probe should be kept.
#' 
#' @return A list with the list of probes passing the minPairwiseCorr and a dataframe
#' with the minPairwiseCorr for each subregion if subregions.annot was provided
#' @export
#'
#' @examples
#' data(betaMatrix_ex1)
#' probes <- colnames(betaMatrix_ex1)
#' minPairwiseCor(betaCluster_mat = betaMatrix_ex1, 
#'                betaToM = FALSE, 
#'                minPairwiseCorr = 0.1,
#'                probes.list = list("chr22:18531243-18531447" = probes),
#'                method = "pearson",
#'                subregions.annot = data.frame(probes,Subregion = rep(1,length(probes))))
#' minPairwiseCor(betaCluster_mat = betaMatrix_ex1, 
#'                betaToM = FALSE, 
#'                minPairwiseCorr = 0.2,
#'                probes.list = list("chr22:18531243-18531447" = probes),
#'                method = "pearson",
#'                subregions.annot = data.frame(probes,Subregion = rep(1,length(probes))))
minPairwiseCor <- function (betaCluster_mat,
                            betaToM = TRUE,
                            minPairwiseCorr = 0.2,
                            probes.list,
                            subregions.annot = NULL,
                            method = c("pearson", "spearman")) {
  
  
  ### Check that betaToM == "TRUE" only if betaCluster_mat has beta values ###
  
  if((min(betaCluster_mat,na.rm = TRUE) < 0 | max(betaCluster_mat > 1,na.rm = TRUE)) & betaToM == "TRUE") {
    message("The input methylation values are not beta values,
         if they are M values, 'betaToM' should be FALSE")
    return(NULL)
  }
  
  if (betaToM == "TRUE") {
    betaCluster_mat <- log2(betaCluster_mat / (1 - betaCluster_mat))
  } 
  
  ### Calculate minPairwiseCor for each CpGsSubregion
  minPairwiseCorlist <- lapply(probes.list, function(probes){
    betaCluster_mat[,probes] %>%
      cor(method = method) %>%
      min(na.rm = TRUE)
  })
  
  keepminPairwiseCor <- minPairwiseCorlist > minPairwiseCorr
  
  if(!is.null(subregions.annot)){
    if(is.null(probes.list)){
      subregions.annot$keepminPairwiseCor <- NA
      subregions.annot$minPairwiseCorlist <- NA
    } else if(all(subregions.annot$Subregion == 0)){
      subregions.annot$keepminPairwiseCor <- as.logical(keepminPairwiseCor)
      subregions.annot$minPairwiseCorlist <- as.numeric(minPairwiseCorlist)
    } else {
      subregions.annot$keepminPairwiseCor <- as.logical(keepminPairwiseCor[subregions.annot$Subregion])
      subregions.annot$minPairwiseCorlist <- as.numeric(minPairwiseCorlist[subregions.annot$Subregion])
    }
  }
  
  if(all(keepminPairwiseCor == FALSE)) {
    probes.list <- NULL
  } else {
    probes.list <- probes.list[keepminPairwiseCor]
  }
  return(list("subregions.annot" = subregions.annot,
              "probes.list.filtered" = probes.list))
}
