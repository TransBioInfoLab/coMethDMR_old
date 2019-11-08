#' Get pre-annotated clusters of probes
#'
#' @param arrayType Type of array, 450k or EPIC
#'
#' @return list of clustered probes
#' @export
#'
#' @examples
#'  cluster <- getPredefinedCluster ("450k")
getPredefinedCluster <- function( arrayType = c("450k","EPIC")){
  arrayType <- match.arg(arrayType)
  files <- system.file("extdata", "", package = "coMethDMR")
  probes.cluster <- dir(files,
      pattern = paste0(arrayType,".*200.rds"),
      full.names = TRUE,
      recursive = TRUE)

  list.of.probes.cluster <- plyr::alply(probes.cluster,1,function(file){readRDS(file)})
  names(list.of.probes.cluster) <- gsub(".rds","",basename(probes.cluster))
  probes.cluster.all <- unlist(list.of.probes.cluster, recursive = FALSE)
  return(probes.cluster.all)
}
