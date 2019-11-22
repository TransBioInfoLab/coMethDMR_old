#' Get pre-annotated clusters of probes
#'
#' @param arrayType Type of array, 450k or EPIC
#' @param clusterType Type of cluster of probes: "gene" Cpgs or "regions" Cpgs
#' clusters.
#' @return list of clustered probes
#' @export
#'
#' @examples
#'  cluster <- getPredefinedCluster ("450k")
getPredefinedCluster <- function(arrayType = c("450k","EPIC"),
                                 clusterType = c("regions","gene")){

  arrayType <- match.arg(arrayType)
  clusterType <- match.arg(clusterType)

  files <- system.file("extdata", "", package = "coMethDMR")

  probes.cluster <- dir(files,
                        pattern = paste0(arrayType,
                                         ifelse(clusterType == "regions",
                                                ".*200.rds",
                                                "*_CpGstoGene_min3CpGs.rds")),
                        full.names = TRUE,
                        recursive = TRUE)

  list.of.probes.cluster <- plyr::alply(probes.cluster,1,function(file){readRDS(file)})
  names(list.of.probes.cluster) <- gsub(".rds","",basename(probes.cluster))
  probes.cluster.all <- unlist(list.of.probes.cluster, recursive = FALSE)
  return(probes.cluster.all)
}
