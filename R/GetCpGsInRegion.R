#' Extract probe IDs for CpGs located in a genomic region
#'
#' @param regionName_char character string with location information for one region in
#'    this format: "chrxx:xxxxxx-xxxxxx"
#' @param arrayType Type of array, 450k or EPIC
#'
#' @return vector of CpG probe IDs mapped to the genomic region
#' @export
#'
#' @importFrom tidyr separate %>%
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges subsetByOverlaps
#'
#' @examples
#'    GetCpGsInRegion(
#'      regionName_char = "chr22:18267969-18268249",
#'      arrayType = "450k"
#'    )
GetCpGsInRegion <- function(regionName_char,
                            arrayType = c("450k","EPIC")){


  arrayType <- match.arg(arrayType)

  switch(arrayType,
         "450k" = {
           CpGlocations_df = IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
         },
         "EPIC" = {
           CpGlocations_df = IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations
         }
  )
  CpGlocations.gr <- CpGlocations_df %>%
    makeGRangesFromDataFrame(end.field = "pos",
                             start.field = "pos",
                             ignore.strand = TRUE)

  ### Split the region name in chr and positions ###
  gr <- regionName_char %>%
    as.data.frame %>%
    separate(col = ".", into = c("seqnames","start","end")) %>%
    makeGRangesFromDataFrame()

  subsetByOverlaps(CpGlocations.gr, gr) %>% sort(ignore.strand = TRUE)  %>% names
}

#' Extract probe IDs for CpGs located in a genomic region
#'
#' @param regionName_char character string with location information for one region in
#'    this format: "chrxx:xxxxxx-xxxxxx"
#' @param arrayType Type of array, 450k or EPIC
#' @param nCores_int Number of computing cores to be used when executing code
#'    in parallel. Defaults to 1 (serial computing).
#' @param ... Dots for additional arguments passed to the cluster constructor.
#'    See \code{\link{CreateParallelWorkers}} for more information.
#'
#' @return vector of CpG probe IDs mapped to the genomic region
#' @export
#'
#' @importFrom tidyr separate %>%
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges subsetByOverlaps
#'
#' @examples
#'    GetCpGsInRegion(
#'      regionName_char = "chr22:18267969-18268249",
#'      arrayType = "450k"
#'    )
GetCpGsInAllRegion <- function(regionName_char,
                               arrayType = c("450k","EPIC"),
                               nCores_int = 1L,
                               ...){


  arrayType <- match.arg(arrayType)

  switch(arrayType,
         "450k" = {
           CpGlocations_df = IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
         },
         "EPIC" = {
           CpGlocations_df = IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations
         }
  )
  CpGlocations.gr <- CpGlocations_df %>%
    makeGRangesFromDataFrame(end.field = "pos",
                             start.field = "pos",
                             ignore.strand = TRUE)

  ### Split the region name in chr and positions ###
  message("Get regions from name")
  regions.df <- regionName_char %>%
    as.data.frame %>%
    separate(col = ".", into = c("seqnames","start","end"))
    regions.gr <- regions.df %>% makeGRangesFromDataFrame()
    regions.list <- regions.df %>% makeGRangesListFromDataFrame()
    names(regions.list) <- regionName_char

    CpGlocations.gr <- subsetByOverlaps(CpGlocations.gr, regions.gr)  %>%
      sort(ignore.strand = TRUE)

  message("Iterating over regions... this may take a while")

  cluster <- CreateParallelWorkers(nCores_int, ...)
  results <- bplapply(
    regions.list,
    function(x) {
      subsetByOverlaps(CpGlocations.gr, x) %>%
        names
    },  BPPARAM = cluster
  )
  results[regionName_char]
}


