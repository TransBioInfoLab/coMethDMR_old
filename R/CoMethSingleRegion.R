#' Wrapper function to find contiguous and comethyalted sub-regions within a
#'   pre-defined genomic region
#'
#' @param CpGs_char vector of CpGs in the inputting pre-defined genomic region.
#' @param dnam matrix (or data frame) of beta values, with row names = CpG ids,
#'    column names = sample ids. This should include the CpGs in \code{CpGs_char},
#'    as well as additional CpGs.
#' @param betaToM indicates if converting methylation beta values mvalues
#' @param method method for computing correlation, can be  "pearson" or "spearman"
#' @param rDropThresh_num thershold for min correlation between a cpg with sum
#'    of the rest of the CpGs
#'
#' @param minPairwiseCorr Calculate minimun Pairs wise correlation and use as a filter ?
#' Default: NULL. To use this filter set a number between 0 and 1.
#'
#' @param minCpGs mininum number of CpGs to be considered a "region".
#'    Only regions with more than \code{minCpGs} will be returned.
#' @param arrayType Type of array, can be "450k" or "EPIC"

#' @param returnAllCpGs When there is not a contiguous comethylated region in
#'    the inputing pre-defined region, \code{returnAllCpGs = 1} indicates
#'    outputting all the CpGs in the input region, while \code{returnAllCpGs = 0}
#'    indicates not returning any CpG.
#'
#' @return A list with two components:
#'   \itemize{
#'     \item{\code{Contiguous_Regions} : }{a data frame with \code{CpG} (CpG ID),
#'       \code{Chr} (chromosome number), \code{MAPINFO} (genomic position),
#'       \code{r_drop} (correlation between the CpG with rest of the CpGs),
#'       \code{keep} (indicator for co-methylated CpG), \code{keep_contiguous}
#'       (index for contiguous comethylated subregion)
#'     }
#'     \item{\code{CpGs_subregions} : }{lists of CpGs in each contiguous
#'       co-methylated subregion
#'     }
#'   }
#'
#' @export
#'
#' @examples
#'    data(betasChr22_df)
#'
#'    CpGsChr22_char <- c(
#'      "cg02953382", "cg12419862", "cg24565820", "cg04234412", "cg04824771",
#'      "cg09033563", "cg10150615", "cg18538332", "cg20007245", "cg23131131",
#'      "cg25703541"
#'    )
#'    CoMethSingleRegion(
#'      CpGs_char = CpGsChr22_char,
#'      dnam = betasChr22_df
#'    )
#'
#'    data(betaMatrix_ex3)
#'    CpGsEx3_char <- c(
#'      "cg14221598", "cg02433884", "cg07372974", "cg13419809", "cg26856676",
#'      "cg25246745"
#'    )
#'    CoMethSingleRegion(
#'      CpGs_char = CpGsEx3_char,
#'      dnam = t(betaMatrix_ex3),
#'      returnAllCpGs = TRUE
#'    )
#'
CoMethSingleRegion <- function(CpGs_char,
                               dnam,
                               betaToM = TRUE,
                               rDropThresh_num = 0.4,
                               minPairwiseCorr = NULL,
                               method = c("pearson", "spearman"),
                               minCpGs = 3,
                               arrayType = c("450k","EPIC"),
                               returnAllCpGs = FALSE){

    arrayType <- match.arg(arrayType)
    method <- match.arg(method)

    ### Order CpGs by genomic location ###
    CpGsOrdered_df <- OrderCpGsByLocation(
        CpGs_char, arrayType, output = "dataframe"
    )

    ### Extract beta matrix for the input CpGs ###
    # take common cpgs in beta matrix and the region first
    commonCpGs_char <- intersect (CpGsOrdered_df$cpg, row.names(dnam))

    if (length(commonCpGs_char) >= minCpGs){

        betaCluster_mat <- dnam[commonCpGs_char, ]

        ### Transpose beta matrix ###
        betaClusterTransp_mat <- t(betaCluster_mat)

        ### Mark comethylated CpGs ###
        keepCpGs_df <- MarkComethylatedCpGs(
            betaCluster_mat = betaClusterTransp_mat,
            method = method,
            betaToM = betaToM,
            rDropThresh_num = rDropThresh_num
        )
        # If MarkComethylatedCpGs was not able to execute
        # (either only NAs or not beta values, function cannot be continued)
        if(is.null(keepCpGs_df)) return(NULL)

        ### Find contiguous comethylated regions ###
        keepContiguousCpGs_df <- FindComethylatedRegions(
            CpGs_df = keepCpGs_df
        )

        ### Split CpG dataframe by Subregion ###
        keepContiguousCpGs_ls <- SplitCpGDFbyRegion(
            keepContiguousCpGs_df, arrayType, returnAllCpGs
        )


        ### calculate minPairwiseCor and subset for each subregion
        ### and add information to df

        if(!is.null(minPairwiseCorr)){
            pairwiseResults <- minPairwiseCor(
                betaCluster_mat = betaClusterTransp_mat,
                method = method,
                betaToM = betaToM,
                minPairwiseCorr = minPairwiseCorr,
                probes.list = keepContiguousCpGs_ls
            )
            keepminPairwiseCor_df <- pairwiseResults$keepminPairwiseCor_df

            # if returnAllCpGs no filtering
            if(!returnAllCpGs){
                keepContiguousCpGs_ls <- pairwiseResults$probes.list.filtered
            }
        } else {
            keepminPairwiseCor_df <- NULL
        }

        ### Create Output Data Frame  ###
        coMethCpGs_df <- CreateOutputDF(
            keepCpGs_df = keepCpGs_df,
            keepContiguousCpGs_df = keepContiguousCpGs_df,
            CpGsOrdered_df =  CpGsOrdered_df,
            returnAllCpGs =  returnAllCpGs,
            keepminPairwiseCor_df = keepminPairwiseCor_df
        )

        ### Create input region & number of co-methylated regions (Lanyu) ###
        if (returnAllCpGs == FALSE & all(keepContiguousCpGs_df$Subregion == 0)){
            inputRegion = NULL
            nCoMethRegions = NULL
        } else {

            input <- as.character(coMethCpGs_df$Region[1])
            nCoMeth <- as.numeric(max(coMethCpGs_df$keep_contiguous))

            if (nCoMeth >= 1){
                inputRegion <- rep(input, times = nCoMeth)
                nCoMethRegions <- rep(nCoMeth, times = nCoMeth)
            } else {
                inputRegion <- input
                nCoMethRegions <- nCoMeth
            }

        }

        ### Create output list of data frame and CpGs by subregion ###
        coMethCpGs_ls <- list(
            contiguousRegions = coMethCpGs_df,
            CpGsSubregions = keepContiguousCpGs_ls,
            inputRegion = inputRegion,
            nCoMethRegions = nCoMethRegions
        )

        coMethCpGs_ls

    } else {
        return(NULL)
    }



}
