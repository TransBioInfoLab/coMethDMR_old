#' Create Output Dataframe
#'
#' @param keepCpGs_df a data frame with CpG = CpG name,
#'    keep = indicator for co-methylated CpG,
#'    r_drop = correlation between the CpG with rest of the CpGs
#' @param keepContiguousCpGs_df a data frame with ProbeID = CpG name,
#'    Subregion = cotiguous comethylated subregion number
#' @param CpGsOrdered_df a data frame of CpG location with
#'    chr = chromosome number, pos = genomic position, cpg = CpG name
#' @param keepminPairwiseCor_df a data frame with
#'    Subregion = cotiguous comethylated subregion number
#' @param returnAllCpGs indicates if outputting all the CpGs in the region
#'    when there is not a contiguous comethylated region or
#'    only the CpGs in the contiguous comethylated regions
#'
#' @return a data frame with CpG = CpG name,
#'    Chr = chromosome number,
#'    MAPINFO = genomic position,
#'    r_drop = correlation between the CpG with rest of the CpGs,
#'    keep = indicator for co-methylated CpG,
#'    keep_contiguous = cotiguous comethylated subregion number
#'
#' @keywords internal
#' @importFrom dplyr left_join
#' @export
#'
#' @examples
#'    data(betasChr22_df)
#'    CpGsChr22_char<-c("cg02953382", "cg12419862", "cg24565820", "cg04234412",
#'       "cg04824771", "cg09033563", "cg10150615", "cg18538332", "cg20007245",
#'       "cg23131131", "cg25703541")
#'    CpGsOrdered_df <- OrderCpGsByLocation(
#'       CpGsChr22_char, arrayType="450k", output = "dataframe"
#'    )
#'    betaCluster_mat <- betasChr22_df[CpGsOrdered_df$cpg,]
#'    betaClusterTransp_mat <- t(betaCluster_mat)
#'    keepCpGs_df <- MarkComethylatedCpGs(betaCluster_mat = betaClusterTransp_mat)
#'    keepContiguousCpGs_df <- FindComethylatedRegions(CpGs_df = keepCpGs_df)
#'    pairwiseCor <- minPairwiseCor(betaClusterTransp_mat,probes.list = list(CpGsChr22_char))
#'    keepminPairwiseCor_df <- pairwiseCor$keepminPairwiseCor.df
#'    CreateOutputDF(keepCpGs_df, keepContiguousCpGs_df, CpGsOrdered_df,keepminPairwiseCor_df)
#'    CreateOutputDF(keepCpGs_df, keepContiguousCpGs_df, CpGsOrdered_df,NULL)
CreateOutputDF <- function(keepCpGs_df,
                           keepContiguousCpGs_df,
                           CpGsOrdered_df,
                           keepminPairwiseCor_df,
                           returnAllCpGs = FALSE){

    if (returnAllCpGs == FALSE & all(keepContiguousCpGs_df$Subregion == 0)){
        NULL
    } else {

        output_df <- merge(
            keepCpGs_df,
            keepContiguousCpGs_df,
            by.x = "CpG",
            by.y = "ProbeID",
            all.x = TRUE
        )
        output_df <- merge(CpGsOrdered_df, output_df, by.x = "cpg", by.y = "CpG")

        print(output_df)
        if (!is.null(keepminPairwiseCor_df)) {
            output_df <- left_join(output_df, keepminPairwiseCor_df, by = "Subregion")
        }

        output_df <-
            output_df[order(output_df$chr, output_df$pos),]
        output_df [is.na(output_df)] <- 0

        coMethCpGs_df <- data.frame(
            Region = NameRegion(CpGsOrdered_df),
            CpG = output_df$cpg,
            Chr = output_df$chr,
            MAPINFO = output_df$pos,
            r_drop = output_df$r_drop,
            keep = output_df$keep,
            keep_contiguous = output_df$Subregion
        )

        if (!is.null(keepminPairwiseCor_df)) {
            coMethCpGs_df <- cbind(
                coMethCpGs_df,
                data.frame(
                    regionMinPairwiseCor = output_df$minPairwiseCor,
                    keep_regionMinPairwiseCor = output_df$keepminPairwiseCor
                )
            )
        }
        coMethCpGs_df
    }
}
