#' Annotate \code{coMethDMR} Pipeline Results
#'
#' @description Given a data frame with regions in the genome, add gene symbols,
#'   UCSC reference gene accession, UCSC reference gene group and
#'   relation to CpG island.
#'
#' @param lmmRes_df A data frame returned by \code{\link{lmmTestAllRegions}}.
#' This data frame must contain the following columns:
#'    \itemize{
#'      \item{\code{chrom} : }{the chromosome the region is on, e.g. ``chr22''}
#'      \item{\code{start} : }{the region start point}
#'      \item{\code{end} : }{the region end point}
#'      }
#'
#' Optionally, the data frame can also has \code{regionType},
#' which is a character string marking the type of genomic region tested.
#'   See details below.
#'
#' @param arrayType Type of array: 450k or EPIC
#'
#' @param nCores_int Number of computing cores to be used when executing code
#'    in parallel. Defaults to 1 (serial computing).
#' @param ... Dots for additional arguments passed to the cluster constructor.
#'    See \code{\link{CreateParallelWorkers}} for more information.
#'
#' @return A data frame with
#'    \itemize{
#'      \item the location of the genomic region's chromosome (\code{chrom}),
#'        start (\code{start}), and end (\code{end});
#'    \item UCSC annotation information (\code{UCSC_RefGene_Group},
#'        \code{UCSC_RefGene_Accession}, and \code{UCSC_RefGene_Name}); and
#'      \item a list of all of the probes in that region (\code{probes}).
#'    }
#'
#' @details The region types include \code{"NSHORE"}, \code{"NSHELF"},
#'    \code{"SSHORE"}, \code{"SSHELF"}, \code{"TSS1500"}, \code{"TSS200"},
#'    \code{"UTR5"}, \code{"EXON1"}, \code{"GENEBODY"}, \code{"UTR3"}, and
#'    \code{"ISLAND"}.
#'
#' @importFrom dplyr bind_rows bind_cols
#' @export
#'
#' @examples
#'    lmmResults_df <- data.frame(
#'      chrom = c("chr22", "chr22", "chr22", "chr22", "chr22"),
#'      start = c("39377790", "50987294", "19746156", "42470063", "43817258"),
#'      end   = c("39377930", "50987527", "19746368", "42470223", "43817384"),
#'      regionType = c("TSS1500", "EXON1", "ISLAND", "TSS200", "ISLAND"),
#'      stringsAsFactors = FALSE
#'    )
#'
#'    AnnotateResults(
#'      lmmRes_df = lmmResults_df,
#'      arrayType = "450k"
#'    )
#'
AnnotateResults <- function(lmmRes_df,
                             arrayType = c("450k","EPIC"),
                             nCores_int = 1L,
                             ...){
    ###  Check Inputs  ###
    stopifnot(
        "data.frame" %in% class(lmmRes_df),
        all(c("chrom", "start", "end") %in% colnames(lmmRes_df))
    )
    arrayType <- match.arg(arrayType)
    checkAnnotationPkg(arrayType)

    ###  Pull Database  ###
    switch(
        arrayType,
        "450k" = {
            UCSCinfo_df  <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other
            IslandsUCSCinfo_df <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Islands.UCSC

        },
        "EPIC" = {
            UCSCinfo_df  <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other
            IslandsUCSCinfo_df <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Islands.UCSC

        }
    )

    # UCSC Gene Info
    UCSCinfo_df <- as.data.frame(UCSCinfo_df)
    interestingColumns_char <- c(
        "UCSC_RefGene_Name",
        "UCSC_RefGene_Accession",
        "UCSC_RefGene_Group"
    )
    UCSCinfo_df <- UCSCinfo_df[, interestingColumns_char]

    # UCSC Island Info
    IslandsUCSCinfo_df <- as.data.frame(IslandsUCSCinfo_df)

    probes.list <- lmmRes_df %>%
        tidyr::unite("region",c("chrom","start","end")) %>%
        dplyr::pull(region) %>%
        GetCpGsInAllRegion

    # memory management
    UCSCinfo_df <- UCSCinfo_df[rownames(UCSCinfo_df) %in% unlist(probes.list),]
    IslandsUCSCinfo_df <- IslandsUCSCinfo_df[rownames(IslandsUCSCinfo_df) %in% unlist(probes.list),]

    ###  Define Wrapper Function  ###
    AnnotateRow <- function(probes_char, info_df, island_df, includeType){

        # Find UCSC Annotation Information for those Probes
        infoOut_df <- info_df[probes_char, ]

        # Find UCSC Relation to Island Information for those Probes
        islandOut_df <- island_df[probes_char, ]

        ###  Wrangle UCSC Annotation  ###
        refGeneGroup_char <-
            sort(unique(unlist(
                strsplit(infoOut_df$UCSC_RefGene_Group, ";")
            )))

        refGeneGroup_char <- sort(unique(refGeneGroup_char))

        refGeneAcc_char <-
            sort(unique(unlist(
                strsplit(infoOut_df$UCSC_RefGene_Accession, ";")
            )))

        refGeneName_char <-
            sort(unique(unlist(
                strsplit(infoOut_df$UCSC_RefGene_Name, ";")
            )))

        refIslandRelation_char <-
            sort(unique(islandOut_df$Relation_to_Island))


        ###  Return Annotated 1-Row Data Frame  ###
        if(includeType){
            # row_df$Relation_to_UCSC_CpG_Island <- ifelse(
            #   test = row_df$regionType %in%
            #     c("NSHELF", "NSHORE", "ISLAND", "SSHORE", "SSHELF"),
            #   yes  = row_df$regionType,
            #   no   = ""
            # )
        }
        tibble::tibble("UCSC_RefGene_Group" = paste0(refGeneGroup_char, collapse = ";"),
                       "UCSC_RefGene_Accession" = paste0(refGeneAcc_char, collapse = ";"),
                       "UCSC_RefGene_Name" = paste0(refGeneName_char, collapse = ";"),
                       "Relation_to_Island" = paste0(refIslandRelation_char, collapse = ";")
        )
    }

    inclType_logi <- !is.null(lmmRes_df$regionType)

    cluster <- CreateParallelWorkers(nCores_int, ...)

    resultsAnno_ls <- bplapply(
        probes.list,
        function(probes){
            AnnotateRow(
                probes_char = probes,
                info_df = UCSCinfo_df,
                island_df = IslandsUCSCinfo_df,
                includeType = inclType_logi
            )},  BPPARAM = cluster
    ) %>% dplyr::bind_rows()
    dplyr::bind_cols(lmmRes_df,resultsAnno_ls)
}
