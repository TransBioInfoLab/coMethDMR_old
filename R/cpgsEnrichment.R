#' @title Calculate probes enrichment
#' @description This function calculates the enrichment of cpgs compared to a background
# using a fisher test and Illumina annotation (relation to island and UCSC RefGene Group).
# It outputs a barplot of frequency of each category within foreground and background,
# a table with the number of counts, frequency, fisher test p-value and OR,
# values in the contigente table
#  a = probes in foregound within the category (i.e. island)
#  b = probes in foregound not within the cateogory (i.e. not island [opensea, shelf, shore])
#  c = probes in backgound within the category (i.e. island)
#  d = probes in backgound not within the cateogory (i.e. not island [opensea, shelf, shore])
#' @param fg.probes foreground probes (i.e. hypo methylated pobes)
#' @param bg.probes backgound probes (i.e. all comethylated probes)
#' @param save.plot create a bar plot of frequncy for foreground and background ? Default TRUE.
#' @param plot.filename filename of the barplot
#' @param fg.label Foreground label
#' @param bg.label Background label
#' @param type of enrichment to caculate: relation to island ("island") or
#' UCSC RefGene Group mapped to the probe ("gene")
#' @param arrayType Type of array, 450k or EPIC
#' @param tab.filename Table file name (csv file). Default no file will be output.
#' @importFrom ggplot2 ggsave position_dodge theme element_blank scale_fill_manual
#' @importFrom purrr reduce
#' @importFrom readr write_csv
#' @importFrom ggpubr ggbarplot
#' @importFrom reshape2 melt
#' @importFrom dplyr full_join %>%
#' @return A list with a ggplot2 barplo, a table and the list of contigency tables
#' @examples
#' data(betasChr22_df)
#' result.list <- cpGsEnrichment(fg.probes = rownames(betasChr22_df)[1:100],
#'                               bg.probes = rownames(betasChr22_df)[-c(1:100)],
#'                               enrichment.type = "island",
#'                               arrayType = "450k",
#'                               save.plot = FALSE)
cpGsEnrichment <- function (fg.probes,
                            bg.probes,
                            fg.label = "foreground",
                            bg.label = "background",
                            tab.filename,
                            arrayType = c("450k","EPIC"),
                            save.plot = TRUE,
                            plot.filename = "barplot.pdf",
                            enrichment.type =  c("island","gene")
){

    arrayType <- match.arg(arrayType)
    enrichment.type <- match.arg(enrichment.type)

    if (missing(fg.probes)) stop("Please, set fg.probes")
    if (missing(bg.probes)) stop("Please, set bg.probes")

    # be sure there is no overlap
    fg.probes <- setdiff(fg.probes, bg.probes)
    bg.probes <- setdiff(bg.probes, fg.probes)

    if (enrichment.type == "island") {

        col.name <- "Relation_to_Island"

        if (arrayType == "450k") {
            annot <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Islands.UCSC
        } else {
            annot <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Islands.UCSC
        }

        annot$Relation_to_Island <-  gsub("N_|S_","",annot$Relation_to_Island)

    } else if (enrichment.type == "gene") {
        col.name <- "UCSC_RefGene_Group_hierarchy"

        if (arrayType == "450k") {
            annot <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other
        } else {
            annot <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other
        }

        annot$UCSC_RefGene_Group_hierarchy <- annot$UCSC_RefGene_Group
        annot$UCSC_RefGene_Group_hierarchy[grep("TSS200",annot$UCSC_RefGene_Group_hierarchy)] <- "TSS200"
        annot$UCSC_RefGene_Group_hierarchy[grep("TSS1500",annot$UCSC_RefGene_Group_hierarchy)] <- "TSS1500"
        annot$UCSC_RefGene_Group_hierarchy[grep("5'UTR",annot$UCSC_RefGene_Group_hierarchy)] <- "5'UTR"
        annot$UCSC_RefGene_Group_hierarchy[grep("1stExon",annot$UCSC_RefGene_Group_hierarchy)] <- "1stExon"
        annot$UCSC_RefGene_Group_hierarchy[grep("Body",annot$UCSC_RefGene_Group_hierarchy)] <- "Body"
        annot$UCSC_RefGene_Group_hierarchy[grep("3'UTR",annot$UCSC_RefGene_Group_hierarchy)] <- "3'UTR"
        annot$UCSC_RefGene_Group_hierarchy[annot$UCSC_RefGene_Group_hierarchy == ""] <- "Intergenic"
        annot$UCSC_RefGene_Group_hierarchy[is.na(annot$UCSC_RefGene_Group_hierarchy)] <- "Intergenic"
    }



    fg.cts <- annot[unique(fg.probes),col.name]
    fg.cts <- plyr::count(fg.cts)
    colnames(fg.cts) <- c(col.name,"foreground")

    bg.cts <- annot[unique(bg.probes),col.name]
    bg.cts <- plyr::count(bg.cts)
    colnames(bg.cts) <- c(col.name,"background")

    cts <- purrr::reduce(.x = list(fg.cts, bg.cts), .f = (full_join))
    cts.freq <- cts
    cts.freq[,2] <- 100 * cts.freq[,2] / sum(cts.freq[,2])
    cts.freq[,3] <- 100 * cts.freq[,3] / sum(cts.freq[,3])

    m.list <- plyr::alply(.data = 1:nrow(cts),
                          .margins = 1,
                          .fun = function(idx){
                              a <- cts[idx,2]
                              b <- cts[-idx,2] %>% sum
                              c <- cts[idx,3] # - a (only if overlapping set)
                              d <- (cts[-idx,3] %>% sum) # - b (only if overlapping set)
                              m <- matrix(c(a,c,b,d),
                                          nrow = 2,
                                          dimnames = list(type = c("Yes","No"), dir = c("Yes","No"))
                              )

                              names(dimnames(m)) <- c(fg.label, cts[idx,1] %>% as.character())
                              return(m)
                          })
    names(m.list) <- cts[,1]  %>% as.character()
    attr(m.list,"split_labels") <- NULL
    attr(m.list,"split_type") <- NULL

    ret <- plyr::ldply(.data = m.list,
                       .fun = function(m){
                           ft <- fisher.test(m,alternative = 'greater')

                           df <- data.frame(
                               "p_value" = ft$p.value,
                               "odds_ratio" = ft$estimate,
                               "a" = m[1,1],
                               "b" = m[1,2],
                               "c" = m[2,1],
                               "d" = m[2,2]
                           )
                           return(df)
                       },.id = NULL)
    ret$X1 <- NULL
    ret[[col.name]] <- cts[,1]  %>% as.character()

    suppressMessages(
        df <- reshape2::melt(cts.freq)
    )
    df$variable <- factor(df$variable, levels = c('background', 'foreground'))
    plot <- ggpubr::ggbarplot(df,
                              y = "value",
                              x = col.name,
                              fill = "variable",
                              color = "white",
                              ylab = "Frequency (% Counts)",
                              xlab = ifelse(col.name == "Relation_to_Island", "Relation to island", "UCSC RefGene Group"),
                              position = position_dodge(0.8)) +
        theme(legend.title = element_blank()) +
        scale_fill_manual(name = "",
                          labels = c(bg.label,
                                     fg.label),
                          values = c("background" = "#1F77B4",
                                     "foreground" = "#FF7F0D"))

    colnames(cts)[2:3] <- paste0("cts_",colnames(cts)[2:3])
    colnames(cts.freq)[2:3] <- paste0("freq_",colnames(cts.freq)[2:3])

    suppressWarnings(
        suppressMessages(
            tab <- purrr::reduce(.x = list(cts, cts.freq, ret), .f = (full_join))
        )
    )

    if(!missing(tab.filename)) {
        message("Saving file as: ", tab.filename)
        readr::write_csv(tab,path = tab.filename)
        return(NULL)
    }

    if(save.plot){
        ggsave(plot.filename, plot = plot)
    }

    return(list("table" =  tab,
                "plot" = plot,
                "m.list" = m.list))
}

