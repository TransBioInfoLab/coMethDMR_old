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
#' UCSC RefGene Group mapped to the probe ("gene"), costumized annotation ("customized")
#' @param annotation.gr Annotation Granges used only if type is set to "customized",
#' first column of the metadata is used to annotate the probes
#' @param arrayType Type of array, 450k or EPIC
#' @param sort.by.or Sort bars by OR descresenting ? Default FALSE
#' @param tab.filename Table file name (csv file). Default no file will be output.
#' @param bar.colors Colors of the bars of each group
#' @param fisher.test.alternative  indicates the alternative hypothesis and must
#' be one of "two.sided" (default), "greater" or "less".
#' @importFrom ggplot2 ggsave position_dodge theme element_blank scale_fill_manual
#' @importFrom purrr reduce
#' @importFrom readr write_csv
#' @importFrom ggpubr ggbarplot
#' @importFrom reshape2 melt
#' @importFrom dplyr full_join %>%
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors values
#' @return A list with a ggplot2 barplo, a table and the list of contigency tables
#' @export
#' @examples
#' data(betasChr22_df)
#' result.list <- cpGsEnrichment(fg.probes = rownames(betasChr22_df)[1:100],
#'                               bg.probes = rownames(betasChr22_df)[-c(1:100)],
#'                               enrichment.type = "island",
#'                               arrayType = "450k",
#'                               save.plot = FALSE)
#'
#' ChmmModels <- readr::read_tsv("https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E073_15_coreMarks_segments.bed",col_names = FALSE)
#' colnames(ChmmModels) <- c("chr","start","end","state")
#' ChmmModels.gr <- makeGRangesFromDataFrame(ChmmModels,keep.extra.columns = TRUE)
#' result.list <- cpGsEnrichment(fg.probes = rownames(betasChr22_df)[1:100],
#'                               bg.probes = rownames(betasChr22_df)[-c(1:100)],
#'                               enrichment.type = "costum",
#'                               annotation.gr = ChmmModels.gr,
#'                               arrayType = "450k",
#'                               save.plot = FALSE)
cpGsEnrichment <- function (fg.probes,
                            bg.probes,
                            fg.label = "foreground",
                            bg.label = "background",
                            tab.filename,
                            fisher.test.alternative = "two.sided",
                            arrayType = c("450k","EPIC"),
                            save.plot = TRUE,
                            plot.filename = "barplot.pdf",
                            bar.colors,
                            sort.by.or = FALSE,
                            annotation.gr,
                            enrichment.type =  c("island","gene","customized")
){

    arrayType <- match.arg(arrayType)
    enrichment.type <- match.arg(enrichment.type)

    if (missing(fg.probes)) stop("Please, set fg.probes")
    if (missing(bg.probes)) stop("Please, set bg.probes")

    # be sure there is no overlap
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
    }  else if (enrichment.type == "customized") {
        probes.gr <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations %>%
            makeGRangesFromDataFrame(start.field = "pos",end.field = "pos")
        hits <- findOverlaps(probes.gr,annotation.gr,select = "first")
        col.name <- names(values(annotation.gr))[1]
        annot <- data.frame(row.names = names(probes.gr),
                            values(annotation.gr)[,1][hits],stringsAsFactors = FALSE)
        colnames(annot) <- col.name
    }


    fg.cts <- annot[unique(fg.probes),col.name]
    fg.cts <- plyr::count(fg.cts)
    colnames(fg.cts) <- c(col.name,"foreground")

    bg.cts <- annot[unique(bg.probes),col.name]
    bg.cts <- plyr::count(bg.cts)
    colnames(bg.cts) <- c(col.name,"background")

    cts <- purrr::reduce(.x = list(fg.cts, bg.cts), .f = (full_join))
    cts[is.na(cts)] <- 0
    cts.freq <- cts
    cts.freq[,2] <- 100 * cts.freq[,2] / sum(cts.freq[,2],na.rm = TRUE)
    cts.freq[,3] <- 100 * cts.freq[,3] / sum(cts.freq[,3],na.rm = TRUE)

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
                           if(any(is.na(m))) {
                               ft <- data.frame(p.value = NA, estimate = NA)
                           } else {
                               ft <- fisher.test(m,alternative = fisher.test.alternative)
                           }
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

    # Adding symbols and pvalues to table
    cuts <- rev(c(Inf, 0.05, 0.01, 0.001, 0.0001, -Inf))
    labs <- rev(c("ns", "*", "**", "***","****"))
    ret$p_value_symbol <- labs[findInterval(ret$p_value,cuts)]

    suppressMessages(
        df <- reshape2::melt(cts.freq)
    )
    df$variable <- factor(df$variable, levels = c('background', 'foreground'))


    if(missing(bar.colors)) {
        bar.colors <- c("background" = "#1F77B4",
                        "foreground" = "#FF7F0D")
    }

    if(sort.by.or){
        order <- unique(ret[order(-ret$odds_ratio),col.name])
    } else if (enrichment.type == "gene") {
        order <- c("TSS1500","TSS200","5'UTR","1stExon","Body","3'UTR","Intergenic")
        df <- df[plyr::alply(order,1,function(x) which(df[[col.name]] == x)) %>% unlist(use.names = F),]
    } else if (enrichment.type == "island") {
        order <- c("Shelf","Shore","Island","OpenSea")
        df <- df[plyr::alply(order,1,function(x) which(df[[col.name]] == x)) %>% unlist(use.names = F),]
    } else {
        # Alphabetically
        order <- sort(unique(df[[col.name]]))
    }

    plot <- ggpubr::ggbarplot(df,
                              y = "value",
                              x = col.name,
                              fill = "variable",
                              color = "white",
                              order = order,
                              palette = bar.colors,
                              ylab = "Frequency (% Counts)",
                              xlab = gsub("[^[:alnum:] ]"," ",col.name),
                              position = position_dodge(0.8)) +
        theme(legend.title = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1)) +
        scale_fill_manual(name = "",
                          labels = c(bg.label,
                                     fg.label),
                          values = bar.colors)

    # Addd pvalue significance
    y.pos <- df %>% dplyr::group_by_(col.name) %>% dplyr::summarise(n = max(value))
    y.pos <- y.pos[match(unique(ret[[col.name]]),y.pos[[col.name]]),]
    y.pos <- y.pos %>% pull(n)
    plot <- plot + geom_text(data = ret,
                             aes_(label = as.name("p_value_symbol"),
                                  x = as.name(col.name),
                                  y = y.pos),
                             position = position_dodge(width = 0.8),
                             vjust = -0.6)

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

#' @title Plot percentage of CpGs for diff. genomic features
#' @description creates barplot for percentage of dff. methylated CpGs
#' for enhancer, intergenic, intron, island, and other.
#' @param probes list of named probes
#' @param save.plot create a bar plot of frequncy ? Default TRUE.
#' @param plot.filename filename of the barplot
#' @param arrayType Type of array, 450k or EPIC
#' @param bar.colors Colors of the bars of each group
#' @param plot.title A string title to the plot
#' @param plot.width Plot width
#' @param plot.hight Plot height
#' @importFrom ggplot2 ggsave position_dodge theme element_blank scale_fill_manual element_text
#' @importFrom purrr reduce
#' @importFrom readr write_csv
#' @importFrom ggpubr ggbarplot
#' @importFrom reshape2 melt
#' @importFrom dplyr full_join %>%
#' @return A list with a ggplot2 barplo, a table and the list of contigency tables
#' @export
#' @examples
#' data(betasChr22_df)
#' probes.list <- list("Hypomethylated CpGs" = rownames(betasChr22_df)[1:100],
#'                     "Hypermethylated CpgS" = rownames(betasChr22_df)[100:200])
#' result.list <- cpGsGenomicFeatures(probes.list = probes.list,
#'                                arrayType = "450k",
#'                                     bar.colors = c("#1F77B4", "#FF7F0D"),
#'                                    save.plot = FALSE)
#' file <- paste0("https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations",
#'                "/ChmmModels/coreMarks/jointModel/final/E073_15_coreMarks_segments.bed")
#' ChmmModels <- readr::read_tsv(file,col_names = FALSE,col_types = readr::cols())
#' colnames(ChmmModels) <- c("chr","start","end","state")
#' states <- readr::read_csv("../../DNAm_RNA/data/chromHMM_labels.csv",col_names = FALSE)
#' states$X1 <- paste0("E",states$X1)
#' ChmmModels$state <- states$X3[match(ChmmModels$state,states$X1)]
#' ChmmModels.gr <- makeGRangesFromDataFrame(ChmmModels,keep.extra.columns = TRUE)
#' regions.and.single.plot <-
#' result.list <- cpGsGenomicFeatures(
#'   probes.list,
#'    annotation.gr = ChmmModels.gr,
#'    plot.width = 12,
#'    enrichment.type = "customized",
#'    save.plot = FALSE
#')
#'
cpGsGenomicFeatures <- function (probes.list,
                                 arrayType = c("450k","EPIC"),
                                 save.plot = TRUE,
                                 bar.colors,
                                 plot.width = 7,
                                 plot.height = 7,
                                 plot.title = NULL,
                                 annotation.gr = NULL,
                                 enrichment.type =  c("island_gene","customized"),
                                 plot.filename = "barplot.pdf"
){
    arrayType <- match.arg(arrayType)
    enrichment.type <- match.arg(enrichment.type)

    if (missing(probes.list)) stop("Please, set fg.probes")

    if(enrichment.type == "island_gene"){
        if (arrayType == "450k") {
            annot <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Islands.UCSC
        } else {
            annot <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Islands.UCSC
        }
        annot$Relation_to_Island <-  gsub("N_|S_","",annot$Relation_to_Island)

        if (arrayType == "450k") {
            annot.other <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other
        } else {
            annot.other <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other
        }

        annot.other$UCSC_RefGene_Group_hierarchy <- annot.other$UCSC_RefGene_Group
        annot.other$UCSC_RefGene_Group_hierarchy[grep("TSS200",annot.other$UCSC_RefGene_Group_hierarchy)] <- "TSS200"
        annot.other$UCSC_RefGene_Group_hierarchy[grep("TSS1500",annot.other$UCSC_RefGene_Group_hierarchy)] <- "TSS1500"
        annot.other$UCSC_RefGene_Group_hierarchy[grep("5'UTR",annot.other$UCSC_RefGene_Group_hierarchy)] <- "5'UTR"
        annot.other$UCSC_RefGene_Group_hierarchy[grep("1stExon",annot.other$UCSC_RefGene_Group_hierarchy)] <- "1stExon"
        annot.other$UCSC_RefGene_Group_hierarchy[grep("Body",annot.other$UCSC_RefGene_Group_hierarchy)] <- "Body"
        annot.other$UCSC_RefGene_Group_hierarchy[grep("3'UTR",annot.other$UCSC_RefGene_Group_hierarchy)] <- "3'UTR"
        annot.other$UCSC_RefGene_Group_hierarchy[annot.other$UCSC_RefGene_Group_hierarchy == ""] <- "Intergenic"
        annot.other$UCSC_RefGene_Group_hierarchy[is.na(annot.other$UCSC_RefGene_Group_hierarchy)] <- "Intergenic"
        if(!is(probes.list,"list")) {
            probes.list <- list("Probes" = probes.list)
        }
        cts.annot <- plyr::ldply(probes.list, function(x){
            cts <- annot.other[x, "UCSC_RefGene_Group_hierarchy"]
            cts <- plyr::count(cts);
            cts2 <- annot[x, "Relation_to_Island"]
            cts2 <- plyr::count(cts2);
            cts <- rbind(cts,
                         cts2
                         #data.frame("x" = "Enhancer",
                         #           freq = sum(x %in% enhancer.probes)
                         #)
            )
            cts$counts <- cts$freq
            cts$freq <- (100 * cts$freq)/length(x)
            return(cts)
        })

        order <-  c(
            "TSS1500",
            "TSS200",
            "5'UTR",
            "1stExon",
            "Body",
            "3'UTR",
            "Intergenic",
            "Shelf",
            "Shore",
            "Island",
            "OpenSea"
        )
    }  else if (enrichment.type == "customized") {
        if(is.null(annotation.gr)){
            stop("Set annotation.gr for customized option")
        }
        probes.gr <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations %>%
            makeGRangesFromDataFrame(start.field = "pos",end.field = "pos")
        hits <- findOverlaps(probes.gr,annotation.gr,select = "first")
        col.name <- names(values(annotation.gr))[1]
        annot <- data.frame(row.names = names(probes.gr),
                            values(annotation.gr)[,1][hits],stringsAsFactors = FALSE)
        colnames(annot) <- col.name

        if(!is(probes.list,"list")) {
            probes.list <- list("Probes" = probes.list)
        }
        cts.annot <- plyr::ldply(probes.list, function(x){
            cts <- annot[x,]
            cts <- plyr::count(cts);
            cts$counts <- cts$freq
            cts$freq <- (100 * cts$freq)/length(x)
            missing <- annot %>% dplyr::pull() %>% unique %>% setdiff(cts$x)
            if(length(missing)) cts <- rbind(cts,data.frame("x" = missing, "freq" = 0, "counts" = 0))
            return(cts)
        })
        # order <- cts.annot$x[order(cts.annot$counts)] %>% unique
        order <- sort(unique(annot[[col.name]]))

    }



    if(missing(bar.colors)) {
        bar.colors <- c("#999999",
                        "#56B4E9", "#E69F00",
                        "#009E73",
                        "#F0E442", "#0072B2",
                        "#D55E00", "#CC79A7")
    }


    plot <- ggpubr::ggbarplot(cts.annot,
                              y = "freq",
                              x = "x",
                              fill = ".id",
                              color = ".id",
                              order = order,
                              palette = bar.colors,
                              ylab = "Percentage CpGs (%)",
                              xlab = "Genomic Feature",
                              position = position_dodge(0.8)) +
        theme(legend.title = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1))
    if(!is.null(plot.title)){
        plot <- plot + ggtitle(plot.title)
    }
    if(save.plot){
        ggsave(plot.filename, plot = plot, width = plot.width,height = plot.height)
    }

    return(list("table" = cts.annot,
                "plot" = plot))
}


