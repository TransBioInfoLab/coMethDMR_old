checkAnnotationPkg <- function(arrayType){
    if(arrayType == "450k"){
        if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
            stop("Package \"IlluminaHumanMethylation450kanno.ilmn12.hg19\" needed for this function to work. Please install it.",
                 call. = FALSE)
        }
    } else {

        if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b2.hg19", quietly = TRUE)) {
            stop("Package \"IlluminaHumanMethylationEPICanno.ilm10b2.hg19\" needed for this function to work. Please install it.",
                 call. = FALSE)
        }
    }
}
