#' make.map.table takes the output from map function (mapping reads onto contigs) and 
#'
#' takes the table which generates the overlap output back to data.frame
#'
#' @param overlaps output from the findOverlaps function from IRanges
#'
#' @return a data.frame object with the names of contigs and the number of reads mapped to them
#'
#' @export
make.map.table <- function(contigs, reads, datype="cDNA"){
    overlapDF = findOverlaps(reads, contigs)
    if(length(overlapDF) > 0){
        countDF            = overlapDF@subjectHits %>% table %>% data.frame %>% setNames(c("contigID", "Freq"))
        contigNames = as.character(seqnames(contigs)@values)[countDF$contigID]
        countDF$contigName = contigNames
        mcolDF = contigs %>% mcols %>% as.data.frame
        cbind(
        countDF,
        mcolDF[countDF$contigID, c("type", "trueLength")]
        ) %>% mutate(type=datype)
    }else{
        message("There were no overlaps")
        data.frame(contigID=factor(), Freq = integer(), contigName = character(), type = character(), trueLength = integer())
    }
}
