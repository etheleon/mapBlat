#'
#' contigStats integrates original assembly data with pAss max diversity analyses
#' 
#'
#' counts number of contigs formed per KO and annotate which of those are captured 
#' within the max diversity region
#'
#' @param contigsDIR folder storing newbler output wrt contigs
#' @param pAssOutputDIR path to directory
#' @param targetFile logFile
#' @param mapping to include mapping information or not to the Max Diversity Region

contigStats <- function(contigsDIR, pAssOutputDIR, pass=F)
{
    tempFile = "mapBlat_tempFile"

    #path to perl Script for processing
    buildDB     = find.package('mapBlat') %>% file.path('perl') %>% file.path('buildDB.pl')
    buildDB.md  = find.package('mapBlat') %>% file.path('perl') %>% file.path('buildDB_md.pl')

    #basic stats of the contigs generated from assembler:newbler
    cmd <- paste('perl',
                     shQuote(buildDB),
                     shQuote(contigsDIR),
                     shQuote(tempFile),
                      sep=" ")
    try(system(cmd))
    allContigs = read.table(targetFile, h=T)
    file.remove(tempFile)

    #basic stats of the max diversity region identified after the MSA
    cmd <- paste('perl',
                     shQuote(buildDB.md),
                     shQuote(pAssOutputDIR),
                     shQuote(targetFile),
                      sep=" ")
    try(system(cmd))
    mdContigs = read.table(targetFile, h=T)
    file.remove(tempFile)

    mdContigs$mx    = 1
    merge(allContigs, mdContigs, by=c("ko", "contig"), all=TRUE)
}
