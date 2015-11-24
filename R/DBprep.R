#' DBprep prepares input for NEO4J database
#' 
#' DBprep takes the output from the mapBlat pipeline , which maps reads to the max diversity region 
#' and prepares it as node, rel directories 
#' for input by the omics package, 
#' 
#' @param dir the directory to store the files
#' @param allKODF data.frameObj which holds all count data
#' 
#' @export
DBprep <- function(theDIR, allKODF){
    allKODF$ko = stringr::str_replace(allKODF$ko, "^(ko:)*", "ko:")
    allKODF$contigID = paste(allKODF$ko, allKODF$contigName, sep = ":")
    nodes <- allKODF                                                                                                         %>%
        mutate(label = 'contigs')                                                                                            %>%
        select(contigID, Freq_cDNA, rpkm_cDNA, Freq_gDNA, rpkm_gDNA, label)                                                  %>%
        setNames(c("contig:string:contigid", "cDNAFreq:long","cDNAFPKM:double","gDNAFreq:long","gDNAFPKM:double","l:label"))

    rels <- allKODF %>%
        mutate(type = 'abundance') %>%
        select(contigID, ko, type) %>%
        setNames(c("contig:string:contigid", "ko:string:koid", "type:string"))

    c("%s", "%s/nodes/", "%s/rels/") %>% sapply(function(path) sprintf(path, theDIR) %>% dir.create) 

#writing
write.table(nodes, file =sprintf("%s/nodes/contigs.node", theDIR), sep="\t", quote=F, row.names=F)
write.table(rels, file=sprintf("%s/rels/contig2ko.rel", theDIR), sep="\t", quote=F, row.names=F)
}
