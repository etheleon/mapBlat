#!/usr/bin/env Rscript
#' ## Introduction

#' We bin reads into functional bins / ortholog groups, KO and assembled them to give KO specific contigs (using Newbler).
#' Assuming each contig is a gene, gene numbers in the max diverse region of their MSA acts as proxy to intra-ortholog diversity
#' within the metagenomic sample
#'
#' In this analysis, we map mRNA reads onto genes in the max diversity regions.
#'
#' ### Identifying contigs within the MAX DIVERSITY window

#' We partition the contigs found in the max diversity regions in the following manner:
#' ![maxDiversity](maxDiversity.jpg)

args <- commandArgs(T) #in this case you should do `ls out/assm.0700.chosenContigWindwsfna.200/*fna | wc -l`

suppressPackageStartupMessages(library(mapBlat))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(water)); mc()

#ko <- args[1] %||% 'K00001'

`%||%` <- function(x, y){
      if (is.na(x)) y else x
}

start         <- args[1] %||% 1
end           <- args[2] %||% 1
finalOutput   <- args[3] %||% sprintf("output-%s-%s.csv", start, end)
passKO        <- as.character(read.table("/export2/home/uesu/reDiamond/data/KOlist_700")$V1)    #kos
msaDIR        <- args[4] %||% "/export2/home/uesu/reDiamond/out/assm.0500"
contigsDIR    <- args[5] %||% "/export2/home/uesu/reDiamond/out/.assm.0200.newbler"
cutContigsDIR <- args[6] %||% "/export2/home/uesu/reDiamond/out/assm.0700.chosenContigWindwsfna.200"
blatOutputDIR <- args[7] %||% "/export2/home/uesu/reDiamond/out/map.0200/blatoutputDIR";
fastaFile     <- args[8] %||% "/export2/home/uesu/reDiamond/out/map.0200/fasta"

dir.create(blatOutputDIR)
dir.create(fastaFile)

allKODF <- passKO[start:end] %>%
mclapply(function(ko){
    message(paste0("Working on ", ko))
    #Loading MSA
    upMSA <-  readDNAStringSet(sprintf("%s/%s.msa", msaDIR, ko))                       #Actual Ulu pandan NOT simulation
    upMSA@ranges@NAMES %<>% strsplit("\\s") %>% sapply(function(x) x[[1]]) #remove names

    #' We convert each contig/gene in the msa, into IRanges obj, then subset and classify those within the max diversity region into spanning, non-spanning, inside.
    #IRanges of contigs and their locations which correspond with the global MSA
    upMSA_int    <-  buildIntervals(stringSet = upMSA, minIntervalSize = 200)
    #the maxdiveristy window
    MDlocation   <-  extractLoc(fileName=sprintf("%s/%s.fna", cutContigsDIR, ko)) 
    #subset contigs which fall in the maxdiveristy regino and assign identity
    contigsSpan  <-  spanningOrNot(upMSA_int, MDlocation$startLoc, MDlocation$endLoc) 

    #' #### Mapping mRNA reads to the contigs
    #' We first perform a high identity alignment (using blat) of the reads to the contigs originally in the msa and later extract out only alignments which fall within the region of interest.

    #+ mRNA-blat
    contigs   <- sprintf("%s/%s/454AllContigs.fna", contigsDIR, ko)
    fastqFile <- 1:2 %>% sapply(function(read)sprintf("%s/%s/input/%s.", contigsDIR, ko, ko) %>% paste0(.,read, ".fq"))
    fastaFile %<>% paste0("/",ko, ".fasta")

    #' We look for transcript alignments (blat output) which overlap max diversity region.

    #Step 1: Mapping reads to their contigs
    blatDF     <- map(fq = fastqFile, fa = fastaFile, blatOutputFile = sprintf("%s/blatoutput_%s.m8", blatOutputDIR, ko), contigs = contigs)
    #Build intervals
    blatRanges <- blatDF %$% GRanges(seqnames=subject, ranges=IRanges(newStart, newEnd), read=query)
    #because the reads were consist of both gDNA and cDNA
    cDNARanges <- mcols(blatRanges)@listData[[1]] %>% as.character %>% grepl(pattern="cDNA", x=.) %>% blatRanges[.]
    gDNARanges <- mcols(blatRanges)@listData[[1]] %>% as.character %>% grepl(pattern="gDNA", x=.) %>% blatRanges[.]

    #Step 2: Built interval range for reads mapping to contigs in the regions of interest
    #99 due to the read size
    ROI <- extendRange(contigsSpan$spanning, contigsSpan$trunc, contigsSpan$inside,99,width(upMSA)[1])
    #Step 3: Find overlap
    countDF_cDNA <-     make.map.table(contigs = ROI, reads = cDNARanges, datype="cDNA")
    countDF_cDNA %<>%   mutate(ko = ko)

    countDF_gDNA <-     make.map.table(contigs = ROI, reads = gDNARanges, datype="gDNA")
    countDF_gDNA %<>%   mutate(ko = ko)

    list(
        gDNA = countDF_gDNA,
        cDNA = countDF_cDNA
        )

combinedDF = merge(
      countDF_cDNA %>% select(-type, -trueLength, -contigID),
      countDF_gDNA %>% select(-type, -trueLength, -contigID),
      by = c("contigName","ko"),
      suffixes = c("_cDNA", "_gDNA"),
      all=T
)
combinedDF[is.na(combinedDF)] = 0
combinedDF
})

write.csv(allKODF, file=finalOutput)
