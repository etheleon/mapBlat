#!/usr/bin/env Rscript

args <- commandArgs(T) #in this case you should do `ls out/assm.0700.chosenContigWindwsfna.200/*fna | wc -l`

suppressPackageStartupMessages(library(mapBlat))
suppressPackageStartupMessages(library(ggplot2))
theme_set(theme_bw())
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(MetamapsDB))
suppressPackageStartupMessages(library(GenomicRanges))
load("~/abu_genes/data/Abundance_Function/pathways.rda")
ko <- args[1] %||% 'K00001'
message(paste0("Working on ", ko))

passKO         =  as.character(read.table("data/KOlist_700")$V1)    #kos
msaDIR         =  "/export2/home/uesu/reDiamond/out/assm.0500"
contigsDIR     =  "/export2/home/uesu/reDiamond/out/.assm.0200.newbler"
cutContigsDIR  =  "/export2/home/uesu/reDiamond/out/assm.0700.chosenContigWindwsfna.200"
blatOutputDIR  =  "out/map.0102_indiv_ko";



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

#Loading MSA
upMSA <-  readDNAStringSet(sprintf("/export2/home/uesu/reDiamond/out/assm.0500/%s.msa",ko))                       #Actual Ulu pandan NOT simulation
upMSA@ranges@NAMES %<>% strsplit("\\s") %>% sapply(function(x) x[[1]]) #remove names

#' We convert each contig/gene in the msa, into IRanges obj, then subset and classify those within the max diversity region into spanning, non-spanning, inside.
#IRanges of contigs and their locations which correspond with the global MSA
upMSA_int    =  buildIntervals(stringSet = upMSA, minIntervalSize = 200) 
#the maxdiveristy window
MDlocation   =  extractLoc(fileName=sprintf("/export2/home/uesu/reDiamond/out/assm.0700.chosenContigWindwsfna.200/%s.fna", ko)) 
#subset contigs which fall in the maxdiveristy regino and assign identity
contigsSpan  =  spanningOrNot(upMSA_int, MDlocation$startLoc, MDlocation$endLoc) 

#' #### Mapping mRNA reads to the contigs
#' We first perform a high identity alignment (using blat) of the reads to the contigs originally in the msa and later extract out only alignments which fall within the region of interest.

#+ mRNA-blat
contigs   <- sprintf("out/.assm.0200.newbler/%s/454AllContigs.fna", ko)
fastqFile <- 1:2 %>% sapply(function(read)sprintf("/export2/home/uesu/reDiamond/out/.assm.0200.newbler/%s/input/%s.", ko, ko) %>% paste0(.,read, ".fq"))
fastaFile <- sprintf("out/map.0100/%s.fasta", ko)

#' We look for transcript alignments (blat output) which overlap max diversity region.

#Step 1: Mapping reads to their contigs
blatDF     <- map(fq = fastqFile, fa = fastaFile, blatOutputFile = sprintf("out/map.0100/blatoutput_%s.m8", ko), contigs = contigs)
#Build intervals
blatRanges <- blatDF %$% GRanges(seqnames=subject, ranges=IRanges(newStart, newEnd), read=query)
#because the reads were consist of both gDNA and cDNA
cDNARanges <- mcols(blatRanges)@listData[[1]] %>% as.character %>% grepl(pattern="cDNA", x=.) %>% blatRanges[.]
gDNARanges <- mcols(blatRanges)@listData[[1]] %>% as.character %>% grepl(pattern="gDNA", x=.) %>% blatRanges[.]
#Step 2: Built interval range for reads mapping to contigs in the regions of interest
regionOfInterest = extendRange(contigsSpan$spanning, contigsSpan$trunc, contigsSpan$inside,99,width(upMSA)[1])
#Step 3: Find overlap
countDF <-  make.map.table(contigs = regionOfInterest, reads = cDNARanges) %>%
            mutate(ko = ko)

