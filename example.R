#!/usr/bin/env Rscript

args <- commandArgs(T) #in this case you should do `ls out/assm.0700.chosenContigWindwsfna.200/*fna | wc -l`

suppressPackageStartupMessages(library(mapBlat))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(MetamapsDB))

scg = read.table("/export2/home/uesu/simulation_fr_the_beginning/data/single.copy.gene") %>% setNames(c("ko", "sym"))
connect(url="192.168.100.253")
scg %<>% cbind(scg$ko %>% 
               lapply(function(koid){ koname(koid)}) %>% 
               do.call(rbind,.) %>% 
               select(-ko.ko)
               )


`%||%` <- function(x, y){
      if (is.na(x)) y else x
}

start         <- args[1] %||% 1
end           <- args[2] %||% 1

finalOutput   <- args[3] %||% sprintf("output-%s.csv", start)
passKO        <- as.character(read.table("/export2/home/uesu/reDiamond/data/KOlist_700")$V1)    #kos
msaDIR        <- args[4] %||% "/export2/home/uesu/reDiamond/out/assm.0500"
contigsDIR    <- args[5] %||% "/export2/home/uesu/reDiamond/out/.assm.0200.newbler"
cutContigsDIR <- args[6] %||% "/export2/home/uesu/reDiamond/out/assm.0700.chosenContigWindwsfna.200"
blatOutputDIR <- args[7] %||% "/export2/home/uesu/reDiamond/out/map.0200/blatoutputDIR";
fastaFile     <- args[8] %||% "/export2/home/uesu/reDiamond/out/map.0200/fasta"

dir.create(blatOutputDIR)
dir.create(fastaFile)

#allKODF <- passKO[start:end] %>%
#mclapply(function(ko){
ko = passKO[as.integer(start)]
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
    print("Step 1: Mapping reads to their contigs")
    blatDF     <- mapper(fq = fastqFile, fa = fastaFile, blatOutputFile = sprintf("%s/blatoutput_%s.m8", blatOutputDIR, ko), contigs = contigs)
    #blatDF_unique = blatDF[!blatDF$query %in% blatDF$query[duplicated(blatDF$query)],]
    #Build intervals

    blatDF %<>%
    blatRanges <- GRanges(seqnames=blatDF$subject, ranges=IRanges(blatDF$newStart, blatDF$newEnd), read=blatDF$query)
    #because the reads were consist of both gDNA and cDNA

    cDNARanges <- mcols(blatRanges)@listData[[1]] %>% as.character %>% grepl(pattern="cDNA", x=.) %>% blatRanges[.]
    gDNARanges <- mcols(blatRanges)@listData[[1]] %>% as.character %>% grepl(pattern="gDNA", x=.) %>% blatRanges[.]

    print("Step 2: Built interval range for reads mapping to contigs in the regions of interest")
    #99 due to the read size
    ROI <- extendRange(contigsSpan$spanning, contigsSpan$trunc, contigsSpan$inside,99,width(upMSA)[1])
    #Step 3: Find overlap
    countDF_cDNA <-     make.map.table(contigs = ROI, reads = cDNARanges, datype="cDNA")
    countDF_cDNA %<>%   mutate(ko = ko)

    countDF_gDNA <-     make.map.table(contigs = ROI, reads = gDNARanges, datype="gDNA")
    countDF_gDNA %<>%   mutate(ko = ko)

#    list(
#        gDNA = countDF_gDNA,
#        cDNA = countDF_cDNA
#        )

combinedDF <- merge(
      countDF_cDNA %>% select(-type, -trueLength, -contigID),
      countDF_gDNA %>% select(-type, -trueLength, -contigID),
      by = c("contigName","ko"),
      suffixes = c("_cDNA", "_gDNA"),
      all=T
)
combinedDF[is.na(combinedDF)] = 0
#combinedDF
#}, mc.cores=10)

combinedDF$rpkm_gDNA = with(combinedDF, rpkm(Freq_gDNA, sum(Freq_gDNA), 200))
combinedDF$rpkm_cDNA = with(combinedDF, rpkm(Freq_cDNA, sum(Freq_cDNA), 200))


#cant use rpkm just on one count
rpkm       =  function(count, allcount, length){
    #change to numeric (bit size too small for integers)
    count     %<>%  as.numeric
    allcount  %<>%  as.numeric
    length    %<>%  as.numeric(length)
    (10^9 * count) / (allcount* length)
}

combinedDF$type    = 'full'
combinedDF_nr$type = 'nr'

cd = rbind(combinedDF, combinedDF_nr)

##################################################
pdf("aplot_scg.pdf")
##################################################

ggplot(cd, aes(rpkm_gDNA,rpkm_cDNA)) +
geom_point(aes(color=type))          +
scale_y_log10()                      +
scale_x_log10()

ggplot(combinedDF, aes(reorder(contigName, rpkm_gDNA),rpkm_gDNA)) +
geom_bar(stat="identity")                                         +
theme(axis.text.x=element_blank())+ ggtitle("gDNA")

ggplot(combinedDF, aes(reorder(contigName, rpkm_gDNA),rpkm_gDNA)) +
geom_bar(stat="identity")                                         +
scale_y_log10()                                                   +
theme(axis.text.x=element_blank())+ ggtitle("gDNA")


ggplot(combinedDF, aes(reorder(contigName, rpkm_cDNA),rpkm_cDNA)) +
geom_bar(stat="identity")                                         +
theme(axis.text.x=element_blank())+ ggtitle("cDNA")

ggplot(combinedDF, aes(reorder(contigName, rpkm_cDNA),rpkm_cDNA)) +
geom_bar(stat="identity")                                         +
scale_y_log10()                                                   +
theme(axis.text.x=element_blank())+ ggtitle("cDNA")


ggplot(combinedDF, aes(reorder(contigName, rpkm_gDNA),rpkm_cDNA)) +
geom_bar(stat="identity")                                         +
theme(axis.text.x=element_blank())+ ggtitle("cDNA")

ggplot(combinedDF, aes(reorder(contigName, rpkm_gDNA),rpkm_cDNA)) +
geom_bar(stat="identity")                                         +
scale_y_log10()                                                   +
theme(axis.text.x=element_blank())+ ggtitle("cDNA")


ggplot(combinedDF, aes(reorder(contigName, rpkm_gDNA),rpkm_cDNA/rpkm_gDNA)) +
geom_bar(stat="identity")                                         +
scale_y_log10()                                                   +
theme(axis.text.x=element_blank())
dev.off()




combinedDF %$% reorder(contigName, rpkm_gDNA) %>% head
combinedDF %$% reorder(contigName, rpkm_gDNA) %>% head


#write.csv(combinedDF, file=finalOutput)
