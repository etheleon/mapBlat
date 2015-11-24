#' Assigns
#'
#' annotateContigs
#' @param blast2lca output file showing taxonomic assignment
#' @export
annotateContigs.taxonomy <- function(blast2lca)
{

}

#' blast2lcaIN
#' reads in and executes perl
#' @param blast2lca output file showing taxonomic assignment
#' @param targetFile path to file which will store the intermediate output from the perl executable
#' @param keepFile to keep file or not
#' @export
blast2lcaIn <- function(blast2lca, targetFile = "blast2lca.csv", keepFile=FALSE){
    perl = findPerl();
    if(targetFile == "blast2lca.csv") targetFile <- paste0(getwd(), '/', targetFile)

    script  <-   find.package('mapBlat') %>% file.path('perl') %>% file.path('mapTaxon.pl')
    ncbimap <-   find.package('mapBlat') %>% file.path('MEGANdata') %>% file.path('ncbi.map')

    cmd <- paste(shQuote(perl),
                     shQuote(script),
                      shQuote(ncbimap),
                      shQuote(blast2lca),
                      shQuote(targetFile),
                      sep=" ")

    results <- try(system(cmd))
    adf = data.table::fread(targetFile) #might be quite big
    if(keepFile != TRUE) file.remove(targetFile)
    adf
}

#' findPerl finds the path to the executable 
#'
#' Code taken from R package gdata
#'
#' @param perl full path to perl executable
#' @param verbose

findPerl <- function(perl, verbose = "FALSE")
{
  errorMsg <- "perl executable not found. Use perl= argument to specify the correct path."

  if (missing(perl))
    {
      perl = "perl"
    }

  perl = Sys.which(perl)
  if (perl=="" || perl=="perl")
    stop(errorMsg)

  if (.Platform$OS == "windows") {
    if (length(grep("rtools", tolower(perl))) > 0) {
      perl.ftype <- shell("ftype perl", intern = TRUE)
      if (length(grep("^perl=", perl.ftype)) > 0) {
        perl <- sub('^perl="([^"]*)".*', "\\1", perl.ftype)
      }
    }
      }
  
  if (verbose) cat("Using perl at", perl, "\n")
  
  perl
}


