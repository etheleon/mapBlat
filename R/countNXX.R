#' Assigns the a NXX to KOs
#'
#' \code{countNXX} takes the RPKM expression value and returns a NXX value describing the number 
#' of contigs required to account for X% of expression eg. 50% of expression
#'
#' @param countDF input data.frame
#' @param perc the percentages given as a sequence of 
#' @param abundanceCol abundance
#'
#' @return a data.frame obj containing the NXX values.
#' We define NXX as metric for describing the number of genes required
#' to satisfy X% total expression within the KO
#'
#' @export
countNXX <- function(countDF, perc=seq(0.1,1,0.05), abundanceCol="rpkm_cDNA"){
    if(length(unique(countDF$ko)) > 1){
        unique(countDF$ko) %>%
        mclapply(function(koID) countDF %>% filter(ko == koID) %>% countInside(abundanceCol, perc)) %>% do.call(rbind,.)
    }else{
        countInside(countDF, abundanceCol, perc)
    }
}

#' counts the rolling abundance scores
#'
#' \code{countInside}
#'
#' @param koID the KO id
#' @param countDF df from the countNXX function
#'
countInside <- function(countDF, abundanceCol, perc=perc){
    abundanceCol %<>% parse(text=.)
    countDF %<>% arrange(desc(eval(abundanceCol, countDF)))
    ko      =   unique(countDF$ko)
    total   <-  sum(eval(abundanceCol, countDF))
   if(total == 0){
       warning("KO has not mappable reads to MD region")
       data.frame(percentage=numeric(), n=integer(), total=numeric(), contigsRequired=integer(),ko=factor())
   }else{
       rollexp         = 0
       roll = 1:nrow(countDF) %>% sapply(function(x) {sumexact(countDF[1:x,as.character(abundanceCol)])})
       perc %>% lapply(function(percentage){
                       NXX             = percentage * total
                       minRoll         = max(which(roll < NXX))
                       n               = nrow(countDF)
                       contigsRequired = ifelse(is.infinite(minRoll), 1, minRoll)
                       data.frame(percentage,
                                  n,
                                  total,
                                  contigsRequired,
                                  ko)
        }) %>% do.call(rbind,.)
   }
}
