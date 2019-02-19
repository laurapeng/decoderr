#' Function for sample weight normalization for PDAC
#'
#' 
#' 
#' @param sampleWeight Sample weights calculated by Decon_single_sample.R
#' @return Normalized sample weights
#' @export



Norm_PDAC_weights <- function(sampleWeight){
  sampleWeight <- sampleWeight[,c(9,5,4,7,2,1,3,6,8)]
  sampleWeightNorm <- data.frame(sweep(sampleWeight, 1, rowSums(sampleWeight), FUN="/"))
  names(sampleWeightNorm) <- c("BasalTumor","ClassicalTumor",
                               "ActivatedStroma","NormalStroma",
                               "Immune","Exocrine","Endocrine","Histone","Olfactory")
  sampleWeightNorm$bcRatio <- (sampleWeightNorm$BasalTumor+0.01)/(sampleWeightNorm$ClassicalTumor+0.01)
  sampleWeightNorm$bcRatio[sampleWeightNorm$bcRatio > 5] <- 5
  sampleWeightNorm$anRatio <- (sampleWeightNorm$ActivatedStroma+0.01)/(sampleWeightNorm$NormalStroma+0.01)
  sampleWeightNorm$anRatio[sampleWeightNorm$anRatio > 5] <- 5
  return(sampleWeightNorm)
}