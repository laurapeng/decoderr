#' Function for sample weight normalization for PDAC
#'
#'
#'
#' @param sampleWeight Sample weights calculated by Decon_single_sample.R
#' @return Normalized sample weights
#' @export



Norm_PDAC_weights <- function(sampleWeight){
  sampleWeight <- sampleWeight[,c(9,5,4,7,2,1,3)]
  sampleWeightNorm <- data.frame(sweep(sampleWeight, 1, rowSums(sampleWeight), FUN="/"))
  names(sampleWeightNorm) <- c("BasalTumor","ClassicalTumor",
                               "ActivatedStroma","NormalStroma",
                               "Immune","Endocrine","Exocrine")

  sampleWeightNorm$bcRatio <- (sampleWeightNorm$BasalTumor+0.01)/(sampleWeightNorm$ClassicalTumor+0.01)
  sampleWeightNorm$bcRatio[sampleWeightNorm$bcRatio > 5] <- 5
  sampleWeightNorm$bcDiff <- sampleWeightNorm$BasalTumor-sampleWeightNorm$ClassicalTumor
  sampleWeightNorm$bcSum <- sampleWeightNorm$BasalTumor+sampleWeightNorm$ClassicalTumor

  sampleWeightNorm$anRatio <- (sampleWeightNorm$ActivatedStroma+0.01)/(sampleWeightNorm$NormalStroma+0.01)
  sampleWeightNorm$anRatio[sampleWeightNorm$anRatio > 5] <- 5
  sampleWeightNorm$anDiff <- sampleWeightNorm$ActivatedStroma-sampleWeightNorm$NormalStroma
  sampleWeightNorm$anSum <- sampleWeightNorm$ActivatedStroma+sampleWeightNorm$NormalStroma

  sampleWeightNorm$aiRatio <- (sampleWeightNorm$ActivatedStroma+0.01)/(sampleWeightNorm$Immune+0.01)
  sampleWeightNorm$aiRatio[sampleWeightNorm$aiRatio > 5] <- 5
  sampleWeightNorm$aiDiff <- sampleWeightNorm$ActivatedStroma-sampleWeightNorm$Immune
  sampleWeightNorm$aiSum <- sampleWeightNorm$ActivatedStroma+sampleWeightNorm$Immune

  sampleWeightNorm$aniRatio <- (sampleWeightNorm$ActivatedStroma+0.01)/(sampleWeightNorm$NormalStroma+sampleWeightNorm$Immune+0.01)
  sampleWeightNorm$aniRatio[sampleWeightNorm$aniRatio > 5] <- 5
  sampleWeightNorm$aniDiff <- sampleWeightNorm$ActivatedStroma-sampleWeightNorm$NormalStroma-sampleWeightNorm$Immune
  sampleWeightNorm$aniSum <- sampleWeightNorm$ActivatedStroma+sampleWeightNorm$NormalStroma+sampleWeightNorm$Immune

  return(sampleWeightNorm)
}
