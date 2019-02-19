#' Function for single sample weights estimation
#'
#' This function allows you to calculate singel sampmle weights
#'
#' @param refSet Label of pre-stored deconvoluted reference dataset
#' @param data Data matrix to be deconvoluted with gene ID as rownames
#' @param geneIDType Choose one among "geneSymbol", "EntrezID" and "TCGA"
#' @return A matrix of the estimated sample compartment weights
#' @export

Decon_single_sample <- function(refSet,
                                data,
                                geneIDType){
  # load reference
  refSet <- get(refSet)
  #load(paste("./data/",refSet,".RData",sep=""))
  indPri <- which(refSet$compInfo$Type == "Primary")
  geneSigRef <- refSet$geneWeight[,indPri]
  dataRef = log2(1+refSet$expr)
  if(geneIDType == "TCGA"){
    geneIDRef = refSet$geneInfo$geneID
  } else if(geneIDType == "EntrezID"){
    geneIDRef = refSet$geneInfo$EntrezID
  } else if(geneIDType == "geneSymbol"){
    geneIDRef = refSet$geneInfo$geneSymbol
  } else{
    warnings("geneIDType not supported...")
  }

  # overlap genes
  geneID <- rownames(data)
  geneOvlp <- geneIDRef[geneIDRef %in% geneID]
  indRef <- match(geneOvlp,geneIDRef)
  dataRef <- dataRef[indRef,]
  geneSigRef <- geneSigRef[indRef,]
  indCur <- match(geneOvlp,geneID)
  data <- data[indCur,]

  # normalize data
  maxRef <- max(dataRef)
  data <- (maxRef/max(data))*data

  # package a nnls function
  callNNLS <- function (A,b){
    res <- nnls::nnls(A,b)$x
    return(res)
  }

  # calculate sample weights
  sampleWeights <- t(apply(data, 2, function(x) callNNLS(as.matrix(geneSigRef),x)))
  colnames(sampleWeights) <- colnames(geneSigRef)

  return(sampleWeights)
}


