#' @title GENIE3
#'
#' @description \code{GENIE3} Infers a gene regulatory network (in the form of a weighted adjacency matrix) from expression data, using ensembles of regression trees.
#'
#' @param exprMatrix Expression matrix (genes x samples). Every row is a gene, every column is a sample.
#' The expression matrix can also be provided as one of the Bioconductor classes:
#' \itemize{
#' \item \link[Biobase]{ExpressionSet}: The matrix will be obtained through exprs(exprMatrix)
#' \item \link[SummarizedExperiment]{RangedSummarizedExperiment}: The matrix will be obtained through assay(exprMatrix), wich will extract the first assay (usually the counts)
#' }
#' @param treeMethod Tree-based method used. Must be either "RF" for Random Forests (default) or "ET" for Extra-Trees.
#' @param K Number of candidate regulators randomly selected at each tree node (for the determination of the best split). Must be either "sqrt" for the square root of the total number of candidate regulators (default), "all" for the total number of candidate regulators, or a stricly positive integer.
#' @param nTrees Number of trees in an ensemble for each target gene. Default: 1000.
#' @param regulators Subset of genes used as candidate regulators. Must be either a vector of indices, e.g. \code{c(1,5,6,7)}, or a vector of gene names, e.g. \code{c("at_12377", "at_10912")}. The default value NULL means that all the genes are used as candidate regulators.
#' @param targets Subset of genes to which potential regulators will be calculated. Must be either a vector of indices, e.g. \code{c(1,5,6,7)}, or a vector of gene names, e.g. \code{c("at_12377", "at_10912")}. If NULL (default), regulators will be calculated for all genes in the input matrix.
#' @param nCores Number of cores to use for parallel computing. Default: 1.
#' @param cluster Cluster for parallel computing on multiple nodes. Default: NULL.
#' @param progress Monitor the progression of the computation. Default: FALSE
#' @param verbose If set to TRUE, a feedback on the progress of the calculations is given. Default: FALSE.
#'
#' @return Weighted adjacency matrix of inferred network. Element w_ij (row i, column j) gives the importance of the link from regulatory gene i to target gene j.
#'
#' @examples
#' ## Generate fake expression matrix
#' exprMatrix <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
#' rownames(exprMatrix) <- paste("Gene", 1:20, sep="")
#' colnames(exprMatrix) <- paste("Sample", 1:5, sep="")
#'
#' ## Run GENIE3
#' set.seed(123) # For reproducibility of results
#' weightMatrix <- GENIE3(exprMatrix, regulators=paste("Gene", 1:5, sep=""))
#'
#' ## Get ranking of edges
#' linkList <- getLinkList(weightMatrix)
#' head(linkList)
#' @export
setGeneric("GENIE3", signature="exprMatrix",
           function(exprMatrix, treeMethod="RF", K="sqrt", nTrees=1000, regulators=NULL, targets=NULL, nCores=1, cluster=NULL, progress=FALSE, verbose=FALSE)
           {
             standardGeneric("GENIE3")
           })

#' @export
setMethod("GENIE3", "matrix",
          function(exprMatrix, treeMethod="RF", K="sqrt", nTrees=1000, regulators=NULL, targets=NULL, nCores=1, cluster=NULL, progress=FALSE, verbose=FALSE)
          {
            .GENIE3(exprMatrix=exprMatrix, treeMethod=treeMethod, K=K, nTrees=nTrees, regulators=regulators, targets=targets, nCores=nCores, cluster=cluster, progress=progress, verbose=verbose)
          })

#' @export
setMethod("GENIE3", "SummarizedExperiment",
          function(exprMatrix, treeMethod="RF", K="sqrt", nTrees=1000, regulators=NULL, targets=NULL, nCores=1, cluster=NULL, progress=FALSE, verbose=FALSE)
          {
            if(length(SummarizedExperiment::assays(exprMatrix))>1) warning("More than 1 assays are available. Only using the first one.")
            exprMatrix <- SummarizedExperiment::assay(exprMatrix)
            .GENIE3(exprMatrix=exprMatrix, treeMethod=treeMethod, K=K, nTrees=nTrees, regulators=regulators, targets=targets, nCores=nCores, cluster=cluster, progress=progress, verbose=verbose)
          })

#' @export
setMethod("GENIE3", "ExpressionSet",
          function(exprMatrix, treeMethod="RF", K="sqrt", nTrees=1000, regulators=NULL, targets=NULL, nCores=1, cluster=NULL, progress=FALSE, verbose=FALSE)
          {
            exprMatrix <- Biobase::exprs(exprMatrix)
            .GENIE3(exprMatrix=exprMatrix, treeMethod=treeMethod, K=K, nTrees=nTrees, regulators=regulators, targets=targets, nCores=nCores, cluster=cluster, progress=progress, verbose=verbose)
          })

.GENIE3 <- function(exprMatrix, treeMethod, K, nTrees, regulators, targets, nCores, cluster, progress, verbose)
{
  .checkArguments(exprMatrix=exprMatrix, treeMethod=treeMethod, K=K, nTrees=nTrees, regulators=regulators, targets=targets, nCores=nCores, cluster=cluster, progress=progress, verbose=verbose)
  
  if(is.numeric(regulators)) regulators <- rownames(exprMatrix)[regulators]
  
  ############################################################
  # transpose expression matrix to (samples x genes)
  exprMatrix <- t(exprMatrix)
  num.samples <- nrow(exprMatrix)
  allGeneNames <- colnames(exprMatrix)
  
  # get names of input genes
  if(is.null(regulators))
  {
    regulatorNames <- allGeneNames
  } else
  {
    # input gene indices given as integers
    if (is.numeric(regulators))
    {
      regulatorNames <- allGeneNames[regulators]
      # input gene indices given as names
    } else
    {
      regulatorNames <- regulators
      # for security, abort if some input gene name is not in gene names
      missingGeneNames <- setdiff(regulatorNames, allGeneNames)
      if (length(missingGeneNames) != 0) stop(paste("Regulator genes missing from the expression matrix:", paste(missingGeneNames, collapse=", ")))
    }
  }
  rm(regulators)
  
  # get names of target genes
  if(is.null(targets))
  {
    targetNames <- allGeneNames
  } else
  {
    # input gene indices given as integers
    if (is.numeric(targets))
    {
      targetNames <- allGeneNames[targets]
      # input gene indices given as names
    } else
    {
      targetNames <- targets
      # for security, abort if some input gene name is not in gene names
      missingGeneNames <- setdiff(targetNames, allGeneNames)
      if (length(missingGeneNames) != 0) stop(paste("Target genes missing from the expression matrix:", paste(missingGeneNames, collapse=", ")))
    }
  }
  
  nGenes <- length(targetNames)
  rm(targets)
  
  # tree method
  if (treeMethod == 'RF')
  {
    RF_randomisation <- 1
    ET_randomisation <- 0
    bootstrap_sampling <- 1
  } else {
    RF_randomisation <- 0
    ET_randomisation <- 1
    bootstrap_sampling <- 0
  }
  
  if (verbose) message(paste("Tree method: ", treeMethod,
                             "\nK: ", K,
                             "\nNumber of trees: ", nTrees, sep=""))
  # other default parameters
  nmin <- 1
  permutation_importance <- 0
  
  # setup weight matrix
  weightMatrix <- matrix(0.0, nrow=length(regulatorNames), ncol=nGenes)
  rownames(weightMatrix) <- regulatorNames
  colnames(weightMatrix) <- targetNames
  
  # compute importances for every target gene
  library(doSNOW)
  library(doRNG)
  library(foreach)
  if(!is.null(cluster)) {
    
    if(cluster$config$type == "raw") {
      # Load utils to build the cluster
      source("https://raw.githubusercontent.com/mase5/r-utils/master/parallel_utils.R")
      # Build the cluster
      message("Building the PS cluster...")
      cluster<-open_PSC(user = cluster$config$def$user, nodes = cluster$config$def$nodes, n.cores = cluster$config$def$n.cores, verbose = cluster$config$def$verbose)
      
    } else if(cluster$config$type == "psc") {
      cluster <- cluster$config$def$cluster
    } else {
      stop("The given cluster object is invalid!")
    }
    
    weightMatrix.reg <- .run(expr.matrix = exprMatrix
                             , weight.matrix = weightMatrix
                             , monitor.progress = progress
                             , target.names = targetNames
                             , regulator.names = regulatorNames
                             , num.samples = num.samples
                             , n.min = nmin
                             , ET.randomisation = ET_randomisation
                             , RF.randomisation = RF_randomisation
                             , K = K
                             , n.trees = nTrees
                             , bootstrap.sampling = bootstrap_sampling
                             , permutation.importance = permutation_importance)
    if(cluster$config$type == "raw") {
      # Close the PS cluster
      close_PSC(cluster = cluster$config$def$cluster, verbose = T)
    }
  } else if(nCores==1)
  {
    # serial computing
    if(verbose) message("Using 1 core.")
    for(targetName in targetNames)
    {
      if(verbose) message(paste("Computing gene ", which(targetNames == targetName), "/", nGenes, ": ",targetName, sep=""))
      
      # remove target gene from input genes
      theseRegulatorNames <- setdiff(regulatorNames, targetName)
      numRegulators <- length(theseRegulatorNames)
      mtry <- .setMtry(K, numRegulators)
      
      x <- exprMatrix[,theseRegulatorNames]
      y <- exprMatrix[,targetName]
      
      im <- .C("BuildTreeEns",as.integer(num.samples),as.integer(numRegulators),
               as.single(c(x)),as.single(c(y)),as.integer(nmin),
               as.integer(ET_randomisation),as.integer(RF_randomisation),
               as.integer(mtry),as.integer(nTrees),
               as.integer(bootstrap_sampling),as.integer(permutation_importance),
               as.double(vector("double",numRegulators)))[[12]]
      
      # normalize variable importances
      im <- im / sum(im)
      weightMatrix[theseRegulatorNames, targetName] <- im
    }
  } else
  {
    # requireNamespace("foreach"); requireNamespace("doRNG"); requireNamespace("doParallel")
    
    # parallel computing
    doParallel::registerDoParallel(); options(cores=nCores)
    if(verbose) message(paste("\nUsing", foreach::getDoParWorkers(), "cores."))
    
    weightMatrix.reg <- .run(expr.matrix = exprMatrix
                             , weight.matrix = weightMatrix
                             , monitor.progress = progress
                             , target.names = targetNames
                             , regulator.names = regulatorNames
                             , num.samples = num.samples
                             , n.min = nmin
                             , ET.randomisation = ET_randomisation
                             , RF.randomisation = RF_randomisation
                             , K = K
                             , n.trees = nTrees
                             , bootstrap.sampling = bootstrap_sampling
                             , permutation.importance = permutation_importance)
  }
  return(weightMatrix.reg)
}

#' Wrapper function tu run GENIE3 with progress bar
#' @author M.D.W
#'
.run <- function(expr.matrix
                 , weight.matrix
                 , monitor.progress = T
                 , target.names
                 , regulator.names
                 , num.samples
                 , n.min
                 , ET.randomisation
                 , RF.randomisation
                 , K
                 , n.trees
                 , bootstrap.sampling
                 , permutation.importance) {
  if(monitor.progress) {
    # Monitoring the progress
    pb <- txtProgressBar(min=1, max=length(target.names), style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
  } else {
    opts <- NULL
  }
  K.arg<-K
  expr.matrix.arg<-expr.matrix
  # weightMatrix.reg <- foreach::foreach(targetName=targetNames, .combine=cbind) %dorng%
  "%dopar%"<- foreach::"%dopar%"
  suppressPackageStartupMessages(weightMatrix.reg <- doRNG::"%dorng%"(foreach::foreach(targetName=target.names, .combine=cbind, .options.snow = opts), 
                                                                      {
                                                                        # remove target gene from input genes
                                                                        theseRegulatorNames <- setdiff(regulator.names, targetName)
                                                                        numRegulators <- length(theseRegulatorNames)
                                                                        mtry <- .setMtry(K.arg, numRegulators)
                                                                        
                                                                        x <- expr.matrix.arg[,theseRegulatorNames]
                                                                        y <- expr.matrix.arg[,targetName]
                                                                        
                                                                        im <- .C("BuildTreeEns", as.integer(num.samples), as.integer(numRegulators),
                                                                                 as.single(c(x)),as.single(c(y)), as.integer(n.min),
                                                                                 as.integer(ET.randomisation), as.integer(RF.randomisation),
                                                                                 as.integer(mtry), as.integer(n.trees),
                                                                                 as.integer(bootstrap.sampling), as.integer(permutation.importance),
                                                                                 as.double(vector("double",numRegulators)))[[12]]
                                                                        
                                                                        # normalize variable importances
                                                                        im <- im / sum(im)
                                                                        
                                                                        c(setNames(0, targetName), setNames(im, theseRegulatorNames))[regulator.names]
                                                                      }))
  attr(weightMatrix.reg, "rng") <- NULL
  weight.matrix[regulator.names,] <- weightMatrix.reg
  
  # Close the monitor progress
  if(monitor.progress) {
    close(pb)
  }
  return (weightMatrix.reg)
}

# mtry <- setMtry(K, numRegulators)
.setMtry <- function(K, numRegulators)
{
  # set mtry
  if (class(K) == "numeric") {
    mtry <- K
  } else if (K == "sqrt") {
    mtry <- round(sqrt(numRegulators))
  } else {
    mtry <- numRegulators
  }
  
  return(mtry)
}

.checkArguments <- function(exprMatrix, treeMethod, K, nTrees, regulators, targets, nCores, cluster, progress, verbose)
{
  ############################################################
  # check input arguments
  if (!is.matrix(exprMatrix) && !is.array(exprMatrix)) {
    stop("Parameter exprMatrix must be a two-dimensional matrix where each row corresponds to a gene and each column corresponds to a condition/sample.")
  }
  
  if (length(dim(exprMatrix)) != 2) {
    stop("Parameter exprMatrix must be a two-dimensional matrix where each row corresponds to a gene and each column corresponds to a condition/sample.")
  }
  
  if (is.null(rownames(exprMatrix))) {
    stop("exprMatrix must specify the names of the genes in rownames(exprMatrix).")
  }
  
  if (treeMethod != "RF" && treeMethod != "ET") {
    stop("Parameter treeMethod must be \"RF\" (Random Forests) or \"ET\" (Extra-Trees).")
  }
  
  if (K != "sqrt" && K != "all" && !is.numeric(K)) {
    stop("Parameter K must be \"sqrt\", or \"all\", or a strictly positive integer.")
  }
  
  if (is.numeric(K) && K<1) {
    stop("Parameter K must be \"sqrt\", or \"all\", or a strictly positive integer.")
  }
  
  if (!is.numeric(nTrees) || nTrees<1) {
    stop("Parameter nTrees should be a stricly positive integer.")
  }
  
  if (!is.null(regulators))
  {
    if(length(regulators)<2) stop("Provide at least 2 potential regulators.")
    
    if (!is.vector(regulators)) {
      stop("Parameter regulators must be either a vector of indices or a vector of gene names.")
    }
    
    if (is.character(regulators) && length(intersect(regulators,rownames(exprMatrix))) == 0) {
      stop("The genes must contain at least one candidate regulator.")
    }
    
    if (is.numeric(regulators) && max(regulators) > nrow(exprMatrix)) {
      stop("At least one index in regulators exceeds the number of genes.")
    }
  }
  
  if (!is.numeric(nCores) || nCores<1)
  {
    stop("Parameter nCores should be a stricly positive integer.")
  }
  
  if (!is.list(cluster))
  {
    stop("Parameter cluster should be either a cluster configuration (list) or NULL.")
  }
  
  if (!is.logical(progress))
  {
    stop("Parameter progress should be a boolean.")
  }
}