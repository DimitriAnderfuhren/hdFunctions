
# Roadmap

# path of experiment folder

# read .fcs
# read channel
# read metadata

# match information from metadata and channel file with flowSet 
  # phenoData

# save raw experiment

# Maybe flowStats

# arcsinh transformation -> arcsinhTransformed flowSet
flowVS::estParamFlowVS(fs, channels)
flowVS::transFlowVS(fs, channels, cofactors)

  # normalization01 for flowSom -> flowSOM flowSet
  # matrixStats
  normC = function(expr, colsToUse = NULL){
    
    if(is.null(colsToUse)){
      colsToUse = 1:ncol(expr)
    }
    
    Xm = as.matrix(expr[,colsToUse])
    
    rng <- colQuantiles(Xm, probs = c(0.01, 0.99))
    Xm <- t((t(Xm) - rng[, 1]) / (rng[, 2] - rng[, 1]))
    Xm[Xm < 0] <- 0
    Xm[Xm > 1] <- 1
    
    expr[,colsToUse] = Xm
    
    
    return(expr)
  }
  # normTsne for Tsne and UMAP -> vis flowSet
  subsampleV = function(expr, balance = "conditionOrSampleID"){
    return(subExpr)
  } 
  
  normV = function(expr, colsToUse = NULL){
    # Adjust data to start from (approximately) 0
    normV_expr = expr
    q = apply(expr[,colsToUse], 2, function(x) quantile(x, 0.01, names = F))
    normV_expr[,colsToUse] = sweep(normV_expr[,colsToUse], 2, q)
    
    
    # Normalize to have everything from 0 to 1
    per = apply(normV_expr[,colsToUse], 2, function(x) quantile(x, 0.99999, names = F))
    normV_expr[,colsToUse] = t(t(normV_expr[,colsToUse]) / as.numeric(per))
    
    # Check whether you adjusted the range approximately from 0 to 1
    return(normV_expr)
  }

# marker alignment
  # Maybe interactive plot with manipulate
  quickDensityPlotly = function(expr){
    return(interactivePlot)
  }

# (flowStats)
  warpExpression = function(expr,sampleIDs, colsToUse = NULL, out = "expr"){
    if(length(sampleIDs) != nrow(expr)){
      stop("sample_id vector not equal length as expression matrix")
    }
    if(is.null(colsToUse)){
      colsToUse = colnames(expr)
    }
    s = split.data.frame(expr,fname)
    ffs = flowCore::flowSet(lapply(s,createFlowFrame))
    ffs_warped = flowStats::warpSet(ffs, stains = colsToUse)
    
    if(out == "fs"){
      return(ffs_warped)
    }
    if(out == "expr"){
      expr_warped = flowCore::fsApply(ffs_warped,exprs)
      return(expr_warped)
    }
  }

# save transformed experiment

# Maybe create a class for plotting object
# Umap
# tSNE

# FlowSOM

# common functions: 
getClusterForCellID = function(){
  return()
}

createFlowFrame = function(expr, colsToUse = NULL){
  
  if(checkColNames(expr,colsToUse)){
    if(is.null(colsToUse)){
      colsToUse = 1:ncol(expr)
    }
    ff = new("flowFrame",as.matrix(expr[,colsToUse]))
    return(ff)
  }
}

extractClusterPerSample = function(expr, clusterMergings = NULL, clusterName = NULL,sampleIDs = NULL, matrixOnly = T){
  
  mergedExpr = data.frame(expr,cluster = clusterMergings,sample_id = sampleIDs)
  
  frames_list = split(mergedExpr,mergedExpr[["sample_id"]])
  
  cluster_list = lapply(frames_list,FUN = function(x){x[x[["cluster"]] == clusterName,]})
  
  if(matrixOnly == T){
    cluster_list = lapply(cluster_list, FUN = function(x){x[,!colnames(x) %in% c("cluster","sample_id")]})
  }
  
  return(cluster_list)
  
}