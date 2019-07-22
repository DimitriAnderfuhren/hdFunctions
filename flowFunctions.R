
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
#flowVS::estParamFlowVS(fs, channels)
#flowVS::transFlowVS(fs, channels, cofactors)

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
  
  
  
##### Visualization
# arcsinh transformation functions

#' asinhTransform applies asinh transformation to the columns of a data frame. It uses the cofactors from "cofactorVector"
#' @param dataFrame A data frame with named columns
#' @param cofactorVector A named integer vector. The names must correspond to the column names of "dataFrame"
#' @return A data frame with transformed values
#' 
asinhTransform = function(dataFrame,cofactorVector){
  # Special case: one vector and single cofactor
  if(is.vector(dataFrame) && length(cofactorVector) == 1){
    t = asinh(dataFrame/cofactorVector)
    return(t)
  }
  
  # Get the names of the cofactors that are above zero
  alteredCofactorNames = names(cofactorValues[cofactorValues>0])
  
  # Check if the alteredCofactorNames are in the columnames of the dataFrame
  # Take out the unfound names
  if(!all(alteredCofactorNames %in% as.character(colnames(dataFrame)))){
    notFound = alteredCofactorNames[!alteredCofactorNames %in% as.character(colnames(dataFrame))]
    warning("Not all cofactors found in the data frame")
    warning(notFound)
    alteredCofactorNames = alteredCofactorNames[!(alteredCofactorNames == notFound)]
  }
  
  # Actual transformation
  tdf = dataFrame
  tdf[,alteredCofactorNames] = asinh(t(t(dataFrame[,alteredCofactorNames])/cofactorValues[alteredCofactorNames]))
  
  return(tdf)
}

#' createCofactorVector creates an empty named vector. The names come from the column names of "dataFrame".
#' @param dataFrame A data frame with named columns
#' @return An empty named vector
createCofactorVector = function(dataFrame, colsToUse = NULL){
  
  if(is.null(colsToUse)){
    colsToUse = 1:ncol(dataFrame)
  }
  
  markerNames = as.character(colnames(dataFrame[,colsToUse]))
  nCol = ncol(dataFrame[,colsToUse])
  cofactorValues = vector(mode = "integer",length = nCol)
  names(cofactorValues) = markerNames
  
  return(cofactorValues)
}

#' interactiveScatterPlot displays a subset of the dataFrame
interactiveScatterPlot = function(dataFrame, markerY, markerX,
                                  asinh_bool = F, histogram_bool = F,
                                  cofactor_int, sampleSize = 5000,
                                  saveValue_bool, transformAll_bool, clearAll_bool,
                                  colsToUse = NULL){
  if(is.null(colsToUse)){
    colsToUse = 1:ncol(dataFrame)
  }
  dataFrame = dataFrame[,colsToUse]
  
  # Store dfs and transform
  untransformed_df = dataFrame
  transformed_df = asinhTransform(dataFrame,cofactorVector = cofactorValues)
  # Take subset of values for performance
  transformedSubset_df = transformed_df[sample(x = 1:nrow(transformed_df),size = sampleSize,replace = F),]
  # X and Y values
  xValues = transformedSubset_df[,markerX]
  yValues = transformedSubset_df[,markerY]
  
  
  # Change xValues: Input from checkbox and slider
  if(asinh_bool == T && cofactor_int>0){
    xValues = asinhTransform(xValues,cofactor_int)
  }
  # Save single cofactor value: Input from button
  if(saveValue_bool == T && asinh_bool == T){
    cofactorValues[markerX] <<- cofactor_int
    print(cofactorValues[markerX])
  }
  # Set current cofactor for all
  if(transformAll_bool == T && asinh_bool == T){
    cofactorValues[1:length(cofactorValues)] <<- cofactor_int
    print(cofactorValues)
  }
  # Clear all
  if(clearAll_bool == T){
    cofactorValues[1:length(cofactorValues)] <<- 0
    print(cofactorValues)
  }
  # Define the label of the yAxis. Add cofactor value to label if there is one.
  yLab = paste(markerY,cofactorValues[[markerY]])
  
  if(cofactor_int>0 && asinh_bool == T){
    xLab = paste(markerX,cofactor_int)
  }else if(cofactorValues[[markerX]]>0){
    xLab = paste(markerX,cofactorValues[[markerX]])
  }else{
    xLab = markerX
  }
  # Choose Color and pch
  if(cofactorValues[[markerX]]>0){
    plotColor = color1
    
  }else{
    plotColor = color2
    
  }
  
  # Plot
  if(histogram_bool == F){
    par(mfrow = c(1,1))
    plot(xValues,yValues,
         xlab = xLab,
         ylab = yLab,
         col = plotColor,
    )
  }else if(histogram_bool == T){
    par(mfrow = c(1,2))
    plot(xValues,yValues,
         xlab = xLab,
         ylab = yLab,
         col = plotColor,
    )
    hist(xValues,
         xlab = xLab,
         main = "",
         col = plotColor)
  }
}
  
  
}