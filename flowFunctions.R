# exampleData
example = list()
example$nVars = 12
example$nObs = 10000
example$nSamples = 10
example$expr = matrix(rnorm(example$nObs*example$nVars),ncol = example$nVars)
example$MarkerNames = c("CD4", "TIGIT", "IFNg", "CD8", "HLADR", "IL10", "CD56","CD45RA", "IL17", "CCR7", "IL2", "TNFa")
example$sampleIDs = rep(letters[1:example$nSamples], rep(example$nObs/example$nSamples,example$nSamples))
example$clustermergings = rep(letters[1:example$nSamples], rep(example$nObs/example$nSamples,example$nSamples))
colnames(example$expr) = example$MarkerNames


# Constants
minCofactor = 0
maxCofactor = 2000
color1 = "cadetblue"
color2 = "chocolate"
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

balancedSubsample = function(expr,sampleIDs, nTot = 1000){
  
  
  sampleIDs = as.factor(sampleIDs)
  
  all_levels = levels(sampleIDs)
  nLevels = length(all_levels)
  nToSample = as.integer(nTot/nLevels)
  
  enoughCells = table(sampleIDs) > nToSample
  
  if(!all(enoughCells)){
    print(names(table(sampleIDs)[!enoughCells]))
    stop("Not enough cells")
  }
  
  s = split.data.frame(expr,sampleIDs)
  l = lapply(s,FUN = function(x){x[sample(1:nrow(x),size = nToSample),]})
  
  expr_out = do.call(rbind,l)
  
  return(expr_out)
}


  
normV = function(expr, colsToUse = NULL){
    # Adjust data to start from (approximately) 0
    if(is.null(colsToUse)){
      colsToUse = 1:ncol(expr)
    }
  
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
    s = split.data.frame(expr,sampleIDs)
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
  colNames = colnames(expr)
  
  mergedExpr = data.frame(expr,cluster = clusterMergings,sample_id = sampleIDs)
  
  frames_list = split(mergedExpr,mergedExpr[["sample_id"]])
  
  cluster_list = lapply(frames_list,FUN = function(x){x[x[["cluster"]] %in% clusterName,]})
  
  if(matrixOnly == T){
    cluster_list = lapply(cluster_list, FUN = function(x){as.matrix(x[,!colnames(x) %in% c("cluster","sample_id")])})
  }
  
  
  return(cluster_list)
}

#' needs flowCore package
readFcsFiles = function(folderPath_string, pattern = "*.fcs$"){
  allFcsNames = list.files(path = folderPath_string,full.names = T, pattern = pattern)
  print(allFcsNames)
  fcsFiles = flowCore::read.flowSet(files = allFcsNames)
  return(fcsFiles)
}

checkColNames = function(expr, colsToCheck){
  colsThereBool = colsToCheck %in% colnames(expr)
  if(all(colsThereBool)){
    return(T)
  }else{
    print("Columns not found in dataframe:")
    print(colsToCheck[!colsThereBool])
    return(F)
  }
}

createColors = function(){
  info = brewer.pal.info
  set_names = rownames(info[info$category == "qual",])
  all_colors = c()
  for(i in 1:length(set_names)){
    all_colors = c(all_colors,brewer.pal(8,set_names[i]) )
  }
  return(all_colors)
}

factorsToColor = function(vector, colors){
  n = nlevels(as.factor(vector))
  color = colors[1:n]
  names(color) = levels(as.factor(vector))
  return(color)
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
  
  sub_inds = sample(x = 1:nrow(transformed_df),size = sampleSize,replace = F)
  
  transformedSubset_df = transformed_df[sub_inds,]
  untransformedSubset_df = untransformed_df[sub_inds,]
  # X and Y values
  xValues = untransformedSubset_df[,markerX]
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
         col = plotColor
    )
    upperB = quantile(xValues,probs = 0.99)
    lowerB = quantile(xValues,probs = 0.01)
    midB = (upperB + lowerB)/2
    abline(v = upperB, lty = 2)
    abline(v = lowerB, lty = 2)
    abline(v = midB, lty = 2)
    
    
    
  }else if(histogram_bool == T){
    par(mfrow = c(1,2))
    plot(xValues,yValues,
         xlab = xLab,
         ylab = yLab,
         col = plotColor
    )
    hist(xValues,
         xlab = xLab,
         main = "",
         col = plotColor)
  }
}
  

plot_clustering_heatmap_wrapper <- function(fcs, expr01, 
                                            cell_clustering, color_clusters, cluster_merging = NULL) {
  
  # Calculate the median expression
  expr_median <- data.frame(fcs, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarize_each(funs(median))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarize_each(funs(median))
  
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  
  # This clustering is based on the markers that were used for the main clustering
  d <- dist(expr_median[, colnames(fcs)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  
  labels_row <- paste0(rownames(expr_heat), " (",
                       round(clustering_table / sum(clustering_table) * 100, 2), "%)")
  labels_col <- colnames(expr_heat)
  
  # Row annotation for the heatmap
  annotation_row <- data.frame(cluster = factor(expr01_median$cell_clustering)) 
  rownames(annotation_row) <- rownames(expr_heat)
  
  color_clusters <- color_clusters[1:nlevels(annotation_row$cluster)]
  names(color_clusters) <- levels(annotation_row$cluster)
  annotation_colors <- list(cluster = color_clusters)
  annotation_legend <- FALSE
  
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$cluster_merging <- cluster_merging$new_cluster
    color_clusters <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters) <- levels(cluster_merging$new_cluster)
    annotation_colors$cluster_merging <- color_clusters
    annotation_legend <- TRUE
  }  
  
  # Colors for the heatmap
  color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  pheatmap(expr_heat, color = color,
           cluster_cols = FALSE, cluster_rows = cluster_rows,
           labels_col = labels_col, labels_row = labels_row,
           display_numbers = TRUE, number_color = "black",
           fontsize = 8, fontsize_number = 0.001,
           annotation_row = annotation_row, annotation_colors = annotation_colors,
           annotation_legend = annotation_legend)
}
  
