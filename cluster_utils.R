##########################################
# Function for labeling clusters based on
# centroid MFI and definition data frame.
##########################################

label_clusters <- function(centroids, defs, outliers=c()) {
  not_outliers <- setdiff(seq(nrow(centroids)), outliers)
  th <- apply(centroids[not_outliers,], 2, get_thresholds)
  
  modality <- lapply(colnames(centroids), function(marker) {
    if_else(centroids[,marker] < th[marker], "lo", "hi")
  }) %>% do.call(what=cbind)
  colnames(modality) <- colnames(centroids)
  
  matches <- match_defs(defs, modality)
  phenos <- c(defs$Phenotype, "Other")
  labels <- phenos[matches]
  labels <- make.unique(labels)
  
  return(labels)
}

##########################################
# Divisive hierarchical clustering (diana)
# of the cluster centroids into "hi" and
# "lo" modalities.
# Independently for all markers.
##########################################
get_thresholds <- function(marker)
{
  diana <- cluster::diana(marker)
  cut <- cutree(as.hclust(diana), k = 2)
  
  if (mean(marker[which(cut == 1)]) < 
      mean(marker[which(cut == 2)])) {
    th <- (max(marker[which(cut == 1)]) + min(marker[which(cut == 2)]))/2
  }
  else {
    th <- (max(marker[which(cut == 2)]) + min(marker[which(cut == 1)]))/2
  }
  return(th)
}


#########################################
# For each phenotype, check for clusters
# whose modality matches the definition,
# and label them.
# IMPORTANT: this process is sequential,
# and matches down the list overwrite
# previous ones. Design the phenotype
# spreadsheet with this in mind.
#########################################
match_defs <- function(defs, modality) {
  ind <- integer(nrow(modality))
  
  for (i in seq(nrow(defs))) {
    markers <- names(defs)[which(defs[i,] %in% c("hi", "lo"))]
    match <- apply(modality[,markers,drop=FALSE], 1, function(x) {
      def <- as.character(defs[i,markers])
      chx <- as.character(x)
      all.equal(chx,def)
    })
    sel <- which(match == "TRUE")
    ind[sel] <- i
  }
  
  ind[which(ind==0)] <- nrow(defs)+1
  return(ind)
}

