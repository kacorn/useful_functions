EdgeLabelTips <- function(tree){
  is_tip <- tree$edge[,2] <= length(tree$tip.label) #since the tips are the first to be labeled
  edges_and_tips <- as.data.frame(which(is_tip == TRUE))
  edges_and_tips[,2] <- tree$tip.label
  colnames(edges_and_tips) <- c("edge_label", "tip_label")
  return(edges_and_tips)
}

#many thanks to this guy: https://stackoverflow.com/questions/34364660/how-to-get-correct-order-of-tip-labels-in-ape-after-calling-ladderize-function
