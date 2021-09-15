EdgeLabelMyTips <- function(tree){
  is_tip <- tree$edge[,2] <= length(tree$tip.label) #since the tips are the first to be labeled
  my_edges_and_tips <- as.data.frame(which(is_tip == TRUE))
  my_edges_and_tips[,2] <- test_tree$tip.label
  colnames(my_edges_and_tips) <- c("edges", "species_name")
  return(my_edges_and_tips)
}

#many thanks to this guy: https://stackoverflow.com/questions/34364660/how-to-get-correct-order-of-tip-labels-in-ape-after-calling-ladderize-function
