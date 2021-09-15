EdgeLabelTips <- function(tree){
  is_tip <- tree$edge[,2] <= length(tree$tip.label) #since the tips are the first to be labeled the first n-your-number-of-tips of edge lengths will be the tip labels
  edges_and_tips <- as.data.frame(which(is_tip == TRUE)) #pull out the tips
  edges_and_tips[,2] <- test_tree$tip.label #add tip labels to a new column
  colnames(edges_and_tips) <- c("edge_label", "species_name") #make names prettier
  return(edges_and_tips)
}

#many thanks to this guy: https://stackoverflow.com/questions/34364660/how-to-get-correct-order-of-tip-labels-in-ape-after-calling-ladderize-function
