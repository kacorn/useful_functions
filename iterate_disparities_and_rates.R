#################################################################################
# relatively quick and repeatable iteration of multivariate & univariate        #
# brownian rate & disparity results in geomorph, including collating results &  # 
# gathering pairwise p values                                                   #
#################################################################################

#Written and last checked April 2024 by Katherine Corn @kacorn on github
#with package versions ape v 5.8, phytools 2.1-1, reshape2 1.4.4, 
#dplyr 1.1.4, tibble 3.2.1 and geomorph 4.0.7

#If there have been major updates to any of these packages, especially geomorph or dplyr.
# use this with extreme caution and double check every result is going to the right place
# I am not liable if one of these packages changes the way it structures results and 
# consequently mismatches your result...so be cautious!! :)

require(ape) #handling a tree
require(phytools) #just using this for a dataset
require(geomorph) #running 
require(reshape2) 
require(dplyr)
require(tibble)

#we're going to use a dataset from phytools
data("sunfish.data")
data("sunfish.tree")

#I assume your data is ordered properly and has a column with the species names
sunfish.data <- sunfish.data %>%
  mutate(tree_name = rownames(.)) %>%
  arrange(match(tree_name, sunfish.tree$tip.label))


###########################################################
# RUN UNIVARIATE DISPARITY FOR MULTIPLE CONTINUOUS TRAITS #
###########################################################

#let's run through disparity and pull out all the important bits

#set up data here
continuous_traits <- sunfish.data %>% 
  dplyr::select(tree_name, gape.width, buccal.length) #select your contx traits 
  #you want here and order them so the species names are first. can be a column named whatever
discrete_regimes <- sunfish.data %>% 
  dplyr::select(tree_name, feeding.mode) #select your discrete trait here 
 #and order them so the sp names are first


discrete <- as.vector(t(discrete_regimes[,2]))
names(discrete) <- t(discrete_regimes[,1])
#set up continuous data as a named data_frame
first_title <- names(continuous_traits[1]) #set of vector to exclude species names column
continuous <- as.data.frame(continuous_traits %>%
                              dplyr::select(everything(), -one_of(first_title)))
rownames(continuous) <- t(continuous_traits[,1]) #adding rownames back in

groups <- levels(as.factor(discrete))
trait_vec <- names(continuous)
#great! our data is all set up and ready to go

#set up data frame to catch disparity results
disparity_results <- data.frame(matrix(ncol = ncol(continuous), nrow = length(groups)))
colnames(disparity_results) <- trait_vec
rownames(disparity_results) <- groups
disparity_results$discrete <- groups

#set up data frame to catch pairwise comparisons significance results. number of pairwise comparisons formula = k*(k-1)/2
num_comparisons <- length(groups)*(length(groups)-1)/2
disparity_sig <- data.frame(matrix(ncol = length(groups) + 2,nrow = 0))
colnames(disparity_sig) <- c("comparison_group", groups, "trait")

sig_res <- data.frame(matrix(ncol = length(groups) + 2, nrow = length(groups)*ncol(continuous)))
colnames(sig_res) <- c("comparison_group", groups, "trait")
sig_res$comparison_group <- rep(groups)

i = 1 #iterates through continuous traits
row = 1 #iterate through rows of pairwise comparisons p-values

#sequential_disparity
stopifnot(!is.na(continuous[,i]))
stopifnot(!is.na(discrete))

for(i in 1:ncol(continuous)){
  
  j = 1 #iterates through groups within disparity analysis
  
  print(paste("Now running disparity for", trait_vec[i]))
  
  run_disparity <- morphol.disparity(f1 = continuous[,i] ~ discrete, groups = ~ discrete, partial = FALSE, 
                                     iter = 10000, seed = NULL, print.progress = F)
  
  #get pairwise disparities
  for (j in 1:length(groups)){
    disparity_results[j,i] <- run_disparity$Procrustes.var[j]
    j <- j + 1
  }
  
  #pairwise p-values
  sig_results <- tibble(as.data.frame(run_disparity$PV.dist.Pval) %>%
                          rownames_to_column(var = "comparison_group")) %>%
    mutate(trait = rep(trait_vec[i]))

  disparity_sig <- rbind(disparity_sig, sig_results)
  
  i <- i + 1
  
}

#returns 2 sets of results
#1) disparity_results, with disparities
#2) disparity_sig, which also has pairwise p-values but formatted differently in case you want that

#clean up disparity_results a little bit before putting it in the list
clean_disparities <- disparity_results %>%
  rename(diet = discrete) %>% #or whatever your group is
  dplyr::select(diet, gape.width:buccal.length)

#save disparity results as a list
univariate_disparity_results <- list(clean_disparities, disparity_sig)
names(univariate_disparity_results) <- c("disparities","pairwise_p_values")

#examples of how to extract the table for premax protrusion pairwise p-values in base r
pairwise <- univariate_disparity_results$pairwise_p_values
pairwise[pairwise$trait == "gape.width", ]


#############################
# RUN UNIVARIATE RATES     #
############################

#set up data here
continuous_traits <- sunfish.data %>% 
  dplyr::select(tree_name, gape.width, buccal.length) #select your contx traits 
#you want here and order them so the species names are first. can be a column named whatever
discrete_regimes <- sunfish.data %>% 
  dplyr::select(tree_name, feeding.mode) #select your discrete trait here 
#and order them so the sp names are first

discrete <- as.vector(t(discrete_regimes[,2]))
names(discrete) <- t(discrete_regimes[,1])
#set up continuous data as a named data_frame
first_title <- names(continuous_traits)[1] #set of vector to exclude species names column
continuous <- as.data.frame(continuous_traits %>%
                              dplyr::select(everything(), -one_of(first_title)))
rownames(continuous) <- t(continuous_traits[,1])

groups <- levels(as.factor(discrete))
trait_vec <- names(continuous)
#great! our data is all set up and ready to go

#set up some data frames to catch our data
rate_results <- data.frame(matrix(ncol = ncol(continuous), nrow = length(groups)))
colnames(rate_results) <- trait_vec
rownames(rate_results) <- groups
rate_results$discrete <- groups

#set up data frame to catch pairwise comparisons significance results. 
num_comparisons <- length(groups)*(length(groups)-1)/2
rate_sig <- data.frame(matrix(ncol = length(groups) + 2,nrow = 0))
colnames(rate_sig) <- c("comparison_group", groups, "trait")

i = 1 #iterates through continuous traits

#sequential_rates
print("starting to run rates")
for(i in 1:ncol(continuous)){
  
  j = 1 #iterates through groups within rate analysis
  
  cont_trait <- continuous[,i]
  names(cont_trait) <- rownames(continuous)
  
  run_rates <- compare.evol.rates(A = cont_trait, phy = sunfish.tree, #YOUR TREE HERE
                                  gp = discrete, iter = 10000, 
                                  method = "permutation", print.progress = T)
  
  rate_sig_res <- tibble(melt(as.matrix(run_rates$pairwise.pvalue), 
                              varnames = c("discrete_1","discrete_2"))) %>%
    mutate(trait = rep(trait_vec[i])) 
  
  for (j in 1:length(groups)){
    rate_results[j,i] <- run_rates$sigma.d.gp[j]
    j <- j + 1
  }
  
  rate_sig <- rbind(rate_sig, rate_sig_res)
}

#returns 2 sets of results
#1) rate_results, with rates
#2) rate_sig, which has pairwise p-values

#clean up and save results
rate_sig_clean <- as.data.frame(rate_sig %>%
                                  rename(pairwise_p = value))

rate_results_clean <- rate_results %>%
  rename(feeding = discrete) %>%
  dplyr::select(feeding, gape.width:buccal.length)

univariate_rate_results <- list(rate_results_clean, rate_sig_clean)
names(univariate_rate_results) <- c("rates","pairwise_p_values")


#################################################
# RUN MULTIVARIATE DISPARITIES & RATES          #
#################################################

#run kinematic rates and disparities
#set up data
continuous_traits <- as.matrix(sunfish.data %>% 
  dplyr::select(gape.width, buccal.length)) #select your contx traits and make a matrix
#you want here and remove the species name column
rownames(continuous_traits) <- sunfish.data$tree_name

discrete_regimes <- as.vector(sunfish.data$feeding.mode) #select your discrete trait here 
names(discrete_regimes) <- sunfish.data$tree_name

#run
multivariate_rates <- compare.evol.rates(A = continuous_traits, 
                                      phy = sunfish.tree, method = "permutation", 
                                      gp = discrete_regimes, iter = 10000)

multivariate_disparities <- morphol.disparity(f1 = continuous_traits ~ discrete_regimes, 
                                           groups = ~ discrete_regimes, partial = FALSE, 
                                           iter = 10000, seed = NULL, print.progress = TRUE)

#also works for super high dimensional stuff like landmark pcs


#################################################
# COLLATE MULTIVARIATE DISPARITIES              #
#################################################

extract_disp <- function(groups, disparity_res){
  res <- data.frame(row.names = groups)
  i = 1
  for (i in 1:length(groups)){
    res$groups[i] <- groups[i]
    res$disparity[i] <- disparity_res$Procrustes.var[i]
  }
  return(tibble(res))
}

these_groups <- groups #just want a vector of your groups here

#put together disparity results

disp_res <- extract_disp(groups = these_groups, disparity_res = multivariate_disparities) 

#get disparities pairwise p-values into a neat data frame
disparities_pvalues <- melt(multivariate_disparities$PV.dist.Pval) %>%
  rename(group_1 = Var1, group_2 = Var2, pairwise_disparity_p = value)

#integrate it into a list
all_multivariate_disparity_results <- list(disp_res, disparities_pvalues)
names(all_multivariate_disparity_results) <- c("disparities","pairwise_p_values")

#################################################
# COLLATE MULTIVARIATE RATES                    #
#################################################

#combine all these data

extract_rates <- function(groups, rate_res){
  res <- data.frame(row.names = groups)
  i = 1
  for (i in 1:length(groups)){
    res$group[i] <- groups[i]
    res$rate[i] <- rate_res$sigma.d.gp[i]
  }
  return(tibble(res))
}

these_groups <- groups

rate_results <- extract_rates(groups = groups, rate_res = multivariate_rates)

#get rate p-values into a single data frame

rate_pvalues <- melt(as.matrix(multivariate_rates$pairwise.pvalue)) %>%
  rename(group_1 = Var1, group_2 = Var2, pairwise_rate_p = value)

#collect overall statistics about multivariate disparities
rate_overall_statistics <- c(multivariate_rates$sigma.d.ratio, 
                                                multivariate_rates$P.value, 
                                                multivariate_rates$Z)
names(rate_overall_statistics) = c("overall_sigma2_ratio", "overall_pvalue", "Z_score")

#integrate it into a list
all_multivariate_rate_results <- list(rate_results, rate_pvalues, rate_overall_statistics)
names(all_multivariate_rate_results) <- c("rates","pairwise_p_values", "overall_statistics")

# hope this worked for you!
# heads up, the pairwise rate labeling will be a little wacky if you only have 2 discrete groups

