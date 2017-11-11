# Recursive Inversion Models

This package was made as part of my UW Statistics PhD prelim exam and is designed for performing inference with recursive inversion models (RIMs). For background on recursive inversion models see the original [paper](https://papers.nips.cc/paper/5579-recursive-inversion-models-for-permutations) or my [project report](https://www.stat.washington.edu/~lucaw/assets/revisiting_recursive_inversion_models_for_permutations.pdf). While the underlying C++ code of the package is well documented, the exposed R functions are not (future work). As a stopgap before I get to writing better documentation, the below specifies how a RIM is defined in this package and gives an example demonstrating the package's functionality.

A RIM tree structure is encoded using a matrix of a special form. In particular the matrix contains must be of dimension (number of nodes in RIM) x 5. Each row of the matrix corresponds to a tree node and its position in the list corresponds to the node's numbering. Each vector should have the form c(left child, right child #, is leaf boolean, theta, rank) where:
* left child # = the number correspond to the left child of the node
* right child # = the number correspond to the right child of the node
* is leaf boolean = 1 if the node is a leaf, 0 otherwise
* theta = the theta value corresponding to the node (for internal nodes)
* rank = the rank of the node (for leaves) left/right child #s and theta are ignored if the node is a leaf and rank is ignored if the node is not a leaf. 
Note that the matrix should be ordered topologically so that PARENTS ALWAYS COME BEFORE CHILDREN. As an example the matrix:

```
matrix(c(1, 2, 0, -.1, 0,
         3, 4, 0, 0.8, 0,
         5, 6, 0, 1.6, 0,
         0, 0, 1, 0.0, 0,
         0, 0, 1, 0.0, 1,
         0, 0, 1, 0.0, 2,
         0, 0, 1, 0.0, 3), ncol = 5, byrow = T)       
```

corresponds to a tree with a root node 0, the root having left and right children as 1,2 respectively. Moreover node 1 has children 3,4 which are both leaves and node 3 has children 5,6 which are also both leaves. If traversing the tree in preorder one visits the leaves in the order 3,4,5,6 which corresponds to the ranking 0,1,2,3.

The below will give an example of how to install the RIM package as well as how to create, sample from, and estimate a RIM from data.

```
# Install the RIM package from github, you’ll need to have already
# installed the devtools package from CRAN.
devtools::install_github("Lucaweihs/RIM")
library(RIM)

# Define a RIM tree structure
treeAsMat = matrix(c(
 1, 2, 0, -.1, 0,
 3, 4, 0, .8, 0,
 5, 6, 0, 1.6, 0,
 0, 0, 1, 0, 0,
 0, 0, 1, 0, 1,
 0, 0, 1, 0, 2,
 0, 0, 1, 0, 3),
 ncol = 5, byrow = T)

# Plot the tree structure to make
# sure it is what you expect
plotTreeMatrix(treeAsMat)

# Generate a matrix of sample
# ranks from the above RIM
set.seed(1)
numSamples = 10000
samples = rRIM(numSamples, treeAsMat)

# Create the average discrepancy
# matrix corresponding to the samples
disMat = averageDiscMatrix(samples)

# Up to a constant scalar, what's the 
# log-probability of this data for RIM?
pRIM(treeAsMat, disMat, log.p = TRUE)

# For a given reference ranking, run
# the dynamic programming algorithm
# to estimate the RIM from the data.
# The makeCanonical parameter just makes
# sure all the theta's are positive
estimate1 = structByDP(disMat, refRanking = c(3,2,1,0), makeCanonical = F)
plotTreeMatrix(estimate1) # Plot the estimated tree

# Run the simulated annealing search
# from the paper to get a, possibly different,
# estimated tree
estimate2 = SASearch(disMat, refRanking = c(3,2,1,0), 
                    inverseTemp = 0.1, maxIter = 1000,
                    makeCanonical = T)
plotTreeMatrix(estimate2) # Plot the estimated tree

# Say you want to, for a given RIM, update it's
# theta parameters so that they are the MLE values
# for a given data set. You can also do that with the
# function
treeWithMLEUpdatedThetas = thetaMLERIM(treeAsMat, disMat)
plotTreeMatrix(estimate2) # Plot the estimated tree
```
