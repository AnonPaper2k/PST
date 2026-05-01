###################################################################################################################################################

# PROJECT: SMMP implementation
# DOC:     Example illustrating the issue with binary covariates held by Party 1
# BY:      X
# DATE:    X
# UPDATE:  --  

# License: https://creativecommons.org/licenses/by-nc-sa/4.0/
# Copyright: X

####################################################################################################################################################

rm(list=ls())

#-------------------------------------------------------------------------------
# Inputs
#-------------------------------------------------------------------------------

# Fix tolerance for recovering binary features
Epsilon <- 1e-10

#-------------------------------------------------------------------------------
# Libraries
#-------------------------------------------------------------------------------

library(RcppAlgos)  # To compute permutation of a given vector

#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------

LoadPathParty2 <- "Outputs/Party2/"

# Party 2 has access to:
n <- read.csv(paste0(LoadPathParty2, "n.csv"))[1,1]
X1X1 <- as.matrix(read.csv(paste0(LoadPathParty2, "X1X1.csv")))

BinaryX1 <- read.csv(paste0(LoadPathParty2, "BinaryX1.csv"))[,1]
Z <- as.matrix(read.csv(paste0(LoadPathParty2, "Z.csv")))

# Extract the number of 1s of the first binary feature of X1
muX1 <- read.csv(paste0(LoadPathParty2, "meansX1.csv"))[1,1]
c <- round(n*muX1,0)

#===============================================================================
# Recovering the first binary feature of X1 from data available at Party 2 after running SMPC
#===============================================================================

# Initialize the list of solutions
ids <- c()
c_ids <- c()

# Centering matrix
Id <- diag(nrow = n)
One <- as.matrix(rep(1, n))
M <- Id - (1/n)*One%*%t(One)

for(j in 1:length(c)){
  
  # Compute number of candidates and batch size to use
  NumberCandidates <- choose(n, c[j])
  BatchSize <- choose(26,13)
  
  # Initialize progress bar 
  pb <- txtProgressBar(min = 0, max = NumberCandidates, style = 3)
  
  # Initialize iterators 
  i <- 0
  Done <- FALSE
  
  while(!Done){
    
    # Select a batch of candidates, and keep those orthogonal to Z, given a tolerence of Epsilon
    p <- permuteGeneral(c(0, 1), freqs = c(n-c[j], c[j]), lower = i*BatchSize+1, upper = min(NumberCandidates,(i+1)*BatchSize), 
                        Parallel = FALSE,
                        FUN = function(p){
                          res <- t(M%*%p)%*%Z
                          return(sum(abs(res))<Epsilon)
                        }, keepResults = TRUE)
    
    # Add new ids to list of solutions, if any
    if(length(which(unlist(p)))>0){
      ids <- c(ids, i*BatchSize + which(unlist(p)))
      c_ids <- c(c_ids, j)
    }
    
    # Update progress bar
    setTxtProgressBar(pb, (i+1)*BatchSize)
    
    # Check if we are 
    if((i+1)*BatchSize>=NumberCandidates){
      Done <- TRUE
    }
    
    # Update iterator
    i <- i + 1
    
    # Clear memory
    gc()
  }
  
  # Update progress bar
  close(pb)
  
}

# Retrieve all solutions
ReconstructedVectors <- matrix(nrow = length(ids), ncol = n)

if(length(ids)>0){
  for(i in 1:length(ids)){
    ReconstructedVectors[i,] <- permuteGeneral(c(0, 1), freqs = c(n-c[c_ids[i]], c[c_ids[i]]), lower = ids[i], upper = ids[i])  
  }  
}

#===============================================================================
# Compare recovered data with true data
#===============================================================================

# Load true data
X1 <- as.matrix(read.csv("Outputs/Party1/X1.csv"))

# Select binary features of X1
X1Bin <- as.matrix(X1[,1])

if(length(ids)==1){
  
  # Counts the number of cells recovered
  n_rec <- sum(t(X1Bin)==ReconstructedVectors)
  
  print(paste0(nrow(ReconstructedVectors), " candidate was found for the first binary column of X1."))
  print(paste0("Expected recovered cells: ", n, "/", n))
  print(paste0("Recovered cells: ", n_rec, "/", n))
  
}else{
  
  # Counts how many cells should be recovered, as several candidates were found
  if(nrow(ReconstructedVectors)>0){
    n_expected <- sum(apply(ReconstructedVectors, 2, function(x) all(x == x[1])))  
  } else{
    n_expected <- 0
  }
  
  
  # Counts the number of cells recovered
  if(nrow(ReconstructedVectors)>0){
    n_rec <- sum(apply(rbind(ReconstructedVectors, X1Bin), 2, function(x) all(x == x[1])))  
  } else{
    n_rec <- 0
  }
  
  
  print(paste0(nrow(ReconstructedVectors), " candidates were found for the first binary column of X1."))
  print(paste0("Expected recovered cells: ", n_expected, "/", n))
  print(paste0("Recovered cells: ", n_rec, "/", n))
}

