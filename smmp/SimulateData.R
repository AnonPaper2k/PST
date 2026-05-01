###############################################################################################################################################################

# PROJECT: SMMP implementation
# DOC:     Example illustrating the quantities known by both parties whilst running SPMC
# BY:      X
# DATE:    X
# UPDATE:  X

# License: https://creativecommons.org/licenses/by-nc-sa/4.0/
# Copyright: X

###############################################################################################################################################################

rm(list=ls())

#-------------------------------------------------------------------------------
# Inputs
#-------------------------------------------------------------------------------

# Reproductibility
set.seed(34)

# Sample size
n_common <- 100

# Party 1
c_to_recover <- 3
p_bin_1 <- 2
p_cont_1 <- 2

# Party2
p_bin_2 <-2
p_cont_2 <- 2

#-------------------------------------------------------------------------------
# Function
#-------------------------------------------------------------------------------

simulate_data <- function(pBin, pCont, n, c = NULL){

  # Initialize the binary columns of X
  Xbin <- matrix(data = 0, nrow = n, ncol = pBin)  
  if(pBin>0){
    # If c is provided, the first column of X will have a specific number of 1s
    if(!is.null(c)){
      Xbin[,1] <- sample(x = c(rep(1, c), rep(0, n-c)), size = n, replace = FALSE)
    }else{
      Xbin[,1] <- rbinom(n = n, size = 1, prob = 1/2)
    }
    
    # For all other binary columns, select random values
    if(pBin>1){
      for(p in 2:pBin){
        Xbin[,p] <- rbinom(n = n, size = 1, prob = 1/2)
      }  
    }
    
  }
  
  # Initialize the continous columns of X
  Xcont <- matrix(data = 0, nrow = n, ncol = pCont)
  if(pCont>0){
    # For all other continous columns, select random values
    for(p in 1:pCont){
      Xcont[,p] <- rnorm(n = n)
    }  
  }
  
  return(cbind(Xbin, Xcont))
  
}

#-------------------------------------------------------------------------------
# Simulate data
#-------------------------------------------------------------------------------

# Simulate matrices
X1 <- simulate_data(pBin = p_bin_1, pCont = p_cont_1, c = c_to_recover, n = n_common)
X2 <- simulate_data(pBin = p_bin_2, pCont = p_cont_2, c = 3, n = n_common)

# Save matrices as .csv
write.csv(X1, file = "Data/DataParty1.csv", row.names = FALSE)
write.csv(X2, file = "Data/DataParty2.csv", row.names = FALSE)

#-------------------------------------------------------------------------------
# Clear all objects in environment
#-------------------------------------------------------------------------------

rm(list=ls())
