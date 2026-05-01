###############################################################################################################################################################

# PROJECT: SMMP implementation
# DOC:     Example illustrating the quantities known by both parties whilst running SPMC
# BY:      X
# DATE:    X
# UPDATE:  --  

# License: https://creativecommons.org/licenses/by-nc-sa/4.0/
# Copyright: X

###############################################################################################################################################################

rm(list=ls())

#-------------------------------------------------------------------------------
# "Save" parameters for each party
#-------------------------------------------------------------------------------

# Should we save all quantities available at Party 1?
# If set to "TRUE", everything will be saved in "Outputs/Party1/"
SaveParty1 <- TRUE
SavePathParty1 <- "Outputs/Party1/"

# Should we save all quantities available at Party 2?
# If set to "TRUE", everything will be saved in "Outputs/Party2/"
SaveParty2 <- TRUE
SavePathParty2 <- "Outputs/Party2/"

#-------------------------------------------------------------------------------
# Libraries
#-------------------------------------------------------------------------------
library(MASS)       # To compute the (left) null space of a matrix
library(Matrix)     # To compute the rank of a matrix

#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------

# Load data from party 1
X1 <- as.matrix(read.csv("Data/DataParty1.csv"))

# Load data from party 2
X2 <- as.matrix(read.csv("Data/DataParty2.csv"))

if(SaveParty1){
  write.csv(X1, file = paste0(SavePathParty1,"X1.csv"), row.names = FALSE)
}

if(SaveParty2){
  write.csv(X2, file = paste0(SavePathParty2,"X2.csv"), row.names = FALSE)
}

#-------------------------------------------------------------------------------
# Shared information
#-------------------------------------------------------------------------------

# Number of rows is known across all parties
n <- nrow(X1)

# It is assumed that parties are willing to share their rank with each other
RankX1 <- rankMatrix(X1)
RankX2 <- rankMatrix(X2)

# It is assumed that parties all know the variables used in the model. As such,
# all parties know which covariates are binary (dummy) variables.
BinaryX1 <- apply(X1, 2, function(col){
  all(col %in% c(0,1))
})
BinaryX2 <- apply(X2, 2, function(col){
  all(col %in% c(0,1))
})

# Center data
Id <- diag(nrow = n)
One <- as.matrix(rep(1, n))
M <- Id - (1/n)*One%*%t(One)
muX1 <- apply(X1, 2, mean)
X1 <- M%*%X1
muX2 <- apply(X2, 2, mean)
X2 <- M%*%X2

# It is assumed that parties will share their local Gram matrix, which is a 
# bloc of the matrix [X^t X]
X1X1 <- t(X1)%*%X1
X2X2 <- t(X2)%*%X2

if(SaveParty1){
  write.csv(X1X1, file = paste0(SavePathParty1,"X1X1.csv"), row.names = FALSE)
  write.csv(X2X2, file = paste0(SavePathParty1,"X2X2.csv"), row.names = FALSE)
  write.csv(BinaryX1, file = paste0(SavePathParty1,"BinaryX1.csv"), row.names = FALSE)
  write.csv(BinaryX2, file = paste0(SavePathParty1,"BinaryX2.csv"), row.names = FALSE)
  write.csv(n, file = paste0(SavePathParty1,"n.csv"), row.names = FALSE)
  write.csv(muX1, file = paste0(SavePathParty1, "meansX1.csv"), row.names = FALSE)
  write.csv(muX2, file = paste0(SavePathParty1, "meansX2.csv"), row.names = FALSE)
}

if(SaveParty2){
  write.csv(X1X1, file = paste0(SavePathParty2,"X1X1.csv"), row.names = FALSE)
  write.csv(X2X2, file = paste0(SavePathParty2,"X2X2.csv"), row.names = FALSE)
  write.csv(BinaryX1, file = paste0(SavePathParty2,"BinaryX1.csv"), row.names = FALSE)
  write.csv(BinaryX2, file = paste0(SavePathParty2,"BinaryX2.csv"), row.names = FALSE)
  write.csv(n, file = paste0(SavePathParty2,"n.csv"), row.names = FALSE)
  write.csv(muX1, file = paste0(SavePathParty2, "meansX1.csv"), row.names = FALSE)
  write.csv(muX2, file = paste0(SavePathParty2, "meansX2.csv"), row.names = FALSE)
}

#===============================================================================
# Party 1: Computes and shares matrix Z
#===============================================================================

# Compute sample size
n <- nrow(X1)

# Compute the left null space of X (orthogonal of the column space of X)
Zfull <- Null(X1)
rankZfull <- n - RankX1

# Choose a subset of Zfull. 
# Compute g according to Karr's equation (8)
g <- round(RankX1/(RankX1+RankX2)*n)

# Choose Z
sampleSize <- min(g, rankZfull)
columns <- sample(1:dim(Zfull)[2], size = sampleSize, replace = FALSE)
Z <- as.matrix(Zfull[, columns])

if(SaveParty1){
  write.csv(g, file = paste0(SavePathParty1,"g.csv"), row.names = FALSE)
  write.csv(Z, file = paste0(SavePathParty1,"Z.csv"), row.names = FALSE)
}

if(SaveParty2){
  write.csv(g, file = paste0(SavePathParty2,"g.csv"), row.names = FALSE)
  write.csv(Z, file = paste0(SavePathParty2,"Z.csv"), row.names = FALSE)
}

#===============================================================================
# Party 2: Computes and shares matrix W
#===============================================================================

# Compute the W matrix
I <- diag(1, nrow(X2))
ZZ <- Z%*%t(Z)
W <- (I - ZZ)%*%X2

if(SaveParty1){
  write.csv(W, file = paste0(SavePathParty1,"W.csv"), row.names = FALSE)
}

if(SaveParty2){
  write.csv(W, file = paste0(SavePathParty2,"W.csv"), row.names = FALSE)
}

#===============================================================================
# Party 1: Computes and shares the secure product t(X1)%*%X2
#===============================================================================

# Compute the matrix product securely
X1_X2 <- t(X1)%*%W

if(SaveParty1){
  write.csv(W, file = paste0(SavePathParty1,"X1X2.csv"), row.names = FALSE)
}

if(SaveParty2){
  write.csv(W, file = paste0(SavePathParty2,"X1X2.csv"), row.names = FALSE)
}

#-------------------------------------------------------------------------------
# Optionnal: Compare SPMC results to matrix multiplication
#-------------------------------------------------------------------------------

# True matrix multiplication
X1X2 <- t(X1)%*%X2

# Compare results
all.equal(X1X2, X1_X2)

#-------------------------------------------------------------------------------
# Clear all objects in environment
#-------------------------------------------------------------------------------

rm(list=ls())
