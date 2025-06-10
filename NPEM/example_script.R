### NPEM Example Script 
###
rm(list=ls())
## Load required packages
library(np)
library(outliers)
source("/home/haoyi/microbiome mediation analysis/NPEM/NPEMv1.R")
## Load in example data
X <- read.csv("/home/haoyi/microbiome mediation analysis/NPEM/geneExp.csv",row.names=1)
M <- read.csv("/home/haoyi/microbiome mediation analysis/NPEM/OTUabnd.csv",row.names=1)
Y <- read.csv("/home/haoyi/microbiome mediation analysis/NPEM/Diagnosis.csv",row.names=1)
class(X)
## Put data in correct format
M <- apply(M,c(1,2),function(s) log(s+1))
Y$Diagnosis <- as.factor(Y$Diagnosis)
X <- X[, 1, drop = FALSE]

## Run NPEM with Univariate Mediator
UV.example <- NPEM(X,M,Y,method="UV") 
UV.example$mediation.p

## Run NPEM with Bivariate Mediator with test after iteration
BVS.example <- NPEM(X,M,Y,method="BV")
BVS.example$mediation.p

