# Sparse-Functional-Linear-Discriminant-Analysis-SFLDA-
This repository contains R implementation of SFLDA based on the papers

1. Sparse Functional Linear Discriminant Aanalysis 
2. Interpretable Classification for Multivariate Gait Analysis of Cerebral Palsy

Source Models.R file and use SFLDA(data, y, tau, lambda) function to train the model. Data should be a list of n x p matrices (list of length 1 for univariate SFLDA) and y should be a vector of classes where order matches with the data. tau and lambda are parameters which control somoothness and sparsity respectively.  After training, use predict_class(model, x_test) function for class assignment. Here, the model is a list object obtained from the SFLDA function.
