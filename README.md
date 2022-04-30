# Application of Matrix Exponential Spatial Specification in Hedonic Price Models

Conventional hedonic models do not take into account the spatial autocorrelation, when applied on spatial datai, which can be a problem. Hence, they were extended to take it into account <i>(Anselin (1998))</i>, mainly two models were specified - spatial lag model and spatial error model. However, these models have some problems of their own, i.e., computation of the inverse of a matrix which might not exists always, and requiring a lot of computational resources if the data is quite large.Hence, <i>Lesage and Pace (2003)</i> created Matrix Exponential Spatial Specification which solves these problems. Therefore, to study sales prices of houses in Baltimore we use a MESS specification for a Hedonic Price Model.

### Note
Most of the matlab code is written by Lesage and Pace. I have only written the following two files:-

1. [MESS/BRAR_code.m](https://github.com/shawnbrar/MESS_Hedonic_Price/blob/master/MESS/BRAR_code.m): code for model results
2. [MESS/obj2mat.R](https://github.com/shawnbrar/MESS_Hedonic_Price/blob/master/MESS/obj2mat.R): code for preparation of data

The following is my project report:-
[BRAR_Sudhakar_DASEE_Project.pdf](https://github.com/shawnbrar/MESS_Hedonic_Price/blob/master/BRAR_Sudhakar_DASEE_Project.pdf): project report
