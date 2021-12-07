
Doubly robust inference for non-probability surveys
author: Ali Rafei
E-mail: arafei@umich.edu
Date  : 12-06-2021
Paper : https://arxiv.org/abs/2101.07456

The following function provides doubly robust (DR) inference for a non-probability survey (S_A) where a parallel probability survey (S_R) is available as the reference survey. In both samples, it is assumed a common set of auxiliary variables, X, are observed, but the outcome, Y, is supposed to be observed only in S_A. On the other hand, S_R comes with a set of sampling weights, while pseudo-weights in S_A are unoberved.

The function AIPW_KH estimates the finite population mean for a continuous or binary outcome [family=c('gaussian', 'binomial')] based on the Augmented Inverse Propensity Weighting (AIPW) idea proposed by Kim & Haziza (2014) that yields DR point and variance estimates by solving the estimating equations associated with the proponsity model and outcome models.

Major limitations:
1. The same set of X should be used in propensity and outcome models.
2. Unique solutions in the joint EE are not guaranteed.
3. The sampling weights of S_R are calculable for units of S_A.

The function uses three different propensity approaches:
1. PMLE proposed by Chen et al (2019)
2. PAPP proposed by Rafei et al (2020)
3. IPSW proposed by Wang et al (2020)
