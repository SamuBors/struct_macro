In the folder “codes” you find:
-“graphs.m”, which produces Figure 1, 2 and 3 of the paper (single plots and variance-covariance estimation)
-“VARs_3R_2019.m”, which contains the analysis with 3 volatility regimes and data up to December 2019 (section 4 of the paper)
-“VARs_3R_2023.m”, which contains the analysis with 3 volatility regimes and data up to June 2023 (section 4 of the paper)
-“VARs_4R.m”, which contains the main analysis, with 4 volatility regimes and data up to June 2023 (section 3 of the paper)

In the folder “functions” for both the subfolders you find:
-“ll_LL” which contains the (negative) log-likelihood function to be minimize for the Lanne and Lukethpol (2008) approach
-“ll_Upper_MS”, “ll_Lower_MS” and “ll_Proposal_MS”  contains the (negative) log-likelihood function to be minimize for the three models in table 2 of the paper (4 volatility regimes)
-“ll_Upper_MS_3” and “ll_Lower_MS_3” contains the (negative) log-likelihood function to be minimize for the three models in table 3 and 4 of the paper (3 volatility regimes)
-“SE_sums” and “SE_sums_2019” contain function to obtain the s.e. of the parameters of the matrices in table 2, 3 and 4
-“to_latex” and “to_latex_sym” allow to export a matrix (and its significance) from Matlab to latex
