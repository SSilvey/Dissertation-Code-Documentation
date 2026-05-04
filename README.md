**Code Block: "Data and Functions.R:"**
This code introduces various functions used to generate Aim 1 results, as well as loading and pre-processing of relevant data sources. 

**An overview of each function is below:**

**Subsample_DF_random:**
Creates multiple random subsamples of row indices (based on row names) from a data frame.
Each subsample contains n unique rows sampled without replacement. The function
returns a list of length 10, where each element is a character vector of sampled row names.

**XGBoost, RF, Lasso, NN Curve Optimal:**
Performs 5-fold cross validation on a dataset of choice, returning AUC and standard error.

**XGBoost, RF, Lasso, NN Curve Evaluate:**
Performs 5-fold cross validation on a dataset and given subsample size of choice, returning AUC at size n. 

**XGBoost, RF, Lasso, NN Curve Random:**
Generates learning curve data over a range of sample sizes, using row indices generates from the Subsample_DF_random function.
Returns a data frame of sample sizes and AUCs at size n over the sequence specified in each function.

**plot_curve_full**:
Allows for plotting of fitted learning curve once generated. The input is the raw data, and users must specify a method of curve fitting : either Power law (estimated c/fixed c) or log-linear model.

**LASSO_Optimal_CoreFeatures**:
Code to extract core linear features from an entire dataset. Returns the raw number of nonzero coefficients from LASSO model.

**Get_N, Get_N_Fixed, Get_N_LOG**:
Various functions used to extract the optimal (tolerance = 0.02) sample size from a fitted learning curve. Which function to use depends on the curve fitting method of choice.

**Get_Fitted_Values, Get_Fitted_Fixed, Get_Fitted_LOG**:
Returns the raw fitted values when applying NLS or linear regression to fit a learning curve. Which function to use depends on the curve fitting method of choice.

**Code Block: "Sample Size ML XGBoost.R:"**
This code applies functions introducted above to the XGBoost algorithm.

**Code Block: "Sample Size ML RandomForest.R:"**
This code applies functions introducted above to the Random Forest algorithm.

**Code Block: "Sample Size ML LASSO Regression.R:"**
This code applies functions introducted above to the logistic regression algorithm.

**Code Block: "Sample Size ML Neural Network.R:"**
This code applies functions introducted above to the Neural Network algorithm.

**Code Block: "Data and Functions Aim 2.R:"**
This code introduces various functions used to generate Aim 2 results, as well as loading and pre-processing of relevant data sources. 

**An overview of each function is below:**

**Simulate_OMICs_Data_BNG_Only**: 
Code used to simulate bayesian network generated data from original sources. 
Parameter "Keep" specifies the DEGs to be modelled.
Parameter "Transform" specifies whether the input data should be log2-transformed prior to simulation.
Parameter "fast" = T will use the Tabu algorithm, otherwise RS2max will be used. (Note the previous 2 parameteters will always be T for this dissertation).
Parameter "Case_Rows" specifies the row indices of the minority class.
Parameter "p" specifies the number of predictors to be modelled.
Parameter "n" specifies the simulated data sample size.

**AvgReadCount**:
Generates average median read count from a dataset. Input should be raw read counts.

**GetDispersion / GetDispersion_Sim**:
Generates median dispersion from a dataset. Input should be raw read counts. UseWhole indicates whether the entire dataset should be used.

**GetCorrelation**
Generates median correlation from a dataset. Input should be normalized read counts. 

**FindModel_XGB/RF/NN**:
Enumerates predictors and summarizes AIC, R2 over various combinations. Used to find best-fitting models.

**The remainder of the functions are various iterations of those introduced in the Aim 1 code and will not be further discussed.**















