# Cross-Validation
Analysis tools for Cross-Validation for model discernment using Python

Cross validation is a statistical methodology that involves the division of the data set into independent subsets used for generation of parameter estimation (training) and validation (also known as ‘testing’). In our implementation, for each fitting trial 20% of the experimental data were randomly assigned to the training set, and the remainder were used for validation. In order to overcome the disadvantage of training a subset of the data, five hundred (500) trials were completed for each of the competing models, and the lowest average sum of the squares of the residuals of the validation fit was considered the “best” model for the data.

The file 20mev_3.75.txt contains RK Metric vs. Time data.  Outputs are the parameters and their statistical uncertainties (One Standard Deviation).
