# Matrix-Variate Regression Model for Multivariate Spatio-Temporal Data

Carlos A. R. Diniz <sup>1</sup>, Victor E. Lachos<sup>2</sup> and  Victor H. Lachos<sup>3</sup>
<sup>1</sup> Statistics Department - Federal University of Sao Carlos, São Carlos, Brazil

<sup>2</sup> Statistics Department - Federal University of Sao Carlos and University of São Paulo, São Carlos, Brazil

<sup>3</sup> Statistics Department - University of Connecticut, Storrs, CT-06269, USA

## Previous installations

- Running the scripts requires the following R packages:matrixNormal, MixMatrix, mvtnorm, MomTrunc, mnormt, and Matrix related to multivariate statistical distributions.

## File description:

The project includes the following files and folders:

- Main.sim.R: This script processes the data simulation and save the results in output_dir. It utilizes one primary function Functions SpatioFile.R.

- Functions SpatioFile.R: contained in the carpet codes, it has the necessary functions to run the model starting for AR(1) correlation, correlation functions,  response sample generator and the Maximum Likelihood for the completed data. Modifications can be made based on the type of dataset to be used, mostly in initial values of parameters, range for rho and phi parameters, number of iterations and numbers of parameters considered in criteria.

-errors_Frobenius.R: It runs the frobenius errors for the Matrix beta and Sigma of the results obtained in the carpet results_test.

- Application.R: it runs the Agriculture dataset. It utilizes one primary load dataset data and changes can be made as well as to the coordinates of the places.

