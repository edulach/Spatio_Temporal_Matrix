# Matrix-Variate Regression Model for Multivariate Spatio-Temporal Data

Carlos A. R. Diniz <sup>1</sup>, Victor E. Lachos<sup>2</sup> and  Victor H. Lachos<sup>3</sup>
<sup>1</sup> Statistics Department - Federal University of Sao Carlos, São Carlos, Brazil

<sup>2</sup> Statistics Department - Federal University of Sao Carlos and University of São Paulo, São Carlos, Brazil

<sup>3</sup> Statistics Department - University of Connecticut, Storrs, CT-06269, USA

## Previous installations

- Running the scripts requires the following R packages:matrixNormal, MixMatrix, mvtnorm, MomTrunc, mnormt, and Matrix related to multivariate statistical distributions.

## File description:

- **Main.sim.R**: This script handles data simulation and saves the results in the `output_dir` folder. It relies primarily on functions defined in `Functions_SpatioFile.R`.

- **Functions_SpatioFile.R** (located in the `codes` folder): Contains all necessary functions to run the model, including:
  - Generation of AR(1) correlation structures,
  - Implementation of spatial correlation functions,
  - Simulation of response samples,
  - Maximum likelihood estimation for complete data.  
  Modifications can be made depending on the dataset, particularly regarding:
  - Initial parameter values,
  - Search ranges for the ρ (rho) and φ (phi) parameters,
  - Number of iterations,
  - Number of parameters used in model selection criteria.

- **errors_Frobenius.R**: Computes Frobenius norm errors for the estimated coefficient matrix (Beta) and covariance matrix (Sigma) using results stored in the `results_test` folder.

- **Application.R**: Applies the methodology to an agricultural dataset. It loads the dataset and uses spatial coordinates of the locations. Both the data and coordinates can be updated or customized as needed.


