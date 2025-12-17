#### Set the path of files ####
base_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(base_dir)


#1. Required packages 
library(matrixNormal)
library(MixMatrix)
library(mvtnorm)
library(MomTrunc)
library(mnormt)
library(Matrix)
library(readxl)

# 2. Load dataset (automatic file detection if only 1 xlsx exists)
data<-"Dados_Completo_Pib_e_Covarivais.xlsx"
data_path <- file.path(base_dir, "dataset", data)

if(length(data_path) > 1) {
  stop("More than one dataset found. Please specify which one to use.")
}
if(length(data_path) == 0) {
  stop("No .xlsx file found in /dataset folder.")
}

data <- read_excel(data_path)
message("\nData loaded: ", basename(data_path))
print(summary(data))

#3. Matrix X and Y
Y<-t(rbind(data[, 3:6]))
X <- t(rbind(data[, 7:9]))

#4. transformations
Y <- log(Y)
X[1, ] <- log(X[1, ])# transformation log for population variable
#4. coordinates
coords <- matrix(c(
  0.00,     0.00,
  -141.25, -100.08,
  -874.23,   20.02,
  -719.07, -518.17,
  -2127.25, 646.04,
  -335.99, 1759.10,
  -1293.68, 1409.95,
  -60.99, 1592.31,
  -1708.86, 780.59,
  -1363.24, 2068.23,
  -46.01, 619.36,
  1305.45, 679.40,
  1007.98, 312.46,
  1004.77, 1341.01,
  388.43, 1473.33,
  1398.55, 962.95,
  1396.41, 859.54,
  548.93, 1188.67,
  1361.09, 1109.73,
  1162.07, 541.52,
  812.16, -504.82,
  426.95, -460.35,
  505.06, -791.71,
  138.04, -863.98,
  -143.39, -1073.03,
  -353.11, -1584.53,
  -66.34, -1314.32
), ncol = 2, byrow = TRUE)

colnames(coords) <- c("x", "y")
coords

#5. Run Spatio-Temporal regression
source(file.path(base_dir, "codes", "Functions SpatioFile.R"))

ResC<-ML.MatrixRegreSPatioT(Y, X, coords, precision=1e-12, MaxIter=100, corr="exponential")
ResC
