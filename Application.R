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
data<-"Data_agriculture_v2.xlsx"
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
Y<-t(rbind(data[, 3:5]))
X_matrix <- t(rbind(data[, 6:8]))

#4. coordinates
coords <- matrix(c(
  0.000,     0.000,
  -51.747,    5.547,
  -415.144,   25.870,
  -662.125,  -55.056,
  -582.202, -105.580,
  -612.205,  -12.316
), ncol = 2, byrow = TRUE)

colnames(coords) <- c("x", "y")
coords

#5. Run Spatio-Temporal regression
source(file.path(base_dir, "codes", "Functions SpatioFile.R"))

ResC<-ML.MatrixRegreSPatioT(Y, X_matrix, coords, precision=0.000001, MaxIter=100, corr="exponential")
ResC
