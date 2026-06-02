library(matrixNormal)
library(MixMatrix)
library(mvtnorm)
library(MomTrunc)
library(mnormt)
library(Matrix)
library(readxl)
library(data.table)
library(dplyr) 
library(tidyr)
library(stringr)

# Defina o caminho base uma vez
base_path <- ""

setwd(base_path)

# Source das funções
source(file.path(base_path, "Functions SpatioFile_b2.R"))

# Leitura dos dados de resposta
dados <- read.csv(file.path(base_path, "dados_pib-2012_2021.csv"))

# Leitura das covariáveis
pnad <- as.data.frame(read_excel(
  file.path(base_path, "pnad_2012_2021.xlsx"),
  sheet = "serie_longa"
))



names(dados)

# Não é necessário filtrar anos nem selecionar colunas —
# a nova base já contém apenas 2012-2021 e as colunas de interesse
dados_filtrados <- dados[, c("ano", "sigla_uf", "va_agropecuaria", 
                              "va_industria", "va_servicos", "va_adespss")]

UF_order <- c("DF", "GO", "MT", "MS", "AC", "AP", "AM", "PA", "RO", "RR", "TO",
              "AL", "BA", "CE", "MA", "PB", "PE", "PI", "RN", "SE", "ES", "MG",
              "RJ", "SP", "PR", "RS", "SC")

dados_filtrados$sigla_uf <- factor(dados_filtrados$sigla_uf, levels = UF_order)



# =============================================================
# Leitura da nova base de covariáveis
# =============================================================


# Ordenar por UF (na ordem do modelo) e depois por ano
pnad <- pnad %>%
  mutate(sigla = factor(sigla, levels = UF_order)) %>%
  arrange(sigla, ano)

pnad <- as.data.frame(pnad)

# =============================================================
# Ordenar dados de resposta da mesma forma
# =============================================================

dados_filtrados <- dados_filtrados[order(dados_filtrados$sigla_uf, dados_filtrados$ano), ]

# =============================================================
# Matriz resposta Y (log-escala, em milhões)
# =============================================================

dados_y <- log(rbind(
  dados_filtrados$va_agropecuaria,
  dados_filtrados$va_industria,
  dados_filtrados$va_servicos,
  dados_filtrados$va_adespss
) / 10^6)

summary(t(dados_y))

# =============================================================
# Matriz de covariáveis X
# Linhas: intercepto, log(populacao), tx_desemprego, tx_participacao
# Colunas: combinações UF x ano
# =============================================================

X <- rbind(
  rep(1, nrow(pnad)),                  # intercepto
  log(pnad$populacao),                 # log-população (escala similar às respostas)
  pnad$tx_desemprego,                  # taxa de desemprego
  pnad$tx_participacao                 # taxa de participação
)

cat("Dimensões de X com intercepto:", dim(X), "\n")

# =============================================================
# Estimação de B e resíduos
# =============================================================

Beta_new <- dados_y %*% t(X) %*% solve(X %*% t(X))
mu_new   <- Beta_new %*% X
E_new    <- dados_y - mu_new

cat("Norma E_hat sem intercepto:", norm(
  dados_y - (dados_y %*% t(X[2:4,]) %*% solve(X[2:4,] %*% t(X[2:4,]))) %*% X[2:4,], "F"
), "\n")
cat("Norma E_hat com intercepto:", norm(E_new, "F"), "\n")

cat("\nPrimeiras colunas mu com intercepto:\n")
print(round(mu_new[, 1:5], 4))





UF <- c("DF", "GO", "MT", "MS", "AC", "AP", "AM", "PA", "RO", "RR", "TO",
            "AL", "BA", "CE", "MA", "PB", "PE", "PI", "RN", "SE", "ES", "MG",
            "RJ", "SP", "PR", "RS", "SC")


# Define coordinates as a named list
coords1 <- matrix(c(
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



coords <- matrix(c(
  -47.86, -15.83,  # DF
  -49.86, -15.98,  # GO
  -55.42, -12.64,  # MT
  -54.54, -20.51,  # MS
  -70.55,  -8.77,  # AC
  -51.77,   1.41,  # AP
  -65.10,  -3.47,  # AM
  -52.48,  -3.79,  # PA
  -63.34, -10.83,  # RO
  -61.33,   1.99,  # RR
  -48.26,  -9.46,  # TO
  -36.82,  -9.62,  # AL
  -41.71, -13.29,  # BA
  -39.53,  -5.20,  # CE
  -45.44,  -5.42,  # MA
  -36.72,  -7.28,  # PB
  -37.86,  -8.38,  # PE
  -42.28,  -6.60,  # PI
  -36.59,  -5.81,  # RN
  -37.45, -10.57,  # SE
  -40.34, -19.19,  # ES
  -44.38, -18.10,  # MG
  -42.66, -22.25,  # RJ
  -48.79, -22.19,  # SP
  -51.55, -24.89,  # PR
  -53.50, -30.17,  # RS
  -50.95, -27.45   # SC
), ncol = 2, byrow = TRUE)*1000




rownames(coords1) <- UF
colnames(coords1) <- c("x", "y")
dim(coords1)

#as.matrix(dist(coords1))/ max(as.matrix(dist(coords1)), na.rm = TRUE)
source(file.path(base_path, "Functions SpatioFile_b2.R"))

# Calibrar intervalo de phi com base nas distâncias reais
D_vals <- as.matrix(dist(coords1))
D_vals_vec <- D_vals[D_vals > 0]
phi_range_correto <- c(min(D_vals_vec) / 10, max(D_vals_vec) * 2)
phi_init <- median(D_vals_vec)

cat("phi_init:", phi_init, "\n")
cat("phi_range: [", phi_range_correto[1], ",", phi_range_correto[2], "]\n")

# Rodar modelos
ResC_exp      <- ML.MatrixRegreSPatioT(dados_y, X, coords1, precision=1e-9, MaxIter=500, corr="exponential", phi_range=phi_range_correto, phi_init=phi_init)
ResC_spherical <- ML.MatrixRegreSPatioT(dados_y, X, coords1, precision=1e-9, MaxIter=50, corr="spherical",   phi_range=phi_range_correto, phi_init=phi_init)
ResC_gaus     <- ML.MatrixRegreSPatioT(dados_y, X, coords1, precision=1e-9, MaxIter=50, corr="gaussian",    phi_range=phi_range_correto, phi_init = phi_range_correto[1] * 10)
ResC_cubic    <- ML.MatrixRegreSPatioT(dados_y, X, coords1, precision=1e-9, MaxIter=50, corr="cubic",       phi_range=phi_range_correto, phi_init=phi_init)
nu <- 1.5
ResC_matern   <- ML.MatrixRegreSPatioT(dados_y, X, coords1, precision=1e-9, MaxIter=50, corr="matern",      phi_range=phi_range_correto, phi_init=phi_init, nu=nu)



library(xtable)
xtable(format(round(ResC_spherical$Sigma, 3),
              scientific = FALSE
              ))
log(2.277632)
#format(round(ResC_matern$Beta, 7), nsmall = 7, scientific = FALSE)
library(xtable)
xtable(format(round(ResC_spherical$Sigma, 7),
       scientific = FALSE,
       nsmall = 7))

ResC_exp$BIC
ResC_gaus$BIC
ResC_cubic$BIC
ResC_spherical$BIC
ResC_matern$BIC




library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
names(dados_filtrados)
# ── 1. Reshape ─────────────────────────────────────────────────────────────────
pib_long <- dados_filtrados %>%
  pivot_longer(
    cols = c(va_agropecuaria, va_industria, va_servicos, va_adespss),
    names_to = "sector",
    values_to = "gdp"
  ) %>%
  mutate(
    sector = factor(sector,
                    levels = c("va_agropecuaria", "va_industria", "va_servicos", "va_adespss"),
                    labels = c("Agriculture", "Industry", "Services", "Public administration")
    )
  )
# ── 2. Build exactly 27 shapes and linetypes ───────────────────────────────────
all_shapes <- rep(c(16, 17, 15, 18, 1, 2, 0, 5, 8), length.out = 27)

all_linetypes <- rep(c("solid","dashed","dotted","dotdash","longdash","twodash"),
                     length.out = 27)

# ── 3. Plot ────────────────────────────────────────────────────────────────────
p <- ggplot(pib_long,
            aes(x        = ano,
                y        = gdp,
                colour   = sigla_uf,
                group    = sigla_uf,
                shape    = sigla_uf,
                linetype = sigla_uf)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 2) +
  facet_wrap(~ sector, scales = "free_y", ncol = 2) +
  scale_colour_hue() +
  scale_shape_manual(values   = all_shapes) +
  scale_linetype_manual(values = all_linetypes) +
  scale_y_continuous(labels = label_comma()) +
  scale_x_continuous(breaks = seq(2010, 2020, by = 1)) +
  labs(
    x        = "Year",
    y        = "GDP (thousands of BRL)",
    colour   = "State",
    shape    = "State",
    linetype = "State"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 12),
    legend.position  = "right",
    legend.key.width = unit(1.5, "cm"),
    legend.key.height= unit(0.45, "cm"),
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 9, face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    panel.background = element_blank()
  ) +
  guides(
    colour   = guide_legend(ncol = 2, override.aes = list(linewidth = 0.9)),
    shape    = guide_legend(ncol = 2),
    linetype = guide_legend(ncol = 2)
  )

print(p)
ggsave("gdp_by_sector.eps", p, width = 12, height = 7, dpi = 300)
##########################################################
##########################################################
##########################################################

names(dfx_long)
dfx_longer <- dfx_long %>%
  pivot_longer(
    cols = c(IDHM_Educacao, IDHM_Longevidade, IDHM_Renda),
    names_to  = "indicator",
    values_to = "value"
  ) %>%
  mutate(
    indicator = factor(indicator,
                       levels = c("IDHM_Educacao", "IDHM_Longevidade", "IDHM_Renda"),
                       labels = c("HDI Education", "HDI Longevity", "HDI Income")
    ),
    UF = factor(UF, levels = c(
      "DF", "GO", "MT", "MS", "AC", "AP", "AM", "PA", "RO", "RR", "TO",
      "AL", "BA", "CE", "MA", "PB", "PE", "PI", "RN", "SE", "ES", "MG",
      "RJ", "SP", "PR", "RS", "SC"
    ))
  )

# ── 2. Named shapes and linetypes (keyed to your exact order) ─────────────────
state_order <- c(
  "DF", "GO", "MT", "MS", "AC", "AP", "AM", "PA", "RO", "RR", "TO",
  "AL", "BA", "CE", "MA", "PB", "PE", "PI", "RN", "SE", "ES", "MG",
  "RJ", "SP", "PR", "RS", "SC"
)

all_shapes <- setNames(
  rep(c(16, 17, 15, 18, 1, 2, 0, 5, 8), length.out = 27),
  state_order
)

all_linetypes <- setNames(
  rep(c("solid","dashed","dotted","dotdash","longdash","twodash"), length.out = 27),
  state_order
)

# ── 3. Plot ────────────────────────────────────────────────────────────────────
p <- ggplot(dfx_longer,
            aes(x        = year,
                y        = value,
                colour   = UF,
                group    = UF,
                shape    = UF,
                linetype = UF)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 2) +
  facet_wrap(~ indicator, scales = "free_y", ncol = 3) +
  scale_colour_hue() +
  scale_shape_manual(values    = all_shapes) +
  scale_linetype_manual(values = all_linetypes) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  scale_x_continuous(breaks = seq(min(dfx_longer$year),
                                  max(dfx_longer$year), by = 1)) +
  labs(
    x        = "Year",
    y        = "IDHM",
    colour   = "State",
    shape    = "State",
    linetype = "State"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 12),
    legend.position   = "right",
    legend.key.width  = unit(1.5, "cm"),
    legend.key.height = unit(0.45, "cm"),
    legend.text       = element_text(size = 7),
    legend.title      = element_text(size = 9, face = "bold"),
    axis.text.x       = element_text(angle = 45, hjust = 1)
  ) +
  guides(
    colour   = guide_legend(ncol = 2, override.aes = list(linewidth = 0.9)),
    shape    = guide_legend(ncol = 2),
    linetype = guide_legend(ncol = 2)
  )

print(p)
ggsave("idhm_by_indicator.eps", p, width = 14, height = 6, dpi = 600)

##########################################################################
###########################################################################


###########################################
###envelope###############################
beta_hat  <- as.matrix(ResC_exp$Beta)
sigma2_hat <- ResC_exp$sigma2
phi_hat    <- ResC_exp$phi
rho_hat    <- ResC_exp$rho
S_hat      <- as.matrix(ResC_exp$Sigma)

D  <- as.matrix(dist(coords1))
TT <- 10
t  <- 1:TT

# CORRETO: sigma2 apenas em PsiS
Psi_spt  <- sigma2_hat * exp(-D / phi_hat)          # sigma2 * C_espacial
Psi_temp <- rho_hat^abs(outer(t, t, "-"))           # apenas C_temporal, SEM sigma2
Psi_hat  <- kronecker(R_spt, R_temp)

E_hat <- dados_y - beta_hat %*% X

row_wise <- function(E, Sigma_hat, Psi_hat) {
  
  r_i <- numeric(nrow(E))
  
  # Loop over columns j = 1 to p
  for (i in 1:nrow(E)) {
    Ei <- E[i,]   # j-th column of E (n x 1 vector)
    Shat_ii<-Sigma_hat[i,i]
    numerator <- t(Ei) %*% solve(Psi_hat) %*% Ei
    denominator <- Shat_ii
    r_i[i] <- numerator / denominator
  }
  
  return(r_i)
}
ri<-row_wise(E_hat, S_hat,Psi_hat)
ri

p_dim <- nrow(E_hat)
r_dim <- ncol(E_hat)

mahalanobis_diagnostic <- function(E, Sigma, Psi) {
  # Check dimensions
  #n <- nrow(E)
  
  # Invert Sigma_hat (use solve for small matrices; consider chol2inv for stability)
  Sigma_inv <- solve(Sigma)
  
  # Initialize result vector
  d_squared <- numeric(ncol(E))
  
  # Loop over columns j = 1 to p
  for (j in 1:ncol(E)) {
    Ej <- E[, j]   # j-th column of E (n x 1 vector)
    Psijj<-Psi[j,j]
    numerator <- t(Ej) %*% Sigma_inv %*% Ej
    denominator <- Psijj
    d_squared[j] <- numerator / denominator
  }
  
  return(d_squared)
}
#t(E[, 1]) %*% solve(S_hat) %*% E[, 1] /Psi_hat[1,1]


#summary(d_squared)
df_val <- p_dim * r_dim


# Mahalanobis diagnostic
d_squared       <- mahalanobis_diagnostic(E_hat, S_hat, Psi_hat)
threshold_value <- qchisq(0.95, df = p_dim)
pval            <- 1 - pchisq(d_squared, df = p_dim)

cat("=== Matrix Normality Test ===\n")
cat(sprintf("Dimensions: p = %d, r = %d, df = %d\n", p_dim, r_dim, p_dim * r_dim))
cat(sprintf("Threshold chi2(0.95, %d) = %.3f\n", p_dim, threshold_value))
cat(sprintf("Points above threshold = %d\n", sum(d_squared > threshold_value)))
cat(sprintf("Decision (alpha=0.05): %s H0\n",
            ifelse(mean(pval) > 0.05, "Fail to reject", "Reject")))

# Data frame para o plot
df_mahal <- data.frame(
  Index           = 1:length(d_squared),
  d_squared       = d_squared,
  above_threshold = d_squared > threshold_value
)


y_max <- threshold_value
# Create data frame
df <- data.frame(
  Index = 1:n,
  d_squared = d_squared,
  above_threshold = d_squared > threshold_value
)

# Plot
library(ggplot2)

# Plot Mahalanobis
mahalanobis_pib <- ggplot(df_mahal, aes(x = Index, y = d_squared)) +
  geom_segment(aes(xend = Index, yend = 0, color = d_squared),
               size = 0.3, alpha = 1) +
  geom_point(aes(color = d_squared), size = 2) +
  geom_vline(xintercept = 10, linetype = "dotted", color = "dodgerblue", size = 0.8) +
  scale_color_gradient(low = "steelblue", high = "red", name = expression(d[j]^2)) +
  labs(
    title    = "Mahalanobis Distance by Column",
    subtitle = paste0("Threshold Chi-square(0.95, ", p_dim, ") = ",
                      round(threshold_value, 3),
                      " - Points above threshold = ",
                      sum(df_mahal$above_threshold)),
    x = "Index j (1-270)",
    y = expression(d[j]^2)
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "gray85", size = 0.2),
    panel.grid.minor = element_line(color = "gray95", size = 0.1),
    axis.line  = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3)
  ) +
  coord_cartesian(ylim = c(0, 50), xlim = c(1, 270)) +
  scale_x_continuous(breaks = seq(0, 270, by = 10), limits = c(1, 270))

ggsave("mahalanobis_pib.eps", plot = mahalanobis_pib, width = 10, height = 6, device = "eps")


### Global Statistic ###
E_hat    <- dados_y - beta_hat %*% X
Sigma_E  <- 1/r_dim * (E_hat %*% solve(Psi_hat) %*% t(E_hat))
Psi_E    <- 1/p_dim * (t(E_hat) %*% solve(S_hat) %*% E_hat)
T_global <- 0.5 * (1/p_dim * sum(diag(solve(S_hat) %*% Sigma_E)) +
                   1/r_dim * sum(diag(solve(Psi_hat) %*% Psi_E)))
T_global

###envelope####

Sigma_inv_sqrt <- function(Sigma, method = c("eigen", "chol"), eps = 1e-8) {
  method <- match.arg(method)
  Sigma  <- (Sigma + t(Sigma))/2  # garantir simetria
  if (method == "eigen") {
    ev  <- eigen(Sigma, symmetric = TRUE)
    lam <- pmax(ev$values, eps)   # evitar divisC#o por zero
    return(ev$vectors %*% diag(1/sqrt(lam)) %*% t(ev$vectors))  # raiz-inversa simC)trica
  } else {
    L <- chol(Sigma)                                   # Sigma = L %*% t(L)
    return(backsolve(L, diag(nrow(L))))                # whitener (nC#o simC)trico)
  }
}



S_inv_sqrt <- Sigma_inv_sqrt(S_hat, method = "eigen")
E_hat<-dados_y -  beta_hat %*% X
E_star <- S_inv_sqrt %*% E_hat %*%  Sigma_inv_sqrt(Psi_hat, method = "eigen")

p_dim <- nrow(E_star)
r_dim <- ncol(E_star)
cat("p_dim:", p_dim, "r_dim:", r_dim, "\n")

test_chisq_pr <- function(E_hat_star, alpha = 0.05) {
  df_val <- length(E_hat_star)   # p_dim * r_dim
  Q      <- sum(E_hat_star^2)
  p_val  <- 1 - pchisq(Q, df_val)
  cat("Degrees of Freedom:", df_val, "\n")
  cat("Observed Q:", Q, "\n")
  cat("P-value:", p_val, "\n")
  cat("Decision:", ifelse(p_val > alpha, "Fail to reject H0", "Reject H0"), "\n")
  return(list(df = df_val, Q = Q, p_value = p_val))
}

z_vec <- as.vector(E_star)
test_chisq_pr(z_vec)


str(E_star)
class(E_star)
dim(E_star)
is.matrix(E_star)
str(E_star)



cell_diagnostic <- function(E) {
  z_ij <- matrix(NA, nrow = nrow(E), ncol = ncol(E))
  for (i in 1:nrow(E))
    for (j in 1:ncol(E))
      z_ij[i, j] <- E[i, j]
  return(z_ij)
}  

zij   <- cell_diagnostic(E_star)   # E_star ainda é matriz
z_vec <- as.vector(zij)     

library(PerformanceAnalytics)
postscript("envelope_pib_2010_2020.eps", 
           width = 8, 
           height = 8, 
           horizontal = FALSE, 
           paper = "special")
chart.QQPlot(z_vec , 
             distribution = "norm", 
             envelope = 0.95,          # 95% simulation envelope
             main = "Envelope Diagnostic Plot for Residuals of the Exponential Model")
dev.off()
#####################################################################
#####################################################################

library(ggplot2)
library(dplyr)
library(tidyr)

# Set row names for the 4 sectors
rownames(E_star) <- c("Agriculture", "Industry", "Services", "PublicAdmin")

# Convert to long format
df <- as.data.frame(t(E_star))
colnames(df) <- c("Agriculture", "Industry", "Services", "PublicAdmin")
df$Column <- 1:ncol(E_star)

df_long <- df %>%
  pivot_longer(
    cols      = -Column,
    names_to  = "Sector",
    values_to = "zii"
  ) %>%
  mutate(
    Sector = factor(Sector, levels = c("Agriculture", "Industry", "Services", "PublicAdmin")),
    Significance = ifelse(abs(zii) > 2, "Significant", "Normal")
  )

# Plot

postscript("cellwise_pib.eps", 
           width = 12, 
           height = 12, 
           horizontal = FALSE, 
           paper = "special")

cellwise_pib<-ggplot(df_long, aes(x = Column, y = zii, color = Significance)) +
  geom_point(size = 0.7, alpha = 1) +
  
  # Reference lines
  geom_hline(yintercept =  0, color = "black",  linewidth = 0.3) +
  geom_hline(yintercept =  2, color = "red", linetype = "dotted", linewidth = 0.5) +
  geom_hline(yintercept = -2, color = "red", linetype = "dotted", linewidth = 0.5) +
  
  facet_wrap(~ Sector, ncol = 1, strip.position = "right", scales = "free_y") +
  
  scale_color_manual(
    values = c("Normal" = "#3355cc", "Significant" = "#cc2200"),
    name   = "Significance Level"
  ) +
  scale_x_continuous(
    breaks = seq(1, 270, by = 30),
    expand = c(0.01, 0)
  )+
  scale_y_continuous(
    breaks = c(-2, 0, 2)
  ) +
  labs(
    x = "Column (Spatio\u2013temporal Condition)",
    y = "Standardized Residual (zii)"
  ) +
  
  theme(
    # Panel appearance — no border, no background
    panel.background   = element_blank(),
    panel.border       = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey75", linetype = "dotted", linewidth = 0.4),
    
    # Axes
    axis.line.x        = element_line(color = "black", linewidth = 0.5),
    axis.line.y        = element_line(color = "black", linewidth = 0.5),
    axis.ticks         = element_line(linewidth = 0.3),
    axis.text          = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    
    # Strip labels: vertical text, no box
    strip.background   = element_blank(),
    strip.text.y.right = element_text(angle = 90, size = 12, vjust = 0.5, margin = margin(l = 4)),
    
    # Spacing between panels
    panel.spacing      = unit(0.15, "lines"),
    
    # Legend
    legend.position    = "right",
    legend.title       = element_text(size = 8),
    legend.text        = element_text(size = 7),
    legend.key         = element_blank(),
    
    plot.margin        = margin(3, 3, 3, 3)
  )
print(cellwise_pib)
dev.off()
ggsave("cellwise_pib.eps", plot = cellwise_pib,
       width = 12, height = 10, device = "eps")
