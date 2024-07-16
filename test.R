rm(list = ls())

library(MASS)
library(stats)
library(graphics)
library(energy)  
library(matrixStats)
library(devtools)

##
source("modima.R")
source("MedOmniTest.R")
# install_github("chanw0/SparseMCMM")
library(SparseMCMM)
#install.packages("ccmm")
library(ccmm)
# install_github("yijuanhu/LDM")
library(LDM)
#install_github( "mkoslovsky/MicroBVS")
library(MicroBVS)
#devtools::install_github("quranwu/MedZIM")
library(MedZIM)
#install.packages("microHIMA_1.0.tar.gz", repos = NULL, type = "source")
library(microHIMA) #install.packages("ncvreg") # install.packages("hommel")



set.seed(1234)
# 读取数据
Y <- read.csv("outcome.csv", header = TRUE)
treatment <- read.csv("treatment.csv", header = TRUE)
T <- treatment$treatment_numeric
M_absolute <- read.csv("mediator(absolute).csv", header = TRUE)
feature <- read.csv("feature_metadata.csv", header = TRUE)

# 准备数据
Y_vector <- as.numeric(Y$x)  # 确保是数值型向量
T_vector <- as.numeric(T)  # 转置并转换为向量
M_a_matrix <- as.matrix(t(M_absolute))  # 确保是数值型矩阵 # ab
colnames(M_a_matrix) <- paste0("taxon_", 1:ncol(M_a_matrix))
M_matrix <- M_a_matrix / rowSums(M_a_matrix) # comp
# pseudo_count <- 1e-5
# temp_m <- M_a_matrix
# temp_m[temp_m == 0] <- pseudo_count
# M_nz_matrix <- temp_m / rowSums(temp_m) #non-zero comp

M_nz_matrix <- apply(M_matrix, 1, function(row) {
  min_value <- min(row[row > 0])
  pseudo_count <- min_value / 100
  row[row == 0] <- pseudo_count
  return(row)
})
M_nz_matrix <- t(M_nz_matrix) #non-zero comp
class(M_nz_matrix)
##################33microHIMA
mhima_fit <- mhima(exposure = T_vector, covariates = NULL, otu.com = M_nz_matrix, outcome = Y_vector)
mhima_fit$ID
print(mhima_fit)

################ MODIMA
T_dist <- dist(T_vector)
Y_dist <- dist(Y_vector)
M_dist <- dist(M_a_matrix)

modima_result <- modima(T_dist, M_dist, Y_dist, nrep=999)
modima_result$p.value
print(modima_result)
## 数据是否compositional有显著影响 #数据是否0值影响不大

################ MEDTEST

m_list <- list(euc = dist(M_a_matrix))

# 运行 MedOmniTest
result <- MedOmniTest(x = T_vector, y = Y_vector, m.list = m_list, z = NULL, nperm = 999)

# 提取 p 值
chen_p <- result$permP

# 打印结果
print(result)
print(paste("MedTest p-value:", chen_p))
## 数据是否compositional有显著影响 #数据是否0值影响不大


############### SparseMCMM
if (any(!is.finite(M_nz_matrix))) {
  stop("M_matrix contains non-finite values!")
}
res = SparseMCMM(T_vector, M_nz_matrix, Y_vector, n.split=1, num.per=200)
ress = res$Test
p_value = ress["OME"]
# res = SparseMCMM(T_vector, M_nz_matrix, Y_vector, n.split=1, num.per=200)
# > print(res)
# $`Esitmated Causal Effects`
# DE         ME         TE 
# -0.5924898 -4.3311540 -4.9236438 
# 
# $`Compontent-wise ME`
# taxon_1     taxon_2     taxon_3     taxon_4     taxon_5 
# 0.05934267  1.24261963  0.01766287  0.29616179 -0.23424203 
# taxon_6     taxon_7     taxon_8     taxon_9    taxon_10 
# 0.34983561 -5.83868936 -0.02322553  0.43861452 -0.63923414 
# 
# $Test
# OME       CME 
# 0.1194030 0.1492537 
# 
# > 
############ CCMM
#help(ccmm)
result <- ccmm(Y_vector, M_nz_matrix, T_vector, x = NULL, w = NULL, method.est.cov = "bootstrap", n.boot = 2000,
     sig.level = 0.05, tol = 1e-06, max.iter = 5000)
CI<-result$TIDE.CI
print(result)

### ldm-med
data <- data.frame(Y = Y_vector, T = T_vector)
result <- ldm(
  formula = M_nz_matrix ~ T + Y,
  data = data,
  seed = 1234,
  test.mediation = TRUE
)
?ldm
print(result$med.detected.otu.omni)

##### microbvs
model_real <- MCMC_Med(trt = T_vector, Y = Y_vector, Z = M_matrix, taxa = 2)
result_real <- Selection_Med1(model = model_real)

result_global <- Selection_Med2(model = model_real)

model_global <- MCMC_Med(trt = T_vector, Y = Y_vector, Z = M_matrix, seed = 1234)

##### MedZim

#?MedZIM_func
#?data_ZIM
taxon_name <- "taxon_"


libsize <- colSums(M_absolute)
dat <- data.frame(Y_vector, T_vector, libsize)
dat <- cbind(dat, M_matrix)
results <- MedZIM_func(
  dat = dat,
  xVar = "T_vector",
  yVar = "Y_vector",
  taxon_name = taxon_name,
  libSize_name = "libsize",
  obs_gt_0 = 2,
  obs_eq_0 = 2,
  inter_x_mg0 = TRUE,
  inter_x_m = FALSE,
  eval.max = 200,
  iter.max = 200,
  x_from = 0,
  x_to = 1,
  type1error = 0.05,
  paraJobs = 2
)
results$fullList$taxon_1
