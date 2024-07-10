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

################ MODIMA
T_dist <- dist(T_vector)
Y_dist <- dist(Y_vector)
M_dist <- dist(M_matrix)

modima_result <- modima(T_dist, M_dist, Y_dist, nrep=999)
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
print(res)

############ CCMM
#help(ccmm)
result <- ccmm(Y_vector, M_nz_matrix, T_vector, x = NULL, w = NULL, method.est.cov = "bootstrap", n.boot = 2000,
     sig.level = 0.05, tol = 1e-06, max.iter = 5000)

print(result)

### ldm-med
data <- data.frame(Y = Y_vector, T = T_vector)
result <- ldm(
  formula = M_nz_matrix ~ T + Y,
  data = data,
  seed = 1234,
  test.mediation = TRUE
)
# 全局中介效应的p值
print(result$med.p.global.omni)
help(ldm)
# 检测到的OTU的中介效应
print(result$med.detected.otu.omni)

##### microbvs
model_real <- MCMC_Med(trt = T_vector, Y = Y_vector, Z = M_matrix, taxa = 1)
result_real <- Selection_Med1(model = model_real)

result_global <- Selection_Med2(model = model_real)

model_global <- MCMC_Med(trt = T_vector, Y = Y_vector, Z = M_matrix, seed = 1234)
