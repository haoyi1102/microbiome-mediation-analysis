rm(list = ls())
library(SparseDOSSA2)
library(caret) 
library(Matrix)
library(matrixcalc)
# 读取生成的50x50数据集
data <- read.csv("generated_50x50_dataset.csv", row.names = 1)
# 标准化数据


set.seed(123)

# 生成50x50的正态分布数据
data <- matrix(rnorm(2500, mean = 100, sd = 100), nrow = 50, ncol = 50)

# 将数据转换为非负的计数数据
data <- round(abs(data))

# 设置列名
colnames(data) <- paste0("X7000000", 0:49)

# 设置行名
rownames(data) <- paste0("k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Taxon_", 0:49)

data <- data[1:50, 1:50]


# 拟合SparseDOSSA2模型
control_params <- control_fit(
  maxit = 200,          # 增加EM算法的最大迭代次数
  rel_tol = 0.01,      # 相对变化收敛阈值
  abs_tol = 0.01,      # 绝对变化收敛阈值
  verbose = TRUE        # 输出详细运行信息
)
fitted <- fit_SparseDOSSA2(data = data, control = control_params)
?fit_SparseDOSSA2
?control_fit
# 显示前几行数据
head(data)
data <- data[1:10, 1:10]
# data("Stool_subset", package = "SparseDOSSA2")
fitted <- ?fit_SparseDOSSA2(data = data,
                           control = list(verbose = TRUE))

# fitted absence probabilities for first two features in zero-inflated model 
fitted$EM_fit$fit$pi0
fitted$EM_fit$fit$pi0 <- rep(0, length(fitted$EM_fit$fit$pi0))
temp = fitted

save(fitted, file = "fitted_model_1.RData")
save(fitted, file = "fitted_model_2.RData")
load("fitted_model_1.RData")
Stool_subset_simulation <- SparseDOSSA2(template = fitted, 
                                        n_sample = 10, 
                                        n_feature = 10,
                                        new_features = TRUE,perc_feature_spiked_metadata = 0.1,
                                        verbose = TRUE)
a = Stool_subset_simulation$simulated_data
sum(a==0)
Stool_subset_simulation$spike_metadata

Stool_subset_simulation <- SparseDOSSA2(
  template = fitted, 
  n_sample = 50, 
  n_feature = 20,
  spike_metadata = "abundance", 
  metadata_effect_size = 8,
  perc_feature_spiked_metadata = 0.1,
  metadata_matrix = metadata_matrix,
  median_read_depth = 100000,  
  verbose = TRUE
)

n_sample = 100, template = "Stool", n_feature = 20, 
metadata_effect_size = 1, perc_feature_spiked_metadata = 0.1, 
median_read_depth = 10000, alpha_0 = 1, alpha_T = 0.5, 
alpha_M_value = 0.2, noise_sd = 0.1

metadata <- data.frame(treatment = rep(c("Control", "Treatment"), each = 25))
metadata$treatment_binary <- ifelse(metadata$treatment == "Control", 0, 1)
metadata_matrix <- as.matrix(metadata$treatment_binary)




# 加载SparseDOSSA2包
library(SparseDOSSA2)

# 创建样本数量
n_samples <- 20
n_features <- 10

# 创建元数据矩阵
set.seed(123)
metadata_matrix <- data.frame(
  Treatment = sample(c(0, 1), n_samples, replace = TRUE)  # 二元治疗变量
)

# 设置spike_metadata数据框，定义元数据与微生物特征的关联
spike_metadata <- data.frame(
  metadata_datum = 1,                # 元数据矩阵中的第1列
  feature_spiked = paste0("Feature", 1:10),  # 关联的微生物特征
  associated_property = "abundance", # 修改特征的丰度
  effect_size = 0.5                  # 关联的效果大小（log fold change）
)

# 使用元数据和配置生成模拟数据
simulated_data <- SparseDOSSA2(
  template = "Stool", 
  n_sample = n_samples, 
  n_feature = n_features,
  new_features = TRUE, 
  spike_metadata = spike_metadata,  # 提供详细的spike-in配置
  metadata_matrix = metadata_matrix, # 提供元数据矩阵
  verbose = TRUE
)

# 查看生成的数据
head(simulated_data$simulated_data)
head(simulated_data$metadata_matrix)









