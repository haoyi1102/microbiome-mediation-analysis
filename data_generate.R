

rm(list = ls())

library(SparseDOSSA2)
library(magrittr)
library(dplyr)
library(ggplot2)

# 设置模拟参数
template <- "Stool"  # 使用预训练的Stool模板
n_sample <- 100  # 生成100个样本
n_feature <- 10  # 每个样本100个特征
metadata_effect_size <- 1  # 设置metadata影响大小
perc_feature_spiked_metadata <- 0.1  # 10%的特征与metadata相关联

# 生成处理组metadata矩阵
set.seed(123)
metadata <- data.frame(treatment = rep(c("Control", "Treatment"), each = 50))

# 将metadata转换为模型矩阵
metadata_matrix <- model.matrix(~ treatment - 1, data = metadata)

# 生成合成数据
sim_data <- SparseDOSSA2(template = template, 
                         n_sample = n_sample, 
                         n_feature = n_feature,
                         spike_metadata = "abundance", 
                         metadata_effect_size = metadata_effect_size,
                         perc_feature_spiked_metadata = perc_feature_spiked_metadata,
                         metadata_matrix = metadata_matrix,  
                         verbose = TRUE)

# 提取生成的计数数据
synthetic_data <- sim_data$simulated_data
generated_metadata <- sim_data$spike_metadata$metadata_matrix
feature_metadata_spike_df <- sim_data$spike_metadata$feature_metadata_spike_df

# 将处理组变量转换为数值型：Control -> 0, Treatment -> 1
metadata$treatment_numeric <- as.numeric(metadata$treatment == "Treatment")

# 假设预定义的回归系数
alpha_0 <- 1
alpha_T <- 0.5  # 处理变量的系数
alpha_M <- rep(0.2, n_feature - 1)  # 微生物类群的系数

# 初始化 outcome Y
outcome_Y <- rep(alpha_0, n_sample)

# 添加处理变量项
outcome_Y <- outcome_Y + metadata$treatment_numeric * alpha_T

# 计算对数对比项并添加到 outcome Y
for (i in 1:n_sample) {
  log_contrasts <- numeric(n_feature - 1)
  for (j in 1:(n_feature - 1)) {
    ratio <- (synthetic_data[j, i] + 1) / (synthetic_data[n_feature, i] + 1)
    if (ratio == 0) {
      log_contrasts[j] <- 0
    } else {
      log_contrasts[j] <- log(ratio)
    }
  }
  outcome_Y[i] <- outcome_Y[i] + sum(alpha_M * log_contrasts)
}

# 查看生成的 outcome Y 值
head(outcome_Y)

write.csv(outcome_Y, file = "outcome.csv", row.names = FALSE)
write.csv(synthetic_data, file = "mediator(absolute).csv", row.names = FALSE)
write.csv(metadata, file = "treatment.csv", row.names = FALSE)
write.csv(feature_metadata_spike_df, file = "feature_metadata.csv", row.names = FALSE)
