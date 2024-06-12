rm(list = ls())

library(SparseDOSSA2)
library(SparseDOSSA2)
library(magrittr)
library(dplyr)
library(ggplot2)

template <- "Stool"  # 使用预训练的Stool模板
n_sample <- 100  # 生成100个样本
n_feature <- 100  # 每个样本100个特征
metadata_effect_size <- 1  # 设置metadata影响大小
perc_feature_spiked_metadata <- 0.1  # 10%的特征与metadata相关联

# 生成处理组metadata矩阵
set.seed(123)
metadata <- data.frame(treatment = rep(c("Control", "Treatment"), each = 50))

# 将metadata转换为模型矩阵，进行哑编码
metadata_matrix <- model.matrix(~ treatment - 1, data = metadata)

# 生成合成数据
sim_data <- SparseDOSSA2(template = template, 
                         n_sample = n_sample, 
                         n_feature = n_feature,
                         spike_metadata = "abundance",  # 关联metadata和丰度
                         metadata_effect_size = metadata_effect_size,
                         perc_feature_spiked_metadata = perc_feature_spiked_metadata,
                         metadata_matrix = metadata_matrix,  # 使用转换后的metadata矩阵
                         verbose = TRUE)

# 提取生成的数据
synthetic_data <- sim_data$simulated_data
generated_metadata <- sim_data$spike_metadata$metadata_matrix

# 查看生成的metadata
head(generated_metadata)




