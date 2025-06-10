original_simulation <- SparseDOSSA2(
  template = "Stool",  # 使用预训练的 Stool 模板
  new_features = TRUE,  # 生成新的特征
  n_sample = 100,  # 模拟100个样本
  n_feature = 1000,  # 模拟1000个特征
  spike_metadata = "none",  # 没有元数据关联
  verbose = F  # 不打印详细信息
)


original_template <- original_simulation$template
original_Ffit <- original_template$F_fit
original_depth <- original_template$depth_fit
original_EMfit <- original_template$EM_fit
original_pi0 <- original_EMfit$fit$pi0 |> unname()  # 提取生成零值的概率
original_mu <- original_EMfit$fit$mu |> unname()  # 提取丰度的均值
original_sigma <- original_EMfit$fit$sigma |> unname()  # 提取丰度的标准差


original_median_abundance <- rep(0, length(original_pi0))
nonzero_mask <- original_pi0 < 0.3  # 只考虑流行率大于50%的特征
truncated_quantile <- qnorm(0.5 - original_pi0[nonzero_mask])
original_median_abundance[nonzero_mask] <- exp(original_mu[nonzero_mask] +
                                                 truncated_quantile * original_sigma[nonzero_mask])
original_median_relabd <- original_median_abundance / sum(original_median_abundance)
subset_taxa <- which(original_median_relabd > 2e-4)  # 过滤掉中位相对丰度小于2e-4的特征



subset_pi0 <- original_pi0[subset_taxa]
subset_mu <- original_mu[subset_taxa]
subset_sigma <- original_sigma[subset_taxa]

subset_EMfit <- list(fit = list(pi0 = subset_pi0, 
                                mu = subset_mu, 
                                sigma = subset_sigma))

subset_template <- list(EM_fit = subset_EMfit,
                        F_fit = original_Ffit,
                        depth_fit = original_depth)


Stool_simulation <- SparseDOSSA2(
  template = subset_template,  # 使用新的模板
  new_features = TRUE,  # 生成新的特征
  n_sample = 100,  # 模拟50个样本
  n_feature = 20,  # 生成1000个特征
  spike_metadata = "none",
  median_read_depth = 5e4,
  verbose = F
)
# 假设已经生成了 Stool_simulation
# 提取模拟的数据矩阵
simulated_data <- Stool_simulation$simulated_data

# 计算总的元素数量
total_elements <- length(simulated_data)

# 计算0值的数量
zero_count <- sum(simulated_data == 0)

# 计算0值的比例
zero_proportion <- zero_count / total_elements

# 输出0值比例
print(zero_proportion)
# 确保你在当前目录中保存
#saveRDS(subset_template, file = "sparsedossa_baseline.rds")





########

rm(list= ls())
original_simulation <- SparseDOSSA2(
  template = "Stool",  # 使用预训练的 Stool 模板
  new_features = TRUE,  # 生成新的特征
  n_sample = 100,  # 模拟100个样本
  n_feature = 1000,  # 模拟1000个特征
  spike_metadata = "none",  # 没有元数据关联
  verbose = F  # 不打印详细信息
)


original_template <- original_simulation$template
original_Ffit <- original_template$F_fit
original_depth <- original_template$depth_fit
original_EMfit <- original_template$EM_fit
original_pi0 <- original_EMfit$fit$pi0 |> unname()  # 提取生成零值的概率
original_mu <- original_EMfit$fit$mu |> unname()  # 提取丰度的均值
original_sigma <- original_EMfit$fit$sigma |> unname()  # 提取丰度的标准差


original_median_abundance <- rep(0, length(original_pi0))
nonzero_mask <- original_pi0 < 0.8  # 只考虑流行率大于50%的特征
truncated_quantile <- qnorm(0.5 - original_pi0[nonzero_mask])
original_median_abundance[nonzero_mask] <- exp(original_mu[nonzero_mask] +
                                                 truncated_quantile * original_sigma[nonzero_mask])
original_median_relabd <- original_median_abundance / sum(original_median_abundance, na.rm = TRUE)
subset_taxa <- which(original_median_relabd > 0)  # 过滤掉中位相对丰度小于2e-4的特征



subset_pi0 <- original_pi0[subset_taxa]*2
subset_mu <- original_mu[subset_taxa]
subset_sigma <- original_sigma[subset_taxa]

subset_EMfit <- list(fit = list(pi0 = subset_pi0, 
                                mu = subset_mu, 
                                sigma = subset_sigma))

subset_template <- list(EM_fit = subset_EMfit,
                        F_fit = original_Ffit,
                        depth_fit = original_depth)


Stool_simulation <- SparseDOSSA2(
  template = subset_template,  # 使用新的模板
  new_features = TRUE,  # 生成新的特征
  n_sample = 100,  # 模拟50个样本
  n_feature = 20,  # 生成1000个特征
  spike_metadata = "none",
  median_read_depth = 5e4,
  verbose = F
)
# 假设已经生成了 Stool_simulation
# 提取模拟的数据矩阵
simulated_data <- Stool_simulation$simulated_data

# 计算总的元素数量
total_elements <- length(simulated_data)

# 计算0值的数量
zero_count <- sum(simulated_data == 0)

# 计算0值的比例
zero_proportion <- zero_count / total_elements

# 输出0值比例
print(zero_proportion)
# 确保你在当前目录中保存
saveRDS(subset_template, file = "sparsedossa_baseline_high.rds")




#######################33


# 第一步：生成新特征数据，不进行元数据关联
sim_initial <- SparseDOSSA2(
  template = "Stool", 
  n_sample = 10, 
  new_features = TRUE, 
  n_feature = 1000
)

# 获取生成的特征名称
generated_features <- rownames(sim_initial$simulated_data)

# 第二步：使用这些生成的特征来设置 spike_metadata
spike_metadata <- data.frame(
  metadata_datum = 1,  # 关联第一列（treatment）
  feature_spiked = generated_features[1:2],  # 关联前两个生成的特征
  associated_property = c("abundance", "abundance"),  # 选择丰度关联
  effect_size = c(2, 1.5)  # 设置效应大小
)

# 重新运行 SparseDOSSA2，加入元数据关联
metadata_matrix_matrix <- model.matrix(~ treatment, data = data.frame(treatment = sample(c(0, 1), 100, replace = TRUE)))

sim <- SparseDOSSA2(
  template = "Stool", 
  n_sample = 100, 
  new_features = TRUE,  # 仍然生成新特征
  spike_metadata = spike_metadata, 
  metadata_matrix = metadata_matrix_matrix  # 使用元数据矩阵
)

# 提取模拟的微生物数据
simulated_data <- sim$simulated_data

