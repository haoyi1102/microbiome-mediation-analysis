library(ggplot2)

# 创建数据框
time_data <- data.frame(
  Method = c("MODIMA", "MEDTEST", "PERMANOVA-MED", "NPEM", "CCMM", "HIMA", "LDM-MED", "Medzim", "Microbvs", "SparseMCMM"),
  Time = c(0.373, 0.481, 0.891, 1.445, 7.002, 13.349, 17.386, 29.162, 182.077, 500)  # SparseMCMM 时间设为最差
)

# 计算性能 (效率的倒数：1 / Time)
time_data$Performance <- 1 / time_data$Time

# 绘制散点图
ggplot(time_data, aes(x = reorder(Method, Performance), y = Performance)) +
  geom_point(size = 3, color = "blue") +
  geom_hline(yintercept = 1 / c(1, 10, 100), linetype = "dashed", color = "red") +  # 添加参考线
  annotate("text", x = 2, y = 1 / 1, label = "1 second", vjust = -0.5, color = "red") +
  annotate("text", x = 2, y = 1 / 10, label = "10 seconds", vjust = -0.5, color = "red") +
  annotate("text", x = 2, y = 1 / 100, label = "100 seconds", vjust = -0.5, color = "red") +
  theme_minimal() +
  labs(
    title = "Time Efficiency Comparison Across Methods",
    x = "Method",
    y = "Time Efficiency"  # 修改纵坐标标签为 Performance
  ) +
  scale_y_log10(
    breaks = c(1 / 100, 1 / 10, 1 / 1, 1 / 0.1)  # 对数比例
  ) +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),  # 居中标题
    axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转X轴标签
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12)
  )


