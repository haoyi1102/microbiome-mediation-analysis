# 安装并加载必要的包
if(!requireNamespace("reactable", quietly = TRUE)) install.packages("reactable")
if(!requireNamespace("webshot", quietly = TRUE)) install.packages("webshot")
if(!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if(!requireNamespace("htmlwidgets", quietly = TRUE)) install.packages("htmlwidgets")
if(!requireNamespace("htmltools", quietly = TRUE)) install.packages("htmltools")
if(!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")

library(reactable)
library(webshot)
library(tidyverse)
library(htmlwidgets)
library(htmltools)
library(RColorBrewer)

# 直接创建一个完整的矩阵形式的performance_data数据框

# sample 50 m1 
performance_data <- data.frame(
  Method = c("MedTest", "MODIMA", "CCMM", "LDM", "PERMANOVA"),
  Power_Low = c(0.14, 0.14, 0.96, 0.42, 0.57),
  Power_Medium = c(0.14, 0.14, 0.89, 0.42, 0.76),
  Power_High = c(0.12, 0.16, 0.40, 0.20, 0.70),
  Type1Error_Low = c(0.0010, 0.0010, 0.1000, 0.0010, 0.0125),
  Type1Error_Medium = c(0.0010, 0.0010, 0.1125, 0.0125, 0.0010),
  Type1Error_High = c(0.0010, 0.0010, 0.0250, 0.0010, 0.0010)
)


# sample 100 m1 
performance_data <- data.frame(
  Method = c("MedTest", "MODIMA", "CCMM", "LDM", "PERMANOVA"),
  Power_Low = c(0.18, 0.24, 0.98, 0.77, 0.88),
  Power_Medium = c(0.18, 0.24, 0.85, 0.62, 0.93),
  Power_High = c(0.17, 0.22, 0.26, 0.38, 0.93),
  Type1Error_Low = c(0.0010, 0.0125, 0.0500, 0.0010, 0.0010),
  Type1Error_Medium = c(0.0010, 0.0125, 0.0625, 0.0010, 0.0010),
  Type1Error_High = c(0.0010, 0.0125, 0.0010, 0.0010, 0.0010)
)


# sample 300 m1 
performance_data <- data.frame(
  Method = c("MedTest", "MODIMA", "CCMM", "LDM", "PERMANOVA"),
  Power_Low = c(0.18, 0.24, 0.98, 0.77, 0.88),
  Power_Medium = c(0.18, 0.24, 0.85, 0.62, 0.93),
  Power_High = c(0.17, 0.22, 0.26, 0.38, 0.93),
  Type1Error_Low = c(0.0010, 0.0125, 0.0500, 0.0010, 0.0010),
  Type1Error_Medium = c(0.0010, 0.0125, 0.0625, 0.0010, 0.0010),
  Type1Error_High = c(0.0010, 0.0125, 0.0010, 0.0010, 0.0010)
)

# sample 50 m2 
performance_data <- data.frame(
  Method = c("MedTest", "MODIMA", "CCMM", "LDM", "PERMANOVA"),
  Power_Low = c(0.14, 0.16, 0.71, 0.96, 0.65),
  Power_Medium = c(0.12, 0.16, 0.65, 0.94, 0.64),
  Power_High = c(0.18, 0.12, 0.73, 0.80, 0.50),
  Type1Error_Low = c(0.0010, 0.0010, 0.0125, 0.0125, 0.0010),
  Type1Error_Medium = c(0.0010, 0.0010, 0.0010, 0.0010, 0.0010),
  Type1Error_High = c(0.0010, 0.0010, 0.0125, 0.0010, 0.0010)
)
# sample 100 m2 
performance_data <- data.frame(
  Method = c("MedTest", "MODIMA", "CCMM", "LDM", "PERMANOVA"),
  Power_Low = c(0.18, 0.24, 0.69, 0.97, 0.86),
  Power_Medium = c(0.17, 0.24, 0.70, 0.96, 0.78),
  Power_High = c(0.18, 0.23, 0.65, 0.91, 0.63),
  Type1Error_Low = c(0.0010, 0.0125, 0.0010, 0.0010, 0.0010),
  Type1Error_Medium = c(0.0010, 0.0125, 0.0010, 0.0010, 0.0010),
  Type1Error_High = c(0.0010, 0.0125, 0.0010, 0.0010, 0.0010)
)
# sample 300 m2 
performance_data <- data.frame(
  Method = c("MedTest", "MODIMA", "CCMM", "LDM", "PERMANOVA"),
  Power_Low = c(0.24, 0.64, 0.85, 0.97, 1.00),
  Power_Medium = c(0.24, 0.56, 0.89, 0.98, 1.00),
  Power_High = c(0.23, 0.55, 0.75, 0.98, 1.00),
  Type1Error_Low = c(0.0010, 0.0010, 0.0010, 0.0010, 0.0010),
  Type1Error_Medium = c(0.0010, 0.0010, 0.0010, 0.0010, 0.0010),
  Type1Error_High = c(0.0010, 0.0010, 0.0010, 0.0010, 0.0010)
)

# sample 50 m3 
performance_data_50 <- data.frame(
  Method = c("MedTest", "MODIMA", "CCMM", "LDM", "PERMANOVA"),
  Power_Low = c(0.140, 0.120, 0.580, 0.150, 0.330),
  Power_Medium = c(0.160, 0.140, 0.020, 0.080, 0.120),
  Power_High = c(0.260, 0.180, 0.001, 0.070, 0.100),
  Type1Error_Low = c(0.0010, 0.0010, 0.0125, 0.0125, 0.0010),
  Type1Error_Medium = c(0.0010, 0.0010, 0.0010, 0.0010, 0.0010),
  Type1Error_High = c(0.0010, 0.0010, 0.0125, 0.0010, 0.0010)
)

# sample 100 m3 
performance_data_100 <- data.frame(
  Method = c("MedTest", "MODIMA", "CCMM", "LDM", "PERMANOVA"),
  Power_Low = c(0.180, 0.230, 0.370, 0.230, 0.480),
  Power_Medium = c(0.180, 0.300, 0.070, 0.050, 0.190),
  Power_High = c(0.140, 0.260, 0.010, 0.020, 0.110),
  Type1Error_Low = c(0.0010, 0.0125, 0.0010, 0.0010, 0.0010),
  Type1Error_Medium = c(0.0010, 0.0125, 0.0010, 0.0010, 0.0010),
  Type1Error_High = c(0.0010, 0.0125, 0.0010, 0.0010, 0.0010)
)

# sample 300 m3 
performance_data_300 <- data.frame(
  Method = c("MedTest", "MODIMA", "CCMM", "LDM", "PERMANOVA"),
  Power_Low = c(0.250, 0.620, 0.380, 0.390, 0.740),
  Power_Medium = c(0.200, 0.620, 0.110, 0.120, 0.510),
  Power_High = c(0.160, 0.720, 0.050, 0.060, 0.410),
  Type1Error_Low = c(0.0010, 0.0010, 0.0010, 0.0010, 0.0010),
  Type1Error_Medium = c(0.0010, 0.0010, 0.0010, 0.0010, 0.0010),
  Type1Error_High = c(0.0010, 0.0010, 0.0010, 0.0010, 0.0010)
)


# 生成虚拟列的名称
NM <- colnames(performance_data)[-1]

# 定义用于生成条形图的函数
bar_chart <- function(value, max_width = 100, height = "10px", fill = "forestgreen", background = "#e0e0e0") {
  width <- paste0(value * max_width, "%")
  bar <- div(style = list(background = fill, width = width, height = height, borderRadius = "4px"))
  chart <- div(style = list(background = background, flexGrow = 1, height = height, borderRadius = "4px"), bar)
  label <- div(style = list(minWidth = "20px", textAlign = "right", color = fill, marginLeft = "1px"), 
               format(round(value, digits = 2), nsmall = 2))
  tagList(
    div(style = list(display = "flex", alignItems = "center", flexDirection = "row"), chart, label)
  )
}

# 定义颜色获取函数
get_color <- function(metric, value){
  if(grepl("Power", metric)){
    if(value >= 0.7){
      col = brewer.pal(8,'Paired')[2]  # 0.7以上 蓝色
    }else if(value >= 0.4){
      col = brewer.pal(8,'Set2')[6]  # 0.4-0.7 黄色
    }else{
      col = brewer.pal(8,'RdGy')[2]  # 0-0.4 红色
    }
  } else {  # Type-1 Error
    if(value <= 0.05){
      col = brewer.pal(8,'Paired')[2]  # 0-0.05 蓝色
    }else if(value <= 0.1){
      col = brewer.pal(8,'Set2')[6]  # 0.05-0.1 黄色
    }else{
      col = brewer.pal(8,'RdGy')[2]  # 0.1以上 红色
    }
  }
  return(col)
}

# 定义coldefs_list函数
coldefs_list <- function(cols, labels){
  coldefs_list <- lapply(seq_along(cols), function(idx){
    col_name <- cols[idx]
    col_metric <- performance_data[[col_name]]
    colors <- sapply(col_metric, get_color, metric = col_name, USE.NAMES = FALSE)
    
    reactable::colDef(
      name = labels[idx],  # 使用自定义标签
      cell = function(value, index) {
        bar_chart(value = value, fill = colors[index])
      },
      style = function(value, index) {
        list(fontWeight = 600, fontSize = 10, color = colors[index])
      },
      align = "center", minWidth = 10
    )
  })
  
  names(coldefs_list) <- cols
  return(coldefs_list)
}

# 定义列和标签
power_cols <- c("Power_Low", "Power_Medium", "Power_High")
type1_cols <- c("Type1Error_Low", "Type1Error_Medium", "Type1Error_High")
labels <- c("Low", "Medium", "High")

# 应用coldefs_list函数
coldefs_list_power <- coldefs_list(power_cols, labels)
coldefs_list_type1 <- coldefs_list(type1_cols, labels)

# 创建reactable表格
performance_table <- reactable(performance_data,
                               columns = c(
                                 list(
                                   Method = colDef(minWidth = 10, align = "center", cell = function(value, index) {  # 调整列宽
                                     div(style = list(fontWeight = 600, fontSize = 15), value)
                                   })
                                 ),
                                 coldefs_list_power,
                                 coldefs_list_type1
                               ),
                               bordered = FALSE,  # 去掉表格边框
                               highlight = FALSE,  # 去掉高亮
                               pagination = FALSE,
                               resizable = FALSE,
                               wrap = FALSE,
                               style = list(fontFamily = "Helvetica", fontSize = "16px"),
                               theme = reactableTheme(
                                 headerStyle = list(
                                   "&:hover[aria-sort]" = list(background = "#f7f7f8"),
                                   "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 98%)"),
                                   borderColor = "#555"
                                 )
                               ),
                               columnGroups = list(
                                 colGroup(name = "Power", columns = power_cols),
                                 colGroup(name = "Type-1 Error", columns = type1_cols)
                               )
)

# 保存表格为临时HTML文件
html_file <- tempfile(fileext = ".html")
saveWidget(performance_table, html_file, selfcontained = TRUE)

# 截图并保存为PNG文件
png_file <- tempfile(fileext = ".png")
webshot(html_file, file = png_file, vwidth = 900, vheight = 300, zoom = 2)  # 调整图片宽度以适应表格内容

# 删除临时HTML文件
unlink(html_file)

# 打开并显示保存的图片
img <- png::readPNG(png_file)
grid::grid.raster(img)

# 删除临时PNG文件
unlink(png_file)


