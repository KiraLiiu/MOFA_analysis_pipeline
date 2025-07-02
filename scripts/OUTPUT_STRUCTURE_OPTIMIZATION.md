# MOFA + MultiGSEA Pipeline - 输出文件结构优化

## 概述
本文档描述了优化后的MOFA + MultiGSEA分析管道的输出文件结构。所有输出文件现在都按照分析步骤和数据类型进行了有组织的分类。

## 主要改进

### 1. 统一的输出基础路径
- 所有输出文件都保存在 `results/mofa_analysis/` 目录下
- 每个分析步骤都有独立的子目录

### 2. 文件命名规范
- 使用数字前缀（01_, 02_, 03_, 04_, 05_）表示分析步骤
- 使用描述性的文件名
- 包含文件类型信息（.csv, .xlsx, .svg, .png, .pdf）

## 详细文件结构

```
results/mofa_analysis/
├── 01_Data_Overview_Plot.pdf                    # 数据概览图
├── 01_MOFA_Data_Structure.rds                   # MOFA数据结构
├── 01_MOFA_Object_Created.rds                   # 创建的MOFA对象
│
├── data_preparation/                            # 步骤01：数据准备
│   ├── elbow_plots/                            # 累积方差分析图
│   │   ├── Transcriptomics_cumulative_variance.pdf
│   │   ├── Proteomics_cumulative_variance.pdf
│   │   └── Metabolomics_cumulative_variance.pdf
│   │
│   ├── distribution_analysis/                   # 分布分析
│   │   ├── 01_Metabolomics_distributions_post_HVF.csv
│   │   ├── 02_Proteomics_distributions_post_HVF.csv
│   │   ├── 03_Transcriptomics_distributions_post_HVF.csv
│   │   ├── Metabolomics_distributions_post_HVF.pdf
│   │   ├── Proteomics_distributions_post_HVF.pdf
│   │   └── Transcriptomics_distributions_post_HVF.pdf
│   │
│   └── feature_selection/                       # 特征选择结果
│       ├── 01_Transcriptomics_Selected_Features.csv
│       ├── 02_Proteomics_Selected_Features.csv
│       ├── 03_Metabolomics_Selected_Features.csv
│       └── 04_Feature_Selection_Summary.csv
│
├── model_training/                              # 步骤02：模型训练
│   ├── 02_MOFA_Model_Trained.hdf5              # 训练好的MOFA模型
│   ├── 02_MOFA_Model_Trained.rds               # R格式的MOFA模型
│   ├── 02_Training_Parameters.csv              # 训练参数
│   └── 02_Model_Information.csv                # 模型信息
│
├── model_explanation/                           # 步骤03：模型解释
│   ├── 03_Analysis_Summary_Report.csv          # 分析总结报告
│   │
│   ├── figures/                                # 模型解释图
│   │   ├── 03_Factor_Correlation_Plot.svg
│   │   ├── 03_Factor_Correlation_Plot.png
│   │   ├── 03_Variance_Explained_Plot.svg
│   │   ├── 03_Variance_Explained_Plot.png
│   │   ├── 03_Factor_Covariate_Correlation.svg
│   │   └── 03_Factor_Covariate_Correlation.png
│   │
│   └── tables/                                 # 模型解释表格
│       ├── Factor1_all_features.xlsx
│       ├── Factor2_all_features.xlsx
│       ├── Factor3_all_features.xlsx
│       ├── 03_Factor_Values_All_Samples.csv
│       ├── 03_Variance_Explained_By_Factors.csv
│       ├── 03_Model_Dimensions.csv
│       └── 03_Convergence_Information.csv
│
├── multigsea_analysis/                          # 步骤04：MultiGSEA分析
│   ├── 04_MultiGSEA_Analysis_Report.csv        # 分析报告
│   │
│   ├── tables/                                 # MultiGSEA结果表格
│   │   ├── 04_Transcriptome_Data_for_MultiGSEA.csv
│   │   ├── 04_Proteome_Data_for_MultiGSEA.csv
│   │   ├── 04_Metabolome_Data_for_MultiGSEA.csv
│   │   ├── 04_MultiGSEA_Complete_Results.csv
│   │   ├── 04_MultiGSEA_Complete_Results.xlsx
│   │   ├── 04_MultiGSEA_Significant_Results.csv
│   │   ├── 04_MultiGSEA_Significant_Results.xlsx
│   │   ├── 04_MultiGSEA_Top_KEGG_Pathways.csv
│   │   ├── 04_MultiGSEA_Top_Reactome_Pathways.csv
│   │   └── 04_MultiGSEA_Summary_Statistics.csv
│   │
│   └── figures/                                # MultiGSEA图形（如果有）
│
└── visualization/                               # 步骤05：可视化
    ├── 05_Visualization_Summary.csv            # 可视化总结
    ├── 05_Generated_Files_List.csv             # 生成文件列表
    ├── 05_Visualization_Analysis_Report.csv    # 可视化分析报告
    │
    ├── mofa_figures/                           # MOFA可视化图
    │   ├── 05_MOFA_Data_Overview.svg
    │   ├── 05_MOFA_Data_Overview.png
    │   ├── 05_MOFA_Factor_Correlation.svg
    │   ├── 05_MOFA_Factor_Correlation.png
    │   ├── 05_MOFA_Variance_Explained.svg
    │   ├── 05_MOFA_Variance_Explained.png
    │   ├── 05_MOFA_Factor_Covariate_Correlation.svg
    │   └── 05_MOFA_Factor_Covariate_Correlation.png
    │
    └── multigsea_figures/                      # MultiGSEA可视化图
        └── pathway_heatmaps/                   # 通路热图
            ├── KEGG_pathway_heatmap.svg
            └── Reactome_pathway_heatmap.svg
```

## 文件命名规则

### 1. 步骤标识
- `01_` - 数据准备
- `02_` - 模型训练
- `03_` - 模型解释
- `04_` - MultiGSEA分析
- `05_` - 可视化

### 2. 文件类型标识
- `_Plot` - 图形文件
- `_Results` - 分析结果
- `_Summary` - 总结报告
- `_Statistics` - 统计信息
- `_Features` - 特征数据
- `_Parameters` - 参数设置

### 3. 数据类型标识
- `MOFA_` - MOFA相关
- `MultiGSEA_` - MultiGSEA相关
- `Transcriptomics_` - 转录组数据
- `Proteomics_` - 蛋白质组数据
- `Metabolomics_` - 代谢组数据

## 主要优势

### 1. 组织性
- 按分析步骤分类
- 按数据类型分类
- 清晰的层次结构

### 2. 可追溯性
- 数字前缀便于排序
- 描述性文件名
- 完整的分析报告

### 3. 易用性
- 统一的文件格式
- 标准化的命名
- 详细的元数据

### 4. 可扩展性
- 模块化的结构
- 易于添加新的分析步骤
- 支持多种输出格式

## 使用建议

1. **运行顺序**：按照数字前缀顺序运行脚本
2. **文件查找**：使用描述性文件名快速定位
3. **结果追踪**：查看各步骤的分析报告
4. **数据备份**：定期备份重要的中间结果

## 注意事项

1. 确保有足够的磁盘空间存储所有输出文件
2. 定期清理临时文件以节省空间
3. 保持文件命名的一致性
4. 及时更新分析报告 