# 胰腺癌代谢性细胞死亡项目交接文件

> 用途：用于在新的对话框中快速恢复项目上下文、角色设定、当前进度、关键文件、目录定位与后续衔接。  
> 项目主题：**胰腺癌代谢性细胞死亡（ferroptosis / cuproptosis / disulfidptosis）系统生物信息学研究**  
> 目标：形成一篇可投稿至 **JCR Q1** 期刊的高质量论文。

---

# 一、建议在新对话中使用的角色设定 Prompt

请将下面整段复制到新的对话框中使用：

## 角色设定 Prompt（可直接复制）

我正在进行一个关于“胰腺癌代谢性细胞死亡（ferroptosis / cuproptosis / disulfidptosis）”的系统生物信息学研究，目标是完成一篇可投稿至 JCR Q1 期刊的高质量论文。请你作为我的长期科研项目总指导，基于我上传的文件、代码、结果、图片、表格和文字记录，先核查我当前进度，再继续指导我推进研究。

你现在是：
1. 生物信息学项目总负责人
2. 医学科研设计与统计顾问
3. R/Python 数据分析与科研绘图助手
4. 论文结果整合与投稿包装顾问

你的目标不是泛泛而谈，而是：
- 审核我已经完成的内容
- 判断是否达到可投稿标准
- 发现技术漏洞、统计问题和逻辑断裂
- 给出下一步最优执行方案
- 必要时直接给出可运行代码、分析流程、绘图方案、结果解读和验收标准

你必须严格执行以下流程：

### 第1步：进度审计
请根据我上传的内容，逐项核查以下模块，并标注：
- 已完成
- 部分完成
- 未完成
- 已完成但需要返工

核查模块包括：
- 项目初始化与目录结构
- 数据下载（TCGA / GTEx / GEO / ICGC / MR / 单细胞）
- 数据预处理
- 基因集构建
- 差异分析
- 泛癌分析
- 聚类分型
- 预后模型构建
- 外部验证
- MR 分析
- 单细胞分析
- 免疫分析
- 药物敏感性 / 药物重定位
- 科研绘图
- 结果写作
- Methods 写作
- 补充材料与投稿验收

### 第2步：成果质量核查
评估我现有结果是否存在以下问题，并分为：
- 严重问题
- 中等问题
- 格式/展示问题

需重点检查：
- 数据来源不统一
- 样本纳入排除标准不清
- 临床结局定义不规范
- 批次效应处理错误
- 特征筛选存在数据泄露
- 生存模型构建不规范
- 外部验证不充分
- MR 工具变量不合格
- 单细胞注释依据不足
- 免疫分析算法单一
- 药物分析缺乏生物学解释
- 图表不符合投稿表达
- 结果与讨论逻辑断裂
- 缺少可重复性与代码管理

### 第3步：判断当前阶段
根据我上传的材料，判断我当前最准确的项目阶段。

### 第4步：制定当前最该做的下一步
只输出当前最该做的下一阶段任务，不要铺陈过多无关内容。要求给出：
1. 当前阶段目标
2. 本阶段任务清单
3. 每个任务的输入文件
4. 每个任务的分析方法
5. 预期输出文件
6. 结果验收标准
7. 常见错误与规避方法

### 第5步：一步一步带我执行
请像项目导师一样推进：
- 直接告诉我先做什么
- 给我可运行的 R 或 Python 代码
- 解释代码目的
- 告诉我运行后应该得到什么结果
- 告诉我如何判断结果是否正确
- 如果结果异常，告诉我怎么排查

工作原则：
1. 优先依据我上传的真实文件判断进度，不要假设我已经做完某些步骤
2. 如果我已经完成某一步，不要重复从头讲，要在我现有基础上继续推进
3. 如果我做错了，请明确指出，并给出返工方案
4. 如果我上传了结果图、表格、代码或报错，请优先进行结果审核与问题诊断
5. 建议必须以“可执行、可复现、可投稿”为标准
6. 同时兼顾：
   - 生物学合理性
   - 统计学规范性
   - 代码可重复性
   - 图表表达质量
   - 投稿逻辑完整性
7. 尽量给出：
   - 可运行代码
   - 文件命名建议
   - 图表设计方案
   - 表格结构
   - 结果描述模板
8. 如果信息不足，请先列出还缺哪些关键材料，再给出基于现有信息的最佳判断
9. 如果我上传多个文件，请先帮我梳理这些文件各自代表哪个分析阶段
10. 任何时候都优先避免数据泄露、训练验证混淆、结局定义错误、批次效应误处理等高风险问题

请使用中文回答，结论明确，优先给结构化审查意见。  
如果你认为我做得对，请明确说“这一步通过”；  
如果你认为需要返工，请明确说“这一步需要返工”，并说明原因。

---

# 二、当前项目进度总览（截至本次会话）

## 当前项目阶段
**已完成 TCGA-PAAD 下载、cohort 构建、表达矩阵清洗与 starter 基因集整理。**  
当前最准确阶段可定义为：

> **“TCGA 训练集已完成 cohort 构建 + 表达矩阵标准化清洗 + starter 代谢性细胞死亡基因集映射，处于 bulk 主分析正式启动前阶段。”**

---

# 三、已完成内容与审查结论

## 1. TCGA-PAAD 下载与初步审计
### 已完成内容
- RNA-seq STAR-counts 已下载
- clinical 已下载
- mutation MAF 已下载成功
- expression-clinical mapping 已完成
- sample audit 已完成

### 关键结果
- 表达样本数：**183**
- 唯一患者数：**178**
- Primary Tumor：**178**
- Solid Tissue Normal：**4**
- Metastatic：**1**
- clinical 患者数：**185**
- mapped samples：**183**
- unmapped samples：**0**
- 选用 assay：**unstranded**

### 审查结论
**这一步通过。**

---

## 2. TCGA-PAAD 主分析 cohort 构建
### 已完成内容
- all-TP 主队列构建完成
- strict-PDAC 敏感性分析队列构建完成
- OS 初步标准化完成

### 关键结果
- all-TP 主队列：**177**
- strict-PDAC 队列：**168**
- all-TP OS events：**93**
- strict-PDAC OS events：**93**
- OS 时间范围：**4–2741 天**

### 当前建议
- **主分析训练集：177 个 all-TP**
- **strict-PDAC：作为敏感性分析**

### 审查结论
**这一步通过。**

---

## 3. 表达矩阵与 cohort 一一对应
### 已完成内容
- 原始表达矩阵与 cohort 匹配
- strict-PDAC 表达矩阵匹配
- ordered cohort 文件生成

### 关键结果
- 原始表达矩阵：**60660 × 183**
- all-TP：**60660 × 177**
- strict-PDAC：**60660 × 168**
- 样本顺序已与 cohort 对齐
- 没有样本丢失

### 审查结论
**这一步通过。**

---

## 4. 基因注释清洗与低表达过滤
### 已完成内容
- Ensembl version 清洗
- gene symbol 映射
- duplicated symbol 处理
- 低表达过滤完成

### 关键结果
- 原始基因行数：**60660**
- valid symbol rows：**60660**
- deduplicated symbol rows：**59427**
- all-TP 过滤后矩阵：**23128 × 177**
- strict-PDAC 过滤后矩阵：**23123 × 168**

### 审查结论
**这一步通过。**

### 已知注意事项
- 当前矩阵不是仅 protein-coding，可能包含 rRNA、pseudogene、lncRNA 等
- 后续具体分析时需明确 gene_type 策略

---

## 5. starter 版代谢性细胞死亡基因集整理与映射
### 已完成内容
采用“**主分析集 + 扩展集**”双轨制 starter 框架。

### starter 结果
- ferroptosis：**31 个**（主分析 19，扩展 12）
- cuproptosis：**16 个**（主分析 9，扩展 7）
- disulfidptosis：**12 个**（主分析 7，扩展 5）

### 映射结果
- 三类基因在 all-TP filtered matrix 中全部检测到
- 三类基因在 strict-PDAC filtered matrix 中全部检测到
- `missing_in_tcga.tsv` 为空
- 当前 main detected 主分析集共 **35 条记录**

### 已知注意事项
- **SLC7A11** 同时属于 ferroptosis 与 disulfidptosis，后续需要统一重复基因使用规则
- 当前 starter gene set 可作为主分析起点，但**不是最终投稿前的终版 Supplementary Table**
- 当前整体评估：**方向正确，不需要推翻重做**

### 审查结论
**这一步通过。**

---

# 四、当前模块状态总表

| 模块 | 当前状态 | 备注 |
|---|---|---|
| 项目初始化与目录结构 | 部分完成 | 基本目录已建立，但还可继续规范 |
| 数据下载（TCGA） | 已完成 | 已完成 bulk 主训练集基础下载 |
| 数据下载（GTEx / GEO / ICGC / MR / 单细胞） | 未完成 | 尚未系统推进 |
| 数据预处理 | 已完成（TCGA bulk 主线） | cohort、矩阵清洗、低表达过滤已完成 |
| 基因集构建 | 部分完成 | starter 版已完成，终版证据表未完成 |
| 差异分析 | 未完成 | 尚未正式开始 |
| 泛癌分析 | 未完成 | 尚未开始 |
| 聚类分型 | 未完成 | 尚未开始 |
| 预后模型构建 | 未完成 | 尚未开始 |
| 外部验证 | 未完成 | 尚未开始 |
| MR 分析 | 未完成 | 尚未开始 |
| 单细胞分析 | 未完成 | 尚未开始 |
| 免疫分析 | 未完成 | 尚未开始 |
| 药物敏感性 / 药物重定位 | 未完成 | 尚未开始 |
| 科研绘图 | 部分完成 | 当前以审计表为主，正式图尚未开始 |
| Results 写作 | 未完成 | 尚未开始 |
| Methods 写作 | 部分完成 | 雏形可根据脚本反推，但未系统成文 |
| 补充材料与投稿验收 | 未完成 | 尚未开始 |

---

# 五、推荐的项目目录（建议长期采用）

> 说明：下面目录是**建议的正式项目目录**。  
> 其中部分目录你已经建立，部分目录建议后续补齐。  
> 根目录示例：`/Users/wmz_mac/Desktop/胰腺癌`

```text
胰腺癌/
├── 00_raw_data/
│   ├── TCGA_PAAD/
│   │   ├── expression/
│   │   ├── clinical/
│   │   └── mutation/
│   ├── GTEx/
│   ├── GEO/
│   ├── ICGC/
│   ├── MR/
│   └── scRNA/
├── 01_metadata/
├── 02_processed_data/
│   ├── TCGA_PAAD/
│   ├── GEO_validation/
│   ├── ICGC_validation/
│   └── scRNA/
├── 03_gene_sets/
├── 04_scripts/
│   ├── 01_download/
│   ├── 02_preprocess/
│   ├── 03_genesets/
│   ├── 04_bulk_analysis/
│   ├── 05_external_validation/
│   ├── 06_MR/
│   ├── 07_single_cell/
│   ├── 08_immune_drug/
│   └── 09_figures_tables/
├── 05_results_tables/
├── 06_figures/
├── 07_models/
├── 08_external_validation/
├── 09_MR/
├── 10_single_cell/
├── 11_logs/
├── 12_manuscript/
│   ├── figures_legends/
│   ├── results_drafts/
│   ├── methods_drafts/
│   └── supplementary_materials/
└── README.md
```

---

# 六、当前已确认生成过的关键文件（按分析阶段整理）

> 说明：以下文件名均为本次会话中已明确出现、并已用于审查的文件名。  
> 若你本地路径与这里不同，请以你本地真实路径为准。

## A. 下载与初步审计阶段
建议位置：
- `00_raw_data/TCGA_PAAD/...`
- `01_metadata/...`
- `11_logs/...`

### 已确认文件
- `tcga_paad_sample_audit.tsv`
- `tcga_paad_sample_type_table.tsv`
- `tcga_paad_clinical_columns.txt`
- `tcga_paad_coldata_from_se.tsv`
- `tcga_paad_sample_clinical_map.tsv`

### 在脚本里曾生成/提及的重要文件名
- `tcga_paad_star_counts_se.rds`
- `tcga_paad_counts_matrix.tsv.gz`
- `tcga_paad_counts_matrix_selected_assay.rds`
- `tcga_paad_assay_names.txt`
- `tcga_paad_gene_annotation.tsv`
- `tcga_paad_clinical_raw.tsv`
- `tcga_paad_clinical_head5.tsv`
- `tcga_paad_coldata_head5.tsv`
- `phase2_tcga_download.log`
- `phase2_sessionInfo.txt`
- `tcga_paad_maf_raw.tsv.gz`

---

## B. cohort 构建阶段
建议位置：
- `02_processed_data/TCGA_PAAD/`

### 已确认文件
- `tcga_paad_cohort_build_audit.tsv`
- `tcga_paad_cohort_master.tsv`
- `tcga_paad_cohort_master_strict_pdac.tsv`
- `tcga_paad_primary_diagnosis_distribution.tsv`

### 后续排序后队列文件
- `tcga_paad_cohort_master_ordered.tsv`

---

## C. 表达矩阵匹配阶段
建议位置：
- `02_processed_data/TCGA_PAAD/`

### 已确认文件
- `tcga_paad_expression_prepare_audit.tsv`
- `tcga_paad_expr_all_tp_counts.rds`
- `tcga_paad_expr_strict_pdac_counts.rds`

---

## D. 基因注释清洗与过滤阶段
建议位置：
- `02_processed_data/TCGA_PAAD/`

### 已确认文件
- `tcga_paad_gene_cleaning_audit.tsv`
- `tcga_paad_gene_annotation_deduplicated.tsv`
- `tcga_paad_expr_all_tp_symbol_counts_filtered.rds`
- `tcga_paad_expr_strict_pdac_symbol_counts_filtered.rds`

### 在脚本中曾导出的相关文件名
- `tcga_paad_expr_all_tp_symbol_counts_unfiltered.rds`
- `tcga_paad_expr_strict_pdac_symbol_counts_unfiltered.rds`
- `tcga_paad_gene_annotation_all_tp_filtered.tsv`
- `tcga_paad_gene_annotation_strict_pdac_filtered.tsv`

---

## E. starter 基因集整理与映射阶段
建议位置：
- `03_gene_sets/`

### 已确认文件
- `metabolic_cell_death_gene_sets_starter_audit.tsv`
- `metabolic_cell_death_gene_sets_mapping_audit.tsv`
- `metabolic_cell_death_gene_sets_main_detected.tsv`
- `metabolic_cell_death_gene_sets_missing_in_tcga.tsv`

### 在脚本中曾导出的相关文件名
- `ferroptosis_genes_starter.csv`
- `cuproptosis_genes_starter.csv`
- `disulfidptosis_genes_starter.csv`
- `metabolic_cell_death_gene_sets_master_starter.tsv`
- `metabolic_cell_death_gene_sets_master_starter_mapped.tsv`
- `metabolic_cell_death_gene_sets_extended_detected.tsv`

---

# 七、当前最重要的项目约束与已知规则

## 1. 主训练集选择
- **主分析训练集：177 个 all-TP**
- **strict-PDAC：作为敏感性分析**

## 2. 当前主终点
- 目前以 **OS** 作为最稳妥主终点

## 3. 当前基因集策略
- 使用 **starter 主分析集 + 扩展集** 双轨制
- 当前 starter 集可用于正式 bulk 主分析起步
- 但不是最终 Supplementary Table 终版

## 4. 已知需要持续注意的问题
- **SLC7A11** 跨死亡类型重复，需要统一使用规则
- 当前 filtered matrix 并非仅 protein-coding，后续分析应明确 gene_type 策略
- 任何后续建模都要避免：
  - 数据泄露
  - 训练验证混淆
  - 终点定义错误
  - 批次效应误处理

---

# 八、建议在新对话里同时上传的优先文件

如果你希望新对话迅速接上当前工作，建议优先上传：

## 若继续 bulk 主分析
- `tcga_paad_expr_all_tp_symbol_counts_filtered.rds`
- `tcga_paad_cohort_master_ordered.tsv`
- `metabolic_cell_death_gene_sets_main_detected.tsv`
- 你最新要运行的脚本或报错文件

## 若需要让新对话先核查全局进度
- 本文件：`PAAD_metabolic_cell_death_project_handover.md`
- `tcga_paad_gene_cleaning_audit.tsv`
- `metabolic_cell_death_gene_sets_mapping_audit.tsv`
- `tcga_paad_cohort_build_audit.tsv`

---

# 九、适合放在新对话开头的简短说明

如果你不想每次都贴很长背景，可以在新对话开头先发这段：

我正在做胰腺癌代谢性细胞死亡（ferroptosis / cuproptosis / disulfidptosis）系统生物信息学研究。  
目前已完成：TCGA-PAAD 下载、cohort 构建、表达矩阵清洗、starter 基因集映射。  
主训练集为 177 个 all-TP，strict-PDAC 168 个作为敏感性分析。  
请先根据我上传的文件做进度核查和质量审查，再继续指导，不要从头讲背景。

---

# 十、当前整体评估

截至本次会话，项目整体评价为：

- **数据基础：良好**
- **队列定义：基本清楚**
- **表达矩阵质量：良好**
- **starter 基因集策略：合理**
- **投稿可行性：有潜力**
- **当前最关键优点：路径清晰，没有明显推倒重来的必要**

结论：
> **当前阶段的 bulk 主线基础已经建立完成，项目具备继续向分型、预后筛选与后续验证推进的条件。**
