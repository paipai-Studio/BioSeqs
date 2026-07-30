# BioSeqs Bioconductor Extensions (MoonBit Implementation)

本扩展在 BioSeqs 项目中用 MoonBit 语言整合实现了 **8 个 Bioconductor 核心包**的功能，用于组学/转录组数据分析流水线中的常见任务：缺失值插补（impute）、方差稳定化归一化（vsn）、基因集管理（GSEABase）、PCA 分析与可视化（PCAtools）、生物数据常量（Bio.Data）、序列近似匹配（Bio.Seq.Approximate）、灵活比对（Bio.Pairwise2）、化合物数据结构（Bio.Compound）。

---

## 已实现模块

| # | 模块文件               | 对应 Bioconductor 包 | 主要用途                                       | 测试数量 | 状态  |
|---|------------------------|----------------------|----------------------------------------------|----------|-------|
| 1 | `src/impute.mbt`       | `impute`             | 缺失值插补（KNN / mean / median / LOCF / NOCB） | 13       | ✅ 通过 |
| 2 | `src/vsn.mbt`          | `vsn`                | 方差稳定化归一化（glog 变换, vsn2, mean-SD）   | 9        | ✅ 通过 |
| 3 | `src/gsea_base.mbt`    | `GSEABase`           | 基因集数据结构 & GMT/GMX 格式解析             | 20       | ✅ 通过 |
| 4 | `src/pcatools.mbt`     | `PCAtools`           | 高级 PCA 分析 (scree/biplot/outliers/correlations) | 11   | ✅ 通过 |
| 5 | `src/data.mbt`         | `Bio.Data`           | IUPAC 碱基、氨基酸映射、密码子表、反向互补     | 15       | ✅ 通过 |
| 6 | `src/seq_approx.mbt`   | `Bio.Seq.Approximate` | 近似字符串匹配（Levenshtein, 错配/插入/缺失）  | 21       | ✅ 通过 |
| 7 | `src/pairwise2.mbt`    | `Bio.Pairwise2`      | 灵活双序列比对（Needleman-Wunsch / Smith-Waterman） | 16     | ✅ 通过 |
| 8 | `src/compound.mbt`     | `Bio.Compound`       | 化合物/反应/代谢通路数据结构、分子式解析        | 22       | ✅ 通过 |

**总计：117 个单元测试全部通过 ✅**

---

## 1. impute.mbt — 缺失值插补（Bioconductor impute）

提供多种微阵列/RNA-seq 表达矩阵常用的缺失值插补方法。

### Structs
- `KNNImputeParam` — KNN 插补超参数 (k, by_row, eps, min_value, max_value)
- `ImputeNAStats` — NA 统计汇总 (total_rows/cols, total_na, na_fraction, row/col_na_counts ...)

### Functions
| 函数签名                              | 说明                                    |
|-------------------------------------|-----------------------------------------|
| `impute_by_row_mean(mat)`           | 按行均值填充                             |
| `impute_by_col_median(mat)`         | 按列中位数填充                           |
| `impute_by_knn(mat, param)`         | KNN 近邻填充（基于欧氏距离加权）           |
| `impute_locf(mat, by_row=true)`     | Last-Observation-Carried-Forward         |
| `impute_nocb(mat, by_row=true)`     | Next-Observation-Carried-Backward        |
| `impute_na_by_zero(mat)`            | 以 0 填充（常用于稀疏矩阵）                |
| `impute_na_stats(mat) -> ImputeNAStats` | 统计缺失值分布                     |
| `impute_na_summary(stats) -> String`    | 人类可读的 NA 汇总报告              |

### Example
```bash
moon run examples/impute_demo
```

---

## 2. vsn.mbt — 方差稳定化归一化（Bioconductor vsn）

实现 `vsn2` / `glog`（广义对数变换），用于芯片/测序计数的方差稳定化，解决 mean-variance 依赖问题。

### Structs
- `VSNColParam`    — 每列拟合的 glog 参数 (a, b)
- `VSNResult`      — vsn 拟合结果，含所有列的参数数组和原始列统计
- `VSNControl`     — 拟合超参数（最大迭代次数、收敛阈值等）
- `MeanSDBin`      — mean-SD 分箱图的桶

### Functions
| 函数签名                               | 说明                                          |
|--------------------------------------|-----------------------------------------------|
| `glog(x, a, b) -> Double`            | 广义对数变换 `asinh` 形式                    |
| `glog_inv(y, a, b) -> Double`        | glog 的反函数                                  |
| `vsn2(mat) -> Array[Array[Double]]`  | 一步完成拟合+变换，返回归一化后的矩阵            |
| `vsn2_with_control(mat, ctrl)`       | 自定义控制参数版本的 vsn2                     |
| `vsn_fit_and_report(mat) -> VSNResult` | 仅拟合不做变换，返回可检查的参数对象          |
| `summarize_vsn_fit(fit) -> String`   | 拟合结果摘要                                    |
| `vsn_denoise(mat) -> Array[Array[Double]]` | vsn2 + 逆变换（glog_inv），实现近似去噪  |
| `mean_sd_bins(mat, n) -> Array[MeanSDBin]` | 按均值分位数分桶，计算 mean-SD 曲线  |
| `mean_sd_ascii(mat, n_bins=...) -> String`   | ASCII 绘制 mean-SD 诊断图                |

### Example
```bash
moon run examples/vsn_demo
```

---

## 3. gsea_base.mbt — 基因集管理（Bioconductor GSEABase）

实现 `GeneSet` / `GeneSetCollection` 数据结构，支持 **GMT** / **GMX** 格式的解析与写出，并提供基因集集合运算（overlap / Jaccard / union / intersect / setdiff）。

### Structs
- `GmtGeneSet`           — 单条基因集 (name, description, gene_ids, collection_type, organism, id)
- `GmtGeneSetCollection` — 基因集集合 (sets[], names[], size)
- `GeneSetCollectionType` — enum: GO_BP/MF/CC, KEGG, Reactome, Hallmark, Canonical, Positional, Motif, Computational, Other

### Functions
#### GmtGeneSet 方法
| 方法                                         | 说明                      |
|----------------------------------------------|---------------------------|
| `GmtGeneSet::new(name, desc, genes)`         | 新建（gene_ids 自动去重）  |
| `GmtGeneSet::with_annotation(name, desc, genes, ct, org, id)` | 带完整注释新建 |
| `GmtGeneSet::size(self) -> Int`              | 基因数量                  |
| `GmtGeneSet::has_gene(self, g) -> Bool`      | 大小写不敏感的成员判断    |

#### 集合运算（自由函数）
| 函数                                            | 返回值            |
|-------------------------------------------------|-------------------|
| `gene_set_overlap(a, b) -> Int`                | 交集大小          |
| `gene_set_overlap_coef(a, b) -> Double`        | 重叠系数 = |∩| / min(|a|,|b|) |
| `gene_set_jaccard(a, b) -> Double`             | Jaccard 指数      |
| `gene_set_union_arr(a, b) -> Array[String]`    | 并集              |
| `gene_set_intersect_arr(a, b) -> Array[String]`| 交集              |
| `gene_set_setdiff(a, b) -> Array[String]`      | 差集 a \ b        |

#### 集合 Collection 方法
| 方法                                                | 说明               |
|-----------------------------------------------------|--------------------|
| `GmtGeneSetCollection::new(sets)`                   | 从基因集数组构建   |
| `GmtGeneSetCollection::push_s(self, gs) -> Self`    | 追加一个基因集     |
| `GmtGeneSetCollection::by_name(self, name) -> Option[GmtGeneSet]` | 按名字查询 |
| `GmtGeneSetCollection::at(self, i) -> Option[GmtGeneSet]`        | 按索引查询 |
| `GmtGeneSetCollection::all_sizes(self) -> Array[Int]`            | 每条基因集的大小 |
| `GmtGeneSetCollection::summary(self) -> String`                  | 摘要 |

#### I/O
| 函数                                           | 说明                     |
|------------------------------------------------|--------------------------|
| `parse_gmt(content) -> GmtGeneSetCollection`   | 解析 GMT 格式字符串      |
| `parse_gmx(content) -> GmtGeneSetCollection`   | 解析 GMX 格式字符串      |
| `write_gmt(col) -> String`                      | 写回 GMT 格式字符串      |

### Example
```bash
moon run examples/gsea_base_demo
```

---

## 4. pcatools.mbt — 高级 PCA 分析（Bioconductor PCAtools）

基于幂迭代的 eigendecomposition，不依赖外部 BLAS，提供完整的 PCA 分析流水线：数据标准化 → 协方差矩阵 → 特征值分解 → 投影；并实现 scree plot、biplot、异常样本检测、变量-主成分相关性分析。

### Structs
- `FullPCAResult`    — 完整 PCA 结果：scores (n×k), loadings (p×k), eigenvalues, variance_explained, cumulative_variance, n_samples/n_variables/n_components, center/scale, used_scaling
- `PCAOutlierResult` — 异常检测结果：indices[]（异常索引）, scores[]（对应卡方距离）, n_tested, cutoff
- `BiplotOptions`    — Biplot 绘图参数

### Functions
| 函数签名                                                  | 说明                                              |
|---------------------------------------------------------|---------------------------------------------------|
| `pcatools_run_pca(data, n_components, scale) -> FullPCAResult` | 主入口；scale=true 时做 z-score 标准化      |
| `pcatools_summary(res) -> String`                       | 人类可读摘要（显示前若干 PC 的解释方差占比）        |
| `scree_plot_ascii(res, height=..., width_per_pc=...) -> String`    | ASCII 碎石图 |
| `pca_biplot_ascii(res, opts=...) -> String`             | ASCII biplot（样本在 PC 平面散点 + 载荷方向叠加）  |
| `find_pca_outliers(res, n_pcs=..., quantile_cutoff=...) -> PCAOutlierResult` | 基于马氏距离分位数的异常点检测 |
| `variable_correlations(data, res, which_pcs=...) -> Array[Array[Double]]` | 每个变量 × 每个 PC 的 Pearson 相关系数矩阵 |

### Example
```bash
moon run examples/pcatools_demo
```

---

## 5. data.mbt — 生物数据常量（Biopython Bio.Data）

提供 IUPAC  ambiguous 碱基映射、氨基酸缩写映射、标准密码子表及反向互补等分子生物学常用数据常量和工具函数。

### Functions
| 函数签名                                        | 说明                                          |
|-----------------------------------------------|-----------------------------------------------|
| `iupac_ambiguous_dna_map() -> Map[String, Array[String]]` | IUPAC DNA 简并碱基映射（11 种） |
| `iupac_ambiguous_rna_map() -> Map[String, Array[String]]` | IUPAC RNA 简并碱基映射（11 种） |
| `is_ambiguous_dna(base) -> Bool`              | 判断是否为 DNA 简并碱基（大小写不敏感）        |
| `is_ambiguous_rna(base) -> Bool`              | 判断是否为 RNA 简并碱基（大小写不敏感）        |
| `expand_ambiguous_dna(base) -> Array[String]` | 扩展 DNA 简并碱基为所有可能碱基                |
| `expand_ambiguous_rna(base) -> Array[String]` | 扩展 RNA 简并碱基为所有可能碱基                |
| `amino_3_to_1() -> Map[String, String]`       | 三字母 → 单字母氨基酸缩写映射（21 种）         |
| `amino_1_to_3() -> Map[String, String]`       | 单字母 → 三字母氨基酸缩写映射（21 种）         |
| `amino3_to_1(code) -> String`                | 三字母→单字母查询，未知返回 `"?"`              |
| `amino1_to_3(code) -> String`                | 单字母→三字母查询，未知返回 `"???"`            |
| `amino_acid_properties(aa) -> AminoAcidProperty` | 获取氨基酸理化性质（疏水性、电荷、等电点等） |
| `standard_codon_table() -> Map[String, String]` | 标准密码子表（64 密码子 → 氨基酸）              |
| `codon_lookup(codon) -> String`              | 密码子查询（大小写不敏感）                     |
| `reverse_codon_table() -> Map[String, Array[String]]` | 反向密码子表（氨基酸 → 密码子列表）          |
| `num_synonymous_codons(aa) -> Int`           | 统计氨基酸的同义密码子数量                     |
| `reverse_complement_iupac(seq) -> String`     | IUPAC 简并碱基反向互补                         |

### Structs
- `AminoAcidProperty` — 氨基酸性质：amino_acid, property, polarity, charge, hydropathy

### Example
```bash
moon run examples/data_demo
```

---

## 6. seq_approx.mbt — 近似序列匹配（Biopython Bio.Seq.Approximate）

实现基于动态规划的近似字符串匹配算法，支持错配、插入和缺失（indels），可用于序列数据库相似性搜索和 CRISPR off-target 检测。

### Structs
- `ApproxMatch` — 近似匹配结果：pattern, query, start, end, mismatches, score, is_match
- `ApproxWordResult` — 单词匹配结果：word, match_start, match_end, edit_distance, is_found

### Functions
| 函数签名                                        | 说明                                          |
|-----------------------------------------------|-----------------------------------------------|
| `count_mismatches(s1, s2) -> Int`            | 等长字符串错配计数（大小写不敏感）              |
| `approx_search(pattern, query, max_mm) -> ApproxMatch` | 无 indels 的近似搜索，返回最佳匹配        |
| `approx_search_with_indels(pattern, query, max_err) -> ApproxMatch` | 允许 indels 的近似搜索       |
| `levenshtein_distance(s1, s2) -> Int`         | Levenshtein 编辑距离                          |
| `approx_find_all(pattern, query, max_mm) -> Array[ApproxMatch]` | 查找所有近似匹配                        |
| `approx_word_search(word, text, max_err) -> ApproxWordResult` | 单词近似搜索                           |
| `approx_best_match(pattern, query, max_err_rate) -> ApproxMatch` | 基于错误率的最佳匹配              |
| `approx_match_new(...) -> ApproxMatch`       | 创建 ApproxMatch 实例                         |

### ApproxMatch 方法
| 方法                                    | 说明               |
|-----------------------------------------|--------------------|
| `get_pattern(self) -> String`           | 获取模式序列       |
| `get_query(self) -> String`             | 获取查询序列       |
| `get_start(self) -> Int`                | 匹配起始位置       |
| `get_end(self) -> Int`                  | 匹配结束位置       |
| `get_mismatches(self) -> Int`           | 错配数             |
| `get_score(self) -> Int`                | 匹配得分           |
| `is_match(self) -> Bool`                | 是否有效匹配       |
| `to_string(self) -> String`             | 字符串表示         |

### Example
```bash
moon run examples/seq_approx_demo
```

---

## 7. pairwise2.mbt — 灵活双序列比对（Biopython Bio.Pairwise2）

实现 Needleman-Wunsch（全局）和 Smith-Waterman（局部）比对算法，支持自定义匹配/错配评分矩阵和开放/延伸空位罚分。

### Structs
- `PairwiseAlignResult` — 比对结果：aligned_seq1, aligned_seq2, score, mode, seq1_start/end, seq2_start/end
- `PairwiseMode` — 比对模式枚举：Global, Local, GlobalMS, LocalMS

### Functions
| 函数签名                                                              | 说明                                    |
|----------------------------------------------------------------------|----------------------------------------|
| `pairwise_globalxx(seq1, seq2, match_sc, mismatch_sc, gap_open, gap_extend) -> PairwiseAlignResult` | 自定义全局比对 |
| `pairwise_localxx(seq1, seq2, match_sc, mismatch_sc, gap_open, gap_extend) -> PairwiseAlignResult`  | 自定义局部比对 |
| `pairwise_global(seq1, seq2) -> PairwiseAlignResult`                 | 便捷全局比对（默认参数）                |
| `pairwise_local(seq1, seq2) -> PairwiseAlignResult`                  | 便捷局部比对（默认参数）                |
| `pairwise_globalms(seq1, seq2, matrix, gap_open, gap_extend) -> PairwiseAlignResult` | 基于替换矩阵的全局比对 |
| `pairwise_localms(seq1, seq2, matrix, gap_open, gap_extend) -> PairwiseAlignResult`  | 基于替换矩阵的局部比对 |
| `simple_score(match_sc, mismatch_sc) -> (String, String) -> Double`   | 创建简单评分函数                        |
| `identity_score() -> (String, String) -> Double`                    | 恒等评分函数（match=1, mismatch=0）      |
| `matrix_score(matrix, default) -> (String, String) -> Double`        | 基于 Map 的替换矩阵评分函数             |
| `dna_matrix(match_sc, mismatch_sc) -> Map[String, Double]`           | 生成 DNA 配对评分矩阵                   |
| `alignment_summary(result) -> String`                                | 比对结果摘要（Identity/Similarity 等）   |

### PairwiseAlignResult 方法
| 方法                                        | 说明               |
|---------------------------------------------|--------------------|
| `get_aligned_seq1(self) -> String`          | 获取比对后序列 1   |
| `get_aligned_seq2(self) -> String`          | 获取比对后序列 2   |
| `get_score(self) -> Double`                 | 获取比对得分       |
| `get_mode(self) -> String`                  | 获取比对模式       |
| `format_alignment(self, line_width) -> String` | 格式化比对结果    |

### Example
```bash
moon run examples/pairwise2_demo
```

---

## 8. compound.mbt — 化合物数据结构（Biopython Bio.Compound）

实现化合物、化学反应和代谢通路的数据结构，支持分子式解析、分子量计算以及预设模板（葡萄糖、丙酮酸、柠檬酸等）。

### Structs
- `Compound` — 化合物：id, name, formula, charge, smiles, aliases[], pathways[], reactions[]
- `ChemicalReaction` — 化学反应：id, name, direction, substrates[], products[], enzymes[]
- `CompoundPathwayMap` — 代谢通路：name, compounds[], reactions[]

### Compound 方法
| 方法                                            | 说明                    |
|-------------------------------------------------|-------------------------|
| `Compound::new(id, name)`                       | 创建空化合物             |
| `Compound::with_chemical(id, name, formula, charge, smiles)` | 创建带化学信息的化合物 |
| `get_id(self) -> String`                        | 获取 ID                 |
| `get_name(self) -> String`                      | 获取名称                |
| `get_formula(self) -> String`                   | 获取分子式              |
| `get_charge(self) -> Int`                       | 获取电荷                |
| `get_smiles(self) -> String`                    | 获取 SMILES 表示        |
| `set_formula(self, formula) -> Compound`        | 设置分子式              |
| `set_charge(self, charge) -> Compound`          | 设置电荷                |
| `set_smiles(self, smiles) -> Compound`          | 设置 SMILES             |
| `add_alias(self, alias) -> Compound`            | 添加别名（自动匹配名称）|
| `has_alias(self, alias) -> Bool`                | 检查是否有别名          |
| `add_pathway(self, pathway) -> Compound`        | 添加代谢通路            |
| `add_reaction(self, reaction) -> Compound`      | 添加反应                |
| `to_string(self) -> String`                     | 字符串表示              |

### ChemicalReaction 方法
| 方法                                            | 说明                    |
|-------------------------------------------------|-------------------------|
| `ChemicalReaction::new(id, name)`               | 创建反应                |
| `get_id(self) / get_name(self) / get_direction(self)` | 获取属性          |
| `add_substrate(self, s) / add_product(self, p)` | 添加底物/产物           |
| `add_enzyme(self, e) -> ChemicalReaction`       | 添加酶                  |
| `set_direction(self, dir) -> ChemicalReaction`  | 设置方向（forward/reverse/reversible） |
| `get_equation(self) -> String`                 | 获取反应方程             |

### 自由函数
| 函数签名                                              | 说明                              |
|-----------------------------------------------------|-----------------------------------|
| `parse_formula(formula) -> Map[String, Int]`        | 解析分子式，返回元素计数          |
| `compound_molecular_weight(formula) -> Double`      | 计算分子式的分子量                |
| `compound_glucose() -> Compound`                    | 葡萄糖模板（C0003）               |
| `compound_fructose() -> Compound`                   | 果糖模板（C0005）                 |
| `compound_pyruvate() -> Compound`                   | 丙酮酸模板（C0024, charge=-1）    |
| `compound_acetate() -> Compound`                    | 乙酸模板（C0009）                 |
| `compound_citrate() -> Compound`                    | 柠檬酸模板（C0015, charge=-3）    |

### CompoundPathwayMap 方法
| 方法                                            | 说明                    |
|-------------------------------------------------|-------------------------|
| `CompoundPathwayMap::new(name)`                 | 创建通路                |
| `add_compound(self, c) -> CompoundPathwayMap`   | 添加化合物              |
| `add_reaction(self, r) -> CompoundPathwayMap`    | 添加反应                |
| `find_compound(self, id) -> Option[Compound]`   | 按 ID 查找化合物        |
| `find_reaction(self, id) -> Option[ChemicalReaction]` | 按 ID 查找反应      |
| `summary(self) -> String`                       | 通路摘要                |

### Example
```bash
moon run examples/compound_demo
```

---

## 测试 & 验证

### 仅运行本扩展相关的 117 个测试
```bash
moon test test/moonbit/impute_test.mbt    # 13 passed
moon test test/moonbit/vsn_test.mbt       # 9 passed
moon test test/moonbit/gsea_base_test.mbt # 20 passed
moon test test/moonbit/pcatools_test.mbt  # 11 passed
moon test test/moonbit/data_test.mbt      # 15 passed
moon test test/moonbit/seq_approx_test.mbt # 21 passed
moon test test/moonbit/pairwise2_test.mbt # 16 passed
moon test test/moonbit/compound_test.mbt  # 22 passed
```

### 完整编译
```bash
moon build   # src + tests + examples 全部编译通过 ✅
```

### 运行所有示例
```bash
moon run examples/impute_demo
moon run examples/vsn_demo
moon run examples/gsea_base_demo
moon run examples/pcatools_demo
moon run examples/data_demo
moon run examples/seq_approx_demo
moon run examples/pairwise2_demo
moon run examples/compound_demo
```

---

## 新增文件清单

### 核心实现（src/）
```
src/impute.mbt       # impute 模块（~220 行）
src/vsn.mbt          # vsn 模块（~480 行）
src/gsea_base.mbt    # GSEABase 模块（~560 行）
src/pcatools.mbt     # PCAtools 模块（~720 行）
src/data.mbt         # Bio.Data 模块（~350 行）
src/seq_approx.mbt   # Bio.Seq.Approximate 模块（~430 行）
src/pairwise2.mbt    # Bio.Pairwise2 模块（~770 行）
src/compound.mbt      # Bio.Compound 模块（~480 行）
```

### 单元测试（test/moonbit/）
```
test/moonbit/impute_test.mbt     # 13 cases
test/moonbit/vsn_test.mbt        # 9 cases
test/moonbit/gsea_base_test.mbt  # 20 cases
test/moonbit/pcatools_test.mbt   # 11 cases
test/moonbit/data_test.mbt       # 15 cases
test/moonbit/seq_approx_test.mbt # 21 cases
test/moonbit/pairwise2_test.mbt  # 16 cases
test/moonbit/compound_test.mbt   # 22 cases
```

### 示例（examples/）
```
examples/impute_demo/main.mbt     + moon.pkg
examples/vsn_demo/main.mbt        + moon.pkg
examples/gsea_base_demo/main.mbt  + moon.pkg
examples/pcatools_demo/main.mbt   + moon.pkg
examples/data_demo/main.mbt       + moon.pkg
examples/seq_approx_demo/main.mbt + moon.pkg
examples/pairwise2_demo/main.mbt  + moon.pkg
examples/compound_demo/main.mbt   + moon.pkg
```

### 文档
```
myREADME.md          # 本文件
```
