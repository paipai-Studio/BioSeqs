# BioSeqs Bioconductor Extensions (MoonBit Implementation)

本扩展在 BioSeqs 项目中用 MoonBit 语言整合实现了 **4 个 Bioconductor 核心包**的功能，用于组学/转录组数据分析流水线中的常见任务：缺失值插补（impute）、方差稳定化归一化（vsn）、基因集管理（GSEABase）、PCA 分析与可视化（PCAtools）。

---

## 已实现模块

| # | 模块文件               | 对应 Bioconductor 包 | 主要用途                                       | 测试数量 | 状态  |
|---|------------------------|----------------------|----------------------------------------------|----------|-------|
| 1 | `src/impute.mbt`       | `impute`             | 缺失值插补（KNN / mean / median / LOCF / NOCB） | 13       | ✅ 通过 |
| 2 | `src/vsn.mbt`          | `vsn`                | 方差稳定化归一化（glog 变换, vsn2, mean-SD）   | 9        | ✅ 通过 |
| 3 | `src/gsea_base.mbt`    | `GSEABase`           | 基因集数据结构 & GMT/GMX 格式解析             | 20       | ✅ 通过 |
| 4 | `src/pcatools.mbt`     | `PCAtools`           | 高级 PCA 分析 (scree/biplot/outliers/correlations) | 11   | ✅ 通过 |

**总计：53 个单元测试全部通过 ✅**

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

## 测试 & 验证

### 仅运行本扩展相关的 53 个测试
```bash
export NODE_OPTIONS="--max-old-space-size=8192"
moon test test/moonbit/impute_test.mbt    # 13 passed
moon test test/moonbit/vsn_test.mbt       # 9 passed
moon test test/moonbit/gsea_base_test.mbt # 20 passed
moon test test/moonbit/pcatools_test.mbt  # 11 passed
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
```

---

## 新增文件清单

### 核心实现（src/）
```
src/impute.mbt       # impute 模块（~220 行）
src/vsn.mbt          # vsn 模块（~480 行）
src/gsea_base.mbt    # GSEABase 模块（~560 行）
src/pcatools.mbt     # PCAtools 模块（~720 行）
```

### 单元测试（test/moonbit/）
```
test/moonbit/impute_test.mbt     # 13 cases
test/moonbit/vsn_test.mbt        # 9 cases
test/moonbit/gsea_base_test.mbt  # 20 cases
test/moonbit/pcatools_test.mbt   # 11 cases
```

### 示例（examples/）
```
examples/impute_demo/main.mbt     + moon.pkg
examples/vsn_demo/main.mbt        + moon.pkg
examples/gsea_base_demo/main.mbt  + moon.pkg
examples/pcatools_demo/main.mbt   + moon.pkg
```

### 文档
```
myREADME.md          # 本文件
```
