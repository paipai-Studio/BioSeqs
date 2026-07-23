# BioSeqs - MoonBit 生物信息学库

> https://github.com/paipai-Studio/BioSeqs
>
> https://gitlink.org.cn/IvanAXu/BioSeqs
> 
> https://mooncakes.io/docs/IvanAXu/BioSeqs
>

## 项目概述

BioSeqs 是一个基于 **MoonBit** 语言开发的生物信息学工具库，旨在复刻主流生物信息学库（Biopython、Bioconductor、scikit-bio 等）的核心功能，并实现高效的序列组装算法。

### 目录

- [核心功能覆盖](#核心功能覆盖)
- [序列组装算法](#序列组装算法)
- [哈希与数据结构](#哈希与数据结构)
- [序列模体与群体遗传学](#序列模体与群体遗传学)
- [架构设计](#架构设计)
- [核心功能实现](#核心功能实现)
- [性能优化](#性能优化)
- [测试验证](#测试验证)

### 核心功能覆盖

| 功能类别 | 对应工具/库 | 功能范围 | 状态 |
| :--- | :--- | :--- | :---: |
| **序列处理** | Biopython `Bio.Seq` | 序列对象、互补、转录、翻译、序列特征 | ✅ |
| **序列 I/O** | Biopython `Bio.SeqIO` | FASTA/FASTQ/GenBank 解析与写入 | ✅ |
| **序列比对** | Biopython / scikit-bio | Needleman-Wunsch、Smith-Waterman、多序列比对、替换矩阵(BLOSUM/PAM) | ✅ |
| **BLAST解析** | Biopython `Bio.Blast` | BLAST结果解析、tabular/xml格式、HSP过滤、最佳匹配 | ✅ |
| **SearchIO** | Biopython `Bio.SearchIO` | 统一搜索结果模型、HMMER3解析、BLAT PSL解析、BLAST转换 | ✅ |
| **系统发育树** | Biopython `Bio.Phylo` | 树结构、Newick 解析、距离计算、可视化 | ✅ |
| **PDB 结构** | Biopython `Bio.PDB` | 原子/残基/链解析、结构操作 | ✅ |
| **SAM/BAM/VCF** | pysam | 比对文件、变异检测、基因型查询 | ✅ |
| **FASTA 索引** | pyfaidx | 快速随机访问、.fai 索引 | ✅ |
| **机器学习特征** | scikit-learn | k-mer 频率、氨基酸组成、理化性质 | ✅ |
| **Biostrings** | Bioconductor Biostrings | IUPAC 支持、RSCU、复杂度、Tm 计算、模式匹配(matchPattern/vmatchPattern)、错配和插入缺失检测、回文序列查找 | ✅ |
| **GenomicRanges** | Bioconductor GenomicRanges | GRanges、区间操作、集合运算、precede/follow、coverage计算、distance_to_nearest | ✅ |
| **DESeq2** | Bioconductor DESeq2 | 差异表达分析、size factors归一化、分散度估计、负二项GLM拟合、Wald检验、LFC收缩 | ✅ |
| **dplyr** | R dplyr | DataFrame 数据操作 | ✅ |
| **enrichplot** | Bioconductor enrichplot | 富集分析结果可视化、dotplot/barplot/heatmap/cnetplot/enrichment map | ✅ |
| **IsoformSwitchAnalyzeR** | Bioconductor IsoformSwitchAnalyzeR | 转录本异构体切换分析、PSI/DPSI/DIF值计算、功能后果预测 | ✅ |
| **VariantFiltering** | Bioconductor VariantFiltering | 变异过滤与遗传模式分析、常染色体显性/隐性、X连锁、复合杂合子模式 |
| **Bio.Alphabet** | Biopython `Bio.Alphabet` | IUPAC字母表定义、DNA/RNA/蛋白质字母表、简化字母表、空位字母表 | ✅ |
| **Bio.Statistics** | scipy/stats | 描述统计、假设检验、相关性分析（Pearson/Spearman）、置信区间、Z-score | ✅ |
| **Bio.FreqAnalysis** | Biopython `Bio.SeqUtils` | 序列频率分析、k-mer计数、密码子使用频率、GC含量、序列复杂度 | ✅ |
| **Bio.Align.analysis** | Biopython `Bio.Align.analysis` | dn/ds计算（Nei-Gojobori）、进化距离（Jukes-Cantor/Kimura）、选择压力分析 | ✅ |
| **SeqUtils 高级功能** | Biopython `Bio.SeqUtils` | GC/AT滑动窗口偏斜、ORF预测、序列相似度、Hamming距离、Levenshtein编辑距离 | ✅ |
| **Motifs 高级功能** | Biopython `Bio.motifs` | 序列Logo生成、 per-position信息含量、总信息含量、模体富集分析、Pearson相关性比较 | ✅ |
| **PDB 结构分析高级功能** | Biopython `Bio.PDB` | SASA计算（Shrake-Rupley算法）、Ramachandran质量评估、Kyte-Doolittle疏水性量表、序列属性距离矩阵 | ✅ |

### 序列组装算法

| 算法 | 对应工具 | 功能说明 | 状态 |
| :--- | :--- | :--- | :---: |
| **De Bruijn Graph** | SPAdes / Velvet | k-mer 节点、欧拉路径、序列组装 | ✅ |
| **Suffix Array & Suffix Tree** | libdivsufsort | 前缀倍增算法、LCP 数组、模式匹配、最长重复子串 | ✅ |
| **Overlap-Layout-Consensus** | Celera Assembler / ARACHNE | 重叠检测、哈密顿路径、一致性序列生成 | ✅ |
| **BWT + FM-index** | Bowtie2 / BWA | 后缀数组、BWT 变换、回溯搜索、精确模式匹配 | ✅ |

### 哈希与数据结构

| 数据结构 | 对应工具 | 功能说明 | 状态 |
| :--- | :--- | :--- | :---: |
| **Bloom Filter** | Jellyfish / khmer | k-mer 计数、成员查询、误判率估算 | ✅ |

### 扩展功能模块

| 功能类别 | 对应工具/库 | 功能说明 | 状态 |
| :--- | :--- | :--- | :---: |
| **序列模体** | Biopython `Bio.motifs` | PWM、MEME格式、模体搜索、信息含量、序列Logo生成、模体富集、Pearson相关性 | ✅ |
| **群体遗传学** | Biopython `Bio.PopGen` | 等位基因频率、FST、哈迪-温伯格检验 | ✅ |
| **edgeR** | Bioconductor edgeR | 差异表达分析、DGEList、精确检验、GLM拟合 | ✅ |
| **limma** | Bioconductor limma | 差异表达分析、线性模型拟合、经验贝叶斯、voom变换 | ✅ |
| **SummarizedExperiment** | Bioconductor SummarizedExperiment | 多维基因组数据容器、Assays、行/列操作 | ✅ |
| **IRanges** | Bioconductor IRanges | 整数区间操作、集合运算、重叠检测、findOverlaps高级类型、nearest、coverage、距离矩阵计算 | ✅ |
| **TxDb** | Bioconductor GenomicFeatures | 转录本数据库、GTF解析、基因/转录本/外显子/CDS提取、UTR/内含子计算、启动子提取 | ✅ |
| **ExPASy** | Biopython `Bio.ExPASy` | 蛋白质分析工具接口、Swiss-Prot条目解析、酶数据库查询、蛋白质参数计算（分子量、等电点、GRAVY、不稳定指数） | ✅ |
| **Prosite** | Biopython `Bio.Prosite` | 蛋白质模体数据库搜索、Prosite模式解析、模体匹配算法、模体得分计算 | ✅ |
| **PAML** | Biopython `Bio.PAML` | 分子进化分析、dN/dS计算（Nei-Gojobori方法）、Jukes-Cantor校正、密码子使用分析 | ✅ |
| **Graphics** | Biopython `Bio.Graphics` | 生物信息学可视化、序列Logo绘制、序列比对可视化、基因组特征绘图 | ✅ |
| **BSgenome** | Bioconductor BSgenome | 基因组序列数据库、染色体序列检索、子序列提取、链特异性基因提取 | ✅ |
| **biomaRt** | Bioconductor biomaRt | 基因ID映射、基因注释查询、批量查询、外部数据库映射 | ✅ |
| **RUVSeq** | Bioconductor RUVSeq | RNA-seq批次效应去除、数据标准化、log2转换、RUVg/RUVs方法 | ✅ |
| **fgsea** | Bioconductor fgsea | 快速基因集富集分析、置换检验、NES/ES计算、Leading Edge基因、BH-FDR校正 | ✅ |
| **sva** | Bioconductor sva | 替代变量分析、ComBat批次校正、经验贝叶斯方法、PCA分析 | ✅ |
| **ballgown** | Bioconductor ballgown | 转录组水平差异表达分析、FPKM计算、t检验、基因/转录本结构表示 | ✅ |
| **限制性内切酶分析** | Biopython `Bio.Restriction` | 内切酶数据库、酶切位点查找、序列酶切、片段分析 | ✅ |
| **序列聚类分析** | Biopython `Bio.Cluster` | 距离矩阵、层次聚类、Newick输出、轮廓系数 | ✅ |
| **比对格式解析** | Biopython `Bio.AlignIO` | ClustalW、FASTA、Stockholm格式解析与写入 | ✅ |
| **进化树格式解析** | Biopython `Bio.TreeIO` | Newick、NHX格式解析、树操作、修剪 | ✅ |
| **蛋白质参数分析** | Biopython `Bio.SeqUtils.ProtParam` | 分子量、不稳定指数、GRAVY评分、等电点、信号肽预测 | ✅ |
| **变异分析** | Biopython `Bio.Variation` | SNP分析、突变检测、氨基酸替换分析、BLOSUM62/Grantham矩阵 | ✅ |
| **DSSP** | Biopython `Bio.PDB.DSSP` | 二级结构预测、溶剂可及表面积(SASA/Shrake-Rupley)、Ramachandran图、结构质量评估 | ✅ |
| **多肽分析** | Biopython `Bio.PDB.Polypeptide` | 氨基酸组成、疏水性分析、跨膜区域预测、等电点计算 | ✅ |
| **基因组轨道格式** | Bioconductor rtracklayer | BED/WIG/BEDGraph/GFF解析与写入、GRanges转换 | ✅ |
| **密码子使用分析** | Biopython `Bio.codonalign` / `Bio.CodonUsage` | CAI、ENC、RSCU、GC3、CBI、Fop、最优密码子检测 | ✅ |
| **比对统计** | Biopython `Bio.Align.AlignInfo` | 一致性序列、保守位点、Shannon熵、成对序列同一性、简约信息位点 | ✅ |
| **密码子比对** | Biopython `Bio.codonalign` | 密码子替换分类、dN/dS选择压力分析（Nei-Gojobori方法）、密码子使用偏好 | ✅ |
| **Entrez数据库** | Biopython `Bio.Entrez` | NCBI数据库搜索、记录获取、PubMed/Gene/Taxonomy解析 | ✅ |
| **GenomicAlignments** | Bioconductor GenomicAlignments | GAlignments对象、coverage计算、summarizeOverlaps、pileup操作 | ✅ |
| **VariantAnnotation** | Bioconductor VariantAnnotation | 变异类型检测、变异定位、编码效应预测、变异汇总 | ✅ |
| **Affy** | Biopython `Bio.Affy` | Affymetrix芯片数据分析、RMA标准化、背景校正、分位数归一化 | ✅ |
| **SVDSuperimposer** | Biopython `Bio.PDB.SVDSuperimposer` | SVD蛋白质结构叠合、旋转矩阵、平移向量、RMSD计算 | ✅ |
| **QCPSuperimposer** | Biopython `Bio.PDB.QCPSuperimposer` | 四元数特征多项式结构叠合、高精度旋转矩阵、平移向量、RMSD计算 | ✅ |
| **ResidueDepth** | Biopython `Bio.PDB.ResidueDepth` | 残基深度计算、溶剂可及表面积(SASA)、表面/核心残基识别 | ✅ |
| **StructureAlignment** | Biopython `Bio.PDB.StructureAlignment` | 多蛋白质结构比对、动态规划比对、RMSD/TM-score计算、渐进式多结构比对 | ✅ |
| **KEGG** | Biopython `Bio.KEGG` | KEGG基因/通路/化合物/酶记录解析、通路分析 | ✅ |
| **Medline** | Biopython `Bio.Medline` | Medline/PubMed文献记录解析、APA引用格式、MeSH过滤 | ✅ |
| **GenomeInfoDb** | Bioconductor GenomeInfoDb | 基因组构建管理、染色体信息、着丝粒位置、染色体臂、标准染色体筛选 | ✅ |
| **InteractionSet** | Bioconductor InteractionSet | Hi-C染色质交互、锚点对、交互矩阵、距离分布、Top交互 | ✅ |
| **MultiAssayExperiment** | Bioconductor MultiAssayExperiment | 多组学数据协调、实验协调、样本映射、跨实验子集 | ✅ |
| **TreeConstruction** | Biopython `Bio.Phylo.TreeConstruction` | 距离矩阵建树、UPGMA/WPGMA/NJ算法、替换模型（Jukes-Cantor、Kimura） | ✅ |
| **NeighborSearch** | Biopython `Bio.PDB.NeighborSearch` | KD树空间搜索、近邻查找、半径搜索、原子对搜索 | ✅ |
| **SwissProt** | Biopython `Bio.SwissProt` | Swiss-Prot/UniProt记录解析、特征注释、参考文献、关键词 | ✅ |
| **mmCIF** | Biopython `Bio.PDB.MMCIFParser` | mmCIF格式解析、数据块、类别、原子位点提取 | ✅ |
| **Nexus** | Biopython `Bio.Nexus` | NEXUS格式解析、数据矩阵、系统发育树、距离矩阵 | ✅ |
| **EMBOSS** | EMBOSS suite | GC偏斜、AT偏斜、分子量、Tm值、ORF查找、距离计算、蛋白质参数 | ✅ |
| **ChIPseeker** | Bioconductor ChIPseeker | ChIP-seq峰注释、基因距离计算、注释分类（启动子/外显子/内含子/UTR/基因间区）、BED格式读取、peak2gene关联分析、多峰值集重叠分析(peakOverlap)、Venn图可视化、饼图可视化、结果汇总与可视化 | ✅ |
| **DOSE** | Bioconductor DOSE | 疾病本体富集分析、超几何检验、多重检验校正（adjusted p-value）、富集结果过滤与汇总 | ✅ |
| **ReactomePA** | Bioconductor ReactomePA | Reactome通路富集分析、通路按顶层术语分组、富集结果可视化、通路注释查询 | ✅ |
| **AnnotationDbi** | Bioconductor AnnotationDbi | 通用注释数据库接口、基因信息存储、ID映射、GO/KEGG注释管理 | ✅ |
| **clusterProfiler** | Bioconductor clusterProfiler | 功能富集分析统一框架、超几何检验、多重检验校正、富集结果可视化 | ✅ |
| **WGCNA** | Bioconductor WGCNA | 加权基因共表达网络分析、邻接矩阵构建、TOM相似度、模块检测 | ✅ |
| **Biobase** | Bioconductor Biobase | ExpressionSet数据结构、AnnotatedDataFrame、数据归一化、log2转换 | ✅ |
| **GEOquery** | Bioconductor GEOquery | GEO数据库数据获取、Series Matrix解析、SOFT格式解析、ExpressionSet转换 | ✅ |
| **tximport** | Bioconductor tximport | 转录本量化数据导入、Salmon quant.sf解析、基因级别汇总、低表达过滤 | ✅ |
| **AnnotationHub** | Bioconductor AnnotationHub | 中心化注释资源访问、资源搜索（类型/提供者/基因组/关键词）、资源管理与缓存 | ✅ |
| **GenomicFeatures** | Bioconductor GenomicFeatures | 基因组注释功能、Gene/Transcript/Exon数据结构、GTF解析、区域查询 | ✅ |
| **graph** | Bioconductor graph | 图数据结构、有向/无向图、最短路径算法、连通分量检测、DOT格式输出 | ✅ |
| **SpatialExperiment** | Bioconductor SpatialExperiment | 空间转录组学数据结构、空间坐标管理、图像数据存储、空间范围查询、spots过滤 | ✅ |
| **phyloseq** | Bioconductor phyloseq | 微生物组分析、OTU丰度计算、分类学过滤、相对丰度计算、稀疏化处理 | ✅ |
| **BiocParallel** | Bioconductor BiocParallel | 并行计算框架、任务分块、并行求和、均值计算、进度追踪 | ✅ |
| **ensembldb** | Bioconductor ensembldb | Ensembl注释数据库接口、基因/转录本/外显子/CDS检索、染色体过滤、biotype过滤、基因长度计算 | ✅ |
| **DropletUtils** | Bioconductor DropletUtils | 空液滴检测、barcode排序、knee点检测、emptyDrops算法、细胞过滤 | ✅ |
| **rhdf5** | Bioconductor rhdf5 | HDF5文件格式支持、数据集读写、组管理、属性操作、文件列表查看 | ✅ |
| **Matrix** | Bioconductor Matrix | 稀疏矩阵操作、CSC/CSR格式、矩阵运算（加法、乘法、转置）、行列统计、范数计算 | ✅ |
| **BiocGenerics** | Bioconductor BiocGenerics | Bioconductor通用函数、NA处理、排序、集合运算、匹配、表统计、序列生成 | ✅ |
| **scran** | Bioconductor scran | 单细胞归一化(sum_factors)、SNN图构建、Leiden聚类、差异标志物分析 | ✅ |
| **monocle3** | Bioconductor monocle3 | 单细胞轨迹分析、PCA/UMAP降维、主图学习、拟时间排序、差异表达分析、分支点检测、分支特异性差异表达 | ✅ |
| **ShortRead** | Bioconductor ShortRead | 短读序列质量控制、QA统计、adapter修剪、质量修剪、读长过滤、FastQC报告生成 | ✅ |
| **scater** | Bioconductor scater | 单细胞质量控制、QC指标计算、细胞/基因过滤、CPM/log-CPM标准化、HVG检测、PCA降维 | ✅ |
| **MAST** | Bioconductor MAST | 单细胞差异表达分析、Hurdle模型、离散/连续组分检验、BH-FDR校正、结果汇总 | ✅ |
| **GenomicFiles** | Bioconductor GenomicFiles | 分布式基因组文件处理、按区间扫描BAM/BED/VCF、批量查询、归约、覆盖度计算 | ✅ |
| **DiffBind** | Bioconductor DiffBind | ChIP-seq差异结合分析、峰值重叠、共识峰识别、TMM归一化、负二项分布检验、对比组定义(contrast)、PCA可视化、热图可视化、火山图可视化、MA图可视化 | ✅ |
| **minfi** | Bioconductor minfi | Illumina 450K/EPIC DNA甲基化分析、NOOB/Illumina/分位数/功能归一化、β/M值计算、DMP/DMR分析 | ✅ |
| **flowCore** | Bioconductor flowCore | 流式细胞术FCS文件处理、数据变换、荧光补偿、门控（矩形/多边形/椭球/四象限） | ✅ |
| **bsseq** | Bioconductor bsseq | 亚硫酸氢盐测序分析、BSmooth平滑、DMR检测、CpG合并、甲基化率计算 | ✅ |
| **SingleCellExperiment** | Bioconductor SingleCellExperiment | 单细胞实验核心数据结构、多assay管理、降维（PCA/t-SNE/UMAP）、size factors、归一化 | ✅ |
| **ComplexHeatmap** | Bioconductor ComplexHeatmap | 复杂热图可视化、行/列聚类、颜色映射、热图注释、多个热图组合、分组拆分 | ✅ |
| **GSVA** | Bioconductor GSVA | 基因集变异分析、单样本通路评分（ssGSEA/zscore/PLAGE）、富集分析、置换检验、富集图可视化(enrichment map)、表型相关性分析(phenotype correlation)、生存分析(survival analysis)、分数分布分析与可视化 | ✅ |
| **ChromVAR** | Bioconductor chromVAR | 染色质变异分析、TF motif富集、GC偏差校正、细胞聚类、变异性分析 | ✅ |
| **DelayedArray** | Bioconductor DelayedArray | 延迟计算数组、懒加载操作、分块处理、行/列聚合、子集操作 | ✅ |
| **AnnotationFilter** | Bioconductor AnnotationFilter | 基因注释过滤、染色体筛选、生物类型过滤、区域重叠检测、符号模式匹配 | ✅ |
| **scDblFinder** | Bioconductor scDblFinder | 单细胞双细胞检测、Doublet评分计算、最近邻搜索、PCA降维、细胞过滤 | ✅ |
| **Batchelor** | Bioconductor batchelor | 单细胞批次校正、rescaleBatches缩放校正、mutual nearest neighbor、fastMNN多批次校正、批次混合评分 | ✅ |
| **Seurat** | Bioconductor Seurat | 单细胞数据分析核心、LogNormalize标准化、高可变基因检测、PCA降维、图聚类、UMAP可视化、差异标志物分析、跨样本整合(FindIntegrationAnchors/IntegrateData) | ✅ |
| **ChIPseeker** | Bioconductor ChIPseeker | ChIP-seq峰值注释、基因组区域分类(启动子/外显子/内含子/UTR/基因间区)、距离TSS分布、BED格式读取、peak2gene关联分析、注释可视化、统计分析 | ✅ |
| **S4Vectors** | Bioconductor S4Vectors | Rle游程编码、DataFrame数据框、Hits匹配数据结构 | ✅ |
| **genefilter** | Bioconductor genefilter | t检验、Wilcoxon秩和检验、方差过滤、CV过滤、分位数过滤 | ✅ |
| **universalmotif** | Bioconductor universalmotif | Motif结构、共识序列计算、motif创建 | ✅ |
| **MeltingTemp** | Biopython `Bio.SeqUtils.MeltingTemp` | DNA熔解温度计算（Wallace规则、GC含量法） | ✅ |
| **PDBList** | Biopython `Bio.PDB.PDBList` | PDB结构下载管理、文件路径生成、过期PDB解析 | ✅ |
| **Application** | Biopython `Bio.Application` | 命令行工具包装、参数管理、命令构建、执行 | ✅ |
| **Align.Applications** | Biopython `Bio.Align.Applications` | 多序列比对工具包装（ClustalW、Clustal Omega、Muscle、MAFFT） | ✅ |
| **Kmer** | Biopython `Bio.Kmer` | k-mer计数、频率分析、Jaccard相似度、Hamming距离、k-mer谱 | ✅ |
| **CheckSum** | Biopython `Bio.SeqUtils.CheckSum` | 序列校验和（GCG校验和、SEGUID） | ✅ |
| **HSExposure** | Biopython `Bio.PDB.HSExposure` | 蛋白质残基半球暴露度计算 | ✅ |
| **Pathway** | Biopython `Bio.Pathway` | 生物化学通路分析（物种、反应、通路数据结构及分析函数） | ✅ |
| **topGO** | Bioconductor topGO | 拓扑GO富集分析（elim算法、weight01算法、Fisher精确检验、GO图构建） | ✅ |
| **DEXSeq** | Bioconductor DEXSeq | 差异外显子使用分析（计数归一化、统计检验、结果过滤） | ✅ |
| **metagenomeSeq** | Bioconductor metagenomeSeq | 零膨胀模型微生物组差异丰度分析（归一化、零膨胀概率计算、差异检验） | ✅ |
| **HilbertCurve** | Bioconductor HilbertCurve | Hilbert曲线坐标映射（编码/解码、距离计算、基因组线性化、网格映射） | ✅ |
| **Taxonomy** | Biopython `Bio.Taxonomy` | NCBI分类数据库解析、分类树操作、谱系查询、共同祖先计算 | ✅ |
| **GFF** | Biopython `Bio.GFF` | GFF3格式解析、基因注释特征提取、属性处理 | ✅ |
| **Phylo.Consensus** | Biopython `Bio.Phylo.Consensus` | 多数规则/严格一致性树、分裂分析、支持度计算 | ✅ |
| **ConsensusClusterPlus** | Bioconductor ConsensusClusterPlus | 共识聚类分析、癌症亚型识别、稳定性评分、最优聚类数选择 | ✅ |
| **SC3** | Bioconductor SC3 | 单细胞共识聚类、PCA降维、k-means聚类、轮廓系数评估、间隙统计量 | ✅ |
| **GENESIS** | Bioconductor GENESIS | 群体结构分析、亲属关系矩阵估计、PCA分析、遗传距离计算（欧氏/曼哈顿/IBS） | ✅ |
| **DSS** | Bioconductor DSS | 离散度收缩估计与差异分析（RNA-seq差异表达、Wald检验、BH-FDR校正、差异甲基化分析DML、差异甲基化区域DMR检测） | ✅ |
| **bamsignals** | Bioconductor bamsignals | ChIP-seq信号提取（计数模式、RPM/RPKM归一化、基因组区域信号分析、染色质状态分析） | ✅ |
| **nucleR** | Bioconductor nucleR | 核小体定位分析（信号平滑、峰值检测、核小体occupancy计算、位置比较、动态变化分析） | ✅ |
| **csaw** | Bioconductor csaw | ChIP-seq窗口差异分析（滑动窗口计数、TMM归一化、窗口过滤、负二项GLM检验、差异区域检测） | ✅ |
| **slingshot** | Bioconductor slingshot | 单细胞轨迹推断（MST构建、主曲线拟合、拟时间计算、分支检测） | ✅ |
| **SCnorm** | Bioconductor SCnorm | 单细胞RNA-seq归一化（分位数回归、深度依赖偏差校正、基因特异性归一化） | ✅ |
| **EDASeq** | Bioconductor EDASeq | RNA-seq探索性分析（GC含量归一化、基因长度校正（Loess）、样本间归一化、RPKM计算） | ✅ |
| **Bio.phenotype** | Biopython `Bio.phenotype` | 表型微阵列分析（PlateRecord/WellRecord、logistic/Gompertz生长曲线拟合、CSV/JSON解析、控制减法） | ✅ |
| **Bio.Blast.Applications** | Biopython `Bio.Blast.Applications` | BLAST命令行工具包装（8种BLAST变体、快速构建器、参数管理、命令构建、验证） | ✅ |
| **Bio.PDB.PSEA** | Biopython `Bio.PDB.PSEA` | PSEA二级结构预测（CA-CA距离、虚拟键角/二面角、H/E/C分配、三态到八态转换） | ✅ |
| **Bio.SeqIO.SffIO** | Biopython `Bio.SeqIO.SffIO` | SFF二进制格式解析（454/Roche流图数据、二进制编码/解码、质量修剪、按名称查找） | ✅ |
| **maftools** | Bioconductor maftools | 癌症基因组学MAF文件解析、突变分类(SNV/Indel)、TMB计算、突变谱分析、共现分析、oncoplot数据生成 | ✅ |
| **CNVkit** | Bioconductor CNVkit | 拷贝数变异检测、CBS（循环二元分割）分割算法、拷贝数状态判定、log2比率平滑、断点检测 | ✅ |
| **destiny** | Bioconductor destiny | 单细胞扩散映射(Diffusion Maps)降维、距离矩阵计算、高斯核构建、特征分解、扩散分量计算 | ✅ |
| **tradeSeq** | Bioconductor tradeSeq | 轨迹差异表达分析、GAM（广义可加模型）拟合、样条基函数、差异表达检验、BH-FDR校正 | ✅ |

项目致力于打造一个完整、高效的生物信息学工具库，覆盖从基础序列处理到高级序列组装的全流程。

## 架构设计

### 项目结构

```
IvanAXu/BioSeqs/
├── moon.mod                    # 模块配置
├── src/                        # 源代码
│   ├── seq.mbt                 # Seq 序列对象
│   ├── seq_record.mbt          # SeqRecord 带注释的序列记录
│   ├── seqfeature.mbt          # SeqFeature 序列特征
│   ├── seqio.mbt               # 统一序列 I/O 接口
│   ├── fasta_io.mbt            # FASTA 格式解析
│   ├── fastq_io.mbt            # FASTQ 格式解析
│   ├── genbank_io.mbt          # GenBank 格式解析
│   ├── align.mbt               # MultipleSeqAlignment 多序列比对
│   ├── alignio.mbt             # 比对文件 I/O
│   ├── clustal_io.mbt          # Clustal 格式
│   ├── phylip_io.mbt           # PHYLIP 格式
│   ├── alignment.mbt           # DNA/RNA/Protein 类型及比对算法
│   ├── align_io.mbt            # AlignIO 比对格式解析 (ClustalW/FASTA/Stockholm)
│   ├── phylo.mbt               # 系统发育树 (Clade/Tree)
│   ├── tree_io.mbt             # 进化树格式解析 (Newick、NHX格式解析、树操作)
│   ├── blast.mbt               # BLAST结果解析 (tabular/xml格式、HSP、Hit、Record、过滤)
│   ├── search_io.mbt           # SearchIO 统一搜索结果模型 (HSPFragment、HSP、Hit、QueryResult、HMMER3/BLAT解析)
│   ├── subsmat.mbt             # 替换矩阵 (BLOSUM62/45、PAM250/30、矩阵解析、分数查询)
│   ├── pdb.mbt                 # PDB 数据类型
│   ├── pdb_io.mbt              # PDB 文件 I/O
│   ├── svd_superimposer.mbt    # SVDSuperimposer SVD蛋白质结构叠合 (旋转矩阵、平移向量、RMSD计算)
│   ├── sequtils.mbt            # 序列工具函数
│   ├── seq_utils.mbt           # 序列工具函数 (分子量、GC含量、Tm值、氨基酸转换、GC/AT滑动窗口偏斜、序列相似度、Hamming距离、Levenshtein编辑距离)
│   ├── complement.mbt          # 互补碱基查找表
│   ├── codon_table.mbt         # 密码子翻译表
│   ├── codon_usage.mbt         # CodonUsage 密码子使用分析 (CAI、ENC、RSCU、GC3、CBI、Fop、最优密码子检测)
│   ├── sam.mbt                 # SAM 格式解析
│   ├── bam.mbt                 # BAM 格式解析
│   ├── bgzf.mbt                # BGZF 解压缩 (支持读取压缩的 BAM 文件)
│   ├── cram.mbt                # CRAM 格式解析 (压缩二进制序列比对格式)
│   ├── vcf.mbt                 # VCF 格式解析
│   ├── variant_annotation.mbt  # VariantAnnotation 变异注释 (变异类型检测、定位、编码效应预测、变异汇总)
│   ├── biostrings.mbt          # Biostrings 序列分析 (IUPAC、k-mer频率、RSCU、复杂度)
│   ├── genomic_ranges.mbt      # GenomicRanges 基因组区间操作 (GRanges、IRanges)
│   ├── iranges.mbt             # IRanges 整数区间操作 (集合运算、重叠检测)
│   ├── genomic_alignments.mbt  # GenomicAlignments 基因组比对分析 (GAlignments、coverage、summarizeOverlaps、pileup)
│   ├── txdb.mbt                # TxDb 转录本数据库 (GTF解析、基因/转录本/外显子/CDS提取、UTR/内含子计算)
│   ├── rtracklayer.mbt         # rtracklayer 基因组轨道格式 (BED/WIG/BEDGraph/GFF解析与写入)
│   ├── deseq2.mbt              # DESeq2 差异表达分析 (size factors归一化、分散度估计、负二项GLM拟合、Wald检验、LFC收缩)
│   ├── edger.mbt               # edgeR 差异表达分析 (DGEList、精确检验、GLM拟合)
│   ├── limma.mbt               # limma 差异表达分析 (线性模型、经验贝叶斯、voom变换、对比矩阵)
│   ├── summarized_experiment.mbt # SummarizedExperiment 多维基因组数据容器
│   ├── dplyr.mbt               # dplyr 数据操作 (DataFrame、filter、select、mutate、arrange、group_by、summarize、join)
│   ├── smith_waterman.mbt      # Smith-Waterman 局部序列比对 (动态规划、自定义打分、回溯矩阵)
│   ├── needleman_wunsch.mbt    # Needleman-Wunsch 全局序列比对 (动态规划、自定义打分、回溯矩阵)
│   ├── bloom_filter.mbt        # Bloom Filter & k-mer 计数 (哈希表、位图、误判率估算、近似去重)
│   ├── bwt_fm.mbt              # BWT + FM-index (后缀数组、BWT变换、回溯搜索、精确模式匹配)
│   ├── de_bruijn.mbt           # De Bruijn Graph (k-mer节点、欧拉路径、序列组装、图简化)
│   ├── suffix_array_tree.mbt   # Suffix Array & Suffix Tree (前缀倍增、LCP数组、模式匹配、最长重复子串)
│   ├── olc.mbt                 # Overlap-Layout-Consensus (重叠检测、哈密顿路径、一致性序列生成)
│   ├── hmm.mbt                 # Hidden Markov Model (前向/后向算法、维特比算法、Baum-Welch训练、基因预测)
│   ├── kmeans.mbt              # K-means Clustering (距离计算、K-means++初始化、聚类、轮廓系数评估)
│   ├── cluster.mbt             # 序列聚类分析 (距离矩阵、层次聚类、轮廓系数)
│   ├── motifs.mbt              # 序列模体识别 (PWM、MEME格式、模体搜索、信息含量、序列Logo生成、模体富集、Pearson相关性)
│   ├── popgen.mbt              # 群体遗传学 (等位基因频率、FST、哈迪-温伯格检验)
│   ├── restriction.mbt         # 限制性内切酶分析 (酶切位点查找、片段分析)
│   ├── protparam.mbt           # ProtParam 蛋白质参数分析 (不稳定指数、等电点、信号肽预测、二级结构倾向)
│   ├── affy.mbt                # Affy Affymetrix芯片数据分析 (RMA标准化、背景校正、分位数归一化)
│   ├── feature_extraction.mbt  # 机器学习特征提取
│   ├── faidx.mbt               # FASTA 快速索引访问 (pyfaidx)
│   ├── go_enrichment.mbt       # GOEnrichment GO功能富集分析 (超几何检验、BH校正、富集结果过滤)
│   ├── single_cell.mbt         # SingleCell 单细胞数据分析 (QC指标、Log标准化、PCA降维、高变异基因)
│   ├── kegg.mbt                # KEGG数据库解析 (基因、通路、化合物、酶记录解析、通路分析)
│   ├── medline.mbt             # Medline/PubMed解析 (文献记录、APA引用格式、MeSH过滤)
│   ├── bsgenome.mbt            # BSgenome 基因组序列数据库 (染色体序列检索、子序列提取、链特异性基因提取)
│   ├── biomart.mbt             # biomaRt 基因ID转换和注释查询 (基因ID映射、注释查询、批量查询)
│   ├── ruvseq.mbt              # RUVSeq RNA-seq批次效应去除 (数据标准化、log2转换、RUVg/RUVs方法)
│   ├── fgsea.mbt               # fgsea 快速基因集富集分析 (置换检验、NES/ES计算、Leading Edge基因、BH-FDR校正)
│   ├── sva.mbt                 # sva 替代变量分析与ComBat批次校正 (经验贝叶斯方法、PCA分析)
│   ├── ballgown.mbt            # ballgown 转录组水平差异表达分析 (FPKM计算、t检验、基因/转录本结构)
│   ├── align_info.mbt          # AlignInfo 比对统计 (一致性序列、保守位点、Shannon熵、成对序列同一性)
│   ├── codon_align.mbt         # CodonAlign 密码子比对 (密码子替换分类、dN/dS选择压力分析、密码子使用偏好)
│   ├── entrez.mbt              # Entrez NCBI数据库访问 (ESearch、EFetch、PubMed/Gene/Taxonomy解析)
│   ├── genome_info_db.mbt      # GenomeInfoDb 基因组信息管理 (染色体信息、着丝粒位置、基因组构建、染色体臂)
│   ├── interaction_set.mbt     # InteractionSet 染色质交互数据 (Hi-C交互、锚点对、交互矩阵、距离分布)
│   ├── multi_assay_experiment.mbt # MultiAssayExperiment 多组学数据协调 (实验协调、样本映射、跨实验子集)
│   ├── tree_construction.mbt   # TreeConstruction 系统发育树构建 (UPGMA/WPGMA/NJ算法、距离计算)
│   ├── neighbor_search.mbt     # NeighborSearch KD树近邻搜索 (空间搜索、半径搜索、最近邻)
│   ├── swissprot.mbt           # SwissProt 蛋白数据库解析 (UniProt/Swiss-Prot记录、特征、参考文献)
│   ├── mmcif.mbt               # mmCIF格式解析 (Bio.PDB.MMCIFParser、数据块、类别、原子位点)
│   ├── nexus.mbt               # Nexus格式解析 (Bio.Nexus、数据矩阵、系统发育树、距离矩阵)
│   ├── emboss.mbt              # EMBOSS工具接口 (GC偏斜、AT偏斜、分子量、Tm值、ORF查找、距离计算)
│   ├── chipseeker.mbt          # ChIPseeker ChIP-seq峰注释分析 (峰-基因距离计算、注释分类(启动子/外显子/内含子/UTR/基因间区)、BED格式读取、peak2gene关联分析、结果汇总与可视化)
│   ├── dose.mbt                # DOSE 疾病本体富集分析 (超几何检验、多重检验校正、结果过滤)
│   ├── reactome_pa.mbt         # ReactomePA Reactome通路分析 (通路富集、按顶层术语分组、结果可视化)
│   ├── annotation_dbi.mbt      # AnnotationDbi 通用注释数据库 (基因信息存储、ID映射、GO/KEGG注释管理)
│   ├── cluster_profiler.mbt    # clusterProfiler 功能富集分析统一框架 (超几何检验、富集结果过滤与可视化)
│   ├── wgcna.mbt               # WGCNA 加权基因共表达网络分析 (邻接矩阵、TOM相似度、模块检测)
│   ├── annotation_hub.mbt      # AnnotationHub 中心化注释资源访问 (资源搜索、类型/提供者/基因组查询、资源管理)
│   ├── genomic_features.mbt    # GenomicFeatures 基因组注释功能 (Gene/Transcript/Exon数据结构、GTF解析、区域查询)
│   ├── graph.mbt               # graph 图数据结构 (有向/无向图、最短路径、连通分量、DOT输出)
│   ├── droplet_utils.mbt       # DropletUtils 空液滴检测 (emptyDrops算法、knee点检测、细胞过滤)
│   ├── scran.mbt               # scran 单细胞归一化与聚类 (sum_factors、SNN图、Leiden聚类、标志物分析)
│   ├── monocle3.mbt            # monocle3 单细胞轨迹分析 (PCA/UMAP降维、主图学习、拟时间排序)
│   ├── short_read.mbt          # ShortRead 短读序列质量控制 (QA统计、adapter修剪、质量修剪、读长过滤、FastQC报告)
│   ├── scater.mbt              # scater 单细胞质量控制 (QC指标计算、细胞/基因过滤、标准化、HVG检测、PCA)
│   ├── mast.mbt                # MAST 单细胞差异表达分析 (Hurdle模型、离散/连续检验、BH-FDR校正)
│   ├── genomic_files.mbt       # GenomicFiles 分布式基因组文件处理 (BAM/BED/VCF扫描、区间查询、归约、覆盖度)
│   ├── diffbind.mbt            # DiffBind ChIP-seq差异结合分析 (峰值重叠、共识峰、TMM归一化、NB检验)
│   ├── minfi.mbt               # minfi DNA甲基化分析 (NOOB/Illumina/分位数/功能归一化、β/M值、DMP/DMR分析)
│   ├── flow_core.mbt           # flowCore 流式细胞术FCS处理 (数据变换、荧光补偿、矩形/多边形/椭球/四象限门控)
│   ├── bsseq.mbt               # bsseq 亚硫酸氢盐测序分析 (BSmooth平滑、DMR检测、CpG合并、甲基化率计算)
│   ├── single_cell_experiment.mbt  # SingleCellExperiment 单细胞核心容器 (多assay、PCA/tSNE/UMAP降维、size factors)
│   ├── complex_heatmap.mbt      # ComplexHeatmap 复杂热图可视化 (行/列聚类、颜色映射、热图注释、分组拆分)
│   ├── gsva.mbt                 # GSVA 基因集变异分析 (ssGSEA/zscore/PLAGE评分、富集分析、置换检验、富集图可视化、表型相关性分析、生存分析、分数分布分析)
│   ├── chromvar.mbt             # ChromVAR 染色质变异分析 (TF motif富集、GC偏差校正、细胞聚类、变异性分析)
│   ├── delayed_array.mbt        # DelayedArray 延迟计算数组 (懒加载操作、分块处理、行/列聚合、子集操作)
│   ├── annotation_filter.mbt    # AnnotationFilter 基因注释过滤 (染色体筛选、生物类型过滤、区域重叠检测、符号模式匹配)
│   ├── sc_dbl_finder.mbt        # scDblFinder 单细胞双细胞检测 (Doublet评分计算、最近邻搜索、PCA降维、细胞过滤)
│   ├── batchelor.mbt            # Batchelor 单细胞批次校正 (rescaleBatches、mutual nearest neighbor、fastMNN、批次混合评分)
│   ├── seurat.mbt               # Seurat 单细胞数据分析核心 (标准化、高可变基因、PCA、聚类、UMAP、差异表达、跨样本整合)
│   ├── chipseeker.mbt           # ChIPseeker ChIP-seq峰值注释 (基因组区域分类(启动子/外显子/内含子/UTR/基因间区)、距离TSS分布、BED格式读取、peak2gene关联分析、注释可视化、统计分析)
│   ├── variation.mbt            # Bio.Variation 变异分析 (SNP分析、突变检测、氨基酸替换分析、BLOSUM62/Grantham矩阵)
│   ├── dssp.mbt                 # Bio.PDB.DSSP 二级结构分析 (DSSP解析、二级结构预测、溶剂可及表面积、Ramachandran图)
│   ├── polypeptide.mbt          # Bio.PDB.Polypeptide 多肽分析 (氨基酸组成、疏水性分析、跨膜区域预测、等电点计算)
│   ├── s4vectors.mbt            # S4Vectors 基础数据结构 (Rle游程编码、S4DataFrame数据框、Hits匹配)
│   ├── genefilter.mbt           # genefilter 基因过滤 (t检验、Wilcoxon秩和检验、方差过滤、CV过滤、分位数过滤)
│   ├── universalmotif.mbt       # universalmotif 模体分析 (S4Motif结构、共识序列计算)
│   ├── melting_temp.mbt         # MeltingTemp DNA熔解温度计算 (Wallace规则、GC含量法)
│   ├── checksum.mbt             # CheckSum 序列校验和 (GCG校验和、SEGUID)
│   ├── hs_exposure.mbt          # HSExposure 蛋白质残基半球暴露度计算
│   ├── pathway.mbt              # Pathway 生物化学通路分析 (物种、反应、通路数据结构及分析函数)
│   ├── topgo.mbt                # topGO 拓扑GO富集分析 (elim算法、weight01算法、Fisher精确检验)
│   ├── dexseq.mbt               # DEXSeq 差异外显子使用分析 (计数归一化、统计检验、结果过滤)
│   ├── metagenomeseq.mbt        # metagenomeSeq 零膨胀模型微生物组差异丰度分析 (归一化、零膨胀概率计算)
│   ├── hilbertcurve.mbt         # HilbertCurve Hilbert曲线坐标映射 (编码/解码、距离计算、基因组线性化)
│   ├── taxonomy.mbt             # Taxonomy 分类学分析 (Taxon/TaxonomyDatabase、谱系查询、共同祖先计算)
│   ├── gff.mbt                  # GFF GFF3格式解析 (GFFFeature/GFFRecord、属性解析、特征提取)
│   ├── phylo_consensus.mbt      # Phylo.Consensus 一致性树构建 (ConsensusNode/Split、多数规则树、支持度计算)
│   ├── consensus_cluster_plus.mbt # ConsensusClusterPlus 共识聚类 (一致性矩阵、k-means聚类、稳定性评分、最优k选择)
│   ├── sc3.mbt                 # SC3 单细胞共识聚类 (PCA降维、共识矩阵、轮廓系数、间隙统计量)
│   ├── genesis.mbt             # GENESIS 群体结构分析 (亲属关系矩阵、PCA分析、遗传距离计算)
│   ├── nucle_r.mbt              # nucleR 核小体定位分析 (NucPosition/NucCallResult、信号平滑、峰值检测、位置比较)
│   ├── alphabet.mbt            # Bio.Alphabet IUPAC字母表定义 (DNA/RNA/蛋白质字母表、简化字母表、空位字母表)
│   ├── statistics.mbt          # Bio.Statistics 统计分析 (描述统计、假设检验、相关性分析)
│   ├── freq_analysis.mbt       # Bio.FreqAnalysis 序列频率分析 (k-mer计数、密码子使用频率、GC含量)
│   ├── align_analysis.mbt      # Bio.Align.analysis 进化分析 (dn/ds计算、Jukes-Cantor距离、Kimura 2-parameter距离)
│   ├── bamsignals.mbt           # bamsignals ChIP-seq信号提取 (BamsigParams/BamsigRegion/BamsigRecord、信号计数与归一化)
│   ├── dss.mbt                  # DSS 离散度收缩与差异分析 (DSSCountData/DSSDEResult/DSSDMResult、Wald检验、DML/DMR检测)
│   ├── phenotype.mbt            # Bio.phenotype 表型微阵列分析 (WellRecord/PlateRecord/PhenFitParams、logistic/Gompertz拟合、CSV/JSON解析)
│   ├── blast_applications.mbt   # Bio.Blast.Applications BLAST命令行工具包装 (8种BLAST变体、快速构建器、参数管理)
│   ├── psea.mbt                 # Bio.PDB.PSEA 二级结构预测 (PseaAtom/PseaResult、CA-CA距离、虚拟二面角、H/E/C分配、三态到八态转换)
│   ├── sff_io.mbt               # Bio.SeqIO.SffIO SFF二进制格式解析 (SffHeader/SffRead/SffFile、二进制编码/解码、质量修剪)
│   ├── csaw.mbt                 # csaw ChIP-seq窗口差异分析 (滑动窗口计数、TMM归一化、窗口过滤、负二项GLM检验、差异区域检测)
│   ├── slingshot.mbt            # slingshot 单细胞轨迹推断 (MST构建、主曲线拟合、拟时间计算、分支检测)
│   ├── scnorm.mbt               # SCnorm 单细胞RNA-seq归一化 (分位数回归、深度依赖偏差校正、基因特异性归一化)
│   ├── edaseq.mbt               # EDASeq RNA-seq探索性分析 (GC含量归一化、基因长度校正Loess、样本间归一化、RPKM计算)
│   ├── searchio.mbt             # Bio.SearchIO 统一搜索结果模型 (BLAST/HMMER解析、QueryResult/Hit/HSP层次结构、E-value过滤)
│   ├── pdb_vectors.mbt          # Bio.PDB.vectors 3D向量与旋转矩阵 (Vector3/RotationMatrix3、叉积、Kabsch叠合、二面角)
│   ├── circ_seq.mbt             # Bio.SeqUtils.CircSeq 环状DNA序列操作 (circ_subseq、circ_rotate、circ_digest、限制性酶切分析)
│   ├── align_abstract.mbt       # Bio.Align.AlignAbstract 抽象比对类型 (AbstractAlignment、一致性序列、Shannon熵、简约信息位点)
│   ├── maftools.mbt             # maftools 癌症基因组学MAF解析与突变分析 (MAFMutation、MAFData、MutationSpectrum、TMB计算、共现分析)
│   ├── cnvkit.mbt               # CNVkit 拷贝数变异检测与CBS分割 (CNVProbe、CNVSegment、CBS分割算法、拷贝数状态判定)
│   ├── destiny.mbt              # destiny 单细胞扩散映射降维 (CellData、距离矩阵、高斯核、特征分解、扩散分量)
│   ├── tradeseq.mbt             # tradeSeq 轨迹差异表达分析 (TrajectoryPoint、GAM拟合、样条基函数、差异检验)
│   └── utils.mbt               # 通用工具函数
├── examples/                   # 示例程序
│   ├── affy_demo/              # Affy Affymetrix芯片数据分析示例 (RMA标准化、背景校正、分位数归一化)
│   ├── align_info_demo/        # AlignInfo 比对统计示例 (一致性序列、保守位点、Shannon熵、成对序列同一性)
│   ├── align_io_demo/          # AlignIO 比对格式解析示例 (ClustalW、FASTA、Stockholm)
│   ├── alignment_demo/         # 序列比对示例
│   ├── ballgown_demo/          # ballgown 转录组水平差异表达分析示例 (FPKM计算、t检验)
│   ├── bam_demo/               # BAM/BGZF 解析示例
│   ├── basic_seq/              # 基础序列操作示例
│   ├── biobase_demo/           # Biobase ExpressionSet示例 (多维基因组数据容器、AnnotatedDataFrame、数据归一化、log2转换)
│   ├── bioconductor_demo/      # Bioconductor模块综合示例 (ChIPseeker峰注释(外显子/内含子/UTR分类、peak2gene关联)、DOSE疾病富集、ReactomePA通路分析、AnnotationDbi注释数据库、clusterProfiler富集框架、WGCNA共表达网络、Batchelor单细胞批次校正、Seurat单细胞分析)
│   ├── biomart_demo/           # biomaRt 基因ID转换和注释查询示例 (基因ID映射、注释查询、批量查询)
│   ├── biostrings_demo/        # Biostrings 序列分析示例
│   ├── blast_demo/             # BLAST结果解析示例 (tabular/xml格式、HSP过滤、E-value/identity过滤、最佳匹配)
│   ├── bloom_filter_demo/      # Bloom Filter & k-mer 计数示例
│   ├── bsgenome_demo/          # BSgenome 基因组序列数据库示例 (染色体序列检索、子序列提取、链特异性基因提取)
│   ├── bwt_fm_demo/            # BWT + FM-index 示例
│   ├── cluster_demo/           # 序列聚类分析示例 (层次聚类、距离矩阵、轮廓系数)
│   ├── codon_align_demo/       # CodonAlign 密码子比对示例 (密码子替换分类、dN/dS选择压力分析、密码子使用偏好)
│   ├── codon_usage_demo/       # CodonUsage 密码子使用分析示例 (CAI、ENC、RSCU、GC3、CBI、Fop、最优密码子检测)
│   ├── cram_demo/              # CRAM 格式解析示例 (压缩二进制序列比对格式、CRAM转BAM、参考序列管理)
│   ├── de_bruijn_demo/         # De Bruijn Graph 序列组装示例
│   ├── deseq2_demo/            # DESeq2 差异表达分析示例
│   ├── dplyr_demo/             # dplyr 数据操作示例
│   ├── edger_demo/             # edgeR 差异表达分析示例
│   ├── emboss_demo/            # EMBOSS工具接口示例 (GC偏斜、AT偏斜、分子量、ORF查找)
│   ├── entrez_demo/            # Entrez NCBI数据库访问示例 (ESearch、EFetch、PubMed/Gene/Taxonomy解析)
│   ├── faidx_demo/             # FASTA 索引示例
│   ├── fgsea_demo/             # fgsea 快速基因集富集分析示例 (置换检验、NES/ES计算、Leading Edge基因)
│   ├── genome_info_db_demo/    # GenomeInfoDb 基因组信息管理示例 (染色体信息、着丝粒位置、染色体臂)
│   ├── genomic_alignments_demo/ # GenomicAlignments 基因组比对分析示例 (GAlignments、coverage、summarizeOverlaps、pileup)
│   ├── genomic_ranges_demo/    # GenomicRanges 基因组区间操作示例
│   ├── geoquery_demo/          # GEOquery GEO数据库示例 (Series Matrix解析、SOFT格式解析、ExpressionSet转换、基因过滤)
│   ├── go_enrichment_demo/     # GOEnrichment GO功能富集分析示例 (超几何检验、BH校正、富集结果过滤)
│   ├── hmm_demo/               # Hidden Markov Model 基因预测示例
│   ├── interaction_set_demo/   # InteractionSet 染色质交互示例 (Hi-C交互、锚点对、交互矩阵、距离分布)
│   ├── iranges_demo/           # IRanges 区间操作示例
│   ├── kegg_demo/              # KEGG数据库解析示例 (基因、通路、化合物、酶记录)
│   ├── kmeans_demo/            # K-means Clustering 聚类分析示例
│   ├── limma_demo/             # limma 差异表达分析示例 (线性模型、经验贝叶斯、voom变换)
│   ├── medline_demo/           # Medline/PubMed解析示例 (文献记录、APA引用、MeSH过滤)
│   ├── ml_features/            # 机器学习特征提取示例
│   ├── mmcif_demo/             # mmCIF格式解析示例 (数据块解析、类别查询、原子位点提取)
│   ├── motifs_demo/            # 序列模体识别示例
│   ├── multi_assay_experiment_demo/ # MultiAssayExperiment 多组学数据协调示例 (实验协调、样本映射)
│   ├── needleman_wunsch_demo/  # Needleman-Wunsch 全局序列比对示例
│   ├── neighbor_search_demo/   # NeighborSearch KD树近邻搜索示例 (半径搜索、最近邻、原子对搜索)
│   ├── nexus_demo/             # Nexus格式解析示例 (数据矩阵、系统发育树、距离矩阵)
│   ├── olc_demo/               # Overlap-Layout-Consensus 序列组装示例
│   ├── pdb_analysis_demo/      # PDB 高级结构分析示例 (主链二面角、氢键检测、二级结构分配、Ramachandran图、SASA计算、疏水性分析)
│   ├── pdb_demo/               # PDB 结构解析示例
│   ├── phylo_demo/             # 系统发育树示例
│   ├── phylo_metrics_demo/     # 系统发育树高级分析示例 (Colless平衡指数、Robinson-Foulds距离、距离矩阵)
│   ├── popgen_demo/            # 群体遗传学分析示例
│   ├── protparam_demo/         # ProtParam 蛋白质参数分析示例 (分子量、不稳定指数、等电点、信号肽预测、二级结构倾向)
│   ├── restriction_demo/       # 限制性内切酶分析示例 (酶切位点查找、片段分析)
│   ├── rtracklayer_demo/       # rtracklayer 基因组轨道格式示例 (BED/WIG/BEDGraph/GFF解析与写入)
│   ├── ruvseq_demo/            # RUVSeq RNA-seq批次效应去除示例 (数据标准化、log2转换、RUVg/RUVs方法)
│   ├── sam_vcf_demo/           # SAM/VCF 解析示例
│   ├── search_io_demo/         # SearchIO 统一搜索结果示例 (HMMER3解析、BLAT PSL解析、BLAST转换)
│   ├── seq_complexity_demo/    # 序列复杂度分析示例 (Shannon熵、语言学复杂度、GC偏斜、混沌游戏表示)
│   ├── seqio_demo/             # 序列 I/O 示例
│   ├── seq_utils_demo/         # 序列工具函数示例
│   ├── single_cell_demo/       # SingleCell 单细胞数据分析示例 (QC指标、Log标准化、PCA降维、高变异基因)
│   ├── smith_waterman_demo/    # Smith-Waterman 局部序列比对示例
│   ├── subsmat_demo/           # 替换矩阵示例 (BLOSUM62/45、PAM250/30矩阵查询、蛋白质比对打分)
│   ├── suffix_array_tree_demo/ # Suffix Array & Suffix Tree 示例
│   ├── summarized_experiment_demo/ # SummarizedExperiment 数据容器示例
│   ├── sva_demo/               # sva 替代变量分析与ComBat批次校正示例 (经验贝叶斯方法、PCA分析)
│   ├── svd_superimposer_demo/  # SVDSuperimposer SVD蛋白质结构叠合示例 (旋转矩阵、平移向量、RMSD计算)
│   ├── swissprot_demo/         # SwissProt 蛋白数据库解析示例 (记录解析、特征提取、参考文献)
│   ├── tree_construction_demo/ # TreeConstruction 系统发育树构建示例 (UPGMA/WPGMA/NJ算法)
│   ├── txdb_demo/              # TxDb 转录本数据库示例 (GTF解析、基因/转录本/外显子/CDS/UTR/内含子提取)
│   ├── tximport_demo/          # tximport转录本量化示例 (Salmon quant.sf解析、转录本到基因汇总、ExpressionSet转换、低表达过滤)
│   ├── variant_annotation_demo/ # VariantAnnotation 变异注释示例 (变异类型检测、定位、编码效应预测、变异汇总)
│   ├── annotation_hub_demo/    # AnnotationHub 中心化注释资源访问示例 (资源搜索、类型/提供者/基因组查询)
│   ├── genomic_features_demo/  # GenomicFeatures 基因组注释示例 (GTF解析、基因/转录本/外显子查询、区域查询)
│   ├── graph_demo/             # graph 图数据结构示例 (有向/无向图构建、最短路径、连通分量、DOT输出)
│   ├── droplet_utils_demo/     # DropletUtils 空液滴检测示例 (emptyDrops算法、knee点检测、细胞过滤)
│   ├── scran_demo/             # scran 单细胞归一化与聚类示例 (sum_factors、SNN图、Leiden聚类、标志物分析)
│   ├── monocle3_demo/          # monocle3 单细胞轨迹分析示例 (PCA/UMAP降维、主图学习、拟时间排序)
│   ├── short_read_demo/        # ShortRead 短读序列质量控制示例 (QA统计、adapter修剪、质量修剪、FastQC报告)
│   ├── scater_demo/            # scater 单细胞质量控制示例 (QC指标计算、细胞/基因过滤、标准化、HVG检测、PCA)
│   ├── mast_demo/              # MAST 单细胞差异表达分析示例 (Hurdle模型、离散/连续检验、BH-FDR校正)
│   ├── genomic_files_demo/     # GenomicFiles 分布式基因组文件处理示例 (BAM/BED/VCF扫描、区间查询、归约、覆盖度)
│   ├── diffbind_demo/          # DiffBind ChIP-seq差异结合分析示例 (峰值重叠、共识峰、TMM归一化、NB检验)
│   ├── minfi_demo/             # minfi DNA甲基化分析示例 (NOOB/Illumina/分位数/功能归一化、β/M值、DMP/DMR分析)
│   ├── flow_core_demo/         # flowCore 流式细胞术示例 (数据变换、荧光补偿、矩形/多边形/椭球/四象限门控)
│   ├── bsseq_demo/             # bsseq 亚硫酸氢盐测序示例 (BSmooth平滑、DMR检测、CpG合并、甲基化率计算)
│   ├── single_cell_experiment_demo/  # SingleCellExperiment 单细胞核心容器示例 (多assay、降维、size factors)
│   ├── complex_heatmap_demo/   # ComplexHeatmap 复杂热图可视化示例 (行/列聚类、颜色映射、热图注释、分组拆分)
│   ├── gsva_demo/              # GSVA 基因集变异分析示例 (ssGSEA/zscore/PLAGE评分、富集分析、置换检验、富集图可视化、表型相关性分析、生存分析、分数分布分析)
│   ├── chromvar_demo/          # ChromVAR 染色质变异分析示例 (TF motif富集、GC偏差校正、细胞聚类、变异性分析)
│   ├── delayed_array_demo/     # DelayedArray 延迟计算数组示例 (懒加载操作、分块处理、行/列聚合、子集操作)
│   ├── annotation_filter_demo/ # AnnotationFilter 基因注释过滤示例 (染色体筛选、生物类型过滤、区域重叠检测、符号模式匹配)
│   ├── sc_dbl_finder_demo/     # scDblFinder 单细胞双细胞检测示例 (Doublet评分计算、最近邻搜索、PCA降维、细胞过滤)
│   ├── batchelor_demo/          # Batchelor 单细胞批次校正示例 (rescaleBatches、fastMNN、mutual nearest neighbor、批次混合评分)
│   ├── seurat_demo/             # Seurat 单细胞数据分析示例 (标准化、高可变基因、PCA、聚类、UMAP、差异标志物分析、跨样本整合)
│   ├── chipseeker_demo/         # ChIPseeker ChIP-seq峰值注释示例 (基因组区域分类(启动子/外显子/内含子/UTR/基因间区)、距离TSS分布、BED格式读取、peak2gene关联分析、注释可视化、统计分析)
│   ├── variation_demo/           # Bio.Variation 变异分析示例 (SNP分析、突变检测、氨基酸替换分析)
│   ├── dssp_demo/                # Bio.PDB.DSSP 二级结构分析示例 (DSSP解析、二级结构预测、溶剂可及表面积)
│   ├── polypeptide_demo/         # Bio.PDB.Polypeptide 多肽分析示例 (氨基酸组成、疏水性分析、跨膜区域预测)
│   ├── s4vectors_demo/           # S4Vectors 基础数据结构示例 (Rle游程编码、S4DataFrame数据框、Hits匹配)
│   ├── genefilter_demo/          # genefilter 基因过滤示例 (t检验、Wilcoxon秩和检验、方差过滤、CV过滤)
│   ├── universalmotif_demo/      # universalmotif 模体分析示例 (S4Motif结构、共识序列计算)
│   ├── melting_temp_demo/        # MeltingTemp DNA熔解温度计算示例 (Wallace规则、GC含量法)
│   ├── checksum_demo/            # CheckSum 序列校验和示例 (GCG校验和、SEGUID)
│   ├── hs_exposure_demo/         # HSExposure 蛋白质残基半球暴露度计算示例
│   ├── pathway_demo/             # Pathway 生物化学通路分析示例 (物种、反应、通路数据结构及分析函数)
│   ├── topgo_demo/               # topGO 拓扑GO富集分析示例 (elim算法、weight01算法、Fisher精确检验)
│   ├── dexseq_demo/              # DEXSeq 差异外显子使用分析示例 (计数归一化、统计检验、结果过滤)
│   ├── metagenomeseq_demo/       # metagenomeSeq 零膨胀模型微生物组差异丰度分析示例 (归一化、零膨胀概率计算)
│   ├── hilbertcurve_demo/        # HilbertCurve Hilbert曲线坐标映射示例 (编码/解码、距离计算、基因组线性化)
│   ├── taxonomy_demo/            # Taxonomy 分类学分析示例 (分类数据库查询、谱系追踪、共同祖先计算)
│   ├── gff_demo/                 # GFF GFF3格式解析示例 (基因结构分析、属性解析、特征提取)
│   ├── phylo_consensus_demo/     # Phylo.Consensus 一致性树构建示例 (多数规则树、分裂分析、支持度计算)
│   ├── alphabet_demo/            # Bio.Alphabet IUPAC字母表示例 (DNA/RNA/蛋白质字母表、字符验证)
│   ├── statistics_demo/          # Bio.Statistics 统计分析示例 (描述统计、假设检验、相关性分析)
│   ├── freq_analysis_demo/       # Bio.FreqAnalysis 序列频率分析示例 (k-mer计数、密码子使用频率、GC含量)
│   ├── nucle_r_demo/             # nucleR 核小体定位分析示例 (信号平滑、峰值检测、位置比较)
│   ├── align_analysis_demo/      # Bio.Align.analysis 进化分析示例 (dn/ds计算、Jukes-Cantor距离、Kimura 2-parameter距离)
│   ├── bamsignals_demo/          # bamsignals ChIP-seq信号提取示例 (信号计数、归一化、染色质状态分析)
│   ├── dss_demo/                 # DSS 离散度收缩与差异分析示例 (RNA-seq差异表达、甲基化分析、DMR检测)
│   ├── phenotype_demo/           # Bio.phenotype 表型微阵列分析示例 (生长曲线拟合、预测、控制减法、CSV/JSON解析)
│   ├── blast_applications_demo/  # Bio.Blast.Applications BLAST命令行工具示例 (8种BLAST变体、快速构建器、参数管理)
│   ├── psea_demo/                # Bio.PDB.PSEA 二级结构预测示例 (CA-CA距离、虚拟二面角、H/E/C分配、三态到八态转换)
│   ├── sff_io_demo/              # Bio.SeqIO.SffIO SFF二进制解析示例 (二进制编码/解码往返、质量修剪、按名称查找)
│   ├── csaw_demo/                # csaw ChIP-seq窗口差异分析示例 (滑动窗口、TMM归一化、差异区域检测)
│   ├── slingshot_demo/           # slingshot 单细胞轨迹推断示例 (MST构建、主曲线、拟时间计算)
│   ├── scnorm_demo/              # SCnorm 单细胞RNA-seq归一化示例 (分位数回归、深度依赖校正)
│   ├── edaseq_demo/              # EDASeq RNA-seq探索性分析示例 (GC归一化、Loess校正、RPKM计算)
│   ├── searchio_demo/            # Bio.SearchIO 统一搜索结果模型示例 (BLAST/HMMER解析、QueryResult/Hit/HSP层次结构、E-value过滤)
│   ├── pdb_vectors_demo/         # Bio.PDB.vectors 3D向量与旋转矩阵示例 (Vector3运算、Kabsch叠合、二面角计算)
│   ├── circ_seq_demo/            # Bio.SeqUtils.CircSeq 环状DNA操作示例 (酶切分析、PCR引物设计、序列旋转)
│   ├── align_abstract_demo/      # Bio.Align.AlignAbstract 抽象比对示例 (一致性序列、Shannon熵、同一性矩阵、简约信息位点)
│   ├── maftools_demo/            # maftools 癌症基因组学示例 (MAF数据创建、突变分类、TMB计算、突变谱分析)
│   ├── cnvkit_demo/              # CNVkit 拷贝数变异示例 (CBS分割、拷贝数判定、断点检测、log2平滑)
│   ├── destiny_demo/             # destiny 扩散映射示例 (距离矩阵、高斯核、特征分解、嵌入坐标)
│   ├── tradeseq_demo/            # tradeSeq 轨迹差异表达示例 (GAM拟合、基因平滑、差异检验、可视化)
├── test/
│   ├── moonbit/                # MoonBit 测试文件
│   │   ├── affy_test.mbt
│   │   ├── align_io_test.mbt
│   │   ├── alignment_test.mbt
│   │   ├── bam_test.mbt
│   │   ├── bgzf_test.mbt
│   │   ├── bio_seq_test.mbt
│   │   ├── bio_seq_wb_test.mbt
│   │   ├── biostrings_test.mbt
│   │   ├── blast_test.mbt
│   │   ├── bloom_filter_test.mbt
│   │   ├── bwt_fm_test.mbt
│   │   ├── cluster_test.mbt
│   │   ├── codon_usage_test.mbt
│   │   ├── de_bruijn_test.mbt
│   │   ├── deseq2_test.mbt
│   │   ├── dplyr_test.mbt
│   │   ├── edger_test.mbt
│   │   ├── faidx_test.mbt
│   │   ├── feature_extraction_test.mbt
│   │   ├── genomic_alignments_test.mbt
│   │   ├── genomic_ranges_test.mbt
│   │   ├── go_enrichment_test.mbt
│   │   ├── hmm_test.mbt
│   │   ├── hmm_wbtest.mbt
│   │   ├── iranges_test.mbt
│   │   ├── kmeans_test.mbt
│   │   ├── kmeans_wbtest.mbt
│   │   ├── limma_test.mbt
│   │   ├── motifs_test.mbt
│   │   ├── needleman_wunsch_test.mbt
│   │   ├── olc_test.mbt
│   │   ├── pdb_test.mbt
│   │   ├── phylo_test.mbt
│   │   ├── popgen_test.mbt
│   │   ├── protparam_test.mbt
│   │   ├── restriction_test.mbt
│   │   ├── rtracklayer_test.mbt
│   │   ├── sam_test.mbt
│   │   ├── search_io_test.mbt
│   │   ├── seqfeature_test.mbt
│   │   ├── seqio_wb_test.mbt
│   │   ├── seq_utils_test.mbt
│   │   ├── sequtils_test.mbt
│   │   ├── single_cell_test.mbt
│   │   ├── smith_waterman_test.mbt
│   │   ├── subsmat_test.mbt
│   │   ├── suffix_array_tree_test.mbt
│   │   ├── suffix_array_tree_wbtest.mbt
│   │   ├── summarized_experiment_test.mbt
│   │   ├── svd_superimposer_test.mbt
│   │   ├── tree_io_test.mbt
│   │   ├── txdb_test.mbt
│   │   ├── variant_annotation_test.mbt
│   │   ├── vcf_test.mbt
│   │   ├── kegg_test.mbt
│   │   ├── medline_test.mbt
│   │   ├── bsgenome_test.mbt
│   │   ├── biomart_test.mbt
│   │   ├── ruvseq_test.mbt
│   │   ├── fgsea_test.mbt
│   │   ├── sva_test.mbt
│   │   ├── ballgown_test.mbt
│   │   ├── align_info_test.mbt
│   │   ├── codon_align_test.mbt
│   │   ├── entrez_test.mbt
│   │   ├── genome_info_db_test.mbt
│   │   ├── interaction_set_test.mbt
│   │   ├── multi_assay_experiment_test.mbt
│   │   ├── tree_construction_test.mbt
│   │   ├── neighbor_search_test.mbt
│   │   ├── swissprot_test.mbt
│   │   ├── mmcif_test.mbt
│       ├── nexus_test.mbt
│       ├── emboss_test.mbt
│       ├── chipseeker_test.mbt
│       ├── dose_test.mbt
│       ├── reactome_pa_test.mbt
│       ├── annotation_dbi_test.mbt
│       ├── cluster_profiler_test.mbt
│       ├── wgcna_test.mbt
│       ├── annotation_hub_test.mbt
│       ├── genomic_features_test.mbt
│       ├── graph_test.mbt
│       ├── droplet_utils_test.mbt
│       ├── scran_test.mbt
│       ├── monocle3_test.mbt
│       ├── short_read_test.mbt
│       ├── scater_test.mbt
│       ├── mast_test.mbt
│       ├── genomic_files_test.mbt
│       ├── diffbind_test.mbt
│       ├── minfi_test.mbt
│       ├── flow_core_test.mbt
│       ├── bsseq_test.mbt
│       ├── single_cell_experiment_test.mbt
│       ├── complex_heatmap_test.mbt
│       ├── gsva_test.mbt
│       ├── chromvar_test.mbt
│       ├── delayed_array_test.mbt
│       ├── annotation_filter_test.mbt
│       ├── sc_dbl_finder_test.mbt
│       ├── batchelor_test.mbt
│       ├── seurat_test.mbt
│       ├── chipseeker_test.mbt
│       ├── variation_test.mbt
│       ├── dssp_test.mbt
│       ├── polypeptide_test.mbt
│       ├── s4vectors_test.mbt
│       ├── genefilter_test.mbt
│       ├── universalmotif_test.mbt
│       ├── melting_temp_test.mbt
│       ├── checksum_test.mbt
│       ├── hs_exposure_test.mbt
│       ├── pathway_test.mbt
│       ├── topgo_test.mbt
│       ├── dexseq_test.mbt
│       ├── consensus_cluster_plus_test.mbt
│       ├── sc3_test.mbt
│       └── genesis_test.mbt
│       ├── metagenomeseq_test.mbt
│       ├── hilbertcurve_test.mbt
│       ├── taxonomy_test.mbt
│       ├── gff_test.mbt
│       ├── phylo_consensus_test.mbt
│       ├── alphabet_test.mbt
│       ├── statistics_test.mbt
│       ├── freq_analysis_test.mbt
│       └── align_analysis_test.mbt
│       ├── bamsignals_test.mbt
│       ├── dss_test.mbt
│       ├── nucle_r_test.mbt
│       ├── phenotype_test.mbt
│       ├── blast_applications_test.mbt
│       ├── psea_test.mbt
│       ├── sff_io_test.mbt
│       ├── csaw_test.mbt
│       ├── slingshot_test.mbt
│       ├── scnorm_test.mbt
│       ├── edaseq_test.mbt
│       ├── searchio_test.mbt
│       ├── pdb_vectors_test.mbt
│       ├── circ_seq_test.mbt
│       ├── align_abstract_test.mbt
│       ├── maftools_test.mbt
│       ├── cnvkit_test.mbt
│       ├── destiny_test.mbt
│       └── tradeseq_test.mbt
│   └── python/                 # Python 参考测试文件
│       ├── python_reference.py
│       ├── python_seqio_reference.py
│       ├── python_alignio_reference.py
│       ├── python_sequtils_reference.py
│       ├── python_seqfeature_reference.py
│       ├── python_phylo_reference.py
│       ├── python_pdb_reference.py
│       ├── python_skbio_alignment_reference.py
│       ├── python_feature_reference.py
│       ├── python_bench.py
│       ├── skbio_pysam_compare.py
│       ├── pyfaidx_compare.py
│       ├── compare.sh
│       └── compare_seqio.sh
└── cmd/                        # 命令行工具
    ├── main/                   # Seq 测试工具
    ├── seqio_main/             # SeqIO 测试工具
    ├── alignio_main/           # AlignIO 测试工具
    └── bench/                  # 性能基准测试
```

### 样例测试
```
moon build                                              # ✅ 成功
moon test --package IvanAXu/BioSeqs/test/moonbit        # ✅ 2229 个测试全部通过
```

### 模块对照表

#### 序列处理与 I/O

| MoonBit 文件 | 对应 Python 库 | 核心功能 |
| :--- | :--- | :--- |
| `seq.mbt` | BioPython `Bio.Seq` | 序列对象、互补、转录、翻译 |
| `seq_record.mbt` | BioPython `Bio.SeqRecord` | 带注释的序列记录 |
| `seqfeature.mbt` | BioPython `Bio.SeqFeature` | 序列特征与位置 |
| `seqio.mbt` | BioPython `Bio.SeqIO` | 统一序列文件 I/O |
| `fasta_io.mbt` | BioPython `Bio.SeqIO.FastaIO` | FASTA 解析 |
| `fastq_io.mbt` | BioPython `Bio.SeqIO.QualityIO` | FASTQ 解析 |
| `genbank_io.mbt` | BioPython `Bio.SeqIO.GenBankIO` | GenBank 解析 |
| `sequtils.mbt` | BioPython `Bio.SeqUtils` | CRC32、GC 含量、Tm 计算 |
| `seq_utils.mbt` | BioPython `Bio.SeqUtils` | 序列工具函数、GC/AT滑动窗口偏斜、ORF预测、Hamming距离、Levenshtein编辑距离 |
| `circ_seq.mbt` | BioPython `Bio.SeqUtils.CircSeq` | 环状DNA序列操作、酶切分析、PCR引物设计 |
| `complement.mbt` | BioPython `Bio.Data.IUPACData` | 互补碱基表 |
| `codon_table.mbt` | BioPython `Bio.Data.CodonTable` | 密码子翻译表 |

#### 比对算法

| MoonBit 文件 | 对应 Python 库 | 核心功能 |
| :--- | :--- | :--- |
| `alignment.mbt` | scikit-bio `skbio.alignment` | Needleman-Wunsch/Smith-Waterman |
| `smith_waterman.mbt` | BioPython `Bio.Align.PairwiseAligner` | Smith-Waterman 局部比对 |
| `needleman_wunsch.mbt` | BioPython `Bio.Align.PairwiseAligner` | Needleman-Wunsch 全局比对 |
| `align.mbt` | BioPython `Bio.Align` | 多序列比对对象 |
| `alignio.mbt` | BioPython `Bio.AlignIO` | 比对文件 I/O |
| `align_io.mbt` | BioPython `Bio.AlignIO` | ClustalW/FASTA/Stockholm 解析 |
| `clustal_io.mbt` | BioPython `Bio.AlignIO.ClustalIO` | Clustal 格式 |
| `phylip_io.mbt` | BioPython `Bio.AlignIO.PhylipIO` | PHYLIP 格式 |
| `subsmat.mbt` | BioPython `Bio.SubsMat` | BLOSUM/PAM 替换矩阵 |
| `align_info.mbt` | BioPython `Bio.Align.AlignInfo` | 比对统计与一致性序列 |
| `align_abstract.mbt` | BioPython `Bio.Align.AlignAbstract` | 抽象比对类型、Shannon熵、同一性矩阵、简约信息位点 |
| `codon_align.mbt` | BioPython `Bio.codonalign` | 密码子比对与 dN/dS 分析 |
| `searchio.mbt` | BioPython `Bio.SearchIO` | 统一搜索结果模型、BLAST/HMMER解析、E-value过滤 |

#### 系统发育树

| MoonBit 文件 | 对应 Python 库 | 核心功能 |
| :--- | :--- | :--- |
| `phylo.mbt` | BioPython `Bio.Phylo` | 系统发育树 (Clade/Tree) |
| `tree_io.mbt` | BioPython `Bio.TreeIO` | Newick/NHX 格式解析 |
| `tree_construction.mbt` | BioPython `Bio.Phylo.TreeConstruction` | UPGMA/WPGMA/NJ 建树算法 |

#### 结构生物学

| MoonBit 文件 | 对应 Python 库 | 核心功能 |
| :--- | :--- | :--- |
| `pdb.mbt` | BioPython `Bio.PDB` | PDB 数据类型 |
| `pdb_io.mbt` | BioPython `Bio.PDB.PDBIO` | PDB 文件 I/O |
| `svd_superimposer.mbt` | BioPython `Bio.PDB.SVDSuperimposer` | SVD 蛋白质结构叠合 |
| `neighbor_search.mbt` | BioPython `Bio.PDB.NeighborSearch` | KD 树近邻搜索 |
| `mmcif.mbt` | BioPython `Bio.PDB.MMCIFParser` | mmCIF 格式解析 |
| `pdb_vectors.mbt` | BioPython `Bio.PDB.vectors` | 3D向量/旋转矩阵、叉积、Kabsch叠合、二面角 |
| `pdb_analysis.mbt` | BioPython `Bio.PDB.StructureAnalysis` | 二面角计算、距离矩阵、接触图、氢键检测、二级结构分配、Ramachandran图、SASA计算(Shrake-Rupley)、结构质量评估、疏水性分析 |

#### 基因组分析

| MoonBit 文件 | 对应 Python 库 | 核心功能 |
| :--- | :--- | :--- |
| `sam.mbt` | pysam | SAM 文件解析 |
| `bam.mbt` | pysam | BAM 文件解析 |
| `bgzf.mbt` | pysam | BGZF 解压缩 |
| `vcf.mbt` | pysam | VCF 文件解析 |
| `cram.mbt` | pysam | CRAM 格式解析 |
| `genomic_ranges.mbt` | Bioconductor GenomicRanges | GRanges 区间操作 |
| `iranges.mbt` | Bioconductor IRanges | 整数区间操作 |
| `genomic_alignments.mbt` | Bioconductor GenomicAlignments | GAlignments 比对分析 |
| `variant_annotation.mbt` | Bioconductor VariantAnnotation | 变异注释 |
| `txdb.mbt` | Bioconductor GenomicFeatures | TxDb 转录本数据库 |
| `rtracklayer.mbt` | Bioconductor rtracklayer | BED/WIG/GFF 解析 |
| `bsgenome.mbt` | Bioconductor BSgenome | 基因组序列数据库 |
| `genome_info_db.mbt` | Bioconductor GenomeInfoDb | 基因组信息管理 |
| `interaction_set.mbt` | Bioconductor InteractionSet | Hi-C 染色质交互 |

#### 差异表达分析

| MoonBit 文件 | 对应 Python 库 | 核心功能 |
| :--- | :--- | :--- |
| `deseq2.mbt` | Bioconductor DESeq2 | size factors归一化、分散度估计、负二项GLM拟合、Wald检验、LFC收缩 |
| `edger.mbt` | Bioconductor edgeR | DGEList 差异表达 |
| `limma.mbt` | Bioconductor limma | 线性模型与 voom 变换 |
| `summarized_experiment.mbt` | Bioconductor SummarizedExperiment | 多维数据容器 |
| `ballgown.mbt` | Bioconductor ballgown | 转录组水平差异表达 |
| `ruvseq.mbt` | Bioconductor RUVSeq | RNA-seq 批次效应去除 |
| `sva.mbt` | Bioconductor sva | 替代变量分析与 ComBat |
| `single_cell.mbt` | Bioconductor SingleCellExperiment | 单细胞数据分析 |
| `csaw.mbt` | Bioconductor csaw | ChIP-seq窗口差异分析、TMM归一化、负二项GLM检验 |
| `slingshot.mbt` | Bioconductor slingshot | 单细胞轨迹推断、MST构建、主曲线拟合、拟时间计算 |
| `scnorm.mbt` | Bioconductor SCnorm | 单细胞RNA-seq归一化、分位数回归、深度依赖偏差校正 |
| `edaseq.mbt` | Bioconductor EDASeq | RNA-seq探索性分析、GC含量归一化、基因长度校正 |

#### 序列组装与数据结构

| MoonBit 文件 | 对应 Python 库 | 核心功能 |
| :--- | :--- | :--- |
| `bloom_filter.mbt` | Jellyfish / khmer | Bloom Filter & k-mer 计数 |
| `bwt_fm.mbt` | Bowtie2 / BWA | BWT + FM-index |
| `de_bruijn.mbt` | SPAdes / Velvet | De Bruijn Graph 组装 |
| `suffix_array_tree.mbt` | libdivsufsort | Suffix Array & Suffix Tree |
| `olc.mbt` | Celera Assembler | Overlap-Layout-Consensus |
| `hmm.mbt` | HMMER / Augustus | Hidden Markov Model |
| `kmeans.mbt` | scikit-learn | K-means 聚类 |

#### 数据库与外部工具

| MoonBit 文件 | 对应 Python 库 | 核心功能 |
| :--- | :--- | :--- |
| `blast.mbt` | BioPython `Bio.Blast` | BLAST 结果解析 |
| `search_io.mbt` | BioPython `Bio.SearchIO` | 统一搜索结果模型 |
| `kegg.mbt` | BioPython `Bio.KEGG` | KEGG 数据库解析 |
| `medline.mbt` | BioPython `Bio.Medline` | Medline/PubMed 解析 |
| `entrez.mbt` | BioPython `Bio.Entrez` | NCBI 数据库访问 |
| `swissprot.mbt` | BioPython `Bio.SwissProt` | UniProt 记录解析 |
| `biomart.mbt` | Bioconductor biomaRt | 基因 ID 映射与注释 |
| `nexus.mbt` | BioPython `Bio.Nexus` | NEXUS 格式解析 |
| `emboss.mbt` | EMBOSS suite | EMBOSS 工具接口 |

#### 扩展功能模块

| MoonBit 文件 | 对应 Python 库 | 核心功能 |
| :--- | :--- | :--- |
| `biostrings.mbt` | Bioconductor Biostrings | IUPAC、k-mer、复杂度 |
| `protparam.mbt` | BioPython `Bio.SeqUtils.ProtParam` | 蛋白质参数分析 |
| `motifs.mbt` | BioPython `Bio.motifs` | PWM 模体识别、信息含量、序列Logo、模体富集、Pearson相关性 |
| `popgen.mbt` | BioPython `Bio.PopGen` | 群体遗传学 |
| `cluster.mbt` | BioPython `Bio.Cluster` | 序列聚类分析 |
| `restriction.mbt` | BioPython `Bio.Restriction` | 限制性内切酶分析 |
| `codon_usage.mbt` | BioPython `Bio.CodonUsage` | 密码子使用分析 |
| `affy.mbt` | BioPython `Bio.Affy` | Affymetrix 芯片分析 |
| `go_enrichment.mbt` | Bioconductor GOstats | GO 功能富集分析 |
| `fgsea.mbt` | Bioconductor fgsea | 基因集富集分析 |
| `faidx.mbt` | pyfaidx | FASTA 快速索引 |
| `feature_extraction.mbt` | 自定义 | 机器学习特征提取 |
| `dplyr.mbt` | R dplyr | DataFrame 数据操作 |
| `multi_assay_experiment.mbt` | Bioconductor MultiAssayExperiment | 多组学数据协调 |
| `utils.mbt` | 自定义 | 通用工具函数 |
| `chipseeker.mbt` | Bioconductor ChIPseeker | ChIP-seq峰注释分析(启动子/外显子/内含子/UTR分类、BED读取、peak2gene关联、多峰值集重叠分析、Venn图、饼图可视化) |
| `variation.mbt` | BioPython `Bio.Variation` | SNP分析、突变检测、氨基酸替换分析、BLOSUM62/Grantham矩阵 |
| `dssp.mbt` | BioPython `Bio.PDB.DSSP` | 二级结构预测、溶剂可及表面积、Ramachandran图 |
| `polypeptide.mbt` | BioPython `Bio.PDB.Polypeptide` | 氨基酸组成、疏水性分析、跨膜区域预测、等电点计算 |
| `dose.mbt` | Bioconductor DOSE | 疾病本体富集分析 |
| `reactome_pa.mbt` | Bioconductor ReactomePA | Reactome通路分析 |
| `annotation_dbi.mbt` | Bioconductor AnnotationDbi | 通用注释数据库接口 |
| `cluster_profiler.mbt` | Bioconductor clusterProfiler | 功能富集分析统一框架 |
| `wgcna.mbt` | Bioconductor WGCNA | 加权基因共表达网络分析 |
| `biobase.mbt` | Bioconductor Biobase | ExpressionSet数据结构、AnnotatedDataFrame |
| `geoquery.mbt` | Bioconductor GEOquery | GEO数据库数据获取、Series Matrix解析 |
| `tximport.mbt` | Bioconductor tximport | 转录本量化数据导入、基因级别汇总 |
| `single_cell_experiment.mbt` | Bioconductor SingleCellExperiment | 单细胞核心容器 (多assay、PCA/tSNE/UMAP降维、size factors) |
| `complex_heatmap.mbt` | Bioconductor ComplexHeatmap | 复杂热图可视化 (行/列聚类、颜色映射、热图注释) |
| `gsva.mbt` | Bioconductor GSVA | 基因集变异分析 (ssGSEA/zscore/PLAGE评分) |
| `chromvar.mbt` | Bioconductor chromVAR | 染色质变异分析 (TF motif富集、GC偏差校正) |
| `delayed_array.mbt` | Bioconductor DelayedArray | 延迟计算数组 (懒加载操作、分块处理、行/列聚合) |
| `annotation_filter.mbt` | Bioconductor AnnotationFilter | 基因注释过滤 (染色体筛选、生物类型过滤、区域重叠检测) |
| `sc_dbl_finder.mbt` | Bioconductor scDblFinder | 单细胞双细胞检测 (Doublet评分、最近邻搜索、PCA降维) |
| `seurat.mbt` | Bioconductor Seurat | 单细胞数据分析核心 (标准化、高可变基因、PCA、聚类、UMAP、差异表达、跨样本整合) |
| `chipseeker.mbt` | Bioconductor ChIPseeker | ChIP-seq峰值注释 (基因组区域分类(启动子/外显子/内含子/UTR/基因间区)、距离TSS分布、BED格式读取、peak2gene关联分析、注释可视化) |
| `topgo.mbt` | Bioconductor topGO | 拓扑GO富集分析 (TopGOTerm/TopGOGraph/TopGOEnrichmentResult数据结构、elim算法、weight01算法、Fisher精确检验、GO图构建) |
| `dexseq.mbt` | Bioconductor DEXSeq | 差异外显子使用分析 (ExonCount/DEXSeqDataSet/DEXSeqResult数据结构、计数归一化、统计检验、结果过滤) |
| `metagenomeseq.mbt` | Bioconductor metagenomeSeq | 零膨胀模型微生物组差异丰度分析 (MRexperiment/MGResult数据结构、中位数归一化、零膨胀概率计算、差异检验) |
| `hilbertcurve.mbt` | Bioconductor HilbertCurve | Hilbert曲线坐标映射 (HilbertCurve数据结构、编码/解码、距离计算、基因组线性化、网格映射) |
| `taxonomy.mbt` | Biopython `Bio.Taxonomy` | 分类学分析 (Taxon/TaxonomyDatabase数据结构、分类数据库解析、谱系查询、共同祖先计算、分类树操作) |
| `gff.mbt` | Biopython `Bio.GFF` | GFF3格式解析 (GFFFeature/GFFRecord数据结构、属性解析、特征提取、基因/转录本/CDS/外显子结构分析) |
| `phylo_consensus.mbt` | Biopython `Bio.Phylo.Consensus` | 一致性树构建 (ConsensusNode/Split/ConsensusTree数据结构、多数规则/严格一致性树、Newick解析、支持度计算) |
| `alphabet.mbt` | Biopython `Bio.Alphabet` | IUPAC字母表定义 (DNA/RNA/蛋白质字母表、简化字母表、空位字母表、字符验证) |
| `statistics.mbt` | scipy/stats | 统计分析 (描述统计、假设检验、相关性分析、置信区间、Z-score) |
| `freq_analysis.mbt` | Biopython `Bio.SeqUtils` | 序列频率分析 (k-mer计数、密码子使用频率、GC含量、序列复杂度、模体查找) |
| `align_analysis.mbt` | Biopython `Bio.Align.analysis` | 进化分析 (dn/ds计算、Jukes-Cantor距离、Kimura 2-parameter距离、选择压力分析) |
| `dss.mbt` | Bioconductor DSS | 离散度收缩估计与差异分析 (DSSCountData/DSSDispResult/DSSDEResult/DSSDMResult/DSSDMRResult数据结构、Wald检验、BH-FDR校正、DML/DMR检测) |
| `bamsignals.mbt` | Bioconductor bamsignals | ChIP-seq信号提取 (BamsigParams/BamsigRegion/BamsigRecord/BamsigSignal数据结构、信号计数、RPM/RPKM归一化、染色质状态分析) |
| `nucle_r.mbt` | Bioconductor nucleR | 核小体定位分析 (NucPosition/NucCallResult/NucCallParams/NucSignalTrack数据结构、信号平滑、峰值检测、位置比较、动态变化分析) |
| `phenotype.mbt` | Biopython `Bio.phenotype` | 表型微阵列分析 (WellRecord/PlateRecord/PhenFitParams/PhenControlSubtracted数据结构、logistic/Gompertz生长曲线拟合、网格搜索优化、CSV/JSON解析、控制减法) |
| `blast_applications.mbt` | Biopython `Bio.Blast.Applications` | BLAST命令行工具包装 (BlastParamSpec/BlastCommandline数据结构、8种BLAST变体(blastp/blastn/blastx/tblastn/tblastx/psiblast/rpsblast/makeblastdb)、快速构建器、命令构建、验证) |
| `psea.mbt` | Biopython `Bio.PDB.PSEA` | PSEA二级结构预测 (PseaAtom/PseaResult/PseaGeometry数据结构、CA-CA距离、虚拟键角/二面角计算、H/E/C三态分配、最小片段长度强制、三态到八态转换) |
| `sff_io.mbt` | Biopython `Bio.SeqIO.SffIO` | SFF二进制格式解析 (SffHeader/SffRead/SffFile数据结构、大端字节序u16/u32/u64读写、二进制编码/解码、质量修剪、均值质量、按名称查找) |
| `csaw.mbt` | Bioconductor csaw | ChIP-seq窗口差异分析 (CswWindow/CswDataSet/CswNormResult/CswResult数据结构、滑动窗口计数、TMM归一化、窗口过滤、负二项GLM检验、BH-FDR校正、差异区域检测) |
| `slingshot.mbt` | Bioconductor slingshot | 单细胞轨迹推断 (SlingshotNode/SlingshotEdge/SlingshotCurve/SlingshotResult数据结构、MST构建、主曲线拟合、拟时间计算、分支检测) |
| `scnorm.mbt` | Bioconductor SCnorm | 单细胞RNA-seq归一化 (SCnormQuantFit/SCnormGeneNormResult/SCnormResult数据结构、分位数回归、深度依赖偏差校正、基因特异性归一化) |
| `edaseq.mbt` | Bioconductor EDASeq | RNA-seq探索性分析 (EDASeqGeneAnno/EDASeqDataSet/EDASeqWithinResult/EDASeqBetweenResult数据结构、GC含量归一化、基因长度Loess校正、样本间归一化、RPKM计算) |
| `maftools.mbt` | Bioconductor maftools | 癌症基因组学MAF分析 (MAFMutation/MAFData/MutationSpectrum/TMBResult数据结构、SNV/Indel分类、TMB计算、突变谱分析、共现分析、Oncoplot数据生成、MAF文件解析) |
| `cnvkit.mbt` | Bioconductor CNVkit | 拷贝数变异检测 (CNVProbe/CNVSegment/CNVDataset/CBSResult数据结构、CBS循环二元分割算法、拷贝数状态判定、log2比率平滑、断点检测) |
| `destiny.mbt` | Bioconductor destiny | 单细胞扩散映射降维 (CellData/DistanceMatrix/KernelMatrix/DiffusionResult数据结构、欧氏距离矩阵、高斯核构建、Markov矩阵归一化、特征分解、扩散分量计算) |
| `tradeseq.mbt` | Bioconductor tradeSeq | 轨迹差异表达分析 (TrajectoryPoint/GeneExpressionData/GAMFit/DifferentialExpressionResult数据结构、GAM广义可加模型拟合、样条基函数、条件效应检验、BH-FDR校正) |

## 核心功能实现

### 1. 序列处理 (Bio.Seq)

提供完整的序列对象支持，包括 DNA、RNA 和蛋白质序列的创建与操作。支持互补、反向互补、转录、反转录和翻译等核心生物信息学操作。翻译功能支持多种模式：普通翻译、翻译到终止密码子、完整 CDS 翻译等，确保满足不同的序列分析需求。

### 2. 序列 I/O (Bio.SeqIO)

实现统一的序列文件解析接口，支持 FASTA、FASTQ 和 GenBank 三种常用格式的解析与写入。通过统一的 API 设计，用户可以轻松切换不同的文件格式，无需关注底层实现细节，极大简化了序列数据的处理流程。

### 3. 比对算法 (scikit-bio)

支持 DNA、RNA 和蛋白质序列的全局与局部比对。实现了 Needleman-Wunsch 全局比对和 Smith-Waterman 局部比对算法，支持自定义打分参数和替换矩阵（如 BLOSUM62）。返回比对结果包括多序列比对对象、比对分数和起始位置信息。

### 4. SAM 文件解析 (pysam)

支持 SAM 格式比对文件的解析，提供丰富的记录访问接口。可以获取读取名、标志位、参考序列名、比对位置、比对质量、CIGAR 数组、序列和质量值等信息。提供便捷方法判断比对状态，如是否配对、是否正确配对、是否未比对、是否反向互补等。

### 5. VCF 文件解析 (pysam)

实现 VCF 变异文件的解析功能，支持变异类型检测（SNP、插入、缺失）、变异定位和基因型查询。可以获取染色体、位置、REF/ALT 等位基因、QUAL 质量值、FILTER 过滤状态和 INFO 字段等信息，满足变异检测和分析需求。

### 6. 系统发育树 (Bio.Phylo)

支持 Newick 格式树的解析与创建，提供丰富的树操作方法。可以计算终端节点数量、节点间距离、寻找共同祖先，并支持 ASCII 树图可视化。树结构基于 Clade 对象构建，支持嵌套子节点和分支长度。

### 7. PDB 结构解析 (Bio.PDB)

实现蛋白质结构文件的解析，支持 Model-Chain-Residue-Atom 的四级结构层次遍历。可以获取原子坐标、计算原子间距离，支持结构叠合和 RMSD 计算。提供丰富的结构操作接口，满足蛋白质结构分析需求。

### 8. FASTA 快速索引访问 (pyfaidx)

提供 FASTA 文件的快速随机访问能力，支持 .fai 索引的构建与使用。可以获取完整序列、快速提取子序列（0-based，左闭右开区间）、获取序列长度和检查序列是否存在。支持从内容直接创建索引或从 .fai 文件加载索引。

### 9. 机器学习特征提取

提供全面的序列特征提取功能，包括 DNA 和蛋白质特征。DNA 特征包括 k-mer 频率、规范 k-mer（考虑反向互补）、核苷酸组成和 42 维特征向量。蛋白质特征包括氨基酸组成、二肽/三肽组成、理化性质（疏水性、电荷、极性、分子量）、二级结构倾向（Chou-Fasman）和 73 维特征向量。

### 10. Biostrings 序列分析 (Bioconductor Biostrings)

实现 IUPAC 核苷酸和氨基酸频率计算、k-mer 频率统计（单/双/三核苷酸）、相对同义密码子使用度（RSCU）、GC 含量（滑动窗口）、Shannon 熵、DUST 复杂度和 IUPAC 模式匹配。支持最近邻法熔解温度计算、序列比对编辑距离、序列一致性百分比和 IUPAC 反向互补操作。

### 11. GenomicRanges 基因组区间操作 (Bioconductor GenomicRanges)

提供 GRanges 数据结构，支持染色体、区间和链信息的存储与操作。支持区间偏移、缩小、调整宽度和侧翼区域获取。集合操作包括合并重叠区间、分割区间、并集、交集和差集。支持重叠检测和距离计算，包括计数重叠数、查找重叠对、计算区间距离和最近邻索引。

### 12. DESeq2 差异表达分析 (Bioconductor DESeq2)

实现完整的 RNA-seq 差异表达分析流程，支持从原始计数到差异表达基因筛选的全流程分析。可以创建 DESeqDataSet 对象管理计数矩阵、样本信息和设计矩阵。支持 size factors 估计（中位数比率法）进行测序深度校正，以及计数矩阵归一化和 log2 CPM 计算。支持分散度估计（parametric fit），结合经验贝叶斯收缩方法。支持负二项 GLM 拟合，通过迭代加权最小二乘法估计回归系数。支持 Wald 检验进行差异表达显著性检验，计算 log2 fold change、标准误、检验统计量和 p 值。支持 Benjamini-Hochberg 多重检验校正。支持 LFC 收缩（apeglm-like 方法），减小低表达基因的 fold change 估计偏差。支持显著基因筛选（按 adjusted p-value 和 LFC 阈值）和获取 top 差异表达基因。适用于 RNA-seq 差异表达分析。

### 13. Suffix Array & Suffix Tree (libdivsufsort)

实现后缀数组和后缀树数据结构，采用前缀倍增算法构建。支持模式匹配（包含、计数、定位）、最长重复子串查找和 LCP 数组计算。后缀树支持高效的字符串搜索和模式定位，适用于序列比对和重复序列分析。

### 14. Overlap-Layout-Consensus 序列组装 (Celera Assembler)

实现基于重叠-布局-一致的序列组装算法。支持重叠检测、构建重叠图、哈密顿路径搜索和一致性序列生成。提供 Graphviz 可视化输出能力，便于分析组装过程和结果。支持 DNA 序列的组装和重叠关系计算。

### 15. Hidden Markov Model 基因预测 (HMMER / Augustus)

实现隐马尔可夫模型，支持前向算法、后向算法和维特比算法。提供基因预测功能，可以从 DNA 序列中预测基因结构并提取外显子。支持 Baum-Welch 参数训练，可用于模型优化和定制。

### 16. TxDb 转录本数据库 (Bioconductor GenomicFeatures)

支持 GTF 文件解析，构建转录本数据库。可以获取基因、转录本、外显子和 CDS 的 GRanges 对象，支持启动子区域提取（上游/下游长度可配置）。支持按转录本分组获取内含子、5' UTR 和 3' UTR，以及基因/转录本/外显子的 ID 列表获取。

### 17. ProtParam 蛋白质参数分析 (Bio.SeqUtils.ProtParam)

提供全面的蛋白质参数分析功能，包括序列长度、分子量（考虑脱水）、氨基酸计数和组成百分比。理化性质分析包括 GRAVY 疏水性评分、芳香性和不稳定指数。支持等电点计算、特定 pH 下的电荷计算、二级结构倾向预测（Chou-Fasman）和信号肽预测。

### 18. rtracklayer 基因组轨道格式 (Bioconductor rtracklayer)

支持多种基因组轨道格式的解析与写入，包括 BED（3-12列）、WIG（variableStep 和 fixedStep）、BEDGraph 和 GFF/GTF。支持格式间转换，如 BED 和 GFF 转 GRanges。提供统一的解析和写入接口，便于基因组数据的处理和交换。

### 19. K-means 聚类分析 (scikit-learn)

实现 K-means 聚类算法，支持 K-means++ 初始化。提供模型训练、新数据预测、质心和标签获取功能。支持惯性计算、轮廓系数评估和最优 k 值搜索。适用于基因表达数据的聚类分析和数据分组。

### 20. SearchIO 统一搜索结果 (Bio.SearchIO)

提供统一的搜索结果模型，支持 HMMER3 tabular 格式和 BLAT PSL 格式的解析。可以获取查询 ID、命中数、top hits（按 E-value 排序）和 HSP 数量统计。支持 BLAST 结果转换为 SearchIO 模型，便于不同搜索工具结果的统一处理。

### 21. BLAST 结果解析 (Bio.Blast)

支持 BLAST tabular 和 XML 格式的解析，提供丰富的结果过滤和访问接口。可以按 E-value 和 identity 过滤 hits，获取最佳匹配和最佳 HSP。支持所有 HSPs 的获取和查询序列长度的访问。

### 22. 替换矩阵 (Bio.SubsMat)

支持多种内置替换矩阵（BLOSUM62、BLOSUM45、PAM250、PAM30），提供分数查询和不区分大小写查询功能。可以获取矩阵名称、尺寸和所有氨基酸列表，支持自定义矩阵的解析和使用。

### 23. 序列模体识别 (Bio.motifs)

实现位置权重矩阵（PWM）的创建和操作，支持序列得分计算、共识序列获取和最可能序列预测。支持 MEME 格式的解析和模体搜索功能，可以在 DNA 序列中搜索特定模体的匹配位置。支持信息含量（IC）计算、序列 Logo 数据生成、模体富集分析和 Pearson 相关性比较。

### 24. 限制性内切酶分析 (Bio.Restriction)

支持限制性内切酶的创建和酶切位点分析。可以查找酶切位点位置、酶切序列生成片段、获取酶信息（名称、识别位点、切割位置）。提供常用酶列表，便于快速访问常用内切酶。

### 25. 序列聚类分析 (Bio.Cluster)

实现距离矩阵计算和层次聚类算法，支持多种距离度量和聚类方法。可以将聚类结果转换为 Newick 格式，便于系统发育树工具的后续分析。支持轮廓系数评估和基因表达数据的聚类分析。

### 26. 群体遗传学 (Bio.PopGen)

提供等位基因频率计算、哈迪-温伯格检验和 FST 计算功能。支持核苷酸多样性（π）和 Tajima's D 的计算，适用于群体遗传学研究和进化分析。

### 27. 密码子使用分析 (Bio.CodonUsage)

实现多种密码子使用指标的计算，包括 CAI（密码子适应指数）、ENC（有效密码子数）、RSCU（相对同义密码子使用度）、GC3 含量、CBI（密码子偏好指数）和 Fop（最优密码子频率）。支持最优密码子检测和密码子翻译功能。

### 28. IRanges 整数区间操作 (Bioconductor IRanges)

提供整数区间的操作和集合运算，支持区间偏移、调整宽度和缩小范围。集合操作包括合并重叠区间、分割区间、并集、交集和差集。支持重叠检测和距离计算，包括计数重叠数、查找重叠对、计算区间距离和最近邻索引。

### 29. AlignIO 比对格式解析 (Bio.AlignIO)

支持 ClustalW、FASTA 和 Stockholm 三种比对格式的解析与写入。提供统一的比对对象接口，便于不同格式比对数据的处理和转换。

### 30. TreeIO 进化树格式解析 (Bio.TreeIO)

支持 Newick 和 NHX（Newick + Extended）格式的解析，提供树操作方法（终端节点计数、节点间距离、共同祖先查找）和 Newick 格式写入功能。

### 31. edgeR 差异表达分析 (Bioconductor edgeR)

实现基于 DGEList 的差异表达分析，支持归一化因子计算、GLM 拟合和差异表达检验（精确检验和似然比检验）。支持 top tags 的获取和显著差异基因的筛选。

### 32. limma 差异表达分析 (Bioconductor limma)

实现基于线性模型的差异表达分析，支持经验贝叶斯校正和 voom 变换（适用于 RNA-seq 数据）。支持设计矩阵创建、线性模型拟合、对比矩阵分析和 top 差异基因的获取。

### 33. SummarizedExperiment 多维数据容器 (Bioconductor SummarizedExperiment)

提供多维基因组数据容器，支持多个 Assays（如表达矩阵、计数矩阵）、行数据（基因信息）和列数据（样本信息）的存储与操作。支持子集操作和合并功能，便于多组学数据的协调管理。

### 34. GenomicAlignments 基因组比对分析 (Bioconductor GenomicAlignments)

提供 GAlignments 对象，支持比对信息的存储和操作。可以计算覆盖度、按特征汇总比对数、执行 pileup 操作。支持按 MAPQ 和 strand 过滤，以及转换为 GRanges 对象。

### 35. VariantAnnotation 变异注释 (Bioconductor VariantAnnotation)

支持变异类型检测（SNP、插入、缺失）、变异定位和编码效应预测（同义/错义/无义变异）。可以获取氨基酸变化信息和变异汇总统计（总变异数、SNP 数量、编码区变异数、受影响基因）。

### 36. Affy Affymetrix芯片数据分析 (Biopython Bio.Affy)

支持 Affymetrix 芯片数据的处理，包括探针集创建、AffyBatch 对象管理、背景校正、PM-MM 差异计算和探针集汇总。支持 RMA 标准化流程（背景校正 + log2 转换 + 分位数归一化 + 中位数汇总）。

### 37. SVDSuperimposer SVD蛋白质结构叠合 (Biopython Bio.PDB.SVDSuperimposer)

实现基于 SVD 的蛋白质结构叠合算法，支持原子坐标的旋转和平移变换。可以计算 RMSD 值、旋转矩阵和平移向量，支持直接叠合并获取结果，以及不进行叠合的 RMSD 直接计算。

### 38. GOEnrichment GO功能富集分析 (Bioconductor GOstats/clusterProfiler)

实现 GO 功能富集分析，支持超几何检验和 BH 校正。可以构建 GO 注释数据库、执行富集分析、过滤结果（按校正后 p 值和富集基因数）和按命名空间筛选（BP、MF、CC）。提供独立的统计检验功能，包括超几何检验、Bonferroni 校正和 BH 校正。

### 39. SingleCell 单细胞数据分析 (Bioconductor SingleCellExperiment)

提供单细胞数据分析功能，包括 QC 指标计算（每个细胞的 UMI 数、检测到的基因数、线粒体基因比例）、细胞过滤、Log 标准化、高变异基因检测和 PCA 降维（幂迭代法）。支持 SingleCellExperiment 对象的创建和管理。

### 40. BAM 文件解析 (pysam)

支持 BAM 格式比对文件的解析，提供丰富的记录访问接口。可以获取读取名、标志位、参考序列名、比对位置、比对质量、CIGAR 数组、序列和质量值等信息。提供便捷方法判断比对状态和计算插入片段长度。

### 41. Bloom Filter & k-mer 计数 (Jellyfish / khmer)

实现 Bloom Filter 概率数据结构，支持 k-mer 的快速成员查询。提供 k-mer 精确计数器，支持近似去重和唯一 k-mer 数量估算。适用于大规模序列数据的快速处理和去重。

### 42. BWT + FM-index (Bowtie2 / BWA)

实现 Burrows-Wheeler Transform 和 FM-index，支持 BWT 正逆变换和高效的模式匹配。可以判断模式是否存在、计数出现次数和定位出现位置，适用于大规模序列数据的快速搜索。

### 43. De Bruijn Graph 序列组装 (SPAdes / Velvet)

实现基于 De Bruijn Graph 的序列组装算法，支持 k-mer 节点构建、欧拉路径查找和序列组装。提供图简化功能（去除气泡），可以获取所有节点和边信息，适用于短读长序列的组装。

### 44. Suffix Array & Suffix Tree (libdivsufsort)

实现后缀数组和后缀树数据结构，支持模式匹配（包含、计数、定位）、最长重复子串查找和 LCP 数组计算。适用于 DNA 序列的重复序列分析和模式搜索。

### 45. Smith-Waterman 局部序列比对

实现 Smith-Waterman 局部序列比对算法，支持自定义打分矩阵（如 BLOSUM62）和空位罚分（空位开放和空位延伸）。返回比对结果包括比对对象、分数和起始位置，支持蛋白质和 DNA 序列的比对。

### 46. Needleman-Wunsch 全局序列比对

实现 Needleman-Wunsch 全局序列比对算法，支持自定义匹配/错配得分和空位罚分。支持蛋白质和 DNA 序列的比对，提供回溯矩阵的获取功能。

### 47. dplyr 数据操作 (R dplyr)

提供 DataFrame 数据操作功能，支持过滤、选择列、添加新列、排序、分组汇总和连接操作。采用链式调用风格，便于数据处理流程的构建和表达。

### 48. KEGG 数据库解析 (Biopython Bio.KEGG)

支持 KEGG Gene、Pathway、Compound 和 Enzyme 记录的解析，提供丰富的字段访问接口。可以获取通路中的基因 ID、计算基因参与的通路数量和检查化合物是否在特定通路中。

### 49. Medline/PubMed 文献解析 (Biopython Bio.Medline)

支持 Medline/PubMed 文献记录的解析，包括单条和多条记录的解析。提供 APA 引用格式生成、按 MeSH 术语过滤记录和按年份统计文献的功能。

### 50. BSgenome 基因组序列数据库 (Bioconductor BSgenome)

提供基因组序列数据库的管理功能，支持染色体序列的添加、检索和子序列提取。支持链特异性基因提取（正链和负链）、染色体长度查询和基因组总长度计算。

### 51. biomaRt 基因ID转换和注释查询 (Bioconductor biomaRt)

提供基因 ID 映射和注释查询功能，支持 Ensembl 基因 ID 到基因名称的映射、批量 ID 映射和基因注释查询（描述、染色体位置、生物类型、外部数据库 ID）。支持批量查询功能，便于大规模基因注释的获取。

### 52. RUVSeq RNA-seq批次效应去除 (Bioconductor RUVSeq)

实现 RNA-seq 数据的批次效应去除，支持数据标准化、log2 转换、RUVg（使用对照基因估计批次效应因子）和 RUVs（使用所有基因估计批次效应因子）方法。支持批次效应的去除和基因索引的获取。

### 53. PDB 高级结构分析 (Bio.PDB.Polypeptide / DSSP)

提供蛋白质结构的高级分析功能，包括主链二面角计算（phi、psi、omega）、四原子二面角计算、CA 原子距离矩阵和接触图生成。支持氢键检测、二级结构分配（DSSP-style）、二级结构统计、回转半径计算和 Ramachandran 图数据生成。支持 SASA（溶剂可及表面积）计算（Shrake-Rupley 算法）、结构质量评估（Ramachandran 区域分类）、Kyte-Doolittle 疏水性分析和序列属性距离矩阵生成。

### 54. 系统发育树高级分析 (Bio.Phylo.TreeMetrics)

提供系统发育树的高级分析功能，包括总分支长度计算、最大深度计算、叶节点名称获取和 Colless 平衡指数计算。支持系统发生距离计算（沿路径的分支长度之和）、距离矩阵生成、内部节点计数、二分体获取和 Robinson-Foulds 距离计算（树拓扑差异度量）。

### 55. 序列复杂度与组成分析 (Bio.SeqUtils.Complexity)

提供序列复杂度和组成的分析功能，包括 Shannon 熵计算、语言学复杂度（observed/max k-mers）、GC/AT 含量计算、GC/AT 偏斜分析和 DUST 低复杂度评分。支持序列签名（k-mer 频率向量）、核苷酸组成分析、混沌游戏表示（CGR）和序列相似度计算（余弦相似度）。

### 56. AlignInfo 比对统计 (Bio.Align.AlignInfo)

提供多序列比对的统计功能，包括一致性序列生成、多数序列生成（>=50% 共识）和严格一致性序列生成（100% 保守）。支持保守性分析、Shannon 熵谱计算、成对序列同一性计算、保守位点查找（>=80% 保守）和可变位点查找（<=50% 保守）。

### 57. CodonAlign 密码子比对与选择压力分析 (Bio.codonalign)

提供密码子比对和选择压力分析功能，包括密码子翻译、密码子替换分类（同义/非同义）、替换数计数和 dN/dS 选择压力分析（Nei-Gojobori 方法）。支持 Jukes-Cantor 多重命中校正、密码子使用偏好（RSCU）、有效密码子数（ENC）计算和同义/非同义位点计数。

### 58. Entrez NCBI 数据库访问 (Bio.Entrez)

提供 NCBI 数据库的访问功能，包括 ESearch（搜索 PubMed、Gene 等数据库）、EFetch（获取完整记录）、EGQuery（全局查询）和 EInfo（数据库信息）。支持 PubMed 文章和 Gene 记录的解析，以及基因信息和分类学信息的获取。

### 59. GenomeInfoDb 基因组信息管理 (Bioconductor GenomeInfoDb)

提供基因组信息的管理功能，支持预构建基因组（hg38、hg19、mm10 等）的加载和访问。可以获取染色体名称、长度、着丝粒位置和染色体臂信息（p 短臂 / q 长臂）。支持标准染色体筛选、基因组构建查询和自定义基因组构建。

### 60. InteractionSet 染色质交互数据 (Bioconductor InteractionSet)

提供染色质交互数据的管理功能，支持锚点对（Anchor）和交互集合（GInteraction）的创建与操作。可以构建 Hi-C 交互矩阵、计算交互距离分布和提取 Top 交互。支持交互数据的子集操作和样本级别的交互分析。

### 61. MultiAssayExperiment 多组学数据协调 (Bioconductor MultiAssayExperiment)

提供多组学实验数据的协调管理功能，支持多个实验（Experiment）的组织和样本映射（SampleMap）。可以建立跨实验的样本对应关系，执行跨实验子集操作和数据整合。适用于多组学数据（如基因组、转录组、蛋白质组）的联合分析。

### 62. TreeConstruction 系统发育树构建 (Bio.Phylo.TreeConstruction)

提供基于距离矩阵的系统发育树构建功能，支持 UPGMA、WPGMA 和 NJ（邻接法）三种算法。支持距离计算（Jukes-Cantor、Kimura 替换模型）和树的构建，适用于分子进化和系统发育分析。

### 63. NeighborSearch KD树近邻搜索 (Bio.PDB.NeighborSearch)

实现基于 KD 树的空间搜索算法，支持半径搜索、最近邻查找和原子对搜索。可以在蛋白质结构中快速查找距离在指定范围内的原子，适用于蛋白质结构分析和分子对接研究。

### 64. SwissProt 蛋白数据库解析 (Bio.SwissProt)

支持 Swiss-Prot/UniProt 蛋白质数据库记录的解析，包括特征注释、参考文献、关键词和序列信息的提取。可以获取蛋白质的完整信息，适用于蛋白质功能注释和数据库查询。

### 65. mmCIF 格式解析 (Bio.PDB.MMCIFParser)

支持 mmCIF（macromolecular Crystallographic Information File）格式的解析，包括数据块（DataBlock）、类别（Category）和原子位点的提取。可以获取蛋白质结构的完整信息，适用于结构生物学研究。

### 66. Nexus 格式解析 (Bio.Nexus)

支持 NEXUS 格式的解析，包括数据矩阵、系统发育树、距离矩阵和分类单元的提取。可以获取系统发育分析所需的完整数据，适用于进化生物学研究。

### 67. EMBOSS 工具接口 (EMBOSS suite)

提供 EMBOSS 生物信息学工具包的接口功能，包括 GC 偏斜、AT 偏斜、分子量计算、Tm 值计算、反向互补、翻译、ORF 查找、Hamming/Levenshtein 距离计算、k-mer 计数、等电点、脂肪族指数、GRAVY 评分、不稳定指数和氨基酸组成分析。

### 68. CRAM 格式解析 (pysam)

支持 CRAM 压缩二进制序列比对格式的解析，包括魔术数字检测、版本识别和压缩类型判断（Gzip、Bzip2、Lzma）。可以解析 CRAM 头部信息（参考序列、流信息、SAM 头部）和记录（比对位置、读取特征），支持 CRAM 到 BAM 的转换，便于与现有比对工具的互操作。

### 69. ChIPseeker ChIP-seq峰注释分析 (Bioconductor ChIPseeker)

实现 ChIP-seq 峰注释功能，支持将峰区域与基因特征进行关联分析。可以计算峰到转录起始位点（TSS）的距离，将峰分类为启动子区、外显子区、内含子区、UTR区（5'UTR/3'UTR）或基因间区。支持 BED 格式峰值文件的读取，提供 peak2gene 峰值到基因关联分析功能。支持多峰值集重叠分析（peak_overlap）和 Venn 图可视化，便于比较不同样本或条件下的峰值重叠情况。提供注释分布饼图（plot_anno_pie）和条形图（plot_annotation_bar）可视化，以及距离 TSS 分布图。提供注释结果汇总统计功能，便于分析峰在基因组上的分布特征。支持自定义启动子区域范围和注释分类标准。

### 70. DOSE 疾病本体富集分析 (Bioconductor DOSE)

实现基于疾病本体（Disease Ontology）的富集分析功能，支持超几何检验计算富集显著性。可以计算原始 p 值和校正后 p 值（adjusted p-value），支持按 p 值阈值过滤富集结果。提供富集分数计算和结果汇总功能，便于识别与基因列表相关的疾病术语。支持自定义背景基因集和疾病数据库。

### 71. ReactomePA Reactome通路分析 (Bioconductor ReactomePA)

实现 Reactome 通路富集分析功能，支持将基因列表与 Reactome 通路数据库进行关联。可以计算通路富集显著性，支持按顶层术语对通路进行分组（如细胞过程、DNA 复制、代谢过程等）。提供富集结果可视化和汇总功能，便于分析基因参与的生物学通路。支持自定义通路数据库和背景基因集。

### 72. AnnotationDbi 通用注释数据库接口 (Bioconductor AnnotationDbi)

提供通用的基因注释数据库接口，支持基因信息的存储和查询。可以创建注释数据库，添加基因信息（基因 ID、符号、描述、染色体位置、链、生物类型），执行 ID 映射（如 symbol 到 ensembl），选择特定列进行查询。支持 GO 和 KEGG 注释的添加与获取，便于基因功能注释的管理和查询。

### 73. clusterProfiler 功能富集分析统一框架 (Bioconductor clusterProfiler)

实现功能富集分析的统一框架，支持通用的富集分析流程。可以构建术语-基因映射数据库，执行超几何检验计算富集显著性，支持多重检验校正（adjusted p-value）。提供结果过滤、top 术语获取、富集图谱和点图可视化功能。适用于各种功能富集分析场景，如 GO、KEGG、DO 等。

### 74. WGCNA 加权基因共表达网络分析 (Bioconductor WGCNA)

实现加权基因共表达网络分析功能，支持表达矩阵的处理和模块检测。可以构建邻接矩阵（基于相关系数的幂函数），计算 TOM（Topological Overlap Measure）相似度，执行模块检测（动态剪切树）。提供软阈值选择、模块特征基因计算和结果汇总功能。适用于基因表达数据的共表达网络构建和模块分析。

### 75. Biobase ExpressionSet基础数据结构 (Bioconductor Biobase)

实现 Bioconductor 的核心数据结构 ExpressionSet，支持多维基因组数据的存储与操作。可以创建包含 assay 数据（表达矩阵）、表型数据（样本信息）、特征数据（基因信息）的完整数据对象。提供维度查询、样本/特征名称获取、子集操作、行过滤（按表达量阈值）、数据归一化（z-score）和 log2 转换功能。支持 AnnotatedDataFrame 数据结构，用于存储带注释的数据框。

### 76. GEOquery GEO数据库数据获取 (Bioconductor GEOquery)

实现 NCBI GEO 数据库数据的解析与获取功能，支持 Series Matrix 格式和 SOFT 格式的解析。可以获取系列 ID、平台 ID、样本 ID、基因名称、表达数据和样本信息。提供基因过滤（按表达量范围）、log2 转换和转换为 ExpressionSet 的功能。支持 GEOSeries、GEOSample 和 GEOPlatform 数据结构的创建与操作，便于 GEO 数据的完整管理和分析。

### 77. tximport 转录本量化数据导入 (Bioconductor tximport)

实现转录本量化数据的导入与处理功能，支持 Salmon 输出的 quant.sf 文件解析。可以解析转录本名称、长度、有效长度、TPM 和读取计数。支持多样本数据合并、转录本到基因级别的汇总（基于 tx2gene 映射）、去除转录本版本号和低表达过滤。提供转换为 ExpressionSet 的功能，支持 counts、TPM 和 length 三种数据类型的导入，便于后续差异表达分析和功能富集分析。

### 78. AnnotationHub 中心化注释资源访问 (Bioconductor AnnotationHub)

实现中心化注释资源的访问与管理功能，支持多种类型的注释资源（GeneModel、Variation、Epigenome、Sequence、Alignment、Expression）和数据提供者（Ensembl、UCSC、NCBI、EMBL、Custom）。可以按类型、提供者、基因组进行资源搜索，支持关键词查询和资源摘要生成。提供资源添加、删除和获取功能，便于注释数据的统一管理和访问。

### 79. GenomicFeatures 基因组注释功能 (Bioconductor GenomicFeatures)

实现基因组注释的核心功能，支持 Gene、Transcript、Exon 数据结构的创建与操作。可以解析 GTF 文件构建基因组注释数据库，支持基因按 ID 或名称查询、按染色体筛选基因、按区域查询基因。提供转录本信息获取、基因生物类型统计和注释摘要生成功能。支持基因/转录本/外显子数量统计和注释信息的完整管理。

### 80. graph 图数据结构 (Bioconductor graph)

实现图数据结构和图算法，支持有向图和无向图的创建与操作。可以添加节点和边，计算节点度和邻居列表，检测边的存在。支持最短路径算法（BFS）、连通分量检测、图摘要生成和 DOT 格式输出。提供示例图和通路图的创建功能，适用于生物学通路分析和网络分析。

### 81. DropletUtils 空液滴检测 (Bioconductor DropletUtils)

实现单细胞 RNA-seq 数据的空液滴检测功能，支持 barcode 排序、knee 点检测和 emptyDrops 算法。可以计算液滴统计指标（总计数、检测基因数），对液滴按总计数排序，找到 knee 点估计细胞数量。支持基于 Monte Carlo 模拟的空液滴检测，计算每个液滴为空的概率和 FDR 值，进行细胞过滤。

### 82. scran 单细胞归一化与聚类 (Bioconductor scran)

实现单细胞 RNA-seq 数据的归一化和聚类分析功能，支持 sum_factors 归一化方法、log 归一化和 SNN（Shared Nearest Neighbors）图构建。可以计算 size factors、构建 SNN 图（基于 Jaccard 相似度）、执行 Leiden 聚类算法（基于社区检测）。提供差异标志物分析功能，支持按聚类结果识别细胞类型特异性基因。

### 83. monocle3 单细胞轨迹分析 (Bioconductor monocle3)

实现单细胞轨迹分析功能，支持 CellDataSet 数据结构的创建与操作。可以执行数据预处理（PCA 降维）、维度约简（PCA 和 UMAP）、主图学习（基于最小生成树）和拟时间排序。支持差异表达分析（按拟时间拟合模型）和结果汇总。新增分支点检测功能（find_branch_points），可以识别主图中的分支节点。支持分支特异性差异表达分析（differential_gene_test_branches），比较不同分支上的基因表达差异。适用于细胞分化轨迹分析和发育过程建模。

### 84. ShortRead 短读序列质量控制 (Bioconductor ShortRead)

实现短读序列质量控制功能，支持 FASTQ 数据的全面质量评估。可以计算 QC 统计指标（总读长、GC 含量、平均质量、重复率、每周期质量分布、N 碱基分布）。支持序列修剪（尾部质量修剪、前端质量修剪、adapter 序列修剪）和读长过滤（最小质量、最小长度、最大 N 碱基）。支持生成 FastQC 格式报告（PASS/WARN/FAIL 状态判定）。适用于测序数据预处理和质量控制。

### 85. scater 单细胞质量控制 (Bioconductor scater)

实现单细胞 RNA-seq 质量控制功能，支持细胞水平和基因水平的 QC 指标计算。可以计算细胞 QC 指标（总计数、检测基因数、线粒体比例、核糖体比例、Top50 基因比例）和基因 QC 指标（平均表达、检测细胞数、检测百分比）。支持细胞过滤（基于计数、基因数、线粒体比例等）和基因过滤（最小检测细胞数）。支持标准化方法（CPM、log2-CPM）、高变基因检测（HVG）和 PCA 降维。适用于单细胞数据预处理和质量控制。

### 86. MAST 单细胞差异表达分析 (Bioconductor MAST)

实现单细胞差异表达分析功能，采用 Hurdle 模型（零膨胀模型）处理单细胞数据的零膨胀特性。模型包含两个组分：离散组分（Fisher 精确检验检测率差异）和连续组分（Welch t 检验表达水平差异）。支持使用卡叉分布合并两个 p 值得到联合检验结果，使用 Benjamini-Hochberg 方法进行 FDR 多重检验校正。支持计算 log2 倍数变化、检测率统计和结果汇总（差异基因计数、Top 基因列表）。适用于单细胞转录组差异表达分析。

### 87. GenomicFiles 分布式基因组文件处理 (Bioconductor GenomicFiles)

实现分布式处理大型基因组文件的功能，支持按染色体/区间查询 BAM、BED、VCF 等文件。可以创建 GenomicFile 对象（支持 BAM、BED、VCF、BigWig 等格式）和 GenomicFilesCollection 集合，按区间进行扫描。支持 BAM 扫描（返回读数和质量值统计）、BED 扫描（返回特征数和分数统计）、VCF 扫描（返回变异数和质量值统计）。支持按区间归约（sum/mean/max/min）、按批次产生 reads、统计各区间 reads 数、计算覆盖度和批量区间查询。支持解析 BED 行、VCF 行和 SAM 行。适用于处理超过内存容量的大型基因组文件。

### 88. DiffBind ChIP-seq差异结合分析 (Bioconductor DiffBind)

实现 ChIP-seq 差异结合分析功能，支持峰值一致性分析和差异结合位点检测。可以创建 DBAConfig 配置对象管理多个样本（DBASample）和峰（DBAPeak）。支持峰值重叠分析（计算重叠矩阵）、共识峰识别（基于最小重叠阈值）、合并重叠峰。支持多种归一化方法（TMM、RPKM、CPM、library size），其中 TMM 归一化包含 M/A 值修剪和加权计算。差异结合分析采用负二项分布模型（类似 DESeq2），包含离散度估计、NB 检验、BH-FDR 多重检验校正。支持报告显著差异峰、生成火山图数据（fold change vs -log10 pvalue）和 MA 图数据（mean expression vs log2 fold change）。适用于 ChIP-seq 实验中识别条件间差异结合位点。

### 89. minfi DNA甲基化分析 (Bioconductor minfi)

实现 Illumina 450K/EPIC DNA 甲基化数组数据分析功能，支持从原始信号到差异甲基化分析的完整流程。可以创建 RGChannelSet（Red/Green 通道信号）、MethylSet（甲基化/未甲基化信号）和 GenomicRatioSet（β值/M值）。支持多种预处理方法：NOOB（Normal-Exponential 卷积模型背景校正）、Illumina（底部 10% 均值背景校正）、分位数归一化、功能归一化（funnorm，结合 NOOB 和 M 值分位数归一化）。支持 β 值计算（M/(M+U+offset)）和 M 值计算（log2(M/U)）相互转换。支持质量控制（QC 报告、检测 p 值、性别预测）、差异甲基化探针（DMP）分析（基于 Welch t 检验）、差异甲基化区域（DMR）分析（基于邻近探针聚合）。支持探针过滤（低质量、SNP 相关）。适用于表观遗传学研究和 DNA 甲基化分析。

### 90. flowCore 流式细胞术FCS文件处理 (Bioconductor flowCore)

实现流式细胞术（Flow Cytometry）FCS 文件处理和分析功能，支持 FCS 文件解析、数据变换、荧光补偿和门控（gating）等核心功能。可以创建 FCSFile（FCS 文件对象）、FlowFrame（单样本数据容器）、FlowSet（多样本集合）。支持多种数据变换方法：对数变换（log）、asinh 变换（反双曲正弦）、logicle 变换（结合对数和线性变换的优点，适用于流式数据）。支持荧光补偿（CompensationMatrix 矩阵乘法）消除荧光通道间溢出。支持多种门控类型：矩形门（RectangleGate，基于参数范围）、多边形门（PolygonGate，点在多边形内判断）、椭球门（EllipsoidGate，点在椭圆内判断）、四象限门（QuadrantGate，四象限分类）。支持归一化（wed-based、range-based）和多个 FlowFrame 合并。适用于临床流式细胞术数据分析和免疫表型研究。

### 91. bsseq 亚硫酸氢盐测序分析 (Bioconductor bsseq)

实现亚硫酸氢盐测序（Bisulfite Sequencing, BS-seq）数据分析功能，支持从甲基化调用到差异甲基化区域（DMR）检测的完整流程。可以创建 BSRecord（单个 CpG 位点记录）和 BSData（多位点数据集合）。支持解析 BED 格式和 Bismark 输出格式的甲基化数据。支持甲基化率计算（meth/(meth+unmeth)）、覆盖度统计、按覆盖度和上下文（CpG/CHG/CHH）过滤。支持 BSmooth 平滑算法（滑动窗口加权平均）处理相邻 CpG 位点的甲基化信号。支持差异甲基化区域检测：基于 t 统计量的 DMR 检测、Welch t 检验、BH-FDR 多重检验校正。支持高甲基化（hyper）和低甲基化（hypo）区域识别、CpG 合并、区域平均甲基化计算。适用于全基因组亚硫酸氢盐测序（WGBS）和简化代表性亚硫酸氢盐测序（RRBS）数据分析。

### 92. SingleCellExperiment 单细胞实验核心数据结构 (Bioconductor SingleCellExperiment)

实现单细胞实验的核心数据结构 SingleCellExperiment (SCE)，作为 scater、scran、monocle3 等单细胞分析包的基础容器。支持多种 assay（表达矩阵）管理（counts、logcounts、normcounts 等），可同时存储原始计数、归一化对数表达等多种数据形式。支持基因/细胞元数据（row_data、col_data）存储。支持降维结果（reduced_dims）管理，可存储 PCA、t-SNE、UMAP 等多种降维结果。支持 Builder 模式流式构建 SCE。支持 size factors 计算（library size 法和 deconvolution 法）和归一化（log2(counts/size_factors + 1)）。支持基因过滤（按最小检测细胞数）和细胞过滤（按最小检测基因数）。支持基因和细胞子集选择、多个 SCE 合并。是单细胞转录组学、单细胞 ATAC-seq、单细胞蛋白质组学等单细胞分析的基础设施。

### 93. ComplexHeatmap 复杂热图可视化 (Bioconductor ComplexHeatmap)

实现复杂热图可视化功能，支持热图注释、行/列聚类、颜色映射和多个热图组合。可以创建 HeatmapData 对象（矩阵、行名、列名、聚类树）和 HeatmapAnnotation 对象（行注释、列注释、单元格注释）。支持颜色映射（ColorMap）和热图选项（HeatmapOptions）配置。支持矩阵标准化（按行/按列）、层次聚类（complete/single/average/ward.D）和多种距离度量（欧氏距离、皮尔逊距离、余弦距离）。支持绘制 ASCII 热图、创建图例、合并多个热图。支持按组拆分行/列和添加热图注释。适用于基因表达数据可视化、差异表达结果展示和多组学数据整合分析。

### 94. GSVA 基因集变异分析 (Bioconductor GSVA)

实现基因集变异分析（Gene Set Variation Analysis）功能，支持单样本通路评分。可以创建 GeneSet 对象（基因集名称、基因列表、描述）和 GSVAData 对象（表达矩阵、基因名、样本名）。支持多种评分方法：ssGSEA（单样本 GSEA）、zscore（Z-score 标准化）、PLAGE（主成分分析）和原始 GSVA 方法。支持基因集过滤（按最小/最大基因数）、基因排序和评分归一化。支持置换检验、差异分析（t 检验）和富集分析（NES/ES 计算、Leading Edge 基因、BH-FDR 校正）。支持统计摘要（均值、标准差、中位数、最小值、最大值）和结果可视化（ASCII 评分图、ASCII 热图）。适用于 RNA-seq 和微阵列数据的基因集富集分析和通路活性评估。

### 95. ChromVAR 染色质变异分析 (Bioconductor chromVAR)

实现单细胞 ATAC-seq 数据分析功能，支持染色质变异分析、TF motif 富集和偏差校正。可以创建 PeakSet 对象（峰值集合、染色体位置）、Motif 对象（TF 模体、PWM 矩阵、共有序列）和 ChromVARData 对象（计数矩阵、峰值集、细胞名称、Motifs、GC 含量）。支持偏差计算（DeviationResult）、GC 偏差校正、Motif 匹配和评分。支持变异性分析（VariabilityResult）、Z-score 转换、差异偏差分析（t 检验）和 Motif 富集分析（MotifEnrichment）。支持细胞聚类（k-means）、高变 Motif 识别、Peak 过滤和归一化（CPM）。支持结果摘要和 ASCII 偏差图可视化。适用于单细胞 ATAC-seq 数据分析和转录因子活性评估。

### 96. DelayedArray 延迟计算数组 (Bioconductor DelayedArray)

实现延迟计算数组数据结构，支持懒加载操作和分块处理。可以创建 DelayedArray 对象（数据矩阵、维度、延迟操作队列、块大小），支持行/列求和（da_row_sums、da_col_sums）、行/列均值（da_row_means、da_col_means）、转置（da_transpose）、行/列子集（da_subset_rows、da_subset_cols）和数据汇总（da_summary）。延迟操作队列允许累积操作而不立即执行，适用于大规模数据的高效处理。支持示例数据生成和数组计算。

### 97. AnnotationFilter 基因注释过滤 (Bioconductor AnnotationFilter)

实现基因注释数据的过滤和查询功能，支持多种过滤条件和统计分析。可以创建 Annotation 对象（基因 ID、符号、染色体、起始/终止位置、链、生物类型、描述）。支持按染色体过滤（af_filter_chromosome）、生物类型过滤（af_filter_biotype）、链过滤（af_filter_strand）、描述包含过滤（af_filter_description_contains）、符号模式匹配（af_filter_symbol_pattern）、生物类型集合过滤（af_filter_biotype_in）和区域重叠检测（af_filter_region）。支持按生物类型和染色体统计（af_count_by_biotype、af_count_by_chromosome）和过滤结果汇总（af_filter_summary）。适用于基因注释数据的筛选和分析。

### 98. scDblFinder 单细胞双细胞检测 (Bioconductor scDblFinder)

实现单细胞 RNA-seq 数据的双细胞（doublet）检测功能，支持 Doublet 评分计算和细胞过滤。可以创建 SingleCellData 对象（计数矩阵、细胞名称、基因名称）和 DoubletScore 对象（细胞名称、评分、是否为双细胞）。支持距离计算（scdf_compute_distance）、最近邻搜索（scdf_find_nearest_neighbors）、Doublet 评分计算（scdf_compute_doublet_score）、双细胞检测（scdf_detect_doublets）、结果汇总（scdf_doublet_summary）和细胞过滤（scdf_filter_doublets）。支持 PCA 降维（scdf_compute_pca）用于降维后距离计算。适用于单细胞数据的质量控制和双细胞去除。

### 99. ChIPseeker ChIP-seq峰值注释 (Bioconductor ChIPseeker)

实现 ChIP-seq 峰值注释和基因组区域分析功能，支持峰值到基因的关联和可视化。可以创建 Peak 对象（染色体、起始/终止位置、峰值 ID）和 PeakAnnotation 对象（注释信息、基因名、到 TSS 的距离、功能类型、链、转录本 ID、外显子 ID、UTR 类型）。支持基于 TxDb 的峰值注释（annotate_peaks），自动查找最近基因并计算到 TSS 的距离。基因组区域分类包括：启动子（<=0bp、0-2kb、2-5kb）、外显子区、内含子区、UTR区（5'UTR/3'UTR）、远端区域（5-10kb、10-50kb）和基因间区域（>50kb）。支持 BED 格式峰值文件读取（read_peak_file）和 peak2gene 峰值到基因关联分析。支持注释分布条形图（plot_annotation_bar）和距离 TSS 分布图（plot_dist_to_tss）的 ASCII 可视化。支持统计分析，包括各区域百分比计算和按功能类型过滤（get_promoter_annotations、get_exon_annotations、get_intron_annotations、get_utr_annotations、get_distal_annotations、get_intergenic_annotations）。支持注释结果汇总（annotation_summary）。适用于 ChIP-seq 数据分析和峰值功能注释。

### 100. Bio.Alphabet IUPAC字母表定义 (Biopython Bio.Alphabet)

实现 IUPAC 字母表的定义和验证功能，支持 DNA、RNA 和蛋白质的各种字母表类型。支持 IUPAC 明确 DNA 字母表（A、C、G、T）、IUPAC 模糊 DNA 字母表（含 R、Y、S、W、K、M、B、D、H、V、N 等简并碱基）、IUPAC 明确 RNA 字母表（A、C、G、U）、IUPAC 蛋白质字母表（20 种标准氨基酸 + B、Z、X 等）、简化 DNA/RNA 字母表（仅含 A、C、G、T/U）以及带空位的字母表（包含 '-' 空位字符）。提供字母验证功能（检查字符是否属于字母表）、字母表属性访问（名称、字母列表、是否带空位）。支持三字母到一字母的氨基酸代码映射（Ala→A、Gly→G 等）。适用于序列数据的有效性验证和标准化。

### 101. Bio.Statistics 统计分析 (scipy/stats)

实现统计分析功能，包括描述统计、假设检验和相关性分析。描述统计支持均值、方差、标准差、中位数、最小值、最大值和求和计算。相关性分析支持 Pearson 相关系数（参数方法，衡量线性相关性）和 Spearman 秩相关系数（非参数方法，衡量单调相关性）。假设检验支持单样本 t 检验（检验样本均值是否等于假设值）和 t 统计量计算。支持置信区间计算（95% 置信水平）和 Z-score 标准化。适用于基因表达数据分析、生物标志物筛选和统计推断。

### 102. Bio.FreqAnalysis 序列频率分析 (Biopython Bio.SeqUtils)

实现序列频率分析功能，包括核苷酸频率、k-mer 计数、密码子使用频率和序列复杂度分析。核苷酸频率支持计算 A、T、G、C 的频率分布。GC 含量和 AT 含量计算适用于序列组成分析。k-mer 分析支持任意 k 值的频率统计，包括二核苷酸和三核苷酸频率。密码子使用频率支持编码序列的密码子计数和频率计算，适用于密码子偏好分析。序列复杂度支持 Shannon 熵计算，衡量序列的信息含量。模体查找支持在序列中搜索特定模式的出现位置。适用于基因组学研究和序列特征分析。

### 103. Bio.Align.analysis 进化分析 (Biopython Bio.Align.analysis)

实现进化分析功能，包括 dn/ds 计算、进化距离和选择压力分析。dn/ds 计算采用 Nei-Gojobori 方法，分别估计同义替换率（ds）和非同义替换率（dn），并计算 dn/ds 比值判断选择压力（>1 为正选择，<1 为净化选择）。支持 Jukes-Cantor 单参数距离校正和 Kimura 双参数距离校正（考虑转换/颠换差异）。支持 p-distance（原始差异比例）计算和距离矩阵生成。支持氨基酸翻译和密码子替换分类。适用于分子进化研究、选择压力分析和系统发育分析。

### 104. csaw ChIP-seq窗口差异分析 (Bioconductor csaw)

实现 ChIP-seq 基于窗口的差异结合分析功能，支持从滑动窗口计数到差异区域检测的完整流程。可以创建 CswWindow（基因组窗口）和 CswDataSet（窗口数据集）进行多样本管理。支持滑动窗口计数（基于 read 位置分配到窗口）和 TMM 归一化（trimmed mean of M-values）进行文库大小校正。支持窗口过滤（基于 abundance 阈值、log2 CPB 转换）和差异结合检验（负二项 GLM 拟合，Wald 检验）。支持 Benjamini-Hochberg 多重检验校正和差异区域检测（合并相邻显著窗口）。提供多种归一化方法选择（TMM、library size、quantile）。适用于 ChIP-seq 数据中条件间差异结合位点的全基因组检测。

### 105. slingshot 单细胞轨迹推断 (Bioconductor slingshot)

实现单细胞轨迹推断功能，支持从细胞聚类到轨迹构建的完整流程。可以创建 SlingshotNode（簇节点）、SlingshotEdge（边）和 SlingshotResult（轨迹结果）进行轨迹管理。支持基于簇质心的最小生成树（MST）构建，识别细胞间的谱系关系。支持主曲线拟合（principal curve fitting），将 M ST 的每条边拟合为平滑曲线，模拟细胞分化路径。支持拟时间（pseudotime）计算，将细胞投影到最近的主曲线上，获得沿谱系的连续时间表示。支持按拟时间排序细胞和分支点检测。适用于单细胞 RNA-seq 数据中细胞分化轨迹和谱系关系的推断。

### 106. SCnorm 单细胞RNA-seq归一化 (Bioconductor SCnorm)

实现单细胞 RNA-seq 归一化功能，采用分位数回归校正深度依赖偏差。可以创建 SCnormQuantFit（分位数拟合结果）、SCnormGeneNormResult（基因归一化结果）和 SCnormResult（完整归一化结果）。支持库大小计算和分位数回归拟合（通过网格搜索优化分位数损失函数）。支持基因特异性归一化，通过对每个基因单独拟合深度-表达关系，校正测序深度导致的偏差。支持多种归一化方法选择（Quantile、Loess、Spline）和基因过滤（表达阈值、细胞比例）。支持完整的归一化流水线，从原始计数到归一化矩阵的全流程处理。适用于单细胞 RNA-seq 数据的深度依赖偏差校正和标准化。

### 107. EDASeq RNA-seq探索性分析 (Bioconductor EDASeq)

实现 RNA-seq 探索性数据分析功能，支持 GC 含量和基因长度归一化。可以创建 EDASeqGeneAnno（基因注释，包含 GC 含量和有效长度）和 EDASeqDataSet（表达数据集）。支持基因注释创建（GC 含量、长度、有效长度）和 RPKM 计算。支持样本内归一化（within-lane）：基于 GC 含量的 Loess 平滑校正和基于基因长度的 Loess 校正。支持样本间归一化（between-lane）：中位数比率法（DESeq2-style）、分位数归一化和文库大小归一化。支持完整的归一化流水线（先内后外）和多种归一化方法选择。适用于 RNA-seq 数据的探索性分析和校正系统偏差。

### 108. SearchIO 统一搜索结果模型 (Bio.SearchIO)

提供统一的搜索结果层次结构，实现了 QueryResult → Hit → HSP 的三级数据模型。支持 BLAST tabular 格式解析、HMMER3 结果解析和 BLAT PSL 格式解析。支持按 E-value、identity 过滤命中，获取 top hits 和 HSP 统计。支持 BLAST 结果转换为 SearchIO 模型，便于不同搜索工具结果的统一处理。

### 109. PDB Vectors 3D向量与旋转矩阵 (Bio.PDB.vectors)

实现 3D 空间向量运算和旋转矩阵计算，支持向量加减法、点积、叉积、范数、归一化、距离和夹角计算。支持绕 X/Y/Z 轴旋转矩阵生成、轴角旋转、向量变换。实现 Kabsch 算法（SVD 结构叠合），支持最优旋转矩阵求解和 RMSD 计算。支持蛋白质结构的二面角计算和质心计算。

### 110. CircSeq 环状DNA序列操作 (Bio.SeqUtils.CircSeq)

实现环状 DNA 序列的完整操作框架，支持环形坐标包装（wrapping）、环形子序列提取（circ_subseq、circ_slice）、序列旋转（circ_rotate）和反向互补。实现限制性酶切位点查找（circ_find、circ_find_all、circ_find_rc）、酶切片段模拟（circ_digest）和常用内切酶数据库。支持 PCR 引物设计（circ_design_primers），包括 Tm 计算和 GC 含量调整。支持序列插入、删除和环状距离计算。

### 111. AlignAbstract 抽象比对类型与统计 (Bio.Align.AlignAbstract)

提供抽象比对类型的验证和统计分析工具。支持比对验证（序列长度一致性、有效字符检查、类型校验）。实现一致性序列计算（基于阈值的保守位点判定）、Shannon 熵多样性分析、每条列的统计量（保守性、多样性、空位比例、唯一字符数）。支持成对同一性矩阵、总体同一性计算、p-distance 距离矩阵。支持可变位点检测、简约信息位点识别、单例位点检测。支持比对覆盖度计算、按空位比例过滤列和比对修剪。

### 112. SeqUtils 高级功能 (Bio.SeqUtils)

提供序列高级分析功能，支持 GC/AT 滑动窗口偏斜分析、ORF 预测、序列相似度计算等。GC 偏斜和 AT 偏斜采用滑动窗口方法，可用于基因组组成分析和复制起点预测。ORF 预测支持多阅读框扫描、自定义起始和终止密码子检测。序列相似度支持 Hamming 距离（逐位比较）和 Levenshtein 编辑距离（动态规划算法），适用于序列比对和相似性搜索。支持从 DNA 序列中提取所有可能的 ORF，包括框架位置、长度和终止密码子信息。

### 113. Motifs 高级功能 (Bio.motifs)

提供序列模体的高级分析功能，包括信息含量（IC）计算、序列 Logo 数据生成和模体富集分析。信息含量计算基于位置权重矩阵，用于评估每个位置的保守程度。序列 Logo 数据生成支持从多序列比对中提取每个位置的频率和信息含量，适用于可视化转录因子结合位点和保守模体。模体富集分析支持在目标序列集合中查找特定模体的出现频率，并与背景模型进行比较，计算富集显著性。支持 Pearson 相关性比较两个模体矩阵的相似性。

### 114. PDB 结构分析高级功能 (Bio.PDB.StructureAnalysis)

提供蛋白质结构的深度分析功能，包括 SASA（溶剂可及表面积）计算、Ramachandran 质量评估和 Kyte-Doolittle 疏水性分析。SASA 计算采用 Shrake-Rupley 算法，通过在原子表面生成采样点来估计可及表面积，支持自定义探针半径和采样点数。Ramachandran 质量评估将二面角分类到允许区域（核心区域、允许区域、 generously allowed 区域和不允许区域），用于评估结构质量。Kyte-Doolittle 疏水性分析支持基于滑动窗口的蛋白质疏水性图谱计算，适用于跨膜区域预测和蛋白质结构域识别。支持基于理化性质的序列属性距离矩阵生成。

### 115. maftools 癌症基因组学MAF分析 (Bioconductor maftools)

实现癌症基因组学 MAF（Mutation Annotation Format）文件的解析与突变分析功能。支持 MAFMutation 数据结构存储突变信息（基因符号、染色体位置、参考/变异等位基因、样本条码、变异分类）。支持突变类型自动分类：SNV（单核苷酸变异）、Indel（插入缺失）和 Complex（复杂变异）。支持转换（Transition: A↔G、C↔T）和颠换（Transversion）判定。支持 TMB（肿瘤突变负荷）计算，包括总突变数、编码区大小和每样本 TMB。支持突变谱分析（统计 SNV/Indel/Complex 数量和转换/颠换计数）。支持基因/样本级别的突变计数。支持共现分析（co-occurrence 和 mutual exclusivity）。支持 Oncoplot 数据生成（突变矩阵、基因和样本列表）。支持 MAF 文件内容解析（TSV 格式）。适用于癌症基因组学研究和体细胞突变分析。

### 116. CNVkit 拷贝数变异检测 (Bioconductor CNVkit)

实现拷贝数变异（Copy Number Variation）检测与分析功能。支持 CNVProbe 数据结构存储探针信息（染色体、位置、log2 比率）。支持 CBS（Circular Binary Segmentation，循环二元分割）算法进行拷贝数分段，通过递归分割找到拷贝数变化的断点。支持 CNVSegment 数据结构存储分段结果（染色体、起止位置、探针数、平均 log2 比率、拷贝数状态）。支持拷贝数状态判定（deletion/neutral/amplification），基于 log2 比率阈值。支持 log2 比率平滑（滑动窗口均值滤波），降低噪声。支持断点检测（基于分段间 log2 比率差异）。支持 CNVDataset 数据集管理和 CBS 结果汇总。适用于癌症基因组拷贝数变异分析和 CNV 区域识别。

### 117. destiny 单细胞扩散映射降维 (Bioconductor destiny)

实现单细胞数据的扩散映射（Diffusion Maps）降维分析功能。支持 CellData 数据结构存储单细胞表达数据。支持距离矩阵计算（欧氏距离），用于衡量细胞间的相似性。支持高斯核（Gaussian Kernel）构建，将距离转换为相似性矩阵。支持核矩阵归一化（Markov 矩阵），使每行和为 1，构建概率转移矩阵。支持特征分解（Eigen Decomposition），计算特征值和特征向量。支持扩散分量（Diffusion Components）计算，将细胞投影到低维空间。支持解释方差（Explained Variance）计算，评估每个分量的信息含量。支持嵌入坐标提取和可视化数据生成。适用于单细胞 RNA-seq 数据的非线性降维和轨迹分析。

### 118. tradeSeq 轨迹差异表达分析 (Bioconductor tradeSeq)

实现基于轨迹的差异表达分析功能，用于识别沿细胞分化轨迹表达变化的基因。支持 TrajectoryPoint 和 GeneExpressionData 数据结构存储拟时间和基因表达数据。支持 GAM（Generalized Additive Model，广义可加模型）拟合，使用样条基函数（Cubic Spline Basis）对基因表达沿拟时间进行平滑拟合。支持最小二乘法系数估计、拟合值计算和残差分析。支持 R-squared 计算评估模型拟合质量。支持条件效应差异检验（trade_test_condition_effect），比较不同条件下基因沿轨迹的表达差异。支持完整分析流程（run_tradeseq_analysis），包含批量基因检验和 Benjamini-Hochberg 多重检验校正。支持基因表达平滑曲线计算（calculate_gene_smooth）。适用于单细胞轨迹分析中差异表达基因的识别。


## 性能优化

### 优化策略

1. **O(1) 互补查找**: 使用 `FixedArray[UInt16]` 实现直接索引的互补碱基查找表，避免 Map 查找开销
2. **无分配字符串操作**: 使用 `unsafe_get` 和 `FixedArray` 避免中间字符串分配
3. **密码子快速查表**: 使用 2-bit 编码（A=0, C=1, G=2, T=3）构建 64 项密码子翻译表，通过位运算 `(b0 << 4) | (b1 << 2) | b2` 计算索引
4. **分支预测优化**: 使用 `combined = b0 | b1 | b2` 合并多条件判断，提升 CPU 分支预测效率
5. **单遍扫描**: 反向互补和翻译操作使用单遍扫描，避免中间字符串分配
6. **预分配数组**: 字符串分割函数先统计分隔符数量，再一次性分配数组
7. **FASTA 快速索引**: 通过 .fai 索引实现 O(1) 随机访问，跳过逐行解析
8. **FixedArray 替代 Array**: BAM 解析中使用 `FixedArray` 存储参考序列和 CIGAR 元素，减少动态内存分配开销

### 性能基准测试

所有测试基于 MoonBit WasmGC 后端。

#### 序列操作性能

| 操作 | 序列长度 | 单次耗时 | 备注 |
| :--- | :---: | :---: | :--- |
| translate | 100,000 | **601 us** | 优化后，提升 32% |
| complement | 100,000 | **836 us** | |
| reverse_complement | 100,000 | **556 us** | |
| transcribe | 100,000 | **485 us** | |
| reverse | 100,000 | **415 us** | |
| upper | 100,000 | **487 us** | |
| lower | 100,000 | **1029 us** | |
| replace A->N | 100,000 | **1453 us** | |
| count AAA | 100,000 | **436 us** | |
| contains GGG | 100,000 | **221 us** | |
| find ACGT | 100,000 | **0.65 us** | |
| slice 0..1000 | 100,000 | **16.6 us** | |

#### SeqIO 性能

| 操作 | 数据规模 | 单次耗时 |
| :--- | :--- | :---: |
| parse_fasta | 100 records × 1000 bp | **1452 us** |
| write_fasta | 100 records × 1000 bp | **1695 us** |
| parse_fastq | 100 records × 1000 bp | **2497 us** |
| write_fastq | 100 records × 1000 bp | **2236 us** |

#### SeqUtils 性能

| 操作 | 序列长度 | 单次耗时 |
| :--- | :---: | :---: |
| gc (GC含量) | 10,000 | **92 us** |
| crc32 | 10,000 | **97 us** |
| seguid | 10,000 | **7370 us** |
| tm_wallace (熔解温度) | 10,000 | **81 us** |
| ProteinAnalysis::count_amino_acids | 99 aa | **24 us** |
| ProteinAnalysis::molecular_weight | 99 aa | **13 us** |

#### Bio.Phylo 性能

| 操作 | 规模 | 单次耗时 |
| :--- | :--- | :---: |
| parse_newick | simple tree | **5.7 us** |
| Tree::count_terminals | simple tree | **0.61 us** |
| Tree::distance | simple tree | **1.85 us** |
| parse_newick | 50 leaves | **54.7 us** |

#### Bio.PDB 性能

| 操作 | 规模 | 单次耗时 |
| :--- | :--- | :---: |
| parse_pdb | 300 atoms | **619 us** |
| write_pdb | 300 atoms | **1804 us** |
| Atom::distance | 2 atoms | **0.30 us** |

#### BAM/CRAM 性能

| 操作 | 规模 | 单次耗时 |
| :--- | :--- | :---: |
| parse_bam_base64 | 1 record | **7.19 us** |
| BamFile::to_sam | 1 record | **4.23 us** |
| BamFile::get_reference_names | 1 reference | **1.37 us** |
| is_cram_magic (positive) | 8 bytes | **0.26 us** |
| is_cram_magic (negative) | 8 bytes | **0.12 us** |

#### Python 对比基准测试

> **注意**: BioPython 的 `complement`、`count`、`replace` 等操作使用 C 扩展实现，因此 MoonBit 在这些操作上较慢是预期的。MoonBit 的优势在于纯 WebAssembly 环境下的便携性和 `translate` 等算法密集型操作的性能。

| 操作 | 序列长度 | Python (BioPython) | MoonBit (WasmGC) | 相对性能 |
| :--- | :---: | :---: | :---: | :---: |
| complement | 100,000 | 106 us | 836 us | ~8× slower |
| reverse_complement | 100,000 | 178 us | 556 us | ~3× slower |
| transcribe | 100,000 | 244 us | 485 us | ~2× slower |
| translate | 100,000 | 17.5 ms | 601 us | **~29× faster** |
| count | 100,000 | 98 us | 436 us | ~4× slower |
| find | 100,000 | 1.1 us | 0.65 us | **~1.7× faster** |
| replace | 100,000 | 234 us | 1453 us | ~6× slower |

## 测试验证

### 测试统计

| 指标 | 数值 |
| :--- | :---: |
| 总测试数 | 2229 |
| 通过数 | 2229 |
| 失败数 | 0 |
| 通过率 | 100% |

### 运行测试

```bash
# 构建项目
moon build

# 运行所有测试
moon test --package IvanAXu/BioSeqs/test/moonbit

# 运行单个模块测试
moon test --package IvanAXu/BioSeqs/test/moonbit --test bio_seq_test

# 更新快照测试
moon test --update
```

### 测试模块分布

| 模块 | 测试文件 | 测试数 |
| :--- | :--- | :---: |
| 序列核心 | `bio_seq_test.mbt` | 28 |
| 序列工具 | `sequtils_test.mbt` | 16 |
| 序列特征 | `seqfeature_test.mbt` | 8 |
| SeqIO | `seqio_wb_test.mbt` | 20 |
| AlignIO | `bio_seq_wb_test.mbt` | 0 |
| 系统发育树 | `phylo_test.mbt` | 22 |
| PDB 结构 | `pdb_test.mbt` | 19 |
| 比对算法 | `alignment_test.mbt` | 13 |
| SAM 解析 | `sam_test.mbt` | 6 |
| VCF 解析 | `vcf_test.mbt` | 7 |
| BAM 解析 | `bam_test.mbt` | 12 |
| BGZF 解压缩 | `bgzf_test.mbt` | 2 |
| FASTA 索引 | `faidx_test.mbt` | 9 |
| 特征提取 | `feature_extraction_test.mbt` | 19 |
| Biostrings | `biostrings_test.mbt` | 21 |
| GenomicRanges | `genomic_ranges_test.mbt` | 22 |
| DESeq2 | `deseq2_test.mbt` | 10 |
| dplyr | `dplyr_test.mbt` | 9 |
| Smith-Waterman | `smith_waterman_test.mbt` | 16 |
| Needleman-Wunsch | `needleman_wunsch_test.mbt` | 20 |
| Bloom Filter & k-mer | `bloom_filter_test.mbt` | 30 |
| BWT + FM-index | `bwt_fm_test.mbt` | 28 |
| De Bruijn Graph | `de_bruijn_test.mbt` | 12 |
| Suffix Array & Tree | `suffix_array_tree_test.mbt` + `suffix_array_tree_wbtest.mbt` | 46 |
| 替换矩阵 | `subsmat_test.mbt` | 13 |
| BLAST解析 | `blast_test.mbt` | 10 |
| SearchIO | `search_io_test.mbt` | 7 |
| Overlap-Layout-Consensus | `olc_test.mbt` | 12 |
| Hidden Markov Model | `hmm_test.mbt` + `hmm_wbtest.mbt` | 30 |
| K-means Clustering | `kmeans_test.mbt` + `kmeans_wbtest.mbt` | 30 |
| CRAM 解析 | `cram_wbtest.mbt` | 1 |
| Sequence Motifs | `motifs_test.mbt` | 20 |
| Sequence Utils | `seq_utils_test.mbt` | 24 |
| Population Genetics | `popgen_test.mbt` | 13 |
| edgeR | `edger_test.mbt` | 7 |
| limma | `limma_test.mbt` | 10 |
| SummarizedExperiment | `summarized_experiment_test.mbt` | 7 |
| IRanges | `iranges_test.mbt` | 14 |
| AlignIO | `align_io_test.mbt` | 12 |
| Cluster | `cluster_test.mbt` | 12 |
| Restriction | `restriction_test.mbt` | 16 |
| TreeIO | `tree_io_test.mbt` | 8 |
| TxDb | `txdb_test.mbt` | 14 |
| ProtParam | `protparam_test.mbt` | 15 |
| rtracklayer | `rtracklayer_test.mbt` | 19 |
| CodonUsage | `codon_usage_test.mbt` | 8 |
| GenomicAlignments | `genomic_alignments_test.mbt` | 7 |
| VariantAnnotation | `variant_annotation_test.mbt` | 6 |
| Affy | `affy_test.mbt` | 8 |
| SVDSuperimposer | `svd_superimposer_test.mbt` | 9 |
| GOEnrichment | `go_enrichment_test.mbt` | 8 |
| SingleCell | `single_cell_test.mbt` | 9 |
| KEGG | `kegg_test.mbt` | 9 |
| Medline | `medline_test.mbt` | 8 |
| BSGenome | `bsgenome_test.mbt` | 6 |
| BioMart | `biomart_test.mbt` | 5 |
| RUVSeq | `ruvseq_test.mbt` | 5 |
| PDB 高级分析 | `pdb_analysis_test.mbt` | 13 |
| Phylo 高级分析 | `phylo_metrics_test.mbt` | 11 |
| 序列复杂度 | `seq_complexity_test.mbt` | 18 |
| GenomeInfoDb | `genome_info_db_test.mbt` | 13 |
| InteractionSet | `interaction_set_test.mbt` | 14 |
| MultiAssayExperiment | `multi_assay_experiment_test.mbt` | 16 |
| TreeConstruction | `tree_construction_test.mbt` | 7 |
| NeighborSearch | `neighbor_search_test.mbt` | 6 |
| SwissProt | `swissprot_test.mbt` | 8 |
| mmCIF | `mmcif_test.mbt` | 2 |
| Nexus | `nexus_test.mbt` | 2 |
| EMBOSS | `emboss_test.mbt` | 15 |
| ChIPseeker | `chipseeker_test.mbt` | 14 |
| DOSE | `dose_test.mbt` | 5 |
| ReactomePA | `reactome_pa_test.mbt` | 5 |
| AnnotationDbi | `annotation_dbi_test.mbt` | 5 |
| clusterProfiler | `cluster_profiler_test.mbt` | 5 |
| WGCNA | `wgcna_test.mbt` | 5 |
| Biobase | `biobase_test.mbt` | 5 |
| GEOquery | `geoquery_test.mbt` | 5 |
| tximport | `tximport_test.mbt` | 5 |
| AnnotationHub | `annotation_hub_test.mbt` | 8 |
| GenomicFeatures | `genomic_features_test.mbt` | 6 |
| graph | `graph_test.mbt` | 8 |
| DropletUtils | `droplet_utils_test.mbt` | 6 |
| scran | `scran_test.mbt` | 8 |
| monocle3 | `monocle3_test.mbt` | 10 |
| ShortRead | `short_read_test.mbt` | 15 |
| scater | `scater_test.mbt` | 17 |
| MAST | `mast_test.mbt` | 12 |
| GenomicFiles | `genomic_files_test.mbt` | 28 |
| DiffBind | `diffbind_test.mbt` | 36 |
| minfi | `minfi_test.mbt` | 40 |
| flowCore | `flow_core_test.mbt` | 35 |
| bsseq | `bsseq_test.mbt` | 54 |
| SingleCellExperiment | `single_cell_experiment_test.mbt` | 58 |
| ComplexHeatmap | `complex_heatmap_test.mbt` | 20 |
| GSVA | `gsva_test.mbt` | 15 |
| ChromVAR | `chromvar_test.mbt` | 12 |
| DelayedArray | `delayed_array_test.mbt` | 10 |
| AnnotationFilter | `annotation_filter_test.mbt` | 10 |
| scDblFinder | `sc_dbl_finder_test.mbt` | 8 |
| ChIPseeker | `chipseeker_test.mbt` | 14 |
| Taxonomy | `taxonomy_test.mbt` | 7 |
| GFF | `gff_test.mbt` | 5 |
| Phylo.Consensus | `phylo_consensus_test.mbt` | 3 |
| maftools | `maftools_test.mbt` | 18 |
| CNVkit | `cnvkit_test.mbt` | 15 |
| destiny | `destiny_test.mbt` | 12 |
| tradeSeq | `tradeseq_test.mbt` | 12 |

### Python 对比测试

运行 Python 参考脚本验证结果一致性：

```bash
# 序列核心功能
python test/python/python_reference.py

# 序列 I/O
python test/python/python_seqio_reference.py

# 比对 I/O
python test/python/python_alignio_reference.py

# 序列工具
python test/python/python_sequtils_reference.py

# 序列特征
python test/python/python_seqfeature_reference.py

# 系统发育树
python test/python/python_phylo_reference.py

# PDB 结构
python test/python/python_pdb_reference.py

# 比对算法
python test/python/python_skbio_alignment_reference.py

# 性能基准测试
python test/python/python_bench.py

# SAM/VCF 对比
python test/python/skbio_pysam_compare.py

# pyfaidx 对比
python test/python/pyfaidx_compare.py

# 机器学习特征提取对比
python test/python/python_feature_reference.py

# 自动对比脚本
bash test/python/compare.sh
bash test/python/compare_seqio.sh
```

## 使用方法

### 构建与测试

```bash
# 构建项目
moon build

# 运行所有测试
moon test --package IvanAXu/BioSeqs/test/moonbit

# 更新接口文件
moon info

# 格式化代码
moon fmt

# 查看测试覆盖率
moon coverage analyze > uncovered.log
```

### 命令行工具

```bash
# Seq 测试工具
moon run cmd/main/main.mbt

# SeqIO 测试工具
moon run cmd/seqio_main/main.mbt

# AlignIO 测试工具
moon run cmd/alignio_main/main.mbt

# 性能基准测试
moon run cmd/bench/main.mbt
```

### 示例程序

项目提供 88 个示例程序，展示各模块的典型用法：

| 示例 | 说明 | 运行命令 |
|------|------|----------|
| basic_seq | 序列操作（互补、转录、翻译、GC含量、熔解温度、校验和） | `moon run examples/basic_seq/main.mbt` |
| seqio_demo | 序列 I/O（FASTA/FASTQ/GenBank 解析与写入、FASTA 索引） | `moon run examples/seqio_demo/main.mbt` |
| alignment_demo | 序列比对（Needleman-Wunsch、Smith-Waterman、Clustal/Phylip） | `moon run examples/alignment_demo/main.mbt` |
| phylo_demo | 系统发育树（Newick 解析、遍历、距离计算、ASCII 可视化） | `moon run examples/phylo_demo/main.mbt` |
| pdb_demo | PDB 结构解析（原子/残基/链访问、距离计算） | `moon run examples/pdb_demo/main.mbt` |
| sam_vcf_demo | SAM/VCF 解析（比对记录、变异检测、基因型查询） | `moon run examples/sam_vcf_demo/main.mbt` |
| faidx_demo | FASTA 索引（pyfaidx 风格随机访问、.fai 序列化） | `moon run examples/faidx_demo/main.mbt` |
| ml_features | 机器学习特征提取（k-mer频率、氨基酸组成、理化性质） | `moon run examples/ml_features/main.mbt` |
| bam_demo | BAM/BGZF 解析（二进制格式、压缩数据读取） | `moon run examples/bam_demo/main.mbt` |
| cram_demo | CRAM 格式解析（压缩二进制序列比对格式、CRAM转BAM、参考序列管理） | `moon run examples/cram_demo/main.mbt` |
| biostrings_demo | Biostrings 序列分析（IUPAC、RSCU、复杂度、Tm） | `moon run examples/biostrings_demo/main.mbt` |
| genomic_ranges_demo | GenomicRanges 基因组区间操作（GRanges、区间运算、集合操作） | `moon run examples/genomic_ranges_demo/main.mbt` |
| deseq2_demo | DESeq2 差异表达分析（size factors归一化、分散度估计、负二项GLM拟合、Wald检验、LFC收缩） | `moon run examples/deseq2_demo/main.mbt` |
| dplyr_demo | dplyr 数据操作（filter、select、mutate、arrange、group_by、summarize、join） | `moon run examples/dplyr_demo/main.mbt` |
| smith_waterman_demo | Smith-Waterman 局部序列比对（DNA/蛋白质比对、自定义打分、得分矩阵） | `moon run examples/smith_waterman_demo/main.mbt` |
| needleman_wunsch_demo | Needleman-Wunsch 全局序列比对（DNA/蛋白质比对、自定义打分、得分矩阵） | `moon run examples/needleman_wunsch_demo/main.mbt` |
| bloom_filter_demo | Bloom Filter & k-mer 计数（成员查询、精确计数、近似去重、FPR 对比） | `moon run examples/bloom_filter_demo/main.mbt` |
| bwt_fm_demo | BWT + FM-index（BWT 变换、逆变换、后缀数组、模式匹配、LF 映射） | `moon run examples/bwt_fm_demo/main.mbt` |
| de_bruijn_demo | De Bruijn Graph 序列组装（图构建、欧拉路径、多重读段组装、分支检测） | `moon run examples/de_bruijn_demo/main.mbt` |
| suffix_array_tree_demo | Suffix Array & Suffix Tree（前缀倍增算法、LCP数组、模式匹配、最长重复子串） | `moon run examples/suffix_array_tree_demo/main.mbt` |
| olc_demo | Overlap-Layout-Consensus 序列组装（重叠检测、哈密顿路径、一致性序列生成、Graphviz 可视化） | `moon run examples/olc_demo/main.mbt` |
| hmm_demo | Hidden Markov Model 基因预测（前向算法、后向算法、维特比算法、Baum-Welch训练、外显子提取） | `moon run examples/hmm_demo/main.mbt` |
| kmeans_demo | K-means 聚类分析（距离计算、K-means++初始化、聚类、轮廓系数评估、最优k值搜索） | `moon run examples/kmeans_demo/main.mbt` |
| motifs_demo | 序列模体识别（PWM构建、模体搜索、一致性序列、信息含量、序列Logo生成、模体富集、Pearson相关性） | `moon run examples/motifs_demo/main.mbt` | 
| seq_utils_demo | 序列工具函数（分子量计算、GC含量、Tm值、序列类型检测、氨基酸转换、k-mer计数、GC/AT滑动窗口偏斜、ORF预测、Hamming距离、Levenshtein编辑距离） | `moon run examples/seq_utils_demo/main.mbt` | 
| popgen_demo | 群体遗传学分析（等位基因频率、基因型频率、哈迪-温伯格检验、FST统计、Watterson's theta、Shannon/Simpson多样性指数、Chao1/ACE丰富度估计） | `moon run examples/popgen_demo/main.mbt` | 
| edger_demo | edgeR 差异表达分析（DGEList创建、归一化因子、分散度估计、精确检验、GLM拟合） | `moon run examples/edger_demo/main.mbt` | 
| limma_demo | limma 差异表达分析（voom变换、线性模型拟合、经验贝叶斯、topTable、对比矩阵） | `moon run examples/limma_demo/main.mbt` | 
| summarized_experiment_demo | SummarizedExperiment 多维数据容器（Assays、行/列操作、合并） | `moon run examples/summarized_experiment_demo/main.mbt` | 
| iranges_demo | IRanges 整数区间操作（shift、resize、reduce、集合运算、重叠检测） | `moon run examples/iranges_demo/main.mbt` |
| align_io_demo | 比对格式解析（ClustalW、FASTA、Stockholm格式解析与写入） | `moon run examples/align_io_demo/main.mbt` |
| cluster_demo | 序列聚类分析（距离矩阵、层次聚类、Newick输出、轮廓系数） | `moon run examples/cluster_demo/main.mbt` |
| restriction_demo | 限制性内切酶分析（酶切位点查找、序列酶切、片段分析） | `moon run examples/restriction_demo/main.mbt` |
| txdb_demo | TxDb 转录本数据库（GTF解析、基因/转录本/外显子/CDS提取、UTR/内含子计算、启动子提取） | `moon run examples/txdb_demo/main.mbt` |
| subsmat_demo | 替换矩阵（BLOSUM62/45、PAM250/30矩阵查询、蛋白质比对打分） | `moon run examples/subsmat_demo/main.mbt` |
| blast_demo | BLAST结果解析（tabular/xml格式、HSP过滤、E-value/identity过滤、最佳匹配） | `moon run examples/blast_demo/main.mbt` |
| search_io_demo | SearchIO 统一搜索结果（HMMER3 tabular解析、BLAT PSL解析、BLAST转换、top_hits、count_hsps） | `moon run examples/search_io_demo/main.mbt` |
| protparam_demo | ProtParam 蛋白质参数分析（分子量、不稳定指数、GRAVY评分、等电点、信号肽预测、二级结构倾向） | `moon run examples/protparam_demo/main.mbt` |
| rtracklayer_demo | rtracklayer 基因组轨道格式（BED/WIG/BEDGraph/GFF解析与写入、GRanges转换） | `moon run examples/rtracklayer_demo/main.mbt` |
| codon_usage_demo | CodonUsage 密码子使用分析（CAI、ENC、RSCU、GC3、CBI、Fop、最优密码子检测、密码子翻译） | `moon run examples/codon_usage_demo/main.mbt` |
| genomic_alignments_demo | GenomicAlignments 基因组比对分析（GAlignments创建、覆盖度计算、summarizeOverlaps、pileup、过滤） | `moon run examples/genomic_alignments_demo/main.mbt` |
| variant_annotation_demo | VariantAnnotation 变异注释（变异类型检测、定位、编码效应预测、变异汇总） | `moon run examples/variant_annotation_demo/main.mbt` |
| affy_demo | Affy Affymetrix芯片数据分析（AffyBatch创建、背景校正、PM-MM差异计算、RMA标准化、探针集汇总） | `moon run examples/affy_demo/main.mbt` |
| svd_superimposer_demo | SVDSuperimposer SVD蛋白质结构叠合（原子坐标创建、结构叠合、旋转矩阵、平移向量、RMSD计算） | `moon run examples/svd_superimposer_demo/main.mbt` |
| go_enrichment_demo | GOEnrichment GO功能富集分析（GOTerm创建、超几何检验、BH校正、富集结果过滤、命名空间筛选） | `moon run examples/go_enrichment_demo/main.mbt` |
| single_cell_demo | SingleCell 单细胞数据分析（QC指标计算、细胞过滤、Log标准化、PCA降维、高变异基因检测） | `moon run examples/single_cell_demo/main.mbt` |
| kegg_demo | KEGG数据库解析（基因、通路、化合物、酶记录解析、通路分析、基因ID提取） | `moon run examples/kegg_demo/main.mbt` |
| medline_demo | Medline/PubMed解析（文献记录、APA引用格式、MeSH过滤、作者过滤、年份统计） | `moon run examples/medline_demo/main.mbt` |
| bsgenome_demo | BSgenome 基因组序列数据库（基因组创建、染色体管理、序列提取、链特异性基因提取） | `moon run examples/bsgenome_demo/main.mbt` |
| biomart_demo | biomaRt 基因ID转换和注释查询（ID映射、基因注释、外部数据库查询、批量查询） | `moon run examples/biomart_demo/main.mbt` |
| ruvseq_demo | RUVSeq RNA-seq批次效应去除（数据标准化、Log2转换、RUVg/RUVs方法、批次效应校正） | `moon run examples/ruvseq_demo/main.mbt` |
| pdb_analysis_demo | PDB高级结构分析（二面角、距离矩阵、接触图、氢键检测、二级结构分配、回转半径、Ramachandran图、SASA计算、疏水性分析、序列属性矩阵） | `moon run examples/pdb_analysis_demo/main.mbt` |
| phylo_metrics_demo | 系统发育树高级分析（总分支长度、最大深度、Colless指数、系统发生距离、Robinson-Foulds距离） | `moon run examples/phylo_metrics_demo/main.mbt` |
| seq_complexity_demo | 序列复杂度与组成分析（Shannon熵、语言学复杂度、DUST评分、CGR、序列相似度） | `moon run examples/seq_complexity_demo/main.mbt` |
| align_info_demo | AlignInfo 比对统计（一致性序列、保守位点、Shannon熵、成对序列同一性） | `moon run examples/align_info_demo/main.mbt` |
| codon_align_demo | CodonAlign 密码子比对（密码子替换分类、dN/dS选择压力分析、密码子使用偏好、ENC） | `moon run examples/codon_align_demo/main.mbt` |
| entrez_demo | Entrez NCBI数据库访问（ESearch、EFetch、PubMed/Gene/Taxonomy解析） | `moon run examples/entrez_demo/main.mbt` |
| genome_info_db_demo | GenomeInfoDb 基因组信息管理（染色体信息、着丝粒位置、染色体臂、基因组构建） | `moon run examples/genome_info_db_demo/main.mbt` |
| interaction_set_demo | InteractionSet 染色质交互（Hi-C交互、锚点对、交互矩阵、距离分布、Top交互） | `moon run examples/interaction_set_demo/main.mbt` |
| multi_assay_experiment_demo | MultiAssayExperiment 多组学数据协调（实验协调、样本映射、跨实验子集） | `moon run examples/multi_assay_experiment_demo/main.mbt` |
| tree_construction_demo | TreeConstruction 系统发育树构建（UPGMA/WPGMA/NJ算法、替换模型、距离矩阵） | `moon run examples/tree_construction_demo/main.mbt` |
| neighbor_search_demo | NeighborSearch KD树近邻搜索（半径搜索、最近邻、原子对搜索） | `moon run examples/neighbor_search_demo/main.mbt` |
| swissprot_demo | SwissProt 蛋白数据库解析（记录解析、特征提取、参考文献） | `moon run examples/swissprot_demo/main.mbt` |
| fgsea_demo | FGSEA 快速基因集富集分析（基因排名、富集分数、NES、p值、Leading Edge基因、BH校正） | `moon run examples/fgsea_demo/main.mbt` |
| sva_demo | SVA 替代变量分析与ComBat批次校正（经验贝叶斯方法、PCA分析、批次效应去除） | `moon run examples/sva_demo/main.mbt` |
| ballgown_demo | Ballgown 转录组水平差异表达分析（FPKM计算、t检验、转录本/基因水平DE分析） | `moon run examples/ballgown_demo/main.mbt` |
| mmcif_demo | mmCIF格式解析（数据块解析、类别查询、原子位点提取、结构信息） | `moon run examples/mmcif_demo/main.mbt` |
| nexus_demo | Nexus格式解析（数据矩阵、系统发育树、距离矩阵、分类单元） | `moon run examples/nexus_demo/main.mbt` |
| emboss_demo | EMBOSS工具接口（GC偏斜、AT偏斜、分子量、Tm值、ORF查找、距离计算） | `moon run examples/emboss_demo/main.mbt` |
| bioconductor_demo | Bioconductor模块综合示例（ChIPseeker峰注释(外显子/内含子/UTR分类、peak2gene关联)、DOSE疾病富集、ReactomePA通路分析、AnnotationDbi注释数据库、clusterProfiler富集框架、WGCNA共表达网络、Batchelor单细胞批次校正、Seurat单细胞分析） | `moon run examples/bioconductor_demo/main.mbt` |
| short_read_demo | ShortRead 短读序列质量控制（QA统计、adapter修剪、质量修剪、读长过滤、FastQC报告生成） | `moon run examples/short_read_demo/main.mbt` |
| scater_demo | scater 单细胞质量控制（QC指标计算、细胞/基因过滤、CPM/log-CPM标准化、HVG检测、PCA降维） | `moon run examples/scater_demo/main.mbt` |
| mast_demo | MAST 单细胞差异表达分析（Hurdle模型、离散/连续检验、BH-FDR校正、结果汇总） | `moon run examples/mast_demo/main.mbt` |
| genomic_files_demo | GenomicFiles 分布式基因组文件处理（BAM/BED/VCF扫描、区间查询、归约、覆盖度计算） | `moon run examples/genomic_files_demo/main.mbt` |
| diffbind_demo | DiffBind ChIP-seq差异结合分析（峰值重叠、共识峰识别、TMM归一化、负二项分布检验） | `moon run examples/diffbind_demo/main.mbt` |
| minfi_demo | minfi DNA甲基化分析（NOOB/Illumina/分位数/功能归一化、β/M值计算、DMP/DMR分析） | `moon run examples/minfi_demo/main.mbt` |
| flow_core_demo | flowCore 流式细胞术（FCS文件处理、数据变换、荧光补偿、矩形/多边形/椭球/四象限门控） | `moon run examples/flow_core_demo/main.mbt` |
| bsseq_demo | bsseq 亚硫酸氢盐测序分析（BSmooth平滑、DMR检测、CpG合并、甲基化率计算） | `moon run examples/bsseq_demo/main.mbt` |
| single_cell_experiment_demo | SingleCellExperiment 单细胞核心容器（多assay管理、PCA/t-SNE/UMAP降维、size factors、归一化） | `moon run examples/single_cell_experiment_demo/main.mbt` |
| complex_heatmap_demo | ComplexHeatmap 复杂热图可视化（行/列聚类、颜色映射、热图注释、分组拆分、ASCII热图） | `moon run examples/complex_heatmap_demo/main.mbt` |
| gsva_demo | GSVA 基因集变异分析（ssGSEA/zscore/PLAGE评分、富集分析、置换检验、统计摘要、富集图可视化、表型相关性分析、生存分析、分数分布分析） | `moon run examples/gsva_demo/main.mbt` |
| chromvar_demo | ChromVAR 染色质变异分析（TF motif富集、GC偏差校正、细胞聚类、变异性分析、偏差图） | `moon run examples/chromvar_demo/main.mbt` |
| delayed_array_demo | DelayedArray 延迟计算数组（懒加载操作、分块处理、行/列聚合、转置、子集操作） | `moon run examples/delayed_array_demo/main.mbt` |
| annotation_filter_demo | AnnotationFilter 基因注释过滤（染色体筛选、生物类型过滤、链过滤、区域重叠检测、符号模式匹配） | `moon run examples/annotation_filter_demo/main.mbt` |
| sc_dbl_finder_demo | scDblFinder 单细胞双细胞检测（Doublet评分计算、最近邻搜索、双细胞检测、PCA降维、细胞过滤） | `moon run examples/sc_dbl_finder_demo/main.mbt` |
| chipseeker_demo | ChIPseeker ChIP-seq峰值注释（基因组区域分类(启动子/外显子/内含子/UTR/基因间区)、距离TSS分布、BED格式读取、peak2gene关联分析、多峰值集重叠分析、Venn图、饼图可视化、统计分析） | `moon run examples/chipseeker_demo/main.mbt` |
| taxonomy_demo | Taxonomy 分类学分析（分类数据库创建、谱系查询、共同祖先计算、分类单元管理） | `moon run examples/taxonomy_demo/main.mbt` |
| gff_demo | GFF GFF3格式解析（基因注释特征提取、属性解析、基因/转录本/CDS/外显子结构分析） | `moon run examples/gff_demo/main.mbt` |
| phylo_consensus_demo | Phylo.Consensus 一致性树构建（Newick解析、分裂分析、多数规则树、支持度计算） | `moon run examples/phylo_consensus_demo/main.mbt` |
| phenotype_demo | Bio.phenotype 表型微阵列分析（生长曲线拟合、预测、控制减法、CSV/JSON解析） | `moon run examples/phenotype_demo/main.mbt` |
| blast_applications_demo | Bio.Blast.Applications BLAST命令行工具（8种BLAST变体、快速构建器、参数管理、命令构建） | `moon run examples/blast_applications_demo/main.mbt` |
| psea_demo | Bio.PDB.PSEA 二级结构预测（CA-CA距离、虚拟二面角、H/E/C分配、三态到八态转换） | `moon run examples/psea_demo/main.mbt` |
| sff_io_demo | Bio.SeqIO.SffIO SFF二进制解析（二进制编码/解码往返、质量修剪、按名称查找） | `moon run examples/sff_io_demo/main.mbt` |

## 技术栈

- **语言**: MoonBit 0.4.x
- **目标平台**: WebAssembly (WasmGC)
- **包管理**: moon.mod / moon.pkg
- **测试框架**: MoonBit 内置测试框架
- **参考实现**: Biopython 1.83+, scikit-bio 0.5.7+, pysam 0.22+, pyfaidx 0.7+

## 开发规范

### 代码风格

- 使用 `///|` 分隔代码块
- 公共 API 使用 `pub` 修饰
- 错误类型使用 `suberror` 定义
- 测试文件使用 `_test.mbt` (黑盒) / `_wbtest.mbt` (白盒) 后缀

### 命名约定

- 类型名使用 PascalCase: `SeqRecord`, `Clade`, `TabularMSA`, `Fasta`
- 函数名使用 snake_case: `parse_fasta`, `reverse_complement`, `build_fai`
- 私有函数不使用 `pub` 修饰
- 文件名使用 snake_case: `fasta_io.mbt`, `alignment.mbt`, `faidx.mbt`
- 测试名称使用 `模块_功能` 格式: `seq_complement`, `alignment_global_pairwise_nucleotide`, `faidx_build_index`

## 未来规划/After 0.1.4

- ✅ 实现 BAM 二进制格式解析
- ✅ 实现 BGZF 解压缩功能，支持读取压缩的 BAM 文件
- ✅ 添加 机器学习特征提取功能
- ✅ 优化 translate 性能（从 880 us → 601 us，提升 32%）
- ✅ 实现 Biostrings 序列分析（IUPAC、RSCU、复杂度、Tm）
- ✅ 实现 DESeq2 差异表达分析（数据集创建、结果分析、显著基因筛选）
- ✅ 实现 GenomicRanges 基因组区间操作（GRanges、区间运算、集合操作）
- ✅ 实现 dplyr 数据操作
- ✅ 实现 Smith-Waterman 局部序列比对（DNA/蛋白质比对、自定义打分、得分矩阵）
- ✅ 实现 Needleman-Wunsch 全局序列比对（DNA/蛋白质比对、自定义打分、得分矩阵）
- ✅ 实现 Bloom Filter & k-mer 计数（成员查询、精确计数、近似去重、FPR 对比）
- ✅ 实现 BWT + FM-index（BWT 变换、逆变换、后缀数组、模式匹配、LF 映射）
- ✅ 实现 De Bruijn Graph 序列组装（图构建、欧拉路径、多重读段组装、分支检测）
- ✅ 实现 Suffix Array & Suffix Tree（前缀倍增算法、LCP数组、模式匹配、最长重复子串）
- ✅ 实现 Overlap-Layout-Consensus 序列组装（重叠检测、哈密顿路径、一致性序列生成、Graphviz 可视化）
- ✅ 实现 Hidden Markov Model 基因预测（前向算法、后向算法、维特比算法、Baum-Welch训练、外显子提取）
- ✅ 实现 kmeans_demo | K-means 聚类分析（距离计算、K-means++初始化、聚类、轮廓系数评估、最优k值搜索）
- ✅ 实现 序列模体识别（PWM构建、模体搜索、一致性序列、信息含量、序列Logo生成、模体富集、Pearson相关性）
- ✅ 实现 序列工具函数（分子量计算、GC含量、Tm值、序列类型检测、氨基酸转换、k-mer计数、GC/AT滑动窗口偏斜、ORF预测、Hamming距离、Levenshtein编辑距离）
- ✅ 实现 群体遗传学分析（等位基因频率、基因型频率、哈迪-温伯格检验、FST统计、Watterson's theta）
- ✅ 实现 edgeR 差异表达分析（DGEList创建、归一化因子、分散度估计、精确检验、GLM拟合）
- ✅ 实现 SummarizedExperiment 多维数据容器（Assays、行/列操作、合并）
- ✅ 实现 IRanges 整数区间操作（shift、resize、reduce、集合运算、重叠检测）
- ✅ 实现 比对格式解析（ClustalW、FASTA、Stockholm格式解析与写入）
- ✅ 实现 序列聚类分析（距离矩阵、层次聚类、Newick输出、轮廓系数）
- ✅ 实现 限制性内切酶分析（酶切位点查找、序列酶切、片段分析）
- ✅ 实现 TxDb 转录本数据库（GTF解析、基因/转录本/外显子/CDS提取、UTR/内含子计算、启动子提取）
- ✅ 实现 替换矩阵、BLAST结果解析、ProtParam 蛋白质参数分析、rtracklayer 基因组轨道格式
- ✅ 实现 CodonUsage 密码子使用分析（CAI、ENC、RSCU、GC3、CBI、Fop、最优密码子检测、密码子翻译）
- ✅ 实现 limma 差异表达分析（voom变换、线性模型拟合、经验贝叶斯、topTable、对比矩阵）
- ✅ 实现 SearchIO 统一搜索结果（HMMER3 tabular解析、BLAT PSL解析、BLAST转换、top_hits、count_hsps）
- ✅ 实现 GenomicAlignments 基因组比对分析（GAlignments创建、覆盖度计算、summarizeOverlaps、pileup、过滤）
- ✅ 实现 VariantAnnotation 变异注释（变异类型检测、定位、编码效应预测、变异汇总）
- ✅ 实现 Affy Affymetrix芯片数据分析（AffyBatch创建、背景校正、PM-MM差异计算、RMA标准化、探针集汇总）
- ✅ 实现 SVDSuperimposer SVD蛋白质结构叠合（原子坐标创建、结构叠合、旋转矩阵、平移向量、RMSD计算）
- ✅ 实现 GOEnrichment GO功能富集分析（GOTerm创建、超几何检验、BH校正、富集结果过滤、命名空间筛选）
- ✅ 实现 SingleCell 单细胞数据分析（QC指标计算、细胞过滤、Log标准化、PCA降维、高变异基因检测）
- ✅ 实现 KEGG数据库解析（基因、通路、化合物、酶记录解析、通路分析、基因ID提取）
- ✅ 实现 Medline/PubMed解析（文献记录、APA引用格式、MeSH过滤、作者过滤、年份统计）
- ✅ 实现 BSgenome 基因组序列数据库（基因组创建、染色体管理、序列提取、链特异性基因提取）
- ✅ 实现 biomaRt 基因ID转换和注释查询（ID映射、基因注释、外部数据库查询、批量查询）
- ✅ 实现 RUVSeq RNA-seq批次效应去除（数据标准化、Log2转换、RUVg/RUVs方法、批次效应校正）
- ✅ 实现 PDB高级结构分析（二面角、距离矩阵、接触图、氢键检测、二级结构分配、回转半径、Ramachandran图、SASA计算、疏水性分析、序列属性矩阵）
- ✅ 实现 phylo_metrics 系统发育树高级分析（总分支长度、最大深度、Colless指数、系统发生距离、Robinson-Foulds距离）
- ✅ 实现 seq_complexity 序列复杂度与组成分析（Shannon熵、语言学复杂度、DUST评分、CGR、序列相似度）
- ✅ 实现 AlignInfo 比对统计（一致性序列、保守位点、Shannon熵、成对序列同一性）
- ✅ 实现 CodonAlign 密码子比对（密码子替换分类、dN/dS选择压力分析、密码子使用偏好、ENC）
- ✅ 实现 Entrez NCBI数据库访问（ESearch、EFetch、PubMed/Gene/Taxonomy解析）
- ✅ 实现 GenomeInfoDb 基因组信息管理（染色体信息、着丝粒位置、染色体臂、基因组构建）
- ✅ 实现 InteractionSet 染色质交互（Hi-C交互、锚点对、交互矩阵、距离分布、Top交互）
- ✅ 实现 MultiAssayExperiment 多组学数据协调（实验协调、样本映射、跨实验子集）
- ✅ 实现 TreeConstruction 系统发育树构建（UPGMA/WPGMA/NJ算法、替换模型、距离矩阵）
- ✅ 实现 NeighborSearch KD树近邻搜索（半径搜索、最近邻、原子对搜索）
- ✅ 实现 SwissProt 蛋白数据库解析（记录解析、特征提取、参考文献）
- ✅ 实现 FGSEA 快速基因集富集分析（基因排名、富集分数、NES、p值、Leading Edge基因、BH校正）
- ✅ 实现 SVA 替代变量分析与ComBat批次校正（经验贝叶斯方法、PCA分析、批次效应去除）
- ✅ 实现 Ballgown 转录组水平差异表达分析（FPKM计算、t检验、转录本/基因水平DE分析）
- ✅ 实现 mmCIF格式解析（数据块解析、类别查询、原子位点提取、结构信息）
- ✅ 实现 Nexus格式解析（数据矩阵、系统发育树、距离矩阵、分类单元）
- ✅ 实现 EMBOSS工具接口（GC偏斜、AT偏斜、分子量、Tm值、ORF查找、距离计算）
- ✅ 实现 群体遗传学分析（等位基因频率、基因型频率、哈迪-温伯格检验、FST统计、Watterson's theta、Shannon/Simpson多样性指数、Chao1/ACE丰富度估计）
- ✅ 实现 CRAM 格式解析（压缩二进制序列比对格式、CRAM转BAM、参考序列管理）
- ✅ 实现 ChIPseeker峰注释、DOSE疾病富集、ReactomePA通路分析
- ✅ 实现 AnnotationDbi注释数据库、clusterProfiler富集框架、WGCNA共表达网络
- ✅ 实现 ShortRead 短读序列质量控制（QA统计、adapter修剪、质量修剪、读长过滤、FastQC报告生成）
- ✅ 实现 scater 单细胞质量控制（QC指标计算、细胞/基因过滤、CPM/log-CPM标准化、HVG检测、PCA降维）
- ✅ 实现 MAST 单细胞差异表达分析（Hurdle模型、离散/连续检验、BH-FDR校正、结果汇总）
- ✅ 实现 GenomicFiles 分布式基因组文件处理（BAM/BED/VCF扫描、区间查询、归约、覆盖度计算）
- ✅ 实现 DiffBind ChIP-seq差异结合分析（峰值重叠、共识峰识别、TMM归一化、负二项分布检验）
- ✅ 实现 minfi DNA甲基化分析（NOOB/Illumina/分位数/功能归一化、β/M值计算、DMP/DMR分析）
- ✅ 实现 flowCore 流式细胞术（FCS文件处理、数据变换、荧光补偿、矩形/多边形/椭球/四象限门控）
- ✅ 实现 bsseq 亚硫酸氢盐测序分析（BSmooth平滑、DMR检测、CpG合并、甲基化率计算）
- ✅ 实现 SingleCellExperiment 单细胞核心容器（多assay管理、PCA/t-SNE/UMAP降维、size factors、归一化）
- ✅ 实现 ComplexHeatmap 复杂热图可视化（行/列聚类、颜色映射、热图注释、分组拆分、ASCII热图）
- ✅ 实现 GSVA 基因集变异分析（ssGSEA/zscore/PLAGE评分、富集分析、置换检验、统计摘要）
- ✅ 实现 ChromVAR 染色质变异分析（TF motif富集、GC偏差校正、细胞聚类、变异性分析、偏差图）
- ✅ 实现 DelayedArray 延迟计算数组（懒加载操作、分块处理、行/列聚合、转置、子集操作）
- ✅ 实现 AnnotationFilter 基因注释过滤（染色体筛选、生物类型过滤、链过滤、区域重叠检测、符号模式匹配）
- ✅ 实现 scDblFinder 单细胞双细胞检测（Doublet评分计算、最近邻搜索、双细胞检测、PCA降维、细胞过滤）
- ✅ 实现 ChIPseeker ChIP-seq峰值注释（基因组区域分类、距离TSS分布、注释可视化、统计分析）
- ✅ 实现 DESeq2 差异表达分析（size factors归一化、分散度估计、负二项GLM拟合、Wald检验、LFC收缩）
- ✅ 实现 ChIPseeker峰注释、DOSE疾病富集、ReactomePA通路分析
- ✅ 实现 AnnotationDbi注释数据库、clusterProfiler富集框架、WGCNA共表达网络
- ✅ 实现 Batchelor单细胞批次校正
- ✅ 实现 ChIPseeker ChIP-seq峰值注释（基因组区域分类(启动子/外显子/内含子/UTR/基因间区)、距离TSS分布、BED格式读取、peak2gene关联分析、注释可视化、统计分析）
- ✅ 实现 ChIPseeker ChIP-seq峰值注释（多峰值集重叠分析、Venn图、饼图可视化）
- ✅ 实现 GSVA 基因集变异分析（ssGSEA/zscore/PLAGE评分、富集分析、置换检验、统计摘要、富集图可视化、表型相关性分析、生存分析、分数分布分析）
- ✅ 实现 Taxonomy 分类学分析（NCBI分类数据库解析、分类树操作、谱系查询、共同祖先计算）
- ✅ 实现 GFF GFF3格式解析（基因注释特征提取、属性解析、基因/转录本/CDS/外显子结构分析）
- ✅ 实现 Phylo.Consensus 一致性树构建（多数规则/严格一致性树、分裂分析、支持度计算）
- ✅ 实现 Bio.phenotype 表型微阵列分析（PlateRecord/WellRecord、logistic/Gompertz生长曲线拟合、CSV/JSON解析、控制减法）
- ✅ 实现 Bio.Blast.Applications BLAST命令行工具包装（8种BLAST变体、快速构建器、参数管理、命令构建、验证）
- ✅ 实现 Bio.PDB.PSEA 二级结构预测（CA-CA距离、虚拟键角/二面角、H/E/C分配、三态到八态转换）
- ✅ 实现 Bio.SeqIO.SffIO SFF二进制格式解析（二进制编码/解码、质量修剪、按名称查找）
- [ ] 添加 SIMD 加速支持

## 许可证

Apache-2.0 License
