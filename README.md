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
| **序列处理** | Biopython `Bio.Seq` | 序列对象、**MutableSeq可变序列**、互补、转录、翻译、序列特征 | ✅ |
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
| **Bio.Statistics** | scipy/stats | 描述统计、假设检验、相关性分析（Pearson/Spearman）、置信区间、Z-score、Wilcoxon/Mann-Whitney/Fisher/KS/chi2/t检验、ANOVA、log-rank生存分析、Benjamini-Hochberg/Yekutieli、Bonferroni/Holm校正 | ✅ |
| **Bio.FreqAnalysis** | Biopython `Bio.SeqUtils` | 序列频率分析、k-mer计数、密码子使用频率、GC含量、**Wooton-Federhen局部组成复杂度(LCC)**、序列复杂度 | ✅ |
| **Bio.Align.analysis** | Biopython `Bio.Align.analysis` | dn/ds计算（Nei-Gojobori）、进化距离（Jukes-Cantor/Kimura）、选择压力分析 | ✅ |
| **SeqUtils 高级功能** | Biopython `Bio.SeqUtils` | GC/AT滑动窗口偏斜、ORF预测、序列相似度、Hamming距离、Levenshtein编辑距离 | ✅ |
| **Motifs 高级功能** | Biopython `Bio.motifs` | 序列Logo生成、 per-position信息含量、总信息含量、模体富集分析、Pearson相关性比较 | ✅ |
| **PDB 结构分析高级功能** | Biopython `Bio.PDB` | SASA计算（Shrake-Rupley算法）、Ramachandran质量评估、Kyte-Doolittle疏水性量表、序列属性距离矩阵 | ✅ |
| ✅ | Bio.Align.AlignClusterer | 渐进式多序列比对 (ClustalW算法)，基于UPGMA向导树和profile-profile比对 |
| ✅ | Bio.SeqUtils (RNA) | RNA二级结构预测，Nussinov动态规划算法，支持发夹环/内部环/凸出环/多环识别 |
| ✅ | Bio.PDB.MAalign | 多蛋白质结构比对，Kabsch算法实现，支持迭代结构对齐和保守性分析 |
| ✅ | Bio.SeqFeature (advanced) | CompoundLocation复合位点、LocationParser位点字符串解析、扩展SeqFeature修饰符 |
| ✅ | Bio.PDB.ParsePDBHeader | PDB头部元数据解析 (HEADER/TITLE/COMPOUND/SOURCE/REMARK/AUTH/DBREF) |
| ✅ | Bio.Phylo.PhyloXML | PhyloXML格式解析、序列化、Newick双向转换、分类单元注释 |
| ✅ | Bio.motifs (advanced) | JASPAR PFM格式解析、TRANSFAC格式解析、模体最优比对、KL/JS散度计算、模体聚类 |
| ✅ | Bio.Phylo.NeXML | NeXML格式解析与序列化、OTUs/Trees/Characters数据模型、Newick转换 |
| ✅ | Bio.PDB.Dice + Selection | PDB结构切割（链/残基/原子/模型提取）、B因子过滤、几何选择、结构统计、序列提取 |
| ✅ | Bio.Align.MAF | MAF (Multiple Alignment Format) 多序列比对格式解析、块处理、百分比一致性计算、统计分析、格式转换 |
| ✅ | Bio.Align.Mauve | Mauve 基因组比对格式解析、LCB(共线性块)检测、倒位检测、断点检测、基因组覆盖率、BED导出 |
| ✅ | Bio.Stockholm | Stockholm 格式解析 (Pfam/Rfam比对格式)、二级结构注释、百分比一致性、保守性分析、FASTA转换 |
| ✅ | Bio.PopGen (advanced) | 高级群体遗传学统计: Tajima's D, Fu & Li's D/F, McDonald-Kreitman检验, 等位基因频率谱, 中性分析 |
| ✅ | Bio.SeqUtils.CodonUsage (advanced) | 高级密码子分析: CAI密码子适应指数、RSCU相对同义密码子使用、ENC有效密码子数、GC3偏斜、最优/稀有密码子检测 |
| ✅ | Bio.PDB.Packing | 蛋白质包装密度分析: 局部包装密度、SASA溶剂可及表面积计算、Lee-Richards算法、低包装区域识别、球体点生成 |
| ✅ | Bio.Statistics.qvalue | Storey's q-value FDR方法: π₀估计、q-value计算、自助法π₀、FDR校正、显著性检验 |
| ✅ | Bio.Statistics.IHW | 独立假设加权(IHW): 协变量加权Bonferroni、局部/全局加权、Storey pi0加权、多协变量支持 |
| ✅ | Bioconductor DelayedMatrixStats | DelayedArray统计层: row/col统计(mean/var/sd/median/min/max/sum)、NA处理、子集操作 |
| ✅ | Bioconductor gcrma | GC校正RMA芯片分析: 背景校正(IdealMM/Express)、GC校正、分位数归一化、探针组汇总 |
| ✅ | Bio.Sequencing.Ace | ACE contig组装格式解析: reads/contigs解析、共有序列、覆盖度分析、GC含量计算 |
| ✅ | Bio.SeqUtils.Proteomics | 蛋白质组学工具: 8种酶切(胰酶/糜酶/胃酶等)、肽段质量计算、同位素分布、b/y碎片离子 |
| ✅ | Bio.PDB.FragmentMapper | PDB片段映射: DSSP二级结构分类、片段分配/合并/过滤、覆盖率分析 |
| ✅ | Bio.Graphics.GenomeDiagram | 基因组图可视化: 数据模型、track/feature管理、SVG生成、样式控制 |
| ✅ | Bio.PDB.internal_coords | 蛋白质内部坐标表示: 键长/键角/扭矩角(phi/psi/omega/chi)、二面角计算、内部坐标到笛卡儿坐标转换、扩展链构建、旋转异构体库、Ramachandran区域 |
| ✅ | Bio.GA | 遗传算法序列优化: 种群初始化、适应度函数、锦标赛/轮盘赌选择、单点/两点交叉、变异、精英保留、世代统计、进化收敛 |
| ✅ | Bio.Graphics.Chromosome | 染色体可视化: 染色体/特征/区域/带数据模型、SVG线性/圆形染色体渲染、G带染色模式、带颜色映射、人类核型图、细菌染色体图 |
| ✅ | Bioconductor Chain文件/liftOver | UCSC Chain格式解析、基因组坐标liftOver转换、链段查找、染色体间坐标映射 |
| ✅ | Bioconductor Biostrings matchPDict | 字典模式匹配(matchPDict/vmatchPattern)、多序列模式计数、错配容忍、最佳匹配查找 |
| ✅ | Bioconductor GenomicRanges gaps/reduce/disjoin | gaps检测、reduce合并、disjoin拆分、setdiff/交集/并集集合运算、coverage计算、promoters提取 |
| ✅ | NGS质量修剪与接头去除 | 质量修剪(滑动窗口)、接头去除、poly-A修剪、长度/GC含量过滤、批量修剪、Fastq解析与序列化 |
| ✅ | Bioconductor Rsubread/featureCounts | 基因特征read计数: 链特异性(Stranded/Reversed/Unstranded)、重叠长度阈值、最小比对质量过滤、多比对read处理(计数/分数计数)、配对末端fragment计数、CPM归一化、多样本矩阵 |
| ✅ | Bioconductor DRIMSeq | 差异转录本使用分析: Dirichlet-multinomial模型、Wald检验、比例计算、counts过滤、BH-FDR校正、显著基因提取、收敛诊断 |
| ✅ | Bioconductor RaggedExperiment | 参差突变数据结构: 突变记录(Missense/Nonsense/Frame_Shift/Splice_Site等)、基因×样本稀疏矩阵、TMB计算(每Mb)、按基因/样本/突变类型过滤、突变计数矩阵、LoF识别 |
| ✅ | Bioconductor HTSFilter | RNA-seq count过滤: CPM归一化、按组最小样本数阈值、keep mask生成、文库大小估计、过滤矩阵导出 |
| ✅ | Bioconductor baySeq | 贝叶斯差异表达分析: 负二项模型、分散度估计、Gamma先验、对数似然比、后验概率、MAP表达估计、差异基因判定 |
| ✅ | Bioconductor CellChat | 细胞间通讯分析: 配体-受体互作对数据库、置换检验、互作评分、细胞类型聚合、显著性检验、FDR校正 |
| ✅ | Bio.File | 智能文件处理: 自动压缩格式检测(gzip/bzip2)、透明压缩读写、文件操作接口 |
| ✅ | Bio.SeqUtils.MolWt | 分子量计算: DNA/RNA/蛋白质分子量、消光系数、吸光度、等电点(pI) |
| ✅ | Bio.Align.Reduced | 简化氨基酸字母表: RAD/Dayhoff/CHARM/SDM12简化字母表、序列比较、简化一致性 |
| ✅ | Hmisc | Bioconductor Hmisc统计工具包: 相关性分析(Pearson/Spearman)、变量聚类、描述性统计、Somers' d统计、缺失值插补 |
| ✅ | rstatix | Bioconductor rstatix tidy统计检验: T检验、Wilcoxon检验、ANOVA、Kruskal-Wallis、Friedman检验、相关性检验、BH-FDR校正、Bonferroni校正 |
| ✅ | VennDiagram | Bioconductor VennDiagram集合分析: 集合运算、Venn区域计算、重叠统计、相似度指标(Jaccard/Dice/Overlap coefficient) |
| ✅ | GENIE3 | Bioconductor GENIE3基因调控网络推断: 回归树特征重要性、方差缩减、加权邻接矩阵、对称化网络 |
| ✅ | decoupleR | Bioconductor decoupleR功能活性推断: WSum/WMean/Norm/ULM/MLM方法、先验知识网络(PKN)、调控子活性评分 |
| ✅ | BayesSpace | Bioconductor BayesSpace空间转录组聚类: t分布混合模型、马尔可夫随机场(MRF)先验、EM算法、六边形/方形网格邻居 |
| ✅ | muscat | Bioconductor muscat单细胞差异状态分析: 伪批量聚合(Sum/Mean/Median)、EdgeR/DESeq2/Limma DS检验、BH-FDR校正、样本QC指标 |
| ✅ | infercnv | Bioconductor infercnv单细胞拷贝数变异推断: 染色体位置排序基因、参考细胞比较、log2FC有界计算、金字塔权重基因组平滑、每细胞中位数中心化+噪声过滤、CNV分数+肿瘤细胞预测 |
| ✅ | SCENIC | Bioconductor SCENIC单细胞调控网络推断与聚类: TF-target共表达模块(GENIE3风格)、Regulon构建(权重剪枝/cisTarget motif排名剪枝)、AUCell活性评分(recovery curve AUC)、二值化阈值(MeanStd/KMeans2/Median)、细胞状态聚类+主控调控因子识别 |

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
| **序列模体** | Biopython `Bio.motifs` | PWM、MEME格式、JASPAR PFM、TRANSFAC、模体搜索、信息含量、序列Logo、模体富集、Pearson相关性、最优比对、KL/JS散度、模体聚类 | ✅ |
| **群体遗传学** | Biopython `Bio.PopGen` | 等位基因频率、FST、哈迪-温伯格检验 | ✅ |
| **edgeR** | Bioconductor edgeR | 差异表达分析、DGEList、精确检验、GLM拟合 | ✅ |
| **limma** | Bioconductor limma | 差异表达分析、线性模型拟合、经验贝叶斯、voom变换、RPKM/CPM/quantile归一化、ComBat/removeBatchEffect批次校正、treat严格检验 | ✅ |
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
| **NeXML格式** | Biopython `Bio.Phylo.NeXML` | NeXML格式解析与序列化、OTUs/Trees/Characters数据模型、Newick转换 | ✅ |
| **蛋白质参数分析** | Biopython `Bio.SeqUtils.ProtParam` | 分子量、不稳定指数、GRAVY评分、等电点、信号肽预测 | ✅ |
| **变异分析** | Biopython `Bio.Variation` | SNP分析、突变检测、氨基酸替换分析、BLOSUM62/Grantham矩阵 | ✅ |
| **DSSP** | Biopython `Bio.PDB.DSSP` | 二级结构预测、溶剂可及表面积(SASA/Shrake-Rupley)、Ramachandran图、结构质量评估 | ✅ |
| **PDB Dice + Selection** | Biopython `Bio.PDB.Dice` | PDB结构切割（链/残基/原子/模型提取）、B因子/occupancy过滤、几何选择、结构统计、序列提取 | ✅ |
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
| **UniProtIO** | Biopython `Bio.SeqIO.UniprotIO` | UniProt XML格式解析、蛋白质条目提取、基因名、物种、序列、功能注释、数据库交叉引用 | ✅ |
| **chem_utils** | Biopython `Bio.PDB.chem_utils` | 化学计算工具：范德华半径、共价半径、键长、键角、二面角、经验式、分子式量、氢键长度 | ✅ |
| **mmCIF** | Biopython `Bio.PDB.MMCIFParser` | mmCIF格式解析、数据块、类别、原子位点提取 | ✅ |
| **Nexus** | Biopython `Bio.Nexus` | NEXUS格式解析、数据矩阵、系统发育树、距离矩阵 | ✅ |
| **EMBOSS** | EMBOSS suite | GC偏斜、AT偏斜、分子量、Tm值、ORF查找、距离计算、蛋白质参数 | ✅ |
| **ChIPseeker** | Bioconductor ChIPseeker | ChIP-seq峰注释、基因距离计算、注释分类（启动子/外显子/内含子/UTR/基因间区）、BED格式读取、peak2gene关联分析、多峰值集重叠分析(peakOverlap)、Venn图可视化、饼图可视化、结果汇总与可视化 | ✅ |
| **DOSE** | Bioconductor DOSE | 疾病本体富集分析、超几何检验、多重检验校正（adjusted p-value）、富集结果过滤与汇总 | ✅ |
| **ReactomePA** | Bioconductor ReactomePA | Reactome通路富集分析、通路按顶层术语分组、富集结果可视化、通路注释查询 | ✅ |
| **AnnotationDbi** | Bioconductor AnnotationDbi | 通用注释数据库接口、基因信息存储、ID映射、GO/KEGG注释管理 | ✅ |
| **clusterProfiler** | Bioconductor clusterProfiler | 功能富集分析统一框架、超几何检验、多重检验校正、富集结果可视化 | ✅ |
| **WGCNA** | Bioconductor WGCNA | 加权基因共表达网络分析、邻接矩阵构建、TOM相似度、模块检测 | ✅ |
| **BiocNeighbors** | Bioconductor BiocNeighbors | 最近邻搜索（KMKNN/Annoy/BruteForce）、欧几里得/曼哈顿/余弦距离、索引构建与查询 | ✅ |
| **Biobase** | Bioconductor Biobase | ExpressionSet数据结构、AnnotatedDataFrame、数据归一化、log2转换 | ✅ |
| **GEOquery** | Bioconductor GEOquery | GEO数据库数据获取、Series Matrix解析、SOFT格式解析、ExpressionSet转换 | ✅ |
| **tximport** | Bioconductor tximport | 转录本量化数据导入、Salmon quant.sf解析、基因级别汇总、低表达过滤 | ✅ |
| **AnnotationHub** | Bioconductor AnnotationHub | 中心化注释资源访问、资源搜索（类型/提供者/基因组/关键词）、资源管理与缓存 | ✅ |
| **GenomicFeatures** | Bioconductor GenomicFeatures | 基因组注释功能、Gene/Transcript/Exon数据结构、GTF解析、区域查询 | ✅ |
| **graph** | Bioconductor graph | 图数据结构、有向/无向图、最短路径算法、连通分量检测、DOT格式输出 | ✅ |
| **SpatialExperiment** | Bioconductor SpatialExperiment | 空间转录组学数据结构、空间坐标管理、图像数据存储、空间范围查询、spots过滤 | ✅ |
| **phyloseq** | Bioconductor phyloseq | 微生物组分析、OTU丰度计算、分类学过滤、相对丰度计算、稀疏化处理 | ✅ |
| **microbiome** | Bioconductor microbiome | 微生物组分析、Alpha多样性（Shannon/Simpson/Chao1/ACE/Fisher/Pielou）、Beta多样性（Bray-Curtis/Jaccard/JSD/weighted/unweighted UniFrac）、PCoA主坐标分析、差异丰度分析（Welch t检验/Wilcoxon秩和检验/BH校正） | ✅ |
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
| **SingleR** | Bioconductor SingleR | 细胞类型注释、Spearman/Pearson相关性、参考图谱匹配、精细调优(Fine-tuning)、Delta score置信度评估 | ✅ |
| **Cyclone** | Bioconductor cyclone | 细胞周期评分、基因对(Gene pairs)比较、G1/S/G2/M期相预测、相别分布统计、平均得分分析 | ✅ |
| **dorothea** | Bioconductor dorothea | 转录因子活性预测、Regulon分析、VIPER算法、置换检验、Z-score评估、Top TF筛选 | ✅ |
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
| **Rtsne** | Bioconductor Rtsne | t-SNE降维算法、成对距离计算、条件概率估计（perplexity优化）、联合概率矩阵构建、梯度下降优化（动量/early exaggeration）、Barnes-Hut近似 | ✅ |
| **uwot** | Bioconductor uwot | UMAP降维算法、k近邻搜索、模糊单纯集构建、局部模糊集并集、低维嵌入优化（SGD/负采样）、min_dist/spread参数控制 | ✅ |
| **tradeSeq** | Bioconductor tradeSeq | 轨迹差异表达分析、GAM（广义可加模型）拟合、样条基函数、差异表达检验、BH-FDR校正 | ✅ |
| **PROGENy** | Bioconductor PROGENy | 通路活性推断、L2正则化线性回归（Ridge回归）、通路基因集权重矩阵、样本通路活性计算 | ✅ |
| **AUCell** | Bioconductor AUCell | 单细胞基因集评分、AUC（曲线下面积）计算、基因排序、min-max归一化、细胞/基因集评分查询 | ✅ |
| **ggtree** | Bioconductor ggtree | 系统发育树可视化布局算法、矩形布局（phylogram）、放射状布局、无根布局、节点坐标映射、边缘/标签数据生成 | ✅ |
| ✅ | BiocSingular | SVD奇异值分解，支持Exact/IRLBA/Randomized三种算法，用于单细胞降维 |
| ✅ | BiocNeighbors | KMKNN和Annoy最近邻搜索，支持欧几里得/曼哈顿/余弦距离 |
| ✅ | mixOmics | 多组学整合方法，包括PLS回归、稀疏PLS (sPLS)、DIABLO多块整合 |
| **MAF格式解析** | Biopython `Bio.Align` | MAF多序列比对格式解析、块操作、百分比一致性、统计分析、选择/过滤/写回 | ✅ |
| **UCSC Chain文件/liftOver** | Bioconductor rtracklayer | Chain格式解析、基因组坐标liftOver转换、链段查找、染色体间坐标映射、位置/区间转换 | ✅ |
| **Biostrings matchPDict** | Bioconductor Biostrings | 字典模式匹配(matchPDict/vmatchPattern)、多序列模式计数(vcountPattern)、错配容忍、最佳匹配查找 | ✅ |
| **GenomicRanges gaps/reduce/disjoin** | Bioconductor GenomicRanges | gaps检测、reduce合并、disjoin拆分、setdiff/交集/并集集合运算、coverage计算、promoters提取、trim | ✅ |
| **NGS质量修剪与接头去除** | Bioconductor ShortRead | 质量修剪(滑动窗口)、接头去除、poly-A修剪、长度/GC含量过滤、批量修剪、Fastq解析与序列化、统计计算 | ✅ |
| **Mauve基因组比对** | Biopython `Bio.Align` | Mauve基因组比对格式、LCB检测、倒位/断点检测、覆盖率分析、BED导出、基因组重排率 | ✅ |
| **Stockholm格式** | Biopython `Bio.Stockholm` | Stockholm/Pfam格式解析、二级结构注释、百分比一致性、保守性、FASTA/Stockholm互转 | ✅ |
| **高级群体遗传学** | Biopython `Bio.PopGen` | Tajima's D中性检验、Fu & Li's D/F、McDonald-Kreitman检验、等位基因频率谱、综合中性分析 | ✅ |
| **高级密码子分析** | Biopython `Bio.SeqUtils.CodonUsage` | CAI密码子适应指数、RSCU相对同义密码子使用、ENC有效密码子数、GC3偏斜、最优/稀有密码子、物种特异性参考表 | ✅ |
| **蛋白质包装密度** | Biopython `Bio.PDB.Packing` | 局部包装密度、SASA计算(Lee-Richards)、范德华半径、低包装区域识别、球体几何分析 | ✅ |
| **Storey q-value** | Bioconductor qvalue | π₀估计、Storey q-value计算、自助法π₀、FDR校正、显著性检验 | ✅ |
| **独立假设加权(IHW)** | Bioconductor IHW | 协变量加权Bonferroni、局部/全局加权、Storey pi0加权、多协变量支持、迭代权重优化 | ✅ |
| **DelayedMatrixStats** | Bioconductor DelayedMatrixStats | DelayedArray统计层、row/col统计(mean/var/sd/median/min/max/sum)、NA处理、子集操作 | ✅ |
| **GC-RMA芯片分析** | Bioconductor gcrma | GC校正RMA、背景校正(IdealMM/Express)、GC查找表、分位数归一化、探针组汇总 | ✅ |
| **ACE contig格式** | Biopython `Bio.Sequencing.Ace` | ACE组装格式解析、reads/contigs提取、共有序列生成、覆盖度分析、GC含量计算、格式化输出 | ✅ |
| **蛋白质组学分析** | Biopython `Bio.SeqUtils.Proteomics` | 8种蛋白酶切(胰酶/糜酶/胃酶/LysC/ArgC/CNBr/GluC/AspN)、单同位素/平均质量计算、同位素分布、b/y碎片离子 | ✅ |
| **PDB片段映射** | Biopython `Bio.PDB.FragmentMapper` | DSSP二级结构分类(H/E/S/L/T/C)、片段分配/合并/过滤、覆盖率分析、边界检测、类型转换 | ✅ |
| **基因组图可视化** | Biopython `Bio.Graphics.GenomeDiagram` | 数据模型(Track/Feature/Diagram)、多种形状(Rectangle/Arrow/Diamond)、SVG生成、样式控制、自动标注、重叠检测 | ✅ |

| **Gviz** | Bioconductor Gviz | 基因组可视化轨道系统、多种轨道类型(AnnotationTrack/GeneRegionTrack/DataTrack/IdeogramTrack/GenomeAxisTrack/SequenceTrack)、ASCII渲染、特征查询与区域可视化 | ✅ |
| **SeqLocation 序列位点** | Biopython `Bio.SeqFeature` | 序列位置类型系统(ExactPosition/BeforePosition/AfterPosition/OneOfPosition/WithinPosition)、SimpleLocation/CompoundLocation、链方向、位置转换(start/end/strand)、位置字符串解析与格式化 | ✅ |
| **BioReference 文献引用** | Biopython `Bio.Reference` / `Bio.Medline` | 文献引用管理(title/authors/journal/year/pubmed_id/doi)、引用位置(locations)、引用类型(journal article/etc)、格式化输出(APA/自定义)、Medline记录解析 | ✅ |
| **ProtDao 蛋白质无序预测** | Biopython `Bio.SeqUtils.ProtDao` | IUPred-like蛋白质无序区域预测、氨基酸无序度评分、能量评分、无序区域检测与分类、阈值可调、ASCII可视化、序列复杂度分析 | ✅ |
| **PairwiseAligner 统一比对器** | Biopython `Bio.Align.PairwiseAligner` | 全局(Needleman-Wunsch)/局部(Smith-Waterman)双序列统一比对、Gotoh仿射空位罚分(open+extend)、独立target/query gap参数、DNA打分/蛋白质BLOSUM62替换矩阵、一致/匹配线、identity%/gaps计数、target_start/end query_start/end定位 | ✅ |
| **NaiveBayes 序列分类器** | Biopython `Bio.NaiveBayes` | 基于k-mer频率的朴素贝叶斯序列分类、Laplace平滑(alpha参数)、对数后验概率、softmax归一化概率、top-K预测、类别先验、准确率计算、AT/GC富集分类示例 | ✅ |
| **Markov 马尔可夫链建模** | Biopython `Bio.Markov` | 1/2/3阶马尔可夫链训练、状态字母表(DNA/蛋白质)、转移概率矩阵+伪计数平滑、序列对数概率打分、单位置log-prob数组、加权采样生成序列、CpG岛log-odds检测、稳态分布迭代求解 | ✅ |

项目致力于打造一个完整、高效的生物信息学工具库，覆盖从基础序列处理到高级序列组装的全流程。

## 架构设计

### 项目结构

```
IvanAXu/BioSeqs/
├── moon.mod                    # 模块配置
├── src/                        # 源代码
│   ├── seq.mbt                 # Seq / MutableSeq 序列对象
│   ├── seq_record.mbt          # SeqRecord 带注释的序列记录
│   ├── seqfeature.mbt          # SeqFeature 序列特征
│   ├── seqfeature_advanced.mbt  # Bio.SeqFeature CompoundLocation与LocationParser
│   ├── seqio.mbt               # 统一序列 I/O 接口
│   ├── fasta_io.mbt            # FASTA 格式解析
│   ├── fastq_io.mbt            # FASTQ 格式解析
│   ├── genbank_io.mbt          # GenBank 格式解析
│   ├── align.mbt               # MultipleSeqAlignment 多序列比对
│   ├── alignio.mbt             # 比对文件 I/O
│   ├── clustal_io.mbt          # Clustal 格式
│   ├── phylip_io.mbt           # PHYLIP 格式
│   ├── alignment.mbt           # DNA/RNA/Protein 类型及比对算法
│   ├── bioc_singular.mbt       # BiocSingular SVD奇异值分解 (Exact/IRLBA/Randomized)
│   ├── application.mbt          # Bio.Application 命令行工具包装
│   ├── align_io.mbt            # AlignIO 比对格式解析 (ClustalW/FASTA/Stockholm)
│   ├── align_applications.mbt   # Align.Applications 比对工具包装 (ClustalW/Clustal Omega/Muscle/MAFFT)
│   ├── phylo.mbt               # 系统发育树 (Clade/Tree)
│   ├── tree_io.mbt             # 进化树格式解析 (Newick、NHX格式解析、树操作)
│   ├── blast.mbt               # BLAST结果解析 (tabular/xml格式、HSP、Hit、Record、过滤)
│   ├── searchio.mbt            # SearchIO 统一搜索结果模型 (HSPFragment、HSP、Hit、QueryResult、HMMER3/BLAT解析)
│   ├── search_io.mbt           # SearchIO 统一搜索结果模型 (HMMER3解析、BLAT PSL解析、BLAST转换)
│   ├── subsmat.mbt             # 替换矩阵 (BLOSUM62/45、PAM250/30、矩阵解析、分数查询)
│   ├── pdb.mbt                 # PDB 数据类型
│   ├── pdb_io.mbt              # PDB 文件 I/O
│   ├── pdb_header.mbt            # Bio.PDB.ParsePDBHeader PDB头部元数据解析
│   ├── pdb_analysis.mbt         # Bio.PDB.StructureAnalysis 结构分析 (二面角、SASA、Ramachandran图)
│   ├── pdb_list.mbt             # Bio.PDB.PDBList PDB结构下载管理
│   ├── pdb_dice.mbt             # Bio.PDB.Dice + Selection PDB结构切割 (链/残基/原子/模型提取、B因子过滤、几何选择、结构统计)
│   ├── svd_superimposer.mbt    # SVDSuperimposer SVD蛋白质结构叠合 (旋转矩阵、平移向量、RMSD计算)
│   ├── structure_alignment.mbt  # Bio.PDB.StructureAlignment 多蛋白质结构比对
│   ├── ma_align.mbt             # Bio.PDB.MAalign 多蛋白质结构比对 (Kabsch算法)
│   ├── sequtils.mbt            # 序列工具函数
│   ├── seq_utils.mbt           # 序列工具函数 (分子量、GC含量、Tm值、氨基酸转换、GC/AT滑动窗口偏斜、序列相似度、Hamming距离、Levenshtein编辑距离)
│   ├── seq_complexity.mbt      # 序列复杂度分析 (Shannon熵、语言学复杂度、**Wooton-Federhen LCC**、DUST分数、GC偏斜、CGR、序列签名、k-mer相似度)
│   ├── rna_structure.mbt       # Bio.SeqUtils RNA二级结构预测 (Nussinov算法)
│   ├── complement.mbt          # 互补碱基查找表
│   ├── codon_table.mbt         # 密码子翻译表
│   ├── codon_usage.mbt         # CodonUsage 密码子使用分析 (CAI、ENC、RSCU、GC3、CBI、Fop、最优密码子检测)
│   ├── sam.mbt                 # SAM 格式解析
│   ├── bam.mbt                 # BAM 格式解析
│   ├── bgzf.mbt                # BGZF 解压缩 (支持读取压缩的 BAM 文件)
│   ├── cram_wbtest.mbt         # CRAM 格式解析与白盒测试 (压缩二进制序列比对格式、内部函数测试、格式验证)
│   ├── vcf.mbt                 # VCF 格式解析
│   ├── variant_annotation.mbt  # VariantAnnotation 变异注释 (变异类型检测、定位、编码效应预测、变异汇总)
│   ├── variant_filtering.mbt    # Bioconductor VariantFiltering 变异过滤与遗传模式分析
│   ├── biostrings.mbt          # Biostrings 序列分析 (IUPAC、k-mer频率、RSCU、复杂度)
│   ├── biobase.mbt             # Bioconductor 基础函数 (常用统计、数据转换、工具函数)
│   ├── bioc_generics.mbt        # Bioconductor 通用函数 (NA处理、排序、集合运算)
│   ├── bioc_parallel.mbt        # Bioconductor 并行计算框架 (任务分块、并行求和)
│   ├── genomic_ranges.mbt      # GenomicRanges 基因组区间操作 (GRanges、IRanges)
│   ├── genomic_ranges_advanced.mbt # GenomicRanges tile/slidingWindows/区间运算
│   ├── iranges.mbt             # IRanges 整数区间操作 (集合运算、重叠检测)
│   ├── genomic_alignments.mbt  # GenomicAlignments 基因组比对分析 (GAlignments、coverage、summarizeOverlaps、pileup)
│   ├── txdb.mbt                # TxDb 转录本数据库 (GTF解析、基因/转录本/外显子/CDS提取、UTR/内含子计算)
│   ├── tximport.mbt             # Bioconductor tximport 转录本量化数据导入
│   ├── rtracklayer.mbt         # rtracklayer 基因组轨道格式 (BED/WIG/BEDGraph/GFF解析与写入)
│   ├── rhdf5.mbt                # Bioconductor rhdf5 HDF5文件格式支持
│   ├── deseq2.mbt              # DESeq2 差异表达分析 (size factors归一化、分散度估计、负二项GLM拟合、Wald检验、LFC收缩)
│   ├── deseq2_advanced.mbt         # DESeq2 VST方差稳定化变换、PCA可视化
│   ├── edger.mbt               # edgeR 差异表达分析 (DGEList、精确检验、GLM拟合)
│   ├── edger_advanced.mbt           # edgeR准似然F检验、camera/roast基因集检验
│   ├── limma.mbt               # limma 差异表达、归一化、批次校正 (线性模型、经验贝叶斯、voom、RPKM/CPM/quantile、ComBat)
│   ├── matrix.mbt               # Bioconductor Matrix 稀疏矩阵操作 (CSC/CSR格式、矩阵运算)
│   ├── bioc_neighbors.mbt      # BiocNeighbors 最近邻搜索 (KMKNN/Annoy)
│   ├── summarized_experiment.mbt # SummarizedExperiment 多维基因组数据容器
│   ├── dplyr.mbt               # dplyr 数据操作 (DataFrame、filter、select、mutate、arrange、group_by、summarize、join)
│   ├── smith_waterman.mbt      # Smith-Waterman 局部序列比对 (动态规划、自定义打分、回溯矩阵)
│   ├── needleman_wunsch.mbt    # Needleman-Wunsch 全局序列比对 (动态规划、自定义打分、回溯矩阵)
│   ├── bloom_filter.mbt        # Bloom Filter & k-mer 计数 (哈希表、位图、误判率估算、近似去重)
│   ├── bwt_fm.mbt              # BWT + FM-index (后缀数组、BWT变换、回溯搜索、精确模式匹配)
│   ├── de_bruijn.mbt           # De Bruijn Graph (k-mer节点、欧拉路径、序列组装、图简化)
│   ├── suffix_array_tree.mbt   # Suffix Array & Suffix Tree (前缀倍增、LCP数组、模式匹配、最长重复子串)
│   ├── olc.mbt                 # Overlap-Layout-Consensus (重叠检测、哈密顿路径、一致性序列生成)
│   ├── paml.mbt                 # Bio.PAML 分子进化分析 (dN/dS、Jukes-Cantor校正)
│   ├── hmm.mbt                 # Hidden Markov Model (前向/后向算法、维特比算法、Baum-Welch训练、基因预测)
│   ├── kmeans.mbt              # K-means Clustering (距离计算、K-means++初始化、聚类、轮廓系数评估)
│   ├── kmer.mbt                 # Bio.Kmer k-mer计数与频率分析
│   ├── cluster.mbt             # 序列聚类分析 (距离矩阵、层次聚类、轮廓系数)
│   ├── motifs.mbt              # 序列模体识别 (PWM、MEME格式、模体搜索、信息含量、序列Logo生成、模体富集、Pearson相关性)
│   ├── motifs_advanced.mbt     # Bio.motifs高级功能 (JASPAR PFM、TRANSFAC格式、模体最优比对、KL/JS散度、模体聚类)
│   ├── popgen.mbt              # 群体遗传学 (等位基因频率、FST、哈迪-温伯格检验)
│   ├── residue_depth.mbt        # 残基深度分析 (晶体学B因子、溶剂可及性、埋藏度分析)
│   ├── restriction.mbt         # 限制性内切酶分析 (酶切位点查找、片段分析)
│   ├── protparam.mbt           # ProtParam 蛋白质参数分析 (不稳定指数、等电点、信号肽预测、二级结构倾向)
│   ├── prosite.mbt              # Bio.Prosite 蛋白质模体数据库搜索
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
│   ├── isoform_switch_analyze_r.mbt # IsoformSwitchAnalyzeR 转录本异构体切换分析
│   ├── multi_assay_experiment.mbt # MultiAssayExperiment 多组学数据协调 (实验协调、样本映射、跨实验子集)
│   ├── tree_construction.mbt   # TreeConstruction 系统发育树构建 (UPGMA/WPGMA/NJ算法、距离计算)
│   ├── neighbor_search.mbt     # NeighborSearch KD树近邻搜索 (空间搜索、半径搜索、最近邻)
│   ├── swissprot.mbt           # SwissProt 蛋白数据库解析 (UniProt/Swiss-Prot记录、特征、参考文献)
│   ├── uniprot_io.mbt          # UniProt XML格式解析 (蛋白质条目、基因名、物种、序列、功能注释)
│   ├── chem_utils.mbt          # 化学计算工具 (范德华半径、共价半径、键长、键角、二面角、分子式量)
│   ├── jaspar.mbt              # JASPAR PFM格式解析 (模体矩阵、PWM转换、共有序列、序列扫描)
│   ├── mmcif.mbt               # mmCIF格式解析 (Bio.PDB.MMCIFParser、数据块、类别、原子位点)
│   ├── nexus.mbt               # Nexus格式解析 (Bio.Nexus、数据矩阵、系统发育树、距离矩阵)
│   ├── emboss.mbt              # EMBOSS工具接口 (GC偏斜、AT偏斜、分子量、Tm值、ORF查找、距离计算)
│   ├── chipseeker.mbt          # ChIPseeker ChIP-seq峰注释分析 (峰-基因距离计算、注释分类(启动子/外显子/内含子/UTR/基因间区)、BED格式读取、peak2gene关联分析、结果汇总与可视化)
│   ├── enrichplot.mbt           # enrichplot 富集分析结果可视化 (dotplot/barplot/heatmap)
│   ├── ensembldb.mbt            # ensembldb Ensembl注释数据库接口
│   ├── expasy.mbt               # ExPASy 蛋白质分析工具接口 (Swiss-Prot解析、蛋白质参数计算)
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
│   ├── seq_quality_trim.mbt    # NGS质量修剪与接头去除 (质量修剪、接头去除、poly-A修剪、长度/GC过滤、批量修剪)
│   ├── scater.mbt              # scater 单细胞质量控制 (QC指标计算、细胞/基因过滤、标准化、HVG检测、PCA)
│   ├── mast.mbt                # MAST 单细胞差异表达分析 (Hurdle模型、离散/连续检验、BH-FDR校正)
│   ├── genomic_files.mbt       # GenomicFiles 分布式基因组文件处理 (BAM/BED/VCF扫描、区间查询、归约、覆盖度)
│   ├── diffbind.mbt            # DiffBind ChIP-seq差异结合分析 (峰值重叠、共识峰、TMM归一化、NB检验)
│   ├── minfi.mbt               # minfi DNA甲基化分析 (NOOB/Illumina/分位数/功能归一化、β/M值、DMP/DMR分析)
│   ├── flow_core.mbt           # flowCore 流式细胞术FCS处理 (数据变换、荧光补偿、矩形/多边形/椭球/四象限门控)
│   ├── bsseq.mbt               # bsseq 亚硫酸氢盐测序分析 (BSmooth平滑、DMR检测、CpG合并、甲基化率计算)
│   ├── single_cell_experiment.mbt  # SingleCellExperiment 单细胞核心容器 (多assay、PCA/tSNE/UMAP降维、size factors)
│   ├── spatial_experiment.mbt   # SpatialExperiment 空间转录组学数据结构
│   ├── complex_heatmap.mbt      # ComplexHeatmap 复杂热图可视化 (行/列聚类、颜色映射、热图注释、分组拆分)
│   ├── gsva.mbt                 # GSVA 基因集变异分析 (ssGSEA/zscore/PLAGE评分、富集分析、置换检验、富集图可视化、表型相关性分析、生存分析、分数分布分析)
│   ├── chromvar.mbt             # ChromVAR 染色质变异分析 (TF motif富集、GC偏差校正、细胞聚类、变异性分析)
│   ├── delayed_array.mbt        # DelayedArray 延迟计算数组 (懒加载操作、分块处理、行/列聚合、子集操作)
│   ├── annotation_filter.mbt    # AnnotationFilter 基因注释过滤 (染色体筛选、生物类型过滤、区域重叠检测、符号模式匹配)
│   ├── sc_dbl_finder.mbt        # scDblFinder 单细胞双细胞检测 (Doublet评分计算、最近邻搜索、PCA降维、细胞过滤)
│   ├── batchelor.mbt            # Batchelor 单细胞批次校正 (rescaleBatches、mutual nearest neighbor、fastMNN、批次混合评分)
│   ├── seurat.mbt               # Seurat 单细胞数据分析核心 (标准化、高可变基因、PCA、聚类、UMAP、差异表达、跨样本整合)
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
│   ├── single_r.mbt             # SingleR 细胞类型注释 (参考图谱、Spearman/Pearson相关性、精细调优)
│   ├── cyclone.mbt              # Cyclone 细胞周期评分 (基因对比较、G1/S/G2/M期相预测)
│   ├── dorothea.mbt             # dorothea 转录因子活性预测 (Regulon、VIPER、置换检验)
│   ├── gff.mbt                  # GFF GFF3格式解析 (GFFFeature/GFFRecord、属性解析、特征提取)
│   ├── graphics.mbt             # Bio.Graphics 生物信息学可视化 (序列Logo、比对可视化)
│   ├── phylo_consensus.mbt      # Phylo.Consensus 一致性树构建 (ConsensusNode/Split、多数规则树、支持度计算)
│   ├── phylo_xml.mbt            # Bio.Phylo.PhyloXML XML系统发育树格式
│   ├── phylo_nexml.mbt          # Bio.Phylo.NeXML NeXML格式解析与序列化 (OTUs/Trees/Characters数据模型、Newick转换)
│   ├── phylo_metrics.mbt        # 系统发育树度量 (Colless平衡指数、Robinson-Foulds距离、距离矩阵)
│   ├── phyloseq.mbt            # phyloseq 微生物组分析 (OTU丰度、分类学过滤)
│   ├── microbiome.mbt          # microbiome 微生物组分析 (Alpha/Beta多样性、PCoA、差异丰度分析)
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
│   ├── qcp_superimposer.mbt     # QCP叠加 (四元数旋转、结构比对、RMSD计算、最优叠加)
│   ├── psea.mbt                 # Bio.PDB.PSEA 二级结构预测 (PseaAtom/PseaResult、CA-CA距离、虚拟二面角、H/E/C分配、三态到八态转换)
│   ├── sff_io.mbt               # Bio.SeqIO.SffIO SFF二进制格式解析 (SffHeader/SffRead/SffFile、二进制编码/解码、质量修剪)
│   ├── seq_complexity.mbt       # 序列复杂度与组成分析 (Shannon熵、GC偏斜、混沌游戏表示)
│   ├── csaw.mbt                 # csaw ChIP-seq窗口差异分析 (滑动窗口计数、TMM归一化、窗口过滤、负二项GLM检验、差异区域检测)
│   ├── slingshot.mbt            # slingshot 单细胞轨迹推断 (MST构建、主曲线拟合、拟时间计算、分支检测)
│   ├── scnorm.mbt               # SCnorm 单细胞RNA-seq归一化 (分位数回归、深度依赖偏差校正、基因特异性归一化)
│   ├── edaseq.mbt               # EDASeq RNA-seq探索性分析 (GC含量归一化、基因长度校正Loess、样本间归一化、RPKM计算)
│   ├── searchio.mbt             # Bio.SearchIO 统一搜索结果模型 (BLAST/HMMER解析、QueryResult/Hit/HSP层次结构、E-value过滤)
│   ├── pdb_vectors.mbt          # Bio.PDB.vectors 3D向量与旋转矩阵 (Vector3/RotationMatrix3、叉积、Kabsch叠合、二面角)
│   ├── circ_seq.mbt             # Bio.SeqUtils.CircSeq 环状DNA序列操作 (circ_subseq、circ_rotate、circ_digest、限制性酶切分析)
│   ├── align_abstract.mbt       # Bio.Align.AlignAbstract 抽象比对类型 (AbstractAlignment、一致性序列、Shannon熵、简约信息位点)
│   ├── align_cluster.mbt        # Bio.Align.AlignClusterer 渐进式多序列比对 (ClustalW算法)
│   ├── maftools.mbt             # maftools 癌症基因组学MAF解析与突变分析 (MAFMutation、MAFData、MutationSpectrum、TMB计算、共现分析)
│   ├── cnvkit.mbt               # CNVkit 拷贝数变异检测与CBS分割 (CNVProbe、CNVSegment、CBS分割算法、拷贝数状态判定)
│   ├── destiny.mbt              # destiny 单细胞扩散映射降维 (CellData、距离矩阵、高斯核、特征分解、扩散分量)
│   ├── rtsne.mbt                # Rtsne t-SNE降维算法 (距离矩阵、条件概率、梯度下降、动量优化)
│   ├── uwot.mbt                 # uwot UMAP降维算法 (k近邻、模糊单纯集、SGD优化、负采样)
│   ├── tradeseq.mbt             # tradeSeq 轨迹差异表达分析 (TrajectoryPoint、GAM拟合、样条基函数、差异检验)
│   ├── progeny.mbt              # PROGENy 通路活性推断 (L2正则化线性回归、Ridge回归、通路基因集权重矩阵)
│   ├── aucell.mbt               # AUCell 单细胞基因集评分 (AUC计算、基因排序、min-max归一化)
│   ├── geoquery.mbt             # GEO数据库查询 (GDS/GSE/GSM解析、数据下载、平台信息)
│   ├── ggtree.mbt               # ggtree 系统发育树可视化布局 (矩形/放射状/无根布局、节点坐标映射)
│   ├── mix_omics.mbt           # mixOmics 多组学整合 (PLS/sPLS/DIABLO)
│   ├── maf.mbt                 # MAF (Multiple Alignment Format) 多序列比对格式解析与分析
│   ├── mauve.mbt               # Mauve 基因组比对格式解析与重排分析
│   ├── stockholm.mbt           # Stockholm 格式解析 (Pfam/Rfam比对) 与二级结构分析
│   ├── popgen_advanced.mbt     # 高级群体遗传学统计 (Tajima's D, Fu & Li's D/F, MK检验)
│   ├── codon_advanced.mbt      # 高级密码子分析 (CAI, RSCU, ENC, GC3)
│   ├── pdb_packing.mbt         # 蛋白质包装密度分析 (SASA, 局部密度, Lee-Richards算法)
│   ├── qvalue.mbt              # Storey's q-value FDR方法 (π₀估计、q-value、自助法、显著性检验)
│   ├── ihw.mbt                 # 独立假设加权IHW (协变量加权Bonferroni、局部/全局加权)
│   ├── delayed_matrix_stats.mbt # DelayedArray统计层 (row/col统计、NA处理、子集操作)
│   ├── gcrma.mbt               # GC-RMA芯片分析 (背景校正、GC校正、分位数归一化、探针组汇总)
│   ├── ace.mbt                 # ACE contig组装格式解析 (reads/contigs、共有序列、覆盖度、GC含量)
│   ├── proteomics.mbt          # 蛋白质组学工具 (8种酶切、肽段质量、同位素分布、b/y碎片离子)
│   ├── fragment_mapper.mbt     # PDB片段映射 (DSSP分类、片段分配/合并/过滤、覆盖率)
│   ├── genome_diagram.mbt      # 基因组图可视化 (数据模型、SVG生成、样式控制、自动标注)
│   ├── impute.mbt              # Bioconductor impute 缺失值插补 (KNN/mean/median/LOCF/NOCB)
│   ├── vsn.mbt                 # Bioconductor vsn 方差稳定化归一化 (glog变换, vsn2)
│   ├── gsea_base.mbt           # Bioconductor GSEABase 基因集管理 (GMT/GMX格式解析、集合运算)
│   ├── pcatools.mbt            # Bioconductor PCAtools 高级PCA分析 (scree/biplot/outliers)
│   ├── data.mbt                # Bio.Data 生物数据常量 (IUPAC、氨基酸映射、密码子表)
│   ├── seq_approx.mbt          # Bio.Seq.Approximate 近似字符串匹配 (Levenshtein、错配/插入/缺失)
│   ├── pairwise2.mbt           # Bio.Pairwise2 灵活双序列比对 (NW/SW算法、自定义评分)
│   ├── compound.mbt            # Bio.Compound 化合物数据结构 (分子式解析、分子量计算)
│   ├── enhanced_volcano.mbt    # EnhancedVolcano 火山图可视化 (差异表达基因分类、ASCII渲染)
│   ├── reporting_tools.mbt     # ReportingTools 报告生成 (文本/表格/图形混合报告)
│   ├── karyoploter.mbt         # karyoploteR 核型可视化 (染色体轨道、数据点、ASCII渲染)
│   ├── system_piper.mbt        # SystemPipeR 流水线编排 (步骤管理、依赖关系、进度追踪)
│   ├── muscat.mbt              # muscat 单细胞差异状态分析 (伪批量聚合、DS检验)
│   ├── infercnv.mbt            # infercnv 单细胞CNV推断 (基因组位置平滑、参考细胞比较、CNV评分)
│   ├── scenic.mbt              # SCENIC 单细胞调控网络推断 (共表达模块、Regulon构建、AUCell活性评分)
│   ├── msstats.mbt             # MSstats 蛋白质显著性分析 (质谱数据归一化、汇总、组间比较)
│   ├── noiseq.mbt              # NOISeq 噪声鲁棒差异表达 (TMM/RPKM/上四分位归一化、NOISeqBio)
│   ├── gviz.mbt                # Gviz 基因组可视化轨道 (注释/数据/核型/序列轨道、ASCII渲染)
│   ├── htsfilter.mbt           # Bioconductor HTSFilter RNA-seq count过滤 (CPM归一化、按组最小样本阈值、keep mask)
│   ├── bayseq.mbt              # Bioconductor baySeq 贝叶斯差异表达分析 (负二项模型、分散度估计、Gamma先验、后验概率)
│   ├── cellchat.mbt            # Bioconductor CellChat 细胞间通讯分析 (配体-受体互作对、置换检验、互作评分、FDR校正)
│   ├── mutational_patterns.mbt # Bioconductor MutationalPatterns 体细胞突变谱分析 (96通道矩阵、三核苷酸上下文、突变签名拟合)
│   ├── gage.mbt                # Bioconductor GAGE 基因集富集分析 (fold change、t检验、BH-FDR校正、配对检验)
│   ├── spia.mbt                # Bioconductor SPIA 信号通路影响分析 (通路图、扰动累积、超几何检验、Fisher合并p值)
│   ├── file.mbt                # Bio.File 智能文件处理 (压缩格式检测、gzip/bzip2透明读写、文件操作接口)
│   ├── mol_wt.mbt              # Bio.SeqUtils.MolWt 分子量计算 (DNA/RNA/蛋白质分子量、消光系数、等电点)
│   ├── reduced.mbt             # Bio.Align.Reduced 简化氨基酸字母表 (RAD/Dayhoff/CHARM/SDM12、序列比较)
│   ├── hmisc.mbt               # Bioconductor Hmisc 统计工具包 (相关性分析、变量聚类、描述性统计、Somers' d、缺失值插补)
│   ├── rstatix.mbt             # Bioconductor rstatix tidy统计检验 (T检验、Wilcoxon、ANOVA、Kruskal-Wallis、Friedman、相关性检验、FDR校正)
│   ├── venn_diagram.mbt        # Bioconductor VennDiagram 集合分析 (集合运算、Venn区域、相似度指标)
│   ├── genie3.mbt              # Bioconductor GENIE3 基因调控网络推断 (回归树、特征重要性、邻接矩阵)
│   ├── decoupler.mbt           # Bioconductor decoupleR 功能活性推断 (WSum/WMean/Norm/ULM/MLM、PKN)
│   ├── bayes_space.mbt         # Bioconductor BayesSpace 空间转录组聚类 (t分布混合模型、MRF先验、EM算法)
│   └── utils.mbt               # 通用工具函数
├── examples/                   # 示例程序
│   ├── affy_demo/              # Affy Affymetrix芯片数据分析示例 (RMA标准化、背景校正、分位数归一化)
│   ├── align_info_demo/        # AlignInfo 比对统计示例 (一致性序列、保守位点、Shannon熵、成对序列同一性)
│   ├── align_io_demo/          # AlignIO 比对格式解析示例 (ClustalW、FASTA、Stockholm)
│   ├── align_applications_demo/ # Align.Applications 比对工具包装示例
│   ├── alignment_demo/         # 序列比对示例
│   ├── bioc_singular_demo/       # BiocSingular SVD演示
│   ├── align_cluster_demo/       # AlignClusterer 渐进式多序列比表示例
│   ├── application_demo/         # Bio.Application 命令行工具示例
│   ├── ballgown_demo/          # ballgown 转录组水平差异表达分析示例 (FPKM计算、t检验)
│   ├── bam_demo/               # BAM/BGZF 解析示例
│   ├── basic_seq/              # 基础序列操作示例 (含MutableSeq可变序列编辑)
│   ├── biobase_demo/           # Biobase ExpressionSet示例 (多维基因组数据容器、AnnotatedDataFrame、数据归一化、log2转换)
│   ├── bioc_generics_demo/      # Bioconductor 通用函数示例
│   ├── bioc_parallel_demo/      # Bioconductor 并行计算示例
│   ├── bioconductor_demo/      # Bioconductor模块综合示例 (ChIPseeker峰注释(外显子/内含子/UTR分类、peak2gene关联)、DOSE疾病富集、ReactomePA通路分析、AnnotationDbi注释数据库、clusterProfiler富集框架、WGCNA共表达网络、Batchelor单细胞批次校正、Seurat单细胞分析)
│   ├── chain_liftover_demo/    # UCSC Chain文件解析与基因组坐标liftOver示例
│   ├── biomart_demo/           # biomaRt 基因ID转换和注释查询示例 (基因ID映射、注释查询、批量查询)
│   ├── biostrings_demo/        # Biostrings 序列分析示例
│   ├── biostrings_matchdict_demo/ # Biostrings 字典模式匹配示例 (matchPDict、vcountPattern)
│   ├── blast_demo/             # BLAST结果解析示例 (tabular/xml格式、HSP过滤、E-value/identity过滤、最佳匹配)
│   ├── bloom_filter_demo/      # Bloom Filter & k-mer 计数示例
│   ├── bsgenome_demo/          # BSgenome 基因组序列数据库示例 (染色体序列检索、子序列提取、链特异性基因提取)
│   ├── bwt_fm_demo/            # BWT + FM-index 示例
│   ├── cluster_demo/           # 序列聚类分析示例 (层次聚类、距离矩阵、轮廓系数)
│   ├── consensus_cluster_plus_demo/ # ConsensusClusterPlus 共识聚类示例
│   ├── cyclone_demo/            # Cyclone 细胞周期评分示例 (基因对比较、G1/S/G2/M期相预测)
│   ├── codon_align_demo/       # CodonAlign 密码子比对示例 (密码子替换分类、dN/dS选择压力分析、密码子使用偏好)
│   ├── codon_usage_demo/       # CodonUsage 密码子使用分析示例 (CAI、ENC、RSCU、GC3、CBI、Fop、最优密码子检测)
│   ├── cram_demo/              # CRAM 格式解析示例 (压缩二进制序列比对格式、CRAM转BAM、参考序列管理)
│   ├── de_bruijn_demo/         # De Bruijn Graph 序列组装示例
│   ├── deseq2_demo/            # DESeq2 差异表达分析示例
│   ├── deseq2_advanced_demo/   # DESeq2 VST方差稳定化变换、PCA可视化示例
│   ├── dorothea_demo/          # dorothea 转录因子活性预测示例 (Regulon、VIPER、置换检验)
│   ├── dplyr_demo/             # dplyr 数据操作示例
│   ├── edger_demo/             # edgeR 差异表达分析示例
│   ├── edger_advanced_demo/    # edgeR准似然F检验、camera/roast基因集检验示例
│   ├── emboss_demo/            # EMBOSS工具接口示例 (GC偏斜、AT偏斜、分子量、ORF查找)
│   ├── enrichplot_demo/        # enrichplot 富集分析结果可视化示例
│   ├── ensembldb_demo/         # ensembldb Ensembl注释数据库接口示例
│   ├── expasy_demo/            # ExPASy 蛋白质分析工具接口示例
│   ├── entrez_demo/            # Entrez NCBI数据库访问示例 (ESearch、EFetch、PubMed/Gene/Taxonomy解析)
│   ├── faidx_demo/             # FASTA 索引示例
│   ├── fgsea_demo/             # fgsea 快速基因集富集分析示例 (置换检验、NES/ES计算、Leading Edge基因)
│   ├── genome_info_db_demo/    # GenomeInfoDb 基因组信息管理示例 (染色体信息、着丝粒位置、染色体臂)
│   ├── genomic_alignments_demo/ # GenomicAlignments 基因组比对分析示例 (GAlignments、coverage、summarizeOverlaps、pileup)
│   ├── genomic_ranges_demo/    # GenomicRanges 基因组区间操作示例
│   ├── genomic_ranges_advanced_demo/ # GenomicRanges tile/slidingWindows/区间运算示例
│   ├── geoquery_demo/          # GEOquery GEO数据库示例 (Series Matrix解析、SOFT格式解析、ExpressionSet转换、基因过滤)
│   ├── go_enrichment_demo/     # GOEnrichment GO功能富集分析示例 (超几何检验、BH校正、富集结果过滤)
│   ├── hmm_demo/               # Hidden Markov Model 基因预测示例
│   ├── interaction_set_demo/   # InteractionSet 染色质交互示例 (Hi-C交互、锚点对、交互矩阵、距离分布)
│   ├── isoform_switch_demo/    # IsoformSwitchAnalyzeR 转录本异构体切换分析示例
│   ├── iranges_demo/           # IRanges 区间操作示例
│   ├── kegg_demo/              # KEGG数据库解析示例 (基因、通路、化合物、酶记录)
│   ├── kmeans_demo/            # K-means Clustering 聚类分析示例
│   ├── kmer_demo/              # Bio.Kmer k-mer计数与频率分析示例
│   ├── limma_demo/             # limma 差异表达分析示例 (线性模型、经验贝叶斯、voom变换)
│   ├── matrix_demo/            # Bioconductor Matrix 稀疏矩阵操作示例
│   ├── bioc_neighbors_demo/      # BiocNeighbors 最近邻搜索演示
│   ├── medline_demo/           # Medline/PubMed解析示例 (文献记录、APA引用、MeSH过滤)
│   ├── ml_features/            # 机器学习特征提取示例
│   ├── mmcif_demo/             # mmCIF格式解析示例 (数据块解析、类别查询、原子位点提取)
│   ├── motifs_demo/            # 序列模体识别示例
│   ├── motifs_advanced_demo/   # 模体高级功能示例 (JASPAR/TRANSFAC解析、模体比对、KL/JS散度、模体聚类)
│   ├── multi_assay_experiment_demo/ # MultiAssayExperiment 多组学数据协调示例 (实验协调、样本映射)
│   ├── needleman_wunsch_demo/  # Needleman-Wunsch 全局序列比对示例
│   ├── neighbor_search_demo/   # NeighborSearch KD树近邻搜索示例 (半径搜索、最近邻、原子对搜索)
│   ├── nexus_demo/             # Nexus格式解析示例 (数据矩阵、系统发育树、距离矩阵)
│   ├── olc_demo/               # Overlap-Layout-Consensus 序列组装示例
│   ├── paml_demo/              # Bio.PAML 分子进化分析示例
│   ├── pdb_analysis_demo/      # PDB 高级结构分析示例 (主链二面角、氢键检测、二级结构分配、Ramachandran图、SASA计算、疏水性分析)
│   ├── pdb_demo/               # PDB 结构解析示例
│   ├── pdb_list_demo/          # Bio.PDB.PDBList PDB结构下载管理示例
│   ├── pdb_header_demo/        # Bio.PDB.ParsePDBHeader PDB头部元数据解析
│   ├── pdb_dice_demo/          # Bio.PDB.Dice + Selection PDB结构切割示例
│   ├── phylo_demo/             # 系统发育树示例
│   ├── phylo_metrics_demo/     # 系统发育树高级分析示例 (Colless平衡指数、Robinson-Foulds距离、距离矩阵)
│   ├── popgen_demo/            # 群体遗传学分析示例
│   ├── protparam_demo/         # ProtParam 蛋白质参数分析示例 (分子量、不稳定指数、等电点、信号肽预测、二级结构倾向)
│   ├── prosite_demo/           # Bio.Prosite 蛋白质模体数据库搜索示例
│   ├── restriction_demo/       # 限制性内切酶分析示例 (酶切位点查找、片段分析)
│   ├── residue_depth_demo/    # Bio.PDB.ResidueDepth 残基深度计算示例
│   ├── rtracklayer_demo/       # rtracklayer 基因组轨道格式示例 (BED/WIG/BEDGraph/GFF解析与写入)
│   ├── rhdf5_demo/             # Bioconductor rhdf5 HDF5文件格式支持示例
│   ├── ruvseq_demo/            # RUVSeq RNA-seq批次效应去除示例 (数据标准化、log2转换、RUVg/RUVs方法)
│   ├── sam_vcf_demo/           # SAM/VCF 解析示例
│   ├── search_io_demo/         # SearchIO 统一搜索结果示例 (HMMER3解析、BLAT PSL解析、BLAST转换)
│   ├── seq_complexity_demo/    # 序列复杂度分析示例 (Shannon熵、语言学复杂度、GC偏斜、**Wooton-Federhen LCC**、混沌游戏表示)
│   ├── seqio_demo/             # 序列 I/O 示例
│   ├── seq_utils_demo/         # 序列工具函数示例
│   ├── seqfeature_advanced_demo/  # Bio.SeqFeature CompoundLocation与LocationParser
│   ├── rna_structure_demo/       # RNA二级结构预测示例
│   ├── single_cell_demo/       # SingleCell 单细胞数据分析示例 (QC指标、Log标准化、PCA降维、高变异基因)
│   ├── single_r_demo/           # SingleR 细胞类型注释示例 (参考图谱、Spearman/Pearson相关性、精细调优)
│   ├── smith_waterman_demo/    # Smith-Waterman 局部序列比对示例
│   ├── subsmat_demo/           # 替换矩阵示例 (BLOSUM62/45、PAM250/30矩阵查询、蛋白质比对打分)
│   ├── suffix_array_tree_demo/ # Suffix Array & Suffix Tree 示例
│   ├── summarized_experiment_demo/ # SummarizedExperiment 数据容器示例
│   ├── sva_demo/               # sva 替代变量分析与ComBat批次校正示例 (经验贝叶斯方法、PCA分析)
│   ├── svd_superimposer_demo/  # SVDSuperimposer SVD蛋白质结构叠合示例 (旋转矩阵、平移向量、RMSD计算)
│   ├── structure_alignment_demo/ # Bio.PDB.StructureAlignment 多蛋白质结构比对示例
│   ├── ma_align_demo/            # 多蛋白质结构比对示例
│   ├── swissprot_demo/         # SwissProt 蛋白数据库解析示例 (记录解析、特征提取、参考文献)
│   ├── uniprot_io_demo/        # UniProt XML格式解析示例 (蛋白质条目解析、功能注释提取)
│   ├── chem_utils_demo/        # 化学计算工具示例 (键长、键角、分子式量计算)
│   ├── jaspar_demo/            # JASPAR PFM格式解析示例 (模体矩阵解析、共有序列、序列扫描)
│   ├── tree_construction_demo/ # TreeConstruction 系统发育树构建示例 (UPGMA/WPGMA/NJ算法)
│   ├── txdb_demo/              # TxDb 转录本数据库示例 (GTF解析、基因/转录本/外显子/CDS/UTR/内含子提取)
│   ├── tximport_demo/          # tximport转录本量化示例 (Salmon quant.sf解析、转录本到基因汇总、ExpressionSet转换、低表达过滤)
│   ├── variant_annotation_demo/ # VariantAnnotation 变异注释示例 (变异类型检测、定位、编码效应预测、变异汇总)
│   ├── variant_filtering_demo/ # Bioconductor VariantFiltering 变异过滤与遗传模式分析示例
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
│   ├── spatial_experiment_demo/ # SpatialExperiment 空间转录组学数据结构示例
│   ├── complex_heatmap_demo/   # ComplexHeatmap 复杂热图可视化示例 (行/列聚类、颜色映射、热图注释、分组拆分)
│   ├── gsva_demo/              # GSVA 基因集变异分析示例 (ssGSEA/zscore/PLAGE评分、富集分析、置换检验、富集图可视化、表型相关性分析、生存分析、分数分布分析)
│   ├── chromvar_demo/          # ChromVAR 染色质变异分析示例 (TF motif富集、GC偏差校正、细胞聚类、变异性分析)
│   ├── delayed_array_demo/     # DelayedArray 延迟计算数组示例 (懒加载操作、分块处理、行/列聚合、子集操作)
│   ├── annotation_filter_demo/ # AnnotationFilter 基因注释过滤示例 (染色体筛选、生物类型过滤、区域重叠检测、符号模式匹配)
│   ├── sc3_demo/                # SC3 单细胞共识聚类示例
│   ├── sc_dbl_finder_demo/     # scDblFinder 单细胞双细胞检测示例 (Doublet评分计算、最近邻搜索、PCA降维、细胞过滤)
│   ├── batchelor_demo/          # Batchelor 单细胞批次校正示例 (rescaleBatches、fastMNN、mutual nearest neighbor、批次混合评分)
│   ├── seurat_demo/             # Seurat 单细胞数据分析示例 (标准化、高可变基因、PCA、聚类、UMAP、差异标志物分析、跨样本整合)
│   ├── chipseeker_demo/         # ChIPseeker ChIP-seq峰值注释示例 (基因组区域分类(启动子/外显子/内含子/UTR/基因间区)、距离TSS分布、BED格式读取、peak2gene关联分析、注释可视化、统计分析)
│   ├── variation_demo/           # Bio.Variation 变异分析示例 (SNP分析、突变检测、氨基酸替换分析)
│   ├── dssp_demo/                # Bio.PDB.DSSP 二级结构分析示例 (DSSP解析、二级结构预测、溶剂可及表面积)
│   ├── polypeptide_demo/         # Bio.PDB.Polypeptide 多肽分析示例 (氨基酸组成、疏水性分析、跨膜区域预测)
│   ├── s4vectors_demo/           # S4Vectors 基础数据结构示例 (Rle游程编码、S4DataFrame数据框、Hits匹配)
│   ├── genefilter_demo/          # genefilter 基因过滤示例 (t检验、Wilcoxon秩和检验、方差过滤、CV过滤)
│   ├── genesis_demo/            # Genesis workflows 示例 (单细胞多组学工作流)
│   ├── graphics_demo/            # Bio.Graphics 生物信息学可视化示例 (序列Logo、比对可视化)
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
│   ├── phylo_xml_demo/          # Bio.Phylo.PhyloXML XML系统发育树格式
│   ├── phylo_nexml_demo/        # Bio.Phylo.NeXML NeXML格式解析与序列化示例
│   ├── phyloseq_demo/            # phyloseq 微生物组分析示例 (OTU丰度、分类学过滤)
│   ├── microbiome_demo/          # microbiome 微生物组分析示例 (Alpha/Beta多样性、PCoA、差异丰度)
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
│   ├── pdb_vectors_demo/         # Bio.PDB.vectors 3D向量与旋转矩阵示例 (Vector3运算、Kabsch叠合、二面角计算)
│   ├── qcp_superimposer_demo/    # Bio.PDB.QCPSuperimposer 四元数结构叠合示例
│   ├── circ_seq_demo/            # Bio.SeqUtils.CircSeq 环状DNA操作示例 (酶切分析、PCR引物设计、序列旋转)
│   ├── align_abstract_demo/      # Bio.Align.AlignAbstract 抽象比对示例 (一致性序列、Shannon熵、同一性矩阵、简约信息位点)
│   ├── maftools_demo/            # maftools 癌症基因组学示例 (MAF数据创建、突变分类、TMB计算、突变谱分析)
│   ├── cnvkit_demo/              # CNVkit 拷贝数变异示例 (CBS分割、拷贝数判定、断点检测、log2平滑)
│   ├── destiny_demo/             # destiny 扩散映射示例 (距离矩阵、高斯核、特征分解、嵌入坐标)
│   ├── rtsne_demo/               # Rtsne t-SNE降维示例 (距离矩阵、条件概率、梯度下降、动量优化)
│   ├── uwot_demo/                # uwot UMAP降维示例 (k近邻、模糊单纯集、SGD优化、负采样)
│   ├── tradeseq_demo/            # tradeSeq 轨迹差异表达示例 (GAM拟合、基因平滑、差异检验、可视化)
│   ├── progeny_demo/             # PROGENy 通路活性推断示例 (L2正则化回归、通路基因集、样本活性计算)
│   ├── aucell_demo/              # AUCell 单细胞基因集评分示例 (AUC计算、归一化、细胞/基因集评分)
│   ├── ggtree_demo/              # ggtree 系统发育树可视化示例 (矩形/放射状/无根布局、节点坐标)
│   ├── mix_omics_demo/          # mixOmics 多组学整合演示
│   ├── maf_demo/                 # MAF 格式解析与分析示例 (解析、统计、选择、过滤、写回)
│   ├── mauve_demo/               # Mauve 基因组比对分析示例 (倒位检测、断点检测、覆盖率、BED导出)
│   ├── stockholm_demo/           # Stockholm 格式解析示例 (Pfam/Rfam格式、二级结构、保守性分析)
│   ├── popgen_advanced_demo/     # 高级群体遗传学示例 (Tajima's D、Fu & Li、MK检验、中性分析)
│   ├── codon_advanced_demo/      # 高级密码子分析示例 (CAI、RSCU、ENC、最优/稀有密码子)
│   ├── pdb_packing_demo/         # 蛋白质包装密度分析示例 (SASA、包装密度、低包装区域识别)
│   ├── qvalue_demo/              # Storey's q-value FDR方法示例 (π₀估计、q-value计算、自助法)
│   ├── ihw_demo/                 # 独立假设加权IHW示例 (协变量加权、局部/全局对比、多协变量)
│   ├── delayed_matrix_stats_demo/ # DelayedArray统计层示例 (row/col统计、NA处理、子集操作)
│   ├── gcrma_demo/               # GC-RMA芯片分析示例 (背景校正、GC校正、分位数归一化)
│   ├── ace_demo/                 # ACE contig格式解析示例 (reads/contigs、共有序列、覆盖度、GC含量)
│   ├── proteomics_demo/          # 蛋白质组学分析示例 (8种酶切、质量计算、同位素分布、碎片离子)
│   ├── fragment_mapper_demo/     # PDB片段映射示例 (DSSP分类、片段分配/合并/过滤、覆盖率)
│   ├── genome_diagram_demo/      # 基因组图可视化示例 (track/feature、SVG生成、样式、自动标注)
│   ├── impute_demo/              # Bioconductor impute 缺失值插补示例
│   ├── vsn_demo/                 # Bioconductor vsn 方差稳定化归一化示例
│   ├── gsea_base_demo/           # Bioconductor GSEABase 基因集管理示例
│   ├── pcatools_demo/            # Bioconductor PCAtools 高级PCA分析示例
│   ├── data_demo/                # Bio.Data 生物数据常量示例
│   ├── seq_approx_demo/          # Bio.Seq.Approximate 近似匹配示例
│   ├── pairwise2_demo/           # Bio.Pairwise2 双序列比对示例
│   ├── compound_demo/            # Bio.Compound 化合物数据结构示例
│   ├── enhanced_volcano_demo/    # EnhancedVolcano 火山图可视化示例
│   ├── reporting_tools_demo/     # ReportingTools 报告生成示例
│   ├── karyoploter_demo/         # karyoploteR 核型可视化示例
│   ├── system_piper_demo/        # SystemPipeR 流水线编排示例
│   ├── muscat_demo/             # muscat 单细胞差异状态分析示例
│   ├── infercnv_demo/           # infercnv 单细胞拷贝数变异推断示例
│   ├── scenic_demo/             # SCENIC 单细胞调控网络推断示例
│   ├── msstats_demo/            # MSstats 蛋白质显著性分析示例
│   ├── noiseq_demo/             # NOISeq 噪声鲁棒差异表达示例
│   ├── gviz_demo/               # Gviz 基因组可视化轨道示例
│   ├── htsfilter_demo/          # HTSFilter RNA-seq count过滤示例 (CPM阈值过滤、保留率分析)
│   ├── bayseq_demo/             # baySeq 贝叶斯差异表达分析示例 (DE基因检测、后验概率)
│   ├── cellchat_demo/           # CellChat 细胞间通讯分析示例 (配体-受体互作、置换检验)
│   ├── mutational_patterns_demo/ # MutationalPatterns 体细胞突变谱分析示例 (96通道矩阵、签名拟合)
│   ├── gage_demo/                # GAGE 基因集富集分析示例 (fold change、t检验、配对检验)
│   ├── spia_demo/                # SPIA 信号通路影响分析示例 (通路扰动、超几何检验、Fisher合并)
│   ├── file_demo/                # Bio.File 智能文件处理示例 (压缩检测、透明读写、文件操作)
│   ├── mol_wt_demo/              # Bio.SeqUtils.MolWt 分子量计算示例 (DNA/RNA/蛋白质分子量、消光系数、等电点)
│   ├── reduced_demo/             # Bio.Align.Reduced 简化氨基酸字母表示例 (RAD/Dayhoff/CHARM字母表、序列比较)
│   ├── hmisc_demo/               # Hmisc 统计工具包示例 (相关性分析、变量聚类、描述性统计、Somers' d、缺失值插补)
│   ├── rstatix_demo/             # rstatix tidy统计检验示例 (T检验、Wilcoxon、ANOVA、Kruskal-Wallis、相关性检验、FDR校正)
│   ├── venn_diagram_demo/        # VennDiagram 集合分析示例 (Venn区域计算、相似度指标、集合运算)
│   ├── genie3_demo/              # GENIE3 基因调控网络推断示例 (回归树、特征重要性、调控边排序)
│   ├── decoupler_demo/           # decoupleR 功能活性推断示例 (WSum/ULM/MLM方法、PKN网络、调控子评分)
│   ├── bayes_space_demo/         # BayesSpace 空间转录组聚类示例 (六边形邻居、EM聚类、聚类可视化)
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
│   │   ├── cyclone_test.mbt
│   │   ├── deseq2_test.mbt
│   │   ├── deseq2_advanced_test.mbt
│   │   ├── dorothea_test.mbt
│   │   ├── dplyr_test.mbt
│   │   ├── edger_test.mbt
│   │   ├── edger_advanced_test.mbt
│   │   ├── faidx_test.mbt
│   │   ├── feature_extraction_test.mbt
│   │   ├── genomic_alignments_test.mbt
│   │   ├── genomic_ranges_test.mbt
│   │   ├── genomic_ranges_advanced_test.mbt
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
│   │   ├── seqfeature_advanced_test.mbt
│   │   ├── seqio_wb_test.mbt
│   │   ├── seq_utils_test.mbt
│   │   ├── sequtils_test.mbt
│   │   ├── single_cell_test.mbt
│   │   ├── single_r_test.mbt
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
│   │   ├── blast_applications_test.mbt
│   │   ├── genome_info_db_test.mbt
│   │   ├── interaction_set_test.mbt
│   │   ├── multi_assay_experiment_test.mbt
│   │   ├── tree_construction_test.mbt
│   │   ├── neighbor_search_test.mbt
│   │   ├── bioc_neighbors_test.mbt
│   │   ├── bioc_singular_test.mbt
│   │   ├── swissprot_test.mbt
│   │   ├── mmcif_test.mbt
│   │   ├── align_abstract_test.mbt
│   │   ├── align_analysis_test.mbt
│   │   ├── align_applications_test.mbt
│   │   ├── alphabet_test.mbt
│   │   ├── annotation_dbi_test.mbt
│   │   ├── annotation_filter_test.mbt
│   │   ├── annotation_hub_test.mbt
│   │   ├── application_test.mbt
│   │   ├── aucell_test.mbt
│   │   ├── bamsignals_test.mbt
│   │   ├── batchelor_test.mbt
│   │   ├── biobase_test.mbt
│   │   ├── bioc_generics_test.mbt
│   │   ├── bioc_parallel_test.mbt
│   │   ├── bsseq_test.mbt
│   │   ├── checksum_test.mbt
│   │   ├── chipseeker_test.mbt
│   │   ├── chromvar_test.mbt
│   │   ├── circ_seq_test.mbt
│   │   ├── cluster_profiler_test.mbt
│   │   ├── cnvkit_test.mbt
│   │   ├── complex_heatmap_test.mbt
│   │   ├── consensus_cluster_plus_test.mbt
│   │   ├── csaw_test.mbt
│   │   ├── delayed_array_test.mbt
│   │   ├── destiny_test.mbt
│   │   ├── rtsne_test.mbt
│   │   ├── uwot_test.mbt
│   │   ├── dexseq_test.mbt
│   │   ├── diffbind_test.mbt
│   │   ├── dose_test.mbt
│   │   ├── droplet_utils_test.mbt
│   │   ├── dss_test.mbt
│   │   ├── dssp_test.mbt
│   │   ├── edaseq_test.mbt
│   │   ├── emboss_test.mbt
│   │   ├── enrichplot_test.mbt
│   │   ├── ensembldb_test.mbt
│   │   ├── expasy_test.mbt
│   │   ├── flow_core_test.mbt
│   │   ├── freq_analysis_test.mbt
│   │   ├── genefilter_test.mbt
│   │   ├── genesis_test.mbt
│   │   ├── genomic_features_test.mbt
│   │   ├── genomic_files_test.mbt
│   │   ├── geoquery_test.mbt
│   │   ├── gff_test.mbt
│   │   ├── ggtree_test.mbt
│   │   ├── graph_test.mbt
│   │   ├── graphics_test.mbt
│   │   ├── gsva_test.mbt
│   │   ├── hilbertcurve_test.mbt
│   │   ├── hs_exposure_test.mbt
│   │   ├── isoform_switch_analyze_r_test.mbt
│   │   ├── kmer_test.mbt
│   │   ├── maftools_test.mbt
│   │   ├── mast_test.mbt
│   │   ├── matrix_test.mbt
│   │   ├── medline_test.mbt
│   │   ├── melting_temp_test.mbt
│   │   ├── metagenomeseq_test.mbt
│   │   ├── motifs_test.mbt
│   │   ├── motifs_advanced_test.mbt
│   │   ├── minfi_test.mbt
│   │   ├── monocle3_test.mbt
│   │   ├── nexus_test.mbt
│   │   ├── nucle_r_test.mbt
│   │   ├── paml_test.mbt
│   │   ├── pathway_test.mbt
│   │   ├── pdb_analysis_test.mbt
│   │   ├── pdb_dice_test.mbt
│   │   ├── pdb_list_test.mbt
│   │   ├── pdb_header_test.mbt
│   │   ├── pdb_vectors_test.mbt
│   │   ├── phenotype_test.mbt
│   │   ├── phylo_consensus_test.mbt
│   │   ├── phylo_metrics_test.mbt
│   │   ├── phylo_nexml_test.mbt
│   │   ├── phylo_test.mbt
│   │   ├── phylo_xml_test.mbt
│   │   ├── phylo_xml_debug_test.mbt
│   │   ├── phyloseq_test.mbt
│   │   ├── microbiome_test.mbt
│   │   ├── polypeptide_test.mbt
│   │   ├── progeny_test.mbt
│   │   ├── prosite_test.mbt
│   │   ├── psea_test.mbt
│   │   ├── qcp_superimposer_test.mbt
│   │   ├── reactome_pa_test.mbt
│   │   ├── residue_depth_test.mbt
│   │   ├── rhdf5_test.mbt
│   │   ├── s4vectors_test.mbt
│   │   ├── sc3_test.mbt
│   │   ├── sc_dbl_finder_test.mbt
│   │   ├── scater_test.mbt
│   │   ├── scnorm_test.mbt
│   │   ├── scran_test.mbt
│   │   ├── search_io_test.mbt
│   │   ├── searchio_new_test.mbt
│   │   ├── seq_complexity_test.mbt
│   │   ├── seurat_test.mbt
│   │   ├── sff_io_test.mbt
│   │   ├── short_read_test.mbt
│   │   ├── seq_quality_trim_test.mbt
│   │   ├── single_cell_experiment_test.mbt
│   │   ├── slingshot_test.mbt
│   │   ├── spatial_experiment_test.mbt
│   │   ├── statistics_test.mbt
│   │   ├── structure_alignment_test.mbt
│   │   ├── taxonomy_test.mbt
│   │   ├── topgo_test.mbt
│   │   ├── tradeseq_test.mbt
│   │   ├── tximport_test.mbt
│   │   ├── universalmotif_test.mbt
│   │   ├── variant_filtering_test.mbt
│   │   ├── variation_test.mbt
│   │   ├── align_cluster_test.mbt
│   │   ├── rna_structure_test.mbt
│   │   ├── mix_omics_test.mbt
│   │   ├── maf_test.mbt
│   │   ├── mauve_test.mbt
│   │   ├── stockholm_test.mbt
│   │   ├── popgen_advanced_test.mbt
│   │   ├── codon_advanced_test.mbt
│   │   ├── pdb_packing_test.mbt
│   │   ├── qvalue_test.mbt
│   │   ├── ihw_test.mbt
│   │   ├── delayed_matrix_stats_test.mbt
│   │   ├── gcrma_test.mbt
│   │   ├── ace_test.mbt
│   │   ├── proteomics_test.mbt
│   │   ├── fragment_mapper_test.mbt
│   │   ├── genome_diagram_test.mbt
│   │   ├── impute_test.mbt
│   │   ├── vsn_test.mbt
│   │   ├── gsea_base_test.mbt
│   │   ├── pcatools_test.mbt
│   │   ├── data_test.mbt
│   │   ├── seq_approx_test.mbt
│   │   ├── pairwise2_test.mbt
│   │   ├── compound_test.mbt
│   │   ├── enhanced_volcano_test.mbt
│   │   ├── reporting_tools_test.mbt
│   │   ├── karyoploter_test.mbt
│   │   ├── system_piper_test.mbt
│   │   ├── muscat_test.mbt
│   │   ├── infercnv_test.mbt
│   │   ├── scenic_test.mbt
│   │   ├── msstats_test.mbt
│   │   ├── noiseq_test.mbt
│   │   ├── gviz_test.mbt
│   │   ├── htsfilter_test.mbt
│   │   ├── bayseq_test.mbt
│   │   ├── cellchat_test.mbt
│   │   ├── mutational_patterns_test.mbt
│   │   ├── gage_test.mbt
│   │   ├── spia_test.mbt
│   │   ├── file_test.mbt
│   │   ├── mol_wt_test.mbt
│   │   ├── reduced_test.mbt
│   │   ├── hmisc_test.mbt
│   │   ├── rstatix_test.mbt
│   │   ├── venn_diagram_test.mbt
│   │   ├── genie3_test.mbt
│   │   ├── decoupler_test.mbt
│   │   ├── bayes_space_test.mbt
│   │   └── ma_align_test.mbt
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
moon test --package IvanAXu/BioSeqs/test/moonbit        # ✅ 4916 个测试全部通过
```

### 模块对照表

#### 序列处理与 I/O

| MoonBit 文件 | 对应 Python 库 | 核心功能 |
| :--- | :--- | :--- |
| `seq.mbt` | BioPython `Bio.Seq` | 序列对象、互补、转录、翻译 |
| `seq_record.mbt` | BioPython `Bio.SeqRecord` | 带注释的序列记录 |
| `seqfeature.mbt` | BioPython `Bio.SeqFeature` | 序列特征与位置 |
| `seqfeature_advanced.mbt` | BioPython `Bio.SeqFeature` | CompoundLocation复合位点、LocationParser位点字符串解析、扩展SeqFeature修饰符 |
| `seqio.mbt` | BioPython `Bio.SeqIO` | 统一序列文件 I/O（seqio_parse / seqio_write / seqio_convert / seqio_read / seqio_to_dict） |
| `fasta_io.mbt` | BioPython `Bio.SeqIO.FastaIO` | FASTA 解析与写入 |
| `fastq_io.mbt` | BioPython `Bio.SeqIO.QualityIO` | FASTQ 解析与写入 |
| `genbank_io.mbt` | BioPython `Bio.SeqIO.GenBankIO` | GenBank 解析与写入 |
| `embl_io.mbt` | BioPython `Bio.SeqIO.EmblIO` | EMBL 格式解析（ID / AC / DE / SQ 字段） |
| `pir_io.mbt` | BioPython `Bio.SeqIO.PdbIO (PIR/NBRF)` | PIR / NBRF 蛋白与核酸序列格式解析 |
| `tab_io.mbt` | BioPython `Bio.SeqIO.TabIO` | Tab 分隔（ID + 序列）解析与写入 |
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
| `phylo_xml.mbt` | BioPython `Bio.Phylo.PhyloXML` | PhyloXML格式解析、序列化、Newick双向转换、分类单元注释 |

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
| `pdb_header.mbt` | BioPython `Bio.PDB.ParsePDBHeader` | PDB头部元数据解析 (HEADER/TITLE/COMPOUND/SOURCE/REMARK/AUTH/DBREF) |

#### 基因组分析

| MoonBit 文件 | 对应 Python 库 | 核心功能 |
| :--- | :--- | :--- |
| `sam.mbt` | pysam | SAM 文件解析 |
| `bam.mbt` | pysam | BAM 文件解析 |
| `bgzf.mbt` | pysam | BGZF 解压缩 |
| `vcf.mbt` | pysam | VCF 文件解析 |
| `cram_wbtest.mbt` | pysam | CRAM 格式解析 |
| `genomic_ranges.mbt` | Bioconductor GenomicRanges | GRanges 区间操作 |
| `genomic_ranges_advanced.mbt` | GenomicRanges Tile/Windows | tile分箱、sliding_windows滑窗、tile_genome基因组覆盖、coverage_by_window覆盖度计算、bin_genome分箱统计、promoters启动子、gaps间隙、subtract区间减法 |
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
| `deseq2_advanced.mbt` | DESeq2 VST/PCA | VST方差稳定化变换、快速VST、PCA主成分分析、高变基因筛选 |
| `edger.mbt` | Bioconductor edgeR | DGEList 差异表达 |
| `edger_advanced.mbt` | edgeR QLF/Camera/Roast | 准似然F检验、QL分散度估计、camera竞争性基因集检验、roast自足基因集检验 |
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
| `uniprot_io.mbt` | BioPython `Bio.SeqIO.UniprotIO` | UniProt XML 格式解析 |
| `chem_utils.mbt` | BioPython `Bio.PDB.chem_utils` | 化学计算工具（键长、键角、二面角、分子式量） |
| `jaspar.mbt` | BioPython `Bio.motifs.Jaspar` | JASPAR PFM 格式解析与模体分析 |
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
| `rtsne.mbt` | Bioconductor Rtsne | t-SNE降维算法 (TsneConfig/TsneResult数据结构、成对距离计算、条件概率估计与perplexity优化、联合概率矩阵构建、梯度下降优化、动量/early exaggeration调度) |
| `uwot.mbt` | Bioconductor uwot | UMAP降维算法 (UmapConfig/UmapResult数据结构、k近邻搜索、模糊单纯集构建、局部模糊集并集、SGD低维嵌入优化、负采样、min_dist/spread参数控制) |
| `tradeseq.mbt` | Bioconductor tradeSeq | 轨迹差异表达分析 (TrajectoryPoint/GeneExpressionData/GAMFit/DifferentialExpressionResult数据结构、GAM广义可加模型拟合、样条基函数、条件效应检验、BH-FDR校正) |
| `maf.mbt` | BioPython `Bio.Align` | MAF多序列比对格式解析与分析 (块处理、百分比一致性、统计分析、选择/过滤) |
| `mauve.mbt` | BioPython `Bio.Align` | Mauve基因组比对格式解析 (LCB检测、倒位/断点、覆盖率、BED导出) |
| `stockholm.mbt` | BioPython `Bio.Stockholm` | Stockholm格式解析与分析 (Pfam/Rfam格式、二级结构、保守性、FASTA转换) |
| `popgen_advanced.mbt` | BioPython `Bio.PopGen` | 高级群体遗传学统计 (Tajima's D、Fu & Li's D/F、MK检验、等位基因频率谱) |
| `codon_advanced.mbt` | BioPython `Bio.SeqUtils.CodonUsage` | 高级密码子分析 (CAI、RSCU、ENC、GC3、物种特异性参考表) |
| `pdb_packing.mbt` | BioPython `Bio.PDB.Packing` | 蛋白质包装密度分析 (SASA计算、Lee-Richards算法、低包装区域识别) |
| `qvalue.mbt` | Bioconductor qvalue | Storey's q-value FDR方法 (π₀估计、自助法、q-value计算、显著性检验) |
| `ihw.mbt` | Bioconductor IHW | 独立假设加权 (协变量加权Bonferroni、局部/全局加权、多协变量支持) |
| `delayed_matrix_stats.mbt` | Bioconductor DelayedMatrixStats | DelayedArray统计层 (row/col统计、NA处理、子集操作) |
| `gcrma.mbt` | Bioconductor gcrma | GC-RMA芯片分析 (背景校正IdealMM/Express、GC查找表、分位数归一化、探针组汇总) |
| `ace.mbt` | Biopython `Bio.Sequencing.Ace` | ACE contig格式解析 (reads/contigs解析、共有序列、覆盖度分析、GC含量) |
| `proteomics.mbt` | Biopython `Bio.SeqUtils.Proteomics` | 蛋白质组学工具 (8种酶切、肽段质量、同位素分布、b/y碎片离子) |
| `fragment_mapper.mbt` | Biopython `Bio.PDB.FragmentMapper` | PDB片段映射 (DSSP分类、片段分配/合并/过滤、覆盖率分析) |
| `genome_diagram.mbt` | Biopython `Bio.Graphics.GenomeDiagram` | 基因组图可视化 (数据模型、SVG生成、样式控制、自动标注) |

#### Bioconductor 新增扩展模块

| MoonBit 文件 | 对应 Bioconductor 包 | 核心功能 |
| :--- | :--- | :--- |
| `impute.mbt` | `impute` | 缺失值插补（KNN / mean / median / LOCF / NOCB） |
| `vsn.mbt` | `vsn` | 方差稳定化归一化（glog 变换, vsn2, mean-SD） |
| `gsea_base.mbt` | `GSEABase` | 基因集数据结构 & GMT/GMX 格式解析、集合运算 |
| `pcatools.mbt` | `PCAtools` | 高级 PCA 分析（scree/biplot/outliers/correlations） |
| `data.mbt` | `Bio.Data` | IUPAC 碱基、氨基酸映射、密码子表、反向互补 |
| `seq_approx.mbt` | `Bio.Seq.Approximate` | 近似字符串匹配（Levenshtein, 错配/插入/缺失） |
| `pairwise2.mbt` | `Bio.Pairwise2` | 灵活双序列比对（Needleman-Wunsch / Smith-Waterman） |
| `compound.mbt` | `Bio.Compound` | 化合物/反应/代谢通路数据结构、分子式解析 |
| `enhanced_volcano.mbt` | `EnhancedVolcano` | 火山图可视化（差异表达基因分类、ASCII 渲染） |
| `reporting_tools.mbt` | `ReportingTools` | 报告生成（文本/表格/图形混合报告、ASCII 表格） |
| `karyoploter.mbt` | `karyoploteR` | 核型可视化（染色体轨道、数据点、ASCII 渲染） |
| `system_piper.mbt` | `SystemPipeR` | 流水线编排（步骤管理、依赖关系、进度追踪） |
| `muscat.mbt` | `muscat` | 单细胞差异状态分析（伪批量聚合、DS 检验、QC） |
| `infercnv.mbt` | `infercnv` | 单细胞拷贝数变异推断（基因组位置排序、参考细胞有界 LFC 计算、金字塔权重平滑、CNV 分数与恶性细胞预测） |
| `scenic.mbt` | `SCENIC` | 单细胞调控网络推断与聚类（TF-target 共表达模块、Regulon 构建、AUCell 活性评分、二值化阈值、细胞状态与主控调控因子） |
| `msstats.mbt` | `MSstats` | 蛋白质显著性分析（质谱归一化、Tukey/Linear 汇总、组间比较） |
| `noiseq.mbt` | `NOISeq` | 噪声鲁棒差异表达（TMM/RPKM/上四分位归一化、NOISeqBio/Sim） |
| `gviz.mbt` | `Gviz` | 基因组可视化轨道（注释/数据/核型/序列轨道、ASCII 渲染） |
| `mutational_patterns.mbt` | `MutationalPatterns` | 体细胞突变谱分析（96通道矩阵、三核苷酸上下文归一化、突变签名拟合、余弦相似度） |
| `gage.mbt` | `gage` | 基因集富集分析（fold change、t检验、BH-FDR校正、配对/非配对检验） |
| `spia.mbt` | `SPIA` | 信号通路影响分析（通路图建模、扰动累积、超几何检验、Fisher合并p值、激活/抑制判定） |
| `file.mbt` | `Bio.File` | 智能文件处理（压缩格式自动检测、gzip/bzip2透明读写、文件操作接口） |
| `mol_wt.mbt` | `Bio.SeqUtils.MolWt` | 分子量计算（DNA/RNA/蛋白质分子量、消光系数、吸光度、等电点） |
| `reduced.mbt` | `Bio.Align.Reduced` | 简化氨基酸字母表（RAD/Dayhoff/CHARM/SDM12字母表、序列比较、简化一致性） |

## 核心功能实现

### 1. 序列处理 (Bio.Seq)

提供完整的序列对象支持，包括 DNA、RNA 和蛋白质序列的创建与操作。支持互补、反向互补、转录、反转录和翻译等核心生物信息学操作。翻译功能支持多种模式：普通翻译、翻译到终止密码子、完整 CDS 翻译等，确保满足不同的序列分析需求。

### 2. 序列 I/O (Bio.SeqIO)

实现统一的序列文件解析接口，支持 **FASTA、FASTQ、GenBank、EMBL、PIR/NBRF、Tab-delimited 六种常用格式的解析与写入（部分格式仅读或仅写）。通过统一的 API 设计，用户可以轻松切换不同的文件格式，无需关注底层实现细节，极大简化了序列数据的处理流程。

- **seqio_parse(content, format)**：统一读取，支持 `fasta` / `fastq` / `fastq-sanger` / `genbank` / `gb` / `embl` / `pir` / `nbrf` / `tab`
- **seqio_write(records, format)**：统一写出，支持 `fasta` / `fastq` / `genbank` / `gb` / `tab`
- **seqio_convert(content, in_format, out_format)**：格式直接互换（如 FASTA → Tab → GenBank 来回转换）
- **seqio_read(content, format)**：读单条记录（多于一条时报错）
- **seqio_to_dict(records)**：按 id 建索引 Map，方便按名字随机访问

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

### 32. limma 差异表达分析与批次校正 (Bioconductor limma)

实现基于线性模型的差异表达分析，支持经验贝叶斯校正和 voom 变换（适用于 RNA-seq 数据）。支持设计矩阵创建、线性模型拟合、对比矩阵分析和 top 差异基因的获取。新增 RPKM/CPM/quantile 归一化、ComBat 批次校正、removeBatchEffect、treat 严格检验。

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

### 118. microbiome 微生物组分析 (Bioconductor microbiome)

实现微生物组数据分析的完整功能，包括 Alpha 多样性、Beta 多样性、排序分析和差异丰度分析。支持 AlphaDiversity 数据结构存储多种多样性指数。Alpha 多样性指数包括：Observed OTUs、Shannon、Simpson、Inverse Simpson、Pielou evenness、Chao1、ACE 和 Fisher's alpha。Beta 多样性度量包括：Bray-Curtis 相异度、Jaccard 距离、Jensen-Shannon 散度、加权和非加权 UniFrac。支持 PCoA（主坐标分析）进行排序，采用经典 MDS（度量多维标度）算法，通过幂迭代计算特征值和特征向量。支持差异丰度分析，包括 Welch's t 检验和 Wilcoxon 秩和检验，以及 Benjamini-Hochberg 多重检验校正。支持 Alpha 多样性表计算和 Beta 多样性距离矩阵生成。适用于 16S rRNA 测序和宏基因组学数据的微生物群落分析。

### 119. Rtsne t-SNE降维算法 (Bioconductor Rtsne)

实现 t-Distributed Stochastic Neighbor Embedding（t-SNE）降维算法，用于高维数据的可视化和探索。支持 TsneConfig 数据结构配置算法参数（perplexity、theta、max_iter、dims、eta、exaggeration_factor、momentum、final_momentum、mom_switch_iter、stop_lying_iter、random_seed）。支持 TsneResult 数据结构存储降维结果（embedding、costs、n_iter）。核心算法包括：成对距离矩阵计算、条件概率估计（使用二分搜索优化 perplexity）、联合概率矩阵 P 构建、梯度下降优化（带动量和 early exaggeration）。支持 Barnes-Hut 近似参数配置。支持伪随机数生成和测试数据创建（3 个聚类的测试数据集）。适用于单细胞 RNA-seq、蛋白质组学等高维生物数据的非线性降维可视化。

### 120. uwot UMAP降维算法 (Bioconductor uwot)

实现 Uniform Manifold Approximation and Projection（UMAP）降维算法，提供快速且保留全局结构的降维方法。支持 UmapConfig 数据结构配置算法参数（n_neighbors、n_components、n_epochs、min_dist、spread、learning_rate、repulsion_strength、negative_sample_rate、random_seed、metric）。支持 UmapResult 数据结构存储降维结果（embedding、n_epochs）。核心算法包括：k 近邻搜索（KNN）、模糊单纯集构建（fuzzy simplicial set）、局部模糊集并集、低维嵌入随机梯度下降（SGD）优化、负采样（negative sampling）。支持 min_dist 和 spread 参数控制嵌入点的分布。支持测试数据创建（4 个聚类的测试数据集）。适用于单细胞转录组学、空间转录组学等大规模生物数据的快速降维和可视化。

### 121. tradeSeq 轨迹差异表达分析 (Bioconductor tradeSeq)

实现基于轨迹的差异表达分析功能，用于识别沿细胞分化轨迹表达变化的基因。支持 TrajectoryPoint 和 GeneExpressionData 数据结构存储拟时间和基因表达数据。支持 GAM（Generalized Additive Model，广义可加模型）拟合，使用样条基函数（Cubic Spline Basis）对基因表达沿拟时间进行平滑拟合。支持最小二乘法系数估计、拟合值计算和残差分析。支持 R-squared 计算评估模型拟合质量。支持条件效应差异检验（trade_test_condition_effect），比较不同条件下基因沿轨迹的表达差异。支持完整分析流程（run_tradeseq_analysis），包含批量基因检验和 Benjamini-Hochberg 多重检验校正。支持基因表达平滑曲线计算（calculate_gene_smooth）。适用于单细胞轨迹分析中差异表达基因的识别。

### 122. MAF 多序列比对格式解析与分析 (Bio.Align.MAF)

实现 MAF（Multiple Alignment Format）格式的完整解析与分析功能，这是基因组多序列比对的标准格式。支持 MAFAliignment、MAFBlock 和 MAFSequence 数据结构存储比对信息。支持从字符串解析 MAF 内容（parse_maf），自动识别版本号、得分、比对块和序列信息。支持 MAF 格式写回（write_maf），可将比对对象序列化为标准 MAF 字符串。支持块级分析：百分比一致性计算（maf_block_percent_identity）、块比对长度统计、块序列数量统计。支持全比对统计：平均百分比一致性、总比对长度、覆盖度分析、块长分布（最小值、最大值、中位数）。支持序列选择（maf_select_seqs）和按长度过滤（maf_filter_by_length），便于从全基因组比对中提取特定子集。支持序列覆盖率计算，评估各序列在基因组上的覆盖情况。适用于全基因组比对分析、保守区域识别和比较基因组学研究。

### 123. Mauve 基因组比对格式解析与重排分析 (Bio.Align.Mauve)

实现 Mauve 基因组比对格式的解析与基因组重排分析功能，用于多基因组共线性比较。支持 MauveAlignment、MauveLCB 和 MauveSequence 数据结构存储比对信息。支持 Mauve 格式解析（parse_mauve），识别 LCB（Locally Collinear Blocks，共线性块）和序列信息。支持基因组重排检测：倒位检测（detect_mauve_inversions）通过比较相邻序列的方向变化识别基因组倒位；断点检测（detect_mauve_breakpoints）通过分析 LCB 边界识别重排断点。支持基因组覆盖率计算（mauve_genome_coverage），量化各基因组在比对中的覆盖比例。支持 Mauve 摘要报告生成（mauve_summary），综合展示比对统计、倒位数量、断点数量和覆盖度信息。支持 BED 格式导出（mauve_to_bed、mauve_inversions_to_bed），便于在基因组浏览器中可视化。支持按得分过滤 LCB（mauve_filter_lcbs）和保守片段分析。支持基因组重排率计算，量化每 Mb 的断点数量。适用于比较基因组学研究、基因组结构变异分析和进化重排研究。

### 124. Stockholm 格式解析与二级结构分析 (Bio.Stockholm)

实现 Stockholm 格式的完整解析与分析功能，这是 Pfam 和 Rfam 数据库使用的标准比对格式。支持 StockholmAlignment、StockholmBlock 和 StockholmSequence 数据结构存储比对及注释信息。支持 Stockholm 格式解析（stockholm_parse），自动识别版本号、块结构、GF/GC/GR/GS 注释行和二级结构信息。支持 Stockholm 格式写回（stockholm_write），保留所有注释和二级结构信息。支持百分比一致性计算（stockholm_percent_identity），排除空位后计算序列间的一致性。支持每列保守性分析（stockholm_conservation），评估比对中每个位置的保守程度。支持 FASTA 格式转换（stockholm_to_fasta），便于与其他分析工具互操作。支持块合并（stockholm_merge_blocks），将多块比对合并为单一块。支持二级结构注释（#=GC SS_cons）解析，提取共识二级结构信息。适用于蛋白质家族分析、RNA 家族比对和功能保守性研究。

### 125. 高级群体遗传学统计与中性检验 (Bio.PopGen Advanced)

实现高级群体遗传学统计分析功能，包括核心中性检验和进化分析方法。支持 PolymorphicSite 和 NeutralityTestResult 等数据结构存储多态位点和检验结果。支持 Tajima's D 检验（popgen_tajima_d），通过比较 θ_π 和 θ_W 检测选择信号或群体人口结构变化：负值表示低频变异过量（正向选择或群体扩张），正值表示低频变异不足（平衡选择或群体收缩）。支持 Fu & Li's D 检验（fu_li_d）和 F 检验（fu_li_f），利用 singletons（单例变异）信息检测近期选择事件。支持 McDonald-Kreitman 检验（mcdonald_kreitman_test），比较种内多态性和种间固定化比率，检测适应性进化。支持 Watterson's θ 估计（popgen_watterson_theta），基于分离位点数估计种群遗传多样性。支持等位基因频率谱计算（calculate_afs），分析单例、双例等不同频率类别的变异数量。支持综合中性分析（run_neutrality_analysis），一次性运行多种检验并生成汇总报告。适用于分子进化研究、选择信号检测和群体遗传学分析。

### 126. 高级密码子使用分析 (Bio.SeqUtils.CodonUsage Advanced)

实现高级密码子使用分析功能，支持多种密码子偏好指标的计算。支持 CodonUsageTable 和 CAIResult 等数据结构存储密码子使用频率和分析结果。支持从 DNA 序列构建密码子使用表（codon_usage_from_sequence、codon_usage_from_sequences），统计各密码子的计数和频率。支持密码子适应指数 CAI 计算（calculate_cai），评估序列的密码子使用偏好与参考密码子表的匹配程度，CAI 值越高表示越适配宿主生物的密码子偏好。支持相对同义密码子使用 RSCU 计算（calculate_rscu），量化各密码子的使用偏离程度。支持有效密码子数 ENC 计算（calculate_enc），衡量密码子使用偏好的强度，ENC 值越低表示偏好越强。支持 GC3 含量计算（密码子第三位置 GC 含量），作为密码子偏好的指标。支持最优密码子检测（optimal_codons）和稀有密码子识别（rare_codons）。支持多种生物物种特异性参考表（高 GC、低 GC、通用）。适用于基因表达优化、密码子偏好研究和进化分析。

### 127. 蛋白质包装密度分析 (Bio.PDB.Packing)

实现蛋白质包装密度分析功能，用于评估蛋白质内部的紧密程度和溶剂可及性。支持 PackingAtom 和 PackingResult 等数据结构存储原子信息和密度分析结果。支持局部包装密度计算（calculate_packing_density），基于 Lee-Richards 算法通过在原子表面生成采样点来估计局部包装紧密程度。支持 per-residue 包装分析（packing_density_per_residue），计算每个残基的包装密度，识别蛋白质的埋藏区域和表面区域。支持 SASA 溶剂可及表面积计算（calculate_packing_sasa），支持自定义探针半径（默认为 1.4 Å，对应水分子）。支持范德华半径查询（get_packing_vdw_radius），覆盖常见原子类型（H、C、N、O、S、P 等）。支持球体点生成（generate_packing_sphere_points），使用黄金螺旋算法在单位球面上均匀分布采样点。支持密度标准化（normalize_packing_density）和低包装区域识别（identify_low_packing），用于发现蛋白质内部的空腔和弱包装区域。支持球体几何计算（球体积、半径互算）。适用于蛋白质结构分析、酶催化位点研究和蛋白质折叠稳定性评估。

### 128. Storey's q-value FDR 方法 (Bio.Statistics.qvalue)

实现 Storey's q-value FDR（错误发现率）校正方法，是多重检验校正的基础工具。支持 QValueResult 数据结构存储 p-values、q-values、π₀、λ 和显著性信息。支持 q-value 主计算（qvalue），自动选择最优 λ 值（通过自助法 MSE 最小化），估计 π₀（原假设比例），并计算 Storey q-values。支持指定 λ 的 q-value 计算（qvalue_with_lambda），允许用户控制 π₀ 估计参数。支持 π₀ 估计（qvalue_estimate_pi0），在给定 λ 下估计原假设比例。支持自助法 π₀ 估计（bootstrap_pi0），通过重复采样提高估计准确性。支持仅 FDR 校正（storey_fdr），返回 q-values 数组。支持显著性判定（qvalue_significance），按阈值获取显著性布尔数组。算法核心：按 p 值升序排序后反向计算 q-value = min(q_(k+1), p_(k) × π₀ × n / (k+1))。适用于基因组学、蛋白质组学等大规模多重检验场景。

### 129. 独立假设加权 IHW (Bio.Statistics.IHW)

实现 IHW（Independent Hypothesis Weighting）方法，利用协变量信息为不同假设分配不同权重，在控制 FDR 的同时提高统计功效。支持 IHWResult 数据结构存储调整后 p-values、权重和协变量信息。支持 IHW 主过程（ihw），使用默认配置（α=0.05，局部加权，最多 10 次迭代）。支持自定义配置 IHW（ihw_with_config），提供 IHWConfig 支持 α、迭代次数和加权类型（global/local）配置。支持多协变量 IHW（ihw_multiple），通过主成分合并多个协变量信息。支持 Bonferroni 基线校正（bonferroni_ihw）作为对比基准。支持 Storey 风格 IHW（storey_ihw），结合协变量加权的 π₀ 估计。支持权重提取（ihw_weights）查看各假设的分配权重。算法核心：加权 Bonferroni α_i = weight_i × α / sum(weights)，通过迭代过程从拒绝集估计最优权重。适用于具有协变量信息的大规模多重检验，如基于基因长度、GC 含量等协变量的差异表达分析。

### 130. DelayedMatrixStats DelayedArray 统计层 (Bioconductor DelayedMatrixStats)

实现 DelayedArray 的高效行/列统计计算层，提供对延迟计算矩阵的统计分析能力。支持 DelayedMatrix 数据结构存储矩阵数据（行名、列名、维度），MatrixStatsResult 存储统计结果。支持通用行统计（row_stats）和列统计（col_stats），涵盖 mean、var、sd、median、min、max、sum、n、nna、nn 等 10 种统计量。支持便捷的行/列均值（dms_row_means、dms_col_means）、中位数（row_medians、col_medians）、方差（dms_row_vars、col_vars）函数。支持 NA 值处理，所有统计函数自动跳过缺失值。支持行求和（row_allsums）和非 NA 值计数（row_n、col_n）。支持矩阵子集操作（delayed_matrix_subset），保留行名和列名信息。采用样本方差（n-1 分母）计算，与 Bioconductor DelayedMatrixStats 保持一致。适用于大规模基因表达矩阵的批量统计分析。

### 131. GC 校正 RMA 芯片分析 (Bioconductor gcrma)

实现 GC-RMA（GC-content adjusted Robust Multi-array Average）方法，用于 Affymetrix 微阵列数据的 GC 含量校正分析。支持 GCRMAConfig 配置（背景方法、归一化开关）、GCRMAResult 结果存储（表达矩阵、探针信息、GC 校正参数）和 ProbeInfo 探针信息结构（ID、GC 计数、序列、亲和力）。支持背景校正（gcrma_background_correction），提供 IdealMM 和 Express 两种方法，均考虑 GC 含量依赖性。支持 GC 校正（gcrma_gc_correction），基于预计算的 GC 查找表对探针强度进行 GC 依赖性调整。支持 GC 查找表计算（compute_gc_lookup_table），存储 GC count 到中位亲和力的映射。支持分位数归一化（gcrma_normalize），使各样本强度分布一致。支持完整 GCRMA 流水线（gcrma_process、gcrma_process_with_config），从背景校正→GC 校正→分位数归一化一站式完成。支持探针组汇总（gcrma_summarize），通过列中位数法汇总探针组到基因水平。支持探针亲和力估计（estimate_affinity）和 GC 含量计算（oligo_gc_count、oligo_gc_fraction）。适用于 Affymetrix 芯片数据的预处理和差异表达分析。

### 132. ACE contig 格式解析 (Bio.Sequencing.Ace)

实现 ACE 组装格式（ACE assembly format）的解析与分析功能。ACE 格式被广泛用于基因组组装项目（如 Celera Assembler、ARACHNE 等）的 contig 数据存储。支持 AceRead（读取记录，含 read_id、sequence、quality、clip_start/end、strand、chemistry、dye）、AceContig（重叠群，含 contig_name、sequence、reads、base_qualities）、AceAlignment（比对条目）和 AceData（顶层容器）等数据结构。支持 ACE 格式解析（ace_parse），自动识别 reads 区域（RD/AQ/QA/AF 关键字）和 contigs 区域（CT/CO/MQ）。支持单独解析 reads（ace_parse_reads）和 contigs（ace_parse_contigs）。支持序列化回 ACE 格式（ace_to_string），保持格式规范。支持查询函数：ace_get_contig 按名称获取 contig、ace_contig_length 获取 contig 长度、ace_contig_reads 获取 contig 中的 reads。支持计算每个位置的读取覆盖度（ace_read_coverage）、生成共有序列（ace_consensus_sequence）和计算 contig 的 GC 含量（ace_contig_gc_content）。适用于基因组组装项目的质控和分析。

### 133. 蛋白质组学工具 (Bio.SeqUtils.Proteomics)

实现蛋白质组学分析的核心功能，覆盖酶切、质量计算和碎片离子分析。支持 8 种蛋白酶：Trypsin（KR，排除 KP/RP）、Chymotrypsin（FYWML，排除 XP）、Pepsin（FL）、LysC（K，排除 KP）、ArgC（R，排除 RP）、CNBr（M）、GluC（DE）、AspN（D 前切割）。支持蛋白质酶切消化（proteomics_digest），可指定酶类型和漏切次数。支持酶切位点查找（proteomics_cleavage_sites），返回所有切割位置。支持肽段质量计算（proteomics_calculate_mass），提供单同位素和平均质量两种模式。支持同位素分布计算（proteomics_isotope_pattern），基于泊松模型生成同位素峰。支持碎片离子计算（proteomics_fragment_ions），生成 b 离子和 y 离子的质量列表。支持列出可用酶（proteomics_available_enzymes）和按名称选择酶（proteomics_select_enzyme）。提供便捷的胰蛋白酶消化函数（proteomics_trypsin_digest）。适用于蛋白质组学质谱分析和肽段鉴定。

### 134. PDB 片段映射 (Bio.PDB.FragmentMapper)

实现蛋白质结构片段的映射和分类功能，基于 DSSP 二级结构分配标准。支持 FragmentType 枚举（Helix、Sheet、Loop、Coil、Bridge、Turn）对片段类型进行分类。支持 Fragment 和 FragmentMapperResult 等数据结构存储片段信息和映射结果。支持从 DSSP 字符串分配片段（fm_assign_fragments），将二级结构字符串解析为结构化片段。支持按类型分类（fm_classify_secondary_structure），将 DSSP 字符（H/G/I→Helix、E→Sheet、L/C/ → Loop/Coil、B/b→Bridge、T/S→Turn）映射到 FragmentType。支持片段查询（fm_get_fragment）、序列提取（fm_fragment_sequence）和边界获取（fm_fragment_boundaries）。支持片段合并（fm_merge_fragments），将相邻同类型片段合并。支持按长度过滤（fm_filter_by_length），保留超过指定长度的片段。支持覆盖率分析（fm_coverage），计算每个残基的片段覆盖率。适用于蛋白质结构分析、二级结构预测和结构域映射。

### 135. 基因组图可视化 (Bio.Graphics.GenomeDiagram)

实现基因组图可视化的数据模型和 SVG 生成功能。支持 FeatureShape 枚举（Rectangle、Arrow、Diamond、CrossedArrow、Terminators）定义图形形状。支持 Diagram、Track、DiagramFeature、TrackFeature 和 DiagramStyle 等数据结构构建基因组图层级数据模型。支持创建图表（gd_create_diagram）并设置起止坐标。支持添加轨道（gd_add_track）和特征（gd_add_feature、gd_add_track_feature），可指定位置、链、标签、颜色和形状。支持样式设置（gd_set_style），包括循环/线性布局、颜色方案等。支持查询特征（gd_get_features）和查找重叠特征（gd_find_overlapping_features）。支持生成 SVG 输出（gd_to_svg_string），将数据模型转换为可渲染的 SVG 字符串。支持特征着色（gd_set_feature_color）和自动标注（gd_label_features），对大于指定阈值的特征添加标签。适用于基因组注释可视化、比较基因组学和教育演示。

### 136. 蛋白质内部坐标 (Bio.PDB.internal_coords)

实现蛋白质结构的内部坐标表示和转换功能。支持 TorsionAngle（扭矩角）、BondLength（键长）、BondAngle（键角）等基本结构单元。支持 InternalCoordResidue 和 InternalCoordChain 表示残基和链的内部坐标。支持计算二面角（ic_dihedral_angle），通过四个原子坐标计算扭矩角。支持主链扭矩角计算：phi（C(i-1)-N(i)-CA(i)-C(i)）、psi（N(i)-CA(i)-C(i)-N(i+1)）、omega（CA(i-1)-C(i-1)-N(i)-CA(i)）。支持从内部坐标到笛卡儿坐标的转换（ic_build_extended_chain），构建扩展多肽链。支持 Rotamer 和 RotamerLibraryEntry 表示侧链旋转异构体库，覆盖 GLU、VAL、PHE、GLY、ALA、SER 等常见氨基酸。支持弧度-角度转换（ic_rad_to_deg、ic_deg_to_rad）。支持 Ramachandran 区域数据，定义允许的 phi/psi 角度区域。适用于蛋白质结构预测、结构比较和构象分析。

### 137. 遗传算法 (Bio.GA)

实现遗传算法用于序列优化问题。支持 GAIndividual 表示个体（包含序列、适应度、代数、ID）和 GAPopulation 表示种群。支持适应度函数作为参数传入，用于评估每个个体的优劣。支持多种选择策略：锦标赛选择（ga_tournament_select）和轮盘赌选择（ga_roulette_select）。支持单点交叉（ga_single_point_crossover）和两点交叉（ga_two_point_crossover）生成后代。支持变异操作（ga_mutate），以一定概率随机替换序列中的碱基。支持精英保留（ga_evolve_generation），每代保留最优个体。支持 GAGenerationStats 统计每代的最优适应度、平均适应度和种群多样性。支持完整的进化循环（ga_evolve），从初始种群开始，经过多代选择、交叉、变异，直到达到收敛条件或最大代数。适用于序列设计、密码子优化和生物分子工程。

### 138. 染色体可视化 (Bio.Graphics.Chromosome)

实现染色体可视化的数据模型和 SVG 渲染功能。支持 ChrFeature（染色体特征）、ChrRegion（染色体区域）、Chromosome（染色体）等数据结构。支持 ChrFeatureType 枚举（Exon、Intron、Promoter、Enhancer、Marker、Other）注释特征类型。支持 ChrOrientation 枚举（Forward、Reverse、None）表示特征方向。支持 ChrBand 表示染色体带型（G带染色模式）。支持添加特征（add_feature）、区域（add_region）和带（add_band）到染色体。支持特征查询：按类型筛选（features_by_type）、区域内查找（features_in_region）。支持 ChrLayout 配置 SVG 布局参数（宽高、标签显示、着丝粒显示、带显示、染色体厚度）。支持生成线性染色体 SVG（chr_chromosome_to_svg）和圆形染色体 SVG（chr_circular_chromosome_to_svg）。支持带颜色映射（chr_band_color），根据染色强度返回 SVG 颜色。支持创建人类核型图（chr_create_human_karyotype）和细菌染色体图（chr_create_bacterial_chromosome）。适用于基因组注释可视化、细胞遗传学研究和教育演示。

### 139. impute 缺失值插补 (Bioconductor impute)

提供多种微阵列/RNA-seq 表达矩阵常用的缺失值插补方法。支持 KNN 插补参数（KNNImputeParam：k、by_row、eps、min_value、max_value）和 NA 统计汇总（ImputeNAStats）。实现按行均值填充（impute_by_row_mean）、按列中位数填充（impute_by_col_median）、KNN 近邻填充（impute_by_knn，基于欧氏距离加权）、LOCF（Last-Observation-Carried-Forward）、NOCB（Next-Observation-Carried-Backward）、以 0 填充（impute_na_by_zero）等方法。支持 NA 分布统计（impute_na_stats）和人类可读的 NA 汇总报告（impute_na_summary）。适用于基因表达数据分析中缺失值的预处理。

### 140. vsn 方差稳定化归一化 (Bioconductor vsn)

实现 vsn2 / glog（广义对数变换），用于芯片/测序计数的方差稳定化，解决 mean-variance 依赖问题。支持每列拟合的 glog 参数（VSNColParam：a、b）、vsn 拟合结果（VSNResult）和拟合超参数（VSNControl）。核心函数包括：glog 变换（asinh 形式）、glog 反函数、vsn2 一步完成拟合+变换、vsn_fit_and_report 仅拟合不做变换、vsn_denoise 去噪（vsn2 + glog_inv）、mean_sd_bins 按均值分位数分桶计算 mean-SD 曲线、mean_sd_ascii 绘制 ASCII mean-SD 诊断图。适用于基因表达数据的方差稳定化预处理。

### 141. GSEABase 基因集管理 (Bioconductor GSEABase)

实现 GeneSet / GeneSetCollection 数据结构，支持 GMT / GMX 格式的解析与写出，并提供基因集集合运算。支持 GmtGeneSet（单条基因集：name、description、gene_ids、collection_type、organism、id）和 GmtGeneSetCollection（基因集集合）。集合类型枚举涵盖 GO_BP/MF/CC、KEGG、Reactome、Hallmark、Canonical 等。支持集合运算：overlap、Jaccard、union、intersect、setdiff。支持 Collection 方法：by_name 查询、at 索引查询、all_sizes 大小列表、summary 摘要。I/O 支持：parse_gmt、parse_gmx、write_gmt。适用于功能富集分析和基因集研究。

### 142. PCAtools 高级 PCA 分析 (Bioconductor PCAtools)

基于幂迭代的 eigendecomposition，提供完整的 PCA 分析流水线。支持 FullPCAResult（scores、loadings、eigenvalues、variance_explained、cumulative_variance）和 PCAOutlierResult。核心功能包括：pcatools_run_pca 主入口（支持 z-score 标准化）、pcatools_summary 摘要、scree_plot_ascii 碎石图、pca_biplot_ascii biplot 可视化、find_pca_outliers 基于马氏距离的异常点检测、variable_correlations 变量-主成分 Pearson 相关系数矩阵。适用于基因表达数据的降维和可视化分析。

### 143. Bio.Data 生物数据常量 (Biopython Bio.Data)

提供 IUPAC 碱基映射、氨基酸缩写映射、标准密码子表及反向互补等分子生物学常用数据常量。支持 IUPAC 模糊碱基映射（DNA/RNA 各 11 种）、氨基酸三字母↔单字母映射（21 种）、氨基酸理化性质查询（疏水性、电荷、等电点）。实现标准密码子表（64 密码子）及反向密码子表、同义密码子数量统计、IUPAC 简并碱基反向互补。适用于序列分析和分子生物学研究的数据查询。

### 144. Bio.Seq.Approximate 近似字符串匹配 (Biopython Bio.Seq.Approximate)

实现基于动态规划的近似字符串匹配算法，支持错配、插入和缺失（indels）。支持 ApproxMatch（pattern、query、start、end、mismatches、score）和 ApproxWordResult 数据结构。核心函数包括：count_mismatches 等长错配计数、approx_search 无 indels 近似搜索、approx_search_with_indels 允许 indels 的近似搜索、levenshtein_distance 编辑距离、approx_find_all 查找所有近似匹配、approx_word_search 单词近似搜索、approx_best_match 基于错误率的最佳匹配。适用于 CRISPR off-target 检测和序列相似性搜索。

### 145. Bio.Pairwise2 灵活双序列比对 (Biopython Bio.Pairwise2)

实现 Needleman-Wunsch（全局）和 Smith-Waterman（局部）比对算法，支持自定义匹配/错配评分矩阵和开放/延伸空位罚分。支持 PairwiseAlignResult（aligned_seq1、aligned_seq2、score、mode、start/end）和 PairwiseMode 枚举。核心函数包括：pairwise_globalxx/pairwise_localxx 自定义参数比对、pairwise_global/pairwise_local 便捷比对、pairwise_globalms/pairwise_localms 基于替换矩阵的比对、simple_score/identity_score/matrix_score 评分函数创建、dna_matrix DNA 配对评分矩阵生成、alignment_summary 比对摘要（Identity/Similarity）。适用于序列比对和进化分析。

### 146. Bio.Compound 化合物数据结构 (Biopython Bio.Compound)

实现化合物、化学反应和代谢通路的数据结构，支持分子式解析、分子量计算以及预设模板。支持 Compound（id、name、formula、charge、smiles、aliases、pathways、reactions）、ChemicalReaction 和 CompoundPathwayMap 数据结构。核心函数包括：parse_formula 分子式解析、compound_molecular_weight 分子量计算、预设模板（葡萄糖、果糖、丙酮酸、乙酸、柠檬酸）。支持化合物属性设置（formula、charge、smiles）、别名添加、通路/反应关联。适用于代谢组学和生物化学研究。

### 147. EnhancedVolcano 火山图可视化 (Bioconductor EnhancedVolcano)

实现差异表达分析结果的火山图可视化功能，支持基因分类（上调/下调/不显著）、p 值与 fold-change 阈值设定以及 ASCII 渲染。支持 VolcanoClassification（Up/Down/NonSig）、VolcanoGene 和 VolcanoResult 数据结构。核心函数包括：volcano_plot 创建火山图结果、volcano_plot_default 使用默认参数、neg_log10_p 计算 -log10(p-value) 并处理 p=0、volcano_sample 创建示例数据集。VolcanoResult 方法包括：get_n_genes/get_n_up/get_n_down/get_n_nonsig 统计、get_up_genes/get_down_genes 基因列表、to_ascii ASCII 渲染、summary 文本摘要。适用于差异表达分析结果的可视化展示。

### 148. ReportingTools 报告生成 (Bioconductor ReportingTools)

实现生物信息学分析报告生成功能，支持文本、表格和图形的混合报告，以及 ASCII 表格渲染。支持 ReportSectionType 枚举（Text/Table/Plot/Section/Header）、ReportSection、ReportColumn、ReportTable 和 ReportDocument 数据结构。ReportDocument 方法包括：创建、set_author、add_text、add_table、add_plot、get_title/get_author/get_n_sections/get_n_tables、render 渲染为文本。ReportTable 方法包括：创建、from_columns、get_column_names/get_n_rows/get_n_columns、to_ascii ASCII 渲染。适用于生物信息学分析报告的自动化生成。

### 149. karyoploteR 核型可视化 (Bioconductor karyoploteR)

实现核型可视化功能，支持染色体轨道、数据点（SNP、GC 含量等）以及 ASCII 渲染。支持 TrackType 枚举（Points/Lines/Bars/Heatmap/Ideogram）、KaryotypeRegion、TrackPoint、KaryotypeTrack、IdeogramBand 和 KaryotypePlot 数据结构。KaryotypePlot 方法包括：创建（支持 hg38/hg19）、get_genome/get_chromosomes/get_chromosome_size、add_track/add_region、to_ascii ASCII 渲染指定染色体。KaryotypeTrack 方法包括：创建、add_point、set_color、set_y_range、get_data/get_n_points、filter_chromosome。适用于基因组数据的染色体级可视化。

### 150. SystemPipeR 流水线编排 (Bioconductor SystemPipeR)

实现生物信息学分析流水线编排功能，支持步骤管理、依赖关系、进度追踪以及 ASCII 可视化。支持 StepStatus 枚举（Pending/Running/Completed/Failed/Skipped）、PipelineStep、Pipeline 和 PipelineConfig 数据结构。Pipeline 方法包括：创建、set_description/set_param/get_param、add_step/get_step、get_id/get_name/get_n_steps/get_steps、get_completed_count/get_failed_count/get_pending_count/get_skipped_count、get_progress、can_run_step、summary、to_ascii。PipelineStep 方法包括：创建、set_description/set_args/add_dependency/add_input/add_output、get_id/get_name/get_command/get_status/get_duration/get_dependencies。PipelineConfig 方法包括：创建、set_cores、get_work_dir/get_input_dir/get_output_dir/get_cores。pipeline_sample 创建示例 RNA-seq 流水线。适用于生物信息学分析工作流的编排与管理。

### 151. muscat 单细胞差异状态分析 (Bioconductor muscat)

实现单细胞 RNA-seq 差异状态（Differential State, DS）分析功能，参考 Bioconductor muscat 包（Crowell HL et al. 2020 Nature Communications），支持伪批量聚合、统计检验和样本级 QC。支持 AggregationMethod 枚举（Sum/Mean/Median 三种聚合方法），辅助构造函数 aggregation_sum()/aggregation_mean()/aggregation_median()。支持 DSMethod 枚举（EdgeR/DESeq2/Limma 三种检验方法），辅助构造函数 ds_method_edger()/ds_method_deseq2()/ds_method_limma()。核心数据结构：SingleCell（cell_id/sample_id/cluster_id/group_id/gene_counts 单细胞计数）、PseudoBulk（sample_id/cluster_id/group_id/n_cells/gene_counts 伪批量样本）、SampleQC（sample_id/cluster_id/n_cells/n_genes/total_counts/median_genes_per_cell QC统计）、DSResult（gene_id/cluster_id/log2fc/p_val/p_adj/mean_ctrl/mean_stim/significant 单基因检验结果）、DSResults（全局结果容器）。核心函数：aggregate_cells() 按 (sample_id, cluster_id) 分组聚合单细胞为伪批量，支持 Sum/Mean/Median 三种统计量；run_ds_analysis() 按聚类分别拆分对照组（ctrl/control）与处理组，计算 log2 fold change、Welch t 检验近似 p 值、全局 Benjamini-Hochberg FDR 校正、根据 fdr_threshold 和 log2fc_threshold 判定显著性；compute_qc() 计算伪批量样本的细胞数、基因数、总计数等 QC 指标；muscat_sample_data() 生成 4 样本 × 2 聚类 × 5 细胞 × 5 基因的示例数据集。DSResults 方法：get_n_results()/get_significant()/get_cluster_results(cluster_id)/get_top_genes(n)/summary()。PseudoBulk 方法：get_count()/total_counts()/n_expressed()。SingleCell 方法：set_count()/get_count()/total_counts()/n_expressed()。适用于单细胞转录组的亚群特异性差异表达分析、伪批量差异状态检验和实验 QC 评估。

### 152. infercnv 单细胞拷贝数变异推断 (Bioconductor infercnv / infercnvpy)

实现单细胞 RNA-seq 的拷贝数变异（CNV）推断功能，参考 Broad Institute inferCNV（Tirosh I et al. 2016 Science）与 infercnvpy 实现，用于识别肿瘤细胞的大片段染色体增减并分离恶性/正常细胞。支持 ReferenceMethod 枚举（GlobalMean / ReferenceCategories / Custom(ref_profile) 三种参考策略），辅助构造函数 ref_method_global_mean() / ref_method_reference_categories() / ref_method_custom(ref_profile)。核心数据结构：GenePosition（gene_id/chromosome/start/end 基因位置）、OrderedGenes（按染色体自然顺序+起点排序的基因列表）、CellAnnotation（cell_id/category 细胞注释）、InferCNVInput（n_cells/n_genes/expression[cell x gene]/annotations/ordered_genes 输入数据）、CNVParams（window_size/lfc_cap/noise_threshold/reference_method/reference_categories 参数集）、CNVResult（cnv_matrix/cell_ids/gene_ids/chromosome_boundaries/per_cell_cnv_score/per_cluster_cnv_score/params 结果）。核心算法流水线 run_infercnv()：① 染色体位置排序 + 排除 chrM/chrX/chrY 等噪声染色体；② 参考表达谱构造：GlobalMean 全体平均、ReferenceCategories 指定正常细胞类型均值、Custom 用户向量；③ 每细胞计算 log fold change（单参考=减法；多参考=有界计算：落在[min,max]区间→0、超过max→cell-max、低于min→cell-min）；④ ±lfc_cap 截断去除极值；⑤ 按染色体独立的金字塔权重滑动窗口平滑（window_size 默认 100）；⑥ 每细胞中位数中心化消除全局偏移；⑦ 噪声阈值（|x|<noise→0）去背景；⑧ 输出 per_cell_cnv_score（行平均|CNV|）与 per_cluster_cnv_score 及恶性细胞预测 predict_tumour_cells(normal_cat, factor=1.5)。辅助函数：ordered_genes() 自然染色体排序（1..22<X<Y<M）、exclude_chromosomes() 过滤、log_normalize_counts() 原始计数→log2(TP100K+1)、chromosome_boundaries() 染色体边界、infercnv_sample_data(n_tumour,n_normal,n_chr,n_genes_per_chr,seed) 合成数据集（chr1 前半段缺失、chr2 后半段扩增）。CNVResult 方法：cluster_score() / chromosome_list() / summary() / at(i,j) / predict_tumour_cells()。适用于单细胞肿瘤微环境分析中的恶性细胞识别、克隆异质性解析、CNV 负荷评估。

### 153. SCENIC 单细胞调控网络推断与聚类 (Bioconductor SCENIC)

实现单细胞 RNA-seq 的基因调控网络推断与细胞状态识别功能，参考 Aibar et al. 2017 Nature Methods 的 SCENIC 三步流水线。核心数据结构：CoExpressionModule（tf_name/targets/weights TF-靶基因共表达模块）、ScenicRegulon（tf_name/targets/weights/n_targets 经剪枝的调控子）、SCENICInput（n_genes/n_cells/expression[gene x cell]/gene_names/cell_names/tf_names 输入数据）、SCENICResult（regulons/auc_matrix/binary_matrix/thresholds/cell_states/master_regulators 结果）、BinarizeMethod 枚举（MeanStd/KMeans2/Median 三种二值化方法），辅助构造函数 binarize_mean_std()/binarize_kmeans2()/binarize_median()。三步算法流水线 run_scenic()：① Step 1 GRN 推断 build_coexpression_modules()：对每个 TF 计算其与所有基因的绝对 Pearson 相关系数作为重要性权重，按权重降序选取 Top-K 靶基因构建共表达模块（min_targets=10, top_k=50）；② Step 2 Regulon 构建 build_regulons()：按权重分位数剪枝共表达模块（默认保留上 75% 权重靶基因），过滤靶基因数 < min_targets 的模块；备选 prune_by_motif_ranking() 提供 cisTarget 风格的 motif 排名剪枝（保留 motif 排名前 5% 的靶基因）；③ Step 3 AUCell 活性评分 compute_regulon_activity()：对每个细胞构建基因表达排名（降序），对每个 regulon 计算其靶基因在细胞排名前 5% 位置的 recovery curve AUC（归一化到 [0,1]），输出 regulon x cell 的 AUC 矩阵；④ Step 4 二值化与细胞状态 binarize_activity() 对每个 regulon 独立计算阈值（MeanStd=均值+0.5σ、KMeans2=1D k-means 双簇中点、Median=中位数），assign_cell_states() 按最高活性 regulon 分配细胞到簇并识别主控调控因子。辅助函数：build_cell_rankings() 构建每细胞基因排名、scenic_abs_correlation() 向量绝对相关系数、scenic_sample_data(n_genes,n_cells,n_tfs,seed) 合成数据集（TF1 活跃于前半细胞、TF2 活跃于后半细胞）。SCENICResult 方法：auc_at(r,c)/binary_at(r,c)/cell_state(c)/summary()/regulon_targets(r)/top_regulons(n=5)。适用于单细胞转录因子调控网络重建、细胞状态分类、主控调控因子识别和调控异质性分析。

### 154. MSstats 蛋白质显著性分析 (Bioconductor MSstats)

实现质谱定量蛋白质组学的蛋白质显著性分析功能，支持数据归一化、肽段到蛋白质汇总、组间差异比较以及样本量设计。支持 MSType 枚举（DDA/DIA/SRM/TMT）、MSNormalization 枚举（None/Median/Quantile/GlobalStandards）、MSSummarization 枚举（Tukey/Linear/LogSum）数据结构。核心函数包括：ms_data_process 数据处理（过滤/Log2转换/归一化）、ms_summarize 肽段汇总到蛋白质（Tukey 中位数平滑/线性模型/LogSum）、ms_group_comparison 组间差异比较（t 检验/BH-FDR 校正）、ms_design_sample_size 样本量计算、msstats_sample_data 示例数据集。MSGroupComparison 方法包括：get_n_results/get_n_significant/get_significant/get_top_proteins/summary。适用于定量质谱数据的蛋白质水平差异丰度分析。

### 155. NOISeq 噪声鲁棒差异表达 (Bioconductor NOISeq)

实现噪声鲁棒的 RNA-seq 差异表达分析功能，支持多种归一化方法和基于噪声分布的概率计算。支持 NOISeqNorm 枚举（RPKM/TMM/UpperQuartile/None）、NOISeqMethod 枚举（NOISeqBio/NOISeqSim）数据结构。核心函数包括：noiseq_normalize 数据归一化（RPKM/TMM/上四分位）、noiseq_run 差异表达分析（M 值/D 统计量/概率计算/显著性排序）、noiseq_qc 质控诊断、noiseq_sample_data 示例数据集。NOISeqSample 包含 sample_id/condition/counts 字段，方法包括 set_count/get_count/library_size/n_expressed。NOISeqResult 包含 gene_id/mean_control/mean_treatment/log2fc/divergence/prob/significant/ranking 字段。NOISeqResults 方法包括：get_top_genes/get_significant/get_up_regulated/get_down_regulated/summary。适用于数据质量感知的 RNA-seq 差异表达分析。

### 156. Gviz 基因组可视化轨道 (Bioconductor Gviz)

实现基因组可视化轨道系统，支持多种轨道类型和 ASCII 渲染。支持 GvizTrackType 枚举（AnnotationTrack/GeneRegionTrack/DataTrack/IdeogramTrack/GenomeAxisTrack/SequenceTrack）、GvizStrand 枚举（Forward/Reverse/Unstranded）数据结构。核心函数包括：gviz_feature 创建基因组特征、gviz_track 创建轨道、gviz_region 创建基因组区域、gviz_plot 创建绘图、gviz_sample_plot 创建示例绘图。GvizTrack 方法包括：set_color/set_label/add_feature/add_data_point/get_n_features/get_n_data_points/get_features_in_region/get_data_in_region。GvizPlot 方法包括：add_track/get_n_tracks/get_region/to_ascii/summary。ASCII 渲染支持坐标轴、特征块、数据点和核型显示。适用于基因组注释和数据的沿坐标可视化。

### 157. SeqLocation 序列位点类型系统 (Bio.SeqFeature)

实现序列位置和位点的完整类型系统，用于基因组特征注释的精确位置表示。支持 5 种位置类型：ExactPosition（精确位置）、BeforePosition（之前位置，如 <100）、AfterPosition（之后位置，如 >100）、OneOfPosition（多选位置，如 100^110）、WithinPosition（内部位置，如 (100)^(110)）。支持 2 种位点类型：SimpleLocation（简单位点，包含起止位置、链方向、类型）和 CompoundLocation（复合位点，包含多个子位点数组，用于外显子/内含子等断裂特征）。提供位点操作方法：start() 获取 0-based 起始位置、end() 获取 0-based 结束位置（不含）、strand() 获取链方向、reverse() 翻转链、is_overlapping() 检查重叠。支持位置字符串的解析与格式化（如 "100..110"、"complement(100..110)"、"join(100..200, 300..400)"）。适用于 GenBank/EMBL 等格式的序列特征注释解析与生成。

### 158. BioReference 文献引用管理 (Bio.Reference / Bio.Medline)

实现文献引用的完整管理功能，用于序列记录的参考文献追踪。支持 BioReference 结构体存储引用信息：标题(title)、作者(authors)、期刊(journal)、年份(year)、PubMed ID(pubmed_id)、DOI(doi)、引用类型(type)、引用位置(locations)、评论(comment)、数据页码(data_pgr)。提供引用操作方法：format() 按 APA 格式格式化输出、get_pubmed_id() 获取 PubMed ID、get_doi() 获取 DOI、get_authors() 获取作者列表。支持从 Medline 记录解析引用信息（parse_medline_record），自动提取标题、作者、期刊、年份、PubMed ID 等字段。支持位置信息存储（如 (1, 100) 表示引用序列的 1-100 位）。适用于 GenBank/EMBL 记录的参考文献注释和文献信息检索。

### 159. ProtDao 蛋白质无序区域预测 (Bio.SeqUtils.ProtDao)

实现基于 IUPred 算法的蛋白质无序区域预测功能，用于识别天然无序蛋白质区域（IDP）。支持 20 种氨基酸的无序度评分（disorder_score）和能量评分（energy_score），采用 IUPred 打分矩阵。提供 DisorderResult 结构体存储预测结果：序列(sequence)、得分数组(scores)、阈值(threshold_disordered/threshold_long)。DisorderResult 方法包括：get_scores() 获取所有位置的无序度得分、get_regions() 获取无序区域列表、get_n_regions() 获取区域数量、get_n_disordered() 获取无序残基数量、get_fraction_disordered() 获取无序残基比例、get_longest_region() 获取最长无序区域、disordered_sequence() 获取无序序列视图（* 表示无序，- 表示有序）、summary() 生成摘要报告、to_ascii() 生成 ASCII 可视化图。支持 DisorderRegion 结构体存储区域信息：起始位置(start)、结束位置(end)、长度(length)、平均得分(avg_score)、区域类型(region_type: disordered/highly disordered)。提供完整预测流程：prot_dao_predict() 一键预测、prot_dao_sample_sequence() 获取示例序列。适用于蛋白质结构预测、功能注释和天然无序区域分析。

### 160. featureCounts 基因特征 read 计数 (Bioconductor Rsubread/featureCounts)

实现 RNA-seq/DNA-seq read 在基因组特征（exon/gene/transcript）上的计数功能，参考 Bioconductor Rsubread 包的 featureCounts 算法。提供 FeatureAnnotation 结构体存储特征注释：染色体(chr)、起止位置(start/end_)、链(strand)、基因ID(gene_id)、转录本ID(transcript_id)、特征类型(feature_type)、唯一标识(feature_id)，支持 length() 计算特征长度、overlaps() 检测区间重叠、overlap_length() 计算重叠长度。提供 ReadAlignment 结构体存储 read 比对信息：read_id、染色体、起止位置、链、比对质量(mapq)、多比对数(n_alignments)、配对末端信息(is_paired/mate_chr/mate_start/fragment_length)，支持 new_paired() 创建配对末端 read、fragment_span() 获取 fragment 跨度。提供 FeatureCountsConfig 配置：最小重叠(min_overlap)、最小比对质量(min_mapq)、链特异性模式(strand_mode: Unstranded/Stranded/Reversed)、多比对计数(multi_count/frac_multi_count)、配对末端 fragment 计数(count_fragments)。核心函数 feature_counts_count() 执行计数：遍历 read、过滤低质量、查找重叠特征、处理歧义(ambiguous)、链匹配检查、分数计数。结果 FeatureCountsResult 提供：get_count() 按特征/样本获取计数、library_size() 计算 library size、cpm() 计算 CPM 归一化、get_feature_total() 跨样本特征总数、summary() 生成报告。提供链模式辅助函数 strand_unstranded()/strand_stranded()/strand_reversed() 和示例数据 feature_counts_sample_data()。适用于 RNA-seq 基因表达定量、外显子计数和多样本计数矩阵生成。

### 161. DRIMSeq 差异转录本使用分析 (Bioconductor DRIMSeq)

实现基于 Dirichlet-multinomial 模型的差异转录本使用（DTU）分析，参考 Bioconductor DRIMSeq 包。提供 TranscriptCount 结构体存储转录本计数：转录本ID(transcript_id)、基因ID(gene_id)、样本ID(sample_id)、条件(condition)、count、gene_count，支持 proportion() 计算转录本占比。提供 DRIMSeqConfig 配置：最小count(min_count)、最小比例(min_proportion)、收敛容差(tolerance)、最大迭代(max_iter)、显著性水平(alpha)、归一化方法(norm_method: None/TMM/Sum)。核心分析流程 drimseq_test_differential()：drimseq_aggregate_by_gene() 按基因/样本聚合计数、drimseq_compute_proportions() 计算转录本比例、drimseq_filter_counts() 过滤低表达转录本、drimseq_wald_test() 执行 Wald 检验（计算 delta、使用方差归一化统计量、卡方分布 p 值）、benjamini_hochberg_correct() 多重检验校正。结果 DRIMSeqResult 提供：get_significant() 获取显著基因、get_top_genes() 获取 top 基因、summary() 生成报告。包含 Dirichlet-multinomial 对数似然计算、Wilson-Hilferty 卡方 p 值近似、正态分布生存函数、收敛诊断。提供 DRIMSeqGeneResult 存储每基因结果：gene_id、p_value、adj_p_value、statistic、converged。提供归一化方法辅助函数 drimseq_norm_none()/drimseq_norm_tmm()/drimseq_norm_sum() 和示例数据 drimseq_sample_data()。适用于 RNA-seq 差异转录本使用分析、可变剪接研究。

### 162. RaggedExperiment 参差突变数据结构 (Bioconductor RaggedExperiment)

实现癌症基因组学中常见的参差（ragged）突变数据结构，参考 Bioconductor RaggedExperiment 包，用于存储基因×样本的稀疏突变矩阵。提供 MutationType 枚举覆盖 11 种突变类型：Missense_Mutation、Nonsense_Mutation、Frame_Shift_Ins/Del、In_Frame_Ins/Del、Splice_Site、Translation_Start_Site/Stop_Site、Multi_Hit_Mutation、Silent_Mutation、Other，支持 to_string() 转换和辅助构造函数 mutation_missense()/mutation_nonsense() 等。提供 MutationRecord 结构体存储突变记录：sample_id、gene_symbol、chrom、pos、ref_allele、alt_allele、mutation_type，支持 mutation_id() 生成唯一标识、is_loss_of_function() 识别 LoF 突变（Nonsense/Frame_Shift）。核心结构 RaggedExperiment 使用三层嵌套数组 data[row][col][records] 存储每个基因×样本的突变记录列表，维护 rownames(基因)、colnames(样本)、counts(突变计数缓存)、sample_tmb(样本 TMB)。方法包括：add_record() 添加突变、n_rows()/n_cols() 维度、get_records() 按基因/样本查询、get_gene_records()/get_sample_records() 按行/列查询、get_tmb() 计算 TMB、get_tmb_per_mb() 计算每 Mb TMB、get_count_matrix() 获取突变计数矩阵、genes_mutated_per_sample() 每样本突变基因数、filter_by_genes()/filter_by_samples()/filter_by_type() 按基因/样本/突变类型过滤、summary() 生成报告。提供 ragged_sample_data() 示例数据。适用于癌症基因组突变数据分析、TMB 计算、突变谱分析。

### 163. PairwiseAligner 统一双序列比对 (Biopython Bio.Align.PairwiseAligner)

实现统一的双序列比对接口，参考 Biopython Bio.Align.PairwiseAligner，整合全局（Needleman-Wunsch）和局部（Smith-Waterman）两种模式。提供 AlignmentMode 枚举：Global（全局对齐，Needleman-Wunsch 回溯到 (0,0)）、Local（局部对齐，Smith-Waterman 零截断 + 最大分值回溯），辅助构造函数 pairaligner_global()/pairaligner_local()。提供 SubstitutionMatrixChoice 枚举：NoSubstMatrix（使用 match_score/mismatch_score）、Blosum62（蛋白质 20 标准氨基酸 BLOSUM62 替换矩阵），辅助构造函数 pairaligner_no_matrix()/pairaligner_blosum62()。核心结构 PairwiseAlignerConfig：11 字段（mode、match_score、mismatch_score、gap_open、gap_extend、target_gap_open/extend、query_gap_open/extend、submatrix、alphabet_type）；default_dna()：Global + DNA + match=1/mismatch=-1/gap_open=-1/gap_extend=-0.5；default_protein()：Global + PROTEIN + BLOSUM62 + gap_open=-10/gap_extend=-0.5；7 个 setter 全部手动列出字段实现不可变更新。核心算法 pairaligner_align() 使用 Gotoh 三矩阵仿射方案：M（匹配/错配）、Ix（target 方向空位）、Iy（query 方向空位）；支持独立的 target/query 空位 open/extend 参数；回溯时构建 match_line（"|" 匹配、"." 错配、" " 空位）。输出结构 PairwiseAlignment：aligned1()/aligned2()、score()、alignment_length()、identities()/identity_pct()、gaps_count()、target_start/end、query_start/end。提供 pairaligner_sample_data() 返回经典蛋白质测试对 HEAGAWGHEE / PAWHEAE。适用于 DNA/蛋白质双序列比对、全局同源性检测、局部 motif 查找、空位模型比较。

### 164. NaiveBayes 生物序列分类器 (Biopython Bio.NaiveBayes)

实现基于 k-mer 频率特征的朴素贝叶斯序列分类器，参考 Biopython Bio.NaiveBayes，用于 DNA/RNA/蛋白质序列的类别预测（AT/GC 富集、蛋白质家族、物种分类等）。结构 NBClassModel：class_label、n_sequences、total_kmers、kmer_counts(Map)、class_prior。主结构 NaiveBayesClassifier：kmer_size(Int, 默认 3)、alpha(Double, 默认 1.0, Laplace 平滑系数)、vocabulary(Array[String], 所有见过的 k-mer 词汇表)、class_labels(Array[String])、models(Map[String, NBClassModel])。::new() 默认 k=3, alpha=1.0；set_kmer_size(val)、set_alpha(val) 手动列出 5 字段不可变更新。核心流程：naive_bayes_extract_kmers(sequence, k) 滑动窗口 k-mer 计数；naive_bayes_train(classifier, sequences, labels) 统计每类 k-mer 计数、总数、先验概率、词汇表；naive_bayes_predict_log_probs() 计算所有类的对数后验 ∝ log(P(class)) + Σ log(P(kmer|class))；P(kmer|class) = (count_kmer + alpha) / (total_kmers + alpha * vocab_size)（Laplace 平滑），未见过的 k-mer 使用默认平滑概率；naive_bayes_predict() 返回 (best_label, best_log_prob) argmax；naive_bayes_predict_proba() 使用 softmax(exp(log_prob - max)) 归一化为概率，总和约为 1.0；naive_bayes_top_k() 返回按概率降序前 k 类；naive_bayes_accuracy() 计算分类准确率。提供 naive_bayes_sample_data() 返回 8 条 DNA：4 条 AT_rich + 4 条 GC_rich，便于示例和快速测试。适用于 DNA/蛋白质序列分类、AT/GC 富集预测、序列家族识别、样本属性判别。

### 165. Markov 高阶马尔可夫链建模 (Biopython Bio.Markov)

实现 1~3 阶高阶马尔可夫链建模，参考 Biopython Bio.Markov，用于生物序列的概率打分、CpG 岛检测、合成序列生成、稳态分析等。枚举 ChainType：FirstOrder / SecondOrder / ThirdOrder，辅助构造函数 markov_first_order()/markov_second_order()/markov_third_order()。结构 MarkovModel：order(Int)、chain_type(ChainType)、states(Array[String])、transition_counts(Map[context, Map[next, count]])、transition_probs(Map[context, Map[next, prob]])、initial_probs(Map[start_kmer, prob])、pseudo_count(Double)。::new() 默认 order=1, states=["A","C","G","T"], pseudo=1.0；default_dna() 相同；default_protein() states=20 氨基酸。setter：set_order(val)、set_states(states~)、set_chain_type(val~)、set_pseudo_count(val)，全部手动列出 7 字段不可变更新。训练 markov_build_model(sequences, order?, states?, pseudo_count?)：对长度 N 的序列从 i=order 到 N-1 扫描 context（seq[i-order:i]）→ next_char，累计转移计数并归一化；初始 k-mer 计数归一化为 initial_probs；states 空自动从数据推断；pseudo_count 用于后续未知 context 平滑。打分 markov_score_sequence()：初始 log P(start_kmer) + Σ_i log P(next|context)；markov_score_per_base() 返回每位置 log 概率（长度 = len - order）。未知转移概率平滑：完全未见过的 context → P = pseudo / (1 + pseudo * N)（高 pseudo 更接近均匀，因此对未知打分更高）；已知 context 但缺少 next_char → 汇总该 row 的 total，P = pseudo / (row_total + pseudo * N)。序列生成 markov_generate_sequence(model, length, seed?)：根据 initial_probs 加权采样起始 context，然后按 transition_probs 加权采样后续字符；seed 驱动简单 LCG RNG (rand = (x*1103515245+12345) & 0x7fffffff)。CpG 岛检测 markov_log_odds(cpg_model, bg_model, sequence)：log P_CpG - log P_bg，正值表示符合 CpG 统计特征。稳态 markov_stationary_distribution()：从均匀分布开始迭代 M100 次或 δ<1e-8 收敛。适用于 CpG 岛检测、DNA 组成分析、序列似然度打分、异常序列检测、合成序列生成。

### 166. HTSFilter RNA-seq count过滤 (Bioconductor HTSFilter)

实现基于 CPM (Count Per Million) 阈值的 RNA-seq count 过滤，参考 Bioconductor HTSFilter 包。核心函数 hts_filter_single_cpm(count, library_size) 将单基因 count 转换为 CPM 值：count / library_size × 1,000,000。hts_filter_library_sizes(counts) 计算每个样本的总文库大小（counts 行求和）。hts_filter_cpm(counts, library_sizes) 应用于整个计数矩阵，生成 CPM 矩阵。hts_filter(counts, groups, cpm_threshold, min_samples_per_group) 主过滤函数：对每个基因计算所有样本的 CPM，然后按组检查是否至少 min_samples_per_group 个样本 CPM ≥ 阈值，若任一组满足条件则保留该基因。返回 HTSFilterResult 结构：cpm_threshold、min_samples_per_group、groups、cpm_matrix、keep 掩码、n_genes_input/kept/removed。辅助函数 hts_filter_apply() 根据 keep 掩码提取过滤后的计数矩阵，hts_filter_apply_names() 同步提取基因名，hts_filter_retention_rate() 计算保留率，hts_filter_summary() 生成中文摘要。hts_filter_sample_data() 提供 20 基因 × 6 样本 (3 control + 3 treatment) 的示例数据，便于测试和演示。适用于 RNA-seq 预处理中去除低表达/噪声基因，为下游差异表达分析提供可靠输入。

### 167. baySeq 贝叶斯差异表达分析 (Bioconductor baySeq)

实现基于负二项模型的贝叶斯差异表达分析，参考 Bioconductor baySeq 包。核心结构 BayesGeneResult 存储每个基因的分析结果：gene_index、gene_name、log_fold_change、dispersion、posterior_prob_de、map_expression_a/map_expression_b、is_de。BayesResult 汇总所有基因结果及全局参数。bayseq_lgamma(x) 使用 Lanczos 近似计算 log-gamma 函数（g=5），支持 (0,1) 反射公式处理小值。bayseq_estimate_dispersion_gene(counts_a, counts_b) 通过基因级别估计离散度 φ：均值 ÷ 方差。bayseq_nb_log_likelihood(counts, mu, phi) 计算负二项分布的对数似然。bayseq_estimate_prior(counts, groups) 估计 Gamma 先验的 shape 和 rate 参数。bayseq_test(counts, groups, gene_names, dispersion_shape?, dispersion_rate?, alpha?) 主函数：对每个基因计算两组的 MAP 表达估计（含 Bayes 正则化项），计算对数似然比 (log_lik_ratio)，近似贝叶斯因子 → 后验概率 P(DE|data)，并根据 alpha 阈值判定显著 DE 基因。辅助函数 bayseq_get_top_de() 按 |logFC| 排序返回前 N 个 DE 基因，bayseq_get_de_genes() 提取所有显著 DE 基因，bayseq_summary() 生成分析摘要。bayseq_sample_data() 提供 20 基因 × 6 样本示例数据。适用于 RNA-seq 差异表达分析，提供基于贝叶斯框架的稳健 DE 判定。

### 168. CellChat 细胞间通讯分析 (Bioconductor CellChat)

实现细胞间通讯分析，参考 Bioconductor CellChat 包。核心结构 InteractionScore 存储每个配体-受体互作的评分：ligand、receptor、source_celltype、target_celltype、score、perm_mean/perm_std、p_value、p_adj、significant。CellChatResult 汇总所有互作评分、细胞类型、基因名及统计参数。LRPair 结构表示配体-受体对（ligand、receptor）。cellchat_lr_database() 返回内置 15 对配体-受体互作数据库（如 TNF→TNFR、VEGF→VEGFR、EGF→EGFR 等）。cellchat_mean_expr(expression, cell_types, celltype, gene_idx) 计算指定细胞类型在某基因上的平均表达。cellchat_analyze(expression, cell_types, gene_names, lr_pairs, n_permutations?, seed?, fdr?) 主分析流程：1) 对每个 LR 对，在所有 source→target 细胞类型组合中计算互作评分（平均配体 × 平均受体）；2) 置换检验：随机打乱细胞类型标签 N 次，生成零分布；3) 计算 p-value = (count_permuted_higher + 1) / (n_perm + 1)；4) BH-FDR 校正；5) 输出 CellChatResult。辅助函数 cellchat_get_top() 按评分排序返回前 N 个互作，cellchat_get_significant() 提取显著互作，cellchat_aggregate() 聚合每个细胞类型对的总通信强度，cellchat_summary() 生成分析摘要。cellchat_sample_data() 提供 30 细胞 × 15 基因的示例数据（3 种细胞类型）。适用于单细胞数据中细胞间通讯信号的探索和可视化。

### 168. MutationalPatterns 体细胞突变谱分析 (Bioconductor MutationalPatterns)

实现体细胞突变模式分析，参考 Bioconductor MutationalPatterns 包。核心枚举 MutationCategory 定义 6 种嘧啶基准突变类型（C>A、C>G、C>T、T>A、T>C、T>G 转换/颠换），并提供构造辅助函数 mutation_c_a() / mutation_c_g() / mutation_c_t() / mutation_t_a() / mutation_t_c() / mutation_t_g() 以支持跨包构造。SomaticMutation 结构记录每条体细胞突变：chromosome、position、ref_base、alt_base、trinucleotide_context。MutationMatrix 存储 96 通道突变计数矩阵（96 行 × n 样本列），含 channel_labels（如 "A[C>A]A"）和 per-sample total_mutations。MutSigProfile 表示突变签名（id + 96 概率分布）。trinucleotide_contexts() 生成 16 种三核苷酸上下文（5'碱基 × 3'碱基）。get_channel_labels() 生成标准顺序的 96 通道标签。get_mutation_type(ref, alt) 通过嘌呤到嘧啶的反向互补转换（complement_base）确保统一使用嘧啶参考基准。normalize_context(trinucleotide) 将中间碱基为嘌呤的三核苷酸上下文反向互补为嘧啶基准。get_channel_index(ref, alt, trinucleotide) 计算 0-95 通道索引（type_idx * 16 + ctx_idx）。build_mutation_matrix(mutations, sample_ids) 主函数：构建样本索引映射、初始化 96×n 计数矩阵、按通道累加、计算 per-sample 总数。MutationMatrix::to_relative_frequencies() 将计数转为相对频率（每样本和为 1）。cosine_similarity(v1, v2) 计算两向量余弦相似度（dot / (||v1|| × ||v2||)）。MutationMatrix::sample_cosine_similarity(i, j) 比较两样本突变谱相似度。known_cosmic_signatures() 返回简化版 SBS1（CpG 位点 C>T 富集）和 SBS5（平坦分布）签名。calculate_exposure(sample_profile, signature) 通过投影计算签名暴露量。fit_signatures(sample_profile, signatures) 通过余弦相似度归一化拟合多签名贡献比例。辅助函数 all_mutation_types() 返回 6 种类型标准顺序，MutationMatrix::spectrum_summary() 生成 6 类型计数摘要，mutational_patterns_sample_data() 提供 9 条示例突变（5 条 S1 + 4 条 S2）。适用于癌症基因组体细胞突变模式分析、突变签名解析、暴露量估计。

### 170. GAGE 基因集富集分析 (Bioconductor gage)

实现基于 fold change 的基因集富集分析，参考 Bioconductor gage 包。核心结构 GageGeneSet（id、name、gene_ids）定义基因集。GageResult 存储每个基因集的检验结果：gene_set_id、n_genes、mean_fc、stat（t/z 统计量）、p_val、p_adj、is_up / is_down、significant。GageResults 汇总所有结果及全局参数（n_gene_sets、n_significant、fdr_threshold、paired）。run_gage(gene_fcs, gene_sets, fdr_threshold?, paired?) 主函数：1) 对每个基因集收集成员基因的 fold change；2) 计算 mean_fc（gage_mean）；3) 计算标准误 se = sqrt(variance / n)，配对检验时方差减半（gage_se）；4) 计算统计量 stat = mean_fc / se；5) 通过正态近似（erfc_approx，使用 Abramowitz & Stegun 7.1.26 公式）计算双侧 p 值（gage_pvalue_from_z）；6) 使用 BH-FDR 校正（gage_bh_adjust：先按 p 值排序，再从大到小回填调整值，确保单调性）。辅助函数 GageResults::get_significant() 返回显著基因集，get_upregulated() / get_downregulated() 分别返回上调/下调基因集，get_top_sets(n) 按调整 p 值排序返回前 N 个，summary() 生成分析摘要。gage_sample_gene_sets() 提供 3 个示例基因集（Pathway A/B/C），gage_sample_fold_changes() 提供 7 个基因的 fold change（含正负值）。适用于 RNA-seq / 微阵列数据的通路富集分析，支持配对/非配对检验。

### 171. SPIA 信号通路影响分析 (Bioconductor SPIA)

实现信号通路影响分析，参考 Bioconductor SPIA 包（Tarca AL et al. 2009）。核心结构 PathwayNode 表示通路节点：gene_id、gene_name、downstream（下游靶基因列表）、upstream（上游调控基因列表）、activation（每条下游互作的激活符号 +1/-1）。SignalingPathway（id、name、nodes）表示完整通路。SpiaResult 存储每通路分析结果：n_genes、n_de、p_or（过表达 p 值）、pert_factor（扰动累积因子）、p_pert（扰动 p 值）、p_combined（合并 p 值）、p_adj、significant、activation_status（+1 激活 / -1 抑制 / 0 中性）。SpiaResults 汇总所有通路结果。PathwayNode::add_downstream(target, is_activation) 添加下游互作并记录激活/抑制符号。run_spia(pathways, de_genes, total_genes, fdr_threshold?) 主分析流程：1) 对每通路统计 DE 基因数；2) 超几何检验 p_or（hypergeometric_pvalue，使用卡方近似 N(0,1)）；3) 计算扰动累积因子 pert_factor（compute_perturbation：每个基因的扰动 = 自身 logFC + Σ 上游 DE 基因的 logFC × 激活符号，最后除以 DE 基因数归一化）；4) 扰动 p 值 p_pert（perturbation_pvalue：使用正态近似 N(0, n_de) 与 erfc 函数）；5) Fisher 合并 p_combined = -2 × (ln p_or + ln p_pert) → 卡方分布 p 值（fisher_combine）；6) BH-FDR 校正（spia_bh_adjust）；7) 根据 pert_factor 符号判定激活/抑制状态。辅助函数 SpiaResults::get_significant() / get_activated() / get_inhibited() / get_top_pathways(n) / summary() 提供结果筛选和摘要。spia_sample_pathways() 提供 2 个示例通路：MAPK Cascade（5 基因激活级联：EGFR→RAS→RAF→MEK→ERK）和 Apoptosis（3 基因含抑制：BAX→CASP3, BCL2 -| CASP3）。spia_sample_de_genes() 提供 4 个 DE 基因（Gene1, Gene2, Gene5, Gene6）。适用于差异表达数据的通路影响分析，结合过表达和扰动信号判定通路激活/抑制状态。

### 172. Bio.File 智能文件处理 (Biopython Bio.File)

实现智能文件处理模块，提供透明压缩格式支持。核心结构 SmartFile 表示文件句柄：path（文件路径）、format（压缩格式：Plain/Gzip/Bzip2）、mode（读写模式）、is_open（文件打开状态）、lines（缓冲行数组，用于读取）、pos（当前读取位置）、written_content（写入内容缓冲）。CompressionFormat 枚举支持三种格式：Plain（无压缩）、Gzip（gzip压缩，.gz扩展名）、Bzip2（bzip2压缩，.bz2扩展名）。自动检测逻辑：SmartFile::new(path, mode?) 根据文件扩展名自动识别压缩格式（.gz → Gzip、.bz2 → Bzip2、其他 → Plain）。支持文件打开/关闭操作（open/close）、逐行读取（read_line）、全量读取（read_all）、字符串写入（write）、追加写入（append）。提供工具函数：detect_compression(path) 检测压缩格式、has_compression_extension(path) 判断是否有压缩扩展名、strip_compression_extension(path) 去除压缩扩展名。适用于生物信息学中常见的压缩序列文件（如 .fasta.gz、.fastq.bz2）的透明读写。

### 173. Bio.SeqUtils.MolWt 分子量计算 (Biopython Bio.SeqUtils.MolWt)

实现分子量计算模块，支持 DNA、RNA 和蛋白质序列的理化性质分析。核心结构 AtomicWeights 存储元素原子量（H/C/N/O/P/S）。支持的计算功能包括：（1）mol_weight_dna(sequence) 计算 DNA 分子量，公式为核苷酸权重之和减去 (n-1)×水分子量（18.01524）；（2）mol_weight_rna(sequence) 计算 RNA 分子量，尿嘧啶替换胸腺嘧啶；（3）mol_weight_protein(sequence) 计算蛋白质分子量，基于氨基酸单同位素质量；（4）extinction_coefficient(sequence) 计算 280nm 消光系数（基于 Trp、Tyr、Cys 残基）；（5）absorbance_280(sequence) 计算 280nm 吸光度；（6）mol_wt_isoelectric_point(sequence) 计算等电点 pI（基于 Lehninger pK 值）。提供示例序列：mol_wt_sample_dna()、mol_wt_sample_protein()。支持序列摘要生成（mol_weight_summary）。适用于蛋白质组学和分子生物学实验中的样品定量与质控。

### 174. Bio.Align.Reduced 简化氨基酸字母表 (Biopython Bio.Align.Reduced)

实现简化氨基酸字母表模块，用于蛋白质序列的简化表示和比较。核心结构 ReducedAlphabet 包含：name（名称）、n_groups（分组数）、mapping（氨基酸到分组字母的映射）、groups（分组定义）。支持四种常用简化字母表：（1）RAD（Reduced Alphabet Database，6组）：(A,G)、(C)、(D,E,N,Q)、(I,L,M,V)、(F,Y,W)、(H,K,R,S,T)，来自 Wang & Wang (1999)；（2）Dayhoff（6组）：(A,G,P,S,T)、(C)、(D,E,N,Q)、(I,L,M,V)、(F,W,Y)、(H,K,R)；（3）CHARM（4组）：(A,C,F,I,L,M,V)、(G,S,T,P)、(D,E,N,Q)、(H,K,R,W,Y)，来自 Li et al. (2003)；（4）SDM12（12组）：基于结构域记忆性的12组简化。支持序列简化（reduce_sequence）、单氨基酸简化（reduce_aa）、简化序列比较（reduced_identity）和字母表摘要（summary）。适用于大规模蛋白质序列比对的预筛选、远缘同源性检测和序列聚类分析。

### 175. GENIE3 基因调控网络推断 (Bioconductor GENIE3)

实现 GENIE3（GEne Network Inference with Ensemble of trees）算法，参考 Huynh-Thu VA et al. (2010) PLoS ONE。核心思想：对每个目标基因，使用其他所有基因作为预测变量构建回归树，通过特征重要性（方差缩减）推断调控关系。核心结构 Genie3TreeNode 表示回归树节点：feature（分裂特征索引，-1 为叶节点）、threshold（分裂阈值）、left_idx/right_idx（子节点索引）、value（叶节点预测值）、is_leaf、importance_gain（该分裂的方差缩减量）。RegressionTree 封装节点数组，支持 predict(sample) 预测和 feature_importance(n_features) 计算特征重要性（累加内部节点的 importance_gain）。genie3_build_tree(features, target, feature_indices, max_depth?, min_samples_split?) 构建单棵回归树：递归选择最佳分裂点（遍历所有特征和阈值，最大化 SSE 减少），支持最大深度和最小样本数约束。genie3_run(expression, gene_names, max_depth?, min_samples_split?, symmetrize?) 主推理函数：1) 对每个目标基因 j 构建回归树（排除 j 自身作为预测变量）；2) 累加各树的特征重要性到权重矩阵 matrix[i][j]；3) 列归一化（每列和为1）；4) 可选对称化（matrix[i][j] 和 matrix[j][i] 取平均）；5) 收集非零边并按权重降序排列。RegulatoryEdge 结构（regulator、target、weight）表示调控边。辅助函数：genie3_top_edges(result, k) 获取前 k 条边、genie3_regulators_of(result, gene) 查询靶基因的调控子、genie3_targets_of(result, gene) 查询调控子的靶基因、genie3_sample_data() 提供示例数据（5 基因、30 样本，G1 调控 G2/G3，G2 调控 G4，G3 调控 G5）。适用于基因调控网络推断、转录因子靶基因预测和共表达网络构建。

### 175. decoupleR 功能活性推断 (Bioconductor decoupleR)

实现 decoupleR 功能活性推断框架，参考 Badia-i-Mompel P et al. (2022) Molecular Systems Biology。核心思想：基于先验知识网络（PKN）从组学数据推断调控子（如转录因子）的活性。PKNEdge 结构（source、target、weight）表示调控边（weight>0 激活、weight<0 抑制）。PriorKnowledgeNetwork（edges、regulators、targets）封装 PKN。pkn_from_edges(edges) 从边列表构建 PKN 并提取唯一调控子和靶基因。支持五种推理方法：（1）WSum（加权求和）：score(r,s) = Σ_t W[r][t] × expr[s][t]；（2）WMean（加权平均）：WSum 除以绝对权重和；（3）Norm（归一化均值）：WMean 的标准化版本；（4）ULM（单变量线性模型）：计算每调控子的加权求和后标准化为 z-score；（5）MLM（多变量线性模型）：对每个调控子，将其靶基因表达对其权重向量做回归，斜率即活性。DecoupleRMethod 枚举封装五种方法，通过 decoupler_run(expression, sample_names, gene_names, pkn, method?) 调用。DecoupleRResult（regulators、samples、matrix、scores、method）存储结果矩阵和逐样本评分。ActivityScore（regulator、sample、score、p_value）表示单个评分记录。辅助函数：decoupler_top_regulators(result, sample, k) 获取某样本前 k 个调控子、decoupler_filter_scores(result, threshold) 按阈值过滤评分、decoupler_sample_data() 提供示例数据（3 TF、11 基因、5 样本、7 条 PKN 边）。适用于转录因子活性推断、信号通路活性分析和单细胞调控网络分析。

### 177. BayesSpace 空间转录组聚类 (Bioconductor BayesSpace)

实现 BayesSpace 空间转录组聚类算法，参考 Zhao E et al. (2021) Nature Biotechnology。核心思想：使用 t 分布混合模型结合马尔可夫随机场（MRF）先验对空间转录组 spot 进行聚类，鼓励相邻 spot 共享聚类标签。SpotCoord 结构（spot_id、row、col、x、y）表示 spot 坐标。bayes_space_hex_neighbors(spots) 构建六边形网格邻居列表（6-连通：6 个候选邻居），bayes_space_square_neighbors(spots) 构建方形网格邻居列表（4-连通）。bayes_space_run(expression, spots, q, neighbors, max_iters?, gamma?, df?, seed?, tol?) 主聚类函数使用 EM 算法：1) k-means++ 初始化聚类中心（第一个中心随机选取，后续选择距已有中心最远的 spot）；2) E-step：计算每个 spot 属于每个簇的后验责任，综合 t 分布似然（bayes_space_log_t_density：多变量 t 分布对数密度）和空间先验（Potts 模型：gamma × Σ 邻居责任）；3) 通过 log-sum-exp 归一化责任；4) M-step：更新聚类中心（加权均值）、尺度（加权标准差）和混合比例；5) 计算对数似然判断收敛。BayesSpaceResult（spot_ids、clusters、q、cluster_centers、responsibilities、n_iterations、log_likelihood、converged）存储完整结果。辅助函数：bayes_space_get_cluster(result, spot_id) 按 ID 查询簇、bayes_space_spots_in_cluster(result, k) 获取簇内 spot、bayes_space_cluster_counts(result) 统计各簇大小、bayes_space_render_clusters(result, spots) 生成文本聚类图、bayes_space_sample_data() 提供 4×4 网格示例数据（左上低表达、右下高表达）。适用于空间转录组数据分析、组织区域识别和空间异质性研究。


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
| 总测试数 | 4916 |
| 通过数 | 4916 |
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
| Motifs Advanced | `motifs_advanced_test.mbt` | 45 |
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
| Phylo NeXML | `phylo_nexml_test.mbt` | 32 |
| TxDb | `txdb_test.mbt` | 14 |
| ProtParam | `protparam_test.mbt` | 15 |
| rtracklayer | `rtracklayer_test.mbt` | 19 |
| CodonUsage | `codon_usage_test.mbt` | 8 |
| GenomicAlignments | `genomic_alignments_test.mbt` | 7 |
| VariantAnnotation | `variant_annotation_test.mbt` | 6 |
| Affy | `affy_test.mbt` | 8 |
| SVDSuperimposer | `svd_superimposer_test.mbt` | 9 |
| PDB Dice | `pdb_dice_test.mbt` | 42 |
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
| BiocNeighbors | `bioc_neighbors_test.mbt` | 61 |
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
| Rtsne | `rtsne_test.mbt` | 10 |
| uwot | `uwot_test.mbt` | 9 |
| microbiome | `microbiome_test.mbt` | 33 |
| tradeSeq | `tradeseq_test.mbt` | 12 |
| QCP叠加 | `qcp_superimposer_test.mbt` | 7 |
| 残基深度 | `residue_depth_test.mbt` | 10 |
| 结构比对 | `structure_alignment_test.mbt` | 8 |
| PDB向量 | `pdb_vectors_test.mbt` | 55 |
| PhyloXML调试 | `phylo_xml_debug_test.mbt` | 8 |
| Seurat | `seurat_test.mbt` | 12 |
| 统计分析 | `statistics_test.mbt` | 16 |
| topGO | `topgo_test.mbt` | 6 |
| 字母表 | `alphabet_test.mbt` | 8 |
| 频率分析 | `freq_analysis_test.mbt` | 13 |
| nucleR | `nucle_r_test.mbt` | 13 |
| 比对分析 | `align_analysis_test.mbt` | 15 |
| bamsignals | `bamsignals_test.mbt` | 10 |
| DSS | `dss_test.mbt` | 12 |
| 表型分析 | `phenotype_test.mbt` | 12 |
| BLAST应用 | `blast_applications_test.mbt` | 19 |
| PSEA | `psea_test.mbt` | 14 |
| SFF_IO | `sff_io_test.mbt` | 16 |
| csaw | `csaw_test.mbt` | 9 |
| slingshot | `slingshot_test.mbt` | 10 |
| SCnorm | `scnorm_test.mbt` | 8 |
| EDASeq | `edaseq_test.mbt` | 10 |
| SearchIO新 | `searchio_new_test.mbt` | 30 |
| SeqFeature高级 | `seqfeature_advanced_test.mbt` | 63 |
| 空间转录组 | `spatial_experiment_test.mbt` | 4 |
| S4Vectors | `s4vectors_test.mbt` | 7 |
| rhdf5 | `rhdf5_test.mbt` | 12 |
| PhyloXML | `phylo_xml_test.mbt` | 39 |
| phyloseq | `phyloseq_test.mbt` | 5 |
| 多肽分析 | `polypeptide_test.mbt` | 8 |
| PROGENy | `progeny_test.mbt` | 13 |
| Prosite | `prosite_test.mbt` | 12 |
| SVA | `sva_test.mbt` | 6 |
| 变异过滤 | `variant_filtering_test.mbt` | 6 |
| 变异分析 | `variation_test.mbt` | 6 |
| universalmotif | `universalmotif_test.mbt` | 3 |
| SC3 | `sc3_test.mbt` | 7 |
| clusterProfiler | `cluster_profiler_test.mbt` | 9 |
| 渐进式比对 | `align_cluster_test.mbt` | 17 |
| RNA结构 | `rna_structure_test.mbt` | 19 |
| mixOmics | `mix_omics_test.mbt` | 42 |
| MAalign | `ma_align_test.mbt` | 17 |
| 内部坐标 | `internal_coords_test.mbt` | 10 |
| 遗传算法 | `ga_test.mbt` | 10 |
| 染色体可视化 | `chromosome_visualization_test.mbt` | 10 |
| impute | `impute_test.mbt` | 13 |
| vsn | `vsn_test.mbt` | 9 |
| GSEABase | `gsea_base_test.mbt` | 20 |
| PCAtools | `pcatools_test.mbt` | 11 |
| Bio.Data | `data_test.mbt` | 15 |
| Seq.Approximate | `seq_approx_test.mbt` | 21 |
| Pairwise2 | `pairwise2_test.mbt` | 16 |
| Compound | `compound_test.mbt` | 22 |
| EnhancedVolcano | `enhanced_volcano_test.mbt` | 15 |
| ReportingTools | `reporting_tools_test.mbt` | 12 |
| karyoploteR | `karyoploter_test.mbt` | 14 |
| SystemPipeR | `system_piper_test.mbt` | 16 |
| muscat | `muscat_test.mbt` | 17 |
| infercnv | `infercnv_test.mbt` | 14 |
| SCENIC | `scenic_test.mbt` | 22 |
| MSstats | `msstats_test.mbt` | 15 |
| NOISeq | `noiseq_test.mbt` | 17 |
| Gviz | `gviz_test.mbt` | 16 |
| HTSFilter | `htsfilter_test.mbt` | 15 |
| baySeq | `bayseq_test.mbt` | 18 |
| CellChat | `cellchat_test.mbt` | 15 |
| MutationalPatterns | `mutational_patterns_test.mbt` | 19 |
| GAGE | `gage_test.mbt` | 11 |
| SPIA | `spia_test.mbt` | 22 |
| Bio.File | `file_test.mbt` | 14 |
| Bio.SeqUtils.MolWt | `mol_wt_test.mbt` | 24 |
| Bio.Align.Reduced | `reduced_test.mbt` | 20 |
| Hmisc | `hmisc_test.mbt` | 68 |
| rstatix | `rstatix_test.mbt` | 71 |
| VennDiagram | `venn_diagram_test.mbt` | 68 |
| GENIE3 | `genie3_test.mbt` | 18 |
| decoupleR | `decoupler_test.mbt` | 17 |
| BayesSpace | `bayes_space_test.mbt` | 17 |

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

项目提供 188 个示例程序，展示各模块的典型用法：

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
| uniprot_io_demo | UniProt XML格式解析（蛋白质条目解析、功能注释提取、序列转换） | `moon run examples/uniprot_io_demo/main.mbt` |
| chem_utils_demo | 化学计算工具（键长、键角、二面角、分子式量、氢键长度） | `moon run examples/chem_utils_demo/main.mbt` |
| jaspar_demo | JASPAR PFM格式解析（模体矩阵解析、共有序列、PWM转换、序列扫描） | `moon run examples/jaspar_demo/main.mbt` |
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
| internal_coords_demo | 蛋白质内部坐标（二面角计算、扩展链构建、Ramachandran区域、旋转异构体库） | `moon run examples/internal_coords_demo/main.mbt` |
| ga_demo | 遗传算法序列优化（个体/种群创建、锦标赛/轮盘赌选择、单点/两点交叉、变异、进化收敛） | `moon run examples/ga_demo/main.mbt` |
| chromosome_visualization_demo | 染色体可视化（特征/区域/带创建、SVG渲染、G带颜色映射、线性/圆形染色体） | `moon run examples/chromosome_visualization_demo/main.mbt` |
| mutational_patterns_demo | MutationalPatterns 体细胞突变谱分析（96通道矩阵、突变类型识别、上下文归一化、签名拟合、余弦相似度） | `moon run examples/mutational_patterns_demo/main.mbt` |
| gage_demo | GAGE 基因集富集分析（fold change、t检验、BH-FDR校正、配对检验、上调/下调分析） | `moon run examples/gage_demo/main.mbt` |
| spia_demo | SPIA 信号通路影响分析（通路图、扰动累积、超几何检验、Fisher合并p值、激活/抑制判定） | `moon run examples/spia_demo/main.mbt` |
| file_demo | Bio.File 智能文件处理（压缩检测、透明读写、文件操作接口） | `moon run examples/file_demo/main.mbt` |
| mol_wt_demo | Bio.SeqUtils.MolWt 分子量计算（DNA/RNA/蛋白质分子量、消光系数、等电点） | `moon run examples/mol_wt_demo/main.mbt` |
| reduced_demo | Bio.Align.Reduced 简化氨基酸字母表（RAD/Dayhoff/CHARM字母表、序列比较） | `moon run examples/reduced_demo/main.mbt` |

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
- ✅ 实现 microbiome 微生物组分析（Alpha多样性、Beta多样性、PCoA排序、差异丰度分析）
- ✅ 实现 Rtsne t-SNE降维算法（距离矩阵、条件概率、梯度下降、动量优化）
- ✅ 实现 uwot UMAP降维算法（k近邻、模糊单纯集、SGD优化、负采样）
- [ ] 添加 SIMD 加速支持

## 许可证

Apache-2.0 License
