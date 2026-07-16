# BioSeqs - MoonBit 生物信息学库

> https://github.com/paipai-Studio/BioSeqs
>
> https://gitlink.org.cn/IvanAXu/BioSeqs
> 
> https://mooncakes.io/docs/IvanAXu/BioSeqs
>

## 项目概述

BioSeqs 是一个基于 **MoonBit** 语言开发的生物信息学工具库，旨在复刻主流生物信息学库的核心功能，并实现高效的序列组装算法。

### 核心功能覆盖

| 功能类别 | 对应工具/库 | 功能范围 | 状态 |
|-----------|-------------|----------|------|
| **序列处理** | Biopython `Bio.Seq` | 序列对象、互补、转录、翻译、序列特征 | ✅ 已实现 |
| **序列 I/O** | Biopython `Bio.SeqIO` | FASTA/FASTQ/GenBank 解析与写入 | ✅ 已实现 |
| **序列比对** | Biopython / scikit-bio | Needleman-Wunsch、Smith-Waterman、多序列比对、替换矩阵(BLOSUM/PAM) | ✅ 已实现 |
| **BLAST解析** | Biopython `Bio.Blast` | BLAST结果解析、tabular/xml格式、HSP过滤、最佳匹配 | ✅ 已实现 |
| **SearchIO** | Biopython `Bio.SearchIO` | 统一搜索结果模型、HMMER3解析、BLAT PSL解析、BLAST转换 | ✅ 已实现 |
| **系统发育树** | Biopython `Bio.Phylo` | 树结构、Newick 解析、距离计算、可视化 | ✅ 已实现 |
| **PDB 结构** | Biopython `Bio.PDB` | 原子/残基/链解析、结构操作 | ✅ 已实现 |
| **SAM/BAM/VCF** | pysam | 比对文件、变异检测、基因型查询 | ✅ 已实现 |
| **FASTA 索引** | pyfaidx | 快速随机访问、.fai 索引 | ✅ 已实现 |
| **机器学习特征** | scikit-learn | k-mer 频率、氨基酸组成、理化性质 | ✅ 已实现 |
| **Biostrings** | Bioconductor Biostrings | IUPAC 支持、RSCU、复杂度、Tm 计算 | ✅ 已实现 |
| **GenomicRanges** | Bioconductor GenomicRanges | GRanges、区间操作、集合运算 | ✅ 已实现 |
| **DESeq2** | Bioconductor DESeq2 | 差异表达分析、负二项分布模型 | ✅ 已实现 |
| **dplyr** | R dplyr | DataFrame 数据操作 | ✅ 已实现 |

### 序列组装算法

| 算法 | 对应工具 | 功能说明 | 状态 |
|------|----------|----------|------|
| **De Bruijn Graph** | SPAdes / Velvet | k-mer 节点、欧拉路径、序列组装 | ✅ 已实现 |
| **Suffix Array & Suffix Tree** | libdivsufsort | 前缀倍增算法、LCP 数组、模式匹配、最长重复子串 | ✅ 已实现 |
| **Overlap-Layout-Consensus** | Celera Assembler / ARACHNE | 重叠检测、哈密顿路径、一致性序列生成 | ✅ 已实现 |
| **BWT + FM-index** | Bowtie2 / BWA | 后缀数组、BWT 变换、回溯搜索、精确模式匹配 | ✅ 已实现 |

### 哈希与数据结构

| 数据结构 | 对应工具 | 功能说明 | 状态 |
|----------|----------|----------|------|
| **Bloom Filter** | Jellyfish / khmer | k-mer 计数、成员查询、误判率估算 | ✅ 已实现 |

### 序列模体与群体遗传学

| 功能类别 | 对应工具 | 功能说明 | 状态 |
|----------|----------|----------|------|
| **序列模体** | Biopython `Bio.motifs` | PWM、MEME格式、模体搜索、信息含量 | ✅ 已实现 |
| **群体遗传学** | Biopython `Bio.PopGen` | 等位基因频率、FST、哈迪-温伯格检验 | ✅ 已实现 |
| **edgeR** | Bioconductor edgeR | 差异表达分析、DGEList、精确检验、GLM拟合 | ✅ 已实现 |
| **limma** | Bioconductor limma | 差异表达分析、线性模型拟合、经验贝叶斯、voom变换 | ✅ 已实现 |
| **SummarizedExperiment** | Bioconductor SummarizedExperiment | 多维基因组数据容器、Assays、行/列操作 | ✅ 已实现 |
| **IRanges** | Bioconductor IRanges | 整数区间操作、集合运算、重叠检测 | ✅ 已实现 |
| **TxDb** | Bioconductor GenomicFeatures | 转录本数据库、GTF解析、基因/转录本/外显子/CDS提取、UTR/内含子计算、启动子提取 | ✅ 已实现 |
| **限制性内切酶分析** | Biopython `Bio.Restriction` | 内切酶数据库、酶切位点查找、序列酶切、片段分析 | ✅ 已实现 |
| **序列聚类分析** | Biopython `Bio.Cluster` | 距离矩阵、层次聚类、Newick输出、轮廓系数 | ✅ 已实现 |
| **比对格式解析** | Biopython `Bio.AlignIO` | ClustalW、FASTA、Stockholm格式解析与写入 | ✅ 已实现 |
| **进化树格式解析** | Biopython `Bio.TreeIO` | Newick、NHX格式解析、树操作、修剪 | ✅ 已实现 |
| **蛋白质参数分析** | Biopython `Bio.SeqUtils.ProtParam` | 分子量、不稳定指数、GRAVY评分、等电点、信号肽预测 | ✅ 已实现 |
| **基因组轨道格式** | Bioconductor rtracklayer | BED/WIG/BEDGraph/GFF解析与写入、GRanges转换 | ✅ 已实现 |
| **密码子使用分析** | Biopython `Bio.codonalign` / `Bio.CodonUsage` | CAI、ENC、RSCU、GC3、CBI、Fop、最优密码子检测 | ✅ 已实现 |
| **GenomicAlignments** | Bioconductor GenomicAlignments | GAlignments对象、coverage计算、summarizeOverlaps、pileup操作 | ✅ 已实现 |
| **VariantAnnotation** | Bioconductor VariantAnnotation | 变异类型检测、变异定位、编码效应预测、变异汇总 | ✅ 已实现 |
| **Affy** | Biopython `Bio.Affy` | Affymetrix芯片数据分析、RMA标准化、背景校正、分位数归一化 | ✅ 已实现 |
| **SVDSuperimposer** | Biopython `Bio.PDB.SVDSuperimposer` | SVD蛋白质结构叠合、旋转矩阵、平移向量、RMSD计算 | ✅ 已实现 |

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
│   ├── seq_utils.mbt           # 序列工具函数 (分子量、GC含量、Tm值、氨基酸转换)
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
│   ├── deseq2.mbt              # DESeq2 差异表达分析 (负二项分布模型、Wald检验、Benjamini-Hochberg校正)
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
│   ├── motifs.mbt              # 序列模体识别 (PWM、MEME格式、模体搜索)
│   ├── popgen.mbt              # 群体遗传学 (等位基因频率、FST、哈迪-温伯格检验)
│   ├── restriction.mbt         # 限制性内切酶分析 (酶切位点查找、片段分析)
│   ├── protparam.mbt           # ProtParam 蛋白质参数分析 (不稳定指数、等电点、信号肽预测、二级结构倾向)
│   ├── affy.mbt                # Affy Affymetrix芯片数据分析 (RMA标准化、背景校正、分位数归一化)
│   ├── feature_extraction.mbt  # 机器学习特征提取
│   ├── faidx.mbt               # FASTA 快速索引访问 (pyfaidx)
│   ├── go_enrichment.mbt       # GOEnrichment GO功能富集分析 (超几何检验、BH校正、富集结果过滤)
│   ├── single_cell.mbt         # SingleCell 单细胞数据分析 (QC指标、Log标准化、PCA降维、高变异基因)
│   └── utils.mbt               # 通用工具函数
├── examples/                   # 示例程序
│   ├── affy_demo/              # Affy Affymetrix芯片数据分析示例 (RMA标准化、背景校正、分位数归一化)
│   ├── align_io_demo/          # AlignIO 比对格式解析示例 (ClustalW、FASTA、Stockholm)
│   ├── alignment_demo/         # 序列比对示例
│   ├── bam_demo/               # BAM/BGZF 解析示例
│   ├── basic_seq/              # 基础序列操作示例
│   ├── biostrings_demo/        # Biostrings 序列分析示例
│   ├── blast_demo/             # BLAST结果解析示例 (tabular/xml格式、HSP过滤、E-value/identity过滤、最佳匹配)
│   ├── bloom_filter_demo/      # Bloom Filter & k-mer 计数示例
│   ├── bwt_fm_demo/            # BWT + FM-index 示例
│   ├── cluster_demo/           # 序列聚类分析示例 (层次聚类、距离矩阵、轮廓系数)
│   ├── codon_usage_demo/       # CodonUsage 密码子使用分析示例 (CAI、ENC、RSCU、GC3、CBI、Fop、最优密码子检测)
│   ├── de_bruijn_demo/         # De Bruijn Graph 序列组装示例
│   ├── deseq2_demo/            # DESeq2 差异表达分析示例
│   ├── dplyr_demo/             # dplyr 数据操作示例
│   ├── edger_demo/             # edgeR 差异表达分析示例
│   ├── faidx_demo/             # FASTA 索引示例
│   ├── genomic_alignments_demo/ # GenomicAlignments 基因组比对分析示例 (GAlignments、coverage、summarizeOverlaps、pileup)
│   ├── genomic_ranges_demo/    # GenomicRanges 基因组区间操作示例
│   ├── go_enrichment_demo/     # GOEnrichment GO功能富集分析示例 (超几何检验、BH校正、富集结果过滤)
│   ├── hmm_demo/               # Hidden Markov Model 基因预测示例
│   ├── iranges_demo/           # IRanges 区间操作示例
│   ├── kmeans_demo/            # K-means Clustering 聚类分析示例
│   ├── limma_demo/             # limma 差异表达分析示例 (线性模型、经验贝叶斯、voom变换)
│   ├── ml_features/            # 机器学习特征提取示例
│   ├── motifs_demo/            # 序列模体识别示例
│   ├── needleman_wunsch_demo/  # Needleman-Wunsch 全局序列比对示例
│   ├── olc_demo/               # Overlap-Layout-Consensus 序列组装示例
│   ├── pdb_demo/               # PDB 结构解析示例
│   ├── phylo_demo/             # 系统发育树示例
│   ├── popgen_demo/            # 群体遗传学分析示例
│   ├── protparam_demo/         # ProtParam 蛋白质参数分析示例 (分子量、不稳定指数、等电点、信号肽预测、二级结构倾向)
│   ├── restriction_demo/       # 限制性内切酶分析示例 (酶切位点查找、片段分析)
│   ├── rtracklayer_demo/       # rtracklayer 基因组轨道格式示例 (BED/WIG/BEDGraph/GFF解析与写入)
│   ├── sam_vcf_demo/           # SAM/VCF 解析示例
│   ├── search_io_demo/         # SearchIO 统一搜索结果示例 (HMMER3解析、BLAT PSL解析、BLAST转换)
│   ├── seqio_demo/             # 序列 I/O 示例
│   ├── seq_utils_demo/         # 序列工具函数示例
│   ├── single_cell_demo/       # SingleCell 单细胞数据分析示例 (QC指标、Log标准化、PCA降维、高变异基因)
│   ├── smith_waterman_demo/    # Smith-Waterman 局部序列比对示例
│   ├── subsmat_demo/           # 替换矩阵示例 (BLOSUM62/45、PAM250/30矩阵查询、蛋白质比对打分)
│   ├── suffix_array_tree_demo/ # Suffix Array & Suffix Tree 示例
│   ├── summarized_experiment_demo/ # SummarizedExperiment 数据容器示例
│   ├── svd_superimposer_demo/  # SVDSuperimposer SVD蛋白质结构叠合示例 (旋转矩阵、平移向量、RMSD计算)
│   ├── txdb_demo/              # TxDb 转录本数据库示例 (GTF解析、基因/转录本/外显子/CDS/UTR/内含子提取)
│   └── variant_annotation_demo/ # VariantAnnotation 变异注释示例 (变异类型检测、定位、编码效应预测、变异汇总)
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
│   │   └── vcf_test.mbt
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
moon test --package IvanAXu/BioSeqs/test/moonbit        # ✅ 707 个测试全部通过
```

### 模块对照表

| MoonBit 文件 | 对应 Python 库 | 核心功能 |
|-------------|---------------|----------|
| `seq.mbt` | BioPython `Bio.Seq` | 序列对象、互补、转录、翻译 |
| `seq_record.mbt` | BioPython `Bio.SeqRecord` | 带注释的序列记录 |
| `seqfeature.mbt` | BioPython `Bio.SeqFeature` | 序列特征与位置 |
| `seqio.mbt` | BioPython `Bio.SeqIO` | 统一序列文件 I/O |
| `fasta_io.mbt` | BioPython `Bio.SeqIO.FastaIO` | FASTA 解析 |
| `fastq_io.mbt` | BioPython `Bio.SeqIO.QualityIO` | FASTQ 解析 |
| `genbank_io.mbt` | BioPython `Bio.SeqIO.GenBankIO` | GenBank 解析 |
| `align.mbt` | BioPython `Bio.Align` | 多序列比对对象 |
| `alignio.mbt` | BioPython `Bio.AlignIO` | 比对文件 I/O |
| `clustal_io.mbt` | BioPython `Bio.AlignIO.ClustalIO` | Clustal 格式 |
| `phylip_io.mbt` | BioPython `Bio.AlignIO.PhylipIO` | PHYLIP 格式 |
| `alignment.mbt` | scikit-bio `skbio.alignment` | DNA/RNA/Protein 类型、Needleman-Wunsch/Smith-Waterman |
| `phylo.mbt` | BioPython `Bio.Phylo` | 系统发育树 (Clade/Tree) |
| `pdb.mbt` | BioPython `Bio.PDB` | PDB 数据类型 |
| `pdb_io.mbt` | BioPython `Bio.PDB.PDBIO` | PDB 文件 I/O |
| `sequtils.mbt` | BioPython `Bio.SeqUtils` | CRC32、GC 含量、Tm 计算、蛋白质分析 |
| `complement.mbt` | BioPython `Bio.Data.IUPACData` | 互补碱基表 |
| `codon_table.mbt` | BioPython `Bio.Data.CodonTable` | 密码子翻译表 |
| `sam.mbt` | pysam | SAM 文件解析 |
| `bam.mbt` | pysam | BAM 文件解析 |
| `bgzf.mbt` | pysam | BGZF 解压缩 |
| `biostrings.mbt` | Bioconductor Biostrings | IUPAC 支持、k-mer 频率、RSCU、复杂度、Tm 计算 |
| `genomic_ranges.mbt` | Bioconductor GenomicRanges | GRanges 数据结构、区间操作、集合操作、重叠检测 |
| `deseq2.mbt` | Bioconductor DESeq2 | 差异表达分析、负二项分布模型、Wald检验、多重检验校正 |
| `dplyr.mbt` | R dplyr | DataFrame 数据操作（filter、select、mutate、arrange、group_by、summarize、join） |
| `smith_waterman.mbt` | BioPython `Bio.Align.PairwiseAligner` | Smith-Waterman 局部序列比对（动态规划、自定义打分、回溯矩阵、得分矩阵导出） |
| `needleman_wunsch.mbt` | BioPython `Bio.Align.PairwiseAligner` | Needleman-Wunsch 全局序列比对（动态规划、自定义打分、回溯矩阵、得分矩阵导出） |
| `bloom_filter.mbt` | Jellyfish / khmer | Bloom Filter 概率成员查询、k-mer 精确计数与近似去重（双哈希、位图、误判率估算） |
| `bwt_fm.mbt` | Bowtie2 / BWA | Burrows-Wheeler Transform + FM-index（后缀数组、BWT正逆变换、回溯搜索、精确模式匹配） |
| `de_bruijn.mbt` | SPAdes / Velvet | De Bruijn Graph 序列组装（k-mer节点、欧拉路径、Hierholzer算法、图简化） |
| `suffix_array_tree.mbt` | libdivsufsort | Suffix Array 和 Suffix Tree（前缀倍增算法、LCP数组、模式匹配、最长重复子串） |
| `olc.mbt` | Celera Assembler / ARACHNE | Overlap-Layout-Consensus 序列组装（重叠检测、哈密顿路径搜索、一致性序列生成） |
| `hmm.mbt` | HMMER / Augustus | Hidden Markov Model（前向算法、后向算法、维特比算法、Baum-Welch训练、基因预测） |
| `kmeans.mbt` | scikit-learn | K-means 聚类（距离计算、K-means++初始化、聚类、轮廓系数评估、最优k值搜索） |
| `vcf.mbt` | pysam | VCF 文件解析 |
| `faidx.mbt` | pyfaidx | FASTA 快速索引访问 |
| `feature_extraction.mbt` | 自定义 | 机器学习特征提取 (k-mer、组成、理化性质) |
| `cram.mbt` | pysam | CRAM 格式解析 (压缩二进制序列比对格式) |
| `motifs.mbt` | Bio.motifs | 序列模体识别 (PWM、MEME格式、模体搜索) |
| `seq_utils.mbt` | Bio.SeqUtils | 序列工具函数 (分子量、GC含量、Tm值、氨基酸转换) |
| `popgen.mbt` | Bio.PopGen | 群体遗传学 (等位基因频率、FST、哈迪-温伯格检验) |
| `edger.mbt` | Bioconductor edgeR | 差异表达分析 (DGEList、精确检验、GLM拟合、P值校正) |
| `limma.mbt` | Bioconductor limma | 差异表达分析 (线性模型拟合、经验贝叶斯、voom变换、对比矩阵、多重检验校正) |
| `summarized_experiment.mbt` | Bioconductor SummarizedExperiment | 多维基因组数据容器 (Assays、行/列操作、合并) |
| `iranges.mbt` | Bioconductor IRanges | 整数区间操作 (shift、resize、reduce、集合运算、重叠检测) |
| `cluster.mbt` | Bio.Cluster | 序列聚类分析 (距离矩阵、层次聚类、轮廓系数) |
| `restriction.mbt` | Bio.Restriction | 限制性内切酶分析 (酶切位点查找、片段分析) |
| `tree_io.mbt` | Bio.TreeIO | 进化树格式解析 (Newick、NHX格式解析、树操作) |
| `align_io.mbt` | Bio.AlignIO | 比对格式解析 (ClustalW、FASTA、Stockholm) |
| `txdb.mbt` | Bioconductor GenomicFeatures | TxDb 转录本数据库 (GTF解析、基因/转录本/外显子/CDS提取、UTR/内含子计算、启动子提取) |
| `subsmat.mbt` | BioPython `Bio.SubsMat` | 替换矩阵 (BLOSUM62/45、PAM250/30、矩阵解析、分数查询、不区分大小写) |
| `blast.mbt` | BioPython `Bio.Blast` | BLAST结果解析 (tabular/xml格式、HSP、Hit、Record、E-value/identity过滤、最佳匹配) |
| `search_io.mbt` | BioPython `Bio.SearchIO` | 统一搜索结果模型 (HSPFragment、HSP、Hit、QueryResult、HMMER3 tabular解析、BLAT PSL解析、BLAST转换、top_hits、count_hsps) |
| `protparam.mbt` | BioPython `Bio.SeqUtils.ProtParam` | ProtParam 蛋白质参数分析 (不稳定指数、等电点、电荷计算、信号肽预测、二级结构倾向) |
| `rtracklayer.mbt` | Bioconductor rtracklayer | rtracklayer 基因组轨道格式 (BED/WIG/BEDGraph/GFF解析与写入、GRanges转换) |
| `codon_usage.mbt` | BioPython `Bio.codonalign` / `Bio.CodonUsage` | CodonUsage 密码子使用分析 (CAI、ENC、RSCU、GC3、CBI、Fop、最优密码子检测、密码子翻译) |
| `genomic_alignments.mbt` | Bioconductor GenomicAlignments | GenomicAlignments 基因组比对分析 (GAlignments对象、coverage计算、summarizeOverlaps、pileup操作、MAPQ/strand过滤、GRanges转换) |
| `variant_annotation.mbt` | Bioconductor VariantAnnotation | VariantAnnotation 变异注释 (变异类型检测、变异定位、编码效应预测、变异汇总、同义/错义/无义变异分类) |
| `affy.mbt` | BioPython `Bio.Affy` | Affy Affymetrix芯片数据分析 (AffyBatch、ProbeSet、RMA标准化、背景校正、分位数归一化、PM-MM差异计算、探针集汇总) |
| `svd_superimposer.mbt` | BioPython `Bio.PDB.SVDSuperimposer` | SVDSuperimposer SVD蛋白质结构叠合 (AtomCoordinate、旋转矩阵、平移向量、RMSD计算、结构叠合、Jacobi特征值分解) |
| `go_enrichment.mbt` | Bioconductor GOstats/clusterProfiler | GO功能富集分析 (GOTerm、超几何检验、Fisher精确检验、Bonferroni/BH校正、富集结果过滤、命名空间筛选) |
| `single_cell.mbt` | Bioconductor SingleCellExperiment | 单细胞数据分析 (SingleCellExperiment、QC指标计算、细胞过滤、Log标准化、PCA降维、高变异基因检测) |
| `utils.mbt` | 自定义 | 通用工具函数 |

## 核心功能实现

### 1. 序列处理 (Bio.Seq)

```moonbit
// 创建序列
let dna = Seq::new("ACGT")
let rna = Seq::new("ACGU")

// 互补操作
dna.complement()          // → Seq("TGCA")
dna.reverse_complement()  // → Seq("ACGT")
rna.complement_rna()      // → Seq("UGCA")

// 转录/反转录
dna.transcribe()          // → Seq("ACGU")
rna.back_transcribe()     // → Seq("ACGT")

// 翻译
dna.translate()           // → Seq("T")
dna.translate(to_stop=true)  // → 翻译到终止密码子
dna.translate(cds=true)      // → 完整 CDS 翻译
```

### 2. 序列 I/O (Bio.SeqIO)

```moonbit
// 解析 FASTA
let records = seqio_parse(fasta_content, "fasta")

// 解析 FASTQ
let records = seqio_parse(fastq_content, "fastq")

// 解析 GenBank
let record = seqio_read(genbank_content, "genbank")

// 写入序列
let text = seqio_write(records, "fasta")
```

### 3. 比对算法 (scikit-bio)

```moonbit
// 创建类型化序列
let dna1 = DNA::new("ACGTACGT")
let dna2 = DNA::new("CGT")

// Needleman-Wunsch 全局比对
let (msa, score, pos) = global_pairwise_align_nucleotide(dna1, dna2)

// Smith-Waterman 局部比对
let (msa, score, pos) = local_pairwise_align_nucleotide(dna1, dna2)

// 蛋白质比对
let prot1 = Protein::new("ACDE")
let prot2 = Protein::new("ACE")
let (msa, score, pos) = global_pairwise_align_protein(prot1, prot2)
```

### 4. SAM 文件解析 (pysam)

```moonbit
// 解析 SAM
let sam = parse_sam(sam_content)

// 访问记录
for record in sam.records {
  record.qname           // 读取名
  record.flag            // 标志位
  record.is_paired()     // 是否配对
  record.is_reverse()    // 是否反向互补
  record.get_cigar()     // CIGAR 数组
}
```

### 5. VCF 文件解析 (pysam)

```moonbit
// 解析 VCF
let vcf = parse_vcf(vcf_content)

// 访问记录
for record in vcf.records {
  record.chrom           // 染色体
  record.pos             // 位置
  record.is_snp()        // 是否 SNP
  record.is_indel()      // 是否插入/缺失
  record.get_info("AC")  // INFO 字段
}
```

### 6. 系统发育树 (Bio.Phylo)

```moonbit
// 解析 Newick
let tree = parse_newick("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")

// 创建树
let clade = Clade::new(name="root", clades=[clade1, clade2])

// 操作
tree.count_terminals()   // → 4
tree.distance("A", "B")  // → 0.3
tree.common_ancestor(["A", "B", "C"])
tree.draw_ascii()        // → ASCII 树图
```

### 7. PDB 结构解析 (Bio.PDB)

```moonbit
// 解析 PDB
let structure = parse_pdb(pdb_content)

// 遍历结构
for model in structure.models {
  for chain in model.chains {
    for residue in chain.residues {
      for atom in residue.atoms {
        atom.get_coord()        // → Vector3
        atom.distance(other)    // → Double
      }
    }
  }
}
```

### 8. FASTA 快速索引访问 (pyfaidx)

```moonbit
// 从内容创建索引
let fa = Fasta::from_content(fasta_content)

// 获取完整序列
let seq = fa.get_seq("chr1")       // → Seq?

// 快速随机访问子序列 (0-based, [start, end))
let sub = fa.fetch("chr1", 1000, 2000)?  // → Seq?

// 获取序列长度
let len = fa.get_length("chr1")    // → Int?

// 检查序列是否存在
fa.contains("chr1")                // → Bool

// 获取所有序列名称
fa.get_names()                     // → Array[String]

// 构建并写入 .fai 索引
let idx = build_fai(fasta_content)
let fai_str = write_fai(idx)

// 从 .fai 索引创建 Fasta
let fa = Fasta::new(fasta_content, fai_str)?
```

### 9. 机器学习特征提取

```moonbit
// DNA 特征提取
let dna_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

// k-mer 频率
let kmer_freq = kmer_frequency(dna_seq, 3, true)

// 规范 k-mer (考虑反向互补)
let canonical = dna_canonical_kmer_frequency(dna_seq, 3, true)

// 核苷酸组成
let comp = nucleotide_composition(dna_seq)
// comp.a, comp.t, comp.g, comp.c

// DNA 特征向量 (42维)
let dna_features = dna_feature_vector(dna_seq)
let dna_names = dna_feature_names()

// 蛋白质特征提取
let prot_seq = "MKKLLLISVLLFLSSAYSR"

// 氨基酸组成
let aa_comp = amino_acid_composition(prot_seq)

// 二肽/三肽组成
let dipep = dipeptide_composition(prot_seq)
let tripep = tripeptide_composition(prot_seq)

// 理化性质
let avg_h = avg_hydrophobicity(prot_seq)  // Kyte-Doolittle 疏水性
let avg_c = avg_charge(prot_seq)           // 电荷
let avg_p = avg_polarity(prot_seq)         // 极性
let mw = molecular_weight(prot_seq)         // 分子量

// 二级结构倾向 (Chou-Fasman)
let (helix, sheet, coil) = secondary_structure_propensity(prot_seq)

// 位置特异性特征
let pos_feat = position_specific_features(prot_seq, 5)

// 蛋白质特征向量 (73维)
let prot_features = protein_feature_vector(prot_seq)
let prot_names = protein_feature_names()
```

### 10. Biostrings 序列分析 (Bioconductor Biostrings)

```moonbit
// IUPAC 核苷酸频率
let freq = letter_frequency_dna("ACGTN")
// freq.a, freq.c, freq.g, freq.t, freq.n, freq.other

// IUPAC 氨基酸频率
let aa_freq = letter_frequency_aa("ACDEFGHIKLMNPQRSTVWY")

// k-mer 频率 (mono-, di-, tri-nucleotide)
let kmer_freq = oligonucleotide_frequency("AATCG", 2, false)

// 相对同义密码子使用度 (RSCU)
let rscu = rscu("ATGTGCTGAATGAA")

// GC 含量 (滑动窗口)
let gc_window = gc_content_by_position("AGCT", 1)

// Shannon 熵
let entropy = shannon_entropy("AGCT")

// DUST 复杂度
let dust_score = dust_complexity(seq)

// IUPAC 模式匹配
let matches = iupac_match("N", "A")

// 最近邻法熔解温度计算
let tm = tm_nn("AGCTAGCT")

// 序列比对编辑距离
let edits = nedit_at("AGCT", "AGTT", true)

// 序列一致性百分比
let pid = percent_identity("AGCT", "AGTT")

// IUPAC 反向互补
let rc = iupac_reverse_complement("ACGTMRWSYKVHDBN")
```

### 11. GenomicRanges 基因组区间操作 (Bioconductor GenomicRanges)

```moonbit
// 创建 GRanges
let gr = granges(
  ["chr1", "chr1", "chr2"],
  [(1, 10), (5, 15), (20, 30)],
  [strand_plus(), strand_minus(), strand_star()]
)

// 区间操作
let shifted = granges_shift(gr, 5)           // 偏移
let narrowed = granges_narrow(gr, 3, 8)      // 缩小范围
let resized = granges_resize(gr, 15, "center")  // 调整宽度
let flanked = granges_flank(gr, 3, true, false) // 侧翼区域

// 集合操作
let reduced = granges_reduce(gr, 0)         // 合并重叠区间
let disjoined = granges_disjoin(gr)         // 分割区间
let union_result = granges_union(gr1, gr2)  // 并集
let intersect_result = granges_intersect(gr1, gr2)  // 交集
let diff_result = granges_setdiff(gr1, gr2) // 差集

// 重叠检测
let counts = count_overlaps(query, subject)  // 计数重叠数
let overlaps = find_overlaps(query, subject) // 查找重叠对

// 距离计算
let distances = granges_distance(gr1, gr2)  // 区间距离
let indices = nearest(gr1, gr2)             // 最近邻索引

// Strand 转换
let strand_char = strand_to_char(strand_plus())  // → "+"
```

### 12. DESeq2 差异表达分析 (Bioconductor DESeq2)

```moonbit
// 创建 DESeqDataSet (原始计数矩阵)
let counts = [
  [100, 120, 200, 220],
  [50, 60, 55, 45],
  [10, 15, 80, 90]
]
let row_names = ["GeneA", "GeneB", "GeneC"]
let col_names = ["Control1", "Control2", "Treat1", "Treat2"]
let design = [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]]

let dds = deseq_dataset(counts, row_names, col_names, design)

// 运行差异表达分析
let results = deseq(dds)  // 或 results(dds)

// 访问结果
results.base_mean           // 标准化后的平均表达量
results.log2_fold_change    // log2 倍数变化
results.p_value             // p 值

// 筛选显著差异基因
let sig_genes = significant_genes(results, 0.05, 1.0)
```

### 13. Suffix Array & Suffix Tree (libdivsufsort)

```moonbit
// 创建 Suffix Array
let sa = SuffixArray::new("banana")

// 模式匹配
sa.contains("ana")           // → true
sa.count("ana")              // → 2
sa.locate("ana")             // → [1, 3]

// 最长重复子串
sa.longest_repeated_substring()  // → "ana"

// LCP 数组
let lcp = LCPArray::new(sa)
lcp.get(0)                   // → LCP 值

// 创建 Suffix Tree
let st = SuffixTree::new("banana")
st.contains("ban")           // → true
st.locate("ana")             // → [1, 3]
st.longest_repeated_substring()  // → "ana"
```

### 14. Overlap-Layout-Consensus 序列组装 (Celera Assembler)

```moonbit
// 简单序列组装
let reads = ["ABCDEF", "DEFXYZ", "XYZ123"]
let assembly = olc_assemble(reads, 2)  // → "ABCDEFXYZ123"

// DNA 序列组装
let dna_reads = [
  "ACGTACGTCCG",
  "CGTCCGATGCA",
  "ATGCATGCTGA"
]
let dna_assembly = olc_assemble(dna_reads, 5)

// 重叠检测
let overlap = compute_overlap("ABCDEF", "DEFXYZ", 2)  // → 3

// 计算所有重叠关系
let overlaps = compute_all_overlaps(reads, 2)

// 构建重叠图
let graph = compute_overlap_graph(reads, 2)

// 生成 Graphviz 可视化
let dot = overlap_graph_to_dot(reads, 2)
```

### 15. Hidden Markov Model 基因预测 (HMMER / Augustus)

```moonbit
// 创建基因预测 HMM
let hmm = create_gene_prediction_hmm()

// 前向算法 - 计算观测序列概率
let prob = hmm.forward(["A", "T", "C", "G"])

// 后向算法 - 计算后向概率
let beta = hmm.backward(["A", "T", "C", "G"])

// 维特比算法 - 最可能状态路径
let path = hmm.viterbi(["A", "T", "C", "G"])

// 基因预测
let dna = "AAATTTCCCGGG"
let prediction = predict_genes(hmm, dna)

// 提取外显子
let exons = extract_exons(prediction, dna)

// Baum-Welch 训练
let trained = hmm.baum_welch(["A", "T", "C", "G"], 10, 1e-6)
```

### 16. TxDb 转录本数据库 (Bioconductor GenomicFeatures)

```moonbit
// 解析 GTF 文件
let gtf_content = "chr1\tensembl\tgene\t100\t500\t.\t+\t.\tgene_id \"GeneA\";\n..."
let txdb = parse_gtf(gtf_content)

// 获取基因/转录本/外显子/CDS 作为 GRanges
let gene_gr = genes(txdb)
let tx_gr = transcripts(txdb)
let exon_gr = exons(txdb)
let cds_gr = cds(txdb)

// 获取启动子区域 (200bp upstream, 50bp downstream)
let prom_gr = promoters(txdb, 200, 50)

// 获取内含子 (按转录本分组)
let introns_map = introns_by_transcript(txdb)

// 获取 5' UTR 和 3' UTR
let five_prime_utr = five_prime_utrs(txdb)
let three_prime_utr = three_prime_utrs(txdb)

// 获取 ID 列表
let gene_ids = gene_ids(txdb)
let tx_ids = transcript_ids(txdb)
let exon_ids = exon_ids(txdb)

// 按基因/转录本分组外显子
let exons_by_gene = exons_by_gene(txdb)
let exons_by_tx = exons_by_transcript(txdb)
let tx_by_gene = transcripts_by_gene(txdb)

// 显示 TxDb 摘要
let summary = txdb_summary(txdb)
```

### 17. ProtParam 蛋白质参数分析 (Bio.SeqUtils.ProtParam)

```moonbit
// 创建蛋白质分析对象
let protein = ProteinAnalysis::new("MKKLLLISVLLFLSSAYSRGVVVDQQCGGNIFRPEQLVSGSEIHARLGVLGSGGGFRLVAVQ")

// 基础属性
protein.length()              // → 序列长度
protein.molecular_weight()    // → 分子量 (考虑脱水)
protein.count_amino_acids()   // → 氨基酸计数 Map

// 氨基酸组成百分比
let percent = protein.get_amino_acids_percent()
percent.get('A')              // → Some(百分比)

// 理化性质
protein.gravy()               // → GRAVY 疏水性评分
protein.aromaticity()         // → 芳香性
protein.instability_index()   // → 不稳定指数 (<40 稳定, >40 不稳定)

// 等电点计算
protein.isoelectric_point()   // → 等电点 pI
protein.charge_at_ph(7.0)     // → pH=7.0 时的电荷

// 二级结构倾向 (Chou-Fasman)
let (helix, sheet, coil) = protein.secondary_structure()
// helix: α-螺旋倾向, sheet: β-折叠倾向, coil: 无规卷曲倾向

// 信号肽预测
let sp_result = protein.predict_signal_peptide()
sp_result.has_signal_peptide  // → 是否有信号肽
sp_result.cleavage_site       // → 剪切位点索引
sp_result.score               // → 预测分数

// 氨基酸组成
let comp = protein.amino_acid_composition()
comp.get('A')                 // → Some(0.5)
```

### 18. rtracklayer 基因组轨道格式 (Bioconductor rtracklayer)

```moonbit
// 解析 BED 格式 (3-12列)
let bed_content = "chr1\t100\t200\tgene1\t500\t+\nchr2\t300\t400\tgene2\t700\t-"
let bed_records = parse_bed(bed_content)

// BED 转 GRanges
let gr_from_bed = bed_to_granges(bed_records)

// 写入 BED 格式
let bed_output = write_bed(bed_records)

// 解析 WIG 格式 (variableStep)
let wig_content = "variableStep chrom=chr1 span=5\n100\t1.5\n105\t2.0\n110\t1.8"
let wig_records = parse_wig(wig_content)

// 解析 WIG 格式 (fixedStep)
let wig_fixed = "fixedStep chrom=chr1 start=100 step=10\n1.5\n2.0\n1.8"
let wig_fixed_records = parse_wig(wig_fixed)

// 写入 WIG 格式
let wig_output = write_wig(wig_records)

// 解析 BEDGraph 格式
let bg_content = "track type=bedGraph\nchr1\t100\t200\t1.5\nchr1\t200\t300\t2.0"
let bg_records = parse_bedgraph(bg_content)

// 写入 BEDGraph 格式
let bg_output = write_bedgraph(bg_records)

// 解析 GFF/GTF 格式
let gff_content = "##gff-version 3\nchr1\tensembl\tgene\t100\t200\t5.5\t+\t0\tID=gene1;Name=GeneA"
let gff_records = parse_gff(gff_content)

// GFF 转 GRanges
let gr_from_gff = gff_to_granges(gff_records)

// 写入 GFF 格式
let gff_output = write_gff(gff_records)
```

### 19. K-means 聚类分析 (scikit-learn)

```moonbit
// 创建 K-means 模型
let kmeans = KMeans::new(3, 100, 1e-6)

// 训练模型 (K-means++ 初始化)
let data = [
  [1.0, 2.0], [2.0, 1.0], [8.0, 7.0], [9.0, 8.0]
]
let fitted = kmeans.fit(data)

// 预测新数据
let predictions = fitted.predict([[1.5, 1.5], [8.5, 8.5]])

// 单个点预测
let cluster = fitted.predict_single([5.0, 5.5])

// 获取质心和标签
let centroids = fitted.get_centroids()
let labels = fitted.get_labels()

// 计算惯性
let inertia = fitted.inertia(data)

// 轮廓系数评估
let score = silhouette_score(data, labels)

// 寻找最优 k 值
let optimal_k = find_optimal_k(data, 2, 5)

// 基因表达聚类
let gene_expr = [[1.2, 3.4, 2.1], [8.5, 7.2, 6.8], ...]
let gene_kmeans = cluster_gene_expression(gene_expr, 3)
```

### 20. SearchIO 统一搜索结果 (Bio.SearchIO)

```moonbit
// 解析 HMMER3 tabular 格式
let hmmer3_tab = "# hmmscan 3.3.2\n" +
  "# query name           target name        accession   E-value  score  bias\n" +
  "PF00001.28            sp|Q9Y2W8|A1BG_HUMAN  Q9Y2W8     1e-50    200    5\n" +
  "PF00001.28            sp|P12345|ABC_HUMAN   P12345     1e-30    150    2\n"
let results = parse_hmmer3_tab(hmmer3_tab)
let qr = results[0]
qr.id                    // → "PF00001.28"
qr.hits.length()         // → 2

// 解析 BLAT PSL 格式
let psl_content = "28\t0\t0\t0\t0\t0\t0\t0\t+\tquery1\t1000\t0\t100\tchr1\t5000\t100\t200\t..."
let psl_results = parse_blat_psl(psl_content)

// 获取 top hits (按 E-value 排序)
let top5 = top_hits(qr, 5)

// 统计 HSP 数量
let total_hsps = count_hsps(qr)

// BLAST 结果转换为 SearchIO 模型
let blast_record = parse_blast_tab(blast_content)
let searchio_qr = blast_to_searchio(blast_record)
```

### 21. BLAST 结果解析 (Bio.Blast)

```moonbit
// 解析 BLAST tabular 格式
let blast_tab = "query1\tsubject1\t98.5\t100\t2\t0\t1\t100\t50\t150\t1e-50\t200\n" +
  "query1\tsubject2\t95.0\t90\t3\t0\t10\t99\t100\t189\t1e-30\t150\n"
let record = parse_blast_tab(blast_tab)
record.query_id          // → "query1"
record.hits.length()     // → 2

// 解析 BLAST XML 格式
let xml_record = parse_blast_xml(blast_xml_content)

// 过滤 hits (E-value < 0.01)
let filtered = filter_hits_by_evalue(record, 0.01)

// 过滤 hits (identity > 90%)
let high_identity = filter_hits_by_identity(record, 90.0)

// 获取最佳匹配
let best = best_hit(record)

// 获取最佳 HSP
let best_hsp = best_hsp(record)

// 获取所有 HSPs
let all_hsps = all_hsps(record)

// 获取查询序列长度
record.query_len         // → 查询序列长度
```

### 22. 替换矩阵 (Bio.SubsMat)

```moonbit
// 使用内置 BLOSUM62 矩阵
let blosum62 = get_matrix("BLOSUM62")

// 查询分数
let score = blosum62.score('A', 'A')       // → 4
let score_mismatch = blosum62.score('A', 'D')  // → -1

// 不区分大小写查询
let score_case = blosum62.score_case_insensitive('a', 'A')  // → 4

// 获取矩阵名称和尺寸
blosum62.name()          // → "BLOSUM62"
blosum62.size()          // → 20

// 获取所有氨基酸
let amino_acids = blosum62.amino_acids()

// 解析矩阵字符串
let matrix_str = "A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V\n" +
  "A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0\n"
let custom_matrix = parse_matrix(matrix_str)

// 支持的矩阵: BLOSUM62, BLOSUM45, PAM250, PAM30
```

### 23. 序列模体识别 (Bio.motifs)

```moonbit
// 创建 PWM (Position Weight Matrix)
let pwm = PWM::new([
  [0.8, 0.1, 0.1, 0.0],
  [0.1, 0.8, 0.1, 0.0],
  [0.0, 0.1, 0.8, 0.1],
  [0.1, 0.0, 0.1, 0.8]
])

// 计算序列得分
let score = pwm.score("ACGT")

// 获取共识序列
let consensus = pwm.consensus()             // → "ACGT"

// 获取最可能序列
let most_probable = pwm.most_probable()     // → "ACGT"

// 解析 MEME 格式
let meme_content = "MEME version 4.11.2\nALPHABET= ACGT\n...\n"
let motifs = parse_meme(meme_content)

// 搜索模体
let hits = search_motif(dna_seq, pwm, 0.8)
```

### 24. 限制性内切酶分析 (Bio.Restriction)

```moonbit
// 创建酶对象
let eco_r1 = RestrictionEnzyme::new("EcoRI", "GAATTC", [1], [1])
let bam_h1 = RestrictionEnzyme::new("BamHI", "GGATCC", [1], [1])

// 查找酶切位点
let sites = find_sites("GAATTCTGAATTC", eco_r1)  // → [0, 7]

// 酶切序列
let fragments = cut_sequence("GAATTCTGAATTC", eco_r1)
// → ["G", "AATTC", "T", "AATTC"]

// 获取酶信息
eco_r1.name              // → "EcoRI"
eco_r1.recognition_site  // → "GAATTC"
eco_r1.cut_positions()   // → [1]

// 常用酶列表
let enzymes = common_enzymes()
```

### 25. 序列聚类分析 (Bio.Cluster)

```moonbit
// 创建距离矩阵
let seqs = ["ACGT", "AGCT", "AAAA"]
let dist_matrix = distance_matrix(seqs, "identity")

// 层次聚类
let tree = hierarchical_clustering(dist_matrix, "average")

// 转换为 Newick 格式
let newick = tree_to_newick(tree)

// 轮廓系数评估
let labels = [0, 0, 1]
let score = silhouette_score(dist_matrix, labels)

// 基因表达聚类
let gene_expr = [[1.2, 3.4], [1.3, 3.5], [8.5, 7.2]]
let clusters = cluster_gene_expression(gene_expr, 2)
```

### 26. 群体遗传学 (Bio.PopGen)

```moonbit
// 计算等位基因频率
let genotypes = ["AA", "Aa", "aa", "AA", "Aa"]
let freq = allele_frequency(genotypes)
// freq.a, freq.A

// 哈迪-温伯格检验
let hw_result = hardy_weinberg_test(genotypes)
hw_result.p_value        // → p 值

// FST 计算
let pop1 = ["AA", "Aa", "aa"]
let pop2 = ["AA", "AA", "Aa"]
let fst = fst(pop1, pop2)

// 核苷酸多样性 (π)
let seqs = ["ACGT", "AGCT", "ATGT"]
let pi = nucleotide_diversity(seqs)

// Tajima's D
let tajima_d = tajima_d(seqs)
```

### 27. 密码子使用分析 (Bio.CodonUsage)

```moonbit
// 计算 CAI (密码子适应指数)
let cai = cai("ATGTGCTGAATGAA", "homo_sapiens")

// 计算 ENC (有效密码子数)
let enc = enc("ATGTGCTGAATGAA")

// 计算 RSCU (相对同义密码子使用度)
let rscu_result = rscu("ATGTGCTGAATGAA")

// GC3 含量
let gc3 = gc3_content("ATGTGCTGAATGAA")

// CBI (密码子偏好指数)
let cbi = cbi("ATGTGCTGAATGAA", "homo_sapiens")

// Fop (最优密码子频率)
let fop = fop("ATGTGCTGAATGAA", "homo_sapiens")

// 检测最优密码子
let optimal = detect_optimal_codons(sequences, "homo_sapiens")

// 密码子翻译
let aa = translate_codon("ATG")              // → 'M'
```

### 28. IRanges 整数区间操作 (Bioconductor IRanges)

```moonbit
// 创建 IRanges
let ir = iranges([1, 5, 20], [10, 15, 30])

// 区间操作
let shifted = iranges_shift(ir, 5)           // 偏移
let resized = iranges_resize(ir, 15, "start")  // 调整宽度
let narrowed = iranges_narrow(ir, 3, 8)      // 缩小范围

// 集合操作
let reduced = iranges_reduce(ir)             // 合并重叠区间
let disjoined = iranges_disjoin(ir)          // 分割区间
let union_result = iranges_union(ir1, ir2)   // 并集
let intersect_result = iranges_intersect(ir1, ir2)  // 交集
let diff_result = iranges_setdiff(ir1, ir2)  // 差集

// 重叠检测
let counts = count_overlaps(ir1, ir2)        // 计数重叠数
let overlaps = find_overlaps(ir1, ir2)       // 查找重叠对

// 距离计算
let distances = iranges_distance(ir1, ir2)   // 区间距离
let indices = nearest(ir1, ir2)              // 最近邻索引
```

### 29. AlignIO 比对格式解析 (Bio.AlignIO)

```moonbit
// 解析 ClustalW 格式
let clustal_content = "CLUSTAL W (1.83) multiple sequence alignment\n\n" +
  "Seq1          ACGTACGT\n" +
  "Seq2          AC--ACGT\n"
let align = parse_clustal(clustal_content)

// 解析 FASTA 比对格式
let fasta_align = parse_fasta_align(fasta_align_content)

// 解析 Stockholm 格式
let stockholm_content = "# STOCKHOLM 1.0\nSeq1 ACGTACGT\nSeq2 AC--ACGT\n//\n"
let stockholm_align = parse_stockholm(stockholm_content)

// 写入比对格式
let output = write_clustal(align)
```

### 30. TreeIO 进化树格式解析 (Bio.TreeIO)

```moonbit
// 解析 Newick 格式
let newick_tree = parse_newick("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")

// 解析 NHX 格式 (Newick + Extended)
let nhx_tree = parse_nhx("(A[&&NHX:conf=0.9]:0.1,B[&&NHX:conf=0.8]:0.2);")

// 树操作
newick_tree.count_terminals()   // → 终端节点数
newick_tree.distance("A", "B")  // → 节点间距离
newick_tree.common_ancestor(["A", "B"])

// 写入 Newick 格式
let newick_str = write_newick(newick_tree)
```

### 31. edgeR 差异表达分析 (Bioconductor edgeR)

```moonbit
// 创建 DGEList
let counts = [[100, 120, 200, 220], [50, 60, 55, 45]]
let dge = dge_list(counts, ["GeneA", "GeneB"], ["Ctrl1", "Ctrl2", "Treat1", "Treat2"])

// 计算归一化因子
let dge_norm = calc_norm_factors(dge)

// 拟合 GLM
let design = [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]]
let fit = glm_fit(dge_norm, design)

// 差异表达检验
let result = glm_lrt(fit, [0, 1])

// 筛选显著差异基因
let sig_genes = top_tags(result, 100)
```

### 32. limma 差异表达分析 (Bioconductor limma)

```moonbit
// 创建设计矩阵
let design = [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]]

// 拟合线性模型
let fit = lm_fit(gene_expr, design)

// 经验贝叶斯
let fit_eb = eBayes(fit)

// 差异表达检验
let contrast = [[0.0, 1.0]]
let result = contrasts.fit(fit_eb, contrast)

// voom 变换 (RNA-seq)
let voom_result = voom(dge_list, design)

// 获取 top 差异基因
let top_genes = top_table(result, n=100)
```

### 33. SummarizedExperiment 多维数据容器 (Bioconductor SummarizedExperiment)

```moonbit
// 创建 SummarizedExperiment
let assays = [[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]]
let row_data = [["GeneA"], ["GeneB"]]
let col_data = [["Sample1"], ["Sample2"]]
let se = summarized_experiment(assays, row_data, col_data)

// 获取数据
se.assays()              // → Assays 列表
se.rows()                // → 行数据
se.cols()                // → 列数据

// 行/列操作
let subset = se_subset(se, [0], [0])        // 子集
let merged = se_merge(se1, se2)             // 合并

// 维度信息
se.nrow()                // → 行数
se.ncol()                // → 列数
se.nassays()             // → Assays 数量
```

### 34. GenomicAlignments 基因组比对分析 (Bioconductor GenomicAlignments)

```moonbit
// 创建 GAlignments 对象
let galn = GAlignments::new(
  ["chr1", "chr1", "chr2"],
  [100, 200, 300],
  [1100, 1200, 1300],
  [Strand::Plus, Strand::Minus, Strand::Plus],
  ["read1", "read2", "read3"],
  [60, 40, 50]
)

// 计算覆盖度
let cov = coverage(galn)
cov.get("chr1")            // → Coverage 对象

// 按特征汇总比对数
let features = granges(["chr1"], [(100, 1200)], [Strand::Plus])
let summary = summarize_overlaps(galn, features, ["gene1"])
summary[0].count           // → 比对数
summary[0].coverage_fraction  // → 覆盖分数

// Pileup 操作
let pileup_result = pileup(galn, "chr1", 500, 600)

// 按 MAPQ 过滤
let filtered = filter_by_mapq(galn, 30)

// 按 strand 过滤
let plus_strand = filter_by_strand(galn, Strand::Plus)

// 转换为 GRanges
let gr = galn_to_granges(galn)
```

### 35. VariantAnnotation 变异注释 (Bioconductor VariantAnnotation)

```moonbit
// 创建 VCF 记录
let records = [VcfRecord {
  chrom: "chr1",
  pos: 100,
  id: ".",
  ref: "A",
  alt: ["T"],
  qual: 99.0,
  filter: ["PASS"],
  info: HashMap::new()
}]

// 检测变异类型
let vtype = detect_variant_type("A", "T")       // → VariantType::SNP
let vtype_ins = detect_variant_type("A", "AT")  // → VariantType::Insertion
let vtype_del = detect_variant_type("AT", "A")  // → VariantType::Deletion

// 变异定位
let loc = locate_variants(records, txdb, "coding")

// 编码效应预测
let predictions = predict_coding(records, txdb, reference)
for pred in predictions {
  pred.variant_type        // → 变异类型
  pred.effect              // → 编码效应 (synonymous/missense/nonsense)
  pred.aa_change           // → 氨基酸变化
}

// 变异汇总
let summary = variant_summary(records, txdb)
summary.total_variants     // → 总变异数
summary.snp_count          // → SNP 数量
summary.coding_variants    // → 编码区变异数
summary.gene_ids           // → 受影响基因
```

### 36. Affy Affymetrix芯片数据分析 (Biopython Bio.Affy)

```moonbit
// 创建探针集
let pm = [[100.0, 150.0], [200.0, 250.0]]
let mm = [[50.0, 75.0], [100.0, 125.0]]
let probe_set = ProbeSet::new("probe_set_1", pm, mm)

// 创建AffyBatch
let probe_sets = [probe_set]
let affy_batch = AffyBatch::new(["probe_set_1"], probe_sets, ["sample_1", "sample_2"])

// 背景校正
let bg_corrected = background_correct_rma(probe_set.pm_intensities)

// PM-MM差异计算
let diff = compute_pm_mm_difference(probe_set)

// 探针集汇总
let summarized = summarize_probeset(probe_set, "median")

// RMA标准化 (背景校正 + log2转换 + 分位数归一化 + 中位数汇总)
let normalized = rma_normalize(affy_batch)

// 批量汇总
let batch_summary = affy_batch_summarize(affy_batch, "mean")
```

### 37. SVDSuperimposer SVD蛋白质结构叠合 (Biopython Bio.PDB.SVDSuperimposer)

```moonbit
// 创建原子坐标
let atom1 = AtomCoordinate::new(1.0, 2.0, 3.0)
let atom2 = AtomCoordinate::new(4.0, 5.0, 6.0)

// 计算原子距离
let dist = atom1.distance(atom2)

// 创建原子集合
let atoms1 = [
  AtomCoordinate::new(1.0, 0.0, 0.0),
  AtomCoordinate::new(0.0, 1.0, 0.0),
  AtomCoordinate::new(0.0, 0.0, 1.0),
]
let atoms2 = [
  AtomCoordinate::new(2.0, 0.0, 0.0),
  AtomCoordinate::new(0.0, 2.0, 0.0),
  AtomCoordinate::new(0.0, 0.0, 2.0),
]

// SVD结构叠合
let superimposer = SVDSuperimposer::set(atoms1, atoms2)
superimposer.rmsd              // → RMSD值
superimposer.get_rotation()    // → 旋转矩阵
superimposer.get_translation() // → 平移向量

// 应用变换
let transformed = superimposer.apply(atoms1)

// 便捷函数: 直接叠合并获取结果
let (aligned, rmsd) = superimpose(atoms1, atoms2)

// 直接计算RMSD (不进行叠合)
let rmsd_direct = calculate_rmsd(atoms1, atoms2)
```

### 38. GOEnrichment GO功能富集分析 (Bioconductor GOstats/clusterProfiler)

```moonbit
// 创建 GO term
let go_term = GOTerm::new(
  "GO:0005515",
  "protein binding",
  "BP",
  ["GeneA", "GeneB", "GeneC"]
)

// 构建 GO 注释数据库
let go_terms = [
  GOTerm::new("GO:0005515", "protein binding", "BP", ["GeneA", "GeneB", "GeneD"]),
  GOTerm::new("GO:0008150", "biological process", "BP", ["GeneE", "GeneF"]),
]

// 差异表达基因列表
let de_genes = ["GeneA", "GeneB", "GeneC"]

// GO 富集分析 (超几何检验 + BH校正)
let results = go_enrich(de_genes, go_terms, 1000)

// 结果访问
results[0].go_id              // → GO ID
results[0].term_name          // → 术语名称
results[0].p_value            // → P值
results[0].adjusted_p_value   // → 校正后P值
results[0].odds_ratio         // → 优势比
results[0].gene_ratio        // → 基因比例
results[0].count             // → 富集基因数

// 过滤结果
let filtered = filter_enrichment(results, 0.05, 2)

// 按命名空间筛选
let bp_terms = go_terms_by_namespace(go_terms, "BP")

// 独立统计检验
let p_hyper = hypergeometric_test(3, 10, 90, 15)
let p_bonferroni = bonferroni_correction([0.01, 0.02, 0.05])
let p_bh = bh_correction([0.01, 0.02, 0.05])
```

### 39. SingleCell 单细胞数据分析 (Bioconductor SingleCellExperiment)

```moonbit
// 创建计数矩阵 (细胞 × 基因)
let counts = [
  [100.0, 200.0, 50.0],
  [120.0, 180.0, 60.0],
]
let gene_names = ["GeneA", "MT-CO1", "GeneB"]

// 计算 QC 指标
let qc = calculate_qc_metrics(counts, gene_names)
qc.n_umi            // → 每个细胞的 UMI 数
qc.n_genes          // → 每个细胞检测到的基因数
qc.mito_percent     // → 线粒体基因比例

// 细胞过滤
let filtered = filter_cells(counts, qc, 100, 10000, 200, 20.0)

// Log 标准化
let normalized = log_normalize(filtered, None)
let normalized_sf = log_normalize(filtered, Some([1.0, 0.8]))

// 寻找高变异基因
let top_genes = find_variable_genes(normalized, 2000)

// 按基因索引筛选
let selected = select_genes_by_indices(normalized, top_genes)

// PCA 降维 (幂迭代法)
let (pca_result, eigenvalues) = pca_power_iteration(selected, 50)

// 创建 SingleCellExperiment 对象
let sce = SingleCellExperiment::new(
  [normalized],
  [["Cell1"], ["Cell2"]],
  [["GeneA"], ["GeneB"], ["GeneC"]]
)
```

### 40. BAM 文件解析 (pysam)

```moonbit
// 解析 BAM 文件内容
let bam = parse_bam(bam_content)

// 访问比对记录
for record in bam.records {
  record.qname              // → 读取名
  record.flag               // → 标志位
  record.rname              // → 参考序列名
  record.pos                // → 比对位置 (1-based)
  record.mapq               // → 比对质量
  record.cigar              // → CIGAR 数组
  record.seq                // → 序列
  record.qual               // → 质量值
}

// 判断比对状态
record.is_paired()         // → 是否配对
record.is_proper_pair()    // → 是否正确配对
record.is_unmapped()       // → 是否未比对
record.is_reverse()        // → 是否反向互补
record.mate_is_unmapped()  // → mate 是否未比对

// 计算插入片段长度
record.tlen                // → 插入片段长度

// 获取标签
record.get_tag("NM")       // → 编辑距离
record.get_tag("MD")       // → MD 字符串
```

### 41. Bloom Filter & k-mer 计数 (Jellyfish / khmer)

```moonbit
// 创建 Bloom Filter (预估100万元素，误判率0.01)
let bf = BloomFilter::new(1_000_000, 0.01)

// 添加 k-mer
bf.add("ACGT")
bf.add("CGTA")

// 查询
bf.contains("ACGT")        // → true
bf.contains("TTTT")        // → false

// k-mer 精确计数
let kmer_counter = KmerCounter::new(21, 10_000_000)
kmer_counter.count("ACGTACGTACGTACGTACGTAC")
let count = kmer_counter.get("ACGTACGTACGTACGTACGTAC")

// 近似去重
let unique_count = kmer_counter.approx_unique_count()

// 获取所有 k-mer
let kmers = kmer_counter.get_kmers()
```

### 42. BWT + FM-index (Bowtie2 / BWA)

```moonbit
// 创建 BWT 索引
let bwt = BWT::new("banana")

// BWT 变换
let transformed = bwt.transform()  // → "annb$aa"

// 反向变换
let original = bwt.inverse(transformed)  // → "banana"

// 创建 FM-index
let fmi = FMIndex::new("banana")

// 模式匹配
fmi.contains("ana")        // → true
fmi.count("ana")           // → 2
fmi.locate("ana")          // → [1, 3]

// 精确模式匹配
let positions = fmi.search("nan")  // → [2]
```

### 43. De Bruijn Graph 序列组装 (SPAdes / Velvet)

```moonbit
// 创建 De Bruijn Graph (k=3)
let dbg = DeBruijnGraph::new(3)

// 添加 reads
dbg.add_read("ACGTCCG")
dbg.add_read("CCGATGC")

// 构建图
dbg.build()

// 简化图 (去除气泡)
dbg.simplify()

// 序列组装
let contigs = dbg.assemble()

// 获取所有节点
let nodes = dbg.get_nodes()

// 获取边
let edges = dbg.get_edges()
```

### 44. Smith-Waterman 局部序列比对

```moonbit
// 创建比对器
let sw = SmithWaterman::new()

// 设置打分矩阵
sw.set_matrix("BLOSUM62")
sw.set_gap_open(-11)
sw.set_gap_extend(-1)

// 蛋白质比对
let (alignment, score, start_pos) = sw.align("ACDEFG", "CDE")

// 获取比对结果
alignment.query           // → 查询序列
alignment.target          // → 目标序列
alignment.query_start     // → 查询起始位置
alignment.score           // → 比对分数

// DNA 比对
let sw_dna = SmithWaterman::new()
let result = sw_dna.align_nucleotide("ACGTACGT", "CGT")
```

### 45. Needleman-Wunsch 全局序列比对

```moonbit
// 创建比对器
let nw = NeedlemanWunsch::new()

// 设置打分参数
nw.set_match(1)
nw.set_mismatch(-1)
nw.set_gap_open(-2)
nw.set_gap_extend(-1)

// DNA 全局比对
let (alignment, score, start_pos) = nw.align("ACGTACGT", "AC--ACGT")

// 蛋白质全局比对
nw.set_matrix("BLOSUM62")
let prot_result = nw.align("ACDEFGH", "ADEFGH")

// 获取回溯矩阵
let backtrack = nw.get_backtrack_matrix()
```

### 46. dplyr 数据操作 (R dplyr)

```moonbit
// 创建 DataFrame
let df = DataFrame::new([
  ["GeneA", "100", "5.2"],
  ["GeneB", "200", "3.1"],
  ["GeneC", "150", "4.5"]
], ["gene", "count", "log2fc"])

// 过滤
let filtered = df.filter(|row| row.get("count").parse_int() > 120)

// 选择列
let selected = df.select(["gene", "log2fc"])

// 添加新列
let mutated = df.mutate("log10count", |row| row.get("count").parse_double().log10())

// 排序
let arranged = df.arrange("count", "desc")

// 分组汇总
let grouped = df.group_by("category").summarize("mean_count", |group| group.mean("count"))

// 连接
let merged = df.inner_join(other_df, "gene")
```

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

### 性能基准测试 (MoonBit WasmGC)

| 操作 | 序列长度 | 单次耗时 | 备注 |
|------|----------|----------|------|
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

### 性能基准测试 (SeqIO)

| 操作 | 数据规模 | 单次耗时 |
|------|----------|----------|
| parse_fasta | 100 records × 1000 bp | **1452 us** |
| write_fasta | 100 records × 1000 bp | **1695 us** |
| parse_fastq | 100 records × 1000 bp | **2497 us** |
| write_fastq | 100 records × 1000 bp | **2236 us** |

### 性能基准测试 (SeqUtils)

| 操作 | 序列长度 | 单次耗时 |
|------|----------|----------|
| gc (GC含量) | 10,000 | **92 us** |
| crc32 | 10,000 | **97 us** |
| seguid | 10,000 | **7370 us** |
| tm_wallace (熔解温度) | 10,000 | **81 us** |
| ProteinAnalysis::count_amino_acids | 99 aa | **24 us** |
| ProteinAnalysis::molecular_weight | 99 aa | **13 us** |

### 性能基准测试 (Bio.Phylo)

| 操作 | 规模 | 单次耗时 |
|------|------|----------|
| parse_newick | simple tree | **5.7 us** |
| Tree::count_terminals | simple tree | **0.61 us** |
| Tree::distance | simple tree | **1.85 us** |
| parse_newick | 50 leaves | **54.7 us** |

### 性能基准测试 (Bio.PDB)

| 操作 | 规模 | 单次耗时 |
|------|------|----------|
| parse_pdb | 300 atoms | **619 us** |
| write_pdb | 300 atoms | **1804 us** |
| Atom::distance | 2 atoms | **0.30 us** |

### 性能基准测试 (BAM)

| 操作 | 规模 | 单次耗时 |
|------|------|----------|
| parse_bam_base64 | 1 record | **7.19 us** |
| BamFile::to_sam | 1 record | **4.23 us** |
| BamFile::get_reference_names | 1 reference | **1.37 us** |
| BamFile::num_references | 1 reference | **0.19 us** |
| BamFile::len | 1 record | **0.28 us** |

### 性能基准测试 (CRAM)

| 操作 | 规模 | 单次耗时 |
|------|------|----------|
| is_cram_magic (positive) | 8 bytes | **0.26 us** |
| is_cram_magic (negative) | 8 bytes | **0.12 us** |

### Python 参考基准测试

> **注意**: BioPython 的 `complement`、`count`、`replace` 等操作使用 C 扩展实现，因此 MoonBit 在这些操作上较慢是预期的。MoonBit 的优势在于纯 WebAssembly 环境下的便携性和 `translate` 等算法密集型操作的性能。

| 操作 | 序列长度 | Python (BioPython) | MoonBit (WasmGC) | 相对性能 |
|------|----------|--------------------|------------------|----------|
| complement | 100,000 | 106 us | 836 us | ~8× slower |
| reverse_complement | 100,000 | 178 us | 556 us | ~3× slower |
| transcribe | 100,000 | 244 us | 485 us | ~2× slower |
| translate | 100,000 | 17.5 ms | 601 us | **~29× faster** |
| count | 100,000 | 98 us | 436 us | ~4× slower |
| find | 100,000 | 1.1 us | 0.65 us | **~1.7× faster** |
| replace | 100,000 | 234 us | 1453 us | ~6× slower |

## 测试验证

### 测试覆盖率

```
Total tests: 673
Passed: 673
Failed: 0
```

### 测试模块分布

| 模块 | 测试文件 | 测试数 |
|------|----------|--------|
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
| DESeq2 | `deseq2_test.mbt` | 3 |
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
| Sequence Motifs | `motifs_test.mbt` | 12 |
| Sequence Utils | `seq_utils_test.mbt` | 16 |
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

项目提供 40 个示例程序，展示各模块的典型用法：

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
| biostrings_demo | Biostrings 序列分析（IUPAC、RSCU、复杂度、Tm） | `moon run examples/biostrings_demo/main.mbt` |
| genomic_ranges_demo | GenomicRanges 基因组区间操作（GRanges、区间运算、集合操作） | `moon run examples/genomic_ranges_demo/main.mbt` |
| deseq2_demo | DESeq2 差异表达分析（数据集创建、结果分析、显著基因筛选） | `moon run examples/deseq2_demo/main.mbt` |
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
| motifs_demo | 序列模体识别（PWM构建、模体搜索、一致性序列、信息含量、MEME格式解析） | `moon run examples/motifs_demo/main.mbt` | 
| seq_utils_demo | 序列工具函数（分子量计算、GC含量、Tm值、序列类型检测、氨基酸转换、k-mer计数） | `moon run examples/seq_utils_demo/main.mbt` | 
| popgen_demo | 群体遗传学分析（等位基因频率、基因型频率、哈迪-温伯格检验、FST统计、Watterson's theta） | `moon run examples/popgen_demo/main.mbt` | 
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
- ✅ 实现 序列模体识别（PWM构建、模体搜索、一致性序列、信息含量、MEME格式解析）
- ✅ 实现 序列工具函数（分子量计算、GC含量、Tm值、序列类型检测、氨基酸转换、k-mer计数）
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
- [ ] 实现 CRAM 格式支持
- [ ] 添加 更多多样性指数计算
- [ ] 添加 SIMD 加速支持

## 许可证

Apache-2.0 License
