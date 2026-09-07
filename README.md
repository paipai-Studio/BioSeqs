# BioSeqs - MoonBit 生物信息学库

> GitHub: https://github.com/paipai-Studio/BioSeqs
>
> GitLink: https://gitlink.org.cn/IvanAXu/BioSeqs
>
> Mooncakes: https://mooncakes.io/docs/IvanAXu/BioSeqs

## 项目概述

BioSeqs 是一个基于 **MoonBit** 语言开发的生物信息学工具库，复刻主流生物信息学库（Biopython、Bioconductor、scikit-bio 等）的核心功能，并实现高效的序列组装算法。项目当前版本：`0.1.9`。

## 目录

- [Biopython 系列功能](#biopython-系列功能)
- [Bioconductor 系列功能](#bioconductor-系列功能)
- [Others 其他功能与算法](#others-其他功能与算法)
- [序列组装算法](#序列组装算法)
- [架构设计](#架构设计)
- [构建与测试](#构建与测试)
- [核心模块速查](#核心模块速查)

---

## Biopython 系列功能

### 序列处理与基础工具

| 功能模块 | 对应 Biopython | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **Seq / MutableSeq** | `Bio.Seq` | 序列对象、可变序列编辑、互补、转录、翻译、环状 DNA 操作 | ✅ |
| **SeqRecord** | `Bio.SeqRecord` | 带注释的序列记录（ID/description/annotations/features） | ✅ |
| **SeqFeature** | `Bio.SeqFeature` | 序列特征与位置系统（Exact/Before/After/OneOf/Within）、SimpleLocation/CompoundLocation、链方向、LocationParser 位点字符串解析、扩展修饰符 | ✅ |
| **SeqReference** | `Bio.Reference` / `Bio.Medline` | 文献引用管理、Medline/PubMed 解析、APA 格式化 | ✅ |
| **Bio.Alphabet** | `Bio.Alphabet` | IUPAC DNA/RNA/蛋白质字母表、简化字母表、空位字母表 | ✅ |
| **Bio.Data** | `Bio.Data` | IUPAC 数据、氨基酸映射、密码子表、互补碱基表 | ✅ |
| **SeqUtils** | `Bio.SeqUtils` | GC 含量（gc_fraction 支持 remove/ignore/weighted 三种歧义碱基模式）、GC/AT 滑动窗口偏斜、GC123 密码子位置 GC、nt_search IUPAC 模糊序列搜索、six_frame_translations 六框翻译与 GC 可视化、seq3/seq1 三字母与单字母转换（IUPAC 扩展码 + custom_map/undef_code）、分子量、Tm 值（Wallace/盐校正）、ORF 预测、序列相似度、Hamming/Levenshtein 距离、CheckSum（GCG/SEGUID） | ✅ |
| **SeqUtils 高级** | `Bio.SeqUtils` / `Bio.SeqUtils.ProtParam` / `MolWt` / `MeltingTemp` | 蛋白质参数（不稳定指数/GRAVY/等电点/信号肽/二级结构倾向）、分子量/消光系数/吸光度、Chou-Fasman 二级结构、IUPred 无序区、COILS 卷曲螺旋、Kolaskar 抗原性、Emini 可及性、Karplus-Schulz 柔性、ProtDao 无序预测、CircSeq 环状 DNA 酶切、**MeltingTemp 完整解链温度**（Tm_Wallace 经验法则、Tm_GC 八套经验公式、Tm_NN 近邻热力学 8 套 DNA/RNA/杂交参数表、盐校正方法 1-7 含 Mg2+/dNTP/Tris、错配与悬挂末端、自互补双链、DMSO/甲酰胺化学校正、自定义热力学表） | ✅ |
| **FreqAnalysis** | `Bio.FreqAnalysis` / `Bio.SeqUtils` | k-mer 计数、密码子使用频率、Shannon 熵、Wooton-Federhen 局部组成复杂度 (LCC)、语言学复杂度、DUST、CGR 混沌游戏表示、序列签名 | ✅ |
| **CodonUsage** | `Bio.SeqUtils.CodonUsage` / `Bio.codonalign` | CAI 密码子适应指数、RSCU 相对同义密码子使用、ENC 有效密码子数、GC3 偏斜、CBI/Fop、最优/稀有密码子检测、物种参考表 | ✅ |
| **Kmer** | `Bio.Kmer` | k-mer 计数与频率分析、Jaccard 相似度、Hamming 距离、k-mer 谱 | ✅ |
| **ApproxMatch** | `Bio.Seq.Approximate` | 近似字符串匹配（Levenshtein、错配/插入/缺失检测） | ✅ |

### 序列 I/O 格式

| 功能模块 | 对应 Biopython | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **SeqIO 统一接口** | `Bio.SeqIO` | seqio_parse / seqio_write / seqio_convert / seqio_read / seqio_to_dict | ✅ |
| **FASTA / FASTQ** | `Bio.SeqIO.FastaIO` / `QualityIO` | FASTA/FASTQ 解析与写入、质量编码处理 | ✅ |
| **GenBank / EMBL** | `Bio.SeqIO.GenBankIO` / `EmblIO` | GenBank/EMBL 严格解析与写入（LOCUS/FEATURES/ORIGIN、ID/AC/DE/SQ） | ✅ |
| **PIR / Tab** | `Bio.SeqIO.PdbIO (PIR/NBRF)` / `TabIO` | PIR/NBRF 蛋白/核酸格式、Tab 分隔（ID+序列） | ✅ |
| **InsdcIO** | `Bio.SeqIO.InsdcIO` | INSDC 国际核苷酸序列数据库协作格式 | ✅ |
| **UniProt XML** | `Bio.SeqIO.UniprotIO` | UniProt XML 解析、蛋白质条目、基因名、物种、功能注释、交叉引用 | ✅ |
| **SwissProt** | `Bio.SwissProt` | Swiss-Prot/UniProt 平面文本解析、特征注释、参考文献、关键词 | ✅ |
| **SeqXml** | `Bio.SeqIO.SeqXmlIO` | SeqXML 序列交换格式、XML 解析/序列化、DNA/RNA/Protein 类型、SeqRecord 互转 | ✅ |
| **2bit / Nib** | `Bio.SeqIO.TwoBitIO` / `NibIO` | 2 位碱基编码、N 块/掩码块、小/大端序列化、子序列提取、反向互补 | ✅ |
| **GFA** | `Bio.SeqIO.GfaIO` | GFA 基因组组装图：Segment/Link/Containment/Path、标签系统（Z/i/f）、图序列化、SeqRecord 转换 | ✅ |
| **XDNA / GCK** | `Bio.SeqIO.XdnaIO` / `GckIO` | Geneious XDNA 2 位编码大端序、GCK_FILE 魔数、特征/序列/拓扑解析、校验和 | ✅ |
| **IMGT / IgIO** | `Bio.SeqIO.ImgtIO` / `IgIO` | IMGT 免疫序列 13 字段管道分隔头部、IntelliGenetics/MASE 注释行、HLA 等位基因解析 | ✅ |
| **SFF** | `Bio.SeqIO.SffIO` | SFF 454/Roche 流图二进制格式、编码/解码、质量修剪、按名称查找 | ✅ |
| **SnapGene** | `Bio.SeqIO.SnapGeneIO` | SnapGene .dna 数据包结构（Comments/Features/Sequence/Primers/Notes）、JSON/XML 解析 | ✅ |
| **ABI / Phd** | `Bio.SeqIO.AbiIO` / `Bio.Sequencing.Phd` | ABI Sanger 测序 trace（ABIF 二进制、4 通道）、Phred 碱基识别（BEGIN_SEQUENCE/碱基-质量-峰位） | ✅ |
| **ACE** | `Bio.Sequencing.Ace` | ACE contig 组装格式：reads/contigs 解析、共有序列、覆盖度、GC 含量 | ✅ |
| **Smart 结构域** | `Bio.Smart` | SMART 蛋白质结构域数据库解析、SMART/Pfam/SignalP/Transmembrane 分类、E 值过滤、GO 注释 | ✅ |
| **ExPASy** | `Bio.ExPASy` | Swiss-Prot 条目解析、酶数据库查询、蛋白质参数计算、Cellosaurus 细胞系数据库平面文本解析 | ✅ |
| **UniGene** | `Bio.UniGene` | NCBI UniGene 固定宽度记录、序列/蛋白相似性/STS/转录本映射、SCOUNT 校验 | ✅ |
| **GOA** | `Bio.UniProt.GOA` | GAF 2.2 基因本体注释 17 列格式、按 GO ID/aspect/evidence/taxon 过滤、序列化往返 | ✅ |
| **GEO SOFT** | `Bio.Geo` | ^PLATFORM/^SAMPLE/^SERIES 实体、!属性行、#列定义、数据表、按类型/编号查询 | ✅ |
| **KEGG** | `Bio.KEGG` | KEGG 基因/通路/化合物/酶记录解析、通路分析、KGML pathway XML 解析与序列化 | ✅ |
| **KEGG REST** | `Bio.KEGG.REST` | REST 风格 KEGG API 客户端：kegg_info/list/find/get/conv/link 六种操作与对应 URL 构造器、多条目 "+" 分隔、化学式/exact_mass/mol_weight 检索选项、可注入 fetch 传输与 now 时钟、内置限速（每秒 3 次查询） | ✅ |
| **SCOP** | `Bio.SCOP` | cla/des/hie 三文件 SCOP 层级（root→class→fold→superfamily→family→protein→species→domain）、PDB 残基片段、SCCS 比较 | ✅ |
| **SCOP Raf / Dom** | `Bio.SCOP.Raf` / `Bio.SCOP.Dom` | ASTRAL RAF 序列图（0.01/0.02 版本头校验、7 列残基字段、3↔1 氨基酸映射 normalize_letters）、seqmap 解析/索引/index/slice/add/extend、scop_parse_residues 链与残基定义、RAF 域切割 get_seqmap、PDB ATOM/HETATM 残基级提取 get_atoms；SCOP Dom 废弃域文件制表符记录解析（sid/pdbid/fragments/hierarchy）、to_string 往返 | ✅ |
| **CAPS** | `Bio.CAPS` | 酶切扩增多态性序列标记：差异酶切位点检测、内切酶切/阻遏分类、多酶并行扫描 | ✅ |
| **GFF** | `Bio.GFF` | GFF3 格式解析、GFFFeature/GFFRecord、属性处理、特征提取 | ✅ |
| **FSSP** | `Bio.FSSP` | FSSP 结构比对数据库：HEADER/TITLE/COMPND、ALIGNMENTS、Z-score/RMSD/PID | ✅ |

### 序列比对

| 功能模块 | 对应 Biopython | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **PairwiseAligner** | `Bio.Align.PairwiseAligner` | Needleman-Wunsch 全局 / Smith-Waterman 局部统一比对、Gotoh 仿射空位罚分（open+extend）、独立 target/query gap 参数、DNA/蛋白质 BLOSUM62 打分、identity%/gaps 计数、target/query 定位 | ✅ |
| **Pairwise2** | `Bio.Pairwise2` | 灵活双序列比对、自定义打分函数 | ✅ |
| **NW / SW 独立实现** | `Bio.Align.PairwiseAligner` | Needleman-Wunsch / Smith-Waterman 动态规划、回溯矩阵、自定义打分 | ✅ |
| **AlignClusterer** | `Bio.Align.AlignClusterer` | 渐进式多序列比对（ClustalW 算法）、UPGMA 向导树、profile-profile 比对 | ✅ |
| **Substitution Matrices** | `Bio.SubsMat` / `Bio.Align.substitution_matrices` | BLOSUM62/45、PAM250/30、ArrayData/SubsMatrix 基础设施、矩阵注册表、频率矩阵、log-odds 打分、Shannon 熵、KL 散度、NCBI 矩阵解析 | ✅ |
| **MultipleSeqAlignment** | `Bio.Align` | 多序列比对对象与操作 | ✅ |
| **AlignInfo / AlignAbstract** | `Bio.Align.AlignInfo` / `AlignAbstract` | 一致性序列、保守位点、Shannon 熵、成对序列同一性、简约信息位点、同一性矩阵 | ✅ |
| **CodonAlign** | `Bio.codonalign` | 密码子替换分类、dN/dS Nei-Gojobori、Z-test 选择检验、Fisher 中性检验、滑窗 dN/dS、BH-FDR、成对 Ka/Ks 表 | ✅ |
| **AlignIO** | `Bio.AlignIO` | 比对文件 I/O（ClustalW/FASTA/Stockholm/PHYLIP） | ✅ |
| **Clustal 现代 API** | `Bio.Align.clustal` | Clustal metadata 与 interleaved block 严格读写、累计残基数、consensus、坐标路径、跨行投影、统计、canonical 写回 | ✅ |
| **PHYLIP 现代 API** | `Bio.Align.phylip` | sequential/interleaved 自动识别、固定宽度名称、坐标统计、canonical writer | ✅ |
| **共享参考比对合并** | `Bio.Align.Alignment.from_alignments_with_same_reference` | 混合 PWA/MSA 输入、首端/内部/末端 insertion 同步、多 query 投影、局部坐标与 metadata 保留 | ✅ |
| **Alignment 坐标组合** | `Bio.Align.Alignment.map/mapall` | alignment path 组合、local clipping、gap/正反链传播、双向坐标查询、PSL、1:1/1:3 codon-aware MSA 投影 | ✅ |
| **Alignment 详细计数评分** | `Bio.Align.Alignment.counts` | 左/内部/右 insertion/deletion、gap open/extend、identity/mismatch/positive、wildcard、替换矩阵与十二类 affine gap 总分 | ✅ |
| **Alignment 频率与替换统计** | `Bio.Align.MultipleSeqAlignment.substitutions` / `Alignment.frequencies` | 每列字母加权计数（`'-'` gap 行）、列归一化频率、全比对背景频率、对称残基替换计数矩阵（序列权重、M=(M+Mᵀ)/2）、观测对频率 q_ij、BLOSUM 式背景 p_i、q_ij/e_ij 相对频率（log-odds 替换矩阵估计） | ✅ |
| **PSL / PSLX** | `Bio.Align.psl` | 21/23 列严格读写、核酸与 translated DNA-protein 3:1 路径、正反链坐标、block/gap 统计、sequence-aware recount、坐标映射 | ✅ |
| **SAM 感知比对** | `Bio.Align.sam` | SAM header 与 typed tag 严格读写、显式 CIGAR path、soft/hard clipping、反向链序列与 PHRED、MD/NM 坐标映射 | ✅ |
| **A2M 状态感知 MSA** | `Bio.Align.a2m` | match/insertion 列状态、大小写与点语义、严格读写、坐标映射、插入槽、pair counts、共识与 match-only 投影 | ✅ |
| **EMBOSS 输出** | `Bio.Align.emboss` | srspair/pair/simple 报告、多 alignment 与多序列、局部/反向坐标、consensus 统计、compact path、canonical 往返 | ✅ |
| **Exonerate 输出** | `Bio.Align.exonerate` | cigar/vulgar 严格读写、完整 M/5/I/3/C/G/N/S/F 操作路径、正反链与 protein strand、3:1 translated 坐标、双向映射 | ✅ |
| **GCG MSF** | `Bio.Align.msf` | AA/NA/PileUp 严格解析、interleaved rows、标准 GCG checksum、gap 规范化、坐标路径、统计、canonical writer | ✅ |
| **NEXUS** | `Bio.Align.nexus` | DATA/CHARACTERS/TAXA block、nested comments、quoted taxa、sequential/interleaved MATRIX、MATCHCHAR、坐标统计 | ✅ |
| **Stockholm** | `Bio.Align.stockholm` | 多记录严格读写、GF/GS/GR/GC、reference/database reference、insertion/deletion 列、all-gap 压缩、坐标统计 | ✅ |
| **UCSC Chain** | `Bio.Align.chain` | 12/13 字段严格读写、连续多记录、正反 target/query 绝对坐标路径、size/dt/dq 块、双向位置/区间映射、反转 | ✅ |
| **MAF / index** | `Bio.Align.maf` | MAF track/header 与 a/s/i/e/q 严格读写、正负链绝对坐标、任意 component 映射、MafIndex 半开查询、多外显子拼接 | ✅ |
| **Mauve XMFA** | `Bio.Align.mauve` | `#SequenceN*` metadata、LCB、正负链坐标、区间索引、跨序列投影、统计、序列重建 | ✅ |
| **BED pairwise** | `Bio.Align.bed` | BED3-BED12 严格读写、正负链 target/query 路径、exon block 重建、双向 residue 映射、半开区间查询、分级 writer | ✅ |
| **BigBed / BigMaf / BigPsl** | `Bio.Align.bigbed` / `bigmaf` / `bigpsl` | BigBed v4 二进制、AutoSQL 扩展、多级 B+ tree/R-tree、zlib DEFLATE、区间/名称查询、bedMaf/bedPsl 语义 | ✅ |
| **Tabular 搜索轨迹** | `Bio.Align.tabular` | BLAST outfmt 7、FASTA 8CB/8CC、BTOP/aln_code traceback、链向与 translated 坐标 | ✅ |

### 搜索结果解析

| 功能模块 | 对应 Biopython | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **BLAST 基础** | `Bio.Blast` (legacy) | 历史 tabular/XML 标签解析、HSP 过滤与最佳匹配 | ✅ |
| **现代 BLAST XML** | `Bio.Blast` | XML1/XML2 严格解析与规范写回、多 query/report、parameters/statistics、description/taxonomy、八类 BLAST 程序的链向/translated HSP 坐标路径 | ✅ |
| **BLAST Applications** | `Bio.Blast.Applications` | 8 种 BLAST 变体命令行包装、快速构建器、参数管理、命令构建与验证 | ✅ |
| **SearchIO 统一模型** | `Bio.SearchIO` | QueryResult / Hit / HSP / Fragment 四层结构、E-value 过滤、BLAST 转换 | ✅ |
| **HMMER3** | `Bio.SearchIO` (HmmerIO) | domtblout 域表、文本格式、domain 聚合、非 verbose 文本与 --noali | ✅ |
| **Infernal** | `Bio.SearchIO.InfernalIO` | cmscan/cmsearch tabular 1/2/3 自动检测、CM/HMM-only、正负链坐标、local-end 多片段 | ✅ |
| **BLAT PSL** | `Bio.SearchIO` (BLAT) | PSL 解析、SearchIO 转换 | ✅ |
| **FASTA 搜索** | `Bio.SearchIO.FastaIO` | -m8 紧凑表格、-m9 注释头、元数据提取 | ✅ |
| **HH-suite HHR** | `Bio.Align.hhr` | HHsearch/HHblits HHR 严格解析、profile-profile 比对、consensus/二级结构/DSSP/confidence 注释、命中筛选、坐标映射、序列化往返 | ✅ |
| **Exonerate vulgar/cigar** | `Bio.SearchIO.ExonerateIO` | vulgar 比对块三元组、cigar 格式、格式自动检测、分数过滤、内含子统计、字符串重建 | ✅ |
| **Exonerate C4 文本** | `Bio.SearchIO.ExonerateIO.exonerate_text` | C4 层次 Document/Query/Hit/HSP/Fragment、3/4/5 行模型、wrapped blocks、intron/NER/split codon/frameshift、翻译与坐标投影 | ✅ |
| **InterProScan** | `Bio.SearchIO.InterproscanIO` | TSV 14 列格式解析（蛋白 ID/数据库/签名/位置/分数/IPR/GO）、按数据库/蛋白过滤、GO 条目提取、按蛋白分组 | ✅ |
| **COMPASS** | `Bio.Compass` | COMPASS profile-profile 比对输出、多记录解析（SW 分数/E 值/百分比一致性/比对序列）、过滤与摘要 | ✅ |

### 系统发育树

| 功能模块 | 对应 Biopython | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **Phylo 基础** | `Bio.Phylo` | 树结构（Clade/Tree）、距离计算、可视化 | ✅ |
| **TreeIO** | `Bio.TreeIO` | Newick / NHX 格式解析与序列化、树操作、修剪 | ✅ |
| **TreeConstruction** | `Bio.Phylo.TreeConstruction` | 距离矩阵建树、UPGMA/WPGMA/NJ 算法、Jukes-Cantor/Kimura 替换模型 | ✅ |
| **PhyloXML** | `Bio.Phylo.PhyloXML` | PhyloXML 格式解析/序列化、Newick 双向转换、分类单元注释 | ✅ |
| **NeXML** | `Bio.Phylo.NeXML` | NeXML 格式解析与序列化、OTUs/Trees/Characters 数据模型、Newick 转换 | ✅ |
| **CDAO** | `Bio.Phylo.CDAO` | CDAO 本体 RDF/XML 格式、Tree/Node/TU/Edge 三元组建模、Newick 双向转换、命名空间处理 | ✅ |
| **Parsimony** | `Bio.Phylo.Parsimony` | Fitch 算法（状态集交集/并集）、Sankoff 算法（动态规划+代价矩阵） | ✅ |
| **Consensus** | `Bio.Phylo.Consensus` | 多数规则/严格一致性树、分裂分析、支持度计算 | ✅ |
| **Phylo Metrics** | `Bio.Phylo` (metrics) | Colless 平衡指数、Robinson-Foulds 距离、距离矩阵 | ✅ |
| **Trie** | `Bio.Phylo.Trie` | 前缀树数据结构、插入/查找/前缀匹配/最长前缀/子串搜索、限制酶位点检测、引物匹配 | ✅ |
| **PopGen 基础** | `Bio.PopGen` | 等位基因频率、FST、哈迪-温伯格检验 | ✅ |
| **PopGen 高级** | `Bio.PopGen` (advanced) | Tajima's D、Fu & Li's D/F、McDonald-Kreitman 检验、等位基因频率谱、中性综合分析 | ✅ |
| **GenePop** | `Bio.PopGen.GenePop` | GenePop 基因型解析、等位基因频率、杂合度、序列化往返 | ✅ |
| **Align.analysis** | `Bio.Align.analysis` | dn/ds Nei-Gojobori、进化距离 Jukes-Cantor/Kimura、选择压力分析 | ✅ |
| **PAML 兼容层** | BioSeqs API | Nei-Gojobori 风格 dN/dS、Jukes-Cantor 校正、密码子使用分析（不解析 PAML 控制文件） | ✅ |
| **PAML codeml** | `Bio.Phylo.PAML.codeml` | 严格 control 读写、CODONML/AAML、多 NSsites/branch-site/clade/free-ratio、pairwise、距离矩阵、BEB/NEB、AIC/BIC/LRT | ✅ |
| **PAML baseml** | `Bio.Phylo.PAML.baseml` | 严格 control 读写、JC69 至 UNRESTu、参数/SE、kappa、REV/UNREST Q 矩阵、rate-class、auto-dGamma、nhomo 节点、AIC/BIC/LRT | ✅ |
| **PAML yn00** | `Bio.Phylo.PAML.yn00` | 严格 control 读写、NG86/YN00/LWL85/LWL85m/LPB93、PAML 4.1-4.9i 兼容、对称矩阵与统计汇总 | ✅ |
| **Phylo Applications** | `Bio.Phylo.Applications` / `Bio.Phylo.Applications._Fasttree` | 系统发育命令行包装：FastTree（-nt/-n/-quote 核苷酸、-boot 支持度、-gtr/-gamma/-spr/-mlnni/-wag/-pseudo）、PhyML（-d 数据类型、-m 替换模型、-b bootstrap、-c/-a/-s SPR、-o tlr、--rand_start/--n_rand_starts/--r_seed）、RAxML（-s/-n/-m 分区模型、-f a 快速 bootstrap、-x/-N 重采样、-q 分区、-t/-r/-g 起始/约束树、-T 线程、-o 外群、-w 工作目录） | ✅ |

### PDB 结构分析

| 功能模块 | 对应 Biopython | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **PDB 数据类型** | `Bio.PDB` | 原子/残基/链/模型解析、结构操作 | ✅ |
| **PDB 结构 I/O** | `Bio.PDB.PDBParser` / `MMCIFParser` | PDB 文件 I/O、mmCIF 数据块/类别/原子位点提取 | ✅ |
| **mmCIF writer** | `Bio.PDB.mmcifio` | Structure 对象序列化、data block/header/atom_site loop 20 列格式化、HETATM 支持、round-trip 验证 | ✅ |
| **BinaryCIF** | `Bio.PDB.binary_cif` | 纯 MoonBit MessagePack 解析、ByteArray/FixedPoint/IntervalQuantization/RunLength/Delta/IntegerPacking/StringArray 七类逆编码、三态缺失值、Structure 转换 | ✅ |
| **MMTF** | `Bio.PDB.mmtf` | MMTF 二进制、IEEE 754 浮点、MsgPack 风格整数、Run-length/Delta/Recursive 编码、层级结构、键/实体/晶体学 | ✅ |
| **NACCESS** | `Bio.PDB.NACCESS` | NACCESS .rsa/.asa 解析、残基/原子级 ASA、绝对/相对/侧链/主链、链总和、Structure 互转 | ✅ |
| **DSSP** | `Bio.PDB.DSSP` | 二级结构预测、溶剂可及表面积 SASA、Ramachandran 图、结构质量评估 | ✅ |
| **SASA** | `Bio.PDB.SASA` | Shrake-Rupley 滚动球（Fibonacci 球面采样）、范德华半径查表、逐原子/残基/链 SASA、骨架/侧链拆分 | ✅ |
| **HSExposure** | `Bio.PDB.HSExposure` | 蛋白质残基半球暴露度计算 | ✅ |
| **Packing** | `Bio.PDB.Packing` | 局部包装密度、SASA Lee-Richards、范德华半径、低包装区域识别、球体点生成 | ✅ |
| **Dice + Selection** | `Bio.PDB.Dice` | PDB 结构切割（链/残基/原子/模型提取）、B 因子/occupancy 过滤、几何选择、结构统计、序列提取 | ✅ |
| **PDBList** | `Bio.PDB.PDBList` | PDB 结构下载管理、文件路径生成、过期 PDB 解析 | ✅ |
| **ParsePDBHeader** | `Bio.PDB.ParsePDBHeader` | PDB 头部元数据解析（HEADER/TITLE/COMPOUND/SOURCE/REMARK/AUTH/DBREF） | ✅ |
| **SVDSuperimposer** | `Bio.PDB.SVDSuperimposer` | SVD 蛋白质结构叠合、旋转矩阵、平移向量、RMSD | ✅ |
| **QCPSuperimposer** | `Bio.PDB.QCPSuperimposer` | 四元数特征多项式叠合、高精度旋转、RMSD | ✅ |
| **CEAligner** | `Bio.PDB.cealign` | CA/C4' 引导原子、AFP 路径搜索、CE Z-score、QCP 刚体叠合、全原子变换 | ✅ |
| **MAalign** | `Bio.PDB.MAalign` | 多蛋白质结构比对、Kabsch 算法迭代对齐、保守性分析 | ✅ |
| **StructureAlignment** | `Bio.PDB.StructureAlignment` | 多结构比对、动态规划、RMSD/TM-score、渐进式多结构比对 | ✅ |
| **ResidueDepth** | `Bio.PDB.ResidueDepth` | 残基深度、晶体学 B 因子、溶剂可及性、埋藏度 | ✅ |
| **PSEA** | `Bio.PDB.PSEA` | 二级结构预测（CA-CA 距离、虚拟键角/二面角、H/E/C 三态→八态） | ✅ |
| **FragmentMapper** | `Bio.PDB.FragmentMapper` | DSSP 二级结构分类 (H/E/S/L/T/C)、片段分配/合并/过滤、覆盖率 | ✅ |
| **internal_coords** | `Bio.PDB.internal_coords` | 键长/键角/扭矩角 (phi/psi/omega/chi)、二面角计算、内部→笛卡儿转换、扩展链构建、旋转异构体库、Ramachandran 区域 | ✅ |
| **Vectors** | `Bio.PDB.vectors` | 3D Vector3、RotationMatrix3、叉积、Kabsch 叠合、二面角 | ✅ |
| **chem_utils** | `Bio.PDB.chem_utils` | 范德华/共价半径、键长/键角/二面角、经验式/分子式量、氢键长度 | ✅ |
| **PDB 结构分析高级** | `Bio.PDB` (advanced) | 主链二面角、氢键检测、二级结构分配、疏水性分析 | ✅ |
| **Polypeptide** | `Bio.PDB.Polypeptide` | 氨基酸组成、疏水性、跨膜预测、等电点 | ✅ |
| **protein_analysis** | `Bio.protein_analysis` | Kyte-Doolittle 疏水性、GOR 二级结构、Hopp-Woods 抗原性、跨膜区段、AA/二肽/三肽组成、Shannon 熵保守性评分 | ✅ |

### 模体、数据库接口与其它

| 功能模块 | 对应 Biopython | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **Motifs 基础** | `Bio.motifs` | PWM、MEME 格式、模体搜索、per-position 信息含量、总信息含量、序列 Logo、模体富集、Pearson 相关性比较 | ✅ |
| **Motifs 高级** | `Bio.motifs` (advanced) | JASPAR PFM、TRANSFAC 格式、模体最优比对、KL/JS 散度、模体聚类 | ✅ |
| **AlignAce / MEME / MAST / Transfac** | `Bio.Motifs.*` | AlignAce PFM 解析、MEME PSPM/E 值/背景、MAST p 值/E 值、TRANSFAC PFM/AC/ID/DE/BF/CC、序列化往返 | ✅ |
| **Motif Scan (FIMO)** | FIMO 风格 | PWM→PSSM log2 似然比、正反向互补扫描、动态规划 p 值卷积、匹配汇总 | ✅ |
| **seqLogo** | Bioconductor seqLogo 风格 | PWM 构建、信息含量 IC、ASCII 艺术 logo 渲染、一致性序列、自定义背景频率 | ✅ |
| **Prosite** | `Bio.Prosite` | Prosite 模体数据库搜索、模式解析、匹配算法、得分计算 | ✅ |
| **Restriction** | `Bio.Restriction` | 内切酶数据库、酶切位点查找、序列酶切、片段分析 | ✅ |
| **ProtParam** | `Bio.SeqUtils.ProtParam` | 分子量、不稳定指数、GRAVY、等电点、信号肽、二级结构倾向 | ✅ |
| **Proteomics** | `Bio.SeqUtils.Proteomics` | 8 种蛋白酶切（胰酶/糜酶/胃酶/LysC/ArgC/CNBr/GluC/AspN）、同位素分布、b/y 碎片离子 | ✅ |
| **Entrez** | `Bio.Entrez` | NCBI 数据库 ESearch/EFetch、PubMed/Gene/Taxonomy 解析 | ✅ |
| **Entrez E-utilities** | `Bio.Entrez` | EPost（WebEnv/QueryKey 会话）、ESummary（DocSum/Item 文档摘要）、ESpell（拼写建议与 Replaced 列表）、ELink（跨库 LinkSetDb/LinkName 关联）、ECitMatch（引文字符串 bdata 检索 PMID）、entrez_quote_plus URL 编码、可注入 fetch/now 传输、内置限速（无 API key 每秒 3 次） | ✅ |
| **Datastore** | BagIt 风格（Bio.Entrez 缓存思路） | MD5（RFC 1321 完整实现）内容寻址 bag：add/manifest-md5.txt/bagit.txt/bag-info.txt（Payload-Oxum）、fetch-through 缓存（hit/miss 统计）、verify 校验、外部清单 load_bag 加载、字典序 names/tags | ✅ |
| **Taxonomy** | `Bio.Taxonomy` | NCBI 分类数据库、分类树操作、谱系查询、共同祖先计算 | ✅ |
| **Medline** | `Bio.Medline` | Medline/PubMed 记录解析、APA 引用、MeSH 过滤 | ✅ |
| **EMBOSS 工具** | EMBOSS suite | GC 偏斜、AT 偏斜、分子量、Tm、ORF 查找、距离计算、蛋白质参数 | ✅ |
| **Primer3** | `Bio.Emboss.Primer3` | Wallace/盐校正 Tm、GC 含量、自互补/交叉二聚体评分、发夹 3' 端检测、正/反向候选流水线 | ✅ |
| **FreqTable** | `Bio.SubsMat.FreqTable` | 计数/频率字典、read_count/read_freq 文本解析、Shannon 熵、KL 散度、Jensen-Shannon 距离、归一化 | ✅ |
| **Affy** | `Bio.Affy` | Affymetrix 芯片、RMA 标准化、背景校正、分位数归一化 | ✅ |
| **Graphics / GenomeDiagram** | `Bio.Graphics` / `Bio.Graphics.GenomeDiagram` | Track/Feature/Diagram 数据模型、Rectangle/Arrow/Diamond 形状、SVG 生成、自动标注 | ✅ |
| **Chromosome** | `Bio.Graphics.Chromosome` | 染色体/特征/区域/带、SVG 线性/圆形渲染、G 带染色、核型图 | ✅ |
| **Wise2** | `Bio.Wise` | GeneWise/ESTwise 输出解析、外显子/内含子/比对列、剪接相位、比特分数 | ✅ |
| **PCD 质谱** | `Bio.PCD` | 质谱 PCD 解析（Scan/RT/PEPMASS）、峰列表、TIC/BPC 色谱图、m/z 过滤、前体离子 | ✅ |
| **NMR** | `Bio.NMR` | NOE 距离约束（XPLOR/CNS）、二面角约束、化学位移表（BMRB-like）、NMRView .xpk、约束违反 | ✅ |
| **Crystal** | `Bio.Crystal` | 小分子 CIF 解析、单胞参数/空间群、原子/化学键、分数↔笛卡尔变换、单胞体积/密度 | ✅ |
| **RNA 二级结构** | `Bio.SeqUtils` (RNA) | Nussinov 动态规划算法、发夹环/内部环/凸出环/多环识别 | ✅ |
| **Pathway** | `Bio.Pathway` | 物种/反应/通路数据结构及分析函数 | ✅ |
| **phenotype** | `Bio.phenotype` | PlateRecord/WellRecord、logistic/Gompertz 生长曲线、CSV/JSON 解析 | ✅ |
| **Cluster 数值聚类** | `Bio.Cluster` (cluster.c) | C 聚类库核心：8 种距离（欧氏[平方均值]/City Block/Pearson/绝对 Pearson/非中心/绝对非中心/Spearman/Kendall）、clusterdistance、层次聚类 treecluster（SLINK 单连接/完全/平均/质心连接）+ Tree cut/sort、kcluster（k-means/k-medians EM）、kmedoids、clustercentroids（均值/中位数/medoid）、somcluster 自组织映射、PCA（Golub-Reinsch SVD）、distancematrix、Record（Eisen Cluster/TreeView 格式读写 + 全套方法） | ✅ |
| **Cluster 序列聚类** | `Bio.Cluster` | 距离矩阵、层次聚类、Newick 输出、轮廓系数 | ✅ |
| **Variation** | `Bio.Variation` | SNP、突变检测、氨基酸替换、BLOSUM62/Grantham 矩阵 | ✅ |
| **Application** | `Bio.Application` / `Align.Applications` | 命令行工具包装（ClustalW/Clustal Omega/Muscle/MAFFT）、参数管理 | ✅ |
| **Bloom Filter** | Jellyfish/khmer 风格 | k-mer 计数、成员查询、误判率估算、近似去重 | ✅ |
| **File 压缩** | `Bio.File` | 自动 gzip/bzip2 检测、透明压缩读写、文件操作接口 | ✅ |
| **Bio.NaiveBayes** | `Bio.NaiveBayes` | k-mer 频率朴素贝叶斯、Laplace 平滑、top-K 预测 | ✅ |
| **Bio.Markov** | `Bio.Markov` | 1/2/3 阶马尔可夫链训练、转移概率、序列对数概率、CpG 岛 log-odds 检测、稳态分布 | ✅ |
| **Bio.LogisticRegression** | `Bio.LogisticRegression` | 二分类、Newton-Raphson、Hessian、操纵子预测示例 | ✅ |
| **Bio.MaxEntropy** | `Bio.MaxEntropy` | IIS 训练、指示特征函数、softmax 归一化 | ✅ |
| **Bio.NeuralNetwork** | `Bio.NeuralNetwork` | 前馈网络、反向传播、sigmoid 激活、XOR/iris 示例 | ✅ |
| **Compound / ChemmineR** | `Bio.Compound` / ChemmineR | 分子式解析/分子量、SDF 解析、分子描述符、原子对指纹、Tanimoto/Dice 相似度、子结构搜索 | ✅ |
| **GA 遗传算法** | `Bio.GA` | 种群初始化、适应度函数、锦标赛/轮盘赌选择、单点/两点交叉、变异、精英保留 | ✅ |
| **Reduced AA** | `Bio.Align.Reduced` | RAD/Dayhoff/CHARM/SDM12 简化氨基酸字母表、序列比较、简化一致性 | ✅ |
| **Statistics 基础** | `Bio.Statistics` | 描述统计、假设检验（Pearson/Spearman 相关、t/Wilcoxon/Mann-Whitney/Fisher/KS/chi2/ANOVA）、Z-score、log-rank、BH/Yekutieli/Bonferroni/Holm 校正 | ✅ |
| **q-value / IHW** | qvalue / IHW | Storey's π₀ 估计（df=3 三次光滑样条 smoother / David Robinson 闭形式 bootstrap）、q-value、pFDR 有限样本校正、fdr.level 显著集、局部 FDR lfdr（probit/logit + 高斯 KDE）、经验 p 值 empPvals（池化/检验特异）、自助法；独立假设加权、协变量加权 Bonferroni、多协变量 | ✅ |
| **Nexus** | `Bio.Nexus` | NEXUS 格式解析、数据矩阵、发育树、距离矩阵 | ✅ |
| **Stockholm 兼容层** | `Bio.Stockholm` (legacy) | Stockholm/Pfam 比对解析、二级结构注释、百分比一致性、保守性分析 | ✅ |
| **MAF 兼容层** | `Bio.Align.MAF` (legacy) | 宽松 MAF 块解析、选择/过滤、百分比一致性、统计分析 | ✅ |
| **Mauve 兼容层** | `Bio.Align.Mauve` (legacy) | 历史 MAF-like 块、LCB 重排摘要、倒位/断点检测、覆盖率、BED 导出 | ✅ |

---

## Bioconductor 系列功能

### 核心基础设施与数据容器

| 功能模块 | 对应 Bioconductor | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **Biobase** | Biobase | ExpressionSet、AnnotatedDataFrame、数据归一化、log2 转换 | ✅ |
| **S4Vectors** | S4Vectors | Rle 游程编码、DataFrame 数据框、Hits 匹配数据结构 | ✅ |
| **BiocGenerics** | Bioconductor BiocGenerics | NA 处理、排序、集合运算、匹配、表统计 | ✅ |
| **BiocParallel** | BiocParallel | 任务分块、并行求和/均值、进度追踪 | ✅ |
| **IRanges** | IRanges | 整数区间操作、集合运算、重叠检测、findOverlaps 高级类型、nearest、coverage、距离矩阵 | ✅ |
| **GenomicRanges** | GenomicRanges | GRanges、复合特征、区间操作、集合运算、precede/follow、coverage、distance_to_nearest、gaps/reduce/disjoin/setdiff/promoters/trim | ✅ |
| **GRangesList** | GenomicRanges | 命名复合基因组特征、split/unlist/relist、逐组区间变换、并行集合运算、特征级重叠/最近邻/覆盖度 | ✅ |
| **plyranges** | plyranges | dplyr-like tidy verbs for GRanges: filter/mutate/select/arrange/rename/group_by+summarise/join_by、metadata 管理 | ✅ |
| **GenomicAlignments** | GenomicAlignments | GAlignments 对象、coverage 计算、summarizeOverlaps、pileup 操作 | ✅ |
| **GenomicFiles** | GenomicFiles | 分布式按区间扫描 BAM/BED/VCF、批量查询、归约、覆盖度 | ✅ |
| **SummarizedExperiment** | SummarizedExperiment | 多维基因组数据容器、Assays、行/列操作 | ✅ |
| **RangedSummarizedExperiment** | SummarizedExperiment | GRanges/GRangesList 行范围、复合特征精确重叠/最近邻、覆盖度、区间变换与协调子集 | ✅ |
| **TreeSummarizedExperiment** | TreeSummarizedExperiment | 行/列树、节点链接、树节点子集、祖先/后代查询、层级聚合 | ✅ |
| **SingleCellExperiment** | SingleCellExperiment | 多 assay 管理、降维（PCA/t-SNE/UMAP）、size factors、归一化 | ✅ |
| **SpatialExperiment** | SpatialExperiment | 空间转录组学数据结构、空间坐标管理、图像数据、spots 过滤 | ✅ |
| **MultiAssayExperiment** | MultiAssayExperiment | 多组学数据协调、实验协调、样本映射、跨实验子集 | ✅ |
| **RaggedExperiment** | RaggedExperiment | 参差突变数据（Missense/Nonsense/Frame_Shift/Splice_Site 等）、基因×样本稀疏矩阵、TMB 每 Mb 计算、按基因/样本/突变过滤、LoF 识别 | ✅ |
| **DelayedArray / DelayedMatrixStats** | DelayedArray | 懒加载操作、分块处理、row/col 统计（mean/var/sd/median/min/max/sum）、NA 处理、子集操作 | ✅ |
| **SparseArray / Matrix** | SparseArray / Matrix | N 维规范化 COO、重复坐标合并、R 列主序、切片/置换/绑定、稀疏算术与矩阵乘法；CSC/CSR 格式、矩阵运算、范数 | ✅ |
| **MatrixGenerics / beachmat** | MatrixGenerics / beachmat | rowMeans/colMeans、rowVars/colVars、rowMedians/colMedians、rowMad/rowCounts/rowAnys/rowAlls；列/行块处理、线性迭代器、子集/转置/绑定 | ✅ |
| **BiocNeighbors / BiocSingular** | BiocNeighbors / BiocSingular | KMKNN/Annoy/BruteForce 最近邻、欧几里得/曼哈顿/余弦；Exact/IRLBA/Randomized SVD | ✅ |
| **GenomeInfoDb** | GenomeInfoDb | 基因组构建、染色体信息、着丝粒位置、染色体臂、标准染色体筛选 | ✅ |

### 注释与基因组数据库

| 功能模块 | 对应 Bioconductor | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **AnnotationDbi** | AnnotationDbi | 通用注释数据库接口、ID 映射、GO/KEGG 注释管理 | ✅ |
| **AnnotationFilter** | AnnotationFilter | 染色体筛选、生物类型过滤、区域重叠检测、符号模式匹配 | ✅ |
| **AnnotationHub** | AnnotationHub | 中心化注释资源、按类型/提供者/基因组/关键词搜索 | ✅ |
| **ensembldb** | ensembldb | Ensembl 注释数据库、基因/转录本/外显子/CDS 检索、biotype 过滤 | ✅ |
| **TxDb / GenomicFeatures** | GenomicFeatures | GTF 解析、Gene/Transcript/Exon/CDS 提取、UTR/内含子计算、启动子提取、区域查询 | ✅ |
| **BSgenome** | BSgenome | 基因组序列数据库、染色体序列检索、子序列提取、链特异性基因提取 | ✅ |
| **biomaRt** | biomaRt | 基因 ID 映射、基因注释查询、批量查询、外部数据库映射 | ✅ |
| **GEOquery** | GEOquery | GEO 数据获取、Series Matrix、SOFT 解析、ExpressionSet 转换 | ✅ |
| **rtracklayer** | rtracklayer | BED/WIG/BEDGraph/GFF 解析与写入、GRanges 转换、Chain liftOver | ✅ |
| **UCSC Chain / liftOver** | rtracklayer | Chain 格式解析、基因组坐标 liftOver 转换、链段查找、染色体间映射 | ✅ |
| **rhdf5** | rhdf5 | HDF5 文件支持、数据集读写、组管理、属性操作 | ✅ |

### 差异表达分析

| 功能模块 | 对应 Bioconductor | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **DESeq2** | DESeq2 | size factors 归一化、分散度估计、负二项 GLM 拟合、Wald 检验、LFC 收缩、VST 方差稳定化变换、PCA 可视化 | ✅ |
| **apeglm** | apeglm | 负二项 GLM、自适应经验贝叶斯 Cauchy/Student-t 先验、MAP、Laplace 后验 SD/区间、FSR/FSOS/s-value | ✅ |
| **edgeR** | edgeR | DGEList、精确检验、GLM 拟合、准似然 F 检验、camera/roast 基因集检验 | ✅ |
| **limma** | limma | 线性模型拟合、经验贝叶斯、voom 变换、RPKM/CPM/quantile 归一化、ComBat/removeBatchEffect 批次校正、treat 严格检验 | ✅ |
| **variancePartition** | variancePartition | 重复测量线性混合模型、ML/REML 方差分量、BLUP、precision weights、dream contrast、Satterthwaite 检验、BH-FDR | ✅ |
| **dreamlet** | dreamlet | sample×cell-type pseudobulk、TMM、cell/sample/gene 过滤、logCPM、Poisson/voom precision weights、分 cell-type 重复测量模型、study-wide FDR | ✅ |
| **DEXSeq** | DEXSeq | 差异外显子使用、计数归一化、统计检验、结果过滤 | ✅ |
| **DRIMSeq** | DRIMSeq | Dirichlet-multinomial 模型、Wald 检验、比例计算、counts 过滤、BH-FDR、显著基因提取 | ✅ |
| **baySeq** | baySeq | 负二项模型、分散度估计、Gamma 先验、对数似然比、后验概率、MAP | ✅ |
| **NOISeq** | NOISeq | TMM/RPKM/上四分位归一化、NOISeqBio 噪声鲁棒差异 | ✅ |
| **glmGamPoi** | glmGamPoi | size factors 估计、伪批量聚合、IWLCS 迭代加权最小二乘、Wald 检验、BH-FDR | ✅ |
| **fishpond (Swish)** | fishpond | Mann-Whitney-Wilcoxon 秩和统计、置换检验、BH-FDR、log2FC 方向判定 | ✅ |
| **ALDEx2** | ALDEx2 | Dirichlet Monte Carlo 组成型差异丰度、六类 denominator、Welch/Wilcoxon、effect/overlap、Aitchison 距离 | ✅ |
| **DirichletMultinomial** | DirichletMultinomial | Dirichlet-multinomial 概率、有限混合 EM/BFGS 聚类、Laplace/AIC/BIC 选 K、生成式分组分类、ROC | ✅ |
| **ANCOMBC** | ANCOMBC | 采样比例估计、对数比偏差校正、参考特征选择、Welch t 检验、ANCOM W 统计量、BH-FDR | ✅ |
| **metagenomeSeq** | metagenomeSeq | 零膨胀模型、归一化、零膨胀概率计算、差异检验 | ✅ |
| **zinbwave** | zinbwave | 零膨胀负二项低维模型、cell/gene 协变量与 offset、latent factors、dispersion shrinkage、observational weights | ✅ |
| **DSS** | DSS | RNA-seq 差异表达 Wald 检验 + BH-FDR；差异甲基化 DML/DMR 检测 | ✅ |
| **tradeSeq** | tradeSeq 1.27.0 | gene × cell 计数合同、cell × lineage 拟时间/权重、多 lineage 负二项 GAM、惩罚 B-spline、五类 Wald 检验、Slingshot 直连 | ✅ |
| **stageR** | stageR | 两阶段检验（筛选+确认）、Simes 聚合、BH-FDR、Holm 步降、OFDR 控制、Dte/Dtu | ✅ |
| **EBSeq** | EBSeq | 经验贝叶斯差异表达、负二项模型、PPDE 后验概率、Laplace 近似边际似然、Stirling lgamma 近似 | ✅ |

### 单细胞分析

| 功能模块 | 对应 Bioconductor | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **scuttle 1.23.1** | scuttle | batch-aware median/MAD 异常值、subset per-feature QC、重叠 feature-set 聚合、精确无放回 count downsampling、batch/block coverage 等化 | ✅ |
| **scrapper** | scrapper | 批次感知 RNA QC、大小因子清洗与居中、log-normalization、LOWESS 方差趋势、HVG 选择、多因子 pseudo-bulk | ✅ |
| **scran** | scran | 单细胞归一化 sum_factors、SNN 图构建、Leiden 聚类、差异标志物 | ✅ |
| **bluster 1.23.0** | bluster | 多起点 K-means、精确 KNN 与加权 SNN 图、Louvain 聚类、two-step 量化聚类、Rand/ARI、silhouette、purity、RMSD、modularity、bootstrap 稳定性 | ✅ |
| **FlowSOM 2.21.0** | FlowSOM | 规则网格 SOM、KWSP/random/PCA 初始化、四类距离、分阶段 MST 邻域、meta-clustering 与自动 elbow、节点 MFI/CV/SD/MAD/阳性率、MAD outlier | ✅ |
| **clusterExperiment** | clusterExperiment | clusterMany 参数网格集成、makeConsensus 共聚类矩阵+一致性聚类、makeDendrogram 聚类树、mergeClusters 显著性合并（Welch t-test+BH-FDR）、RSEC 流水线 | ✅ |
| **SC3** | SC3 | PCA 降维、k-means 聚类、轮廓系数评估、间隙统计量、共识聚类 | ✅ |
| **ConsensusClusterPlus** | ConsensusClusterPlus | 癌症亚型识别、稳定性评分、最优聚类数选择 | ✅ |
| **apcluster 1.4.15** | apcluster | Affinity Propagation 亲和传播（Frey & Dueck 2007）：responsibility/availability 消息传递 + 阻尼迭代、exemplar 代表点识别、negDistMat（Minkowski^r/平方欧氏）/expSimMat（高斯 RBF）/linSimMat/corSimMat（Pearson/Spearman）相似度矩阵、中位数偏好、per-point 偏好、apclusterK 偏好二分定 K + 最近 exemplar 剪枝合并、确定性噪声 runs | ✅ |
| **Seurat** | Seurat | LogNormalize、高可变基因、PCA、图聚类、UMAP、差异标志物、FindIntegrationAnchors/IntegrateData 跨样本整合 | ✅ |
| **Batchelor** | Batchelor | rescaleBatches 缩放校正、mutual nearest neighbor、fastMNN 多批次校正、批次混合评分 | ✅ |
| **scDblFinder 1.27.6** | scDblFinder | top-variable 特征、library/PCA、随机或跨 cluster 人工 doublet、kNN/cxds 特征、迭代正则化 logistic 分类、capture 分层、homotypic 修正 | ✅ |
| **DropletUtils 1.33.0** | DropletUtils | Simple Good-Turing ambient profile、knee/inflection 曲线追踪、multinomial/Dirichlet-multinomial 概率、alpha 估计、Phipson–Smyth p 值、emptyDrops | ✅ |
| **decontX** | decontX | cluster-native/contaminant 多项式混合、Beta/Dirichlet 先验 EM、empty-droplet background、计数分解 | ✅ |
| **celda** | celda `celda_CG` | 细胞群与基因模块联合聚类、分层 Dirichlet-multinomial、collapsed likelihood、EM/Gibbs、多链、K/L 网格选择 | ✅ |
| **miloR** | miloR | 精确 KNN 图、精炼重叠邻域、邻域×样本计数、NB-GLM/Wald 检验、BH 与四种 graph spatial FDR | ✅ |
| **muscat 1.27.4** | muscat | gene × cell 严格合同、五类 cluster-sample 伪批量、任意满秩设计/多 contrast、负二项 IRLS、经验贝叶斯 dispersion、DS/DD 与局部/全局 BH-FDR | ✅ |
| **monocle3** | monocle3 | PCA/UMAP 降维、主图学习、拟时间排序、差异表达、分支点检测 | ✅ |
| **slingshot 2.21.0** | slingshot | MST 构建、主曲线拟合、拟时间计算、soft membership、协方差缩放距离、start/end 约束、omega forest、拟时间、新数据映射 | ✅ |
| **velociraptor** | velociraptor | 稳态线性回归（gamma/beta）、EM 动力学模型（alpha/beta/gamma）、velocity 向量、KNN 嵌入投影、根细胞识别 | ✅ |
| **SCENIC** | SCENIC | GENIE3 TF-target 共表达模块、Regulon 构建（权重剪枝/cisTarget motif 排名剪枝）、AUCell 活性评分、二值化阈值（MeanStd/KMeans2/Median）、细胞状态聚类 | ✅ |
| **infercnv** | infercnv | 染色体位置排序基因、参考细胞比较、log2FC 有界计算、金字塔权重基因组平滑、每细胞中位数中心化+噪声过滤、CNV 分数+肿瘤细胞预测 | ✅ |
| **scmap** | scmap | scmap-cluster 质心 Spearman 相关、scmap-cell k 近邻投票、阈值过滤、参考图谱投影 | ✅ |
| **SingleR 2.15.2** | SingleR | 严格 gene × sample/cell 模型、成对 classic markers、ties-aware Spearman、标签内相关分位数、fine-tuning、delta/MAD 剪枝、cluster 注释 | ✅ |
| **MAST 1.39.0** | MAST | 严格 feature × cell 模型、任意设计矩阵与 CDR、Bayesian logistic/Gaussian hurdle GLM、嵌套模型 LRT、经验贝叶斯方差收缩、边际 logFC | ✅ |
| **cyclone** | cyclone | 细胞周期评分、基因对比较、G1/S/G2/M 期相预测 | ✅ |
| **scater** | scater | QC 指标、细胞/基因过滤、CPM/log-CPM、HVG 检测、PCA 降维 | ✅ |
| **SCnorm** | SCnorm | 分位数回归、深度依赖偏差校正、基因特异性归一化 | ✅ |
| **ShortRead** | ShortRead | QA 统计、adapter 修剪、质量修剪、读长过滤、FastQC 报告 | ✅ |
| **dorothea** | dorothea | Regulon、VIPER 算法、置换检验、Z-score、Top TF | ✅ |
| **PROGENy** | PROGENy | L2 正则化 Ridge 回归、通路基因集权重矩阵、样本通路活性 | ✅ |
| **AUCell** | AUCell | AUC 计算、基因排序、min-max 归一化、细胞/基因集评分查询 | ✅ |
| **scDblFinder / decontX / celda / milo 等容器接入** | 所有上面模块 | 统一 SingleCellExperiment 不可变写回与跨模块集成 | ✅ |

### 空间转录组

| 功能模块 | 对应 Bioconductor | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **nnSVG** | nnSVG | nearest-neighbor Gaussian process、空间变异基因检验、gene-specific length scale、协变量设计、空间方差占比、BH-FDR | ✅ |
| **Banksy** | Banksy | H0 邻域均值与 H1+方位 harmonic、六类空间核、lambda 联合特征、分组标准化、PCA、多起点 k-means、标签平滑、参数扫描 | ✅ |
| **Voyager** | Voyager | kNN/distance-band/inverse-distance 空间权重（W/B/C/S 编码）、全局 Moran's I 与 Geary's c、局部 Moran's I (LISA)、局部 Geary's c、Getis-Ord Gi/Gi*、Lee's L、多元局部 Geary、经验变差函数与 spherical/exponential/gaussian 拟合、Moran correlogram | ✅ |
| **spicyR** | spicyR | 有序细胞类型对 cross-L 曲线、矩形窗口边界校正、图像级共定位统计、precision weights、重复受试者随机截距、条件对比 | ✅ |
| **lisaClust** | lisaClust | 每细胞多类型 local-K/centered local-L 曲线、Gaussian KDE 强度校正、矩形/凸包窗口、圆盘边界修正、确定性多起点 k-means、silhouette、区域富集 | ✅ |
| **SpatialDecon** | SpatialDecon | 背景感知加权 log-normal 非负回归、两阶段异常点重拟合、Hessian 不确定度、细胞丰度/比例/计数尺度、cell-type collapse、reverse deconvolution、负探针背景 | ✅ |
| **BayesSpace** | BayesSpace | t 分布混合模型、马尔可夫随机场 (MRF) 先验、EM 算法、六边形/方形网格邻居 | ✅ |

### 表观基因组与甲基化

| 功能模块 | 对应 Bioconductor | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **minfi** | minfi | Illumina 450K/EPIC DNA 甲基化分析、NOOB/Illumina/分位数/功能归一化、β/M 值计算、DMP/DMR 分析 | ✅ |
| **missMethyl** | missMethyl | Beta/M 值转换、limma t 检验、BH-FDR、探针-基因映射偏倚校正、GO 富集 | ✅ |
| **bsseq** | bsseq | 亚硫酸氢盐测序、BSmooth 平滑、DMR 检测、CpG 合并、甲基化率 | ✅ |
| **methylKit** | methylKit | 亚硫酸氢盐、甲基化胞嘧啶统计、覆盖率过滤/归一化、Fisher 精确检验差异甲基化、BH-FDR、DMR、BED 导出 | ✅ |
| **methylSeekR** | methylSeekR | 基因组分 tile 计算甲基化水平、UMR/LMR/PMD/FMR 分类、滑窗 PMD 检测、相邻区域合并 | ✅ |
| **bumphunter** | bumphunter | t 统计量、滑动均值平滑、候选 bump 识别、置换检验 p 值、BH-FDR | ✅ |
| **ChIPseeker** | ChIPseeker | 峰-TSS 距离、基因组特征分配（Promoter/5'UTR/3'UTR/Exon/Intron/Downstream/Distal Intergenic）、最近基因查找、peakOverlap、Venn/饼图可视化 | ✅ |
| **DiffBind** | DiffBind | ChIP-seq 差异结合、峰值重叠、共识峰、TMM 归一化、负二项检验、PCA/热图/火山/MA 图 | ✅ |
| **MACS2 peak_calling** | MACS2 风格 | 滑动窗口扫描、局部 lambda 背景估计、Poisson 显著性（不完全 Gamma）、峰合并、BH-FDR、fold enrichment | ✅ |
| **csaw** | csaw | ChIP-seq 窗口差异、滑动窗口计数、TMM 归一化、窗口过滤、负二项 GLM 检验、差异区域 | ✅ |
| **nucleR** | nucleR | 核小体定位、信号平滑、峰值检测、核小体 occupancy 计算、位置比较 | ✅ |
| **bamsignals** | bamsignals | ChIP-seq 信号提取、计数模式、RPM/RPKM 归一化、染色质状态分析 | ✅ |
| **ChromVAR** | chromVAR | 染色质变异、TF motif 富集、GC 偏差校正、细胞聚类、变异性分析 | ✅ |
| **DNAshapeR** | DNAshapeR | 二核苷酸查找表（MGW/Roll/ProT/HelT/EP）、k-mer 滑动窗口、A-tract 窄小沟、GC-rich 宽小沟、反向互补 | ✅ |
| **HMMcopy** | HMMcopy | GC 偏差校正 LOESS 分箱、Viterbi 解码、高斯发射、前向算法对数似然、片段合并、CN0/CN3 | ✅ |
| **SGSeq** | SGSeq | 剪接图构建（外显子/连接）、SE/A5SS/A3SS/MXE/RI 五类事件、PSI 量化、STAR SJ 格式解析 | ✅ |
| **StructuralVariantAnnotation** | StructuralVariantAnnotation | VCF BND 断裂端解析（4 种 ALT 格式）、SV 类型推断（DEL/DUP/INV/INS/TRA/BND）、伴侣配对、基因区域重叠注释 | ✅ |
| **HiC-DC+** | HiC-DC+ | 负二项 GLM 背景建模（IRLS）、z-score 显著性检验、BH-FDR、方向性指数 TAD、PCA A/B compartment | ✅ |
| **CNVkit** | CNVkit | CBS（循环二元分割）分割、拷贝数状态判定、log2 比率平滑、断点检测 | ✅ |
| **DMRcate** | DMRcate | 差异甲基化区域检测、limma 调节 t 统计量、指数空间核 exp(-d/λ)、CpG 聚类成 DMR | ✅ |

### 基因集 / 通路 / 富集

| 功能模块 | 对应 Bioconductor | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **clusterProfiler** | clusterProfiler | 功能富集统一框架、超几何检验、多重检验校正、结果过滤与可视化 | ✅ |
| **DOSE** | DOSE | 疾病本体富集、超几何检验、adjusted p-value | ✅ |
| **ReactomePA** | ReactomePA | Reactome 通路富集、按顶层术语分组、通路注释查询 | ✅ |
| **topGO** | topGO | 拓扑 GO 富集（elim/weight01 算法、Fisher 精确检验、GO 图） | ✅ |
| **enrichplot** | enrichplot | 富集结果可视化：dotplot/barplot/heatmap/cnetplot/enrichment map | ✅ |
| **fgsea** | fgsea | 快速基因集富集、置换检验、NES/ES、Leading Edge 基因、BH-FDR | ✅ |
| **GSVA** | GSVA | 单样本通路评分（ssGSEA/zscore/PLAGE）、富集分析、置换检验、enrichment map、phenotype correlation、survival analysis | ✅ |
| **GAGE** | GAGE | fold change、t 检验、BH-FDR 校正、配对检验 | ✅ |
| **SPIA** | SPIA | 信号通路影响分析、通路图扰动累积、超几何检验、Fisher 合并 p 值 | ✅ |
| **decoupleR** | decoupleR | WSum/WMean/Norm/ULM/MLM 方法、先验知识网络 (PKN)、调控子活性评分 | ✅ |
| **GSEABase** | GSEABase | GMT/GMX 格式解析、基因集管理与集合运算 | ✅ |
| **GOstats** | GOstats | GO 超几何检验富集分析、条件检验（移除子项基因）、log-sum-exp 数值稳定性、lgamma log-choose | ✅ |

### 微生物组与免疫

| 功能模块 | 对应 Bioconductor | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **phyloseq** | phyloseq | OTU 丰度计算、分类学过滤、相对丰度、稀疏化处理 | ✅ |
| **microbiome** | microbiome | Alpha 多样性（Shannon/Simpson/Chao1/ACE/Fisher/Pielou）、Beta 多样性（Bray-Curtis/Jaccard/JSD/weighted/unweighted UniFrac）、PCoA、差异丰度（Welch t/Wilcoxon/BH） | ✅ |
| **CellChat** | CellChat | 配体-受体互作对数据库、置换检验、互作评分、细胞类型聚合、显著性检验、FDR 校正 | ✅ |
| **CIBERSORT** | CIBERSORT 风格 | NNLS 非负最小二乘求解细胞类型分数、投影梯度下降、LM22 风格特征矩阵、Pearson 拟合优度+RMSE、分数归一化（Σ=1.0） | ✅ |
| **MCPcounter** | MCPcounter | 标记基因几何均值表达评分（Becht 2016）、10 种细胞种群默认标记集、对数域稳健计算 | ✅ |
| **ESTIMATE** | ESTIMATE | ssGSEA-style 富集打分（Yoshihara 2013）、50 基质基因+30 免疫基因特征集、ECDF 跨样本标准化、肿瘤纯度推断（cos 公式） | ✅ |
| **xCell** | xCell | 单样本 GSEA (ssGSEA) 富集打分、64 种细胞类型默认签名集（T/B/NK/髓系/树突/内皮/成纤维）、自定义签名集支持、按细胞类别汇总 | ✅ |
| **MutationalPatterns** | MutationalPatterns | 体细胞突变谱、96 通道矩阵、三核苷酸上下文、突变签名拟合 | ✅ |
| **maftools** | maftools | 癌症基因组学 MAF 解析、突变分类（SNV/Indel）、TMB 计算、突变谱、共现分析、oncoplot | ✅ |
| **LEfSe** | LEfSe | Kruskal-Wallis 检验、Wilcoxon 秩和检验、LDA 效应量（log10 类均值比）、双阈值过滤 | ✅ |

### 可视化

| 功能模块 | 对应 Bioconductor | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **pheatmap** | pheatmap | 增强型热图：层次聚类（complete/average/ward）、距离矩阵（euclidean/manhattan/correlation）、行/列注释、颜色方案、聚类间隙 | ✅ |
| **ComplexHeatmap** | ComplexHeatmap | 行/列聚类、颜色映射、热图注释、多个热图组合、分组拆分 | ✅ |
| **EnrichedHeatmap** | EnrichedHeatmap | 基因组信号归一化、目标区域窗口化、四种均值模式、行平滑、百分位裁剪、链方向处理 | ✅ |
| **EnhancedVolcano** | EnhancedVolcano | 差异表达基因分类、ASCII 渲染 | ✅ |
| **Gviz** | Gviz | 基因组可视化轨道系统（AnnotationTrack/GeneRegionTrack/DataTrack/IdeogramTrack/GenomeAxisTrack/SequenceTrack）、ASCII 渲染、特征查询与区域可视化 | ✅ |
| **GenomeDiagram** | GenomeDiagram | 数据模型 Track/Feature/Diagram、多种形状、SVG 生成、样式控制、自动标注、重叠检测 | ✅ |
| **karyoploteR** | karyoploteR | 染色体轨道、数据点、ASCII 渲染、核型可视化 | ✅ |
| **ggtree** | ggtree | 系统发育树可视化布局（矩形/放射状/无根）、节点坐标映射 | ✅ |
| **VennDiagram** | VennDiagram | 集合运算、Venn 区域计算、重叠统计、相似度（Jaccard/Dice/Overlap） | ✅ |
| **seqLogo** | seqLogo | PWM、信息含量 IC、ASCII logo、一致性序列 | ✅ |
| **graphics / reporting** | Bio.Graphics / ReportingTools | 序列 Logo、比对可视化、文本/表格/图形混合报告 | ✅ |

### 芯片、预处理与杂项

| 功能模块 | 对应 Bioconductor | 核心功能 | 状态 |
| :--- | :--- | :--- | :---: |
| **affy / gcrma / preprocessCore** | affy / gcrma / preprocessCore | RMA 标准化、背景校正（IdealMM/Express）、GC 校正、分位数归一化；quantile/invariant set/cyclic loess/contrast/percentile shift 归一化、log2 变换组合、归一化质量统计 | ✅ |
| **vsn** | vsn | 方差稳定化归一化 glog 变换、vsn2 | ✅ |
| **EDASeq** | EDASeq | RNA-seq 探索性分析、GC 含量归一化、基因长度校正 Loess、样本间归一化、RPKM | ✅ |
| **HTSFilter** | HTSFilter | RNA-seq count 过滤、CPM 归一化、按组最小样本数阈值、keep mask、文库大小 | ✅ |
| **RUVSeq / sva / ComBat** | RUVSeq / sva | RNA-seq 批次效应去除、log2 转换、RUVg/RUVs；替代变量分析、经验贝叶斯方法、PCA、ComBat/removeBatchEffect | ✅ |
| **NanoString** | NanoString | nCounter 数字条码：阳性/阴性对照归一化、管家基因归一化、content 归一化、LOQ 过滤、成像 QC、线性度、差异表达 | ✅ |
| **QFeatures** | QFeatures | 多 assay 层级容器（PSM/peptide/protein）、跨层级特征链接、聚合（sum/mean/median/max）、过滤/归一化/缺失值插补（KNN/mean/zero） | ✅ |
| **MSnbase / MsCoreUtils / MSstats** | MSnbase / MsCoreUtils / MSstats | Spectrum/Chromatogram/MSnSet、peak 查找、TIC 归一化、MA 平滑、SNIP 基线校正、SNR centroiding、refineCentroids、localMaxima、joinPeaks m/z 匹配、Savitzky-Golay 平滑、KNN 插补、medianPolish 聚合、肽段→蛋白汇总、样本 QC、质谱数据归一化与组间比较 | ✅ |
| **survival** | survival | Kaplan-Meier 估计器（Greenwood SE）、log-rank 检验、Cox 比例风险（Newton-Raphson + Breslow ties）、卡方 p 值、中位生存期 | ✅ |
| **VariantAnnotation** | VariantAnnotation | 变异类型检测、编码效应预测、变异汇总 | ✅ |
| **VariantFiltering** | VariantFiltering | 变异过滤与遗传模式：常染色体显性/隐性、X 连锁、复合杂合子 | ✅ |
| **Rsubread / featureCounts** | Rsubread/featureCounts | 链特异性（Stranded/Reversed/Unstranded）、重叠长度阈值、最小比对质量、多比对 read 处理（计数/分数计数）、配对末端 fragment 计数、CPM 归一化、多样本矩阵 | ✅ |
| **ballgown** | ballgown | 转录组水平差异、FPKM、t 检验、基因/转录本结构 | ✅ |
| **IsoformSwitchAnalyzeR** | IsoformSwitchAnalyzeR | 转录本异构体切换、PSI/DPSI/DIF、功能后果预测 | ✅ |
| **tximport** | tximport | Salmon quant.sf 解析、转录本→基因级别汇总、低表达过滤、ExpressionSet 转换 | ✅ |
| **impute** | impute | KNN/mean/median/LOCF/NOCB 缺失值插补 | ✅ |
| **Hmisc / rstatix** | Hmisc / rstatix | 相关性（Pearson/Spearman）、变量聚类、描述性统计、Somers' d、缺失值插补；T/Wilcoxon/ANOVA/Kruskal-Wallis/Friedman、BH-FDR/Bonferroni | ✅ |
| **PCAtools / factoextra** | PCAtools / factoextra | scree/biplot/outliers PCA 工具；特征值计算、方差解释率、个体/变量坐标、cos2 质量、贡献度评分、维度描述 | ✅ |
| **WGCNA** | WGCNA | 加权基因共表达网络、邻接矩阵、TOM 相似度、模块检测 | ✅ |
| **GENIE3** | GENIE3 | 回归树特征重要性、方差缩减、加权邻接矩阵、对称化网络 | ✅ |
| **minet** | minet | 互信息基因共表达网络：equalfreq/equalwidth 离散化、四种熵估计（empirical/Miller-Madow/Hausser-Strimmer shrinkage/Schürmann-Grassberger）、MIM 构建、ARACNE 数据处理不等式剪枝、CLR 行/列背景 z-score 校正、MRNET mRMR 前向选择、validate 精度/召回曲线与 AUROC/AUPR/maxF 评分 | ✅ |
| **genefilter** | genefilter | t/Wilcoxon 检验、方差过滤、CV 过滤、分位数过滤 | ✅ |
| **mixOmics** | mixOmics | 多组学整合：PLS 回归、稀疏 PLS (sPLS)、DIABLO 多块整合 | ✅ |
| **destiny / Rtsne / uwot** | destiny / Rtsne / uwot | 扩散映射（距离矩阵/高斯核/特征分解）；t-SNE（条件概率/梯度下降/Barnes-Hut）；UMAP（kNN/模糊单纯集/SGD/负采样） | ✅ |
| **GENESIS** | GENESIS | 亲属关系矩阵估计、PCA、遗传距离（欧氏/曼哈顿/IBS）、群体结构 | ✅ |
| **HilbertCurve** | HilbertCurve | Hilbert 曲线编码/解码、距离计算、基因组线性化、网格映射 | ✅ |
| **flowCore / openCyto / FlowSOM / diffcyt** | flowCore / openCyto / FlowSOM / diffcyt | FCS 处理、荧光补偿、门控（矩形/多边形/椭球/四象限/mindensity KDD/tailgate/flowClust t 混合 EM/rangeGate）、模板驱动门控流水线；高维流式：FlowSOM 聚类、负二项 GLM DA 检验、经验贝叶斯调节 t 检验 DS、BH-FDR | ✅ |
| **universalmotif** | universalmotif | Motif 结构、共识序列计算 | ✅ |
| **SystemPipeR** | SystemPipeR | 流水线编排、步骤管理、依赖关系、进度追踪 | ✅ |
| **pqsfinder** | pqsfinder | G-四链体（PQS）检测：G-run 穷举搜索、四分体加分/bulge/错配罚分评分系统、loop 长度惩罚、正负链搜索、非重叠/全重叠（overlapping）导出、deep 搜索逐位置 density 与 maxScores 向量、PQS 序列提取 | ✅ |
| **DNAcopy** | DNAcopy | 循环二分分割（CBS）拷贝数分割：置换检验显著性（alpha/nperm）、trimmed variance、max 统计量、多染色体独立分割、sdundo 相邻片段合并 | ✅ |
| **regioneR** | regioneR | 基因组区域置换检验 permTest：randomizeRegions 随机化（保宽度/染色体、避开 mask 掩蔽区）、circularRandomizeRegions 环状平移、numOverlaps/meanDistance 统计、z-score 与经验 p 值（greater/less/two.sided）、可复现种子 | ✅ |
| **MatrixEQTL** | MatrixEQTL | 快速 eQTL 分析、线性回归 SNP-基因关联、Benjamini-Hochberg FDR 校正、t 统计量、正则化不完全 Beta 函数 p 值 | ✅ |

---

## Others 其他功能与算法

### 序列组装算法

| 算法 | 对应工具 | 功能说明 | 状态 |
| :--- | :--- | :--- | :---: |
| **De Bruijn Graph** | SPAdes / Velvet | k-mer 节点、欧拉路径、序列组装、图简化 | ✅ |
| **Suffix Array & Suffix Tree** | libdivsufsort | 前缀倍增算法、LCP 数组、模式匹配、最长重复子串 | ✅ |
| **Overlap-Layout-Consensus** | Celera Assembler / ARACHNE | 重叠检测、哈密顿路径、一致性序列生成 | ✅ |
| **BWT + FM-index** | Bowtie2 / BWA | 后缀数组、BWT 变换、回溯搜索、精确模式匹配 | ✅ |
| **Bloom Filter** | Jellyfish / khmer | k-mer 计数、成员查询、误判率估算、近似去重 | ✅ |

### 机器学习 / 统计算法（独立实现）

| 模块 | 功能说明 | 状态 |
| :--- | :--- | :---: |
| **k-means++** | 距离计算、K-means++ 初始化、轮廓系数评估 | ✅ |
| **Hidden Markov Model** | 前向/后向算法、Viterbi 解码、Baum-Welch 训练、基因预测 | ✅ |
| **朴素贝叶斯** | 基于 k-mer 频率的序列分类、Laplace 平滑、top-K 预测 | ✅ |
| **马尔可夫链** | 1/2/3 阶训练、转移概率矩阵、CpG 岛 log-odds、加权采样生成序列 | ✅ |
| **逻辑回归** | Newton-Raphson 优化、sigmoid、对数似然收敛、操纵子预测示例 | ✅ |
| **最大熵分类器** | IIS 改进迭代尺度训练、指示特征函数、softmax 归一化 | ✅ |
| **神经网络** | 前馈 + 反向传播、sigmoid 激活、学习率/动量、XOR/iris 示例 | ✅ |
| **遗传算法** | 种群初始化、锦标赛/轮盘赌选择、单点/两点交叉、变异、精英保留、世代统计 | ✅ |

### 降维 / 可视化算法

| 模块 | 功能说明 | 状态 |
| :--- | :--- | :---: |
| **PCA / factoextra / PCAtools** | SVD 特征分解、方差解释率、个体/变量坐标、cos2、贡献度、scree/biplot/outliers | ✅ |
| **t-SNE (Rtsne)** | Barnes-Hut 近似、成对距离、条件概率、梯度下降、动量/early exaggeration | ✅ |
| **UMAP (uwot)** | k 近邻搜索、模糊单纯集、局部模糊集并集、SGD/负采样、min_dist/spread | ✅ |
| **扩散映射 (destiny)** | 距离矩阵计算、高斯核构建、特征分解、扩散分量 | ✅ |

---

## 架构设计

### 项目结构

```
IvanAXu/BioSeqs/
├── moon.mod                          # 模块配置 (name="IvanAXu/BioSeqs", version=0.1.9)
├── src/                              # 源代码（约 600 个 .mbt 模块）
│   ├── seq.mbt                       # Bio.Seq 序列对象
│   ├── seqio.mbt / fasta_io.mbt / ... # 序列 I/O（30+ 种格式）
│   ├── alignment.mbt / align_*.mbt   # 比对算法 + 20+ 种比对格式严格 API
│   ├── searchio.mbt / blast_*.mbt / ... # 搜索结果解析（BLAST/HMMER/Infernal/...）
│   ├── phylo.mbt / tree_*.mbt / paml_*.mbt # 发育树 + PAML 分子进化
│   ├── pdb*.mbt / mmcif*.mbt / binary_cif.mbt / ... # 结构分析与格式
│   ├── sam.mbt / bam.mbt / vcf.mbt / cram_wbtest.mbt # NGS 文件
│   ├── genomic_ranges.mbt / iranges.mbt / plyranges.mbt # 区间操作
│   ├── deseq2.mbt / edger.mbt / limma.mbt / seurat.mbt / ... # Bioconductor 分析套件
│   ├── de_bruijn.mbt / suffix_array_tree.mbt / olc.mbt / bwt_fm.mbt # 序列组装四大算法
│   ├── statistics.mbt / kmeans.mbt / hmm.mbt / neural_network.mbt / ... # ML 与统计
│   └── utils.mbt / data.mbt / ...    # 通用工具与常量
├── examples/                         # 示例程序（约 250 个演示 demo）
│   ├── basic_seq/                    # 基础序列操作
│   ├── seqcode_demo/ / seq_utils_demo/ # SeqUtils 序列工具（gc_fraction/GC123/nt_search/六框翻译/seq3/seq1）
│   ├── pqs_demo/                     # pqsfinder G-四链体检测（评分系统/双链/deep density 与 maxScores/overlapping/序列提取）
│   ├── melting_temp_demo/            # MeltingTemp 解链温度（Tm_GC 八套公式/Tm_NN 近邻热力学/盐与 Mg2+/错配/DMSO）
│   ├── dnacopy_demo/                 # DNAcopy CBS 拷贝数分割（sdundo 合并/多染色体）
│   ├── regione_r_demo/               # regioneR 区域置换检验（随机化/mask/permTest）
│   ├── minet_demo/                   # minet 互信息网络（离散化/MIM/ARACNE/CLR/MRNET/AUROC 验证）
│   ├── qvalue_demo/                  # qvalue FDR 全流程（π₀ smoother/bootstrap、q-value/pFDR、lfdr、empPvals、BH 对比）
│   ├── alignment_frequencies_demo/   # Bio.Align MSA 统计（替换计数矩阵/列频率/保守位点/序列权重/log-odds）
│   ├── bio_cluster_demo/             # Bio.Cluster 数值聚类（8 种距离/层次聚类/k-means/k-medoids/SOM/PCA/Record）
│   ├── hmm_markov_model_demo/        # 经典 Bio.HMM（构建器/Viterbi/前向-后向/Baum-Welch/KnownState 训练）
│   ├── togows_demo/                  # Bio.TogoWS（percent quoting/URL 构造/entry/search/convert/限速）
│   ├── apcluster_demo/               # apcluster 亲和传播（消息传递/exemplar/4 种相似度/apclusterK 定 K/噪声 runs）
│   ├── kegg_rest_demo/               # Bio.KEGG.REST（info/list/find/get/conv/link URL 构造/多条目/限速）
│   ├── scop_raf_dom_demo/            # Bio.SCOP.Raf/Dom（RAF 解析/索引/域切割/ATOM 提取/DOM 往返）
│   ├── entrez_eutils_demo/           # Bio.Entrez E-utilities（EPost/ESummary/ESpell/ELink/ECitMatch/请求日志）
│   ├── phylo_applications_demo/      # Bio.Phylo.Applications（FastTree/PhyML/RAxML 命令行构建）
│   ├── datastore_demo/               # Bio.Datastore（MD5/BagIt 清单/标签/fetch-through 缓存/校验）
│   ├── pdb_demo/ / phylo_demo/       # 结构与发育树
│   ├── deseq2_demo/ / edger_demo/ / limma_demo/ # 差异表达
│   ├── ebseq_demo/ / matrixeqtl_demo/ / dmrcate_demo/ / gostats_demo/ / lefse_demo/ # EBSeq/eQTL/DMR/GO/LEfSe
│   ├── seurat_demo/ / milo_demo/ / monocle3_demo/ # 单细胞
│   ├── de_bruijn_demo/ / olc_demo/   # 序列组装算法
│   └── ... （更多 examples/*_demo/）
├── test/
│   ├── moonbit/                      # MoonBit 单元测试（约 500 个测试文件，12551 用例）
│   │   ├── bio_seq_test.mbt / seqio_wb_test.mbt / ...
│   │   ├── alignment_test.mbt / pdb_test.mbt / phylo_test.mbt / ...
│   │   ├── deseq2_test.mbt / seurat_test.mbt / scran_test.mbt / ...
│   │   └── ...
│   └── python/                       # Python 参考实现与对比脚本
│       ├── python_reference.py / python_seqio_reference.py / ...
│       ├── compare.sh / compare_seqio.sh
│       └── python_bench.py
└── cmd/                              # 命令行工具
    ├── main/                         # Seq 基础测试工具
    ├── seqio_main/                   # SeqIO 测试工具
    ├── alignio_main/                 # AlignIO 测试工具
    └── bench/                        # 性能基准测试
```

### 样例测试

```bash
moon build                                              # ✅ 成功
moon test                                               # ✅ 12451 个测试全部通过
```

---

## 构建与测试

```bash
moon build      # 构建项目
moon test       # 运行全部测试 (12551 个测试用例)
```

---

## 核心模块速查

| MoonBit 模块文件 | 对应功能 |
| :--- | :--- |
| `seq.mbt` | `Bio.Seq` 序列对象与基础操作 |
| `seq_record.mbt` / `seqfeature.mbt` / `seqfeature_advanced.mbt` | 序列记录与特征（位置/CompoundLocation/修饰符） |
| `seqio.mbt` / `fasta_io.mbt` / `fastq_io.mbt` / `genbank_io.mbt` | 序列 I/O 统一接口（FASTA/FASTQ/GenBank 等 30+ 格式） |
| `sequtils.mbt` / `seq_utils.mbt` / `sequtils_advanced.mbt` / `seq_complexity.mbt` | GC 含量与 gc_fraction 歧义模式、GC123、nt_search、six_frame_translations、seq3/seq1、Tm/ORF/分子量/LCC 复杂度/Hamming 距离 |
| `melting_temp.mbt` / `melting_temp_advanced.mbt` | 解链温度：Tm_Wallace 经验法则、Tm_GC 八套经验公式、Tm_NN 近邻热力学（DNA_NN1-4/RNA_NN1-3/R_DNA_NN1 八套表 + 错配/末端/悬挂端表）、盐校正 1-7（Na+/K+/Tris/Mg2+/dNTP，Owczarzy 2008）、自互补、DMSO/甲酰胺化学校正、自定义热力学表 |
| `codon_usage.mbt` / `codon_align.mbt` / `codon_align_advanced.mbt` | 密码子使用分析与 dN/dS 选择压力检验 |
| `alignment.mbt` / `pairaligner.mbt` / `smith_waterman.mbt` / `needleman_wunsch.mbt` | 全局/局部比对算法 + PairwiseAligner |
| `align_*.mbt` (clustal/phylip/stockholm/msf/nexus/a2m/emboss/exonerate/maf/mauve/psl/sam/chain/bed/bigbed/...) | 各类比对格式严格读写 + 坐标映射 + 统计 + canonical 往返 |
| `alignment_map.mbt` / `alignment_counts.mbt` / `shared_reference_alignment.mbt` | Alignment map/mapall 坐标路径、counts gap/composition 评分、共享参考 PWA/MSA 合并 |
| `alignment_frequencies.mbt` | Bio.Align MSA 统计：每列字母计数/频率（`'-'` gap 行）、对称替换计数矩阵（序列权重、M=(M+Mᵀ)/2）、观测对频率与背景频率、q_ij/e_ij 相对频率/log-odds |
| `searchio.mbt` / `blast.mbt` / `blast_xml_advanced.mbt` | SearchIO 统一模型 + BLAST tabular/XML1/XML2 |
| `hmmer_io.mbt` / `infernal_io.mbt` / `hhr.mbt` / `exonerate_text.mbt` / `interproscan.mbt` | HMMER3/Infernal/HH-suite HHR/Exonerate C4/InterProScan 解析 |
| `phylo.mbt` / `tree_io.mbt` / `tree_construction.mbt` | 系统发育树解析 + UPGMA/NJ/WPGMA 建树 |
| `phylo_xml.mbt` / `phylo_nexml.mbt` / `phylo_cdao.mbt` / `parsimony.mbt` / `phylo_consensus.mbt` | PhyloXML/NeXML/CDAO、Fitch/Sankoff 简约性、共识树 |
| `pdb.mbt` / `pdb_io.mbt` / `mmcif.mbt` / `binary_cif.mbt` / `mmtf.mbt` | PDB/mmCIF/BinaryCIF/MMTF 结构解析与转换 |
| `cealign.mbt` / `qcp_superimposer.mbt` / `svd_superimposer.mbt` / `structure_alignment.mbt` / `ma_align.mbt` | CE/SVD/QCP 结构叠合 + 多结构比对 |
| `dssp.mbt` / `sasa.mbt` / `pdb_packing.mbt` / `internal_coords.mbt` | DSSP 二级结构、SASA (Shrake-Rupley/Lee-Richards)、包装密度、内部坐标 |
| `sam.mbt` / `bam.mbt` / `bgzf.mbt` / `vcf.mbt` / `variant_annotation.mbt` / `structural_variant.mbt` | NGS SAM/BAM/BGZF/VCF/SV 解析与注释 |
| `faidx.mbt` | pyfaidx 风格 FASTA 索引 |
| `genomic_ranges.mbt` / `granges_list.mbt` / `iranges.mbt` / `plyranges.mbt` | GenomicRanges/IRanges/plyranges 区间与 tidy 操作 |
| `pqs_finder.mbt` | pqsfinder G-四链体（PQS）检测：G-run 穷举搜索、四分体加分/bulge/错配/loop 罚分评分、正负链坐标映射、非重叠解析（revised non-overlapping storage）、`overlapping=true` 全构象导出（按 (start,end) 去重保留最高分）、`pqs_search` deep 搜索逐位置 `density`（构象重叠计数）与 `max_scores`（位置最高分）向量、`Pqs::sequence` 提取命中序列 |
| `dnacopy.mbt` | DNAcopy 循环二分分割（CBS）：置换检验 alpha/nperm、trimmed variance、多染色体独立分割、sdundo 合并 |
| `regione_r.mbt` | regioneR 区域置换检验：randomizeRegions/circularRandomizeRegions（含 mask）、numOverlaps/meanDistance、permTest z-score 与经验 p 值 |
| `qvalue.mbt` | Bioconductor qvalue FDR 分析：Storey π₀ 估计（Reinsch df=3 三次光滑样条 smoother、David Robinson 闭形式 MSE bootstrap、λ 路径 π₀(λ)）、q-value 与 pFDR 有限样本校正、fdr.level 显著集、局部 FDR lfdr（probit/logit 变换 + nrd0 高斯 KDE、cummax 单调化、截断）、经验 p 值 empPvals（池化合并排序 / 检验特异逐行矩阵，平局观测优先、下界 1/m₀） |
| `bio_cluster.mbt` | Bio.Cluster 数值聚类核心（cluster.c 逐行移植）：8 种距离（欧氏平方均值/City Block/Pearson/绝对/非中心/绝对非中心/Spearman/Kendall）、clusterdistance、treecluster 层次聚类（单 SLINK/完全/平均/质心连接 + Tree cut/sort）、treecluster_from_distance、kcluster（k-means/k-medians，多次 npass 取最优）、kmedoids、clustercentroids、somcluster 自组织映射（Kohonen 训练+分配）、cluster_pca（Golub-Reinsch SVD，均值中心化 + 奇异值排序）、distancematrix、Record Eisen Cluster/TreeView 格式解析（NAME/GWEIGHT/GORDER/EWEIGHT/EORDER）及全套聚类方法、可复现 cluster_seed |
| `hmm_markov_model.mbt` | 经典 Bio.HMM 移植（Biopython 1.84 退役前版本）：MarkovModelBuilder（状态/发射字母表、allow_all_transitions、转移/发射概率与伪计数、随机初始化 hmm_gen_random_array、set_initial_probabilities 残余概率均分）、HiddenMarkovModel（Viterbi 对数空间最优路径解码 + 回溯）、ScaledDp 缩放前向/后向算法（Durbin et al. p78，对数空间缩放防下溢）、BaumWelchTrainer（EM 迭代，change/迭代次数早停回调）、KnownStateTrainer（有标注路径最大似然估计）、TrainingSequence |
| `togows.mbt` | Bio.TogoWS REST 客户端：togows_entry/togows_search/togows_search_count/togows_search_iter/togows_convert 与对应 *_url 构造器、togows_quote percent quoting（urllib.parse.quote 语义，safe 默认 "/"）、数据库/字段/格式/转换白名单校验（错误消息字典序排序）、TogoWSClient 可注入 fetch 与 now 时钟、内置限速 TOGOWS_DELAY=1/3 秒（每秒 3 次查询） |
| `apcluster.mbt` | Bioconductor apcluster 亲和传播聚类：ap_neg_dist_mat（负 Minkowski^r，默认平方欧氏）/ap_exp_sim_mat（高斯 RBF，w 默认距离中位数）/ap_lin_sim_mat（线性，w 默认最大距离）/ap_cor_sim_mat（Pearson/Spearman，Spearman 为平局平均秩的 Pearson）；apcluster 消息传递（responsibility top-two 技巧 + availability 列更新、lam=0.9 阻尼、convits=100 exemplar 集合稳定判停）、ApResult（exemplars/clusters/labels/dpsim/expref/netsim）、标量与逐点偏好（默认有限相似度中位数）、apclusterK 偏好二分定 K + 最近 exemplar 剪枝合并、exact=false 容错、确定性噪声 runs（apcluster_seed） |
| `kegg_rest.mbt` | Bio.KEGG.REST REST 客户端：kegg_rest_info/list/list_entries/find/find_multi/get/get_multi/conv/conv_multi/link/link_multi 与 kegg_rest_url 构造器（info/list/find/get/conv/link 六种操作）、多条目 "+" 分隔、find 的 formula/exact_mass/mol_weight 选项、get 的 aaseq/ntseq/mol/kcf/image/kgml/json 输出选项、KeggRestClient 可注入 fetch 与 now 时钟、内置 1/3 秒限速（每秒 3 次查询） |
| `scop_raf_dom.mbt` | Bio.SCOP.Raf / Bio.SCOP.Dom：scop_raf_normalize_letters（3↔1 氨基酸映射）、scop_raf_seqmap_parse（38 列定宽头 + pdbid/chainid/version/flags + 7 列残基字段，0.01/0.02 版本校验）、scop_raf_parse/scop_raf_index（pdbid+chainid 键）、ScopSeqMap（length/get/residues/index/slice/add/extend）、scop_raf_get_seqmap（SCOP 链与残基定义切域）、scop_raf_get_atoms（PDB ATOM/HETATM 残基级提取）、scop_parse_residues；DomRecord 制表符记录解析（sid/pdbid/fragments/hierarchy）与 to_string 往返 |
| `entrez_eutils.mbt` | Bio.Entrez E-utilities：epost（WebEnv/QueryKey 会话）、esummary + parse_esummary（DocSum/Item 文档摘要）、espell + parse_espell（CorrectedQuery/Replaced）、elink + parse_elink（LinkSet/LinkSetDb/LinkName 跨库关联）、ecitmatch + parse_ecitmatch（journal|year|vol|page|author 引文 bdata 检索 PMID）、entrez_quote_plus URL 编码、EntrezClient 可注入 fetch 与 now 时钟、内置限速（无 API key 1/3 秒、有 key 1/11 秒） |
| `datastore.mbt` | Bio.Datastore MD5 键控 BagIt 风格 bag：md5_hex（RFC 1321 完整实现：K 表 floor(|sin(i+1)|·2³²)、四轮 64 步、大端输出）、Datastore::add/get/contains/remove/checksum/entry、manifest_md5_txt/manifest_lines/bagit_txt/bag_info_txt（Payload-Oxum 字节·文件数）、verify 逐条目校验、fetch-through 缓存（下载回调 + hit/miss 统计）、load_bag 外部清单加载（一致/过期清单判定）、字典序 names/tags |
| `phylo_applications.mbt` | Bio.Phylo.Applications 系统发育命令行包装（可变 builder 模式）：FastTreeCommandline（-nt/-n/-quote、-boot、-gtr/-gamma、-spr/-mlnni、-wag/-pseudo/-mlacc）、PhymlCommandline（-i/-d/-m/-b/-c/-a/-s、-o tlr、--rand_start/--n_rand_starts/--r_seed、--quiet）、RaxmlCommandline（-s/-n/-m、-f a、-x/-N、-q 分区、-t/-r/-g 起始/约束树、-e/-c、-k/-T/-o/-w、RAPID bootstrap 种子）与 build() 命令行构造 |
| `wgcna.mbt` / `genie3.mbt` / `minet.mbt` | 基因共表达/调控网络：WGCNA 加权共表达与 TOM 模块、GENIE3 回归树特征重要性；minet 互信息网络（equalfreq/equalwidth 离散化 + empirical/Miller-Madow/Hausser-Strimmer/Schürmann-Grassberger 四种熵估计，ARACNE DPI 剪枝、CLR 背景 z-score、MRNET mRMR 前向选择，validate precision/recall/AUROC/AUPR/maxF 评分） |
| `summarized_experiment.mbt` / `single_cell_experiment.mbt` / `spatial_experiment.mbt` / `multi_assay_experiment.mbt` / `tree_summarized_experiment.mbt` / `ragged_experiment.mbt` | Bioconductor 数据容器家族 |
| `deseq2.mbt` + `deseq2_advanced.mbt` / `edger.mbt` + `edger_advanced.mbt` / `limma.mbt` / `apeglm.mbt` | 差异表达三大套件 + apeglm LFC 收缩 |
| `ebseq.mbt` / `matrixeqtl.mbt` / `dmrcate.mbt` / `gostats.mbt` / `lefse.mbt` | EBSeq 经验贝叶斯 DE / MatrixEQTL eQTL 线性回归 / DMRcate 差异甲基化区域 / GOstats GO 富集 / LEfSe 微生物组生物标志物 |
| `seurat.mbt` / `scran.mbt` / `scuttle.mbt` / `scrapper.mbt` / `bluster.mbt` / `monocle3.mbt` / `slingshot.mbt` / `tradeSeq.mbt` / `velociraptor.mbt` / `scenic.mbt` / `infercnv.mbt` / `milo.mbt` / `muscat.mbt` / `zinbwave.mbt` / `celda.mbt` / `decontx.mbt` / `batchelor.mbt` / `sc_dbl_finder.mbt` / `droplet_utils_advanced.mbt` / `single_r_advanced.mbt` / `mast_advanced.mbt` | 单细胞 / 空间组学全栈分析套件 |
| `minfi.mbt` / `bsseq.mbt` / `methylkit.mbt` / `chipseeker.mbt` / `diffbind.mbt` / `peak_calling.mbt` / `bumphunter.mbt` | 甲基化与 ChIP-seq 分析 |
| `cluster_profiler.mbt` / `fgsea.mbt` / `gsva.mbt` / `gage.mbt` / `spia.mbt` / `enrichplot.mbt` | 富集分析统一框架与可视化 |
| `complex_heatmap.mbt` / `pheatmap.mbt` / `gviz.mbt` / `genome_diagram.mbt` / `enhanced_volcano.mbt` | 热图 / 基因组轨道 / 火山图可视化 |
| `de_bruijn.mbt` / `suffix_array_tree.mbt` / `olc.mbt` / `bwt_fm.mbt` / `bloom_filter.mbt` | 序列组装四大算法 + Bloom Filter |
| `motifs.mbt` / `motifs_advanced.mbt` / `jaspar.mbt` / `transfac.mbt` / `meme.mbt` / `motif_scan.mbt` / `seqlogo.mbt` | 序列模体识别 + JASPAR/TRANSFAC/MEME + FIMO 扫描 + Logo |
| `kmeans.mbt` / `hmm.mbt` / `logistic_regression.mbt` / `markov.mbt` / `neural_network.mbt` / `ga.mbt` | 独立实现 ML / 统计算法（k-means/HMM/逻辑回归/马尔可夫/神经网络/遗传算法） |
