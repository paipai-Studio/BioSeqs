# BioSeqs - MoonBit 生物信息学库

> https://github.com/paipai-Studio/BioSeqs
>
> https://gitlink.org.cn/IvanAXu/BioSeqs
> 
> https://mooncakes.io/docs/IvanAXu/BioSeqs
>

## 项目概述

BioSeqs 是一个基于 **MoonBit** 语言开发的生物信息学工具库，旨在复刻以下 Python 生物信息学库的核心功能：

| Python 库 | 功能范围 | 状态 |
|-----------|----------|------|
| **Biopython** | 序列处理、序列 I/O、比对、系统发育树、PDB 结构 | ✅ 已实现 |
| **scikit-bio** | 序列类型、比对算法 | ✅ 已实现 |
| **pysam** | SAM/BAM/VCF 文件解析 | ✅ 已实现 |
| **pyfaidx** | FASTA 快速索引访问 | ✅ 已实现 |
| **scikit-learn** | 机器学习特征提取功能 | ✅ 已实现 |

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
│   ├── phylo.mbt               # 系统发育树 (Clade/Tree)
│   ├── newick_io.mbt           # Newick 格式解析
│   ├── pdb.mbt                 # PDB 数据类型
│   ├── pdb_io.mbt              # PDB 文件 I/O
│   ├── sequtils.mbt            # 序列工具函数
│   ├── complement.mbt          # 互补碱基查找表
│   ├── codon_table.mbt         # 密码子翻译表
│   ├── sam.mbt                 # SAM 格式解析
│   ├── bam.mbt                 # BAM 格式解析
│   ├── bgzf.mbt               # BGZF 解压缩 (支持读取压缩的 BAM 文件)
│   ├── biostrings.mbt         # Biostrings 序列分析 (IUPAC、k-mer频率、RSCU、复杂度)
│   ├── vcf.mbt                 # VCF 格式解析
│   ├── faidx.mbt               # FASTA 快速索引访问 (pyfaidx)
│   ├── feature_extraction.mbt   # 机器学习特征提取
│   └── utils.mbt               # 通用工具函数
├── examples/                   # 示例程序
│   ├── basic_seq/              # 基础序列操作示例
│   ├── seqio_demo/             # 序列 I/O 示例
│   ├── alignment_demo/         # 序列比对示例
│   ├── phylo_demo/             # 系统发育树示例
│   ├── pdb_demo/               # PDB 结构解析示例
│   ├── sam_vcf_demo/           # SAM/VCF 解析示例
│   ├── faidx_demo/             # FASTA 索引示例
│   ├── ml_features/            # 机器学习特征提取示例
│   ├── bam_demo/               # BAM/BGZF 解析示例
│   ├── biostrings_demo/        # Biostrings 序列分析示例
│   └── deploy/                 # 部署配置
│       ├── ci.yml              # GitHub Actions CI/CD
│       └── Makefile            # 便捷操作命令
├── test/
│   ├── moonbit/                # MoonBit 测试文件
│   │   ├── bio_seq_test.mbt
│   │   ├── bio_seq_wb_test.mbt
│   │   ├── sequtils_test.mbt
│   │   ├── seqfeature_test.mbt
│   │   ├── seqio_wb_test.mbt
│   │   ├── phylo_test.mbt
│   │   ├── pdb_test.mbt
│   │   ├── alignment_test.mbt
│   │   ├── sam_test.mbt
│   │   ├── vcf_test.mbt
│   │   ├── faidx_test.mbt
│   │   ├── feature_extraction_test.mbt
│   │   ├── bgzf_test.mbt
│   │   ├── bam_test.mbt
│   │   └── biostrings_test.mbt
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
moon test --package IvanAXu/BioSeqs/test/moonbit        # ✅ 214 个测试全部通过
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
| `newick_io.mbt` | BioPython `Bio.Phylo.NewickIO` | Newick 格式 |
| `pdb.mbt` | BioPython `Bio.PDB` | PDB 数据类型 |
| `pdb_io.mbt` | BioPython `Bio.PDB.PDBIO` | PDB 文件 I/O |
| `sequtils.mbt` | BioPython `Bio.SeqUtils` | CRC32、GC 含量、Tm 计算、蛋白质分析 |
| `complement.mbt` | BioPython `Bio.Data.IUPACData` | 互补碱基表 |
| `codon_table.mbt` | BioPython `Bio.Data.CodonTable` | 密码子翻译表 |
| `sam.mbt` | pysam | SAM 文件解析 |
| `bam.mbt` | pysam | BAM 文件解析 |
| `bgzf.mbt` | pysam | BGZF 解压缩 |
| `biostrings.mbt` | Bioconductor Biostrings | IUPAC 支持、k-mer 频率、RSCU、复杂度、Tm 计算 |
| `vcf.mbt` | pysam | VCF 文件解析 |
| `faidx.mbt` | pyfaidx | FASTA 快速索引访问 |
| `feature_extraction.mbt` | 自定义 | 机器学习特征提取 (k-mer、组成、理化性质) |

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
Total tests: 214
Passed: 214
Failed: 0
```

### 测试模块分布

| 模块 | 测试文件 | 测试数 |
|------|----------|--------|
| 序列核心 | `bio_seq_test.mbt` | 27 |
| 序列工具 | `sequtils_test.mbt` | 14 |
| 序列特征 | `seqfeature_test.mbt` | 7 |
| SeqIO | `seqio_wb_test.mbt` | 18 |
| AlignIO | `bio_seq_wb_test.mbt` | 1 |
| 系统发育树 | `phylo_test.mbt` | 13 |
| PDB 结构 | `pdb_test.mbt` | 13 |
| 比对算法 | `alignment_test.mbt` | 14 |
| SAM 解析 | `sam_test.mbt` | 6 |
| VCF 解析 | `vcf_test.mbt` | 7 |
| BAM 解析 | `bam_test.mbt` | 11 |
| BGZF 解压缩 | `bgzf_test.mbt` | 2 |
| FASTA 索引 | `faidx_test.mbt` | 9 |
| 特征提取 | `feature_extraction_test.mbt` | 18 |
| Biostrings | `biostrings_test.mbt` | 21 |

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

项目提供 8 个示例程序，展示各模块的典型用法：

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

## 未来规划

- ✅ 实现 BAM 二进制格式解析
- ✅ 实现 BGZF 解压缩功能，支持读取压缩的 BAM 文件
- [ ] 实现 CRAM 格式支持
- [ ] 添加更多多样性指数计算
- ✅ 添加机器学习特征提取功能
- ✅ 优化 translate 性能（从 880 us → 601 us，提升 32%）
- [ ] 添加 SIMD 加速支持

## 许可证

Apache-2.0 License
