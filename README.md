# BioSeqs - 基于MoonBit的生物计算库

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

## 架构设计

### 项目结构

```
biolab/bio_seq/
├── moon.mod                    # 模块配置
├── moon.pkg                    # 包配置
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
│   ├── vcf.mbt                 # VCF 格式解析
│   ├── faidx.mbt               # FASTA 快速索引访问 (pyfaidx)
│   └── utils.mbt               # 通用工具函数
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
│   │   └── faidx_test.mbt
│   └── python/                 # Python 参考测试文件
│       ├── python_reference.py
│       ├── python_seqio_reference.py
│       ├── python_alignio_reference.py
│       ├── python_sequtils_reference.py
│       ├── python_seqfeature_reference.py
│       ├── python_phylo_reference.py
│       ├── python_pdb_reference.py
│       ├── python_skbio_alignment_reference.py
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
moon test --package IvanAXu/BioSeqs/test/moonbit        # ✅ 148 个测试全部通过
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
| `vcf.mbt` | pysam | VCF 文件解析 |
| `faidx.mbt` | pyfaidx | FASTA 快速索引访问 |

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

## 性能优化

### 优化策略

1. **O(1) 互补查找**: 使用 `FixedArray[UInt16]` 实现直接索引的互补碱基查找表，避免 Map 查找开销
2. **无分配字符串操作**: 使用 `unsafe_get` 和 `FixedArray` 避免中间字符串分配
3. **密码子快速查表**: 使用整数编码（c0*65536 + c1*256 + c2）作为键的密码子翻译表
4. **单遍扫描**: 反向互补操作使用单遍扫描，避免中间字符串分配
5. **预分配数组**: 字符串分割函数先统计分隔符数量，再一次性分配数组
6. **FASTA 快速索引**: 通过 .fai 索引实现 O(1) 随机访问，跳过逐行解析

### 性能基准测试 (Python 参考)

| 操作 | 序列长度 | 单次耗时 |
|------|----------|----------|
| complement | 100,000 | 106 us |
| reverse_complement | 100,000 | 178 us |
| transcribe | 100,000 | 244 us |
| translate | 100,000 | 17.5 ms |
| count | 100,000 | 98 us |
| find | 100,000 | 1.1 us |
| replace | 100,000 | 234 us |

## 测试验证

### 测试覆盖率

```
Total tests: 148
Passed: 148
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
| FASTA 索引 | `faidx_test.mbt` | 9 |

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
moon test --package biolab/bio_seq/test/moonbit

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

- [ ] 实现 BAM 二进制格式解析
- [ ] 实现 CRAM 格式支持
- [ ] 添加更多多样性指数计算
- [ ] 添加机器学习特征提取功能
- [ ] 优化 translate 性能（当前热点）
- [ ] 添加 SIMD 加速支持

## 许可证

MIT License
