///|
/// Bio.SeqUtils Codon Analysis — MoonBit 实现
/// 对应 Biopython Bio.SeqUtils + Bio.Data.CodonTable 核心算法
///
/// 密码子表覆盖: 27 个 NCBI 标准表 (1-6, 9-16, 21-33)
/// 分析函数: GC123, 六框翻译, CAI, Nc/ENC, stop/start codons
///
/// 注意: 与 sequtils_advanced.mbt / codon_align.mbt 重名的函数已加 seq_ 前缀区分

// ============ NCBI 密码子表数据 ============


/// NCBI Table 1: Standard
let _codon_table_1 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "*"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 2: Vertebrate Mitochondrial
let _codon_table_2 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "*"), ("AGC", "S"), ("AGG", "*"), ("AGT", "S"), ("ATA", "M"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "W"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 3: Yeast Mitochondrial
let _codon_table_3 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "M"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "T"), ("CTC", "T"), ("CTG", "T"), ("CTT", "T"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "W"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 4: Mold Mitochondrial
let _codon_table_4 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "W"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 5: Invertebrate Mitochondrial
let _codon_table_5 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "S"), ("AGC", "S"), ("AGG", "S"), ("AGT", "S"), ("ATA", "M"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "W"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 6: Ciliate Nuclear
let _codon_table_6 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "Q"), ("TAC", "Y"), ("TAG", "Q"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "*"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 9: Echinoderm Mitochondrial
let _codon_table_9 : Map[String, String] = Map(
  [("AAA", "N"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "S"), ("AGC", "S"), ("AGG", "S"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "W"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 10: Euplotid Nuclear
let _codon_table_10 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "C"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 11: Bacterial
let _codon_table_11 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "*"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 12: Alternative Yeast Nuclear
let _codon_table_12 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "S"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "*"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 13: Ascidian Mitochondrial
let _codon_table_13 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "G"), ("AGC", "S"), ("AGG", "G"), ("AGT", "S"), ("ATA", "M"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "W"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 14: Alternative Flatworm Mitochondrial
let _codon_table_14 : Map[String, String] = Map(
  [("AAA", "N"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "S"), ("AGC", "S"), ("AGG", "S"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "Y"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "W"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 15: Blepharisma Macronuclear
let _codon_table_15 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "Q"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "*"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 16: Chlorophycean Mitochondrial
let _codon_table_16 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "L"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "*"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 21: Trematode Mitochondrial
let _codon_table_21 : Map[String, String] = Map(
  [("AAA", "N"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "S"), ("AGC", "S"), ("AGG", "S"), ("AGT", "S"), ("ATA", "M"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "W"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 22: Scenedesmus obliquus Mitochondrial
let _codon_table_22 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "L"), ("TAT", "Y"), ("TCA", "*"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "*"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 23: Thraustochytrium Mitochondrial
let _codon_table_23 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "*"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "*"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 24: Pterobranchia Mitochondrial
let _codon_table_24 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "S"), ("AGC", "S"), ("AGG", "K"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "W"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 25: Candidate Division SR1
let _codon_table_25 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "G"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 26: Pachysolen tannophilus Nuclear
let _codon_table_26 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "A"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "*"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 27: Karyorelict Nuclear
let _codon_table_27 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "Q"), ("TAC", "Y"), ("TAG", "Q"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "W"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 28: Condylostoma Nuclear
let _codon_table_28 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "Q"), ("TAC", "Y"), ("TAG", "Q"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "W"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 29: Mesodinium Nuclear
let _codon_table_29 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "Y"), ("TAC", "Y"), ("TAG", "Y"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "*"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 30: Peritrich Nuclear
let _codon_table_30 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "E"), ("TAC", "Y"), ("TAG", "E"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "*"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 31: Blastocrithidia Nuclear
let _codon_table_31 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "E"), ("TAC", "Y"), ("TAG", "E"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "W"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 32: Balanophoraceae Plastid
let _codon_table_32 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "R"), ("AGC", "S"), ("AGG", "R"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "*"), ("TAC", "Y"), ("TAG", "W"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "*"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

/// NCBI Table 33: Cephalodiscidae Mitochondrial
let _codon_table_33 : Map[String, String] = Map(
  [("AAA", "K"), ("AAC", "N"), ("AAG", "K"), ("AAT", "N"), ("ACA", "T"), ("ACC", "T"), ("ACG", "T"), ("ACT", "T"), ("AGA", "S"), ("AGC", "S"), ("AGG", "K"), ("AGT", "S"), ("ATA", "I"), ("ATC", "I"), ("ATG", "M"), ("ATT", "I"), ("CAA", "Q"), ("CAC", "H"), ("CAG", "Q"), ("CAT", "H"), ("CCA", "P"), ("CCC", "P"), ("CCG", "P"), ("CCT", "P"), ("CGA", "R"), ("CGC", "R"), ("CGG", "R"), ("CGT", "R"), ("CTA", "L"), ("CTC", "L"), ("CTG", "L"), ("CTT", "L"), ("GAA", "E"), ("GAC", "D"), ("GAG", "E"), ("GAT", "D"), ("GCA", "A"), ("GCC", "A"), ("GCG", "A"), ("GCT", "A"), ("GGA", "G"), ("GGC", "G"), ("GGG", "G"), ("GGT", "G"), ("GTA", "V"), ("GTC", "V"), ("GTG", "V"), ("GTT", "V"), ("TAA", "Y"), ("TAC", "Y"), ("TAG", "*"), ("TAT", "Y"), ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"), ("TGA", "W"), ("TGC", "C"), ("TGG", "W"), ("TGT", "C"), ("TTA", "L"), ("TTC", "F"), ("TTG", "L"), ("TTT", "F")], capacity=64)

// ============ NCBI 密码子表元数据 ============

/// 密码子表名称 (按 NCBI ID 索引)
let _codon_table_names : Map[Int, String] = Map([
  (1, "Standard"),
  (2, "Vertebrate Mitochondrial"),
  (3, "Yeast Mitochondrial"),
  (4, "Mold Mitochondrial"),
  (5, "Invertebrate Mitochondrial"),
  (6, "Ciliate Nuclear"),
  (9, "Echinoderm Mitochondrial"),
  (10, "Euplotid Nuclear"),
  (11, "Bacterial"),
  (12, "Alternative Yeast Nuclear"),
  (13, "Ascidian Mitochondrial"),
  (14, "Alternative Flatworm Mitochondrial"),
  (15, "Blepharisma Macronuclear"),
  (16, "Chlorophycean Mitochondrial"),
  (21, "Trematode Mitochondrial"),
  (22, "Scenedesmus obliquus Mitochondrial"),
  (23, "Thraustochytrium Mitochondrial"),
  (24, "Pterobranchia Mitochondrial"),
  (25, "Candidate Division SR1"),
  (26, "Pachysolen tannophilus Nuclear"),
  (27, "Karyorelict Nuclear"),
  (28, "Condylostoma Nuclear"),
  (29, "Mesodinium Nuclear"),
  (30, "Peritrich Nuclear"),
  (31, "Blastocrithidia Nuclear"),
  (32, "Balanophoraceae Plastid"),
  (33, "Cephalodiscidae Mitochondrial"),
], capacity=33)

/// 终止密码子 (按 NCBI ID 索引)
let _codon_table_stops : Map[Int, Array[String]] = Map([
  (1, ["TAA", "TAG", "TGA"]),
  (2, ["TAA", "TAG", "AGA", "AGG"]),
  (3, ["TAA", "TAG"]),
  (4, ["TAA", "TAG"]),
  (5, ["TAA", "TAG"]),
  (6, ["TGA"]),
  (9, ["TAA", "TAG"]),
  (10, ["TAA", "TAG"]),
  (11, ["TAA", "TAG", "TGA"]),
  (12, ["TAA", "TAG", "TGA"]),
  (13, ["TAA", "TAG"]),
  (14, ["TAG"]),
  (15, ["TAA", "TGA"]),
  (16, ["TAA", "TGA"]),
  (21, ["TAA", "TAG"]),
  (22, ["TCA", "TAA", "TGA"]),
  (23, ["TTA", "TAA", "TAG", "TGA"]),
  (24, ["TAA", "TAG"]),
  (25, ["TAA", "TAG"]),
  (26, ["TAA", "TAG", "TGA"]),
  (27, ["TGA"]),
  (28, ["TAA", "TAG", "TGA"]),
  (29, ["TGA"]),
  (30, ["TGA"]),
  (31, ["TAA", "TAG"]),
  (32, ["TAA", "TGA"]),
  (33, ["TAG"]),
], capacity=33)

/// 起始密码子 (按 NCBI ID 索引)
let _codon_table_starts : Map[Int, Array[String]] = Map([
  (1, ["TTG", "CTG", "ATG"]),
  (2, ["ATT", "ATC", "ATA", "ATG", "GTG"]),
  (3, ["ATA", "ATG", "GTG"]),
  (4, ["TTA", "TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"]),
  (5, ["TTG", "ATT", "ATC", "ATA", "ATG", "GTG"]),
  (6, ["ATG"]),
  (9, ["ATG", "GTG"]),
  (10, ["ATG"]),
  (11, ["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"]),
  (12, ["CTG", "ATG"]),
  (13, ["TTG", "ATA", "ATG", "GTG"]),
  (14, ["ATG"]),
  (15, ["ATG"]),
  (16, ["ATG"]),
  (21, ["ATG", "GTG"]),
  (22, ["ATG"]),
  (23, ["ATT", "ATG", "GTG"]),
  (24, ["TTG", "CTG", "ATG", "GTG"]),
  (25, ["TTG", "ATG", "GTG"]),
  (26, ["CTG", "ATG"]),
  (27, ["ATG"]),
  (28, ["ATG"]),
  (29, ["ATG"]),
  (30, ["ATG"]),
  (31, ["ATG"]),
  (32, ["ATA", "ATC", "ATG", "ATT", "CTG", "GTG", "TTG"]),
  (33, ["ATG", "CTG", "GTG", "TTG"]),
], capacity=33)

/// 有效 NCBI 密码子表 ID 列表
let _valid_codon_table_ids : Array[Int] = [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33]

///|
/// 获取 NCBI 密码子表 (Biopython Bio.Data.CodonTable.unambiguous_dna_by_id)
/// 支持 27 个 NCBI 标准表 (1-6, 9-16, 21-33)，无效 ID 默认返回 Table 1
pub fn codon_table(table_id : Int) -> Map[String, String] {
  match table_id {

    1 => _codon_table_1
    2 => _codon_table_2
    3 => _codon_table_3
    4 => _codon_table_4
    5 => _codon_table_5
    6 => _codon_table_6
    9 => _codon_table_9
    10 => _codon_table_10
    11 => _codon_table_11
    12 => _codon_table_12
    13 => _codon_table_13
    14 => _codon_table_14
    15 => _codon_table_15
    16 => _codon_table_16
    21 => _codon_table_21
    22 => _codon_table_22
    23 => _codon_table_23
    24 => _codon_table_24
    25 => _codon_table_25
    26 => _codon_table_26
    27 => _codon_table_27
    28 => _codon_table_28
    29 => _codon_table_29
    30 => _codon_table_30
    31 => _codon_table_31
    32 => _codon_table_32
    33 => _codon_table_33
    _ => _codon_table_1
  }
}

///|
/// 六框翻译结果
pub struct SixFrameTranslation {
  frame_id : String      // "Fwd-1", "Fwd-2", "Fwd-3", "Rev-1", ...
  strand : String        // "+" 或 "-"
  frame_offset : Int      // 0, 1, 2
  protein : String       // 翻译后的氨基酸序列
}

// ============ 内部辅助函数 ============

///|
/// 单框翻译（内部）
fn _cu_translate_frame(seq : String, frame : Int, codon_map : Map[String, String]) -> String {
  let n = seq.length()
  let buf = StringBuilder()
  let mut i = frame
  while i + 3 <= n {
    let codon = seq[i:i + 3].to_owned().to_upper()
    let aa = match codon_map.get(codon) {
      Some(a) => a
      None => "X"
    }
    buf.write_string(aa)
    i = i + 3
  }
  buf.to_string()
}

///|
/// 密码子频率统计（内部）
fn _cu_codon_counts(sequences : Array[String], codon_map : Map[String, String]) -> Map[String, Int] {
  let counts : Map[String, Int] = Map([], capacity=64)
  let codon_keys = codon_map.keys().to_array()
  let n = codon_keys.length()
  let mut c = 0
  while c < n {
    counts[codon_keys[c]] = 0
    c = c + 1
  }

  let n_seqs = sequences.length()
  let mut s = 0
  while s < n_seqs {
    let seq = sequences[s].to_upper()
    let mut i = 0
    while i + 3 <= seq.length() {
      let codon = seq[i:i + 3].to_owned().to_upper()
      if counts.contains(codon) {
        counts[codon] = counts[codon] + 1
      }
      i = i + 3
    }
    s = s + 1
  }
  counts
}

///|
/// 按氨基酸分组密码子（内部）
fn _cu_aa_groups(codon_map : Map[String, String]) -> Map[String, Array[String]] {
  let groups : Map[String, Array[String]] = Map([], capacity=20)
  let codon_keys = codon_map.keys().to_array()
  let n = codon_keys.length()
  let mut k = 0
  while k < n {
    let codon = codon_keys[k]
    match codon_map.get(codon) {
      Some(aa) => {
        if aa != "*" {
          if !groups.contains(aa) {
            groups[aa] = Array::new()
          }
          groups[aa].push(codon)
        }
      }
      None => ()
    }
    k = k + 1
  }
  groups
}


// ============ 公开 API ============

///|
/// **GC123** — 密码子三个位置的 GC 含量
/// 对应 Biopython Bio.SeqUtils.GC123
/// 返回 [GC_pos1, GC_pos2, GC_pos3]（比例值 0.0-1.0）
///
/// 注意: 本函数返回 Array[Double]，与 sequtils_advanced.mbt 的 gc123 返回 GC123Result 不同
///
/// 例: seq_gc123("ATGGCCATTGTA") => [0.5, 0.5, 0.333333]
pub fn seq_gc123(seq : String) -> Array[Double] {
  let seq_upper = seq.to_upper()
  let gc_counts = [0, 0, 0]
  let total_counts = [0, 0, 0]

  let mut i = 0
  while i + 3 <= seq_upper.length() {
    let pos0 = seq_upper[i:i + 1].to_owned()
    let pos1 = seq_upper[i + 1:i + 2].to_owned()
    let pos2 = seq_upper[i + 2:i + 3].to_owned()
    let bases = [pos0, pos1, pos2]

    let mut j = 0
    while j < 3 {
      let b = bases[j]
      if b == "G" || b == "C" {
        gc_counts[j] = gc_counts[j] + 1
      }
      if b == "A" || b == "T" || b == "G" || b == "C" {
        total_counts[j] = total_counts[j] + 1
      }
      j = j + 1
    }
    i = i + 3
  }

  let result : Array[Double] = Array::new()
  let mut j = 0
  while j < 3 {
    if total_counts[j] == 0 {
      result.push(0.0)
    } else {
      result.push(gc_counts[j].to_double() / total_counts[j].to_double())
    }
    j = j + 1
  }
  result
}

///|
/// **六框翻译** — 对 DNA 序列进行所有 6 种阅读框的翻译
/// 对应 Biopython Bio.SeqUtils.six_frame_translations
///
/// 支持全部 25 个 NCBI 密码子表
///
/// 参数:
///   - seq: DNA 序列
///   - table: NCBI 密码子表编号 (默认 1=Standard)
///
/// 例:
///   seq_six_frame_translations("ATGGCCATTGTAATGGGCCGCTGAA", table=1)
///   返回 Fwd-1/2/3 + Rev-1/2/3 共 6 个框架
pub fn seq_six_frame_translations(seq : String, table? : Int = 1) -> Array[SixFrameTranslation] {
  let seq_upper = seq.to_upper()
  let codon_map = codon_table(table)
  let result : Array[SixFrameTranslation] = Array::new()

  let mut frame = 0
  while frame < 3 {
    let protein = _cu_translate_frame(seq_upper, frame, codon_map)
    result.push(SixFrameTranslation::{
      frame_id: "Fwd-" + (frame + 1).to_string(),
      strand: "+",
      frame_offset: frame,
      protein: protein,
    })
    frame = frame + 1
  }

  let rc_seq = reverse_complement(seq_upper)
  let mut frame2 = 0
  while frame2 < 3 {
    let protein = _cu_translate_frame(rc_seq, frame2, codon_map)
    result.push(SixFrameTranslation::{
      frame_id: "Rev-" + (frame2 + 1).to_string(),
      strand: "-",
      frame_offset: frame2,
      protein: protein,
    })
    frame2 = frame2 + 1
  }

  result
}

///|
/// **Codon Adaptation Index (CAI) 权重** — 计算密码子使用偏倚度权重
/// 对应 Biopython Bio.SeqUtils.CodonAdaptationIndex
///
/// 算法 (Sharp & Li 1987):
///   w(codon) = freq(codon) / max_freq(同氨基酸密码子组)
///
/// 参数:
///   - sequences: 参考序列集合（通常是高表达基因）
///   - table: NCBI 密码子表编号
///
/// 返回: Map[String, Double] — 每个密码子的 CAI 权重 (0.0-1.0)
///   高 w 值表示该密码子在参考集中是"偏好密码子"
///
/// 例:
///   let w = cai_weights(["ATGGCCAAGG...", "ATGGCCTACAT..."])
///   // w["ATG"] = 1.0 (Met 只有一个密码子)
///   // w["CTG"] = 1.0, w["CTT"] = 0.45 (若 CTG 是最常用 Leu 密码子)
pub fn cai_weights(sequences : Array[String], table? : Int = 1) -> Map[String, Double] {
  let codon_map = codon_table(table)
  let counts = _cu_codon_counts(sequences, codon_map)
  let codon_keys = codon_map.keys().to_array()

  let mut total = 0
  let n_codons = codon_keys.length()
  let mut k = 0
  while k < n_codons {
    match counts.get(codon_keys[k]) {
      Some(c) => total = total + c
      None => ()
    }
    k = k + 1
  }

  let freqs : Map[String, Double] = Map([], capacity=64)
  let mut k2 = 0
  while k2 < n_codons {
    let codon = codon_keys[k2]
    if total > 0 {
      match counts.get(codon) {
        Some(c) => freqs[codon] = c.to_double() / total.to_double()
        None => freqs[codon] = 0.0
      }
    } else {
      freqs[codon] = 0.0
    }
    k2 = k2 + 1
  }

  let aa_groups = _cu_aa_groups(codon_map)
  let aa_keys = aa_groups.keys().to_array()
  let n_aa = aa_keys.length()

  let weights : Map[String, Double] = Map([], capacity=64)
  let mut a = 0
  while a < n_aa {
    let aa = aa_keys[a]
    let aa_codons = match aa_groups.get(aa) {
      Some(c) => c
      None => Array::new()
    }

    let mut max_freq = 0.0
    let n_aa_codons = aa_codons.length()
    let mut ac = 0
    while ac < n_aa_codons {
      let codon = aa_codons[ac]
      match freqs.get(codon) {
        Some(f) => if f > max_freq { max_freq = f }
        None => ()
      }
      ac = ac + 1
    }

    let mut ac2 = 0
    while ac2 < n_aa_codons {
      let codon = aa_codons[ac2]
      if max_freq > 0.0 {
        match freqs.get(codon) {
          Some(f) => weights[codon] = f / max_freq
          None => weights[codon] = 1.0
        }
      } else {
        weights[codon] = 1.0
      }
      ac2 = ac2 + 1
    }
    a = a + 1
  }

  weights
}

///|
/// **Codon Adaptation Index (CAI)** — 计算单条序列的密码子适应指数
/// 对应 Biopython Bio.SeqUtils.CodonAdaptationIndex
///
/// 算法 (Sharp & Li 1987):
///   CAI = exp(mean(ln(w_i)))
///
/// 参数:
///   - seq: 编码 DNA 序列 (必须是 3 的倍数)
///   - weights: 由 cai_weights() 计算的密码子权重
///   - table: NCBI 密码子表编号
///
/// 返回: CAI 值 (0.0 - 1.0)
///   - 接近 1.0: 密码子使用高度优化
///   - 低值: 密码子使用与参考集差异大
pub fn codon_adaptation_index(
  seq : String,
  weights : Map[String, Double],
  table? : Int = 1,
) -> Double {
  let seq_upper = seq.to_upper()
  let codon_map = codon_table(table)

  let ln_w_values : Array[Double] = Array::new()
  let mut i = 0
  while i + 3 <= seq_upper.length() {
    let codon = seq_upper[i:i + 3].to_owned().to_upper()
    let is_valid = match codon_map.get(codon) {
      Some(aa) => aa != "*"
      None => false
    }
    if is_valid {
      let w_opt = match weights.get(codon) {
        Some(w) => w
        None => 0.0
      }
      if w_opt > 0.0 {
        ln_w_values.push(@math.ln(w_opt))
      }
    }
    i = i + 3
  }

  if ln_w_values.length() == 0 {
    return 0.0
  }

  let mut sum = 0.0
  let n = ln_w_values.length()
  let mut j = 0
  while j < n {
    sum = sum + ln_w_values[j]
    j = j + 1
  }

  @math.exp(sum / n.to_double())
}

///|
/// **Effective Number of Codons (Nc/ENC)** — Wright 1990 精确算法
/// 对应 Biopython Bio.SeqUtils
///
/// 算法 (Wright 1990):
///   F-hat = (n * sum(p_i^2) - 1) / (n - 1)
///   Nc_aa = 1 / F-hat
///   最终 Nc = sum(Nc_aa for 每个氨基酸)
///   Met 和 Trp 各只有 1 密码子，跳过不计入
///
/// 参数:
///   - seq: 编码 DNA 序列 (3 的倍数)
///   - table: NCBI 密码子表编号
///
/// 返回: Nc 值 (约 20.0 - 61.0)
///   - ~20: 密码子使用极度偏倚
///   - ~61: 密码子使用完全均匀
pub fn seq_enc_wright(seq : String, table? : Int = 1) -> Double {
  let seq_upper = seq.to_upper()
  let codon_map = codon_table(table)

  let codon_counts = _cu_codon_counts([seq_upper], codon_map)
  let aa_groups = _cu_aa_groups(codon_map)
  let aa_keys = aa_groups.keys().to_array()
  let n_aa = aa_keys.length()

  let mut total_nc = 0.0

  let mut a = 0
  while a < n_aa {
    let aa = aa_keys[a]
    let aa_codons = match aa_groups.get(aa) {
      Some(c) => c
      None => Array::new()
    }

    if aa_codons.length() <= 1 {
      a = a + 1
      continue
    }

    let mut total_c = 0
    let n_aa_c = aa_codons.length()
    let mut ac = 0
    while ac < n_aa_c {
      match codon_counts.get(aa_codons[ac]) {
        Some(c) => total_c = total_c + c
        None => ()
      }
      ac = ac + 1
    }

    if total_c <= 1 {
      total_nc = total_nc + aa_codons.length().to_double()
    } else {
      let mut sum_p2 = 0.0
      let mut ac2 = 0
      while ac2 < n_aa_c {
        let codon = aa_codons[ac2]
        let count = match codon_counts.get(codon) {
          Some(c) => c
          None => 0
        }
        let p = count.to_double() / total_c.to_double()
        sum_p2 = sum_p2 + p * p
        ac2 = ac2 + 1
      }

      let f_hat = (total_c.to_double() * sum_p2 - 1.0) / (total_c.to_double() - 1.0)
      if f_hat > 0.0 {
        total_nc = total_nc + 1.0 / f_hat
      } else {
        total_nc = total_nc + aa_codons.length().to_double()
      }
    }

    a = a + 1
  }

  total_nc
}


// ============ 新增便利 API (Biopython 对齐) ============

///|
/// 查询密码子表的官方名称
/// 对应 Biopython: CodonTable.unambiguous_dna_by_id[table_id].names[0]
///
/// 例: codon_table_name(2) => "Vertebrate Mitochondrial"
pub fn codon_table_name(table_id : Int) -> String {
  match _codon_table_names.get(table_id) {
    Some(n) => n
    None => "Unknown"
  }
}

///|
/// 查询密码子表的终止密码子
/// 对应 Biopython: CodonTable.unambiguous_dna_by_id[table_id].stop_codons
///
/// 例: codon_table_stop_codons(1) => ["TAA", "TAG", "TGA"]
pub fn codon_table_stop_codons(table_id : Int) -> Array[String] {
  match _codon_table_stops.get(table_id) {
    Some(s) => s.copy()
    None => ["TAA", "TAG", "TGA"]
  }
}

///|
/// 查询密码子表的起始密码子
/// 对应 Biopython: CodonTable.unambiguous_dna_by_id[table_id].start_codons
///
/// 例: codon_table_start_codons(1) => ["TTG", "CTG", "ATG"]
pub fn codon_table_start_codons(table_id : Int) -> Array[String] {
  match _codon_table_starts.get(table_id) {
    Some(s) => s.copy()
    None => ["ATG"]
  }
}

///|
/// 检查密码子表 ID 是否为有效 NCBI 表
///
/// 例: codon_table_is_valid(1) => true, codon_table_is_valid(7) => false
pub fn codon_table_is_valid(table_id : Int) -> Bool {
  _valid_codon_table_ids.contains(table_id)
}

///|
/// 返回所有可用的 NCBI 密码子表 ID
/// 对应 Biopython: list(CodonTable.unambiguous_dna_by_id.keys())
pub fn codon_table_available_ids() -> Array[Int] {
  _valid_codon_table_ids.copy()
}

///|
/// 获取反向密码子表 (AA -> Array[codon])
/// 对应 Biopython: table.back_table (近似)
///
/// 参数:
///   - table_id: NCBI 密码子表编号
///
/// 返回: Map[String, Array[String]] — 每个氨基酸对应的密码子列表
///   终止密码子 (*) 不包含在结果中
///
/// 例: codon_table_reverse(1)["M"] => ["ATG"]
pub fn codon_table_reverse(table_id? : Int = 1) -> Map[String, Array[String]] {
  let codon_map = codon_table(table_id)
  let reverse : Map[String, Array[String]] = Map([], capacity=20)
  let codon_keys = codon_map.keys().to_array()
  let n = codon_keys.length()
  let mut k = 0
  while k < n {
    let codon = codon_keys[k]
    match codon_map.get(codon) {
      Some(aa) => {
        if aa != "*" {
          if !reverse.contains(aa) {
            reverse[aa] = Array::new()
          }
          reverse[aa].push(codon)
        }
      }
      None => ()
    }
    k = k + 1
  }
  reverse
}
