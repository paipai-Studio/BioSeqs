///|
/// Bio.SeqUtils Codon Analysis — Demo
/// 密码子使用分析：GC123、六框翻译、CAI、Nc
///
/// 对应 Biopython Bio.SeqUtils 核心功能

fn main {
  let globin = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG"

  // 1. GC123 — 密码子三个位置的 GC 含量
  let gc = @src.seq_gc123(globin)
  println("=== GC123 (globin-like) ===")
  println("GC1(pos1): " + gc[0].to_string())
  println("GC2(pos2): " + gc[1].to_string())
  println("GC3(pos3): " + gc[2].to_string())
  println("(密码子第三位 GC3 通常最高，因为简并密码子差异多在 wobble 位)")
  println("")

  // 2. 六框翻译
  println("=== 六框翻译 (globin-like) ===")
  let frames = @src.seq_six_frame_translations(globin)
  let mut i = 0
  while i < 6 {
    let f = frames[i]
    println(f.frame_id + " [" + f.strand + "]: " + f.protein)
    i = i + 1
  }
  println("")

  // 3. 密码子表切换 (Vertebrate Mitochondrial)
  println("=== Vertebrate Mito 密码子表 (table=2) ===")
  println("Table 1: ATA=" + @src.codon_table(1)["ATA"] + ", AGA=" + @src.codon_table(1)["AGA"] + ", TGA=" + @src.codon_table(1)["TGA"])
  println("Table 2: ATA=" + @src.codon_table(2)["ATA"] + ", AGA=" + @src.codon_table(2)["AGA"] + ", TGA=" + @src.codon_table(2)["TGA"])
  println("(Mito: ATA=M, AGA/AGG=*, TGA=W — 与标准表不同)")
  println("")

  // 4. CAI — 密码子适应指数
  println("=== Codon Adaptation Index (CAI) ===")
  let ref_set : Array[String] = [
    "ATGGCCAAGCACGTGAGCGTGGAGAAGGTGGCCGAGTGCAAGAGCACCGGCACCCTGAACGGCACCGAGAAGGTGGAGGTGCTGGACAC",
    "ATGGCCTACATCGGCTGCGTGCTGAACAAGCACGAGAAGGACGACGTGATGGTCGTCAACGACCCGTGCGTGACCCACAAGCCCGAG",
    "ATGGCCGACAAGACCAAGGGCTGCGTGAAGCTGATCGTGAACGACGTGTTTGACCGGCTGATCGGCTTCGAGCACTGCTGCGCTGAC",
  ]
  let weights = @src.cai_weights(ref_set)
  let test_gene = "ATGGCCAAGATGATCGGCTACGTGCAGCGCGGCGTCGACGCCACCCTGGCCGAGATGCGCAAGGCCGGCATCACCGTGGCC"
  let cai = @src.codon_adaptation_index(test_gene, weights)
  println("测试基因 CAI (3-gene ref_set): " + cai.to_string())
  println("CAI ≈ 1.0 → 密码子使用高度优化")
  println("")

  // 5. Nc — 有效密码子数
  println("=== Effective Number of Codons (Nc/ENC) ===")
  let nc = @src.seq_enc_wright(globin)
  let nc_test = @src.seq_enc_wright(test_gene)
  println("globin Nc:   " + nc.to_string() + " (20-61)")
  println("test_gene Nc: " + nc_test.to_string() + " (20-61)")
  println("Nc ≈ 20 → 密码子极度偏倚, Nc ≈ 61 → 完全均匀")
}
