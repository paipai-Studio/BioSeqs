///|
/// 密码子优化演示: 宿主 CAI 权重 + GC/motif/均聚物约束.
fn main raise {
  println("=== 密码子优化 (宿主偏好 + 序列约束, 本地完成) ===")

  let protein = "MKKAVLAGGGEEEKKK*"
  // 模拟高表达基因参考集, 用于估计宿主密码子偏好.
  let reference = [
    "ATGAAAAAAGCGGTGCTGGCGGGTGGTGGTGAAGAAGAAAAAAAATAA", "ATGAAAAAAGCGGTGCTGGCGGGCGGCGGTGAAGAAGAAAAAAAATAA",
  ]
  let weights = @src.cai_weights(reference, table=1)
  let forbidden = @src.codon_optimizer_sites_for_enzymes(["EcoRI", "BamHI"])
  let config = @src.codon_optimizer_config(
    table=1,
    gc_min=0.45,
    gc_max=0.6,
    max_run=4,
    avoid_cpg=true,
    forbidden_motifs=forbidden,
    seed=2026UL,
  )

  println("\n蛋白序列: \{protein}")
  let dna = @src.codon_optimize(protein, weights, config~)
  println("优化后 DNA (长度 \{dna.length()}):")
  println("  \{dna}")

  // 回译验证.
  let table = @src.codon_table(1)
  let sb = StringBuilder()
  for i = 0; i + 3 <= dna.length(); i = i + 3 {
    sb.write_string(table.get(dna[i:i + 3].to_owned()).unwrap_or("X"))
  }
  println("\n回译验证: \{sb.to_string()}")
  println("与输入一致: \{sb.to_string() == protein}")

  let report = @src.analyze_codon_optimization(
    dna,
    weights,
    table=1,
    forbidden_motifs=forbidden,
  )
  println("\n质量报告:")
  println("  GC = \{report.gc}")
  println("  CAI = \{report.cai}")
  println("  最长均聚物 = \{report.max_run}")
  println("  禁止 motif 命中 = \{report.forbidden_count}")

  println("\n完成: 密码子设计全程本地完成, 数据不出本机.")
}
