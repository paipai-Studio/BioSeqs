// ===========================================================================
// codon_optimizer_test.mbt — 密码子优化 / 约束修复测试
// ===========================================================================

///|
fn codon_opt_test_weights(
  codons : Array[String],
  w : Double,
) -> Map[String, Double] {
  let m : Map[String, Double] = Map([], capacity=codons.length())
  for c in codons {
    m[c] = w
  }
  m
}

///|
fn codon_opt_test_translate(dna : String) -> String {
  let table = @src.codon_table(1)
  let sb = StringBuilder()
  for i = 0; i + 3 <= dna.length(); i = i + 3 {
    let codon = dna[i:i + 3].to_owned().to_upper()
    sb.write_string(table.get(codon).unwrap_or("X"))
  }
  sb.to_string()
}

///|
test "codon optimizer config defaults and validation" {
  let c = @src.codon_optimizer_config()
  assert_eq(c.table, 1)
  assert_eq(c.gc_min, 0.0)
  assert_eq(c.gc_max, 1.0)

  let mut raised = false
  try {
    let bad = @src.codon_optimizer_config(table=99)
    ignore(bad)
  } catch {
    _ => raised = true
  }
  assert_true(raised)

  raised = false
  try {
    let bad = @src.codon_optimizer_config(gc_min=0.8, gc_max=0.2)
    ignore(bad)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "codon optimizer revcomp and iupac matching" {
  assert_eq(@src.codon_optimizer_revcomp("GAATTC"), "GAATTC")
  assert_eq(@src.codon_optimizer_revcomp("ACGTN"), "NACGT")
  assert_true(@src.codon_optimizer_motif_matches("ACGT", 0, "NNNN"))
  assert_true(@src.codon_optimizer_motif_matches("AGCT", 0, "ARYW"))
  assert_true(!@src.codon_optimizer_motif_matches("ACGT", 0, "GAAT"))
  assert_true(!@src.codon_optimizer_motif_matches("ACGT", 2, "ACGTT"))
}

///|
test "codon optimizer forbidden sites from enzymes" {
  let sites = @src.codon_optimizer_sites_for_enzymes([
    "EcoRI", "BamHI", "NotAnEnzyme",
  ])
  assert_true(sites.contains("GAATTC"))
  assert_true(sites.contains("GGATCC"))
  assert_eq(sites.length(), 2)
}

///|
test "codon optimizer reverse translates protein exactly" {
  let weights = codon_opt_test_weights(
    ["ATG", "AAA", "AAG", "TAA", "TAG", "TGA"],
    1.0,
  )
  let config = @src.codon_optimizer_config(seed=7UL)
  let dna = @src.codon_optimize("MKK*", weights, config~)
  assert_eq(dna.length(), 12)
  assert_eq(codon_opt_test_translate(dna), "MKK*")
  assert_true(dna.has_prefix("ATG"))
}

///|
test "codon optimizer enforces GC bounds" {
  // Ala 的同义密码子为 GCN; 高 GC 目标迫使选择 GCC/GCG.
  let weights = codon_opt_test_weights(["GCA", "GCC", "GCG", "GCT"], 1.0)
  let config = @src.codon_optimizer_config(gc_min=0.9, seed=3UL)
  let dna = @src.codon_optimize("AAAAA", weights, config~)
  assert_eq(codon_opt_test_translate(dna), "AAAAA")
  let report = @src.analyze_codon_optimization(dna, weights)
  assert_true(report.gc >= 0.9)

  let config_low = @src.codon_optimizer_config(gc_max=0.7, seed=3UL)
  let dna_low = @src.codon_optimize("AAAAA", weights, config=config_low)
  let report_low = @src.analyze_codon_optimization(dna_low, weights)
  assert_true(report_low.gc <= 0.7)
}

///|
test "codon optimizer caps homopolymer runs" {
  let weights = codon_opt_test_weights(["AAA", "AAG"], 1.0)
  let config = @src.codon_optimizer_config(max_run=4, seed=11UL)
  let dna = @src.codon_optimize("KKKKKKKK", weights, config~)
  assert_eq(codon_opt_test_translate(dna), "KKKKKKKK")
  let report = @src.analyze_codon_optimization(dna, weights)
  assert_true(report.max_run <= 4)
}

///|
test "codon optimizer avoids forbidden motif" {
  // 禁止 AAA: 需在 Lys 的 AAA/AAG 中选择 AAG.
  let weights = codon_opt_test_weights(["AAA", "AAG"], 1.0)
  let config = @src.codon_optimizer_config(forbidden_motifs=["AAA"], seed=5UL)
  let dna = @src.codon_optimize("KKK", weights, config~)
  assert_true(!dna.contains("AAA"))
  let report = @src.analyze_codon_optimization(dna, weights, forbidden_motifs=[
    "AAA",
  ])
  assert_eq(report.forbidden_count, 0)
}

///|
test "codon optimizer rejects unknown residues" {
  let mut raised = false
  try {
    let empty : Map[String, Double] = Map([], capacity=0)
    let dna = @src.codon_optimize(
      "MZ",
      empty,
      config=@src.codon_optimizer_config(),
    )
    ignore(dna)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "codon optimizer avoids cross-codon CpG" {
  // Gly 的 GGC(末位 C) + Glu 的 GAG/GAA(首位 G) 拼接产生跨界 CpG(GGCGAG).
  // 权重偏好这两个密码子, 逼迫优化器避开跨界组合.
  let weights : Map[String, Double] = Map([], capacity=0)
  weights["GGC"] = 2.0
  weights["GAG"] = 2.0
  let config = @src.codon_optimizer_config(avoid_cpg=true, seed=5UL)
  let dna = @src.codon_optimize("GE", weights, config~)
  assert_eq(codon_opt_test_translate(dna), "GE")
  assert_true(!dna.to_upper().contains("CG"))
  // 更长的蛋白质同样保持无 CpG.
  let dna_long = @src.codon_optimize("GEGEGE", weights, config~)
  assert_eq(codon_opt_test_translate(dna_long), "GEGEGE")
  assert_true(!dna_long.to_upper().contains("CG"))
}

///|
test "codon optimizer caps cross-codon homopolymer runs" {
  // Glu 的 GAA(尾部 AA, run=2) + Thr 的 ACN(首位 A) 拼接产生跨界 AAA run=3.
  let weights = codon_opt_test_weights(
    ["GAA", "GAG", "ACT", "ACC", "ACA", "ACG"],
    1.0,
  )
  let config = @src.codon_optimizer_config(max_run=2, seed=5UL)
  let dna = @src.codon_optimize("ET", weights, config~)
  assert_eq(codon_opt_test_translate(dna), "ET")
  let report = @src.analyze_codon_optimization(dna, weights)
  assert_true(report.max_run <= 2)
}
