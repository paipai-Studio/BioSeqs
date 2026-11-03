// ===========================================================================
// vdj_test.mbt — 免疫组库分析测试(V(D)J 指认 / CDR3 / 克隆型 / 多样性 / TCRdist)
// ===========================================================================

///|
fn vdj_test_db() -> @src.VdjGeneDb {
  let v : Array[@src.VdjGene] = [
    {
      name: "TRBV01",
      locus: "TRB",
      kind: @src.VGene,
      sequence: "AAAACTGGGGTTTACCTGT",
    },
    {
      name: "TRBV02",
      locus: "TRB",
      kind: @src.VGene,
      sequence: "CCCCATCGGGGAACGATGT",
    },
  ]
  let d : Array[@src.VdjGene] = [
    {
      name: "TRBD1",
      locus: "TRB",
      kind: @src.DGene,
      sequence: "GCGCGATTATCAGGCA",
    },
    {
      name: "TRBD2",
      locus: "TRB",
      kind: @src.DGene,
      sequence: "TTATCGGATTAAGCCT",
    },
  ]
  let j : Array[@src.VdjGene] = [
    {
      name: "TRBJ01",
      locus: "TRB",
      kind: @src.JGene,
      sequence: "TGGGGCCATTTGAAAGTGG",
    },
    {
      name: "TRBJ02",
      locus: "TRB",
      kind: @src.JGene,
      sequence: "TGGGGTTAAACCGTTAGG",
    },
  ]
  @src.vdj_build_gene_db(v, d, j, k=6)
}

///|
test "vdj assign full productive read" {
  let db = vdj_test_db()
  // read = TRBV01 全段 + N 区(ACG) + TRBJ01 全段
  // CDR3 = TGT + ACG + TGG -> CTW
  let read = "AAAACTGGGGTTTACCTGT" + "ACG" + "TGGGGCCATTTGAAAGTGG"
  let chain = @src.vdj_assign(read, db, sequence_id="r1")
  assert_eq(chain.v_call, "TRBV01")
  assert_eq(chain.j_call, "TRBJ01")
  assert_eq(chain.cdr3_nt, "TGTACGTGG")
  assert_eq(chain.cdr3_aa, "CTW")
  assert_true(chain.productive)
  assert_true(chain.in_frame)
  assert_false(chain.stop_codon)
}

///|
test "vdj assign with D segment in N region" {
  let db = vdj_test_db()
  // N 区 = D 段(16nt) + 2nt P 核苷酸 -> CDR3 = 3+16+2+3 = 24nt 框内
  let read = "AAAACTGGGGTTTACCTGT" +
    "GCGCGATTATCAGGCA" +
    "AC" +
    "TGGGGCCATTTGAAAGTGG"
  let chain = @src.vdj_assign(read, db, sequence_id="r2")
  assert_eq(chain.v_call, "TRBV01")
  assert_eq(chain.j_call, "TRBJ01")
  assert_eq(chain.d_call, "TRBD1")
  assert_eq(chain.cdr3_nt, "TGTGCGCGATTATCAGGCAACTGG")
  assert_eq(chain.cdr3_aa, "CARLSGNW")
  assert_eq(chain.cdr3_aa.length(), 8)
  assert_true(chain.productive)
}

///|
test "vdj assign second V and J genes" {
  let db = vdj_test_db()
  let read = "CCCCATCGGGGAACGATGT" + "TTA" + "TGGGGTTAAACCGTTAGG"
  let chain = @src.vdj_assign(read, db, sequence_id="r3")
  assert_eq(chain.v_call, "TRBV02")
  assert_eq(chain.j_call, "TRBJ02")
  assert_eq(chain.cdr3_nt, "TGTTTATGG")
  assert_eq(chain.cdr3_aa, "CLW")
}

///|
test "vdj extract cdr3 standalone" {
  // 无数据库模式: 直接从读序提取
  let read = "GACACAGCGGCCACCTGGTGTACGTGGGGCCATTTGAAAG"
  let cdr3 = @src.vdj_extract_cdr3_nt(read).unwrap()
  assert_eq(cdr3, "TGTACGTGG")
  let (aa, in_frame, stop, prod) = @src.vdj_cdr3_aa(cdr3)
  assert_eq(aa, "CTW")
  assert_true(in_frame)
  assert_false(stop)
  assert_true(prod)
}

///|
test "vdj unproductive out of frame" {
  let db = vdj_test_db()
  // N 区 4nt -> CDR3 长度 10, 非框内
  let read = "AAAACTGGGGTTTACCTGT" + "ACGT" + "TGGGGCCATTTGAAAGTGG"
  let chain = @src.vdj_assign(read, db, sequence_id="r4")
  assert_false(chain.in_frame)
  assert_false(chain.productive)
  assert_eq(chain.cdr3_nt, "TGTACGTTGG")
}

///|
test "vdj unproductive stop codon" {
  let db = vdj_test_db()
  // N 区 TAA -> 终止密码子
  let read = "AAAACTGGGGTTTACCTGT" + "TAA" + "TGGGGCCATTTGAAAGTGG"
  let chain = @src.vdj_assign(read, db, sequence_id="r5")
  assert_true(chain.stop_codon)
  assert_false(chain.productive)
}

///|
test "vdj gene family strips allele" {
  assert_eq(@src.vdj_gene_family("TRBV05*01"), "TRBV05")
  assert_eq(@src.vdj_gene_family("TRBJ02"), "TRBJ02")
  assert_eq(@src.vdj_gene_family(""), "")
}

///|
test "vdj clonotypes grouping and ordering" {
  let db = vdj_test_db()
  let read_a = "AAAACTGGGGTTTACCTGT" + "ACG" + "TGGGGCCATTTGAAAGTGG"
  let read_b = "CCCCATCGGGGAACGATGT" + "TTA" + "TGGGGTTAAACCGTTAGG"
  // read_c 与 read_a 同克隆型(相同 CDR3+V+J)
  let read_c = read_a
  // read_d: 不同 CDR3
  let read_d = "AAAACTGGGGTTTACCTGT" + "GCA" + "TGGGGCCATTTGAAAGTGG"
  let reads : Array[(String, String)] = [
    (read_a, "r1"),
    (read_b, "r2"),
    (read_c, "r3"),
    (read_d, "r4"),
  ]
  let chains = @src.vdj_assign_all(reads, db)
  let clonotypes = @src.vdj_clonotypes(chains)
  assert_eq(clonotypes.length(), 3)
  // 最大克隆型是 CTW|TRBV01|TRBJ01, 计数 2
  let top = clonotypes[0]
  assert_eq(top.clonotype_id, 1)
  assert_eq(top.cdr3_aa, "CTW")
  assert_eq(top.count, 2)
  assert_eq(top.member_ids, ["r1", "r3"])
  let freq = top.frequency
  assert_true(freq > 0.49 && freq < 0.51)
}

///|
test "vdj diversity metrics" {
  let db = vdj_test_db()
  let read_a = "AAAACTGGGGTTTACCTGT" + "ACG" + "TGGGGCCATTTGAAAGTGG"
  let read_b = "CCCCATCGGGGAACGATGT" + "TTA" + "TGGGGTTAAACCGTTAGG"
  let read_d = "AAAACTGGGGTTTACCTGT" + "GCA" + "TGGGGCCATTTGAAAGTGG"
  // 单克隆型扩张: 4 条同 read_a
  let reads_expanded : Array[(String, String)] = [
    (read_a, "e1"),
    (read_a, "e2"),
    (read_a, "e3"),
    (read_a, "e4"),
  ]
  // 均匀分布: 三种各一条
  let reads_even : Array[(String, String)] = [
    (read_a, "u1"),
    (read_b, "u2"),
    (read_d, "u3"),
  ]
  let chains_exp = @src.vdj_assign_all(reads_expanded, db)
  let chains_even = @src.vdj_assign_all(reads_even, db)
  let div_exp = @src.vdj_diversity(chains_exp)
  let div_even = @src.vdj_diversity(chains_even)
  // 扩张样本: 1 个克隆型, clonality=1(极值), shannon=0
  assert_eq(div_exp.total_clonotypes, 1)
  assert_true(div_exp.shannon < 1.0e-9)
  assert_true(div_exp.clonality > 0.99)
  assert_true(div_exp.simpson < 1.0e-9)
  assert_eq(div_exp.d50, 1)
  // 均匀样本: 3 克隆型
  assert_eq(div_even.total_clonotypes, 3)
  let ln3 = 1.0986122886681098
  assert_true(div_even.shannon > ln3 - 0.001 && div_even.shannon < ln3 + 0.001)
  assert_true(div_even.clonality < 0.001)
  assert_true(div_even.simpson > 0.66 && div_even.simpson < 0.67)
  assert_eq(div_even.d50, 2)
  // Gini: 扩张(单一)=0, 均匀=0
  assert_true(div_exp.gini < 1.0e-9)
  assert_true(div_even.gini < 1.0e-9)
}

///|
test "vdj cdr3 distance" {
  // 相同序列距离为 0
  assert_true(@src.vdj_cdr3_distance("CASSLFW", "CASSLFW") < 1.0e-9)
  // 相似替换(F->Y 同芳香组)距离小于非相似替换
  let d_sim = @src.vdj_cdr3_distance("CASFLFW", "CASYLFW")
  let d_diff = @src.vdj_cdr3_distance("CASFLFW", "CASKLFW")
  assert_true(d_sim < d_diff)
  // gap 罚分: 更长的差异距离更大
  let d_gap1 = @src.vdj_cdr3_distance("CASS", "CASSA")
  let d_gap2 = @src.vdj_cdr3_distance("CASS", "CASSAA")
  assert_true(d_gap2 > d_gap1)
  // 空序列
  assert_true(@src.vdj_cdr3_distance("", "") < 1.0e-9)
}

///|
test "vdj tcrdist matrix symmetric" {
  let aas : Array[String] = ["CASSLFW", "CASSYFW", "CATKLFW"]
  let m = @src.vdj_tcrdist_matrix(aas)
  assert_eq(m.length(), 3)
  assert_eq(m[0].length(), 3)
  assert_true(m[0][0] < 1.0e-9)
  assert_true(m[1][1] < 1.0e-9)
  for i = 0; i < 3; i = i + 1 {
    for j = 0; j < 3; j = j + 1 {
      let a = m[i][j]
      let b = m[j][i]
      let diff = if a > b { a - b } else { b - a }
      assert_true(diff < 1.0e-9)
    }
  }
}

///|
test "vdj default databases buildable" {
  let trb = @src.vdj_default_trb_db()
  let igh = @src.vdj_default_igh_db()
  ignore(trb)
  ignore(igh)
  assert_true(true)
}

///|
test "vdj assign reads with mismatches tolerated" {
  let db = vdj_test_db()
  // 引入 2 个错配(非种子区): V 段中部 A->T, T->G
  let read = "AAATCTGGGGTTTGCCTGT" + "ACG" + "TGGGGCCATTTGAAAGTGG"
  let chain = @src.vdj_assign(read, db, sequence_id="mm")
  assert_eq(chain.v_call, "TRBV01")
  assert_true(chain.v_identity > 0.85 && chain.v_identity < 1.0)
}

///|
test "vdj nonproductive excluded from clonotypes" {
  let db = vdj_test_db()
  let good = "AAAACTGGGGTTTACCTGT" + "ACG" + "TGGGGCCATTTGAAAGTGG"
  // 非框内读序
  let bad = "AAAACTGGGGTTTACCTGT" + "ACGT" + "TGGGGCCATTTGAAAGTGG"
  let reads : Array[(String, String)] = [(good, "g1"), (bad, "b1"), (bad, "b2")]
  let chains = @src.vdj_assign_all(reads, db)
  let clonotypes = @src.vdj_clonotypes(chains)
  // 仅 1 个 productive 克隆型, 频率按 3 条链总量计 = 1/3
  assert_eq(clonotypes.length(), 1)
  assert_eq(clonotypes[0].count, 1)
  let freq = clonotypes[0].frequency
  assert_true(freq > 0.32 && freq < 0.34)
}
