///|
/// codon_utils.mbt — NCBI 密码子表扩展测试
/// frozen values 由 Python Bio.Data.CodonTable 1.88 权威生成并验证
/// 覆盖: 27 个 NCBI 表 (1-6, 9-16, 21-33)

// ============ 表可用性 ============

test "codon_table_valid_ids_count" {
  let ids = @src.codon_table_available_ids()
  assert_true(ids.length() == 27)
}

test "codon_table_is_valid_true" {
  assert_true(@src.codon_table_is_valid(1))
  assert_true(@src.codon_table_is_valid(2))
  assert_true(@src.codon_table_is_valid(6))
  assert_true(@src.codon_table_is_valid(12))
  assert_true(@src.codon_table_is_valid(22))
  assert_true(@src.codon_table_is_valid(33))
}

test "codon_table_is_valid_false" {
  assert_false(@src.codon_table_is_valid(7))
  assert_false(@src.codon_table_is_valid(8))
  assert_false(@src.codon_table_is_valid(17))
  assert_false(@src.codon_table_is_valid(99))
}

test "codon_table_default_for_invalid_id" {
  // 无效 ID 应返回 Table 1
  let t1 = @src.codon_table(1)
  let t99 = @src.codon_table(99)
  // Table 1 中 TGA 是终止密码子
  let v1 = t1["TGA"]
  let v99 = t99["TGA"]
  assert_true(v1 == "*")
  assert_true(v99 == "*")
}

// ============ 表 1: Standard ============

test "table1_basic_codons" {
  let t = @src.codon_table(1)
  assert_true(t["ATG"] == "M")
  assert_true(t["TTT"] == "F")
  assert_true(t["CTG"] == "L")
  assert_true(t["TGG"] == "W")
}

test "table1_stop_codons" {
  let stops = @src.codon_table_stop_codons(1)
  assert_true(stops.length() == 3)
  assert_true(stops.contains("TAA"))
  assert_true(stops.contains("TAG"))
  assert_true(stops.contains("TGA"))
}

test "table1_start_codons" {
  let starts = @src.codon_table_start_codons(1)
  assert_true(starts.length() == 3)
  assert_true(starts.contains("ATG"))
  assert_true(starts.contains("TTG"))
  assert_true(starts.contains("CTG"))
}

test "table1_name" {
  assert_true(@src.codon_table_name(1) == "Standard")
}

test "table1_reverse_met" {
  let rev = @src.codon_table_reverse(table_id=1)
  // Met 只有 ATG
  let met_codons = rev["M"]
  assert_true(met_codons.length() == 1)
  assert_true(met_codons.contains("ATG"))
}

// ============ 表 2: Vertebrate Mitochondrial ============

test "table2_variants" {
  let t = @src.codon_table(2)
  assert_true(t["ATA"] == "M")   // Ile -> Met
  assert_true(t["AGA"] == "*")   // Arg -> Stop
  assert_true(t["AGG"] == "*")   // Arg -> Stop
  assert_true(t["TGA"] == "W")   // Stop -> Trp
}

test "table2_stop_codons" {
  let stops = @src.codon_table_stop_codons(2)
  assert_true(stops.length() == 4)
  assert_true(stops.contains("AGA"))
  assert_true(stops.contains("AGG"))
  assert_true(stops.contains("TAA"))
  assert_true(stops.contains("TAG"))
}

test "table2_start_codons" {
  let starts = @src.codon_table_start_codons(2)
  assert_true(starts.length() == 5)
  assert_true(starts.contains("ATA"))
  assert_true(starts.contains("GTG"))
}

// ============ 表 3: Yeast Mitochondrial ============

test "table3_variants" {
  let t = @src.codon_table(3)
  assert_true(t["ATA"] == "M")
  assert_true(t["CTA"] == "T")   // Leu -> Thr
  assert_true(t["CTC"] == "T")
  assert_true(t["CTG"] == "T")
  assert_true(t["CTT"] == "T")
  assert_true(t["TGA"] == "W")
}

// ============ 表 5: Invertebrate Mitochondrial ============

test "table5_variants" {
  let t = @src.codon_table(5)
  assert_true(t["AGA"] == "S")   // Arg -> Ser
  assert_true(t["AGG"] == "S")
  assert_true(t["ATA"] == "M")
  assert_true(t["TGA"] == "W")
}

// ============ 表 6: Ciliate Nuclear (只有 TGA 终止) ============

test "table6_variants" {
  let t = @src.codon_table(6)
  assert_true(t["TAA"] == "Q")   // Stop -> Gln
  assert_true(t["TAG"] == "Q")   // Stop -> Gln
  assert_true(t["TGA"] == "*")   // 唯一终止
}

test "table6_stop_codons" {
  let stops = @src.codon_table_stop_codons(6)
  assert_true(stops.length() == 1)
  assert_true(stops.contains("TGA"))
}

test "table6_start_codons" {
  let starts = @src.codon_table_start_codons(6)
  assert_true(starts.length() == 1)
  assert_true(starts.contains("ATG"))
}

// ============ 表 10: Euplotid Nuclear ============

test "table10_variants" {
  let t = @src.codon_table(10)
  assert_true(t["TGA"] == "C")   // Stop -> Cys
}

// ============ 表 12: Alternative Yeast Nuclear ============

test "table12_variants" {
  let t = @src.codon_table(12)
  assert_true(t["CTG"] == "S")   // Leu -> Ser
}

test "table12_start_codons" {
  let starts = @src.codon_table_start_codons(12)
  assert_true(starts.contains("CTG"))
  assert_true(starts.contains("ATG"))
}

// ============ 表 14: Alternative Flatworm Mitochondrial ============

test "table14_variants" {
  let t = @src.codon_table(14)
  assert_true(t["AAA"] == "N")   // Lys -> Asn
  assert_true(t["AGA"] == "S")
  assert_true(t["AGG"] == "S")
  assert_true(t["TAA"] == "Y")   // Stop -> Tyr
  assert_true(t["TAG"] == "*")   // 唯一终止
  assert_true(t["TGA"] == "W")
}

test "table14_stop_codons" {
  let stops = @src.codon_table_stop_codons(14)
  assert_true(stops.length() == 1)
  assert_true(stops.contains("TAG"))
}

// ============ 表 15: Blepharisma Macronuclear ============

test "table15_variants" {
  let t = @src.codon_table(15)
  assert_true(t["TAG"] == "Q")   // Stop -> Gln
}

// ============ 表 16: Chlorophycean Mitochondrial ============

test "table16_variants" {
  let t = @src.codon_table(16)
  assert_true(t["TAG"] == "L")   // Stop -> Leu
}

// ============ 表 22: Scenedesmus obliquus Mitochondrial ============

test "table22_variants" {
  let t = @src.codon_table(22)
  assert_true(t["TAG"] == "L")   // Stop -> Leu
  assert_true(t["TCA"] == "*")   // Ser -> Stop (独特!)
}

test "table22_name" {
  assert_true(@src.codon_table_name(22) == "Scenedesmus obliquus Mitochondrial")
}

test "table22_stop_codons" {
  let stops = @src.codon_table_stop_codons(22)
  assert_true(stops.contains("TCA"))
  assert_true(stops.contains("TAA"))
  assert_true(stops.contains("TGA"))
}

// ============ 表 23: Thraustochytrium Mitochondrial ============

test "table23_variants" {
  let t = @src.codon_table(23)
  assert_true(t["TTA"] == "*")   // Leu -> Stop (独特!)
}

test "table23_stop_codons" {
  let stops = @src.codon_table_stop_codons(23)
  assert_true(stops.length() == 4)
  assert_true(stops.contains("TTA"))
  assert_true(stops.contains("TAA"))
  assert_true(stops.contains("TAG"))
  assert_true(stops.contains("TGA"))
}

// ============ 表 25: Candidate Division SR1 ============

test "table25_variants" {
  let t = @src.codon_table(25)
  assert_true(t["TGA"] == "G")   // Stop -> Gly
}

test "table25_stop_codons" {
  let stops = @src.codon_table_stop_codons(25)
  assert_true(stops.length() == 2)
  assert_true(stops.contains("TAA"))
  assert_true(stops.contains("TAG"))
  assert_false(stops.contains("TGA"))
}

// ============ 表 26: Pachysolen tannophilus Nuclear ============

test "table26_variants" {
  let t = @src.codon_table(26)
  assert_true(t["CTG"] == "A")   // Leu -> Ala
}

// ============ 表 27/28: Karyorelict/Condylostoma Nuclear ============

test "table27_variants" {
  let t = @src.codon_table(27)
  assert_true(t["TAA"] == "Q")
  assert_true(t["TAG"] == "Q")
  assert_true(t["TGA"] == "W")
}

test "table28_variants" {
  let t = @src.codon_table(28)
  assert_true(t["TAA"] == "Q")
  assert_true(t["TAG"] == "Q")
}

// ============ 表 29: Mesodinium Nuclear ============

test "table29_variants" {
  let t = @src.codon_table(29)
  assert_true(t["TAA"] == "Y")   // Stop -> Tyr
  assert_true(t["TAG"] == "Y")   // Stop -> Tyr
  assert_true(t["TGA"] == "*")
}

test "table29_stop_codons" {
  let stops = @src.codon_table_stop_codons(29)
  assert_true(stops.length() == 1)
  assert_true(stops.contains("TGA"))
}

// ============ 表 30: Peritrich Nuclear ============

test "table30_variants" {
  let t = @src.codon_table(30)
  assert_true(t["TAA"] == "E")   // Stop -> Glu
  assert_true(t["TAG"] == "E")   // Stop -> Glu
  assert_true(t["TGA"] == "*")
}

// ============ 表 31: Blastocrithidia Nuclear ============

test "table31_variants" {
  let t = @src.codon_table(31)
  assert_true(t["TAA"] == "E")
  assert_true(t["TAG"] == "E")
  assert_true(t["TGA"] == "W")   // Stop -> Trp
}

// ============ 表 32: Balanophoraceae Plastid (新增) ============

test "table32_variants" {
  let t = @src.codon_table(32)
  assert_true(t["TAG"] == "W")   // Stop -> Trp
  assert_true(t["TAA"] == "*")
  assert_true(t["TGA"] == "*")
}

test "table32_stops" {
  let stops = @src.codon_table_stop_codons(32)
  assert_true(stops.length() == 2)
  assert_true(stops.contains("TAA"))
  assert_true(stops.contains("TGA"))
  assert_false(stops.contains("TAG"))
}

test "table32_starts" {
  let starts = @src.codon_table_start_codons(32)
  assert_true(starts.length() == 7)
  assert_true(starts.contains("ATA"))
  assert_true(starts.contains("GTG"))
}

// ============ 表 33: Cephalodiscidae Mitochondrial (新增) ============

test "table33_variants" {
  let t = @src.codon_table(33)
  assert_true(t["AGA"] == "S")
  assert_true(t["AGG"] == "K")   // Arg -> Lys
  assert_true(t["TAA"] == "Y")   // Stop -> Tyr
  assert_true(t["TAG"] == "*")   // 唯一终止
  assert_true(t["TGA"] == "W")
}

test "table33_stops" {
  let stops = @src.codon_table_stop_codons(33)
  assert_true(stops.length() == 1)
  assert_true(stops.contains("TAG"))
}

// ============ 反向密码子表 ============

test "reverse_table1_all_amino_acids_present" {
  let rev = @src.codon_table_reverse(table_id=1)
  // 20 种标准氨基酸
  let amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                     "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
  let mut i = 0
  while i < amino_acids.length() {
    let aa = amino_acids[i]
    assert_true(rev.contains(aa))
    i = i + 1
  }
}

test "reverse_table1_no_stop_codons" {
  let rev = @src.codon_table_reverse(table_id=1)
  assert_false(rev.contains("*"))
}

test "reverse_table6_contains_gln_from_stops" {
  // Table 6 的 TAA/TAG 编码 Q (Gln)
  let rev = @src.codon_table_reverse(table_id=6)
  let gln_codons = rev["Q"]
  assert_true(gln_codons.contains("TAA"))
  assert_true(gln_codons.contains("TAG"))
}

test "reverse_table25_no_tga_stop" {
  // Table 25 的 TGA 编码 G
  let rev = @src.codon_table_reverse(table_id=25)
  let gly_codons = rev["G"]
  assert_true(gly_codons.contains("TGA"))
}

// ============ 六框翻译 + 多密码子表 ============

test "six_frame_table2_mitochondrial" {
  // 用 Table 2，ATA 应翻译为 M (而不是 I)
  let frames = @src.seq_six_frame_translations("ATGATACGT", table=2)
  assert_true(frames.length() == 6)
  // Fwd-1 起始 ATG -> M, 然后 ATA -> M (Table 2), CGT -> R
  assert_true(frames[0].protein == "MMR")
}

test "six_frame_table6_ciliate" {
  // 用 Table 6，TGA 是终止密码子但 TAA/TAG 是 Q
  let frames = @src.seq_six_frame_translations("ATGTGATAA", table=6)
  assert_true(frames[0].protein == "M*Q")
}

// ============ 全部 27 个表可正常访问 ============

test "all_valid_tables_accessible" {
  let ids = @src.codon_table_available_ids()
  let mut i = 0
  while i < ids.length() {
    let tid = ids[i]
    let t = @src.codon_table(tid)
    // 每个表必须有 64 个密码子
    assert_true(t.keys().length() == 64)
    // 必须有 name
    let name = @src.codon_table_name(tid)
    assert_true(name.length() > 0)
    // 必须有 stops
    let stops = @src.codon_table_stop_codons(tid)
    assert_true(stops.length() > 0)
    i = i + 1
  }
}
