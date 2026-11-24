///|
/// codon_table_custom.mbt — 统一遗传密码表对象测试
/// 覆盖: 27 表构造 / 名称 / 非法 ID / 双向索引 / 起始与终止密码子 /
///       IUPAC 简并解析 / 翻译入口 (translate, translate_frame, translate_cds)

// ============================================================
// 构造与元数据
// ============================================================

///|
test "ct_default_is_table1" {
  let t = @src.codon_table_default()
  assert_eq(t.id(), 1)
  assert_eq(t.name(), "Standard")
}

///|
test "ct_named_presets" {
  assert_eq(@src.codon_table_bacterial().id(), 11)
  assert_eq(@src.codon_table_alt_yeast().id(), 12)
  assert_eq(@src.codon_table_ciliate().id(), 6)
  assert_eq(@src.codon_table_vertebrate_mito().id(), 2)
}

///|
test "ct_all_27_tables_constructible" {
  let ids = @src.codon_table_available_ids()
  assert_eq(ids.length(), 27)
  let mut i = 0
  while i < ids.length() {
    let t = @src.codon_table_by_id(ids[i])
    assert_eq(t.id(), ids[i])
    assert_true(t.name().length() > 0)
    i = i + 1
  }
}

///|
test "ct_invalid_id_raises" {
  let raised = try {
    ignore(@src.codon_table_by_id(7))
    false
  } catch {
    CodonTableError(_) => true
  }
  assert_true(raised)
}

///|
test "ct_or_default_falls_back_to_table1" {
  let t = @src.codon_table_or_default(99)
  assert_eq(t.id(), 1)
  let t2 = @src.codon_table_or_default(12)
  assert_eq(t2.id(), 12)
}

///|
test "ct_from_map_is_custom_id_zero" {
  let m : Map[String, String] = Map([])
  m["ATG"] = "M"
  m["TAA"] = "*"
  let t = @src.codon_table_from_map(m, name="MyTable", start_codons=["ATG"], stop_codons=[
    "TAA",
  ])
  assert_eq(t.id(), 0)
  assert_eq(t.name(), "MyTable")
  assert_eq(t.codon_to_aa("ATG"), 'M')
  assert_eq(t.codon_to_aa("TAA"), '*')
}

// ============================================================
// 正向 / 反向索引
// ============================================================

///|
test "ct_forward_table_lookup" {
  let t = @src.codon_table_default()
  let f = t.forward_table()
  assert_eq(f["ATG"], "M")
  assert_eq(f["TAA"], "*")
  assert_eq(f.length(), 64)
}

///|
test "ct_forward_table_is_copy" {
  let t = @src.codon_table_default()
  let f = t.forward_table()
  f["ATG"] = "X"
  assert_eq(t.forward_table()["ATG"], "M")
}

///|
test "ct_backward_table_excludes_stop" {
  let t = @src.codon_table_default()
  let b = t.backward_table()
  assert_eq(b["M"].length(), 1)
  assert_true(b["M"].contains("ATG"))
  assert_eq(b["L"].length(), 6)
  assert_false(b.contains("*"))
}

///|
test "ct_backward_table_is_copy" {
  let t = @src.codon_table_default()
  let b = t.backward_table()
  match b.get("M") {
    Some(list) => list.push("AAA")
    None => ()
  }
  assert_eq(t.aa_to_codons("M").length(), 1)
}

///|
test "ct_aa_to_codons" {
  let t = @src.codon_table_default()
  assert_true(t.aa_to_codons("W").contains("TGG"))
  assert_eq(t.aa_to_codons("Z").length(), 0)
}

///|
test "ct_start_and_stop_codons" {
  let t1 = @src.codon_table_default()
  assert_eq(t1.stop_codons().length(), 3)
  assert_true(t1.start_codons().contains("ATG"))
  assert_true(t1.start_codons().contains("TTG"))
  assert_true(t1.start_codons().contains("CTG"))
  let t6 = @src.codon_table_ciliate()
  assert_eq(t6.stop_codons().length(), 1)
  assert_true(t6.stop_codons().contains("TGA"))
}

///|
test "ct_is_start_and_stop_codon" {
  let t = @src.codon_table_default()
  assert_true(t.is_start_codon("atg"))
  assert_true(t.is_stop_codon("TAA"))
  assert_false(t.is_stop_codon("AAA"))
}

///|
test "ct_canonical_codon_order" {
  let order = @src.canonical_codon_order()
  assert_eq(order.length(), 64)
  assert_eq(order[0], "TTT")
  assert_eq(order[63], "GGG")
}

// ============================================================
// 密码子解析 (含 IUPAC 简并)
// ============================================================

///|
test "ct_codon_to_aa_basic" {
  let t = @src.codon_table_default()
  assert_eq(t.codon_to_aa("ATG"), 'M')
  assert_eq(t.codon_to_aa("atg"), 'M')
  assert_eq(t.codon_to_aa("TAA"), '*')
  assert_eq(t.codon_to_aa("XXX"), 'X')
}

///|
test "ct_codon_to_aa_iupac_consistent" {
  let t = @src.codon_table_default()
  assert_eq(t.codon_to_aa("GCN"), 'A')
  assert_eq(t.codon_to_aa("ATH"), 'I')
  assert_eq(t.codon_to_aa("GAR"), 'E')
  assert_eq(t.codon_to_aa("TAR"), '*')
}

///|
test "ct_codon_to_aa_iupac_inconsistent" {
  let t = @src.codon_table_default()
  assert_eq(t.codon_to_aa("NNN"), 'X')
  assert_eq(t.codon_to_aa("NN"), 'X')
  assert_eq(t.codon_to_aa("1AT"), 'X')
}

///|
test "ct_codon_to_aa_uses_target_table" {
  // 表 2 中 TGA = W, 表 1 中 TGA = *
  assert_eq(@src.codon_table_vertebrate_mito().codon_to_aa("TGA"), 'W')
  assert_eq(@src.codon_table_default().codon_to_aa("TGA"), '*')
  // 表 12 中 CTG = S
  assert_eq(@src.codon_table_alt_yeast().codon_to_aa("CTG"), 'S')
}

// ============================================================
// 翻译
// ============================================================

///|
test "ct_translate" {
  let t = @src.codon_table_default()
  assert_eq(t.translate("ATGGCCGCTGAA"), "MAAE")
  assert_eq(t.translate("atggccgctgaa"), "MAAE")
  assert_eq(t.translate("ATGGCCGCGTAA"), "MAA*")
}

///|
test "ct_translate_frame" {
  let t = @src.codon_table_default()
  // 第 0 框: AAT GGC CGC TGA
  assert_eq(t.translate_frame("AATGGCCGCTGAA", 0), "NGR*")
  // 第 1 框: ATG GCC GCT GAA
  assert_eq(t.translate_frame("AATGGCCGCTGAA", 1), "MAAE")
}

///|
test "ct_translate_cds_stops_at_first_stop" {
  let t = @src.codon_table_default()
  assert_eq(t.translate_cds("ATGGCCGCGTAA"), "MAA*")
  assert_eq(t.translate_cds("ATGGCC"), "MA")
}

///|
test "ct_translate_cds_start_codon_to_met" {
  let t = @src.codon_table_default()
  // TTG 是表 1 的起始密码子, CDS 模式下统一翻译为 M
  assert_eq(t.translate_cds("TTGGCCTGCTAA"), "MAC*")
}

///|
test "ct_translate_with_table" {
  let t = @src.codon_table_vertebrate_mito()
  assert_eq(@src.translate_with_table("ATGTGA", t), "MW")
  assert_eq(
    @src.translate_with_table("ATGTGA", @src.codon_table_default()),
    "M*",
  )
}

///|
test "ct_translate_with_table_id" {
  assert_eq(@src.translate_with_table_id("ATGTGA", 2), "MW")
  assert_eq(@src.translate_with_table_id("ATGTGA", 1), "M*")
}

///|
test "ct_translate_with_table_id_invalid_raises" {
  let raised = try {
    ignore(@src.translate_with_table_id("ATG", 7))
    false
  } catch {
    CodonTableError(_) => true
  }
  assert_true(raised)
}

///|
test "ct_translate_with_codon_map" {
  let m : Map[String, String] = Map([])
  m["ATG"] = "M"
  m["TAA"] = "*"
  m["TGA"] = "W"
  assert_eq(@src.translate_with_codon_map("ATGTGA", m), "MW")
}

///|
test "ct_translate_cds_with_table_id" {
  assert_eq(@src.translate_cds_with_table_id("ATGGCCGCGTAA", 1), "MAA*")
  assert_eq(@src.translate_cds_with_table_id("ATGTGA", 2), "MW")
}
