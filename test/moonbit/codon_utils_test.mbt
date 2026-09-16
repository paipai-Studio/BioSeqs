///|
/// codon_utils.mbt 测试
/// frozen values 由 Python 精确生成并验证

fn _cu_eq(a : Double, b : Double) -> Bool {
  let tol = 0.0001
  let d = a - b
  if d < 0.0 { -d <= tol } else { d <= tol }
}

// ============ 测试序列 ============
let _GLOBIN : String = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG"
let _SIMPLE : String = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
let _REF_SET : Array[String] = [
  "ATGGCCAAGCACGTGAGCGTGGAGAAGGTGGCCGAGTGCAAGAGCACCGGCACCCTGAACGGCACCGAGAAGGTGGAGGTGCTGGACAC",
  "ATGGCCTACATCGGCTGCGTGCTGAACAAGCACGAGAAGGACGACGTGATGGTCGTCAACGACCCGTGCGTGACCCACAAGCCCGAG",
  "ATGGCCGACAAGACCAAGGGCTGCGTGAAGCTGATCGTGAACGACGTGTTTGACCGGCTGATCGGCTTCGAGCACTGCTGCGCTGAC",
]
let _TEST_GENE : String = "ATGGCCAAGATGATCGGCTACGTGCAGCGCGGCGTCGACGCCACCCTGGCCGAGATGCGCAAGGCCGGCATCACCGTGGCC"

// ============ reverse_complement 验证 ============

test "rc_basic" {
  assert_true(@src.reverse_complement("ATGGCC") == "GGCCAT")
  assert_true(@src.reverse_complement("AAAA") == "TTTT")
  assert_true(@src.reverse_complement("ATCG") == "CGAT")
}

// ============ GC123 ============

test "gc123_globin" {
  let gc = @src.seq_gc123(_GLOBIN)
  assert_true(_cu_eq(gc[0], 0.733333))
  assert_true(_cu_eq(gc[1], 0.400000))
  assert_true(_cu_eq(gc[2], 0.633333))
}

test "gc123_simple" {
  let gc = @src.seq_gc123(_SIMPLE)
  assert_true(_cu_eq(gc[0], 0.538462))
  assert_true(_cu_eq(gc[1], 0.538462))
  assert_true(_cu_eq(gc[2], 0.615385))
}

test "gc123_test_gene" {
  let gc = @src.seq_gc123(_TEST_GENE)
  assert_true(_cu_eq(gc[0], 0.629630))
  assert_true(_cu_eq(gc[1], 0.444444))
  assert_true(_cu_eq(gc[2], 1.000000))
}

test "gc123_empty" {
  let gc = @src.seq_gc123("")
  assert_true(gc[0] == 0.0 && gc[1] == 0.0 && gc[2] == 0.0)
}

test "gc123_short" {
  let gc = @src.seq_gc123("AT")
  assert_true(gc[0] == 0.0 && gc[1] == 0.0 && gc[2] == 0.0)
}

test "gc123_single_codon" {
  let gc = @src.seq_gc123("ATG")
  assert_true(gc[0] == 0.0)  // A
  assert_true(gc[1] == 0.0)  // T
  assert_true(gc[2] == 1.0)  // G
}

// ============ 六框翻译 ============

test "six_frame_globin_fwd" {
  let frames = @src.seq_six_frame_translations(_GLOBIN)
  assert_true(frames[0].frame_id == "Fwd-1")
  assert_true(frames[0].strand == "+")
  assert_true(frames[0].protein == "MVHLTPEEKSAVTALWGKVNVDEVGGEALG")

  assert_true(frames[1].frame_id == "Fwd-2")
  assert_true(frames[1].protein == "WCI*LLRRSLPLLPCGAR*TWMKLVVRPWA")

  assert_true(frames[2].frame_id == "Fwd-3")
  assert_true(frames[2].protein == "GASDS*GEVCRYCPVGQGERG*SWW*GPGQ")
}

test "six_frame_globin_rev" {
  let frames = @src.seq_six_frame_translations(_GLOBIN)
  assert_true(frames[3].frame_id == "Rev-1")
  assert_true(frames[3].strand == "-")
  assert_true(frames[3].protein == "LPRASPPTSSTFTLPHRAVTADFSSGVRCT")

  assert_true(frames[4].frame_id == "Rev-2")
  assert_true(frames[4].protein == "CPGPHHQLHPRSPCPTGQ*RQTSPQESDAP")

  assert_true(frames[5].frame_id == "Rev-3")
  assert_true(frames[5].protein == "AQGLTTNFIHVHLAPQGSNGRLLLRSQMHH")
}

test "six_frame_simple_fwd1" {
  let frames = @src.seq_six_frame_translations(_SIMPLE)
  assert_true(frames[0].protein == "MAIVMGR*KGAR*")
}

test "six_frame_test_gene_all6" {
  let frames = @src.seq_six_frame_translations(_TEST_GENE)
  assert_true(frames.length() == 6)
  assert_true(frames[0].frame_id == "Fwd-1")
  assert_true(frames[0].protein == "MAKMIGYVQRGVDATLAEMRKAGITVA")
  assert_true(frames[3].frame_id == "Rev-1")
  assert_true(frames[3].protein == "GHGDAGLAHLGQGGVDAALHVADHLGH")
}

// ============ 密码子表切换 ============

test "codon_table_mito_table2" {
  let seq_mito = "ATGATACGATGAGCCACTGA"
  let f1 = @src.seq_six_frame_translations(seq_mito, table=1)
  assert_true(f1[0].protein == "MIR*AT")
  let f2 = @src.seq_six_frame_translations(seq_mito, table=2)
  assert_true(f2[0].protein == "MMRWAT")
}

test "codon_table_yeast_mito_table3" {
  let f1 = @src.seq_six_frame_translations("ATGCTG", table=1)
  assert_true(f1[0].protein == "ML")
  let f3 = @src.seq_six_frame_translations("ATGCTG", table=3)
  assert_true(f3[0].protein == "MT")
}

test "codon_table_accessor" {
  let t1 = @src.codon_table(1)
  assert_true(t1["ATG"] == "M")
  assert_true(t1["TAA"] == "*")
  assert_true(t1["CTG"] == "L")

  let t2 = @src.codon_table(2)
  assert_true(t2["ATA"] == "M")
  assert_true(t2["AGA"] == "*")
  assert_true(t2["TGA"] == "W")

  let t3 = @src.codon_table(3)
  assert_true(t3["CTG"] == "T")
  assert_true(t3["ATA"] == "M")
}

// ============ CAI 权重 ============

test "cai_weights_globin_self" {
  let w = @src.cai_weights([_GLOBIN])
  assert_true(_cu_eq(w["ATG"], 1.0))
  assert_true(w.length() > 0)
}

test "cai_weights_ref_set" {
  let w = @src.cai_weights(_REF_SET)
  assert_true(_cu_eq(w["ATG"], 1.0))
  let keys = w.keys().to_array()
  let n = keys.length()
  let mut i = 0
  while i < n {
    let v = w[keys[i]]
    assert_true(v >= 0.0 && v <= 1.0 + 0.000001)
    i = i + 1
  }
}

test "cai_weights_empty_ref" {
  let w = @src.cai_weights([])
  assert_true(w["ATG"] == 1.0)
  assert_true(w["TTT"] == 1.0)
}

// ============ CAI ============

test "cai_globin_self" {
  let w = @src.cai_weights([_GLOBIN])
  let cai = @src.codon_adaptation_index(_GLOBIN, w)
  assert_true(_cu_eq(cai, 0.938332))
}

test "cai_simple_self" {
  let w = @src.cai_weights([_SIMPLE])
  let cai = @src.codon_adaptation_index(_SIMPLE, w)
  assert_true(_cu_eq(cai, 1.0))
}

test "cai_test_gene_with_ref" {
  let w = @src.cai_weights(_REF_SET)
  let cai = @src.codon_adaptation_index(_TEST_GENE, w)
  assert_true(_cu_eq(cai, 0.934083))
}

test "cai_empty" {
  let w = @src.cai_weights(["ATG"])
  let cai = @src.codon_adaptation_index("", w)
  assert_true(cai == 0.0)
}

test "cai_short" {
  let w = @src.cai_weights(["ATG"])
  let cai = @src.codon_adaptation_index("A", w)
  assert_true(cai == 0.0)
}

// ============ Nc/ENC ============

test "enc_globin" {
  let nc = @src.seq_enc_wright(_GLOBIN)
  assert_true(_cu_eq(nc, 44.5))
}

test "enc_simple" {
  let nc = @src.seq_enc_wright(_SIMPLE)
  assert_true(_cu_eq(nc, 56.0))
}

test "enc_test_gene" {
  let nc = @src.seq_enc_wright(_TEST_GENE)
  assert_true(_cu_eq(nc, 41.0))
}

test "enc_diverse_sequence" {
  let diverse = "ATGGCCAAGATGATCGGCTACGTGCAGCGCGGCGTCGACGCCACCCTGGCCGAGATGCGCAAGGCCGGCATCACCGTGGCCAAAGTTTGGGCA"
  let nc = @src.seq_enc_wright(diverse)
  assert_true(nc > 20.0 && nc < 61.0)
}

// ============ 六框结构验证 ============

test "six_frame_structure" {
  let frames = @src.seq_six_frame_translations("ATG")
  assert_true(frames.length() == 6)
  assert_true(frames[0].strand == "+")
  assert_true(frames[0].frame_offset == 0)
  assert_true(frames[3].strand == "-")
  assert_true(frames[3].frame_offset == 0)
  assert_true(frames[0].frame_id == "Fwd-1")
  assert_true(frames[5].frame_id == "Rev-3")
}

test "six_frame_short_seq" {
  let frames = @src.seq_six_frame_translations("ATG")
  assert_true(frames[0].protein == "M")
  assert_true(frames[1].protein == "")
  assert_true(frames[2].protein == "")
}

// ============ CAI 边界 ============

test "cai_single_codon_ref" {
  let w = @src.cai_weights(["ATGGCC"])
  assert_true(_cu_eq(w["ATG"], 1.0))
  assert_true(_cu_eq(w["GCC"], 1.0))
  assert_true(w.contains("TTT"))
}
