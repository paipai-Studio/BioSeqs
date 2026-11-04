// ===========================================================================
// motif_discovery_test.mbt — de novo motif 发现测试(EM/Gibbs)
// ===========================================================================

///|
/// 构造含植入 motif 的序列集(两侧随机背景).
fn motif_disco_test_seqs(motif : String, n : Int32) -> Array[String] {
  let rng = @src.Rng::from_int(12345)
  let bg_chars : Array[Char] = ['A', 'C', 'G', 'T']
  let seqs : Array[String] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let sb = StringBuilder()
    // 左侧 8-12 背景
    let left_len = 8 + rng.next_int(5)
    for j = 0; j < left_len; j = j + 1 {
      sb.write_char(bg_chars[rng.next_int(4)])
    }
    sb.write_string(motif)
    // 右侧 8-12 背景
    let right_len = 8 + rng.next_int(5)
    for j = 0; j < right_len; j = j + 1 {
      sb.write_char(bg_chars[rng.next_int(4)])
    }
    ignore(i)
    seqs.push(sb.to_string())
  }
  seqs
}

///|
test "motif em recovers planted motif" {
  let seqs = motif_disco_test_seqs("ACGTTGCA", 10)
  let result = @src.motif_em(seqs, width=8, iterations=50)
  // consensus 应恢复植入 motif(容许少量误差)
  let cons = result.consensus()
  let mut match_chars = 0
  for i = 0; i < 8; i = i + 1 {
    if cons.unsafe_get(i) == "ACGTTGCA".unsafe_get(i) {
      match_chars = match_chars + 1
    }
  }
  assert_true(match_chars >= 6)
  // 高 IC(强 motif)
  assert_true(result.information_content > 10.0)
  assert_eq(result.n_sequences, 10)
  assert_eq(result.width, 8)
}

///|
test "motif em site positions recovered" {
  let seqs = motif_disco_test_seqs("GGATCCAA", 8)
  let result = @src.motif_em(seqs, width=8, iterations=50)
  // 多数位点应被正确指认(位置在 8-12 区间, 双链均可)
  let mut correct = 0
  for p in result.site_positions {
    if p.0 >= 8 && p.0 <= 12 {
      correct = correct + 1
    }
  }
  assert_true(correct >= 6)
}

///|
test "motif gibbs recovers planted motif" {
  let seqs = motif_disco_test_seqs("TTAGGCAT", 10)
  let result = @src.motif_gibbs(seqs, width=8, restarts=3, sweeps=40, seed=42)
  // 双链语义: consensus 可为 motif 或其反向互补
  let cons = result.consensus()
  let cons_rc = @src.reverse_complement(cons)
  let planted = "TTAGGCAT"
  let mut match_chars = 0
  let mut match_rc = 0
  for i = 0; i < 8; i = i + 1 {
    if cons.unsafe_get(i) == planted.unsafe_get(i) {
      match_chars = match_chars + 1
    }
    if cons_rc.unsafe_get(i) == planted.unsafe_get(i) {
      match_rc = match_rc + 1
    }
  }
  assert_true(match_chars >= 6 || match_rc >= 6)
  assert_true(result.information_content > 8.0)
}

///|
test "motif gibbs deterministic with same seed" {
  let seqs = motif_disco_test_seqs("ACGTTGCA", 6)
  let r1 = @src.motif_gibbs(seqs, width=8, restarts=2, sweeps=20, seed=7)
  let r2 = @src.motif_gibbs(seqs, width=8, restarts=2, sweeps=20, seed=7)
  // 同种子 -> 同 IC 与同位点
  assert_eq(r1.information_content, r2.information_content)
  assert_eq(r1.consensus(), r2.consensus())
}

///|
test "motif em weak background no false structure" {
  // 纯随机序列: 发现的低 IC motif
  let rng = @src.Rng::from_int(99)
  let seqs : Array[String] = Array(capacity=0)
  let bg : Array[Char] = ['A', 'C', 'G', 'T']
  for i = 0; i < 10; i = i + 1 {
    let sb = StringBuilder()
    for j = 0; j < 30; j = j + 1 {
      sb.write_char(bg[rng.next_int(4)])
    }
    ignore(i)
    seqs.push(sb.to_string())
  }
  let result = @src.motif_em(seqs, width=8, iterations=30)
  assert_true(result.information_content < 12.0)
}

///|
test "motif em empty raises" {
  let mut raised = false
  try {
    let r = @src.motif_em([])
    ignore(r)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "motif em short sequence raises" {
  let mut raised = false
  try {
    let r = @src.motif_em(["ACG"], width=8)
    ignore(r)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "motif gibbs empty raises" {
  let mut raised = false
  try {
    let r = @src.motif_gibbs([])
    ignore(r)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "motif discovered pwm row sums to one" {
  let seqs = motif_disco_test_seqs("ACGTTGCA", 8)
  let result = @src.motif_em(seqs, width=8, iterations=30)
  for row in result.pwm {
    let mut s = 0.0
    for v in row {
      s = s + v
    }
    assert_true(s > 0.999 && s < 1.001)
  }
}

///|
test "motif discovered to motif conversion" {
  let seqs = motif_disco_test_seqs("ACGTTGCA", 6)
  let result = @src.motif_em(seqs, width=8, iterations=30)
  let m = result.to_motif()
  assert_eq(m.name, "discovered-EM")
  assert_eq(m.pwm.length(), 8)
  assert_eq(m.alphabet, ["A", "C", "G", "T"])
  assert_eq(m.sites.length(), 6)
}

///|
test "motif discovered summary" {
  let seqs = motif_disco_test_seqs("ACGTTGCA", 6)
  let result = @src.motif_em(seqs, width=8, iterations=30)
  let s = result.summary()
  assert_true(s.contains("motif-discovery"))
  assert_true(s.contains("IC="))
  assert_true(s.contains("consensus="))
}

///|
test "motif em reverse strand motif" {
  // 反向链植入: 植入 motif 的反向互补出现在序列中
  let planted = "ACGTTGCA"
  let rc = @src.reverse_complement(planted)
  let rng = @src.Rng::from_int(7)
  let bg : Array[Char] = ['A', 'C', 'G', 'T']
  let seqs : Array[String] = Array(capacity=0)
  for i = 0; i < 8; i = i + 1 {
    let sb = StringBuilder()
    for j = 0; j < 10; j = j + 1 {
      sb.write_char(bg[rng.next_int(4)])
    }
    sb.write_string(rc)
    for j = 0; j < 10; j = j + 1 {
      sb.write_char(bg[rng.next_int(4)])
    }
    ignore(i)
    seqs.push(sb.to_string())
  }
  let result = @src.motif_em(seqs, width=8, iterations=50)
  // consensus 应为原 motif(双链搜索语义)或其反向互补
  let cons = result.consensus()
  let cons_rc = @src.reverse_complement(cons)
  let mut match_direct = 0
  let mut match_rc = 0
  for i = 0; i < 8; i = i + 1 {
    if cons.unsafe_get(i) == planted.unsafe_get(i) {
      match_direct = match_direct + 1
    }
    if cons_rc.unsafe_get(i) == planted.unsafe_get(i) {
      match_rc = match_rc + 1
    }
  }
  assert_true(match_direct >= 6 || match_rc >= 6)
}
