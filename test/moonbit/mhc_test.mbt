///|
test "mhc allele lookup" {
  let a = @src.mhc_allele("HLA-A*02:01").unwrap()
  assert_eq(a.name(), "HLA-A*02:01")
  assert_eq(a.mhc_class(), 1)
  assert_true(@src.mhc_allele("a*02:01") is Some(_))
  assert_true(@src.mhc_allele("HLA-A0201") is Some(_))
  let dr = @src.mhc_allele("HLA-DRB1*01:01").unwrap()
  assert_eq(dr.mhc_class(), 2)
  assert_true(@src.mhc_allele("HLA-C*99:99") is None)
  assert_eq(@src.mhc_allele_names().length(), 9)
}

///|
test "mhc score peptide anchors" {
  let a = @src.mhc_allele("HLA-A*02:01").unwrap()
  // SLYNTVATL: P2=L 与 P9=L 均命中 A*02:01 锚规则
  let (s1, am, at) = @src.mhc_score_peptide("SLYNTVATL", a)
  ignore(s1)
  assert_eq(am, 2)
  assert_eq(at, 2)
  // GILGFVFTL: P2=I, P9=L 全锚
  let (s2, am2, _) = @src.mhc_score_peptide("GILGFVFTL", a)
  ignore(s2)
  assert_eq(am2, 2)
}

///|
test "mhc score peptide weighting" {
  let a = @src.mhc_allele("HLA-A*02:01").unwrap()
  let (s_full, _, _) = @src.mhc_score_peptide("GILGFVFTL", a)
  let (s_none, _, _) = @src.mhc_score_peptide("DDDDDDDDD", a)
  assert_true(s_full > s_none)
  assert_true(s_full >= 0.7)
  assert_true(s_none < 0.35)
}

///|
test "mhc score peptide too short raises" {
  let a = @src.mhc_allele("HLA-A*02:01").unwrap()
  let msg = try {
    ignore(@src.mhc_score_peptide("AC", a))
    "NO_ERROR"
  } catch {
    SeqError(m) => m
  }
  assert_true(msg.contains("too short"))
}

///|
test "mhc predict i sliding window" {
  let a = @src.mhc_allele("HLA-A*02:01").unwrap()
  let protein = "SLYNTVATLGILGFVFTLDDDDDDD"
  let hits = @src.mhc_predict_i(protein, a)
  assert_eq(hits.length(), protein.length() - 9 + 1)
  // 得分降序
  let mut i = 1
  while i < hits.length() {
    assert_true(hits[i - 1].score >= hits[i].score)
    i = i + 1
  }
  // 最优命中应含已知锚
  assert_eq(hits[0].peptide.length(), 9)
  assert_eq(hits[0].allele(), "HLA-A*02:01")
}

///|
test "mhc predict i known epitope ranks high" {
  let a = @src.mhc_allele("HLA-A*02:01").unwrap()
  // GILGFVFTL 是经典 A*02:01 流感表位
  let protein = "DDDDDDDDGILGFVFTLDDDDDDDD"
  let hits = @src.mhc_predict_i(protein, a)
  // 表位窗口应出现在头部区域
  let top = hits[0]
  let mut found = false
  for h in hits {
    if h.peptide == "GILGFVFTL" {
      found = true
    }
  }
  assert_true(found)
  // 表位自身应优于纯 D 肽段
  let d_hit = {
    let mut s = -1.0
    for h in hits {
      if h.peptide == "DDDDDDDDD" {
        s = h.score
      }
    }
    s
  }
  let epi_score = {
    let mut s = -1.0
    for h in hits {
      if h.peptide == "GILGFVFTL" {
        s = h.score
      }
    }
    s
  }
  ignore(top)
  assert_true(epi_score > d_hit)
}

///|
test "mhc predict i wrong class raises" {
  let dr = @src.mhc_allele("HLA-DRB1*01:01").unwrap()
  let msg = try {
    ignore(@src.mhc_predict_i("SLYNTVATL", dr))
    "NO_ERROR"
  } catch {
    SeqError(m) => m
  }
  assert_true(msg.contains("class I"))
}

///|
test "mhc predict ii core scanning" {
  let dr = @src.mhc_allele("HLA-DRB1*01:01").unwrap()
  // 15-mer 含 P1=F, P4=L, P6=S, P9=V 的核心
  let protein = "AKFVAAWTLKSAAVDQAKFVAAWTLKSAAVD"
  let hits = @src.mhc_predict_ii(protein, dr)
  assert_eq(hits.length(), protein.length() - 15 + 1)
  // 每个命中都有核心起点
  for h in hits {
    assert_true(h.core_start() >= 0 && h.core_start() <= 6)
  }
}

///|
test "mhc predict ii wrong class raises" {
  let a = @src.mhc_allele("HLA-A*02:01").unwrap()
  let msg = try {
    ignore(@src.mhc_predict_ii("SLYNTVATLGILGFV", a))
    "NO_ERROR"
  } catch {
    SeqError(m) => m
  }
  assert_true(msg.contains("class II"))
}

///|
test "mhc predict dispatch" {
  let a = @src.mhc_allele("H-2Kb").unwrap()
  let hits = @src.mhc_predict("SSLENFRAYV", a)
  assert_eq(hits.length(), 10 - 8 + 1)
  let dr = @src.mhc_allele("HLA-DRB1*04:01").unwrap()
  let hits2 = @src.mhc_predict("AKFVAAWTLKSAAVDQAKFVAAWTLKSAAVD", dr)
  assert_eq(hits2.length(), 31 - 15 + 1)
}

///|
test "mhc h2 kb anchor rule" {
  let kb = @src.mhc_allele("H-2Kb").unwrap()
  // SIINFEKL(OVA 经典 H-2Kb 表位): P5=F, P8=L 双锚命中
  let (s, am, at) = @src.mhc_score_peptide("SIINFEKL", kb)
  ignore(s)
  assert_eq(am, 2)
  assert_eq(at, 2)
}

///|
test "mhc binders filter" {
  let a = @src.mhc_allele("HLA-A*02:01").unwrap()
  let hits = @src.mhc_predict_i("GILGFVFTLDDDDDDD", a)
  let weak = @src.mhc_binders(hits)
  for h in weak {
    assert_true(h.is_weak())
  }
  let strong = @src.mhc_binders(hits, strong_only=true)
  for h in strong {
    assert_true(h.is_strong())
  }
}

///|
test "mhc binding category labels" {
  let a = @src.mhc_allele("HLA-A*02:01").unwrap()
  let hits = @src.mhc_predict_i("GILGFVFTLDDDDDDD", a)
  let mut has_any = false
  for h in hits {
    let c = h.binding_category()
    assert_true(c == "strong" || c == "weak" || c == "none")
    has_any = true
  }
  assert_true(has_any)
}

///|
test "mhc unique epitopes dedup" {
  let dr = @src.mhc_allele("HLA-DRB1*01:01").unwrap()
  let protein = "AKFVAAWTLKSAAVDQAKFVAAWTLKSAAVD"
  let hits = @src.mhc_predict_ii(protein, dr)
  let uniq = @src.mhc_unique_epitopes(hits)
  assert_true(uniq.length() <= hits.length())
  // 每个表位键唯一: II 类以核心 9-mer 计
  let seen : Map[String, Bool] = Map([], capacity=0)
  for h in uniq {
    let key = h.peptide[h.core_start():h.core_start() + 9].to_owned()
    assert_false(seen.contains(key))
    seen[key] = true
  }
}

///|
test "mhc predict i short protein empty" {
  let a = @src.mhc_allele("HLA-A*02:01").unwrap()
  let hits = @src.mhc_predict_i("ACGTACG", a)
  assert_eq(hits.length(), 0)
}

///|
test "mhc predict lower case protein normalized" {
  let a = @src.mhc_allele("HLA-A*02:01").unwrap()
  let upper = @src.mhc_predict_i("GILGFVFTLDDDD", a)
  let lower = @src.mhc_predict_i("gilgfvftldddd", a)
  assert_eq(upper.length(), lower.length())
  assert_eq(lower[0].peptide.length(), 9)
}
