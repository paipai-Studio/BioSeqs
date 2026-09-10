///|
/// Tests for the multtest port (Bioconductor multtest::mt.rawp2adjp,
/// mt.reject, and the Westfall-Young maxT/minP permutation engines).
///
/// Raw-p adjustment vectors, sort indices and adaptive h0 estimates are
/// asserted against the reference R implementation (multtest 2.8.0,
/// R/mt.basic.R sourced directly); maxT/minP vectors against a numpy mirror
/// of the same step-down algorithms.
fn mt_vec_eq(got : Array[Double], want : Array[Double]) -> Bool {
  if got.length() != want.length() {
    return false
  }
  let mut i = 0
  while i < got.length() {
    if want[i].is_nan() {
      if !got[i].is_nan() {
        return false
      }
    } else if (got[i] - want[i]).abs() > 1.0e-12 {
      return false
    }
    i = i + 1
  }
  true
}

///|
test "mt_rawp2adjp grid: nine procedures against reference R" {
  let rawp = [0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.4, 0.7, 0.99]
  let r = @src.mt_rawp2adjp(rawp~)
  assert_eq(r.index, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
  assert_true(r.h0_abh == 5.0)
  assert_eq(r.h0_tsbh, [(0.05, 10.0)])
  assert_true(mt_vec_eq(r.rawp, rawp))
  assert_true(
    mt_vec_eq(r.adjp.get("Bonferroni").unwrap(), [
      0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("Holm").unwrap(), [
      0.1, 0.18, 0.24, 0.28, 0.3, 0.5, 0.8, 1.0, 1.0, 1.0,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("Hochberg").unwrap(), [
      0.1, 0.18, 0.24, 0.28, 0.3, 0.5, 0.8, 0.99, 0.99, 0.99,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("SidakSS").unwrap(), [
      0.0956179249912, 0.182927193112, 0.262575873105, 0.335167364008, 0.401263060762,
      0.6513215599, 0.8926258176, 0.9939533824, 0.9999940951, 1.0,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("SidakSD").unwrap(), [
      0.0956179249912, 0.16625223787, 0.216256640562, 0.248552521892, 0.264908109375,
      0.40951, 0.5904, 0.784, 0.91, 0.99,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("BH").unwrap(), [
      0.1, 0.1, 0.1, 0.1, 0.1, 0.166666666667, 0.285714285714, 0.5, 0.777777777778,
      0.99,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("BY").unwrap(), [
      0.292896825397, 0.292896825397, 0.292896825397, 0.292896825397, 0.292896825397,
      0.488161375661, 0.836848072562, 1.0, 1.0, 1.0,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("ABH").unwrap(), [
      0.05, 0.05, 0.05, 0.05, 0.05, 0.0833333333333, 0.142857142857, 0.25, 0.388888888889,
      0.495,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("TSBH_0.05").unwrap(), [
      0.1, 0.1, 0.1, 0.1, 0.1, 0.166666666667, 0.285714285714, 0.5, 0.777777777778,
      0.99,
    ]),
  )
  // mt.reject counts at alpha = 0.05 for columns in reference order
  let names = [
    "rawp", "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH",
    "TSBH_0.05",
  ]
  let counts = [5, 0, 0, 0, 0, 0, 0, 0, 5, 0]
  for i = 0; i < names.length(); i = i + 1 {
    let (rr, _) = @src.mt_reject(adjp=r.adjp.get(names[i]).unwrap())
    assert_eq(rr, counts[i])
  }
}

///|
test "mt_rawp2adjp ties share the same adjusted values" {
  let rawp = [0.01, 0.01, 0.05, 0.05, 0.05, 0.2, 0.2, 0.8]
  let r = @src.mt_rawp2adjp(rawp~)
  assert_true(r.h0_abh == 3.0)
  assert_eq(r.h0_tsbh, [(0.05, 6.0)])
  assert_true(
    mt_vec_eq(r.adjp.get("Bonferroni").unwrap(), [
      0.08, 0.08, 0.4, 0.4, 0.4, 1.0, 1.0, 1.0,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("Holm").unwrap(), [
      0.08, 0.08, 0.3, 0.3, 0.3, 0.6, 0.6, 0.8,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("Hochberg").unwrap(), [
      0.07, 0.07, 0.2, 0.2, 0.2, 0.4, 0.4, 0.8,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("SidakSS").unwrap(), [
      0.0772553055721, 0.0772553055721, 0.336579568711, 0.336579568711, 0.336579568711,
      0.83222784, 0.83222784, 0.99999744,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("SidakSD").unwrap(), [
      0.0772553055721, 0.0772553055721, 0.264908109375, 0.264908109375, 0.264908109375,
      0.488, 0.488, 0.8,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("BH").unwrap(), [
      0.04, 0.04, 0.08, 0.08, 0.08, 0.228571428571, 0.228571428571, 0.8,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("BY").unwrap(), [
      0.108714285714, 0.108714285714, 0.217428571429, 0.217428571429, 0.217428571429,
      0.621224489796, 0.621224489796, 1.0,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("ABH").unwrap(), [
      0.015, 0.015, 0.03, 0.03, 0.03, 0.0857142857143, 0.0857142857143, 0.3,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("TSBH_0.05").unwrap(), [
      0.03, 0.03, 0.06, 0.06, 0.06, 0.171428571429, 0.171428571429, 0.6,
    ]),
  )
  let names = [
    "rawp", "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH",
    "TSBH_0.05",
  ]
  let counts = [5, 0, 0, 0, 0, 0, 2, 0, 5, 2]
  for i = 0; i < names.length(); i = i + 1 {
    let (rr, _) = @src.mt_reject(adjp=r.adjp.get(names[i]).unwrap())
    assert_eq(rr, counts[i])
  }
}

///|
test "mt_rawp2adjp boundary p-values 0 and 1" {
  let rawp = [0.0, 0.0, 1.0e-6, 0.5, 1.0, 1.0, 1.0]
  let r = @src.mt_rawp2adjp(rawp~)
  assert_true(r.h0_abh == 6.0)
  assert_eq(r.h0_tsbh, [(0.05, 4.0)])
  assert_true(
    mt_vec_eq(r.adjp.get("Bonferroni").unwrap(), [
      0.0, 0.0, 7.0e-6, 1.0, 1.0, 1.0, 1.0,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("Holm").unwrap(), [
      0.0, 0.0, 5.0e-6, 1.0, 1.0, 1.0, 1.0,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("Hochberg").unwrap(), [
      0.0, 0.0, 5.0e-6, 1.0, 1.0, 1.0, 1.0,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("SidakSS").unwrap(), [
      0.0, 0.0, 6.99997900022e-6, 0.9921875, 1.0, 1.0, 1.0,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("SidakSD").unwrap(), [
      0.0, 0.0, 4.99999000014e-6, 0.9375, 1.0, 1.0, 1.0,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("BH").unwrap(), [
      0.0, 0.0, 2.33333333333e-6, 0.875, 1.0, 1.0, 1.0,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("BY").unwrap(), [0.0, 0.0, 6.05e-6, 1.0, 1.0, 1.0, 1.0]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("ABH").unwrap(), [
      0.0, 0.0, 2.0e-6, 0.75, 0.857142857143, 0.857142857143, 0.857142857143,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("TSBH_0.05").unwrap(), [
      0.0, 0.0, 1.33333333333e-6, 0.5, 0.571428571429, 0.571428571429, 0.571428571429,
    ]),
  )
  // every column rejects exactly the three smallest p-values
  for name in @src.mt_procs {
    let key = if name == "TSBH" { "TSBH_0.05" } else { name }
    let (rr, _) = @src.mt_reject(adjp=r.adjp.get(key).unwrap())
    assert_eq(rr, 3)
  }
  let (rr, _) = @src.mt_reject(adjp=r.rawp)
  assert_eq(rr, 3)
}

///|
test "mt_rawp2adjp sort order and stable tie placement" {
  let rawp = [0.5, 0.001, 0.02, 0.03, 0.04, 0.5, 0.6, 0.7]
  let r = @src.mt_rawp2adjp(rawp~)
  assert_eq(r.index, [1, 2, 3, 4, 0, 5, 6, 7])
  assert_true(r.h0_abh == 6.0)
  assert_eq(r.h0_tsbh, [(0.05, 7.0)])
  assert_true(
    mt_vec_eq(r.adjp.get("Bonferroni").unwrap(), [
      0.008, 0.16, 0.24, 0.32, 1.0, 1.0, 1.0, 1.0,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("Hochberg").unwrap(), [
      0.008, 0.14, 0.18, 0.2, 0.7, 0.7, 0.7, 0.7,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("SidakSS").unwrap(), [
      0.00797205593006, 0.149236977418, 0.216256640562, 0.278610421016, 0.99609375,
      0.99609375, 0.99934464, 0.99993439,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("SidakSD").unwrap(), [
      0.00797205593006, 0.131874466753, 0.167027995071, 0.1846273024, 0.9375, 0.9375,
      0.9375, 0.9375,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("BH").unwrap(), [
      0.008, 0.08, 0.08, 0.08, 0.666666666667, 0.666666666667, 0.685714285714, 0.7,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("BY").unwrap(), [
      0.0217428571429, 0.217428571429, 0.217428571429, 0.217428571429, 1.0, 1.0,
      1.0, 1.0,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("ABH").unwrap(), [
      0.006, 0.06, 0.06, 0.06, 0.5, 0.5, 0.514285714286, 0.525,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("TSBH_0.05").unwrap(), [
      0.007, 0.07, 0.07, 0.07, 0.583333333333, 0.583333333333, 0.6, 0.6125,
    ]),
  )
}

///|
test "mt_rawp2adjp NaN values are sorted last and propagated" {
  let nan = 0.0 / 0.0
  let rawp = [0.01, nan, 0.05, nan, 0.5]
  // na.rm = FALSE: mgood stays m = 5
  let r = @src.mt_rawp2adjp(rawp~)
  assert_eq(r.index, [0, 2, 4, 1, 3])
  assert_true(r.h0_abh == 5.0)
  assert_eq(r.h0_tsbh, [(0.05, 5.0)])
  assert_true(mt_vec_eq(r.rawp, [0.01, 0.05, 0.5, nan, nan]))
  assert_true(
    mt_vec_eq(r.adjp.get("Bonferroni").unwrap(), [0.05, 0.25, 1.0, nan, nan]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("SidakSS").unwrap(), [
      0.0490099501, 0.2262190625, 0.96875, nan, nan,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("SidakSD").unwrap(), [
      0.0490099501, 0.18549375, 0.875, nan, nan,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("BH").unwrap(), [0.05, 0.125, 0.833333333333, nan, nan]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("BY").unwrap(), [
      0.114166666667, 0.285416666667, 1.0, nan, nan,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("ABH").unwrap(), [
      0.05, 0.125, 0.833333333333, nan, nan,
    ]),
  )
  assert_true(
    mt_vec_eq(r.adjp.get("TSBH_0.05").unwrap(), [
      0.05, 0.125, 0.833333333333, nan, nan,
    ]),
  )
  // R returns NA for the rejection total with NaN comparisons
  let (rr, which) = @src.mt_reject(adjp=r.adjp.get("BH").unwrap())
  assert_eq(rr, 0)
  assert_eq(which, [Some(true), Some(false), Some(false), None, None])

  // na.rm = TRUE: mgood = 3 changes every multiplier
  let r2 = @src.mt_rawp2adjp(rawp~, na_rm=true)
  assert_eq(r2.index, [0, 2, 4, 1, 3])
  assert_true(r2.h0_abh.is_nan())
  assert_eq(r2.h0_tsbh, [(0.05, 2.0)])
  assert_true(
    mt_vec_eq(r2.adjp.get("Bonferroni").unwrap(), [0.03, 0.15, 1.0, nan, nan]),
  )
  assert_true(
    mt_vec_eq(r2.adjp.get("Holm").unwrap(), [0.03, 0.1, 0.5, nan, nan]),
  )
  assert_true(
    mt_vec_eq(r2.adjp.get("SidakSS").unwrap(), [
      0.029701, 0.142625, 0.875, nan, nan,
    ]),
  )
  assert_true(
    mt_vec_eq(r2.adjp.get("SidakSD").unwrap(), [0.029701, 0.0975, 0.5, 0.5, nan]),
  )
  assert_true(
    mt_vec_eq(r2.adjp.get("BH").unwrap(), [0.03, 0.075, 0.5, nan, nan]),
  )
  assert_true(
    mt_vec_eq(r2.adjp.get("BY").unwrap(), [
      0.055, 0.1375, 0.916666666667, nan, nan,
    ]),
  )
  assert_true(mt_vec_eq(r2.adjp.get("ABH").unwrap(), [nan, nan, nan, nan, nan]))
  assert_true(
    mt_vec_eq(r2.adjp.get("TSBH_0.05").unwrap(), [
      0.02, 0.05, 0.333333333333, nan, nan,
    ]),
  )
}

///|
test "mt_rawp2adjp all-NaN with na_rm and mgood zero" {
  let nan = 0.0 / 0.0
  let r = @src.mt_rawp2adjp(rawp=[nan, nan], na_rm=true)
  assert_true(r.h0_abh.is_nan())
  assert_eq(r.h0_tsbh, [(0.05, 0.0)])
  assert_true(mt_vec_eq(r.adjp.get("Bonferroni").unwrap(), [nan, nan]))
  assert_true(mt_vec_eq(r.adjp.get("Holm").unwrap(), [nan, nan]))
  assert_true(mt_vec_eq(r.adjp.get("Hochberg").unwrap(), [nan, nan]))
  // 1 - (1 - NA)^0 == 0 in R
  assert_true(mt_vec_eq(r.adjp.get("SidakSS").unwrap(), [0.0, 0.0]))
  assert_true(mt_vec_eq(r.adjp.get("SidakSD").unwrap(), [0.0, nan]))
  assert_true(mt_vec_eq(r.adjp.get("BH").unwrap(), [nan, nan]))
  assert_true(mt_vec_eq(r.adjp.get("BY").unwrap(), [nan, nan]))
  assert_true(mt_vec_eq(r.adjp.get("ABH").unwrap(), [nan, nan]))
  assert_true(mt_vec_eq(r.adjp.get("TSBH_0.05").unwrap(), [nan, nan]))
}

///|
test "mt_rawp2adjp multiple nominal levels for TSBH" {
  let rawp = [0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.4, 0.7, 0.99]
  let r = @src.mt_rawp2adjp(rawp~, procs=["BH", "TSBH"], alpha=[0.05, 0.1, 0.2])
  assert_eq(r.h0_tsbh, [(0.05, 10.0), (0.1, 10.0), (0.2, 5.0)])
  // only rawp, BH and the three TSBH columns exist
  assert_eq(r.adjp.length(), 5)
  let bh = [
    0.1, 0.1, 0.1, 0.1, 0.1, 0.166666666667, 0.285714285714, 0.5, 0.777777777778,
    0.99,
  ]
  assert_true(mt_vec_eq(r.adjp.get("BH").unwrap(), bh))
  assert_true(mt_vec_eq(r.adjp.get("TSBH_0.05").unwrap(), bh))
  assert_true(mt_vec_eq(r.adjp.get("TSBH_0.1").unwrap(), bh))
  assert_true(
    mt_vec_eq(r.adjp.get("TSBH_0.2").unwrap(), [
      0.05, 0.05, 0.05, 0.05, 0.05, 0.0833333333333, 0.142857142857, 0.25, 0.388888888889,
      0.495,
    ]),
  )
}

///|
test "mt_rawp2adjp input validation" {
  let raised_empty = try {
    let _ = @src.mt_rawp2adjp(rawp=[])
    false
  } catch {
    _ => true
  }
  assert_true(raised_empty)
  let raised_range = try {
    let _ = @src.mt_rawp2adjp(rawp=[0.01, 1.2])
    false
  } catch {
    _ => true
  }
  assert_true(raised_range)
  let raised_neg = try {
    let _ = @src.mt_rawp2adjp(rawp=[0.01, -0.5])
    false
  } catch {
    _ => true
  }
  assert_true(raised_neg)
  // unknown procedure names are silently ignored, as in R
  let r = @src.mt_rawp2adjp(rawp=[0.01, 0.04], procs=["BH", "nope"])
  assert_eq(r.adjp.length(), 2)
  assert_true(mt_vec_eq(r.adjp.get("BH").unwrap(), [0.02, 0.04]))
}

///|
test "mt_reject vector flags and alpha argument" {
  let adjp = [0.01, 0.04, 0.049, 0.05, 0.2]
  let (r05, w05) = @src.mt_reject(adjp~)
  assert_eq(r05, 4)
  assert_eq(w05, [Some(true), Some(true), Some(true), Some(true), Some(false)])
  let (r10, _) = @src.mt_reject(adjp~, alpha=0.1)
  assert_eq(r10, 4)
  let (r20, _) = @src.mt_reject(adjp~, alpha=0.2)
  assert_eq(r20, 5)
  let raised = try {
    let _ = @src.mt_reject(adjp=[])
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "mt_max_t case B ties and small B" {
  let v = [
    [0.5, -0.5, 0.0, 1.0, -1.0, 0.2, -0.2, 0.9, -0.9, 0.1],
    [0.0, 0.0, 0.0, 2.0, -2.0, 0.1, -0.1, 0.3, -0.3, 0.2],
    [1.5, -1.5, 0.0, 0.0, 0.0, 1.2, -1.2, 0.8, -0.8, 0.4],
  ]
  let adj = @src.mt_max_t(observed=[0.5, 2.0, 1.5], null_stats=v)
  assert_true(mt_vec_eq(adj, [0.8, 0.2, 0.4]))
}

///|
test "mt_max_t single hypothesis" {
  let v = [[0.1, -0.2, 0.3, -0.4, 0.5]]
  let adj = @src.mt_max_t(observed=[0.25], null_stats=v)
  assert_true(mt_vec_eq(adj, [0.6]))
}

///|
test "mt_max_t multi-gene step-down case D" {
  let v = [
    [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
    [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2],
    [0.05, 0.15, 0.25, 0.95, 0.05, 0.15, 0.25, 0.35],
    [0.2, 0.2, 0.2, 0.2, 0.2, 0.9, 0.2, 0.2],
  ]
  let adj = @src.mt_max_t(observed=[0.75, 0.85, 0.9, 0.88], null_stats=v)
  assert_true(mt_vec_eq(adj, [0.625, 0.375, 0.375, 0.375]))
}

///|
test "mt_min_p case B ties and small B" {
  let v = [
    [0.5, -0.5, 0.0, 1.0, -1.0, 0.2, -0.2, 0.9, -0.9, 0.1],
    [0.0, 0.0, 0.0, 2.0, -2.0, 0.1, -0.1, 0.3, -0.3, 0.2],
    [1.5, -1.5, 0.0, 0.0, 0.0, 1.2, -1.2, 0.8, -0.8, 0.4],
  ]
  let adj = @src.mt_min_p(observed=[0.5, 2.0, 1.5], null_stats=v)
  assert_true(mt_vec_eq(adj, [0.9, 0.4, 0.4]))
}

///|
test "mt_min_p single hypothesis" {
  let v = [[0.1, -0.2, 0.3, -0.4, 0.5]]
  let adj = @src.mt_min_p(observed=[0.25], null_stats=v)
  assert_true(mt_vec_eq(adj, [0.6]))
}

///|
test "mt_min_p multi-gene step-down case D" {
  let v = [
    [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
    [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2],
    [0.05, 0.15, 0.25, 0.95, 0.05, 0.15, 0.25, 0.35],
    [0.2, 0.2, 0.2, 0.2, 0.2, 0.9, 0.2, 0.2],
  ]
  let adj = @src.mt_min_p(observed=[0.75, 0.85, 0.9, 0.88], null_stats=v)
  assert_true(mt_vec_eq(adj, [0.5, 0.5, 0.5, 0.5]))
}

///|
test "mt_max_t and mt_min_p matrix validation" {
  let bad_rows = try {
    let _ = @src.mt_min_p(observed=[0.1, 0.2], null_stats=[[0.1, 0.2]])
    false
  } catch {
    _ => true
  }
  assert_true(bad_rows)
  let ragged = try {
    let _ = @src.mt_max_t(observed=[0.1, 0.2], null_stats=[[0.1, 0.2], [0.3]])
    false
  } catch {
    _ => true
  }
  assert_true(ragged)
  let no_perms = try {
    let _ = @src.mt_max_t(observed=[0.1], null_stats=[[]])
    false
  } catch {
    _ => true
  }
  assert_true(no_perms)
  let empty_obs = try {
    let _ = @src.mt_max_t(observed=[], null_stats=[])
    false
  } catch {
    _ => true
  }
  assert_true(empty_obs)
}
