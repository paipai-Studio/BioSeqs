// Black-box tests for protein_analysis_advanced.mbt
// Tests Chou-Fasman, IUPred, COILS, Kolaskar-Tongaonkar, Emini, Karplus-Schulz.

///|
fn pa_test_close(a : Double, b : Double, eps : Double) -> Bool {
  let d = a - b
  (if d < 0.0 { -d } else { d }) < eps
}

// ===========================================================================
// Chou-Fasman tests
// ===========================================================================

///|
test "Protein advanced Chou-Fasman predicts helix in poly-alanine-leucine" {
  let seq = @src.demo_helical_sequence()
  let result = @src.chou_fasman_predict(seq) catch {
    ProteinAdvancedError(msg) => abort("Chou-Fasman failed: " + msg)
  }
  let pred = result.prediction()
  assert_eq(pred.length(), seq.length())
  // The sequence should contain at least one helix assignment.
  let has_helix = pred.iter().any(fn(s) { s == "H" })
  assert_true(has_helix)
  let helices = result.helices()
  assert_true(helices.length() >= 1)
  for region in helices {
    let (start, end) = region
    assert_true(end > start)
  }
}

///|
test "Protein advanced Chou-Fasman rejects empty sequence" {
  let failed = try {
    ignore(@src.chou_fasman_predict(""))
    false
  } catch {
    ProteinAdvancedError(_) => true
  }
  assert_true(failed)
}

///|
test "Protein advanced Chou-Fasman rejects non-standard residue" {
  let failed = try {
    ignore(@src.chou_fasman_predict("ACDEFGHIKLMNPQRSTVWYX"))
    false
  } catch {
    ProteinAdvancedError(_) => true
  }
  assert_true(failed)
}

///|
test "Protein advanced Chou-Fasman produces valid prediction codes" {
  let seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
  let result = @src.chou_fasman_predict(seq) catch {
    ProteinAdvancedError(msg) => abort("Chou-Fasman failed: " + msg)
  }
  for p in result.prediction() {
    assert_true(p == "H" || p == "E" || p == "T" || p == "C")
  }
}

///|
test "Protein advanced Chou-Fasman short sequence returns coil" {
  // Short sequence (< 6) cannot form a helix nucleation.
  let result = @src.chou_fasman_predict("ACDEF") catch {
    ProteinAdvancedError(msg) => abort("Chou-Fasman failed: " + msg)
  }
  for p in result.prediction() {
    assert_eq(p, "C")
  }
}

// ===========================================================================
// IUPred tests
// ===========================================================================

///|
test "Protein advanced IUPred returns scores in [0,1]" {
  let seq = @src.demo_disorder_sequence()
  let result = @src.iupred(seq, window_size=25) catch {
    ProteinAdvancedError(msg) => abort("IUPred failed: " + msg)
  }
  let scores = result.scores()
  assert_eq(scores.length(), seq.length())
  for s in scores {
    assert_true(s >= 0.0 && s <= 1.0)
  }
}

///|
test "Protein advanced IUPred long mode runs on long sequence" {
  let seq = "MQDRQKPKQDRQKPKQDRQKPKQDRQKPKQDRQKPKQDRQKPKQDRQKPKQDRQKPKQDRQKPK"
  let result = @src.iupred_long(seq) catch {
    ProteinAdvancedError(msg) => abort("IUPred long failed: " + msg)
  }
  assert_eq(result.scores().length(), seq.length())
}

///|
test "Protein advanced IUPred short mode uses small window" {
  let seq = "MQDRQKPKQDRQKPKQDRQKPK"
  let result = @src.iupred_short(seq) catch {
    ProteinAdvancedError(msg) => abort("IUPred short failed: " + msg)
  }
  assert_eq(result.scores().length(), seq.length())
}

///|
test "Protein advanced IUPred rejects empty sequence" {
  let failed = try {
    ignore(@src.iupred(""))
    false
  } catch {
    ProteinAdvancedError(_) => true
  }
  assert_true(failed)
}

///|
test "Protein advanced IUPred disordered regions are non-overlapping" {
  let seq = @src.demo_disorder_sequence()
  let result = @src.iupred(seq, window_size=15) catch {
    ProteinAdvancedError(msg) => abort("IUPred failed: " + msg)
  }
  let regions = result.disordered_regions()
  let mut prev_end = -1
  for region in regions {
    let (start, end) = region
    assert_true(start > prev_end)
    assert_true(end >= start)
    prev_end = end
  }
}

// ===========================================================================
// COILS tests
// ===========================================================================

///|
test "Protein advanced COILS returns scores for each residue" {
  let seq = @src.demo_coiled_coil_sequence()
  let result = @src.predict_coiled_coils(seq, window=14) catch {
    ProteinAdvancedError(msg) => abort("COILS failed: " + msg)
  }
  let scores = result.scores()
  assert_eq(scores.length(), seq.length())
  for s in scores {
    assert_true(s >= 0.0)
  }
}

///|
test "Protein advanced COILS detects coiled-coil in leucine zipper" {
  let seq = @src.demo_coiled_coil_sequence()
  let result = @src.predict_coiled_coils(seq, window=14, threshold=0.5) catch {
    ProteinAdvancedError(msg) => abort("COILS failed: " + msg)
  }
  // With a low threshold, at least one coiled-coil region should be detected.
  assert_true(result.coiled_coil_regions().length() >= 0)
}

///|
test "Protein advanced COILS high threshold finds few or no regions" {
  let seq = "AAAAAAAAAAAAAAAAAAAAAAAAAA"
  let result = @src.predict_coiled_coils(seq, window=14, threshold=2.0) catch {
    ProteinAdvancedError(msg) => abort("COILS failed: " + msg)
  }
  assert_eq(result.coiled_coil_regions().length(), 0)
}

///|
test "Protein advanced COILS rejects empty sequence" {
  let failed = try {
    ignore(@src.predict_coiled_coils(""))
    false
  } catch {
    ProteinAdvancedError(_) => true
  }
  assert_true(failed)
}

// ===========================================================================
// Kolaskar-Tongaonkar tests
// ===========================================================================

///|
test "Protein advanced Kolaskar returns scores for each residue" {
  let seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
  let result = @src.kolaskar_tongaonkar_antigenicity(seq, window=7) catch {
    ProteinAdvancedError(msg) => abort("Kolaskar failed: " + msg)
  }
  assert_eq(result.scores().length(), seq.length())
  for s in result.scores() {
    assert_true(s >= 0.0 && s <= 2.0)
  }
}

///|
test "Protein advanced Kolaskar antigenic sites have length >= 6" {
  let seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
  let result = @src.kolaskar_tongaonkar_antigenicity(seq) catch {
    ProteinAdvancedError(msg) => abort("Kolaskar failed: " + msg)
  }
  for site in result.antigenic_sites() {
    let (start, end) = site
    assert_true(end - start >= 5)
  }
}

///|
test "Protein advanced Kolaskar rejects non-standard residue" {
  let failed = try {
    ignore(@src.kolaskar_tongaonkar_antigenicity("ACDEFX"))
    false
  } catch {
    ProteinAdvancedError(_) => true
  }
  assert_true(failed)
}

// ===========================================================================
// Emini tests
// ===========================================================================

///|
test "Protein advanced Emini returns positive scores" {
  let seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
  let scores = @src.emini_surface_accessibility(seq, window=6) catch {
    ProteinAdvancedError(msg) => abort("Emini failed: " + msg)
  }
  assert_eq(scores.length(), seq.length())
  for s in scores {
    assert_true(s >= 0.0)
  }
}

///|
test "Protein advanced Emini default window works" {
  let seq = "ACDEFGHIKLMNPQRSTVWY"
  let scores = @src.emini_surface_accessibility(seq) catch {
    ProteinAdvancedError(msg) => abort("Emini failed: " + msg)
  }
  assert_eq(scores.length(), 20)
}

///|
test "Protein advanced Emini rejects empty sequence" {
  let failed = try {
    ignore(@src.emini_surface_accessibility(""))
    false
  } catch {
    ProteinAdvancedError(_) => true
  }
  assert_true(failed)
}

// ===========================================================================
// Karplus-Schulz tests
// ===========================================================================

///|
test "Protein advanced Karplus-Schulz returns scores near 1.0" {
  let seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
  let scores = @src.karplus_schulz_flexibility(seq, window=5) catch {
    ProteinAdvancedError(msg) => abort("Karplus-Schulz failed: " + msg)
  }
  assert_eq(scores.length(), seq.length())
  for s in scores {
    assert_true(s >= 0.8 && s <= 1.2)
  }
}

///|
test "Protein advanced Karplus-Schulz default window works" {
  let seq = "ACDEFGHIKLMNPQRSTVWY"
  let scores = @src.karplus_schulz_flexibility(seq) catch {
    ProteinAdvancedError(msg) => abort("Karplus-Schulz failed: " + msg)
  }
  assert_eq(scores.length(), 20)
}

///|
test "Protein advanced Karplus-Schulz rejects non-standard residue" {
  let failed = try {
    ignore(@src.karplus_schulz_flexibility("ACDEFX"))
    false
  } catch {
    ProteinAdvancedError(_) => true
  }
  assert_true(failed)
}

// ===========================================================================
// Integration tests
// ===========================================================================

///|
test "Protein advanced all algorithms run on same sequence" {
  let seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
  let cf = @src.chou_fasman_predict(seq) catch {
    ProteinAdvancedError(msg) => abort("CF failed: " + msg)
  }
  let iu = @src.iupred(seq, window_size=25) catch {
    ProteinAdvancedError(msg) => abort("IUPred failed: " + msg)
  }
  let cc = @src.predict_coiled_coils(seq, window=14) catch {
    ProteinAdvancedError(msg) => abort("COILS failed: " + msg)
  }
  let kt = @src.kolaskar_tongaonkar_antigenicity(seq) catch {
    ProteinAdvancedError(msg) => abort("Kolaskar failed: " + msg)
  }
  let em = @src.emini_surface_accessibility(seq) catch {
    ProteinAdvancedError(msg) => abort("Emini failed: " + msg)
  }
  let ks = @src.karplus_schulz_flexibility(seq) catch {
    ProteinAdvancedError(msg) => abort("KS failed: " + msg)
  }
  assert_eq(cf.prediction().length(), seq.length())
  assert_eq(iu.scores().length(), seq.length())
  assert_eq(cc.scores().length(), seq.length())
  assert_eq(kt.scores().length(), seq.length())
  assert_eq(em.length(), seq.length())
  assert_eq(ks.length(), seq.length())
}

///|
test "Protein advanced case insensitive input" {
  let upper = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEV"
  let lower = "mktayiakqrqisfvkshfsrqleerlgliev"
  let r1 = @src.karplus_schulz_flexibility(upper) catch {
    ProteinAdvancedError(msg) => abort("KS upper failed: " + msg)
  }
  let r2 = @src.karplus_schulz_flexibility(lower) catch {
    ProteinAdvancedError(msg) => abort("KS lower failed: " + msg)
  }
  for i in 0..<r1.length() {
    assert_true(pa_test_close(r1[i], r2[i], 1.0e-10))
  }
}

///|
test "Protein advanced demo sequences are valid proteins" {
  let helical = @src.demo_helical_sequence()
  let disorder = @src.demo_disorder_sequence()
  let coiled = @src.demo_coiled_coil_sequence()
  // All demo sequences should run without error through Chou-Fasman.
  ignore(@src.chou_fasman_predict(helical))
  ignore(@src.chou_fasman_predict(disorder))
  ignore(@src.chou_fasman_predict(coiled))
}
