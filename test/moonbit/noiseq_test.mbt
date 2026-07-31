///|
/// Test file for NOISeq module.

test "noiseq_sample_creation" {
  let s = @src.NOISeqSample::new("sample1", "control")
  assert_eq(s.sample_id, "sample1")
  assert_eq(s.condition, "control")
}

test "noiseq_sample_set_get_count" {
  let s = @src.NOISeqSample::new("sample1", "control")
  s.set_count("Gene1", 100.0)
  s.set_count("Gene2", 50.0)
  assert_eq(s.get_count("Gene1"), 100.0)
  assert_eq(s.get_count("Gene2"), 50.0)
  assert_eq(s.get_count("Gene3"), 0.0) // Non-existent
}

test "noiseq_sample_library_size" {
  let s = @src.NOISeqSample::new("sample1", "control")
  s.set_count("Gene1", 100.0)
  s.set_count("Gene2", 200.0)
  s.set_count("Gene3", 50.0)
  assert_eq(s.library_size(), 350.0)
}

test "noiseq_sample_n_expressed" {
  let s = @src.NOISeqSample::new("sample1", "control")
  s.set_count("Gene1", 100.0)
  s.set_count("Gene2", 0.0)
  s.set_count("Gene3", 50.0)
  assert_eq(s.n_expressed(), 2)
}

test "noiseq_normalize_tmm" {
  let samples : Array[@src.NOISeqSample] = Array::new()
  let s1 = @src.NOISeqSample::new("s1", "ctrl")
  s1.set_count("Gene1", 100.0)
  s1.set_count("Gene2", 200.0)
  samples.push(s1)
  let s2 = @src.NOISeqSample::new("s2", "ctrl")
  s2.set_count("Gene1", 200.0)
  s2.set_count("Gene2", 400.0)
  samples.push(s2)

  let norm = @src.noiseq_normalize(samples, norm_method=@src.noiseq_norm_tmm())
  // TMM normalizes to ref sample library size
  // s1 lib=300, s2 lib=600, ref=s1 lib=300
  // s2 scale = 300/600 = 0.5, so Gene1=100, Gene2=200
  assert_true((norm[1].get_count("Gene1") - 100.0).abs() < 0.01)
}

test "noiseq_normalize_rpkm" {
  let samples : Array[@src.NOISeqSample] = Array::new()
  let s1 = @src.NOISeqSample::new("s1", "ctrl")
  s1.set_count("Gene1", 500000.0)
  samples.push(s1)

  let norm = @src.noiseq_normalize(samples, norm_method=@src.noiseq_norm_rpkm())
  // CPM: 500000 * 1e6 / 500000 = 1e6
  assert_true((norm[0].get_count("Gene1") - 1000000.0).abs() < 0.01)
}

test "noiseq_normalize_none" {
  let samples : Array[@src.NOISeqSample] = Array::new()
  let s1 = @src.NOISeqSample::new("s1", "ctrl")
  s1.set_count("Gene1", 100.0)
  samples.push(s1)

  let norm = @src.noiseq_normalize(samples, norm_method=@src.noiseq_norm_none())
  assert_eq(norm[0].get_count("Gene1"), 100.0)
}

test "noiseq_run_basic" {
  let (ctrl, trt) = @src.noiseq_sample_data()
  let results = @src.noiseq_run(ctrl, trt)
  assert_eq(results.get_n_results(), 10) // 10 genes
  assert_eq(results.n_genes, 10)
}

test "noiseq_run_significant" {
  let (ctrl, trt) = @src.noiseq_sample_data()
  let results = @src.noiseq_run(ctrl, trt, prob_threshold=0.3)
  let sig = results.get_significant()
  assert_true(sig.length() > 0)
}

test "noiseq_run_top_genes" {
  let (ctrl, trt) = @src.noiseq_sample_data()
  let results = @src.noiseq_run(ctrl, trt)
  let top = results.get_top_genes(3)
  assert_eq(top.length(), 3)
  // Results should be sorted by probability (descending)
  assert_true(top[0].prob >= top[1].prob)
}

test "noiseq_run_up_down_regulated" {
  let (ctrl, trt) = @src.noiseq_sample_data()
  let results = @src.noiseq_run(ctrl, trt, prob_threshold=0.3)
  let up = results.get_up_regulated()
  let down = results.get_down_regulated()
  // We expect some up-regulated (first 3 genes) and down-regulated (next 3)
  assert_true(up.length() + down.length() > 0)
}

test "noiseq_run_summary" {
  let (ctrl, trt) = @src.noiseq_sample_data()
  let results = @src.noiseq_run(ctrl, trt)
  let s = results.summary()
  assert_true(s.contains("NOISeq"))
  assert_true(s.contains("Genes tested"))
}

test "noiseq_qc" {
  let (ctrl, trt) = @src.noiseq_sample_data()
  let all : Array[@src.NOISeqSample] = Array::new()
  let mut i = 0
  while i < ctrl.length() {
    all.push(ctrl[i])
    i = i + 1
  }
  let mut j = 0
  while j < trt.length() {
    all.push(trt[j])
    j = j + 1
  }
  let qc = @src.noiseq_qc(all)
  assert_true(qc.biotype_counts.size() > 0)
}

test "noiseq_norm_to_string" {
  assert_eq(@src.noiseq_norm_rpkm().to_string(), "RPKM")
  assert_eq(@src.noiseq_norm_tmm().to_string(), "TMM")
  assert_eq(@src.noiseq_norm_upper().to_string(), "UpperQuartile")
  assert_eq(@src.noiseq_norm_none().to_string(), "none")
}

test "noiseq_method_to_string" {
  assert_eq(@src.noiseq_method_bio().to_string(), "NOISeqBio")
  assert_eq(@src.noiseq_method_sim().to_string(), "NOISeqSim")
}

test "noiseq_sample_data" {
  let (ctrl, trt) = @src.noiseq_sample_data()
  assert_eq(ctrl.length(), 3)
  assert_eq(trt.length(), 3)
  assert_eq(ctrl[0].condition, "control")
  assert_eq(trt[0].condition, "treatment")
}

test "noiseq_result_ranking" {
  let (ctrl, trt) = @src.noiseq_sample_data()
  let results = @src.noiseq_run(ctrl, trt)
  // Top gene should have ranking 1
  assert_eq(results.results[0].ranking, 1)
  assert_eq(results.results[1].ranking, 2)
}
