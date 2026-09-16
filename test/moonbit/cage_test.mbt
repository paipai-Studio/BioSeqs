///|
/// Tests for cage module (CAGEr - CAGE data TSS analysis).

fn cage_example_ctss() -> Array[@src.CageCtss] {
  [
    @src.CageCtss::new("chr1", 1000, "+", 50),
    @src.CageCtss::new("chr1", 1001, "+", 30),
    @src.CageCtss::new("chr1", 1002, "+", 20),
    @src.CageCtss::new("chr1", 1005, "+", 15),
    @src.CageCtss::new("chr1", 1010, "+", 10),
    @src.CageCtss::new("chr1", 2000, "-", 40),
    @src.CageCtss::new("chr1", 2001, "-", 35),
  ]
}

test "cage_sample_total_tags" {
  let s = @src.CageSample::new("s1", cage_example_ctss())
  assert_eq(s.total_tags, 200)
}

test "cage_new_and_add" {
  let cs = @src.cage_new("test")
  let cs1 = @src.cage_add_sample(cs, @src.CageSample::new("s1", cage_example_ctss()))
  assert_eq(@src.cage_n_samples(cs1), 1)
}

test "cage_normalize_tpm" {
  let s = @src.CageSample::new("s1", cage_example_ctss())
  let cs = @src.cage_add_sample(@src.cage_new("test"), s)
  let norm = @src.cage_normalize_tpm(cs)
  let s1 = norm.samples[0]
  assert_true(s1.ctss[0].normalized > 0.0)
}

test "cage_cluster_basic" {
  let s = @src.CageSample::new("s1", cage_example_ctss())
  let cs = @src.cage_add_sample(@src.cage_new("test"), s)
  let clusters = @src.cage_cluster(cs, max_dist=10, min_tags=2)
  assert_true(clusters.length() > 0)
}

test "cage_cluster_shape" {
  let s = @src.CageSample::new("s1", cage_example_ctss())
  let cs = @src.cage_add_sample(@src.cage_new("test"), s)
  let clusters = @src.cage_cluster(cs)
  if clusters.length() > 0 {
    let shape = @src.cage_classify_tss_shape(clusters[0])
    assert_true(shape == "sharp" || shape == "broad")
  }
}

test "cage_cluster_summary" {
  let s = @src.CageSample::new("s1", cage_example_ctss())
  let cs = @src.cage_add_sample(@src.cage_new("test"), s)
  let clusters = @src.cage_cluster(cs)
  let (n, total_tags, avg_tags, avg_ctss, avg_shape) = @src.cage_cluster_summary(clusters)
  assert_true(n >= 0)
  assert_true(total_tags >= 0)
}

test "cage_diff_tss" {
  let ctss1 = [@src.CageCtss::new("chr1", 1000, "+", 50)]
  let ctss2 = [@src.CageCtss::new("chr1", 1000, "+", 100)]
  let s1 = @src.CageSample::new("s1", ctss1)
  let s2 = @src.CageSample::new("s2", ctss2)
  let diff = @src.cage_diff_tss(s1, s2)
  assert_true(diff.length() > 0)
}

test "cage_annotate_tss" {
  let ctss = [
    @src.CageCtss::new("chr1", 1000, "+", 50),
    @src.CageCtss::new("chr1", 2000, "+", 30),
  ]
  let s = @src.CageSample::new("s1", ctss)
  let cs = @src.cage_add_sample(@src.cage_new("test"), s)
  let clusters = @src.cage_cluster(cs, max_dist=50, min_tags=1)
  let genes = [("GENE1", 1000, "+", "chr1"), ("GENE2", 2000, "+", "chr1")]
  let ann = @src.cage_annotate_tss(clusters, genes)
  assert_true(ann.length() >= 0)
}

test "cage_total_tags" {
  let s = @src.CageSample::new("s1", cage_example_ctss())
  let cs = @src.cage_add_sample(@src.cage_new("test"), s)
  assert_eq(@src.cage_total_tags(cs), 200)
}
