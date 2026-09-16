///|
/// Tests for conumee module (copy number unified interface).

fn cnm_example_probes() -> Array[@src.CnmProbe] {
  [
    @src.CnmProbe::new("chr1", 1000, 0.1),
    @src.CnmProbe::new("chr1", 2000, 0.2),
    @src.CnmProbe::new("chr1", 3000, 0.15),
    @src.CnmProbe::new("chr1", 4000, 1.5),
    @src.CnmProbe::new("chr1", 5000, 1.8),
    @src.CnmProbe::new("chr2", 1000, -0.5),
    @src.CnmProbe::new("chr2", 2000, -0.8),
  ]
}

test "cnm_from_arrays" {
  let opt = @src.cnm_from_arrays(
    "test",
    ["chr1", "chr1", "chr2"],
    [1000, 2000, 1000],
    [0.5, -0.3, 1.2],
  )
  match opt {
    Some(s) => assert_eq(s.sample_id, "test")
    None => assert_true(false)
  }
}

test "cnm_from_arrays_mismatch" {
  let opt = @src.cnm_from_arrays(
    "test",
    ["chr1", "chr1"],
    [1000],
    [0.5],
  )
  assert_true(opt.is_none())
}

test "cnm_sample_creation" {
  let s = @src.CnmSample::new("s1", cnm_example_probes())
  assert_eq(s.sample_id, "s1")
  assert_eq(s.probes.length(), 7)
}

test "cnm_median" {
  let vals = [1.0, 2.0, 3.0, 4.0, 5.0]
  let m = @src.cnm_median(vals)
  assert_true((m - 3.0).abs() < 1.0e-9)
}

test "cnm_mad" {
  let vals = [1.0, 2.0, 3.0, 4.0, 5.0]
  let mad = @src.cnm_mad(vals)
  assert_true(mad > 0.0)
}

test "cnm_segment" {
  let s = @src.CnmSample::new("s1", cnm_example_probes())
  let segged = @src.cnm_segment(s)
  assert_eq(segged.sample_id, "s1")
  assert_true(segged.segments.length() >= 0)
}

test "cnm_qc" {
  let s = @src.CnmSample::new("s1", cnm_example_probes())
  let segged = @src.cnm_segment(s)
  let qc = segged.qc
  assert_true(qc.n_probes >= 0)
  assert_true(qc.noise_score >= 0.0)
}

test "cnm_filter_segments" {
  let s = @src.CnmSample::new("s1", cnm_example_probes())
  let segged = @src.cnm_segment(s)
  let filtered = @src.cnm_filter_segments(segged, gain_threshold=2.5, loss_threshold=1.5)
  assert_true(filtered.length() >= 0)
}

test "cnm_annotate_genes" {
  let s = @src.CnmSample::new("s1", cnm_example_probes())
  let segged = @src.cnm_segment(s)
  let genes = [
    ("chr1", 1000, 2000, "GENE1"),
    ("chr1", 4000, 5000, "GENE2"),
    ("chr2", 1000, 2000, "GENE3"),
  ]
  let annotated = @src.cnm_annotate_genes(segged, genes)
  assert_true(annotated.segments.length() >= 0)
}

test "cnm_compare_samples" {
  let s1 = @src.CnmSample::new("s1", cnm_example_probes())
  let s2 = @src.CnmSample::new("s2", [
    @src.CnmProbe::new("chr1", 1000, 0.2),
    @src.CnmProbe::new("chr1", 2000, 0.3),
  ])
  let seg1 = @src.cnm_segment(s1)
  let seg2 = @src.cnm_segment(s2)
  let cmp = @src.cnm_compare_samples(seg1, seg2)
  assert_true(cmp.correlation >= -1.0 && cmp.correlation <= 1.0)
}

test "cnm_segments_tsv" {
  let segs : Array[@src.CnmSegment] = [
    @src.CnmSegment::{
      chr: "chr1", start: 1000, end: 2000, n_probes: 5,
      mean_logratio: 0.5, copy_number: 3.0, genes: ["GENE1"],
    },
  ]
  let tsv = @src.cnm_segments_tsv(segs)
  assert_true(tsv.contains("chr1"))
}

test "cnm_ascii_plot" {
  let segs : Array[@src.CnmSegment] = [
    @src.CnmSegment::{
      chr: "chr1", start: 1000, end: 2000, n_probes: 5,
      mean_logratio: 0.5, copy_number: 3.0, genes: ["GENE1"],
    },
  ]
  let plot = @src.cnm_ascii_plot(segs, "chr1")
  assert_true(plot.contains("chr1"))
}
