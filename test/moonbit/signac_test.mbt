///|
/// Tests for Signac (scATAC) core: feature matrix, TSS enrichment,
/// nucleosome signal.

///|
fn siv(chrom : String, start : Int, end : Int) -> @src.BedInterval {
  @src.BedInterval::new(chrom, start, end)
}

///|
fn sfrag(chrom : String, start : Int, end : Int, cell : String) -> @src.SignacFragment {
  @src.SignacFragment::new(chrom, start, end, cell)
}

///|
test "signac_feature_matrix" {
  let peaks = [siv("chr1", 0, 100), siv("chr1", 200, 300)]
  let cells = ["cell1", "cell2"]
  let fragments = [
    sfrag("chr1", 10, 50, "cell1"),
    sfrag("chr1", 30, 80, "cell1"),
    sfrag("chr1", 150, 250, "cell2"),
    sfrag("chr1", 210, 220, "cell1"),
  ]
  let fm = @src.signac_feature_matrix(fragments, peaks, cells)
  assert_eq(fm.peaks.length(), 2)
  assert_eq(fm.cells.length(), 2)
  assert_eq(fm.counts[0][0], 2)
  assert_eq(fm.counts[0][1], 0)
  assert_eq(fm.counts[1][0], 1)
  assert_eq(fm.counts[1][1], 1)
}

///|
test "signac_feature_matrix_ignores_unknown_cell" {
  let peaks = [siv("chr1", 0, 100)]
  let cells = ["cell1"]
  let fragments = [sfrag("chr1", 10, 20, "cellX"), sfrag("chr1", 30, 40, "cell1")]
  let fm = @src.signac_feature_matrix(fragments, peaks, cells)
  assert_eq(fm.counts[0][0], 1)
}

///|
test "signac_nucleosome_signal" {
  // 100,120 < 147 (nuc-free); 150,160,200,293 in [147,294) (mono); 300 ignored
  let fragments = [
    sfrag("chr1", 0, 100, "c"),
    sfrag("chr1", 0, 120, "c"),
    sfrag("chr1", 0, 150, "c"),
    sfrag("chr1", 0, 160, "c"),
    sfrag("chr1", 0, 200, "c"),
    sfrag("chr1", 0, 293, "c"),
    sfrag("chr1", 0, 300, "c"),
  ]
  assert_eq(@src.signac_nucleosome_signal(fragments), 2.0)
}

///|
test "signac_tss_enrichment" {
  let tss = [@src.TssLoc::new("chr1", 1000)]
  // 6 fragments in the TSS window [900,1100); 2 left, 2 right flanks
  let fragments = [
    sfrag("chr1", 910, 920, "c"),
    sfrag("chr1", 930, 940, "c"),
    sfrag("chr1", 950, 960, "c"),
    sfrag("chr1", 970, 980, "c"),
    sfrag("chr1", 990, 1000, "c"),
    sfrag("chr1", 1010, 1020, "c"),
    sfrag("chr1", 0, 50, "c"),
    sfrag("chr1", 100, 150, "c"),
    sfrag("chr1", 1200, 1250, "c"),
    sfrag("chr1", 1500, 1550, "c"),
  ]
  assert_eq(@src.signac_tss_enrichment(fragments, tss, 100, 1000), 3.0)
}

// ============================================================
// Signac 深化：peak calling + gene activity
// ============================================================

test "signac_call_peaks" {
  // chr1 两个重叠 fragment [0,100) 与 [50,150)：coverage>=2 仅在 [50,100)
  let frags = [
    sfrag("chr1", 0, 100, "c1"),
    sfrag("chr1", 50, 150, "c1"),
    sfrag("chr2", 0, 10, "c1"),
  ]
  let peaks = @src.signac_call_peaks(frags, 2, 10)
  assert_eq(peaks.length(), 1)
  assert_eq(peaks[0].chrom, "chr1")
  assert_eq(peaks[0].start, 50)
  assert_eq(peaks[0].end, 100)
}

test "signac_call_peaks_chroms_separate" {
  // 每个染色体各有一个独立峰
  let frags = [
    sfrag("chr1", 0, 100, "c1"),
    sfrag("chr1", 50, 150, "c1"),
    sfrag("chr2", 0, 100, "c1"),
    sfrag("chr2", 50, 150, "c1"),
  ]
  let peaks = @src.signac_call_peaks(frags, 2, 10)
  assert_eq(peaks.length(), 2)
}

test "signac_gene_activity" {
  let peaks = [siv("chr1", 0, 100), siv("chr1", 200, 300)]
  let cells = ["c1", "c2"]
  let fragments = [
    sfrag("chr1", 10, 50, "c1"),
    sfrag("chr1", 30, 80, "c1"),
    sfrag("chr1", 210, 220, "c2"),
  ]
  let fm = @src.signac_feature_matrix(fragments, peaks, cells)
  // gene1 覆盖 peak0（含 upstream promoter），gene2 覆盖 peak1
  let genes = [siv("chr1", 0, 100), siv("chr1", 200, 300)]
  let ga = @src.signac_gene_activity(fm, genes, ["g1", "g2"], 50)
  // peak0 → c1 计数 2，peak1 → c2 计数 1
  assert_eq(ga.activity[0][0], 2)
  assert_eq(ga.activity[0][1], 0)
  assert_eq(ga.activity[1][0], 0)
  assert_eq(ga.activity[1][1], 1)
}

test "signac_gene_activity_promoter_upstream" {
  // peak 落在 gene 的 upstream 区（不在 body 内）也应计入
  let peaks = [siv("chr1", 150, 180)]
  let cells = ["c1"]
  let fragments = [sfrag("chr1", 160, 170, "c1")]
  let fm = @src.signac_feature_matrix(fragments, peaks, cells)
  // gene body [200,300]，upstream 100 → promoter 范围 [100,300]，peak[150,180] 在内
  let genes = [siv("chr1", 200, 300)]
  let ga = @src.signac_gene_activity(fm, genes, ["g_up"], 100)
  assert_eq(ga.activity[0][0], 1)
}
