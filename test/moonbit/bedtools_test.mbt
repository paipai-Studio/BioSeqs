///|
/// Tests for bedtools-style interval arithmetic, cross-validated against a
/// Python reference implementation of the bedtools `merge`/`intersect` default
/// semantics, plus a brute-force check of the IntervalTree.

///|
fn sort_ints(a : Array[Int]) -> Array[Int] {
  let out = a.copy()
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && out[j] > key {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

///|
fn ints_eq(a : Array[Int], b : Array[Int]) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i = 0; i < a.length(); i = i + 1 {
    if a[i] != b[i] {
      return false
    }
  }
  true
}

///|
fn iv(chrom : String, start : Int, end : Int) -> @src.BedInterval {
  @src.BedInterval::new(chrom, start, end)
}

///|
test "bedtools_merge" {
  let merged = @src.bedtools_merge([
    iv("chr1", 0, 10),
    iv("chr1", 5, 20),
    iv("chr1", 20, 30),
    iv("chr1", 40, 50),
    iv("chr2", 0, 5),
    iv("chr2", 5, 10),
  ])
  assert_eq(merged.length(), 3)
  assert_eq(merged[0].chrom, "chr1")
  assert_eq(merged[0].start, 0)
  assert_eq(merged[0].end, 30)
  assert_eq(merged[1].start, 40)
  assert_eq(merged[1].end, 50)
  assert_eq(merged[2].chrom, "chr2")
  assert_eq(merged[2].start, 0)
  assert_eq(merged[2].end, 10)
}

///|
test "bedtools_merge_unsorted_input" {
  let merged = @src.bedtools_merge([iv("chr1", 10, 20), iv("chr1", 30, 40), iv("chr1", 25, 35)])
  assert_eq(merged.length(), 2)
  assert_eq(merged[0].start, 10)
  assert_eq(merged[0].end, 20)
  assert_eq(merged[1].start, 25)
  assert_eq(merged[1].end, 40)
}

///|
test "bedtools_intersect" {
  let a = [iv("chr1", 0, 10), iv("chr1", 15, 25), iv("chr1", 100, 110)]
  let b = [iv("chr1", 5, 20), iv("chr1", 200, 300)]
  let hit = @src.bedtools_intersect(a, b)
  assert_eq(hit.length(), 2)
  assert_eq(hit[0].start, 0)
  assert_eq(hit[1].start, 15)
  // no overlap
  let none = @src.bedtools_intersect([iv("chr1", 0, 10), iv("chr1", 50, 60)], [iv("chr1", 20, 30)])
  assert_eq(none.length(), 0)
}

///|
test "interval_tree_bruteforce" {
  let n = 200
  let starts = Array::make(n, 0)
  let ends = Array::make(n, 0)
  let mut rng = 12345
  for i = 0; i < n; i = i + 1 {
    rng = (rng * 1103515245 + 12345) & 0x7FFFFFFF
    let s = rng % 1000
    rng = (rng * 1103515245 + 12345) & 0x7FFFFFFF
    let l = 1 + rng % 50
    starts[i] = s
    ends[i] = s + l
  }
  let tree = @src.interval_tree_new(starts, ends)
  for q = 0; q < 300; q = q + 1 {
    rng = (rng * 1103515245 + 12345) & 0x7FFFFFFF
    let a = rng % 1000
    rng = (rng * 1103515245 + 12345) & 0x7FFFFFFF
    let b = a + 1 + rng % 100
    let got = sort_ints(tree.query(a, b))
    let want : Array[Int] = Array::new()
    for i = 0; i < n; i = i + 1 {
      if starts[i] < b && ends[i] > a {
        want.push(i)
      }
    }
    assert_true(ints_eq(got, want))
  }
}
///|
test "bedtools_closest" {
  let a = [iv("chr1", 100, 200), iv("chr1", 500, 600), iv("chr1", 900, 1000)]
  let b = [iv("chr1", 450, 480), iv("chr1", 250, 300), iv("chr1", 700, 750)]
  let r = @src.bedtools_closest(a, b)
  assert_eq(r.length(), 3)
  // A(100,200) -> B(250,300) 下游 +51（gap+1）
  assert_eq(r[0].b.start, 250)
  assert_eq(r[0].distance, 51)
  // A(500,600) -> B(450,480) 上游 -21
  assert_eq(r[1].b.start, 450)
  assert_eq(r[1].distance, -21)
  // A(900,1000) -> B(700,750) 上游 -151
  assert_eq(r[2].b.start, 700)
  assert_eq(r[2].distance, -151)
  // 重叠 distance 0
  let r2 = @src.bedtools_closest([iv("chr1", 100, 200)], [iv("chr1", 150, 160)])
  assert_eq(r2[0].distance, 0)
}

test "bedtools_window" {
  let a = [iv("chr1", 100, 200), iv("chr1", 500, 600)]
  let b = [iv("chr1", 130, 140), iv("chr1", 230, 240), iv("chr1", 480, 490)]
  let r = @src.bedtools_window(a, b, 50)
  assert_eq(r.length(), 2)
  assert_eq(r[0].start, 100)
  assert_eq(r[1].start, 500)
  // w=10：fudge [90,210) 命中 130-140；[490,610) 不命中 480-490（490 是边界）
  let r2 = @src.bedtools_window(a, b, 10)
  assert_eq(r2.length(), 1)
  assert_eq(r2[0].start, 100)
  // A.end 下游附近命中（旧「围绕 start」语义会漏）
  let r3 = @src.bedtools_window([iv("chr1", 100, 200)], [iv("chr1", 205, 210)], 10)
  assert_eq(r3.length(), 1)
}

///|
test "bedtools_coverage" {
  let a = [iv("chr1", 0, 10), iv("chr1", 20, 30)]
  let b = [iv("chr1", 3, 7), iv("chr1", 6, 9), iv("chr1", 25, 28)]
  let r = @src.bedtools_coverage(a, b)
  assert_eq(r.length(), 2)
  // A(0,10): count=2, covered=6, len=10, frac=0.6
  assert_eq(r[0].count, 2)
  assert_eq(r[0].bases_covered, 6)
  assert_eq(r[0].length, 10)
  assert_true(assert_close(r[0].fraction, 0.6, 1.0e-9))
  // A(20,30): count=1, covered=3, len=10, frac=0.3
  assert_eq(r[1].count, 1)
  assert_eq(r[1].bases_covered, 3)
  assert_true(assert_close(r[1].fraction, 0.3, 1.0e-9))
}

///|
test "bedtools_genomecov_bg" {
  let ivs = [iv("chr1", 3, 7), iv("chr1", 6, 9), iv("chr1", 25, 28)]
  let r = @src.bedtools_genomecov_bg(ivs)
  assert_eq(r.length(), 4)
  assert_eq(r[0].start, 3)
  assert_eq(r[0].end, 6)
  assert_eq(r[0].depth, 1)
  assert_eq(r[1].start, 6)
  assert_eq(r[1].end, 7)
  assert_eq(r[1].depth, 2)
  assert_eq(r[2].start, 7)
  assert_eq(r[2].end, 9)
  assert_eq(r[2].depth, 1)
  assert_eq(r[3].start, 25)
  assert_eq(r[3].end, 28)
  assert_eq(r[3].depth, 1)
}

///|
test "bedtools_subtract_complement" {
  let sizes = [("chr1", 1000), ("chr2", 500)]
  // subtract
  let s = @src.bedtools_subtract([iv("chr1", 0, 100)], [iv("chr1", 30, 50), iv("chr1", 70, 80)])
  assert_eq(s.length(), 3)
  assert_eq(s[0].start, 0)
  assert_eq(s[0].end, 30)
  assert_eq(s[1].start, 50)
  assert_eq(s[1].end, 70)
  assert_eq(s[2].start, 80)
  assert_eq(s[2].end, 100)
  // 全被覆盖 → 空
  assert_eq(@src.bedtools_subtract([iv("chr1", 0, 100)], [iv("chr1", 0, 100)]).length(), 0)
  // complement
  let c = @src.bedtools_complement([iv("chr1", 100, 200), iv("chr1", 300, 400)], sizes)
  assert_eq(c.length(), 4)
  assert_eq(c[0].start, 0)
  assert_eq(c[0].end, 100)
  assert_eq(c[1].start, 200)
  assert_eq(c[1].end, 300)
  assert_eq(c[2].start, 400)
  assert_eq(c[2].end, 1000)
  assert_eq(c[3].chrom, "chr2")
  assert_eq(c[3].start, 0)
  assert_eq(c[3].end, 500)
}

///|
test "bedtools_flank_shift_slop" {
  let sizes = [("chr1", 1000)]
  // flank
  let f = @src.bedtools_flank([iv("chr1", 100, 200)], sizes, 50, 50)
  assert_eq(f.length(), 2)
  assert_eq(f[0].start, 50)
  assert_eq(f[0].end, 100)
  assert_eq(f[1].start, 200)
  assert_eq(f[1].end, 250)
  // shift
  let sh1 = @src.bedtools_shift([iv("chr1", 100, 200)], sizes, 50)
  assert_eq(sh1[0].start, 150)
  assert_eq(sh1[0].end, 250)
  let sh2 = @src.bedtools_shift([iv("chr1", 100, 200)], sizes, -50)
  assert_eq(sh2[0].start, 50)
  assert_eq(sh2[0].end, 150)
  // slop
  let sl = @src.bedtools_slop([iv("chr1", 100, 200)], sizes, 50, 50)
  assert_eq(sl[0].start, 50)
  assert_eq(sl[0].end, 250)
  let sl2 = @src.bedtools_slop([iv("chr1", 10, 20)], sizes, 50, 50)
  assert_eq(sl2[0].start, 0)
  assert_eq(sl2[0].end, 70)
}

///|
test "bedtools_intersect_frac_stranded" {
  let a = [iv("chr1", 0, 100)]
  // overlay 50/100 = 0.5
  let frac_hit = @src.bedtools_intersect_frac(a, [iv("chr1", 50, 100)], 0.5)
  assert_eq(frac_hit.length(), 1)
  let frac_miss = @src.bedtools_intersect_frac(a, [iv("chr1", 50, 100)], 0.6)
  assert_eq(frac_miss.length(), 0)
  // strand：A '+'，B 一条 '+' 一条 '-'
  let sa = ["+"]
  let sb = ["+", "-"]
  let hit = @src.bedtools_intersect_stranded([iv("chr1", 0, 10)], sa, [iv("chr1", 0, 10), iv("chr1", 0, 10)], sb)
  assert_eq(hit.length(), 1)
}

///|
test "bedtools_map" {
  let a = [iv("chr1", 0, 100)]
  let b = [iv("chr1", 0, 50), iv("chr1", 25, 75), iv("chr1", 200, 300)]
  let bv = [10.0, 20.0, 5.0]
  let cnt = @src.bedtools_map_count(a, b)
  assert_eq(cnt[0], 2)
  let sum = @src.bedtools_map_sum(a, b, bv)
  assert_true(assert_close(sum[0], 30.0, 1.0e-9))
  let mean = @src.bedtools_map_mean(a, b, bv)
  assert_true(assert_close(mean[0], 15.0, 1.0e-9))
  let mn = @src.bedtools_map_min(a, b, bv)
  assert_true(assert_close(mn[0], 10.0, 1.0e-9))
  let mx = @src.bedtools_map_max(a, b, bv)
  assert_true(assert_close(mx[0], 20.0, 1.0e-9))
}

///|
test "bedtools_groupby" {
  let records = [
    ["chr1", "100", "200", "10"],
    ["chr1", "300", "400", "20"],
    ["chr1", "350", "450", "30"],
    ["chr2", "50", "60", "5"],
  ]
  // sum by chrom (col 1), aggregate col 4
  let sum = @src.bedtools_groupby(records, [1], 4, "sum")
  assert_eq(sum.length(), 2)
  assert_eq(sum[0].group[0], "chr1")
  assert_true(assert_close(sum[0].value, 60.0, 1.0e-9))
  assert_eq(sum[1].group[0], "chr2")
  assert_true(assert_close(sum[1].value, 5.0, 1.0e-9))
  // count
  let cnt = @src.bedtools_groupby(records, [1], 4, "count")
  assert_true(assert_close(cnt[0].value, 3.0, 1.0e-9))
  assert_true(assert_close(cnt[1].value, 1.0, 1.0e-9))
  // mean / min / max
  let mean = @src.bedtools_groupby(records, [1], 4, "mean")
  assert_true(assert_close(mean[0].value, 20.0, 1.0e-9))
  let mn = @src.bedtools_groupby(records, [1], 4, "min")
  assert_true(assert_close(mn[0].value, 10.0, 1.0e-9))
  let mx = @src.bedtools_groupby(records, [1], 4, "max")
  assert_true(assert_close(mx[0].value, 30.0, 1.0e-9))
}
///|
test "bedtools_map_text_ops" {
  let a = [iv("chr1", 0, 100)]
  let b = [iv("chr1", 0, 50), iv("chr1", 25, 75), iv("chr1", 50, 60)]
  let bs = ["a", "b", "a"]
  assert_eq(@src.bedtools_map_mode(a, b, bs)[0], "a")
  assert_eq(@src.bedtools_map_antimode(a, b, bs)[0], "b")
  assert_eq(@src.bedtools_map_distinct(a, b, bs)[0], "a,b")
  assert_eq(@src.bedtools_map_distinct_only(a, b, bs)[0], "b")
  assert_eq(@src.bedtools_map_collapse(a, b, bs)[0], "a,b,a")
  assert_eq(@src.bedtools_map_concat(a, b, bs)[0], "aba")
  assert_eq(@src.bedtools_map_first(a, b, bs)[0], "a")
  assert_eq(@src.bedtools_map_last(a, b, bs)[0], "a")
  assert_eq(@src.bedtools_map_count_distinct(a, b, bs)[0], 2)
}

///|
test "bedtools_map_numeric_extended" {
  let a = [iv("chr1", 0, 100)]
  let b = [iv("chr1", 0, 50), iv("chr1", 25, 75), iv("chr1", 200, 300)]
  let bv = [10.0, 20.0, 5.0]
  // overlap values: [10, 20]
  assert_true(assert_close(@src.bedtools_map_median(a, b, bv)[0], 15.0, 1.0e-9))
  // population stddev (÷N)
  assert_true(assert_close(@src.bedtools_map_stdev(a, b, bv)[0], 5.0, 1.0e-9))
  // sample stddev (÷N-1)
  assert_true(assert_close(@src.bedtools_map_sstdev(a, b, bv)[0], 7.0710678118654755, 1.0e-9))
  assert_true(assert_close(@src.bedtools_map_absmin(a, b, bv)[0], 10.0, 1.0e-9))
  assert_true(assert_close(@src.bedtools_map_absmax(a, b, bv)[0], 20.0, 1.0e-9))
}

///|
test "bedtools_groupby_extended" {
  let records = [
    ["chr1", "100", "200", "10"],
    ["chr1", "300", "400", "20"],
    ["chr1", "350", "450", "30"],
    ["chr2", "50", "60", "5"],
  ]
  let median = @src.bedtools_groupby(records, [1], 4, "median")
  assert_true(assert_close(median[0].value, 20.0, 1.0e-9))
  assert_true(assert_close(median[1].value, 5.0, 1.0e-9))
  let stdev = @src.bedtools_groupby(records, [1], 4, "stdev")
  assert_true(assert_close(stdev[0].value, 8.16496580927726, 1.0e-9))
  assert_true(assert_close(stdev[1].value, 0.0, 1.0e-9))
  let sstdev = @src.bedtools_groupby(records, [1], 4, "sstdev")
  assert_true(assert_close(sstdev[0].value, 10.0, 1.0e-9))
  assert_true(assert_close(sstdev[1].value, 0.0, 1.0e-9))
  let absmin = @src.bedtools_groupby(records, [1], 4, "absmin")
  assert_true(assert_close(absmin[0].value, 10.0, 1.0e-9))
  let absmax = @src.bedtools_groupby(records, [1], 4, "absmax")
  assert_true(assert_close(absmax[0].value, 30.0, 1.0e-9))
  let cdist = @src.bedtools_groupby(records, [1], 4, "count_distinct")
  assert_true(assert_close(cdist[0].value, 3.0, 1.0e-9))
  assert_true(assert_close(cdist[1].value, 1.0, 1.0e-9))
}

///|
test "bedtools_groupby_text" {
  let records = [
    ["chr1", "x", "a"],
    ["chr1", "y", "b"],
    ["chr1", "z", "a"],
    ["chr2", "w", "c"],
  ]
  let mode = @src.bedtools_groupby_text(records, [1], 3, "mode")
  assert_eq(mode[0].group[0], "chr1")
  assert_eq(mode[0].value, "a")
  assert_eq(mode[1].value, "c")
  assert_eq(@src.bedtools_groupby_text(records, [1], 3, "antimode")[0].value, "b")
  let distinct = @src.bedtools_groupby_text(records, [1], 3, "distinct")
  assert_eq(distinct[0].value, "a,b")
  assert_eq(distinct[1].value, "c")
  let distinct_only = @src.bedtools_groupby_text(records, [1], 3, "distinct_only")
  assert_eq(distinct_only[0].value, "b")
  assert_eq(@src.bedtools_groupby_text(records, [1], 3, "collapse")[0].value, "a,b,a")
  assert_eq(@src.bedtools_groupby_text(records, [1], 3, "concat")[0].value, "aba")
  let first = @src.bedtools_groupby_text(records, [1], 3, "first")
  assert_eq(first[0].value, "a")
  assert_eq(first[1].value, "c")
  let last = @src.bedtools_groupby_text(records, [1], 3, "last")
  assert_eq(last[0].value, "a")
  assert_eq(last[1].value, "c")
}///|
test "bedtools_intersect_opposite_strand" {
  let a = [iv("chr1", 0, 10)]
  let sa = ["+"]
  // B has one "+" and one "-"; only "-" is opposite.
  let b = [iv("chr1", 0, 10), iv("chr1", 0, 10)]
  let sb = ["+", "-"]
  assert_eq(@src.bedtools_intersect_opposite_strand(a, sa, b, sb).length(), 1)
  // both same strand → no opposite hit
  let b2 = [iv("chr1", 0, 10), iv("chr1", 0, 10)]
  let sb2 = ["+", "+"]
  assert_eq(@src.bedtools_intersect_opposite_strand(a, sa, b2, sb2).length(), 0)
  // unknown strand never satisfies -S
  let b3 = [iv("chr1", 0, 10)]
  let sb3 = ["."]
  assert_eq(@src.bedtools_intersect_opposite_strand(a, sa, b3, sb3).length(), 0)
}
///|
test "bedtools_intersect_frac_b_reciprocal" {
  // A(0,100) vs B(90,110): overlap=10, B.len=20 → overlap/B = 0.5
  let a = [iv("chr1", 0, 100)]
  let frac_b_hit = @src.bedtools_intersect_frac_b(a, [iv("chr1", 90, 110)], 0.5)
  assert_eq(frac_b_hit.length(), 1)
  let frac_b_miss = @src.bedtools_intersect_frac_b(a, [iv("chr1", 90, 110)], 0.6)
  assert_eq(frac_b_miss.length(), 0)
  // reciprocal: A(0,100) vs B(80,200): overlap=20, A.len=100, B.len=120
  // overlap/A=0.2, overlap/B≈0.1667
  let rec_hit = @src.bedtools_intersect_reciprocal(a, [iv("chr1", 80, 200)], 0.1)
  assert_eq(rec_hit.length(), 1)
  let rec_miss = @src.bedtools_intersect_reciprocal(a, [iv("chr1", 80, 200)], 0.2)
  assert_eq(rec_miss.length(), 0)
}///|
test "bedtools_jaccard" {
  // a=[0,100), b=[50,150): a_len=100, b_len=100, overlap=50
  let r = @src.bedtools_jaccard([iv("chr1", 0, 100)], [iv("chr1", 50, 150)])
  assert_eq(r.intersection, 50)
  assert_eq(r.union, 150)
  assert_eq(r.n_intersections, 1)
  assert_true(assert_close(r.jaccard, 1.0 / 3.0, 1.0e-9))
  // disjoint -> intersection 0, jaccard 0
  let r2 = @src.bedtools_jaccard([iv("chr1", 0, 100)], [iv("chr1", 200, 300)])
  assert_eq(r2.intersection, 0)
  assert_eq(r2.n_intersections, 0)
  assert_true(assert_close(r2.jaccard, 0.0, 1.0e-12))
}

///|
test "bedtools_makewindows" {
  // chr1:10, size=3, step=3 -> [0,3),[3,6),[6,9),[9,10)
  let w = @src.bedtools_makewindows([("chr1", 10)], 3, 3)
  assert_eq(w.length(), 4)
  assert_eq((w[0].start, w[0].end), (0, 3))
  assert_eq((w[1].start, w[1].end), (3, 6))
  assert_eq((w[2].start, w[2].end), (6, 9))
  assert_eq((w[3].start, w[3].end), (9, 10))
  // overlapping windows when step < size
  let w2 = @src.bedtools_makewindows([("chr1", 10)], 4, 3)
  assert_eq(w2.length(), 4)
  assert_eq((w2[2].start, w2[2].end), (6, 10))
}

///|
test "bedtools_makewindows_count" {
  // chr1:10, count=3 -> [0,3),[3,6),[6,10)
  let w = @src.bedtools_makewindows_count([("chr1", 10)], 3)
  assert_eq(w.length(), 3)
  assert_eq((w[0].start, w[0].end), (0, 3))
  assert_eq((w[1].start, w[1].end), (3, 6))
  assert_eq((w[2].start, w[2].end), (6, 10))
  // count larger than length -> skipped
  let w2 = @src.bedtools_makewindows_count([("chr1", 3)], 5)
  assert_eq(w2.length(), 0)
}///|
test "bedtools_reldist" {
  // B midpoints: [25, 75]
  let b = [iv("chr1", 0, 50), iv("chr1", 50, 100)]
  // A(0,100) midpoint 50 -> min(25,25)/(75-25) = 0.5
  // A(0,50)  midpoint 25 -> coincides with left B midpoint -> 0.0
  let a = [iv("chr1", 0, 100), iv("chr1", 0, 50)]
  let r = @src.bedtools_reldist(a, b)
  assert_eq(r.length(), 2)
  assert_eq(r[0].rel_dist, 0.0)
  assert_eq(r[0].count, 1)
  assert_eq(r[0].total, 2)
  assert_true(assert_close(r[0].fraction, 0.5, 1.0e-9))
  assert_eq(r[1].rel_dist, 0.5)
  assert_eq(r[1].count, 1)
}///|
test "bedtools_getfasta_nuc" {
  // chr1 = "AACCCGGGTNNX" (A2 C3 G3 T1 N2 X1, len 12)
  let fasta = @src.Fasta::from_content(">chr1\nAACCCGGGTNNX\n")
  // getfasta
  let f = @src.bedtools_getfasta(fasta, [iv("chr1", 0, 4), iv("chr1", 5, 9)])
  assert_eq(f.length(), 2)
  assert_eq(f[0].header, ">chr1:0-4")
  assert_eq(f[0].seq, "AACC")
  assert_eq(f[1].header, ">chr1:5-9")
  assert_eq(f[1].seq, "GGGT")
  // nuc
  let r = @src.bedtools_nuc(fasta, [iv("chr1", 0, 12)])
  assert_eq(r.length(), 1)
  assert_eq(r[0].a, 2)
  assert_eq(r[0].c, 3)
  assert_eq(r[0].g, 3)
  assert_eq(r[0].t, 1)
  assert_eq(r[0].n, 2)
  assert_eq(r[0].other, 1)
  assert_eq(r[0].length, 12)
  assert_true(assert_close(r[0].pct_at, 0.25, 1.0e-9))
  assert_true(assert_close(r[0].pct_gc, 0.5, 1.0e-9))
  // out-of-range interval is skipped
  let r2 = @src.bedtools_nuc(fasta, [iv("chr1", 0, 100)])
  assert_eq(r2.length(), 0)
}