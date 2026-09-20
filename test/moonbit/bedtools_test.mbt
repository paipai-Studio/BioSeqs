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
