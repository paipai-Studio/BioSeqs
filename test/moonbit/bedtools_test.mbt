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