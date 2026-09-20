///|
/// bedtools-style interval arithmetic demo: merge, intersect, and interval tree.

///|
fn main {
  // merge book-ended / overlapping intervals
  let merged = @src.bedtools_merge([
    @src.BedInterval::new("chr1", 0, 10),
    @src.BedInterval::new("chr1", 5, 20),
    @src.BedInterval::new("chr1", 20, 30),
    @src.BedInterval::new("chr2", 5, 10),
  ])
  println("=== merge ===")
  for iv in merged {
    println("  \{iv.chrom}:\{iv.start}-\{iv.end}")
  }

  // intersect A against B
  let hit = @src.bedtools_intersect(
    [@src.BedInterval::new("chr1", 0, 10), @src.BedInterval::new("chr1", 15, 25)],
    [@src.BedInterval::new("chr1", 5, 20)],
  )
  println("=== intersect (A overlapping B) ===")
  for iv in hit {
    println("  \{iv.chrom}:\{iv.start}-\{iv.end}")
  }

  // interval tree query
  let tree = @src.interval_tree_new([0, 10, 20], [5, 15, 25])
  println("=== interval tree query [8, 12) → \{tree.query(8, 12)}")

  // closest / window
  let cl = @src.bedtools_closest(
    [@src.BedInterval::new("chr1", 100, 200), @src.BedInterval::new("chr1", 500, 600)],
    [@src.BedInterval::new("chr1", 250, 300), @src.BedInterval::new("chr1", 450, 480)],
  )
  for r in cl {
    println("  closest: \{r.a.start}-\{r.a.end} -> \{r.b.start}-\{r.b.end} (d=\{r.distance})")
  }
}
