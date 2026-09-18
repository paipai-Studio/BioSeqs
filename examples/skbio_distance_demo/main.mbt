///|
/// scikit-bio distance-based statistical tests demo.
///
/// Runs ANOSIM and PERMANOVA on a small distance matrix with two groups.

///|
fn main {
  let dm = [[0.0, 1.0, 1.0, 4.0], [1.0, 0.0, 3.0, 2.0], [1.0, 3.0, 0.0, 3.0], [4.0, 2.0, 3.0, 0.0]]
  let g = [0, 0, 1, 1]

  let a = @src.skbio_anosim(dm, g, 999, 42)
  println("ANOSIM    R   = \{a.statistic}")
  println("          p   = \{a.p_value}")
  println("          n   = \{a.n}, groups = \{a.num_groups}")

  let p = @src.skbio_permanova(dm, g, 999, 42)
  println("PERMANOVA F   = \{p.statistic}")
  println("          p   = \{p.p_value}")
  println("          n   = \{p.n}, groups = \{p.num_groups}")
}