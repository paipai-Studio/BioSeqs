///|
/// scikit-bio distance-based statistical tests demo.
///
/// Runs ANOSIM, PERMDISP and PERMANOVA on a small Euclidean distance matrix
/// with two groups (a tight cluster vs a spread cluster).

///|
fn main {
  let dm = [
    [0.0, 0.223607, 0.111803, 2.828427, 3.162278, 3.590265],
    [0.223607, 0.0, 0.304138, 2.617250, 3.008322, 3.373426],
    [0.111803, 0.304138, 0.0, 2.865746, 3.148412, 3.677295],
    [2.828427, 2.617250, 2.865746, 0.0, 1.414214, 1.920937],
    [3.162278, 3.008322, 3.148412, 1.414214, 0.0, 3.330165],
    [3.590265, 3.373426, 3.677295, 1.920937, 3.330165, 0.0],
  ]
  let g = [0, 0, 0, 1, 1, 1]

  let a = @src.skbio_anosim(dm, g, 999, 42)
  println("ANOSIM    R   = \{a.statistic}")
  println("          p   = \{a.p_value}")

  let pd = @src.skbio_permdisp(dm, g, 999, 42)
  println("PERMDISP  F   = \{pd.statistic}")
  println("          p   = \{pd.p_value}")

  let p = @src.skbio_permanova(dm, g, 999, 42)
  println("PERMANOVA F   = \{p.statistic}")
  println("          p   = \{p.p_value}")
  println("          n = \{p.n}, groups = \{p.num_groups}")
}