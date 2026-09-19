///|
/// scikit-bio distance-based statistical tests demo.
///
/// Runs ANOSIM, PERMDISP and PERMANOVA on a small Euclidean distance matrix
/// with two groups (a tight cluster vs a spread cluster).

///|
fn main {
  println("=== BioEnv (best environmental-variable subset) ===")
  let bdm = [
    [0.0, 0.5, 0.25, 1.0, 0.8, 0.6],
    [0.5, 0.0, 0.4, 0.9, 0.7, 0.55],
    [0.25, 0.4, 0.0, 1.1, 0.9, 0.5],
    [1.0, 0.9, 1.1, 0.0, 0.3, 0.75],
    [0.8, 0.7, 0.9, 0.3, 0.0, 0.65],
    [0.6, 0.55, 0.5, 0.75, 0.65, 0.0],
  ]
  let bdata = [
    [7.0, 6.0, 400.0],
    [8.0, 6.5, 380.0],
    [6.5, 6.2, 420.0],
    [10.0, 7.5, 100.0],
    [11.0, 7.8, 80.0],
    [9.0, 7.0, 250.0],
  ]
  let be = @src.skbio_bioenv(bdm, bdata)
  println("best subset correlations: \{be.correlations}")

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