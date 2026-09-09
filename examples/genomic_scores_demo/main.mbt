///|
/// Bioconductor GenomicScores Demo - position-specific genomic scores.
///
/// Demonstrates the port of the Bioconductor `GenomicScores` package
/// (Puigdevall & Ramos, 2018): installing score populations
/// (phastCons / phyloP style) per chromosome, single-position queries
/// with step or linear interpolation, multi-base range summaries
/// (mean/min/max/median) and GRanges-style batch fetch.

///|
fn main raise {
  println("=== GScores: install two populations on chr1 ===")
  let gss = @src.GScores::new("phastCons7way", "hg38")
  gss.set_scores(
    "chr1",
    "phastCons",
    [100, 200, 300, 400],
    [0.1, 0.5, 0.9, 0.4],
  )
  gss.set_scores(
    "chr1",
    "phyloP",
    [100, 200, 300, 400],
    [0.9, 0.7, 0.2, 0.2],
  )
  gss.set_scores(
    "chr2",
    "phastCons",
    [150, 250, 350],
    [0.6, 0.2, 0.8],
  )
  println("populations      : \{gss.populations()}")
  println("chr1 stored bases: \{gss.n_stored_positions("chr1")}")
  println("chr2 stored bases: \{gss.n_stored_positions("chr2")}")

  println("\n=== Per-position scores (step interpolation, default) ===")
  println("phastCons chr1:100 = \{gss.score_at("chr1", 100, population="phastCons")}")
  println("phastCons chr1:150 = \{gss.score_at("chr1", 150, population="phastCons")} (step: keeps 0.1)")
  println("phastCons chr1:250 = \{gss.score_at("chr1", 250, population="phastCons")} (step: keeps 0.5)")
  println("phyloP   chr1:250 = \{gss.score_at("chr1", 250, population="phyloP")}")
  println("phastCons chr1:50  = \{gss.score_at("chr1", 50, population="phastCons")} (before first base: NaN)")

  println("\n=== Linear interpolation between stored bases ===")
  gss.set_interpolation(@src.GsInterpolation::GsLinear)
  println("phastCons chr1:150 = \{gss.score_at("chr1", 150, population="phastCons")} (0.1 -> 0.5 midpoint)")
  println("phastCons chr1:250 = \{gss.score_at("chr1", 250, population="phastCons")}")
  println("phastCons chr1:350 = \{gss.score_at("chr1", 350, population="phastCons")}")
  println("phastCons chr1:450 = \{gss.score_at("chr1", 450, population="phastCons")} (clamped to last)")

  println("\n=== Range summaries (mean / min / max / median) ===")
  println(
    "phastCons chr1:100-400 mean   = \{gss.range_score("chr1", 100, 400, summary="mean", population="phastCons")}",
  )
  println(
    "phastCons chr1:100-400 min    = \{gss.range_score("chr1", 100, 400, summary="min", population="phastCons")}",
  )
  println(
    "phastCons chr1:100-400 max    = \{gss.range_score("chr1", 100, 400, summary="max", population="phastCons")}",
  )
  println(
    "phastCons chr1:100-400 median = \{gss.range_score("chr1", 100, 400, summary="median", population="phastCons")}",
  )

  println("\n=== Default population switch ===")
  gss.set_default_population("phyloP")
  println("default is now phyloP; chr1:300 = \{gss.score_at("chr1", 300)}")
  gss.set_default_population("phastCons")

  println("\n=== Batch fetch over GRanges-like coordinates ===")
  let seqs = ["chr1", "chr1", "chr2", "chr2"]
  let starts = [100, 100, 150, 300]
  let ends = [200, 400, 250, 350]
  let fetched = @src.gscores_fetch(gss, seqs, starts, ends)
  for i in 0..<fetched.length() {
    println("fetch \{seqs[i]}:\{starts[i]}-\{ends[i]} mean = \{fetched[i]}")
  }

  println("\n=== Extending a chromosome with new bases ===")
  gss.set_scores("chr2", "phastCons", [300], [0.4])
  println("chr2 stored bases: \{gss.n_stored_positions("chr2")}")
  println("chr2@300 = \{gss.score_at("chr2", 300)} (new base)")
  println("chr2@250 = \{gss.score_at("chr2", 250)} (previous base kept)")
}
