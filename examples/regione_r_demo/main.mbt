///|
/// regioneR demo: permutation testing for genomic region sets, ported from
/// Bioconductor regioneR. Shows randomization strategies, association
/// statistics and permTest with empirical p-values and z-scores.
///
/// Run: moon run examples/regione_r_demo
fn main {
  println("=== regioneR Demo (region permutation testing) ===\n")

  // A tiny genome: chr1 = 10 kb, chr2 = 5 kb
  let genome = @bio.GenomeLayout::new([("chr1", 10000), ("chr2", 5000)]) catch {
    _ => abort("bad genome")
  }

  // b: 20 probes of width 51 spaced 450 bases apart on chr1
  let probe_ranges : Array[(Int, Int)] = []
  let mut k = 0
  while k < 20 {
    probe_ranges.push((250 + k * 450, 300 + k * 450))
    k = k + 1
  }
  let probes = mk("chr1", probe_ranges)

  // a: the first 5 probes -> strongly associated with b
  let a_ranges : Array[(Int, Int)] = []
  let mut j = 0
  while j < 5 {
    a_ranges.push(probe_ranges[j])
    j = j + 1
  }
  let a = mk("chr1", a_ranges)

  // 1. Association statistics
  println("1. Statistics")
  println("   numOverlaps(a, b)        = \{@bio.stat_num_overlaps(a, probes)}")
  println(
    "   numOverlaps(a, b, once)  = \{@bio.stat_num_overlaps(a, probes, count_once=true)}",
  )
  println(
    "   meanDistance(probes, a)  = \{@bio.stat_mean_distance(probes, a)}",
  )

  // 2. Randomization strategies
  println("\n2. Randomization")
  let randomized = @bio.randomize_regions(a, genome, options=@bio.RegionPermOptions::of(seed=42)) catch {
    _ => abort("randomize failed")
  }
  println("   randomizeRegions (per chromosome, mask-free):")
  let mut i = 0
  while i < @bio.granges_length(randomized) {
    println(
      "     \{randomized.seqnames[i]}:\{randomized.starts[i]}-\{randomized.ends[i]} (width \{randomized.widths[i]})",
    )
    i = i + 1
  }
  let circ = @bio.circular_randomize_regions(a, genome, options=@bio.RegionPermOptions::of(seed=42)) catch {
    _ => abort("circular randomize failed")
  }
  println("   circularRandomizeRegions keeps widths and chromosomes:")
  let mut i2 = 0
  while i2 < @bio.granges_length(circ) {
    println(
      "     \{circ.seqnames[i2]}:\{circ.starts[i2]}-\{circ.ends[i2]} (width \{circ.widths[i2]})",
    )
    i2 = i2 + 1
  }
  // a mask keeps randomized regions out of the first 9 kb
  let mask = mk("chr1", [(1, 9000)])
  let masked = @bio.randomize_regions(
    a,
    genome,
    mask=mask,
    options=@bio.RegionPermOptions::of(seed=42),
  ) catch { _ => abort("masked randomize failed") }
  println("   with mask chr1:1-9000 all regions move past 9000:")
  let mut i3 = 0
  while i3 < @bio.granges_length(masked) {
    println("     \{masked.seqnames[i3]}:\{masked.starts[i3]}-\{masked.ends[i3]}")
    i3 = i3 + 1
  }

  // 3. Overlap permutation test: a coincides with 5 of the 20 probes
  println("\n3. overlapPermTest(a, probes)")
  let res = @bio.overlap_perm_test(
    a,
    probes,
    genome,
    options=@bio.RegionPermOptions::of(ntimes=200, seed=1),
  ) catch { _ => abort("permTest failed") }
  println("   observed numOverlaps = \{res.statistic}")
  println("   permutations         = \{res.ntimes}")
  println("   z-score              = \{res.z_score}")
  println("   p-value              = \{res.p_value}")
  let mut ge = 0
  for p in res.permuted {
    if p >= res.statistic {
      ge = ge + 1
    }
  }
  println(
    "   permuted >= observed : \{ge}/\{res.permuted.length()} -> strongly associated",
  )

  // 4. A control set far from the probes is not associated
  println("\n4. Control region set (scattered probes, no association)")
  let ctrl = mk("chr1", [(9100, 9150), (9200, 9250), (9300, 9350)])
  let res2 = @bio.overlap_perm_test(
    ctrl,
    probes,
    genome,
    options=@bio.RegionPermOptions::of(ntimes=200, seed=1),
  ) catch { _ => abort("permTest failed") }
  println("   observed numOverlaps = \{res2.statistic}")
  println("   z-score              = \{res2.z_score}")
  println("   p-value              = \{res2.p_value}")

  println("\nDone.")
}

///|
/// Build a single-chromosome region set with star strands.
fn mk(seqname : String, ranges : Array[(Int, Int)]) -> @bio.GRanges {
  let n = ranges.length()
  @bio.granges(Array::make(n, seqname), ranges, Array::make(n, @bio.strand_star()))
}
