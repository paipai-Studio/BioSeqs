///|
/// DNAcopy demo: Circular Binary Segmentation (CBS) of copy-number data,
/// ported from Bioconductor DNAcopy::segment(). Shows profile construction,
/// segmentation with default and merged options, and multi-chromosome input.
///
/// Run: moon run examples/dnacopy_demo
fn main {
  println("=== DNAcopy Demo (CBS copy-number segmentation) ===\n")

  // 1. A synthetic copy-number profile: 60 probes on chr1 with
  //    normal(0) x15 -> gain(+1) x15 -> loss(-1) x15 -> normal(0) x15,
  //    plus a little deterministic wiggle.
  println("1. Four-level synthetic profile (60 probes, 1 kb apart)")
  let vals : Array[Double] = []
  let wiggle = [0.1, -0.1, 0.05, -0.05, 0.0]
  let mut i = 0
  while i < 60 {
    let level = if i < 15 {
      0.0
    } else if i < 30 {
      1.0
    } else if i < 45 {
      -1.0
    } else {
      0.0
    }
    vals.push(level + wiggle[i % wiggle.length()])
    i = i + 1
  }
  let maploc : Array[Int] = []
  let chrom : Array[String] = Array::make(60, "chr1")
  let mut k = 0
  while k < 60 {
    maploc.push((k + 1) * 1000)
    k = k + 1
  }
  let cna = @bio.Cna::new(chrom, maploc, vals) catch { _ => abort("bad CNA") }
  let res = @bio.segment(cna, options=@bio.DnaCopyOptions::of(nperm=1000)) catch {
    _ => abort("segment failed")
  }
  println("   \{res.segments.length()} segments found:")
  for s in res.segments {
    println(
      "   \{s.chrom} \{s.start}-\{s.end} num.mark=\{s.num_mark} seg.mean=\{s.seg_mean}",
    )
  }

  // 2. sdundo: merge adjacent segments with nearly equal levels
  println("\n2. undo.splits=\"sdundo\" merges similar levels")
  // 30 probes: 0 x10, +0.2 x10, +1 x10 - with a generous undo_sd the small
  // 0.2 step is merged back into the 0 level
  let vals2 : Array[Double] = []
  let mut j = 0
  while j < 30 {
    vals2.push(if j < 10 {
      0.0
    } else if j < 20 {
      0.2
    } else {
      1.0
    })
    j = j + 1
  }
  let maploc2 : Array[Int] = []
  let mut k2 = 0
  while k2 < 30 {
    maploc2.push((k2 + 1) * 1000)
    k2 = k2 + 1
  }
  let cna2 = @bio.Cna::new(Array::make(30, "chr2"), maploc2, vals2) catch {
    _ => abort("bad CNA")
  }
  let plain = @bio.segment(cna2, options=@bio.DnaCopyOptions::of(nperm=1000)) catch {
    _ => abort("segment failed")
  }
  let merged = @bio.segment(
    cna2,
    options=@bio.DnaCopyOptions::of(nperm=1000, undo_splits="sdundo", undo_sd=10.0),
  ) catch { _ => abort("segment failed") }
  println("   without undo: \{plain.segments.length()} segments")
  for s in plain.segments {
    println("     \{s.start}-\{s.end} mean=\{s.seg_mean} (n=\{s.num_mark})")
  }
  println("   with sdundo:   \{merged.segments.length()} segments")
  for s in merged.segments {
    println("     \{s.start}-\{s.end} mean=\{s.seg_mean} (n=\{s.num_mark})")
  }

  // 3. Two chromosomes are segmented independently
  println("\n3. Multi-chromosome input")
  let chrom3 : Array[String] = []
  let maploc3 : Array[Int] = []
  let vals3 : Array[Double] = []
  // chrA: flat 0
  let mut m = 0
  while m < 20 {
    chrom3.push("chrA")
    maploc3.push((m + 1) * 1000)
    vals3.push(0.0)
    m = m + 1
  }
  // chrB: gain in the middle
  let mut m2 = 0
  while m2 < 20 {
    chrom3.push("chrB")
    maploc3.push((m2 + 1) * 1000)
    vals3.push(if m2 >= 8 && m2 < 14 { 1.0 } else { 0.0 })
    m2 = m2 + 1
  }
  let cna3 = @bio.Cna::new(chrom3, maploc3, vals3) catch { _ => abort("bad CNA") }
  let res3 = @bio.segment(cna3, options=@bio.DnaCopyOptions::of(nperm=1000)) catch {
    _ => abort("segment failed")
  }
  println("   \{res3.segments.length()} segments:")
  for s in res3.segments {
    println(
      "   \{s.chrom} \{s.start}-\{s.end} num.mark=\{s.num_mark} seg.mean=\{s.seg_mean}",
    )
  }

  println("\nDone.")
}
