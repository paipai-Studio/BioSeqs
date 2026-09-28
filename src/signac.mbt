///|
/// Signac (single-cell ATAC-seq) core: fragment → peak×cell count matrix,
/// TSS enrichment and nucleosome-signal QC metrics.

///|
/// A single Tn5 insertion fragment (BED-like 0-based half-open + cell barcode).
pub struct SignacFragment {
  chrom : String
  start : Int
  end : Int
  cell : String
}

///|
pub fn SignacFragment::new(chrom : String, start : Int, end : Int, cell : String) -> SignacFragment {
  SignacFragment::{ chrom, start, end, cell }
}

///|
/// A TSS location (0-based position).
pub struct TssLoc {
  chrom : String
  pos : Int
}

///|
pub fn TssLoc::new(chrom : String, pos : Int) -> TssLoc {
  TssLoc::{ chrom, pos }
}

///|
/// Peak × cell integer count matrix produced by `signac_feature_matrix`.
pub struct FeatureMatrix {
  peaks : Array[BedInterval]
  cells : Array[String]
  counts : Array[Array[Int]]
}

///|
fn index_of(items : Array[String], key : String) -> Int {
  for i = 0; i < items.length(); i = i + 1 {
    if items[i] == key {
      return i
    }
  }
  -1
}

///|
/// Build a peak × cell count matrix from Tn5 fragments (`Signac::FeatureMatrix`).
/// A fragment contributes one count to every peak it overlaps (same chromosome,
/// half-open overlap) in its cell. Returns counts indexed `[peak][cell]` in the
/// order of `peaks` / `cells`.
pub fn signac_feature_matrix(
  fragments : Array[SignacFragment],
  peaks : Array[BedInterval],
  cells : Array[String],
) -> FeatureMatrix {
  let npeaks = peaks.length()
  let ncells = cells.length()
  let counts : Array[Array[Int]] = Array::make(npeaks, Array(capacity=0))
  for pi = 0; pi < npeaks; pi = pi + 1 {
    counts[pi] = Array::make(ncells, 0)
  }
  for f in fragments {
    let ci = index_of(cells, f.cell)
    if ci < 0 {
      continue
    }
    for pi = 0; pi < npeaks; pi = pi + 1 {
      let pk = peaks[pi]
      if f.chrom == pk.chrom && f.start < pk.end && pk.start < f.end {
        counts[pi][ci] = counts[pi][ci] + 1
      }
    }
  }
  FeatureMatrix::{ peaks, cells, counts }
}

///|
fn frag_overlaps_region(f : SignacFragment, chrom : String, s : Int, e : Int) -> Bool {
  f.chrom == chrom && f.start < e && s < f.end
}

///|
/// TSS enrichment: total fragments overlapping a centered window of
/// ±`tss_halfwin` bp around each TSS, divided by the average of the two
/// flanking windows (each `flank` bp wide, offset by `tss_halfwin`). This is
/// Signac's promoter-relative enrichment QC (NaN when the flanks are empty).
pub fn signac_tss_enrichment(
  fragments : Array[SignacFragment],
  tss : Array[TssLoc],
  tss_halfwin : Int,
  flank : Int,
) -> Double {
  let mut tss_count = 0
  let mut left_count = 0
  let mut right_count = 0
  for t in tss {
    for f in fragments {
      if f.chrom != t.chrom {
        continue
      }
      if frag_overlaps_region(f, t.chrom, t.pos - tss_halfwin, t.pos + tss_halfwin) {
        tss_count = tss_count + 1
      }
      if frag_overlaps_region(f, t.chrom, t.pos - tss_halfwin - flank, t.pos - tss_halfwin) {
        left_count = left_count + 1
      }
      if frag_overlaps_region(f, t.chrom, t.pos + tss_halfwin, t.pos + tss_halfwin + flank) {
        right_count = right_count + 1
      }
    }
  }
  let avg_flank = (left_count.to_double() + right_count.to_double()) / 2.0
  if avg_flank == 0.0 {
    0.0 / 0.0
  } else {
    tss_count.to_double() / avg_flank
  }
}

///|
/// Nucleosome signal: ratio of mononucleosome fragments (147–294 bp) to
/// nucleosome-free fragments (< 147 bp), the standard fragment-length
/// periodicity QC (`Signac::NucleosomeSignal`). NaN when there are no
/// nucleosome-free fragments.
pub fn signac_nucleosome_signal(fragments : Array[SignacFragment]) -> Double {
  let mut nuc_free = 0
  let mut mono = 0
  for f in fragments {
    let len = f.end - f.start
    if len < 147 {
      nuc_free = nuc_free + 1
    } else if len < 294 {
      mono = mono + 1
    }
  }
  if nuc_free == 0 {
    0.0 / 0.0
  } else {
    mono.to_double() / nuc_free.to_double()
  }
}

// ============================================================================
// Signac 深化：peak calling + gene activity
// ============================================================================

///|
/// Insertion-sort 事件列表（按位置升序，delta 不参与比较）。
fn signac_sort_events(events : Array[(Int, Int)]) -> Array[(Int, Int)] {
  let out = events.copy()
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && out[j].0 > key.0 {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

///|
/// Call peaks from Tn5 fragments by fragment-coverage thresholding (a
/// simplified `Signac::CallPeaks` / MACS-style scan). For each chromosome the
/// per-base coverage is encoded as `+1` at a fragment start and `-1` at its
/// end; contiguous regions with coverage `>= threshold` and length `>=
/// min_length` become peaks (0-based half-open).
pub fn signac_call_peaks(
  fragments : Array[SignacFragment],
  threshold : Int,
  min_length : Int,
) -> Array[BedInterval] {
  // 去重染色体
  let chroms : Array[String] = Array(capacity=0)
  for f in fragments {
    if index_of(chroms, f.chrom) < 0 {
      chroms.push(f.chrom)
    }
  }

  let result : Array[BedInterval] = Array(capacity=0)
  for chrom in chroms {
    let events : Array[(Int, Int)] = Array(capacity=0)
    for f in fragments {
      if f.chrom != chrom {
        continue
      }
      events.push((f.start, 1))
      events.push((f.end, -1))
    }
    if events.length() == 0 {
      continue
    }
    let sorted = signac_sort_events(events)
    let mut coverage = 0
    let mut peak_start = -1
    let mut i = 0
    while i < sorted.length() {
      let pos = sorted[i].0
      // 累加同位置的所有 delta
      let mut j = i
      while j < sorted.length() && sorted[j].0 == pos {
        coverage = coverage + sorted[j].1
        j = j + 1
      }
      if coverage >= threshold {
        if peak_start < 0 {
          peak_start = pos
        }
      } else {
        if peak_start >= 0 {
          if pos - peak_start >= min_length {
            result.push(BedInterval::new(chrom, peak_start, pos))
          }
          peak_start = -1
        }
      }
      i = j
    }
    // 末尾未闭合的峰（coverage 直到最后仍达标）
    if peak_start >= 0 {
      let last = sorted[sorted.length() - 1].0
      if last - peak_start >= min_length {
        result.push(BedInterval::new(chrom, peak_start, last))
      }
    }
  }
  result
}

///|
/// Gene-activity matrix produced by `signac_gene_activity` (gene × cell).
pub struct GeneActivityMatrix {
  genes : Array[String]
  cells : Array[String]
  activity : Array[Array[Int]]
}

///|
/// Aggregate a peak × cell matrix into gene activity (`Signac::GeneActivity`):
/// every peak overlapping a gene's interval extended `upstream` bp (the
/// promoter) contributes its per-cell counts to that gene. Returns activity
/// indexed `[gene][cell]`.
pub fn signac_gene_activity(
  fm : FeatureMatrix,
  genes : Array[BedInterval],
  gene_names : Array[String],
  upstream : Int,
) -> GeneActivityMatrix {
  let ngenes = genes.length()
  let ncells = fm.cells.length()
  let activity : Array[Array[Int]] = Array::make(ngenes, Array(capacity=0))
  for gi = 0; gi < ngenes; gi = gi + 1 {
    activity[gi] = Array::make(ncells, 0)
  }
  for pi = 0; pi < fm.peaks.length(); pi = pi + 1 {
    let pk = fm.peaks[pi]
    for gi = 0; gi < ngenes; gi = gi + 1 {
      let g = genes[gi]
      let promoter_start = g.start - upstream
      if pk.chrom == g.chrom && pk.start < g.end && promoter_start < pk.end {
        for ci = 0; ci < ncells; ci = ci + 1 {
          activity[gi][ci] = activity[gi][ci] + fm.counts[pi][ci]
        }
      }
    }
  }
  GeneActivityMatrix::{ genes: gene_names, cells: fm.cells, activity }
}
