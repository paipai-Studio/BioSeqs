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