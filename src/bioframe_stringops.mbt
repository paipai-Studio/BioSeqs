/// bioframe_stringops.mbt — String utilities for genomic regions.
///
/// Bioframe equivalent: `bioframe.parse_region_string`,
/// `bioframe.to_ucsc_string`, `bioframe.is_complete_ucsc_string`,
/// and the natural-chromosome sort helpers.

// ================================================================
// 1. to_ucsc_string
// ================================================================

///|
/// Format a genomic coordinate as a UCSC-style region string.
///
/// Input coordinates are **0-based half-open** (BioSeqs internal
/// convention, matches BAM/BED).  The output uses **1-based inclusive**
/// ranges (the UCSC convention).  Conversion rules:
///
///   - start is bumped up by one: 0-based s → 1-based s+1
///   - end stays the same:         0-based half-open `[s, e)` covers
///                                 positions s..e-1 (0-based), which are
///                                 1-based positions s+1..e
///
///     bf_to_ucsc_string("chr1", 0, 10)   // "chr1:1-10"
///     bf_to_ucsc_string("chr1", 100, 200)  // "chr1:101-200"
///
pub fn bf_to_ucsc_string(chrom : String, start : Int32, end : Int32) -> String {
  let s = start + 1
  let e = end
  chrom + ":" + s.to_string() + "-" + e.to_string()
}

///|
/// Version that takes end position directly — also 0-based half-open.
pub fn bf_to_ucsc_region(iv : BedInterval) -> String {
  bf_to_ucsc_string(iv.chrom, iv.start, iv.end)
}

// ================================================================
// 2. is_complete_ucsc_string
// ================================================================

///|
/// True if `s` matches the UCSC region pattern `chrom:start-end`.
///
/// Checks:
///   - exactly one colon ':'
///   - exactly one dash '-' after the colon
///   - start and end are both non-negative integers with no leading zeros
///     (optional for small numbers)
///   - chrom part is non-empty
///
/// Note: MoonBit has no regex stdlib, so this is a hand-written scanner.
/// It doesn't try to be perfect — the goal is to reject obvious junk so
/// that callers can safely `parse_region` / `bf_parse_region_string`.
pub fn bf_is_complete_ucsc_string(s : String) -> Bool {
  let n = s.length()
  if n < 5 {
    return false // minimum: "c:1-2"
  }
  // Find colon
  let mut colon = -1
  let mut i = 0
  while i < n {
    if s.unsafe_get(i).to_int() == ':'.to_int() {
      colon = i
      break
    }
    i = i + 1
  }
  if colon <= 0 || colon >= n - 3 {
    return false // chrom empty, or nothing after ':'
  }
  // Find dash AFTER colon
  let mut dash = -1
  i = colon + 1
  while i < n {
    if s.unsafe_get(i).to_int() == '-'.to_int() {
      dash = i
      break
    }
    i = i + 1
  }
  if dash < colon + 1 || dash >= n - 1 {
    return false
  }
  // Make sure there's no *second* colon or dash after the first colon
  i = colon + 1
  while i < n {
    if i != dash {
      let c = s.unsafe_get(i).to_int()
      if c == ':'.to_int() {
        return false
      }
      // digits only after colon (dash is the separator we already found)
      if c != '-'.to_int() && (c < '0'.to_int() || c > '9'.to_int()) {
        return false
      }
    }
    i = i + 1
  }
  true
}

// ================================================================
// 3. bf_parse_region_string — enhancement over bam.mbt parse_region
// ================================================================

///|
/// Parse a UCSC-style region string into (canonical chrom, 0-based start,
/// 0-based end).  Mirrors `bam.mbt::parse_region` but adds assembly
/// awareness: when `assembly` is passed, alias resolution is applied.
///
/// Supported formats:
///   "chr1"            → ("chr1", None, None)
///   "chr1:100"        → ("chr1", 99, 100)
///   "chr1:100-200"    → ("chr1", 99, 200)
///   "chr1:100..200"   → ("chr1", 99, 200)
///   "chr1:100+50"     → ("chr1", 99, 149)   (length, 1-based)
///
/// With alias resolution (assembly ~= None):
///   "MT:1-100"        with assembly="hg38" → ("chrM", 0, 100)
///   "1:100-200"       with assembly="hg38" → ("chr1", 99, 200)
///   "X"               with assembly="hg38" → ("chrX", None, None)
///
pub fn bf_parse_region_string(
  region : String,
  assembly~ : String? = None,
) -> (String, Int32?, Int32?) {
  // reuse existing parse_region logic for the interval part
  // but handle alias resolution on the chrom part
  let (contig_raw, s, e) = parse_region(region)

  // Apply alias resolution if an assembly is supplied
  let contig = match assembly {
    None => contig_raw
    Some(asm_name) =>
      match bf_resolve_alias(asm_name, contig_raw) {
        Some(c) => c
        None => contig_raw
      }
  }
  (contig, s, e)
}

// ================================================================
// 4. Natural sort (arrownatsort equivalent)
// ================================================================

///|
/// Chromosome-name comparison key for natural sort.
///
/// Splits a chromosome string into a text prefix and a numeric suffix,
/// then compares lexicographically first on the prefix, numerically on
/// the suffix.  This produces the familiar biological order:
///
///     chr1, chr2, chr3, ..., chr10, chr11, ..., chr22, chrX, chrY, chrM
///
/// instead of a purely lexicographic
///
///     chr1, chr10, chr11, ..., chr2, chr20, chr21, ...
///
pub struct BfNatKey {
  text : String
  /// -1 means "no numeric suffix"
  num : Int64
}

///|
/// Extract the natural-sort key from a chromosome name.
///
/// Rules:
///   - leading text (e.g. "chr") is taken verbatim
///   - trailing digits form the numeric suffix
///   - if no trailing digits, num stays -1 and comparison falls back
///     to raw string
///   - special case: "chrX" / "chrY" / "chrM" → no numeric suffix,
///     so they sort alphabetically after the numbered autosomes
pub fn bf_nat_key(name : String) -> BfNatKey {
  let mut num_start = name.length()
  let mut j = name.length() - 1
  while j >= 0 {
    let c = name.unsafe_get(j).to_int()
    if c >= '0'.to_int() && c <= '9'.to_int() {
      num_start = j
    } else {
      break
    }
    j = j - 1
  }
  if num_start == name.length() {
    // No numeric suffix
    BfNatKey::{ text: name.to_lower(), num: -1, }
  } else {
    let prefix = name[0:num_start].to_owned().to_lower()
    let digits = name[num_start:name.length()].to_owned()
    let n = parse_int(digits).to_int64()
    BfNatKey::{ text: prefix, num: n, }
  }
}

///|
/// Compare two natural keys, useful as a sort_by closure.
pub fn bf_nat_cmp(a : BfNatKey, b : BfNatKey) -> Int32 {
  if a.text < b.text {
    -1
  } else if a.text > b.text {
    1
  } else if a.num < b.num {
    -1
  } else if a.num > b.num {
    1
  } else {
    0
  }
}

///|
/// Sort an array of chromosome names using natural ordering.
/// Returns a *new* sorted array; input is not mutated.
pub fn bf_natsort_chroms(names : Array[String]) -> Array[String] {
  let out = names.copy()
  out.sort_by(fn(a, b) { bf_nat_cmp(bf_nat_key(a), bf_nat_key(b)) })
  out
}

///|
/// Sort BedIntervals primarily by natural chromosome order, secondarily
/// by start coordinate (then end).  Returns a new array.
pub fn bf_natsort_intervals(
  intervals : Array[BedInterval],
) -> Array[BedInterval] {
  let out = intervals.copy()
  out.sort_by(fn(a, b) {
    let ka = bf_nat_key(a.chrom)
    let kb = bf_nat_key(b.chrom)
    let c = bf_nat_cmp(ka, kb)
    if c != 0 {
      return c
    }
    if a.start < b.start {
      return -1
    } else if a.start > b.start {
      return 1
    }
    if a.end < b.end {
      return -1
    } else if a.end > b.end {
      return 1
    }
    0
  })
  out
}
