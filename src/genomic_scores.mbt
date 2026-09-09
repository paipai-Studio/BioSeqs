///|
/// GenomicScores - position-specific genomic scores.
/// Inspired by the Bioconductor GenomicScores package (Puigdevall & Ramos,
/// 2018). Stores per-chromosome score series (one per score population,
/// e.g. phastCons/phyloP), retrieves scores at single positions with a
/// step (constant) or linear interpolation, and summarises scores over
/// multi-base genomic ranges.

///|
/// Interpolation mode for positions between stored points.
pub(all) enum GsInterpolation {
  /// Score of the nearest stored position at or before the query (step
  /// function), matching the per-position storage of conservation scores.
  GsConstant
  /// Linear interpolation between the bracketing stored positions.
  GsLinear
}

///|
/// Per-chromosome score series. Positions are sorted ascending; each
/// population holds one score per stored position.
pub struct GsChromosomeData {
  mut positions : Array[Int]
  populations : Map[String, Array[Double]]
}

///|
/// A genomic scores container, mirroring the GScores class.
pub struct GScores {
  provider : String
  organism : String
  mut default_population : String
  mut interpolation : GsInterpolation
  chromosomes : Map[String, GsChromosomeData]
}

///|
/// Create an empty scores container for the given provider/organism.
pub fn GScores::new(provider : String, organism : String) -> GScores {
  GScores::{
    provider,
    organism,
    default_population: "default",
    interpolation: GsInterpolation::GsConstant,
    chromosomes: Map([], capacity=0),
  }
}

///|
/// Provider name.
pub fn GScores::provider_name(self : GScores) -> String {
  self.provider
}

///|
/// Organism name.
pub fn GScores::organism_name(self : GScores) -> String {
  self.organism
}

///|
/// Install one score population on a chromosome. Positions must be
/// sorted ascending; unsorted input is sorted together with the scores.
pub fn GScores::set_scores(
  self : GScores,
  seqname : String,
  population : String,
  positions : Array[Int],
  scores : Array[Double],
) -> Unit raise SeqIOError {
  if positions.length() != scores.length() {
    raise SeqIOError("positions and scores length mismatch")
  }
  let chrom_data = match self.chromosomes.get(seqname) {
    Some(existing) => existing
    None => {
      let fresh = GsChromosomeData::{ positions: [], populations: Map([], capacity=0) }
      self.chromosomes.set(seqname, fresh)
      fresh
    }
  }
  let n = positions.length()
  let order : Array[Int] = Array::new(capacity=n)
  for i in 0..<n {
    order.push(i)
  }
  // simple stable insertion sort by position
  let mut i = 1
  while i < n {
    let key = order[i]
    let mut j = i - 1
    while j >= 0 && positions[order[j]] > positions[key] {
      order[j + 1] = order[j]
      j = j - 1
    }
    order[j + 1] = key
    i = i + 1
  }
  let sorted_positions : Array[Int] = Array::new(capacity=n)
  let sorted_scores : Array[Double] = Array::new(capacity=n)
  for idx in order {
    sorted_positions.push(positions[idx])
    sorted_scores.push(scores[idx])
  }
  if chrom_data.positions.length() == 0 {
    chrom_data.positions = sorted_positions
    // align the new population to the freshly installed positions
    let aligned = gs_align_series(
      sorted_positions,
      sorted_scores,
      chrom_data.positions,
      chrom_data.positions,
    )
    chrom_data.populations.set(population, aligned)
  } else {
    // merge new positions with existing ones (they may be distinct)
    let merged : Array[Int] = []
    let mut a = 0
    let mut b = 0
    let old = chrom_data.positions
    while a < old.length() || b < n {
      if a >= old.length() || (b < n && sorted_positions[b] < old[a]) {
        merged.push(sorted_positions[b])
        b = b + 1
      } else if b >= n || old[a] < sorted_positions[b] {
        merged.push(old[a])
        a = a + 1
      } else {
        merged.push(old[a])
        a = a + 1
        b = b + 1
      }
    }
    chrom_data.positions = merged
    // re-align every existing population to merged positions (missing ->
    // NaN) and fill gaps at the newly added positions with the new scores
    let fresh = gs_align_series(sorted_positions, sorted_scores, merged, merged)
    let exists = match chrom_data.populations.get(population) {
      Some(_) => true
      None => false
    }
    for _name, series in chrom_data.populations {
      let realigned = gs_align_series(old, series, merged, merged)
      if _name == population {
        for i in 0..<merged.length() {
          if realigned[i].is_nan() && !fresh[i].is_nan() {
            realigned[i] = fresh[i]
          }
        }
      }
      series.clear()
      for v in realigned {
        series.push(v)
      }
    }
    if !exists {
      chrom_data.populations.set(population, fresh)
    }
  }
}

///|
/// Align a score series defined on `src_positions` onto `dst_positions`;
/// entries with no source position become NaN.
fn gs_align_series(
  src_positions : Array[Int],
  src_scores : Array[Double],
  dst_positions : Array[Int],
  _dst_mirror : Array[Int],
) -> Array[Double] {
  let _ = _dst_mirror
  let out : Array[Double] = Array::new(capacity=dst_positions.length())
  let mut s = 0
  for d in dst_positions {
    let mut value = 0.0 / 0.0
    while s < src_positions.length() && src_positions[s] < d {
      s = s + 1
    }
    if s < src_positions.length() && src_positions[s] == d {
      value = src_scores[s]
    }
    out.push(value)
  }
  out
}

///|
/// Available score population names for a chromosome (union across
/// chromosomes).
pub fn GScores::populations(self : GScores) -> Array[String] {
  let names : Array[String] = []
  for _seq, chrom_data in self.chromosomes {
    for pop, _series in chrom_data.populations {
      let mut found = false
      for existing in names {
        if existing == pop {
          found = true
          break
        }
      }
      if !found {
        names.push(pop)
      }
    }
  }
  names
}

///|
/// Select the population used when none is specified.
pub fn GScores::set_default_population(self : GScores, population : String) -> Unit {
  self.default_population = population
}

///|
/// Select the interpolation mode for queries between stored positions.
pub fn GScores::set_interpolation(self : GScores, mode : GsInterpolation) -> Unit {
  self.interpolation = mode
}

///|
/// Score at a single genomic position. Returns NaN when the position
/// precedes the first stored position or the chromosome is unknown.
pub fn GScores::score_at(
  self : GScores,
  seqname : String,
  position : Int,
  population? : String = "",
) -> Double {
  let pop = if population.length() > 0 { population } else { self.default_population }
  match self.chromosomes.get(seqname) {
    Some(chrom_data) =>
      match chrom_data.populations.get(pop) {
        Some(series) => {
          // index of the last stored position <= query (or -1)
          let mut lo = 0
          let mut hi = chrom_data.positions.length() - 1
          let mut idx = -1
          while lo <= hi {
            let mid = (lo + hi) / 2
            if chrom_data.positions[mid] <= position {
              idx = mid
              lo = mid + 1
            } else {
              hi = mid - 1
            }
          }
          if idx < 0 {
            return 0.0 / 0.0
          }
          if self.interpolation is GsConstant || idx + 1 >= chrom_data.positions.length() {
            return series[idx]
          }
          let p1 = chrom_data.positions[idx]
          let p2 = chrom_data.positions[idx + 1]
          let s1 = series[idx]
          let s2 = series[idx + 1]
          if p2 == p1 {
            return s1
          }
          s1 + (s2 - s1) * ((position - p1).to_double() / (p2 - p1).to_double())
        }
        None => 0.0 / 0.0
      }
    None => 0.0 / 0.0
  }
}

///|
/// Summarise scores across the closed [start, end] range of single-base
/// positions. Supported summaries: "mean" (default), "min", "max",
/// "median". Only stored positions inside the range contribute; when the
/// range holds no stored position the step/linear interpolated scores at
/// both endpoints are used, and a range without any data yields NaN.
pub fn GScores::range_score(
  self : GScores,
  seqname : String,
  start : Int,
  end : Int,
  summary? : String = "mean",
  population? : String = "",
) -> Double {
  let pop = if population.length() > 0 { population } else { self.default_population }
  let values : Array[Double] = []
  match self.chromosomes.get(seqname) {
    Some(chrom_data) =>
      match chrom_data.populations.get(pop) {
        Some(series) => {
          // lower bound: first stored position >= start
          let mut lo = 0
          let mut hi = chrom_data.positions.length()
          while lo < hi {
            let mid = (lo + hi) / 2
            if chrom_data.positions[mid] < start {
              lo = mid + 1
            } else {
              hi = mid
            }
          }
          let mut idx = lo
          while idx < chrom_data.positions.length() && chrom_data.positions[idx] <= end {
            let v = series[idx]
            if !v.is_nan() {
              values.push(v)
            }
            idx = idx + 1
          }
          if values.length() == 0 {
            // fall back to the interpolated endpoint scores
            let v1 = self.score_at(seqname, start, population=pop)
            if !v1.is_nan() {
              values.push(v1)
            }
            if end > start {
              let v2 = self.score_at(seqname, end, population=pop)
              if !v2.is_nan() {
                values.push(v2)
              }
            }
          }
        }
        None => ()
      }
    None => ()
  }
  if values.length() == 0 {
    return 0.0 / 0.0
  }
  match summary {
    "min" => {
      let mut best = values[0]
      for v in values {
        if v < best {
          best = v
        }
      }
      best
    }
    "max" => {
      let mut best = values[0]
      for v in values {
        if v > best {
          best = v
        }
      }
      best
    }
    "median" => gs_median(values)
    _ => gs_mean(values)
  }
}

///|
/// Fetch scores for a batch of ranges given parallel arrays (GRanges-like).
pub fn gscores_fetch(
  gss : GScores,
  seqnames : Array[String],
  starts : Array[Int],
  ends : Array[Int],
  summary? : String = "mean",
  population? : String = "",
) -> Array[Double] raise SeqIOError {
  if seqnames.length() != starts.length() || starts.length() != ends.length() {
    raise SeqIOError("gscores_fetch: seqnames/starts/ends length mismatch")
  }
  let out : Array[Double] = Array::new(capacity=seqnames.length())
  for i in 0..<seqnames.length() {
    out.push(gss.range_score(seqnames[i], starts[i], ends[i], summary~, population~))
  }
  out
}

///|
fn gs_mean(values : Array[Double]) -> Double {
  let mut total = 0.0
  for v in values {
    total = total + v
  }
  total / values.length().to_double()
}

///|
fn gs_median(values : Array[Double]) -> Double {
  let sorted = values.copy()
  // insertion sort (score series are small in tests/examples)
  let mut i = 1
  while i < sorted.length() {
    let key = sorted[i]
    let mut j = i - 1
    while j >= 0 && sorted[j] > key {
      sorted[j + 1] = sorted[j]
      j = j - 1
    }
    sorted[j + 1] = key
    i = i + 1
  }
  let n = sorted.length()
  if n % 2 == 1 {
    sorted[n / 2]
  } else {
    0.5 * (sorted[n / 2 - 1] + sorted[n / 2])
  }
}

///|
/// Number of stored positions on a chromosome (0 if unknown).
pub fn GScores::n_stored_positions(self : GScores, seqname : String) -> Int {
  match self.chromosomes.get(seqname) {
    Some(chrom_data) => chrom_data.positions.length()
    None => 0
  }
}
