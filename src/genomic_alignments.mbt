///|
/// GenomicAlignments module, inspired by Bioconductor GenomicAlignments.

///|
pub struct GAlignments {
  seqnames : Array[String]
  starts : Array[Int]
  ends : Array[Int]
  strands : Array[Strand]
  widths : Array[Int]
  qnames : Array[String]
  mapqs : Array[Int]
  flags : Array[Int]
}

///|
pub fn GAlignments::new(
  seqnames : Array[String],
  starts : Array[Int],
  ends : Array[Int],
  strands : Array[Strand],
  qnames : Array[String],
  mapqs : Array[Int],
  flags : Array[Int],
) -> GAlignments {
  let n = seqnames.length()
  let widths : Array[Int] = Array::make(n, 0)
  let mut i = 0
  while i < n {
    widths[i] = ends[i] - starts[i] + 1
    i = i + 1
  }
  GAlignments::{ seqnames, starts, ends, strands, widths, qnames, mapqs, flags }
}

///|
pub fn GAlignments::length(self : GAlignments) -> Int {
  self.seqnames.length()
}

///|
pub fn GAlignments::filter_by_mapq(
  self : GAlignments,
  min_mapq : Int,
) -> GAlignments {
  let seqnames : Array[String] = Array::new()
  let starts : Array[Int] = Array::new()
  let ends : Array[Int] = Array::new()
  let strands : Array[Strand] = Array::new()
  let qnames : Array[String] = Array::new()
  let mapqs : Array[Int] = Array::new()
  let flags : Array[Int] = Array::new()

  let mut i = 0
  while i < self.length() {
    if self.mapqs[i] >= min_mapq {
      seqnames.push(self.seqnames[i])
      starts.push(self.starts[i])
      ends.push(self.ends[i])
      strands.push(self.strands[i])
      qnames.push(self.qnames[i])
      mapqs.push(self.mapqs[i])
      flags.push(self.flags[i])
    }
    i = i + 1
  }

  GAlignments::new(seqnames, starts, ends, strands, qnames, mapqs, flags)
}

///|
pub fn GAlignments::filter_by_strand(
  self : GAlignments,
  strand : Strand,
) -> GAlignments {
  let seqnames : Array[String] = Array::new()
  let starts : Array[Int] = Array::new()
  let ends : Array[Int] = Array::new()
  let strands : Array[Strand] = Array::new()
  let qnames : Array[String] = Array::new()
  let mapqs : Array[Int] = Array::new()
  let flags : Array[Int] = Array::new()

  let mut i = 0
  while i < self.length() {
    if self.strands[i] == strand {
      seqnames.push(self.seqnames[i])
      starts.push(self.starts[i])
      ends.push(self.ends[i])
      strands.push(self.strands[i])
      qnames.push(self.qnames[i])
      mapqs.push(self.mapqs[i])
      flags.push(self.flags[i])
    }
    i = i + 1
  }

  GAlignments::new(seqnames, starts, ends, strands, qnames, mapqs, flags)
}

///|
pub fn GAlignments::to_granges(self : GAlignments) -> GRanges {
  let n = self.length()
  let ranges : Array[(Int, Int)] = Array::make(n, (0, 0))
  let mut i = 0
  while i < n {
    ranges[i] = (self.starts[i], self.ends[i])
    i = i + 1
  }
  granges(self.seqnames, ranges, self.strands)
}

///|
pub struct Coverage {
  seqname : String
  positions : Array[Int]
  depths : Array[Int]
}

///|
pub fn coverage(galn : GAlignments) -> Array[Coverage] {
  let seqnames : Array[String] = Array::new()
  let all_positions : Array[Array[Int]] = Array::new()
  let all_depths : Array[Array[Int]] = Array::new()

  let mut i = 0
  while i < galn.length() {
    let seqname = galn.seqnames[i]
    let start = galn.starts[i]
    let end = galn.ends[i]

    let mut seq_idx = -1
    let mut j = 0
    while j < seqnames.length() {
      if seqnames[j] == seqname {
        seq_idx = j
      }
      j = j + 1
    }

    if seq_idx == -1 {
      seqnames.push(seqname)
      all_positions.push(Array::new())
      all_depths.push(Array::new())
      seq_idx = seqnames.length() - 1
    }

    let mut pos = start
    while pos <= end {
      let mut found = false
      let mut k = 0
      while k < all_positions[seq_idx].length() {
        if all_positions[seq_idx][k] == pos {
          all_depths[seq_idx][k] = all_depths[seq_idx][k] + 1
          found = true
          k = all_positions[seq_idx].length()
        }
        k = k + 1
      }
      if !found {
        all_positions[seq_idx].push(pos)
        all_depths[seq_idx].push(1)
      }
      pos = pos + 1
    }
    i = i + 1
  }

  let result : Array[Coverage] = Array::new()
  let mut l = 0
  while l < seqnames.length() {
    let n_pos = all_positions[l].length()
    let mut m = 0
    while m < n_pos - 1 {
      let mut n = m + 1
      while n < n_pos {
        if all_positions[l][m] > all_positions[l][n] {
          let temp = all_positions[l][m]
          all_positions[l][m] = all_positions[l][n]
          all_positions[l][n] = temp
          let temp_d = all_depths[l][m]
          all_depths[l][m] = all_depths[l][n]
          all_depths[l][n] = temp_d
        }
        n = n + 1
      }
      m = m + 1
    }
    result.push(Coverage::{
      seqname: seqnames[l],
      positions: all_positions[l],
      depths: all_depths[l],
    })
    l = l + 1
  }

  result
}

///|
pub struct OverlapSummary {
  feature_id : String
  mut count : Int
  mut sum_width : Int
  mut mean_width : Double
  mut coverage_fraction : Double
}

///|
pub fn OverlapSummary::new(
  feature_id : String,
  count : Int,
  sum_width : Int,
  mean_width : Double,
  coverage_fraction : Double,
) -> OverlapSummary {
  OverlapSummary::{
    feature_id,
    count,
    sum_width,
    mean_width,
    coverage_fraction,
  }
}

///|
fn ranges_overlap(start1 : Int, end1 : Int, start2 : Int, end2 : Int) -> Bool {
  start1 <= end2 && start2 <= end1
}

///|
pub fn summarize_overlaps(
  galn : GAlignments,
  features : GRanges,
  feature_ids : Array[String],
) -> Array[OverlapSummary] {
  let n_features = granges_length(features)
  let result : Array[OverlapSummary] = Array::make(
    n_features,
    OverlapSummary::new("", 0, 0, 0.0, 0.0),
  )

  let mut i = 0
  while i < n_features {
    result[i] = OverlapSummary::new(feature_ids[i], 0, 0, 0.0, 0.0)
    i = i + 1
  }

  let mut j = 0
  while j < galn.length() {
    let aln_seq = galn.seqnames[j]
    let aln_start = galn.starts[j]
    let aln_end = galn.ends[j]
    let aln_strand = galn.strands[j]

    let mut k = 0
    while k < n_features {
      if features.seqnames[k] == aln_seq &&
        ranges_overlap(aln_start, aln_end, features.starts[k], features.ends[k]) &&
        (
          features.strands[k] == Strand::Star ||
          features.strands[k] == aln_strand
        ) {
        result[k].count = result[k].count + 1
        let overlap_start = if aln_start > features.starts[k] {
          aln_start
        } else {
          features.starts[k]
        }
        let overlap_end = if aln_end < features.ends[k] {
          aln_end
        } else {
          features.ends[k]
        }
        result[k].sum_width = result[k].sum_width +
          (overlap_end - overlap_start + 1)
      }
      k = k + 1
    }
    j = j + 1
  }

  let mut l = 0
  while l < n_features {
    if result[l].count > 0 {
      result[l].mean_width = result[l].sum_width.to_double() /
        result[l].count.to_double()
      result[l].coverage_fraction = result[l].sum_width.to_double() /
        features.widths[l].to_double()
    }
    l = l + 1
  }

  result
}

///|
pub struct PileupElement {
  seqname : String
  position : Int
  mut depth : Int
  ref_base : String
  alts : Array[String]
  alt_counts : Array[Int]
}

///|
pub fn pileup(
  galn : GAlignments,
  reference : @hashmap.HashMap[String, String],
) -> Array[PileupElement] {
  let result : Array[PileupElement] = Array::new()

  let mut i = 0
  while i < galn.length() {
    let seqname = galn.seqnames[i]
    let start = galn.starts[i]
    let seq = galn.qnames[i]
    let ref_seq = match reference.get(seqname) {
      Some(s) => s
      None => ""
    }

    let mut pos = start
    let mut seq_pos = 0
    while pos <= galn.ends[i] && seq_pos < seq.length() {
      let base = seq[seq_pos:seq_pos + 1].to_string()

      let mut found = false
      let mut j = 0
      while j < result.length() {
        if result[j].seqname == seqname && result[j].position == pos {
          let mut k = 0
          while k < result[j].alts.length() {
            if result[j].alts[k] == base && base != result[j].ref_base {
              result[j].alt_counts[k] = result[j].alt_counts[k] + 1
              found = true
              k = result[j].alts.length()
            }
            k = k + 1
          }
          if !found && base != result[j].ref_base {
            result[j].alts.push(base)
            result[j].alt_counts.push(1)
            found = true
          }
          result[j].depth = result[j].depth + 1
          found = true
          j = result.length()
        }
        j = j + 1
      }

      if !found {
        let ref_base = if pos > 0 && pos <= ref_seq.length() {
          ref_seq[pos - 1:pos].to_string()
        } else {
          "N"
        }
        let alts = if base != ref_base { [base] } else { [] }
        let alt_counts = if base != ref_base { [1] } else { [] }
        result.push(PileupElement::{
          seqname,
          position: pos,
          depth: 1,
          ref_base,
          alts,
          alt_counts,
        })
      }

      pos = pos + 1
      seq_pos = seq_pos + 1
    }
    i = i + 1
  }

  let n = result.length()
  let mut k = 0
  while k < n - 1 {
    let mut l = k + 1
    while l < n {
      if result[k].seqname > result[l].seqname ||
        (
          result[k].seqname == result[l].seqname &&
          result[k].position > result[l].position
        ) {
        let temp = result[k]
        result[k] = result[l]
        result[l] = temp
      }
      l = l + 1
    }
    k = k + 1
  }

  result
}
