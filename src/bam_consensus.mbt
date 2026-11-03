/// Consensus sequence generation from BAM pileup.
///
/// pysam equivalent: pileup iteration + consensus calling
/// (samtools consensus / pysam.AlignmentFile.pileup with majority vote)

///|
/// Consensus call result for a single position
pub struct ConsensusPosition {
  pos : Int32
  base : Char
  coverage : Int32
  consensus_score : Double
}

///|
pub fn ConsensusPosition::pos(self : ConsensusPosition) -> Int32 {
  self.pos
}

///|
pub fn ConsensusPosition::base(self : ConsensusPosition) -> Char {
  self.base
}

///|
pub fn ConsensusPosition::coverage(self : ConsensusPosition) -> Int32 {
  self.coverage
}

///|
pub fn ConsensusPosition::consensus_score(self : ConsensusPosition) -> Double {
  self.consensus_score
}

///|
/// Consensus result
pub struct ConsensusResult {
  contig : String
  sequence : String
  positions : Array[ConsensusPosition]
  mean_coverage : Double
}

///|
pub fn ConsensusResult::contig(self : ConsensusResult) -> String {
  self.contig
}

///|
pub fn ConsensusResult::sequence(self : ConsensusResult) -> String {
  self.sequence
}

///|
pub fn ConsensusResult::positions(
  self : ConsensusResult,
) -> Array[ConsensusPosition] {
  self.positions
}

///|
pub fn ConsensusResult::mean_coverage(self : ConsensusResult) -> Double {
  self.mean_coverage
}

///|
/// Generate consensus sequence from a BamFile for a given contig.
///
/// min_coverage: positions with fewer reads are output as 'N'
/// min_frac: base must appear in at least this fraction of reads to be called
/// use_qualities: weight base votes by Phred quality
pub fn bam_consensus(
  bam : BamFile,
  contig : String,
  min_coverage~ : Int32 = 1,
  min_frac~ : Double = 0.5,
  use_qualities~ : Bool = false,
) -> ConsensusResult {
  let ref_len = match bam.get_reference_length(contig) {
    Some(l) => l
    None => return { contig, sequence: "", positions: [], mean_coverage: 0.0, }
  }
  // base_counts[pos] = Map[base_char_code, weight]
  let counts : Array[Map[Int32, Double]] = Array::make(ref_len, Map([]))
  let mut total_cov = 0.0
  let mut n_positions = 0
  for r in bam.records {
    let ref_name = r.get_reference_name(bam.header)
    if ref_name != contig || r.is_unmapped() {
      continue
    }
    let seq = r.seq
    let qual = r.qual
    let mut qpos : Int32 = 0
    let mut rpos = r.pos
    for el in r.cigar {
      let oplen = el.length
      match el.op {
        BamCigarOp::Match | BamCigarOp::Equal | BamCigarOp::Diff =>
          for k = 0; k < oplen; k = k + 1 {
            if rpos >= 0 && rpos < ref_len && qpos < seq.length() {
              let base_code = seq.unsafe_get(qpos).to_int()
              let base_upper = if base_code >= 97 && base_code <= 122 {
                base_code - 32
              } else {
                base_code
              }
              if base_upper >= 65 && base_upper <= 90 {
                let weight = if use_qualities && qpos < qual.length() {
                  let q = qual.unsafe_get(qpos).to_int()
                  let phred = q - 33
                  @math.pow(10.0, phred.to_double() / -10.0)
                } else {
                  1.0
                }
                let cur = match counts[rpos].get(base_upper) {
                  Some(v) => v
                  None => 0.0
                }
                counts[rpos][base_upper] = cur + weight
              }
            }
            qpos = qpos + 1
            rpos = rpos + 1
          }
        BamCigarOp::Insert | BamCigarOp::SoftClip => qpos = qpos + oplen
        BamCigarOp::Delete | BamCigarOp::Skip => rpos = rpos + oplen
        BamCigarOp::HardClip | BamCigarOp::Pad => ()
      }
    }
  }
  // Build consensus
  let buf = StringBuilder()
  let positions : Array[ConsensusPosition] = []
  for pi = 0; pi < ref_len; pi = pi + 1 {
    let total = counts[pi]
      .values()
      .to_array()
      .fold(init=0.0, fn(acc : Double, v : Double) -> Double { acc + v })
    let cov_int = counts[pi].length()
    if cov_int < min_coverage || total == 0.0 {
      buf.write_string("N")
      positions.push({
        pos: pi,
        base: 'N',
        coverage: cov_int,
        consensus_score: 0.0,
      })
    } else {
      let keys = counts[pi].keys().to_array()
      let mut best_key : Int32 = -1
      let mut best_weight = 0.0
      for ki = 0; ki < keys.length(); ki = ki + 1 {
        let w = counts[pi][keys[ki]]
        if w > best_weight {
          best_weight = w
          best_key = keys[ki]
        }
      }
      let frac = best_weight / total
      if frac >= min_frac && best_key >= 0 {
        let base_char = Int32::unsafe_to_char(best_key)
        buf.write_string(Char::to_string(base_char))
        positions.push({
          pos: pi,
          base: base_char,
          coverage: cov_int,
          consensus_score: frac,
        })
      } else {
        buf.write_string("N")
        positions.push({
          pos: pi,
          base: 'N',
          coverage: cov_int,
          consensus_score: frac,
        })
      }
    }
    total_cov = total_cov + cov_int.to_double()
    n_positions = n_positions + 1
  }
  let mean_cov = if n_positions > 0 {
    total_cov / n_positions.to_double()
  } else {
    0.0
  }
  { contig, sequence: buf.to_string(), positions, mean_coverage: mean_cov, }
}

///|
/// Generate consensus for all contigs in the BAM header.
pub fn bam_consensus_all(
  bam : BamFile,
  min_coverage~ : Int32 = 1,
  min_frac~ : Double = 0.5,
  use_qualities~ : Bool = false,
) -> Array[ConsensusResult] {
  let results : Array[ConsensusResult] = []
  let ref_names = bam.get_reference_names()
  for i = 0; i < ref_names.length(); i = i + 1 {
    results.push(
      bam_consensus(bam, ref_names[i], min_coverage~, min_frac~, use_qualities~),
    )
  }
  results
}
