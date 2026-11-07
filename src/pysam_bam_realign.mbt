/// pysam_bam_realign.mbt — Local realignment around indels
/// Inspired by pysam/GATK IndelRealigner (simplified)

///|
/// A target interval for realignment.
pub struct PysamRealignTarget {
  contig : String
  start : Int32
  end : Int32
}

///|
/// Identify realignment targets from SAM.
/// Finds regions with indels that may need realignment.
pub fn pysam_find_realign_targets(
  sam_text : String,
  window~ : Int32 = 50,
) -> Array[PysamRealignTarget] {
  let lines = sam_text.split("\n").to_array()
  let targets : Array[PysamRealignTarget] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 6 {
      continue
    }
    let rname = fields[2].to_owned()
    let pos = parse_int(fields[3].to_owned()) - 1
    let cigar = fields[5].to_owned()
    // Check for indels in CIGAR
    let mut ref_pos = pos
    let mut num = 0
    for c = 0; c < cigar.length(); c = c + 1 {
      let ch = cigar[c].to_int()
      if ch >= 48 && ch <= 57 {
        num = num * 10 + (ch - 48)
      } else {
        if ch == 73 || ch == 68 { // I or D
          // Indel at ref_pos
          let t_start = if ref_pos - window > 0 { ref_pos - window } else { 0 }
          let t_end = ref_pos + window
          targets.push({ contig: rname, start: t_start, end: t_end, })
        }
        if ch == 77 || ch == 68 || ch == 78 || ch == 61 || ch == 88 {
          ref_pos = ref_pos + num
        }
        num = 0
      }
    }
  }
  // Merge overlapping targets
  pysam_merge_targets(targets)
}

///|
/// Merge overlapping targets.
fn pysam_merge_targets(
  targets : Array[PysamRealignTarget],
) -> Array[PysamRealignTarget] {
  if targets.length() <= 1 {
    return targets
  }
  // Sort by (contig, start)
  let sorted = targets.copy()
  sorted.sort_by(fn(a, b) {
    if a.contig < b.contig {
      -1
    } else if a.contig > b.contig {
      1
    } else if a.start < b.start {
      -1
    } else if a.start > b.start {
      1
    } else {
      0
    }
  })
  let merged : Array[PysamRealignTarget] = Array(capacity=0)
  merged.push(sorted[0])
  for i = 1; i < sorted.length(); i = i + 1 {
    let last = merged[merged.length() - 1]
    let curr = sorted[i]
    if curr.contig == last.contig && curr.start <= last.end {
      // Merge
      merged[merged.length() - 1] = {
        contig: last.contig,
        start: last.start,
        end: if curr.end > last.end {
          curr.end
        } else {
          last.end
        },
      }
    } else {
      merged.push(curr)
    }
  }
  merged
}

///|
/// Check if a position is within any realignment target.
pub fn pysam_in_realign_target(
  targets : Array[PysamRealignTarget],
  contig : String,
  pos : Int32,
) -> Bool {
  for i = 0; i < targets.length(); i = i + 1 {
    if targets[i].contig == contig &&
      pos >= targets[i].start &&
      pos < targets[i].end {
      return true
    }
  }
  false
}

///|
/// Count reads overlapping realignment targets.
pub fn pysam_count_reads_in_targets(
  sam_text : String,
  targets : Array[PysamRealignTarget],
) -> Int32 {
  let lines = sam_text.split("\n").to_array()
  let mut count = 0
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 4 {
      continue
    }
    let rname = fields[2].to_owned()
    let pos = parse_int(fields[3].to_owned()) - 1
    if pysam_in_realign_target(targets, rname, pos) {
      count = count + 1
    }
  }
  count
}
