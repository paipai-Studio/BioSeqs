///|
/// ORFik-style ORF detection and classification.

pub(all) enum OrffFrame {
  OrffFrameF0
  OrffFrameF1
  OrffFrameF2
  OrffFrameR0
  OrffFrameR1
  OrffFrameR2
}

pub(all) enum OrffClass {
  OrffCDS
  OrffUTR5
  OrffUTR3
  OrffUORF
  OrffDORF
  OrffOverlapCDS
  OrffIntergenic
}

pub(all) struct OrffORF {
  start : Int
  stop : Int
  frame : OrffFrame
  length : Int
  sequence : String
  protein : String
  kozak_score : Double
  orf_class : OrffClass
}

pub fn orff_frame_label(f : OrffFrame) -> String {
  match f {
    OrffFrameF0 => "F0"
    OrffFrameF1 => "F1"
    OrffFrameF2 => "F2"
    OrffFrameR0 => "R0"
    OrffFrameR1 => "R1"
    OrffFrameR2 => "R2"
  }
}

pub fn orff_class_label(c : OrffClass) -> String {
  match c {
    OrffCDS => "CDS"
    OrffUTR5 => "5UTR"
    OrffUTR3 => "3UTR"
    OrffUORF => "uORF"
    OrffDORF => "dORF"
    OrffOverlapCDS => "overlapCDS"
    OrffIntergenic => "intergenic"
  }
}

fn orff_is_reverse(f : OrffFrame) -> Bool {
  match f {
    OrffFrameR0 | OrffFrameR1 | OrffFrameR2 => true
    _ => false
  }
}

fn orff_frame_offset(f : OrffFrame) -> Int {
  match f {
    OrffFrameF0 => 0
    OrffFrameF1 => 1
    OrffFrameF2 => 2
    OrffFrameR0 => 0
    OrffFrameR1 => 1
    OrffFrameR2 => 2
  }
}

fn orff_revcomp(seq : String) -> String {
  let n = seq.length()
  let mut result = ""
  let mut i = n - 1
  while i >= 0 {
    let c = seq.unsafe_get(i).unsafe_to_char()
    let comp = match c {
      'A' => 'T'
      'T' | 'U' => 'A'
      'G' => 'C'
      'C' => 'G'
      'N' => 'N'
      _ => 'N'
    }
    result = result + comp.to_string()
    i = i - 1
  }
  result
}

pub fn orff_translate(seq : String) -> String {
  let n = seq.length()
  let mut protein = ""
  let mut i = 0
  while i + 2 < n {
    let codon = substring(seq, i, 3)
    let aa = match codon {
      "TTT" | "TTC" => "F"
      "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => "L"
      "ATT" | "ATC" | "ATA" => "I"
      "ATG" => "M"
      "GTT" | "GTC" | "GTA" | "GTG" => "V"
      "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => "S"
      "CCT" | "CCC" | "CCA" | "CCG" => "P"
      "ACT" | "ACC" | "ACA" | "ACG" => "T"
      "GCT" | "GCC" | "GCA" | "GCG" => "A"
      "TAT" | "TAC" => "Y"
      "TAA" | "TAG" | "TGA" => "*"
      "CAT" | "CAC" => "H"
      "CAA" | "CAG" => "Q"
      "AAT" | "AAC" => "N"
      "AAA" | "AAG" => "K"
      "GAT" | "GAC" => "D"
      "GAA" | "GAG" => "E"
      "TGT" | "TGC" => "C"
      "TGG" => "W"
      "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => "R"
      "GGT" | "GGC" | "GGA" | "GGG" => "G"
      _ => "X"
    }
    if aa == "*" { break }
    protein = protein + aa
    i = i + 3
  }
  protein
}

pub fn orff_kozak_score(seq : String, start : Int) -> Double {
  let s = seq.to_upper()
  let n = s.length()
  let mut score = 0.0
  if start >= 3 {
    let p = s.unsafe_get(start - 3).unsafe_to_char()
    if p == 'A' || p == 'G' { score = score + 0.4 }
  }
  if start + 3 < n {
    let p = s.unsafe_get(start + 3).unsafe_to_char()
    if p == 'G' { score = score + 0.3 }
  }
  if start >= 2 {
    let p = s.unsafe_get(start - 2).unsafe_to_char()
    if p == 'C' || p == 'T' { score = score + 0.15 }
  }
  if start >= 1 {
    let p = s.unsafe_get(start - 1).unsafe_to_char()
    if p == 'C' || p == 'T' { score = score + 0.15 }
  }
  score
}

fn orff_scan_frame(
  seq : String, frame : OrffFrame, min_length : Int,
  start_codons : Array[String], stop_codons : Array[String],
) -> Array[OrffORF] {
  let orfs : Array[OrffORF] = Array(capacity=0)
  let n = seq.length()
  let offset = orff_frame_offset(frame)
  let is_rev = orff_is_reverse(frame)
  let work_seq = if is_rev { orff_revcomp(seq) } else { seq }
  let work_n = work_seq.length()
  let mut i = offset
  while i < work_n - 2 {
    let codon = substring(work_seq, i, 3)
    let mut is_start = false
    for sc in start_codons {
      if codon == sc { is_start = true; break }
    }
    if is_start {
      let mut j = i + 3
      let mut found_stop = false
      while j < work_n - 2 {
        let stop_codon = substring(work_seq, j, 3)
        let mut is_stop = false
        for stc in stop_codons {
          if stop_codon == stc { is_stop = true; break }
        }
        if is_stop { found_stop = true; break }
        j = j + 3
      }
      if found_stop {
        let orf_len = j - i + 3
        if orf_len >= min_length {
          let orf_seq = substring(work_seq, i, orf_len)
          let protein = orff_translate(orf_seq)
          let (orig_start, orig_stop) = if is_rev {
            (work_n - (j + 3), work_n - i)
          } else {
            (i, j + 3)
          }
          orfs.push(OrffORF::{
            start: orig_start, stop: orig_stop, frame: frame,
            length: orf_len, sequence: orf_seq, protein: protein,
            kozak_score: orff_kozak_score(work_seq, i),
            orf_class: OrffIntergenic,
          })
        }
        i = j + 3
      } else {
        i = i + 3
      }
    } else {
      i = i + 3
    }
  }
  orfs
}

pub fn orff_find_all_orfs(
  seq : String,
  min_length? : Int = 30,
  start_codons? : Array[String] = ["ATG"],
  stop_codons? : Array[String] = ["TAA", "TAG", "TGA"],
) -> Array[OrffORF] {
  let s = seq.to_upper()
  let all_orfs : Array[OrffORF] = Array(capacity=0)
  let frames = [OrffFrameF0, OrffFrameF1, OrffFrameF2,
    OrffFrameR0, OrffFrameR1, OrffFrameR2]
  for f in frames {
    let frame_orfs = orff_scan_frame(s, f, min_length, start_codons, stop_codons)
    for orf in frame_orfs { all_orfs.push(orf) }
  }
  all_orfs.sort_by(fn(a, b) { a.start.compare(b.start) })
  all_orfs
}

pub fn orff_classify_orf(orf : OrffORF, cds_start : Int, cds_stop : Int) -> OrffClass {
  let overlaps_cds = orf.start < cds_stop && orf.stop > cds_start
  if overlaps_cds {
    if orf.start == cds_start && orf.stop == cds_stop { return OrffCDS }
    return OrffOverlapCDS
  }
  if orf.stop <= cds_start { return OrffUORF }
  if orf.start >= cds_stop { return OrffDORF }
  OrffIntergenic
}

pub fn orff_classify_all(
  orfs : Array[OrffORF], cds_start : Int, cds_stop : Int,
) -> Array[OrffORF] {
  let result : Array[OrffORF] = Array(capacity=0)
  for orf in orfs {
    let cls = orff_classify_orf(orf, cds_start, cds_stop)
    result.push(OrffORF::{
      start: orf.start, stop: orf.stop, frame: orf.frame,
      length: orf.length, sequence: orf.sequence, protein: orf.protein,
      kozak_score: orf.kozak_score, orf_class: cls,
    })
  }
  result
}

pub fn orff_find_uorfs(seq : String, cds_start : Int, min_length? : Int = 30) -> Array[OrffORF] {
  let all = orff_find_all_orfs(seq, min_length=min_length)
  let uorfs : Array[OrffORF] = Array(capacity=0)
  for orf in all {
    if orf.stop <= cds_start {
      uorfs.push(OrffORF::{
        start: orf.start, stop: orf.stop, frame: orf.frame,
        length: orf.length, sequence: orf.sequence, protein: orf.protein,
        kozak_score: orf.kozak_score, orf_class: OrffUORF,
      })
    }
  }
  uorfs
}

pub fn orff_find_dorfs(seq : String, cds_stop : Int, min_length? : Int = 30) -> Array[OrffORF] {
  let all = orff_find_all_orfs(seq, min_length=min_length)
  let dorfs : Array[OrffORF] = Array(capacity=0)
  for orf in all {
    if orf.start >= cds_stop {
      dorfs.push(OrffORF::{
        start: orf.start, stop: orf.stop, frame: orf.frame,
        length: orf.length, sequence: orf.sequence, protein: orf.protein,
        kozak_score: orf.kozak_score, orf_class: OrffDORF,
      })
    }
  }
  dorfs
}

pub fn orff_overlapping_orfs(orfs : Array[OrffORF]) -> Array[(OrffORF, OrffORF)] {
  let pairs : Array[(OrffORF, OrffORF)] = Array(capacity=0)
  let n = orfs.length()
  for i in 0..<n {
    let j_start = i + 1
    for j in j_start..<n {
      let a = orfs[i]
      let b = orfs[j]
      if a.start < b.stop && b.start < a.stop {
        pairs.push((a, b))
      }
    }
  }
  pairs
}

pub fn orff_filter_by_kozak(orfs : Array[OrffORF], min_kozak : Double) -> Array[OrffORF] {
  let result : Array[OrffORF] = Array(capacity=0)
  for orf in orfs {
    if orf.kozak_score >= min_kozak { result.push(orf) }
  }
  result
}

pub fn orff_filter_by_aa_length(orfs : Array[OrffORF], min_aa : Int) -> Array[OrffORF] {
  let result : Array[OrffORF] = Array(capacity=0)
  for orf in orfs {
    if orf.protein.length() >= min_aa { result.push(orf) }
  }
  result
}

pub fn orff_export_orfs(orfs : Array[OrffORF]) -> String {
  let sb = StringBuilder()
  sb.write_string("start\tstop\tframe\tlength\tkozak_score\tclass\tprotein\n")
  for orf in orfs {
    sb.write_string(
      orf.start.to_string() + "\t" + orf.stop.to_string() + "\t" +
      orff_frame_label(orf.frame) + "\t" + orf.length.to_string() + "\t" +
      orf.kozak_score.to_string() + "\t" + orff_class_label(orf.orf_class) + "\t" +
      orf.protein + "\n"
    )
  }
  sb.to_string()
}

pub fn orff_count_by_frame(orfs : Array[OrffORF]) -> Map[String, Int] {
  let counts : Map[String, Int] = Map([], capacity=6)
  counts.set("F0", 0); counts.set("F1", 0); counts.set("F2", 0)
  counts.set("R0", 0); counts.set("R1", 0); counts.set("R2", 0)
  for orf in orfs {
    let label = orff_frame_label(orf.frame)
    counts.set(label, counts.get(label).unwrap_or(0) + 1)
  }
  counts
}

pub fn orff_count_by_class(orfs : Array[OrffORF]) -> Map[String, Int] {
  let counts : Map[String, Int] = Map([], capacity=7)
  counts.set("CDS", 0); counts.set("5UTR", 0); counts.set("3UTR", 0)
  counts.set("uORF", 0); counts.set("dORF", 0)
  counts.set("overlapCDS", 0); counts.set("intergenic", 0)
  for orf in orfs {
    let label = orff_class_label(orf.orf_class)
    counts.set(label, counts.get(label).unwrap_or(0) + 1)
  }
  counts
}
