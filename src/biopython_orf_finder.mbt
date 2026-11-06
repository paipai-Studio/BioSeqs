/// biopython_orf_finder.mbt — Open Reading Frame finder
/// Inspired by Biopython Bio.SeqUtils: ORF detection in all 6 frames

///|
/// An ORF (Open Reading Frame) record.
pub struct BiopythonORF {
  start : Int32 // start position in original sequence (0-based)
  end : Int32 // end position (exclusive)
  frame : Int32 // frame 1-6 (1-3 forward, 4-6 reverse)
  length : Int32 // length in nucleotides (including stop codon)
  protein : String // translated protein sequence
  has_stop : Bool // whether ORF ends with a stop codon
}

///|
/// Standard codon table for translation.
fn biopython_orf_translate_codon(codon : String) -> String {
  match codon {
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
}

///|
/// Reverse complement of a DNA sequence.
fn biopython_orf_revcomp(seq : String) -> String {
  let n = seq.length()
  let result = StringBuilder()
  for i = n - 1; i >= 0; i = i - 1 {
    let ch = match seq[i].to_int() {
      65 | 97 => 'T' // A -> T
      84 | 116 => 'A' // T -> A
      71 | 103 => 'C' // G -> C
      67 | 99 => 'G' // C -> G
      _ => 'N'
    }
    let _ = result.write_char(ch)
  }
  result.to_string()
}

///|
/// Find all ORFs in a sequence.
/// Searches all 6 reading frames.
pub fn biopython_find_all_orfs(
  sequence : String,
  min_length~ : Int32 = 30, // minimum ORF length in nucleotides
  start_codon~ : String = "ATG",
) -> Array[BiopythonORF] {
  let seq = sequence.to_upper()
  let orfs : Array[BiopythonORF] = Array(capacity=0)
  // Forward frames 1-3
  for frame = 0; frame < 3; frame = frame + 1 {
    biopython_find_orfs_in_frame(
      seq,
      frame,
      frame + 1,
      min_length,
      start_codon,
      orfs,
    )
  }
  // Reverse frames 4-6
  let revcomp = biopython_orf_revcomp(seq)
  for frame = 0; frame < 3; frame = frame + 1 {
    biopython_find_orfs_in_frame(
      revcomp,
      frame,
      frame + 4,
      min_length,
      start_codon,
      orfs,
    )
  }
  // Sort by length descending
  orfs.sort_by(fn(a, b) {
    if b.length < a.length {
      -1
    } else if b.length > a.length {
      1
    } else {
      0
    }
  })
  orfs
}

///|
/// Find ORFs in a single frame.
fn biopython_find_orfs_in_frame(
  seq : String,
  frame : Int32,
  frame_num : Int32,
  min_length : Int32,
  start_codon : String,
  orfs : Array[BiopythonORF],
) -> Unit {
  let n = seq.length()
  let mut i = frame
  while i + 3 <= n {
    let codon = seq[i:i + 3].to_owned()
    if codon == start_codon {
      // Found start codon, look for stop
      let mut j = i + 3
      let mut found_stop = false
      let protein_sb = StringBuilder()
      let _ = protein_sb.write("M")
      while j + 3 <= n {
        let c = seq[j:j + 3].to_owned()
        let aa = biopython_orf_translate_codon(c)
        if aa == "*" {
          found_stop = true
          break
        }
        let _ = protein_sb.write(aa)
        j = j + 3
      }
      if found_stop {
        let orf_len = j + 3 - i
        if orf_len >= min_length {
          orfs.push({
            start: i,
            end: j + 3,
            frame: frame_num,
            length: orf_len,
            protein: protein_sb.to_string(),
            has_stop: true,
          })
        }
        i = j + 3
      } else {
        // No stop codon found, skip to next codon
        i = i + 3
      }
    } else {
      i = i + 3
    }
  }
}

///|
/// Find the longest ORF in a sequence.
pub fn biopython_longest_orf(
  sequence : String,
  min_length~ : Int32 = 30,
) -> BiopythonORF? {
  let orfs = biopython_find_all_orfs(sequence, min_length~)
  if orfs.length() == 0 {
    None
  } else {
    Some(orfs[0])
  }
}

///|
/// Count total ORFs across all frames.
pub fn biopython_count_orfs(
  sequence : String,
  min_length~ : Int32 = 30,
) -> Int32 {
  biopython_find_all_orfs(sequence, min_length~).length()
}
