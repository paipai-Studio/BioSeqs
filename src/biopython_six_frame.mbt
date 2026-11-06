/// biopython_six_frame.mbt — Six-frame translation
/// Inspired by Biopython Bio.Seq.translate with frame parameter

///|
/// Standard genetic code table (codon -> amino acid).
fn biopython_codon_to_aa(codon : String) -> String {
  match codon {
    // Phenylalanine
    "TTT" | "TTC" => "F"
    // Leucine
    "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => "L"
    // Isoleucine
    "ATT" | "ATC" | "ATA" => "I"
    // Methionine / Start
    "ATG" => "M"
    // Valine
    "GTT" | "GTC" | "GTA" | "GTG" => "V"
    // Serine
    "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => "S"
    // Proline
    "CCT" | "CCC" | "CCA" | "CCG" => "P"
    // Threonine
    "ACT" | "ACC" | "ACA" | "ACG" => "T"
    // Alanine
    "GCT" | "GCC" | "GCA" | "GCG" => "A"
    // Tyrosine
    "TAT" | "TAC" => "Y"
    // Stop
    "TAA" | "TAG" | "TGA" => "*"
    // Histidine
    "CAT" | "CAC" => "H"
    // Glutamine
    "CAA" | "CAG" => "Q"
    // Asparagine
    "AAT" | "AAC" => "N"
    // Lysine
    "AAA" | "AAG" => "K"
    // Aspartate
    "GAT" | "GAC" => "D"
    // Glutamate
    "GAA" | "GAG" => "E"
    // Cysteine
    "TGT" | "TGC" => "C"
    // Tryptophan
    "TGG" => "W"
    // Arginine
    "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => "R"
    // Glycine
    "GGT" | "GGC" | "GGA" | "GGG" => "G"
    _ => "X" // Unknown
  }
}

///|
/// DNA complement.
fn biopython_dna_complement(seq : String) -> String {
  let sb = StringBuilder()
  for i = 0; i < seq.length(); i = i + 1 {
    let ch = seq[i].to_int()
    let comp = match ch {
      65 => 84 // A->T
      84 => 65 // T->A
      71 => 67 // G->C
      67 => 71 // C->G
      97 => 116 // a->t
      116 => 97 // t->a
      103 => 99 // g->c
      99 => 103 // c->g
      _ => ch
    }
    let _ = sb.write_char(comp.unsafe_to_char())
  }
  sb.to_string()
}

///|
/// Reverse a string.
fn biopython_reverse(seq : String) -> String {
  let sb = StringBuilder()
  for i = seq.length() - 1; i >= 0; i = i - 1 {
    let _ = sb.write_char(seq[i].unsafe_to_char())
  }
  sb.to_string()
}

///|
/// Translate a single frame.
/// frame: 0, 1, 2 = forward frames starting at position 0, 1, 2
/// frame: 3, 4, 5 = reverse complement frames
fn biopython_translate_frame(sequence : String, frame : Int32) -> String {
  let seq = if frame < 3 {
    sequence[frame:]
  } else {
    let rc = biopython_dna_complement(sequence)
    biopython_reverse(rc)
  }
  let start = if frame < 3 { 0 } else { frame - 3 }
  let work = seq[start:]
  let sb = StringBuilder()
  let n = work.length()
  let mut i = 0
  while i + 3 <= n {
    let codon = work[i:i + 3].to_owned().to_upper()
    let _ = sb.write(biopython_codon_to_aa(codon))
    i = i + 3
  }
  sb.to_string()
}

///|
/// Six-frame translation: returns all 6 reading frames.
/// Frames 0-2 are forward, frames 3-5 are reverse complement.
pub fn biopython_six_frame_translation(sequence : String) -> Array[String] {
  let result : Array[String] = Array::make(6, "")
  for f = 0; f < 6; f = f + 1 {
    result[f] = biopython_translate_frame(sequence, f)
  }
  result
}

///|
/// Three-frame forward translation.
pub fn biopython_three_frame_forward(sequence : String) -> Array[String] {
  [
    biopython_translate_frame(sequence, 0),
    biopython_translate_frame(sequence, 1),
    biopython_translate_frame(sequence, 2),
  ]
}

///|
/// Three-frame reverse complement translation.
pub fn biopython_three_frame_reverse(sequence : String) -> Array[String] {
  [
    biopython_translate_frame(sequence, 3),
    biopython_translate_frame(sequence, 4),
    biopython_translate_frame(sequence, 5),
  ]
}

///|
/// Find all ORFs in a given frame.
/// Returns array of (start_pos, end_pos, protein_sequence).
pub fn biopython_find_orfs(
  sequence : String,
  frame : Int32,
  min_length~ : Int32 = 30,
) -> Array[(Int32, Int32, String)] {
  let protein = biopython_translate_frame(sequence, frame)
  let orfs : Array[(Int32, Int32, String)] = Array(capacity=0)
  let mut i = 0
  let n = protein.length()
  while i < n {
    // Find start codon (M)
    if protein[i:i + 1].to_owned() == "M" {
      let mut j = i + 1
      let mut found_stop = false
      while j < n {
        let ch = protein[j:j + 1].to_owned()
        if ch == "*" {
          found_stop = true
          break
        }
        j = j + 1
      }
      if found_stop && j - i >= min_length {
        orfs.push((i, j, protein[i:j].to_owned()))
      }
      i = j + 1
    } else {
      i = i + 1
    }
  }
  orfs
}
