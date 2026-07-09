// Complement lookup tables, mirroring BioPython's
// Bio.Data.IUPACData.ambiguous_dna_complement / ambiguous_rna_complement.
//
// BioPython's Seq.complement() treats any U as a T (i.e. it uses the DNA
// complement table augmented with U->A), while Seq.complement_rna() treats
// any T as a U. See Bio/Seq.py lines 53-60 and 1654/1703.

///|
/// Mapping from a nucleotide character to its DNA complement.
/// Includes U -> A so that Seq.complement() works on RNA input just like
/// BioPython (which adds "U": "T"'s complement to the DNA table).
let dna_complement_map : Map[Char, Char] = {
  let m = Map::new()
  m['A'] = 'T'
  m['C'] = 'G'
  m['G'] = 'C'
  m['T'] = 'A'
  m['U'] = 'A'
  // IUPAC ambiguous codes:
  m['M'] = 'K'
  m['R'] = 'Y'
  m['W'] = 'W'
  m['S'] = 'S'
  m['Y'] = 'R'
  m['K'] = 'M'
  m['V'] = 'B'
  m['H'] = 'D'
  m['D'] = 'H'
  m['B'] = 'V'
  m['X'] = 'X'
  m['N'] = 'N'
  m
}

///|
/// Mapping from a nucleotide character to its RNA complement.
/// Includes T -> A so that Seq.complement_rna() works on DNA input just like
/// BioPython (which adds "T": "U"'s complement to the RNA table).
let rna_complement_map : Map[Char, Char] = {
  let m = Map::new()
  m['A'] = 'U'
  m['C'] = 'G'
  m['G'] = 'C'
  m['U'] = 'A'
  m['T'] = 'A'
  m['M'] = 'K'
  m['R'] = 'Y'
  m['W'] = 'W'
  m['S'] = 'S'
  m['Y'] = 'R'
  m['K'] = 'M'
  m['V'] = 'B'
  m['H'] = 'D'
  m['D'] = 'H'
  m['B'] = 'V'
  m['X'] = 'X'
  m['N'] = 'N'
  m
}

///|
/// Complement every character in `s` using `table`. Lowercase input is
/// complemented to lowercase output, matching BioPython's _maketrans.
fn complement_chars(s : String, table : Map[Char, Char]) -> String {
  let buf = StringBuilder::new()
  for c in s.iter() {
    if c.is_ascii_lowercase() {
      let up = c.to_ascii_uppercase()
      match table.get(up) {
        Some(v) => buf.write_char(v.to_ascii_lowercase())
        None => buf.write_char(c)
      }
    } else if c.is_ascii_uppercase() {
      match table.get(c) {
        Some(v) => buf.write_char(v)
        None => buf.write_char(c)
      }
    } else {
      buf.write_char(c)
    }
  }
  buf.to_string()
}
