/// seq_variant_library.mbt — combinatorial variant libraries
/// (virtual saturation mutagenesis, item B4)
///
/// Site-saturation over an alphabet, degenerate-codon schemes (NNK/NNS/
/// NDT/...), free per-position combination libraries, and reverse
/// translation of proteins into degenerate DNA sequences.

///|
/// Nucleotide set covered by an IUPAC ambiguity code.
pub fn seq_var_iupac_bases(code : Char) -> Array[Char] {
  match code {
    'A' => ['A']
    'C' => ['C']
    'G' => ['G']
    'T' => ['T']
    'U' => ['U']
    'R' => ['A', 'G']
    'Y' => ['C', 'T']
    'K' => ['G', 'T']
    'M' => ['A', 'C']
    'S' => ['C', 'G']
    'W' => ['A', 'T']
    'B' => ['C', 'G', 'T']
    'D' => ['A', 'G', 'T']
    'H' => ['A', 'C', 'T']
    'V' => ['A', 'C', 'G']
    'N' => ['A', 'C', 'G', 'T']
    _ => []
  }
}

///|
/// Map a set of nucleotides to the tightest IUPAC code.
pub fn seq_var_iupac_from_set(bases : Set[Char]) -> Char {
  let has = (c : Char) => bases.contains(c)
  let n = bases.length()
  if n == 1 {
    if has('A') { 'A' } else if has('C') {
      'C'
    } else if has('G') {
      'G'
    } else {
      'T'
    }
  } else if n == 2 {
    if has('A') && has('G') {
      'R'
    } else if has('C') && has('T') {
      'Y'
    } else if has('G') && has('T') {
      'K'
    } else if has('A') && has('C') {
      'M'
    } else if has('C') && has('G') {
      'S'
    } else {
      'W'
    }
  } else if n == 3 {
    if !has('A') {
      'B'
    } else if !has('C') {
      'D'
    } else if !has('G') {
      'H'
    } else {
      'V'
    }
  } else {
    'N'
  }
}

///|
fn seq_var_string_to_seq(s : String) -> Seq {
  Seq::new(s)
}

///|
/// Enumerate all single-position substitutions of `pos` over `alphabet`
/// (e.g. "ACGT" for DNA saturation, the 20 amino acids for proteins).
/// The wild-type sequence is excluded from the result.
pub fn seq_saturate_site(
  seq : Seq,
  pos : Int32,
  alphabet~ : String = "ACGT",
  max_variants~ : Int32 = 256,
) -> Array[Seq] {
  let out : Array[Seq] = Array(capacity=0)
  let s = seq.to_string()
  let n = s.length()
  if pos < 0 || pos >= n {
    return out
  }
  let wild = s.unsafe_get(pos)
  for i = 0; i < alphabet.length(); i = i + 1 {
    if alphabet.unsafe_get(i) == wild {
      continue
    }
    if out.length() >= max_variants {
      break
    }
    let sb = StringBuilder()
    for j = 0; j < n; j = j + 1 {
      if j == pos {
        let _ = sb.write_char(alphabet.unsafe_get(i).unsafe_to_char())
      } else {
        let _ = sb.write_char(s.unsafe_get(j).unsafe_to_char())
      }
    }
    out.push(seq_var_string_to_seq(sb.to_string()))
  }
  out
}

///|
/// Write a degenerate-codon scheme (e.g. "NNK", "NNS", "NDT") over the
/// given codon sites (0-based amino-acid positions; DNA positions
/// 3*i .. 3*i+3 are replaced by the scheme's three codes).
pub fn seq_degenerate_sequence(
  seq : Seq,
  codon_sites : Array[Int32],
  scheme~ : String = "NNK",
) -> Seq {
  let s = seq.to_string()
  let n = s.length()
  let sb = StringBuilder()
  for j = 0; j < n; j = j + 1 {
    let _ = sb.write_char(s.unsafe_get(j).unsafe_to_char())
  }
  let chars = sb.to_string()
  let out = StringBuilder()
  let mut i = 0
  while i < n {
    let site = i / 3
    let offset = i % 3
    let mut in_site = false
    for cs in codon_sites {
      if cs == site {
        in_site = true
      }
    }
    if in_site && offset < scheme.length() {
      let _ = out.write_char(scheme.unsafe_get(offset).unsafe_to_char())
    } else {
      let _ = out.write_char(chars.unsafe_get(i).unsafe_to_char())
    }
    i = i + 1
  }
  seq_var_string_to_seq(out.to_string())
}

///|
/// Number of concrete nucleotide combinations covered by the scheme over
/// the given codon sites (product of per-position alphabet sizes).
pub fn seq_degenerate_variant_count(
  codon_sites : Array[Int32],
  scheme~ : String = "NNK",
) -> Int {
  let mut total = 1
  for _site in codon_sites {
    for i = 0; i < scheme.length(); i = i + 1 {
      total = total * seq_var_iupac_bases(scheme.unsafe_get(i).unsafe_to_char()).length()
    }
  }
  total
}

///|
/// Expand a degenerate scheme into concrete DNA sequences (Cartesian
/// product over the per-position nucleotide sets), capped by
/// `max_variants`.
pub fn seq_degenerate_variants(
  seq : Seq,
  codon_sites : Array[Int32],
  scheme~ : String = "NNK",
  max_variants~ : Int32 = 1024,
) -> Array[Seq] {
  let degenerate = seq_degenerate_sequence(seq, codon_sites, scheme~)
  let base = degenerate.to_string()
  // collect (position, choices) for every position carrying an ambiguity code
  let positions : Array[Int32] = Array(capacity=0)
  let choices : Array[Array[Char]] = Array(capacity=0)
  let n = base.length()
  for i = 0; i < n; i = i + 1 {
    let c = base.unsafe_get(i).unsafe_to_char()
    let opts = seq_var_iupac_bases(c)
    if opts.length() > 1 {
      positions.push(i)
      choices.push(opts)
    }
  }
  let out : Array[Seq] = Array(capacity=0)
  if positions.length() == 0 {
    out.push(seq_var_string_to_seq(base))
    return out
  }
  // odometer over the ambiguity positions
  let counters : Array[Int32] = Array::make(positions.length(), 0)
  let mut done = false
  while !done {
    if out.length() >= max_variants {
      break
    }
    let sb = StringBuilder()
    let mut pi = 0
    for i = 0; i < n; i = i + 1 {
      if pi < positions.length() && positions[pi] == i {
        let _ = sb.write_char(choices[pi][counters[pi]])
        pi = pi + 1
      } else {
        let _ = sb.write_char(base.unsafe_get(i).unsafe_to_char())
      }
    }
    out.push(seq_var_string_to_seq(sb.to_string()))
    // increment odometer
    let mut carry = positions.length() - 1
    let mut overflow = true
    while carry >= 0 {
      counters[carry] = counters[carry] + 1
      if counters[carry] < choices[carry].length() {
        overflow = false
        break
      }
      counters[carry] = 0
      carry = carry - 1
    }
    if overflow {
      done = true
    }
  }
  out
}

///|
/// Free combination library: `variants` maps a position to the list of
/// allowed characters there; every combination is materialized
/// (Cartesian product, capped at `max_combos`).
pub fn seq_combo_library(
  seq : Seq,
  variants : Map[Int32, Array[String]],
  max_combos~ : Int32 = 4096,
) -> Array[Seq] {
  let out : Array[Seq] = Array(capacity=0)
  let s = seq.to_string()
  let n = s.length()
  // deterministic position order
  let positions : Array[Int32] = Array(capacity=0)
  for pos, _opts in variants {
    positions.push(pos)
  }
  for i = 1; i < positions.length(); i = i + 1 {
    let key = positions[i]
    let mut j = i - 1
    while j >= 0 && positions[j] > key {
      positions[j + 1] = positions[j]
      j = j - 1
    }
    positions[j + 1] = key
  }
  let choice_sets : Array[Array[String]] = Array(capacity=0)
  for pos in positions {
    choice_sets.push(variants[pos])
  }
  let counters : Array[Int32] = Array::make(positions.length(), 0)
  let mut done = positions.length() == 0
  while !done {
    if out.length() >= max_combos {
      break
    }
    let sb = StringBuilder()
    let mut vi = 0
    for i = 0; i < n; i = i + 1 {
      if vi < positions.length() && positions[vi] == i {
        let _ = sb.write(choice_sets[vi][counters[vi]])
        vi = vi + 1
      } else {
        let _ = sb.write_char(s.unsafe_get(i).unsafe_to_char())
      }
    }
    out.push(seq_var_string_to_seq(sb.to_string()))
    let mut carry = positions.length() - 1
    let mut overflow = true
    while carry >= 0 {
      counters[carry] = counters[carry] + 1
      if counters[carry] < choice_sets[carry].length() {
        overflow = false
        break
      }
      counters[carry] = 0
      carry = carry - 1
    }
    if overflow {
      done = true
    }
  }
  out
}

///|
/// Reverse-translate a protein into a degenerate DNA sequence: every
/// amino acid is expanded to its codon family, and each codon position
/// is aggregated into the tightest IUPAC code. Returns None when an
/// amino acid has no codons in the table (e.g. X, B, Z).
pub fn seq_reverse_translate_degenerate(
  protein : Seq,
  table_id~ : Int32 = 1,
) -> Option[String] {
  let s = protein.to_string()
  let n = s.length()
  if n == 0 {
    return None
  }
  let rev = codon_table_reverse(table_id~)
  let sb = StringBuilder()
  for i = 0; i < n; i = i + 1 {
    let aa = s.unsafe_get(i).unsafe_to_char().to_string()
    match rev.get(aa) {
      None => return None
      Some(codons) => {
        if codons.length() == 0 {
          return None
        }
        // aggregate bases per codon position
        for offset = 0; offset < 3; offset = offset + 1 {
          let bases : Set[Char] = Set([])
          for cd in codons {
            if cd.length() == 3 {
              let _ = bases.add(cd.unsafe_get(offset).unsafe_to_char())
            }
          }
          let _ = sb.write_char(seq_var_iupac_from_set(bases))
        }
      }
    }
  }
  Some(sb.to_string())
}
