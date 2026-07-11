// Complement lookup tables, mirroring BioPython's
// Bio.Data.IUPACData.ambiguous_dna_complement / ambiguous_rna_complement.
//
// BioPython's Seq.complement() treats any U as a T (i.e. it uses the DNA
// complement table augmented with U->A), while Seq.complement_rna() treats
// any T as a U. See Bio/Seq.py lines 53-60 and 1654/1703.
//
// Implemented with FixedArray[UInt16] of size 128 for O(1) direct-index
// lookup instead of Map[Char, Char]. 0 means "no complement" (pass-through).
// Using UInt16 so values can be directly written to a FixedArray[UInt16]
// output buffer without conversion.

///|
/// Intrinsic: create a String directly from a FixedArray[UInt16] of UTF-16
/// code units, without copying. This is much faster than StringBuilder for
/// building strings from code units.
fn unsafe_from_uint16_fixedarray(buffer : FixedArray[UInt16]) -> String = "%string.unsafe_from_uint16_fixedarray"

///|
/// DNA complement table indexed by ASCII code. Entries for BOTH uppercase
/// and lowercase are pre-computed so the inner loop has zero branches.
/// 0 means "pass through unchanged" (for non-nucleotide characters).
/// Includes U->A so that Seq.complement() works on RNA input.
let dna_complement_arr : FixedArray[UInt16] = {
  let t : FixedArray[UInt16] = FixedArray::make(128, 0)
  // Uppercase
  t['A'.to_int()] = 'T'.to_int().to_uint16()
  t['C'.to_int()] = 'G'.to_int().to_uint16()
  t['G'.to_int()] = 'C'.to_int().to_uint16()
  t['T'.to_int()] = 'A'.to_int().to_uint16()
  t['U'.to_int()] = 'A'.to_int().to_uint16()
  t['M'.to_int()] = 'K'.to_int().to_uint16()
  t['R'.to_int()] = 'Y'.to_int().to_uint16()
  t['W'.to_int()] = 'W'.to_int().to_uint16()
  t['S'.to_int()] = 'S'.to_int().to_uint16()
  t['Y'.to_int()] = 'R'.to_int().to_uint16()
  t['K'.to_int()] = 'M'.to_int().to_uint16()
  t['V'.to_int()] = 'B'.to_int().to_uint16()
  t['H'.to_int()] = 'D'.to_int().to_uint16()
  t['D'.to_int()] = 'H'.to_int().to_uint16()
  t['B'.to_int()] = 'V'.to_int().to_uint16()
  t['X'.to_int()] = 'X'.to_int().to_uint16()
  t['N'.to_int()] = 'N'.to_int().to_uint16()
  // Lowercase (pre-computed to eliminate branches in inner loop)
  t['a'.to_int()] = 't'.to_int().to_uint16()
  t['c'.to_int()] = 'g'.to_int().to_uint16()
  t['g'.to_int()] = 'c'.to_int().to_uint16()
  t['t'.to_int()] = 'a'.to_int().to_uint16()
  t['u'.to_int()] = 'a'.to_int().to_uint16()
  t['m'.to_int()] = 'k'.to_int().to_uint16()
  t['r'.to_int()] = 'y'.to_int().to_uint16()
  t['w'.to_int()] = 'w'.to_int().to_uint16()
  t['s'.to_int()] = 's'.to_int().to_uint16()
  t['y'.to_int()] = 'r'.to_int().to_uint16()
  t['k'.to_int()] = 'm'.to_int().to_uint16()
  t['v'.to_int()] = 'b'.to_int().to_uint16()
  t['h'.to_int()] = 'd'.to_int().to_uint16()
  t['d'.to_int()] = 'h'.to_int().to_uint16()
  t['b'.to_int()] = 'v'.to_int().to_uint16()
  t['x'.to_int()] = 'x'.to_int().to_uint16()
  t['n'.to_int()] = 'n'.to_int().to_uint16()
  t
}

///|
/// RNA complement table indexed by ASCII code. Includes BOTH uppercase and
/// lowercase entries. Includes T->A so that Seq.complement_rna() works on
/// DNA input.
let rna_complement_arr : FixedArray[UInt16] = {
  let t : FixedArray[UInt16] = FixedArray::make(128, 0)
  // Uppercase
  t['A'.to_int()] = 'U'.to_int().to_uint16()
  t['C'.to_int()] = 'G'.to_int().to_uint16()
  t['G'.to_int()] = 'C'.to_int().to_uint16()
  t['U'.to_int()] = 'A'.to_int().to_uint16()
  t['T'.to_int()] = 'A'.to_int().to_uint16()
  t['M'.to_int()] = 'K'.to_int().to_uint16()
  t['R'.to_int()] = 'Y'.to_int().to_uint16()
  t['W'.to_int()] = 'W'.to_int().to_uint16()
  t['S'.to_int()] = 'S'.to_int().to_uint16()
  t['Y'.to_int()] = 'R'.to_int().to_uint16()
  t['K'.to_int()] = 'M'.to_int().to_uint16()
  t['V'.to_int()] = 'B'.to_int().to_uint16()
  t['H'.to_int()] = 'D'.to_int().to_uint16()
  t['D'.to_int()] = 'H'.to_int().to_uint16()
  t['B'.to_int()] = 'V'.to_int().to_uint16()
  t['X'.to_int()] = 'X'.to_int().to_uint16()
  t['N'.to_int()] = 'N'.to_int().to_uint16()
  // Lowercase
  t['a'.to_int()] = 'u'.to_int().to_uint16()
  t['c'.to_int()] = 'g'.to_int().to_uint16()
  t['g'.to_int()] = 'c'.to_int().to_uint16()
  t['u'.to_int()] = 'a'.to_int().to_uint16()
  t['t'.to_int()] = 'a'.to_int().to_uint16()
  t['m'.to_int()] = 'k'.to_int().to_uint16()
  t['r'.to_int()] = 'y'.to_int().to_uint16()
  t['w'.to_int()] = 'w'.to_int().to_uint16()
  t['s'.to_int()] = 's'.to_int().to_uint16()
  t['y'.to_int()] = 'r'.to_int().to_uint16()
  t['k'.to_int()] = 'm'.to_int().to_uint16()
  t['v'.to_int()] = 'b'.to_int().to_uint16()
  t['h'.to_int()] = 'd'.to_int().to_uint16()
  t['d'.to_int()] = 'h'.to_int().to_uint16()
  t['b'.to_int()] = 'v'.to_int().to_uint16()
  t['x'.to_int()] = 'x'.to_int().to_uint16()
  t['n'.to_int()] = 'n'.to_int().to_uint16()
  t
}

///|
/// Complement every character in `s` using the given FixedArray lookup table.
/// The table includes both uppercase and lowercase entries, so the inner loop
/// is branch-free: just one array lookup per character. Builds the result
/// directly in a FixedArray[UInt16] and converts to String via intrinsic.
fn complement_chars(s : String, table : FixedArray[UInt16]) -> String {
  let n = s.length()
  let buf : FixedArray[UInt16] = FixedArray::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    let ci = s.unsafe_get(i).to_int()
    if ci < 128 {
      let comp = table[ci]
      if comp != 0 {
        buf[i] = comp
      } else {
        buf[i] = s.unsafe_get(i)
      }
    } else {
      buf[i] = s.unsafe_get(i)
    }
  }
  unsafe_from_uint16_fixedarray(buf)
}

///|
/// Reverse-complement in a single pass: read from the end of `s`, complement
/// each character, and write to the buffer. Avoids the intermediate reversed
/// string allocation of `complement_chars(s.rev(), table)`.
fn reverse_complement_chars(s : String, table : FixedArray[UInt16]) -> String {
  let n = s.length()
  let buf : FixedArray[UInt16] = FixedArray::make(n, 0)
  for i = n - 1; i >= 0; i = i - 1 {
    let ci = s.unsafe_get(i).to_int()
    let out_idx = n - 1 - i
    if ci < 128 {
      let comp = table[ci]
      if comp != 0 {
        buf[out_idx] = comp
      } else {
        buf[out_idx] = s.unsafe_get(i)
      }
    } else {
      buf[out_idx] = s.unsafe_get(i)
    }
  }
  unsafe_from_uint16_fixedarray(buf)
}
