/// seq_descriptors.mbt — batch sequence descriptors
/// (CalcMolDescriptors analogue, item B1)
///
/// One-call descriptor maps for proteins (composition, molecular weight,
/// GRAVY, aromaticity, isoelectric point, charge at pH) and nucleotides
/// (base counts, GC/AT content, purine fraction, CpG count, Shannon
/// entropy), mirroring the rdkit_calc_descriptors batch pattern.

///|
/// Kyte-Doolittle hydropathy values.
fn seq_desc_kd(a : Char) -> Double {
  match a {
    'A' => 1.8
    'R' => -4.5
    'N' => -3.5
    'D' => -3.5
    'C' => 2.5
    'Q' => -3.5
    'E' => -3.5
    'G' => -0.4
    'H' => -3.2
    'I' => 4.5
    'L' => 3.8
    'K' => -3.9
    'M' => 1.9
    'F' => 2.8
    'P' => -1.6
    'S' => -0.8
    'T' => -0.7
    'W' => -0.9
    'Y' => -1.3
    'V' => 4.2
    _ => 0.0
  }
}

///|
/// Average residue masses (without water).
fn seq_desc_residue_mass(a : Char) -> Double {
  match a {
    'A' => 71.0788
    'R' => 156.1875
    'N' => 114.1038
    'D' => 115.0886
    'C' => 103.1388
    'E' => 129.1155
    'Q' => 128.1307
    'G' => 57.0519
    'H' => 137.1411
    'I' => 113.1594
    'L' => 113.1594
    'K' => 128.1741
    'M' => 131.1926
    'F' => 147.1766
    'P' => 97.1167
    'S' => 87.0782
    'T' => 101.1051
    'W' => 186.2132
    'Y' => 163.1760
    'V' => 99.1326
    _ => 0.0
  }
}

///|
/// Net protein charge at a given pH (Sillero pKa set).
pub fn seq_protein_charge_at_ph(seq : Seq, ph : Double) -> Double {
  let s = seq.to_string()
  let n = s.length()
  if n == 0 {
    return 0.0
  }
  // positive group counts
  let mut n_k = 0
  let mut n_r = 0
  let mut n_h = 0
  // negative group counts
  let mut n_d = 0
  let mut n_e = 0
  let mut n_c = 0
  let mut n_y = 0
  for i = 0; i < n; i = i + 1 {
    let c = s.unsafe_get(i).unsafe_to_char()
    match c {
      'K' => n_k = n_k + 1
      'R' => n_r = n_r + 1
      'H' => n_h = n_h + 1
      'D' => n_d = n_d + 1
      'E' => n_e = n_e + 1
      'C' => n_c = n_c + 1
      'Y' => n_y = n_y + 1
      _ => ()
    }
  }
  // termini
  let pos_part = n_k.to_double() / (1.0 + @math.pow(10.0, ph - 10.0)) +
    n_r.to_double() / (1.0 + @math.pow(10.0, ph - 12.0)) +
    n_h.to_double() / (1.0 + @math.pow(10.0, ph - 5.98)) +
    1.0 / (1.0 + @math.pow(10.0, ph - 7.5))
  let neg_part = n_d.to_double() / (1.0 + @math.pow(10.0, 4.05 - ph)) +
    n_e.to_double() / (1.0 + @math.pow(10.0, 4.45 - ph)) +
    n_c.to_double() / (1.0 + @math.pow(10.0, 9.0 - ph)) +
    n_y.to_double() / (1.0 + @math.pow(10.0, 10.0 - ph)) +
    1.0 / (1.0 + @math.pow(10.0, 3.55 - ph))
  pos_part - neg_part
}

///|
/// Isoelectric point via bisection on the net charge curve.
pub fn seq_protein_isoelectric_point(seq : Seq) -> Double {
  let n = seq.length()
  if n == 0 {
    return 7.0
  }
  let mut lo = 0.0
  let mut hi = 14.0
  // charge decreases monotonically with pH
  for _iter = 0; _iter < 80; _iter = _iter + 1 {
    let mid = (lo + hi) / 2.0
    if seq_protein_charge_at_ph(seq, mid) > 0.0 {
      lo = mid
    } else {
      hi = mid
    }
  }
  (lo + hi) / 2.0
}

///|
/// Batch protein descriptors: per-residue composition, molecular weight
/// (average residue masses + water), GRAVY, aromaticity (F+W+Y
/// fraction), isoelectric point and charge at pH 7.
pub fn seq_calc_protein_descriptors(seq : Seq) -> Map[String, Double] {
  let out : Map[String, Double] = Map([], capacity=32)
  let s = seq.to_string()
  let n = s.length()
  out["length"] = n.to_double()
  if n == 0 {
    return out
  }
  let letters = "ACDEFGHIKLMNPQRSTVWY"
  let counts : Map[Char, Int32] = Map([], capacity=24)
  for ci = 0; ci < letters.length(); ci = ci + 1 {
    counts[letters.unsafe_get(ci).unsafe_to_char()] = 0
  }
  let mut mass = 18.01528 // water
  let mut gravy_sum = 0.0
  let mut aromatic = 0
  let mut counted = 0
  for i = 0; i < n; i = i + 1 {
    let c = s.unsafe_get(i).unsafe_to_char()
    match counts.get(c) {
      Some(v) => {
        counts[c] = v + 1
        counted = counted + 1
      }
      None => ()
    }
    mass = mass + seq_desc_residue_mass(c)
    gravy_sum = gravy_sum + seq_desc_kd(c)
    match c {
      'F' | 'W' | 'Y' => aromatic = aromatic + 1
      _ => ()
    }
  }
  for a, cnt in counts {
    out["A_" + a.to_string()] = cnt.to_double() / n.to_double()
  }
  out["counted"] = counted.to_double()
  out["mol_weight"] = mass
  out["gravy"] = gravy_sum / n.to_double()
  out["aromaticity"] = aromatic.to_double() / n.to_double()
  out["charge_ph7"] = seq_protein_charge_at_ph(seq, 7.0)
  out["isoelectric_point"] = seq_protein_isoelectric_point(seq)
  out
}

///|
/// Batch nucleotide descriptors: base counts (fractions), GC content,
/// AT(U) content, purine fraction, CpG dinucleotide count (per kb) and
/// Shannon entropy in bits.
pub fn seq_calc_nucleotide_descriptors(seq : Seq) -> Map[String, Double] {
  let out : Map[String, Double] = Map([], capacity=24)
  let s = seq.to_string()
  let n = s.length()
  out["length"] = n.to_double()
  if n == 0 {
    return out
  }
  let mut a = 0
  let mut c = 0
  let mut g = 0
  let mut t = 0
  let mut u = 0
  let mut other = 0
  for i = 0; i < n; i = i + 1 {
    let ch = s.unsafe_get(i).unsafe_to_char()
    match ch {
      'A' => a = a + 1
      'C' => c = c + 1
      'G' => g = g + 1
      'T' => t = t + 1
      'U' => u = u + 1
      _ => other = other + 1
    }
  }
  let nf = n.to_double()
  out["A"] = a.to_double() / nf
  out["C"] = c.to_double() / nf
  out["G"] = g.to_double() / nf
  out["T"] = t.to_double() / nf
  out["U"] = u.to_double() / nf
  out["other"] = other.to_double() / nf
  let gc = (c + g).to_double() / nf
  let at = (a + t + u).to_double() / nf
  out["GC"] = gc
  out["AT"] = at
  out["purine"] = (a + g).to_double() / nf
  // CpG count (overlapping), normalized per kilobase
  let mut cpg = 0
  for i = 0; i < n - 1; i = i + 1 {
    if s.unsafe_get(i).unsafe_to_char() == 'C' &&
      s.unsafe_get(i + 1).unsafe_to_char() == 'G' {
      cpg = cpg + 1
    }
  }
  out["cpg_count"] = cpg.to_double()
  out["cpg_per_kb"] = cpg.to_double() / (nf / 1000.0)
  // Shannon entropy over the observed alphabet
  out["shannon_bits"] = shannon_entropy_bits(s)
  out
}
