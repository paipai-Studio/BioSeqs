///|
/// MeltingTemp demo: oligonucleotide melting temperatures, ported from
/// Biopython's Bio.SeqUtils.MeltingTemp. Shows the Wallace rule of thumb,
/// empirical GC-content formulas (Tm_GC) and full nearest-neighbor
/// thermodynamics (Tm_NN) with eight published parameter tables, salt /
/// Mg2+ / dNTP corrections, mismatch and dangling-end duplexes,
/// self-complementary strands, custom tables and DMSO/formamide corrections.
///
/// Run: moon run examples/melting_temp_demo

///|
fn p2(x : Double) -> Double {
  (x * 100.0).round() / 100.0
}

///|
fn main {
  println("=== MeltingTemp Demo (oligo melting temperatures) ===\n")

  // A 29-mer primer, the same example used throughout Biopython's
  // MeltingTemp documentation.
  let s = "CGTTCCAAAGATGTGGGCATGAGCTTAC"

  // 1. Wallace rule of thumb: Tm = 4*(G+C) + 2*(A+T), extended with
  //    expected contributions from ambiguous bases.
  println("1. Wallace rule of thumb")
  let wallace = @bio.tm_wallace_full(s) catch {
    _ => abort("wallace failed")
  }
  println("   Tm_Wallace(primer) = \{p2(wallace)} C")
  println("   (non-base characters such as digits/spaces are ignored)\n")

  // 2. Empirical GC-content formulas; valueset 7 is the Primer3Plus
  //    formula (default), sets 1-8 cover Marmur-Doty, QuikChange,
  //    Wetmur RNA / RNA-DNA and von Ahsen variants.
  println("2. Tm_GC empirical formulas")
  for vs in [1, 2, 3, 4, 5, 6, 7, 8] {
    let tm = @bio.tm_gc_full(s, valueset=vs) catch {
      _ => abort("tm_gc failed")
    }
    println("   valueset \{vs}: Tm = \{p2(tm)} C")
  }
  // 'X' marks a mismatch position (QuikChange formula)
  let tm_x = @bio.tm_gc_full("CTGCTGATXGCACGAGGTTATGG", valueset=2) catch {
    _ => abort("tm_gc mismatch failed")
  }
  println("   QuikChange set with one 'X' mismatch: Tm = \{p2(tm_x)} C\n")

  // 3. Nearest-neighbor thermodynamics (Allawi & SantaLucia 1997 default)
  println("3. Tm_NN nearest-neighbor thermodynamics")
  let tm_nn = @bio.tm_nn_full(s) catch {
    _ => abort("tm_nn failed")
  }
  println("   DNA_NN3 (default), 50 mM Na+: Tm = \{p2(tm_nn)} C")
  for table in ["DNA_NN1", "DNA_NN2", "DNA_NN4"] {
    let tm = @bio.tm_nn_full(s, options=@bio.TmOptions::of(nn_table=table)) catch {
      _ => abort("tm_nn table failed")
    }
    println("   \{table}: Tm = \{p2(tm)} C")
  }
  let r = "CGUUCCAAAGAUGUGGGCAUGAGCUUAC"
  for rtable in ["RNA_NN1", "RNA_NN2", "RNA_NN3"] {
    let tm = @bio.tm_nn_full(r, options=@bio.TmOptions::of(nn_table=rtable)) catch {
      _ => abort("tm_nn rna table failed")
    }
    println("   \{rtable} (RNA): Tm = \{p2(tm)} C")
  }
  let tm_hybrid = @bio.tm_nn_full(s, options=@bio.TmOptions::of(nn_table="R_DNA_NN1")) catch {
    _ => abort("tm_nn hybrid failed")
  }
  println("   R_DNA_NN1 (RNA/DNA hybrid): Tm = \{p2(tm_hybrid)} C\n")

  // 4. Salt conditions: Na+ equivalents, Tris, Mg2+ and dNTP binding.
  println("4. Salt, Mg2+ and dNTP effects (salt corrections 1-7)")
  for m in [1, 2, 3, 4, 5, 6, 7] {
    let tm = @bio.tm_nn_full(s, options=@bio.TmOptions::of(saltcorr=m)) catch {
      _ => abort("tm_nn saltcorr failed")
    }
    println("   saltcorr \{m}: Tm = \{p2(tm)} C")
  }
  let tm_tris = @bio.tm_nn_full(s, options=@bio.TmOptions::of(tris=10.0)) catch {
    _ => abort("tm_nn tris failed")
  }
  println("   + 10 mM Tris: Tm = \{p2(tm_tris)} C")
  let tm_mg = @bio.tm_nn_full(s, options=@bio.TmOptions::of(tris=10.0, mg=1.5, saltcorr=7)) catch {
    _ => abort("tm_nn mg failed")
  }
  println("   + 1.5 mM Mg2+ (Owczarzy 2008): Tm = \{p2(tm_mg)} C")
  let tm_dntp = @bio.tm_nn_full(
    s,
    options=@bio.TmOptions::of(tris=10.0, mg=1.5, dntps=0.6, saltcorr=7),
  ) catch {
    _ => abort("tm_nn dntp failed")
  }
  println("   + 0.6 mM dNTPs (sequester Mg2+): Tm = \{p2(tm_dntp)} C\n")

  // 5. Mismatches and dangling ends: c_seq is the template in 3'->5'
  //    direction; a lowercase base marks the mismatch site.
  println("5. Mismatch and dangling-end corrections")
  let template = "GCAAGGcTTCTACACCCGTACTCGAATG"
  let tm_mis = @bio.tm_nn_full(
    s,
    options=@bio.TmOptions::of(c_seq=template),
  ) catch {
    _ => abort("tm_nn mismatch failed")
  }
  println("   one internal mismatch: Tm = \{p2(tm_mis)} C")
  let tm_de = @bio.tm_nn_full(
    s,
    options=@bio.TmOptions::of(shift=1, c_seq="TGCAAGGcTTCTACACCCGTACTCGAATGC"),
  ) catch {
    _ => abort("tm_nn dangling failed")
  }
  println("   mismatch with one dangling end each side: Tm = \{p2(tm_de)} C")
  let tm_self = @bio.tm_nn_full(
    "CGCGAATTCGCG",
    options=@bio.TmOptions::of(selfcomp=true),
  ) catch {
    _ => abort("tm_nn selfcomp failed")
  }
  println("   self-complementary CGCGAATTCGCG: Tm = \{p2(tm_self)} C\n")

  // 6. Chemical additives: DMSO and formamide lower the melting temp.
  println("6. DMSO / formamide chemical corrections")
  let tm_dmso = @bio.tm_chem_correction(70.0, dmso=3.0) catch {
    _ => abort("dmso failed")
  }
  println("   70 C with 3% DMSO: Tm = \{p2(tm_dmso)} C")
  let tm_fmd = @bio.tm_chem_correction(70.0, fmd=5.0) catch {
    _ => abort("fmd failed")
  }
  println("   70 C with 5% formamide: Tm = \{p2(tm_fmd)} C")
}
