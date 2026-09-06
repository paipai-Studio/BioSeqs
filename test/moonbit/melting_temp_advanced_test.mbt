// Tests for the full Biopython MeltingTemp port
// (src/melting_temp_advanced.mbt): Tm_Wallace / Tm_GC / Tm_NN with the eight
// NN thermodynamic tables, salt corrections 1-7, dangling ends, terminal and
// internal mismatches, self-complementary duplexes, custom tables and
// DMSO/formamide corrections. Reference values were produced with
// Biopython 1.88 Bio.SeqUtils.MeltingTemp.

///|
fn mt_approx(actual : Double, expected : Double) -> Bool {
  (actual - expected).abs() < 0.02
}

///|
test "mt_advanced_wallace" {
  let s = "CGTTCCAAAGATGTGGGCATGAGCTTAC"
  assert_true(mt_approx(@src.tm_wallace_full(s), 84.0))
  // non-DNA characters and whitespace are ignored
  assert_true(mt_approx(@src.tm_wallace_full("1ACGT 2TGCA\n3ATGC\t4CGTA"), 48.0))
  // ambiguous bases with expected intermediate values
  assert_true(mt_approx(@src.tm_wallace_full("W"), 2.0))
  assert_true(mt_approx(@src.tm_wallace_full("S"), 4.0))
  assert_true(mt_approx(
    @src.tm_wallace_full("N", options=@src.TmOptions::of(strict=false)),
    3.0,
  ))
  // strict mode rejects ambiguous bases
  let raised = try {
    ignore(@src.tm_wallace_full("ACGTN"))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(raised)
}

///|
test "mt_advanced_gc_valuesets" {
  let s = "CGTTCCAAAGATGTGGGCATGAGCTTAC"
  // valueset 7 (Primer3Plus) is the default
  assert_true(mt_approx(@src.tm_gc_full(s), 58.9743))
  // Biopython 1.88 references for valuesets 1-8 (Na=50 mM defaults)
  assert_true(mt_approx(@src.tm_gc_full(s, valueset=1), 66.5857))
  assert_true(mt_approx(@src.tm_gc_full(s, valueset=2), 77.8929))
  assert_true(mt_approx(@src.tm_gc_full(s, valueset=3), 56.2958))
  assert_true(mt_approx(@src.tm_gc_full(s, valueset=4), 62.2977))
  assert_true(mt_approx(@src.tm_gc_full(s, valueset=5), 73.2977))
  assert_true(mt_approx(@src.tm_gc_full(s, valueset=6), 67.2977))
  assert_true(mt_approx(@src.tm_gc_full(s, valueset=7), 58.9743))
  assert_true(mt_approx(@src.tm_gc_full(s, valueset=8), 63.5208))
  // QuikChange formula with an X mismatch position
  assert_true(mt_approx(
    @src.tm_gc_full("CTGCTGATXGCACGAGGTTATGG", valueset=2),
    69.1957,
  ))
  // custom constant set, equivalent to the valueset 7 constants
  let u = Some((81.5, 0.41, 600.0, 1.0))
  assert_true(mt_approx(
    @src.tm_gc_full(s, userset=u, options=@src.TmOptions::of(saltcorr=1)),
    58.9743,
  ))
  // invalid valueset
  let bad_valueset = try {
    ignore(@src.tm_gc_full(s, valueset=9))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(bad_valueset)
  // strict mode rejects ambiguous bases
  let ambiguous = try {
    ignore(@src.tm_gc_full("ACGTKMNRYBVDH"))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(ambiguous)
}

///|
test "mt_advanced_nn_default_and_tables" {
  let s = "CGTTCCAAAGATGTGGGCATGAGCTTAC"
  // Allawi & SantaLucia (1997) DNA_NN3 is the default table
  assert_true(mt_approx(@src.tm_nn_full(s), 60.3209))
  // the other DNA/DNA tables
  assert_true(mt_approx(
    @src.tm_nn_full(s, options=@src.TmOptions::of(nn_table="DNA_NN1")),
    72.1950,
  ))
  assert_true(mt_approx(
    @src.tm_nn_full(s, options=@src.TmOptions::of(nn_table="DNA_NN2")),
    65.4693,
  ))
  assert_true(mt_approx(
    @src.tm_nn_full(s, options=@src.TmOptions::of(nn_table="DNA_NN4")),
    60.2715,
  ))
  // RNA/RNA tables (sequences are back-transcribed internally)
  let r = "CGUUCCAAAGAUGUGGGCAUGAGCUUAC"
  assert_true(mt_approx(
    @src.tm_nn_full(r, options=@src.TmOptions::of(nn_table="RNA_NN1")),
    73.3504,
  ))
  assert_true(mt_approx(
    @src.tm_nn_full(r, options=@src.TmOptions::of(nn_table="RNA_NN2")),
    76.4178,
  ))
  assert_true(mt_approx(
    @src.tm_nn_full(r, options=@src.TmOptions::of(nn_table="RNA_NN3")),
    75.1915,
  ))
  // RNA/DNA hybrid (Sugimoto 1995); seq may be given as DNA here
  assert_true(mt_approx(
    @src.tm_nn_full(s, options=@src.TmOptions::of(nn_table="R_DNA_NN1")),
    58.4496,
  ))
  // unknown table name raises
  let unknown = try {
    ignore(@src.tm_nn_full(s, options=@src.TmOptions::of(nn_table="NOPE")))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(unknown)
}

///|
test "mt_advanced_nn_salt_corrections" {
  let s = "CGTTCCAAAGATGTGGGCATGAGCTTAC"
  // salt corrections 1-7 (Na=50 mM default)
  assert_true(mt_approx(
    @src.tm_nn_full(s, options=@src.TmOptions::of(saltcorr=1)),
    54.2657,
  ))
  assert_true(mt_approx(
    @src.tm_nn_full(s, options=@src.TmOptions::of(saltcorr=2)),
    54.0177,
  ))
  assert_true(mt_approx(
    @src.tm_nn_full(s, options=@src.TmOptions::of(saltcorr=3)),
    59.5999,
  ))
  assert_true(mt_approx(
    @src.tm_nn_full(s, options=@src.TmOptions::of(saltcorr=4)),
    60.6408,
  ))
  assert_true(mt_approx(
    @src.tm_nn_full(s, options=@src.TmOptions::of(saltcorr=5)),
    60.3209,
  ))
  assert_true(mt_approx(
    @src.tm_nn_full(s, options=@src.TmOptions::of(saltcorr=6)),
    59.7775,
  ))
  assert_true(mt_approx(
    @src.tm_nn_full(s, options=@src.TmOptions::of(saltcorr=7)),
    59.7775,
  ))
  // other cations: Tris, Mg2+, dNTPs
  assert_true(mt_approx(
    @src.tm_nn_full(s, options=@src.TmOptions::of(tris=10.0)),
    60.7940,
  ))
  assert_true(mt_approx(
    @src.tm_nn_full(s, options=@src.TmOptions::of(tris=10.0, mg=1.5)),
    67.3880,
  ))
  assert_true(mt_approx(
    @src.tm_nn_full(s,
      options=@src.TmOptions::of(tris=10.0, mg=1.5, saltcorr=7)),
    66.8062,
  ))
  assert_true(mt_approx(
    @src.tm_nn_full(s,
      options=@src.TmOptions::of(tris=10.0, mg=1.5, dntps=0.6, saltcorr=7)),
    66.0371,
  ))
}

///|
test "mt_advanced_salt_correction_raw" {
  let s = "CGTTCCAAAGATGTGGGCATGAGCTTAC"
  assert_true(mt_approx(
    @src.tm_salt_correction(na=50.0, salt_method=1),
    -21.5971,
  ))
  assert_true(mt_approx(
    @src.tm_salt_correction(na=50.0, salt_method=2),
    -21.8451,
  ))
  assert_true(mt_approx(
    @src.tm_salt_correction(na=100.0, tris=20.0, salt_method=2, seq=s),
    -16.4477,
  ))
  assert_true(mt_approx(
    @src.tm_salt_correction(
      na=100.0, tris=20.0, mg=1.5, salt_method=2, seq=s,
    ),
    -10.9885,
  ))
  // method > 7 raises
  let bad_method = try {
    ignore(@src.tm_salt_correction(na=50.0, salt_method=8))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(bad_method)
  // methods 5/6/7 need the sequence
  let missing_seq = try {
    ignore(@src.tm_salt_correction(na=50.0, salt_method=5))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(missing_seq)
}

///|
test "mt_advanced_chem_correction" {
  assert_true(mt_approx(@src.tm_chem_correction(70.0), 70.0))
  assert_true(mt_approx(@src.tm_chem_correction(70.0, dmso=3.0), 67.75))
  assert_true(mt_approx(@src.tm_chem_correction(70.0, fmd=5.0), 66.75))
  assert_true(mt_approx(
    @src.tm_chem_correction(70.0, fmd=1.25, fmd_method=2, gc=50.0),
    66.6831,
  ))
  // fmdmethod=2 needs GC content
  let no_gc = try {
    ignore(@src.tm_chem_correction(70.0, fmd=1.0, fmd_method=2))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(no_gc)
}

///|
test "mt_advanced_nn_mismatches_and_dangling_ends" {
  let s = "CGTTCCAAAGATGTGGGCATGAGCTTAC"
  // internal mismatch flagged by a lowercase base in c_seq
  assert_true(mt_approx(
    @src.tm_nn_full(s,
      options=@src.TmOptions::of(c_seq="GCAAGGcTTCTACACCCGTACTCGAATG")),
    55.3892,
  ))
  // dangling ends: c_seq longer by 1 on each side plus shift
  assert_true(mt_approx(
    @src.tm_nn_full(s,
      options=@src.TmOptions::of(
        shift=1, c_seq="TGCAAGGcTTCTACACCCGTACTCGAATGC",
      )),
    55.6878,
  ))
  // RNA duplex with mismatch and dangling ends (RNA_NN3 + RNA_DE1); the
  // default DNA TMM/IMM tables are kept, matching the Biopython doctest
  let r = "CGUUCCAAAGAUGUGGGCAUGAGCUUAC"
  assert_true(mt_approx(
    @src.tm_nn_full(r,
      options=@src.TmOptions::of(
        shift=1, nn_table="RNA_NN3", de_table="RNA_DE1",
        c_seq="UGCAAGGcUUCUACACCCGUACUCGAAUGC",
      )),
    72.9977,
  ))
}

///|
test "mt_advanced_nn_self_complementary" {
  let s = "CGCGAATTCGCG"
  assert_true(mt_approx(
    @src.tm_nn_full(s, options=@src.TmOptions::of(selfcomp=true)),
    42.4019,
  ))
  assert_true(mt_approx(@src.tm_nn_full(s), 42.4242))
}

///|
test "mt_advanced_make_table_custom" {
  // replicate the Biopython make_table doctest: custom NN keys for a
  // modified base '1' require check=false so the base is not discarded
  let values : @src.ThermoTable = Map([
    ("A1/T1", (-6.608, -17.235)),
    ("1A/1T", (-6.893, -15.923)),
  ])
  let mytable = @src.tm_make_table(base="DNA_NN3", values=values)
  assert_true(mt_approx(
    @src.tm_nn_full("CGTTCCA1AGATGTGGGCATGAGCTTAC",
      options=@src.TmOptions::of(nn_custom=mytable, check=false)),
    62.5271,
  ))
  // blank table via tm_make_table
  let blank = @src.tm_make_table()
  assert_eq(blank.length(), 17)
}
