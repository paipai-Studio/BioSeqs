///|
/// Bio.SeqUtils advanced function tests, validated against Biopython 1.88.
/// Reference values were produced with:
///   python3 -c "from Bio.SeqUtils import ..."

///|
test "gc_fraction_default_remove" {
  assert_eq(@src.gc_fraction("ACTG"), 0.5)
  assert_eq(@src.gc_fraction("GGAUCUUCGGAUCU"), 0.5)
  assert_eq(@src.gc_fraction("ACTGSSSS"), 0.75)
  assert_eq(@src.gc_fraction("ACTGN"), 0.5)
  assert_eq(@src.gc_fraction("GDVV"), 1.0)
  // Mixed case is handled
  assert_eq(@src.gc_fraction("actg"), 0.5)
}

///|
test "gc_fraction_ignore" {
  assert_eq(@src.gc_fraction("ACTGN", ambiguous="ignore"), 0.4)
  assert_eq(@src.gc_fraction("ACTGSSSS", ambiguous="ignore"), 0.75)
  assert_eq(@src.gc_fraction("GDVV", ambiguous="ignore"), 0.25)
}

///|
test "gc_fraction_weighted" {
  assert_eq(@src.gc_fraction("ACTGN", ambiguous="weighted"), 0.5)
  assert_eq(@src.gc_fraction("ACTGSSSS", ambiguous="weighted"), 0.75)
  assert_eq(@src.gc_fraction("GDVV", ambiguous="weighted"), 0.6666666666666666)
}

///|
test "gc_fraction_empty" {
  assert_eq(@src.gc_fraction(""), 0.0)
}

///|
test "gc_fraction_invalid_mode_raises" {
  let result = try {
    ignore(@src.gc_fraction("ACTG", ambiguous="bogus"))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(result)
}

///|
test "gc123_biopython_doctest" {
  let r = @src.gc123("ACTGTN")
  assert_eq(r.gc_all, 40.0)
  assert_eq(r.gc1, 50.0)
  assert_eq(r.gc2, 50.0)
  assert_eq(r.gc3, 0.0)
}

///|
test "gc123_mixed_case" {
  let r = @src.gc123("actgtn")
  assert_eq(r.gc_all, 40.0)
  assert_eq(r.gc1, 50.0)
  assert_eq(r.gc2, 50.0)
  assert_eq(r.gc3, 0.0)
}

///|
test "gc123_incomplete_codon" {
  let r = @src.gc123("ACTG")
  assert_eq(r.gc_all, 50.0)
  assert_eq(r.gc1, 50.0)
  assert_eq(r.gc2, 100.0)
  assert_eq(r.gc3, 0.0)
}

///|
test "gc123_coding_sequence" {
  let r = @src.gc123("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  assert_eq(r.gc_all, 56.41025641025641)
  assert_eq(r.gc1, 53.84615384615385)
  assert_eq(r.gc2, 53.84615384615385)
  assert_eq(r.gc3, 61.53846153846154)
}

///|
test "gc123_empty" {
  let r = @src.gc123("")
  assert_eq(r.gc_all, 0.0)
  assert_eq(r.gc1, 0.0)
  assert_eq(r.gc2, 0.0)
  assert_eq(r.gc3, 0.0)
}

///|
test "nt_search_exact_overlap" {
  let r = @src.nt_search("AAAA", "AA")
  assert_eq(r.pattern, "AA")
  assert_eq(r.positions, [0, 1, 2])
}

///|
test "nt_search_ambiguous" {
  let r = @src.nt_search("GATCGATCGA", "NGA")
  assert_eq(r.pattern, "[GATC]GA")
  assert_eq(r.positions, [3, 7])
}

///|
test "nt_search_no_match" {
  let r = @src.nt_search("GATCGATCGA", "RGAWS")
  assert_eq(r.pattern, "[AG]GA[AT][CG]")
  assert_eq(r.positions, [])
}

///|
test "nt_search_case_sensitive_on_seq" {
  // Biopython builds an uppercase pattern and matches it literally,
  // so a lowercase target sequence produces no hits.
  let r = @src.nt_search("gatcgatcga", "NGA")
  assert_eq(r.pattern, "[GATC]GA")
  assert_eq(r.positions, [])
}

///|
test "nt_search_x_code" {
  // X behaves like N (any nucleotide). Biopython 1.88:
  // nt_search("GATCGA", "XGA") == ['[GATC]GA', 3] -- "GGA" only occurs at 3.
  let r = @src.nt_search("GATCGA", "XGA")
  assert_eq(r.pattern, "[GATC]GA")
  assert_eq(r.positions, [3])
}

///|
test "nt_search_invalid_letter_raises" {
  let result = try {
    ignore(@src.nt_search("GATC", "Z"))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(result)
}

///|
test "seq3_doctest" {
  assert_eq(
    @src.seq3("MAIVMGRWKGAR*"),
    "MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer",
  )
}

///|
test "seq3_custom_map" {
  assert_eq(
    @src.seq3("MAIVMGRWKGAR*", custom_map=Map([("*", "***")])),
    "MetAlaIleValMetGlyArgTrpLysGlyAlaArg***",
  )
  // A provided custom map replaces the default terminator mapping
  assert_eq(@src.seq3("M*", custom_map=Map([])), "MetXaa")
}

///|
test "seq3_undef_code" {
  assert_eq(
    @src.seq3("MAIVMGRWKGA--R*", undef_code="---"),
    "MetAlaIleValMetGlyArgTrpLysGlyAla------ArgTer",
  )
  assert_eq(@src.seq3("MAIVMGRWKGA--R*"), "MetAlaIleValMetGlyArgTrpLysGlyAlaXaaXaaArgTer")
}

///|
test "seq3_extended_iupac" {
  assert_eq(@src.seq3("BXZJUO*"), "AsxXaaGlxXleSecPylTer")
  assert_eq(@src.seq3(""), "")
  // lowercase one-letter codes are not mapped (Biopython parity)
  assert_eq(@src.seq3("j"), "Xaa")
}

///|
test "seq1_doctest" {
  assert_eq(
    @src.seq1("MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer"),
    "MAIVMGRWKGAR*",
  )
}

///|
test "seq1_case_insensitive" {
  assert_eq(
    @src.seq1("METalaIlEValMetGLYArgtRplysGlyAlaARGTer"),
    "MAIVMGRWKGAR*",
  )
}

///|
test "seq1_custom_map" {
  assert_eq(
    @src.seq1("MetAlaIleValMetGlyArgTrpLysGlyAla***", custom_map=Map([("***", "*")])),
    "MAIVMGRWKGA*",
  )
  // A provided custom map replaces the default "Ter" -> "*" mapping
  assert_eq(@src.seq1("MetTer", custom_map=Map([("***", "*")])), "MX")
}

///|
test "seq1_undef_code" {
  assert_eq(
    @src.seq1("MetAlaIleValMetGlyArgTrpLysGlyAla------ArgTer"),
    "MAIVMGRWKGAXXR*",
  )
  assert_eq(
    @src.seq1("MetAlaIleValMetGlyArgTrpLysGlyAla------ArgTer", undef_code="?"),
    "MAIVMGRWKGA??R*",
  )
}

///|
test "seq1_extended_iupac" {
  assert_eq(@src.seq1("AsxGlxXleSecPylTer"), "BZJUO*")
  // Trailing characters that do not form a full triplet are dropped
  assert_eq(@src.seq1("AlaG"), "A")
}

///|
test "six_frame_translations_doctest" {
  // Exact output of Biopython's six_frame_translations doctest
  let expected = "GC_Frame: a:5 t:0 g:8 c:5\nSequence: auggccauug ... gggccgcuga, 24 nt, 54.17 %GC\n\n\n1/1\n  G  H  C  N  G  P  L\n W  P  L  *  W  A  A\nM  A  I  V  M  G  R  *\nauggccauuguaaugggccgcuga   54 %\nuaccgguaacauuacccggcgacu\nA  M  T  I  P  R  Q\n H  G  N  Y  H  A  A  S\n  P  W  Q  L  P  G  S\n\n"
  assert_eq(@src.six_frame_translations("AUGGCCAUUGUAAUGGGCCGCUGA"), expected)
}

///|
test "six_frame_translations_dna_39nt" {
  let expected = "GC_Frame: a:9 t:8 g:14 c:8\nSequence: atggccattg ... tgcccgatag, 39 nt, 56.41 %GC\n\n\n1/1\n  G  H  C  N  G  P  L  K  G  C  P  I\n W  P  L  *  W  A  A  E  R  V  P  D\nM  A  I  V  M  G  R  *  K  G  A  R  *\natggccattgtaatgggccgctgaaagggtgcccgatag   56 %\ntaccggtaacattacccggcgactttcccacgggctatc\nA  M  T  I  P  R  Q  F  P  A  R  Y\n H  G  N  Y  H  A  A  S  L  T  G  S  L\n  P  W  Q  L  P  G  S  F  P  H  G  I\n\n"
  assert_eq(@src.six_frame_translations("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"), expected)
}

///|
test "six_frame_translations_multiblock" {
  // 154 nt sequence wraps into two 60 nt blocks plus a 34 nt block
  let unit = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGTAGGCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
  let out = @src.six_frame_translations(unit + unit)
  // Header
  assert_true(out.has_prefix("GC_Frame: a:36 t:32 g:56 c:30\n"))
  // Block positions: 1/1, 61/21, 121/41
  assert_true(out.contains("\n61/21\n"))
  assert_true(out.contains("\n121/41\n"))
  // Each block shows the nucleotide line followed by the reverse complement
  assert_true(out.contains("atggccattgtaatgggccgctgaaagggtgcccgatagtaggcattgtaatgggccgct   55 %\n"))
  assert_true(out.contains("taccggtaacattacccggcgactttcccacgggctatcatccgtaacattacccggcga\n"))
}

///|
test "six_frame_translations_rna_vs_dna" {
  // Any 'u' (even mixed with T) switches to the RNA complement tables
  // (A <-> U), mirroring Biopython 1.88's `if "u" in seq.lower()` rule.
  let rna = @src.six_frame_translations("AUGGCCAUUGUAAUGGGCCGCUGA")
  assert_true(rna.contains("GC_Frame: a:5 t:0 g:8 c:5\n"))
  assert_true(rna.contains("auggccauuguaaugggccgcuga   54 %\n"))
  assert_true(rna.contains("uaccgguaacauuacccggcgacu\n"))
  // A mixed T/U sequence still uses the RNA tables; only the counts differ.
  let mixed = @src.six_frame_translations("ATGGCCAUUGUAAUGGGCCGCUGA")
  assert_true(mixed.contains("GC_Frame: a:5 t:1 g:8 c:5\n"))
  assert_true(mixed.contains("atggccauuguaaugggccgcuga   54 %\n"))
  assert_true(mixed.contains("uaccgguaacauuacccggcgacu\n"))
}

///|
test "six_frame_translations_short_seq" {
  // Very short sequence: the +2/+3 and -2/-3 frames are empty but the
  // layout is preserved (exact output verified against Biopython 1.88).
  let out = @src.six_frame_translations("ATG")
  assert_eq(
    out,
    "GC_Frame: a:1 t:1 g:1 c:0\nSequence: atg, 3 nt, 33.33 %GC\n\n\n1/1\n  \n \nM\natg   33 %\ntac\n\n H\n  \n\n",
  )
}

///|
test "six_frame_translations_bad_genetic_code_raises" {
  let result = try {
    ignore(@src.six_frame_translations("ATG", genetic_code=11))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(result)
}
