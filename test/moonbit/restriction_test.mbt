///|
test "restriction_builtin_enzymes" {
  let enzymes = @src.builtin_enzymes()

  assert_true(enzymes.keys().collect().length() > 0)
}

///|
test "restriction_get_enzyme_by_name" {
  let enzyme = @src.get_enzyme_by_name("EcoRI")

  assert_eq(enzyme.is_some(), true)
}

///|
test "restriction_find_cut_sites" {
  let enzyme = @src.restriction_enzyme("EcoRI", "GAATTC", 1, "AATT", "")
  let sites = @src.find_cut_sites("GAATTCGAATTC", enzyme)

  assert_eq(sites.length(), 2)
}

///|
test "restriction_digest" {
  let enzyme = @src.restriction_enzyme("EcoRI", "GAATTC", 1, "AATT", "")
  let fragments = @src.digest("GAATTCGAATTC", enzyme)

  assert_eq(fragments.length(), 3)
}

///|
test "restriction_digest_with_multiple" {
  let enzyme1 = @src.restriction_enzyme("EcoRI", "GAATTC", 1, "AATT", "")
  let enzyme2 = @src.restriction_enzyme("HindIII", "AAGCTT", 1, "AGCT", "")
  let fragments = @src.digest_with_multiple("GAATTCAAGCTT", [enzyme1, enzyme2])

  assert_eq(fragments.length(), 3)
}

///|
test "restriction_fragment_sizes" {
  let fragments = ["ACGT", "ACAT"]
  let sizes = @src.fragment_sizes(fragments)

  assert_eq(sizes[0], 4)
}

///|
test "restriction_is_compatible_overhang" {
  let enzyme1 = @src.restriction_enzyme("EcoRI", "GAATTC", 1, "AATT", "")
  let enzyme2 = @src.restriction_enzyme("EcoRI", "GAATTC", 1, "AATT", "")

  assert_true(@src.is_compatible_overhang(enzyme1, enzyme2))
}

///|
test "restriction_list_enzymes" {
  let list = @src.list_enzymes()

  assert_true(list.length() > 0)
}

///|
test "restriction_enzymes_by_overhang" {
  let enzymes = @src.enzymes_by_overhang("AATT")

  assert_true(enzymes.length() > 0)
}

///|
test "restriction_enzymes_by_site_length" {
  let enzymes = @src.enzymes_by_site_length(6, 6)

  assert_true(enzymes.length() > 0)
}

///|
fn br_enzyme(name : String) -> @src.BioRestrictionEnzyme {
  @src.get_bio_restriction_enzyme(name).unwrap()
}

///|
fn br_sites(m : Map[String, Array[Int]], key : String) -> Array[Int] {
  m.get(key).unwrap()
}

///|
fn br_search(
  e : @src.BioRestrictionEnzyme,
  seq : String,
  linear : Bool,
) -> Array[Int] {
  e.search(seq, linear~) catch {
    _ => abort("search should not raise here")
  }
}

///|
fn br_catalyse(
  e : @src.BioRestrictionEnzyme,
  seq : String,
  linear : Bool,
) -> Array[String] {
  e.catalyse(seq, linear~) catch {
    _ => abort("catalyse should not raise here")
  }
}

///|
test "bio_restriction_database_count" {
  assert_eq(@src.bio_restriction_enzyme_count(), 1088)
  assert_eq(@src.restriction_db_size, 1088)
  assert_true(@src.get_bio_restriction_enzyme("EcoRI").is_some())
  assert_true(@src.get_bio_restriction_enzyme("Nope123").is_none())
  let raised = try {
    ignore(@src.require_bio_restriction_enzyme("Nope123"))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(raised)
}

///|
test "bio_restriction_search_basic" {
  let e0 = br_enzyme("EcoRI")
  assert_eq(br_search(e0, "GAATTC", true), [2])
  let e1 = br_enzyme("EcoRI")
  assert_eq(br_search(e1, "AAGAATTCGGGAATTCAA", true), [4, 12])
  let e2 = br_enzyme("SmaI")
  assert_eq(br_search(e2, "CCCGGG", true), [4])
  let e3 = br_enzyme("KpnI")
  assert_eq(br_search(e3, "GGTACC", true), [6])
  let e4 = br_enzyme("BbsI")
  assert_eq(br_search(e4, "AAGAAGACAAGGGGGGGG", true), [11])
  let e5 = br_enzyme("BspMI")
  assert_eq(br_search(e5, "AAACCTGCAAGGGGGGG", true), [13])
  let e6 = br_enzyme("BstXI")
  assert_eq(br_search(e6, "CCAAAAAAATGG", true), [9])
  let e7 = br_enzyme("BstYI")
  assert_eq(br_search(e7, "AGATCT", true), [2])
  let e8 = br_enzyme("BstYI")
  assert_eq(br_search(e8, "GGATCC", true), [2])
  let e9 = br_enzyme("AvaI")
  assert_eq(br_search(e9, "CTCGAG", true), [2])
  let e10 = br_enzyme("AvaI")
  assert_eq(br_search(e10, "CCCGGG", true), [2])
  assert_eq(br_search(br_enzyme("EcoRI"), "gaattc", true), [2])
  assert_eq(br_search(br_enzyme("EcoRI"), "G AA TTC", true), [2])
  assert_eq(br_search(br_enzyme("EcoRI"), "G2AATTC", true), [2])
  assert_eq(br_search(br_enzyme("EcoRI"), "GAATXC", true), [])
  let raised = try {
    ignore(br_enzyme("EcoRI").search("GAA.TC"))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(raised)
}

///|
test "bio_restriction_search_reverse_strand" {
  // reverse-complement hit with cut inside
  assert_eq(br_search(br_enzyme("BbsI"), "AAAAAAAAAAGTCTTCAA", true), [5])
  // reverse-complement site near the start wraps only on circular
  assert_eq(br_search(br_enzyme("BbsI"), "TTGTCTTCGGGGG", true), [])
  assert_eq(br_search(br_enzyme("BbsI"), "TTGTCTTCGGGGG", false), [10])
  assert_eq(br_search(br_enzyme("BbsI"), "GAAGACAAGG", true), [])
  assert_eq(br_search(br_enzyme("BbsI"), "GAAGACAAGG", false), [9])
}

///|
test "bio_restriction_search_circular" {
  assert_eq(br_search(br_enzyme("BbsI"), "GTCTTCAAAAAAAAAAAAAAAAAAAA", false), [
    21,
  ])
  assert_eq(br_search(br_enzyme("BbsI"), "GTCTTCAAAAAAAAAAAAAAAAAAAA", true), [])
  assert_eq(br_search(br_enzyme("BaeI"), "AACAAAAGTACCA", false), [5, 12])
  assert_eq(br_search(br_enzyme("AbaSI"), "AAAACAAAAAAAAAAA", false), [1])
  assert_eq(br_search(br_enzyme("AbaSI"), "AAAACAAAAAAAAAAA", true), [])
}

///|
test "bio_restriction_search_twocuts" {
  assert_eq(
    br_search(br_enzyme("BcgI"), "AAAAAAAAAAAACGAAAAAAATGCAAAAAAAAAAAAA", true),
    [37],
  )
  assert_eq(
    br_search(br_enzyme("BcgI"), "AAAAAAAAAAAACGAAAAAAATGCAAAAAAAAAAAAA", false),
    [3, 37],
  )
  assert_eq(
    br_search(br_enzyme("BaeI"), "AAAAAAAAAAAAACAAAAGTACCAAAAAAAAAAAAAA", true),
    [36],
  )
  assert_true(br_enzyme("BcgI").cut_twice())
  assert_true(br_enzyme("EcoRI").cut_once())
}

///|
test "bio_restriction_catalyse" {
  let f0 = br_catalyse(br_enzyme("EcoRI"), "GAATTC", true)
  assert_eq(f0, ["G", "AATTC"])
  let f1 = br_catalyse(br_enzyme("EcoRI"), "AAGAATTCGGGAATTCAA", true)
  assert_eq(f1, ["AAG", "AATTCGGG", "AATTCAA"])
  let f2 = br_catalyse(br_enzyme("SmaI"), "CCCGGG", true)
  assert_eq(f2, ["CCC", "GGG"])
  let f3 = br_catalyse(br_enzyme("KpnI"), "GGTACC", true)
  assert_eq(f3, ["GGTAC", "C"])
  let f4 = br_catalyse(br_enzyme("BbsI"), "AAGAAGACAAGGGGGGGG", true)
  assert_eq(f4, ["AAGAAGACAA", "GGGGGGGG"])
  let f5 = br_catalyse(br_enzyme("EcoRI"), "GAATTCAAGAATTCAA", false)
  assert_eq(f5, ["AATTCAAG", "AATTCAAG"])
  let f6 = br_catalyse(br_enzyme("BbsI"), "GTCTTCAAAAAAAAAAAAAAAAAAAA", false)
  assert_eq(f6, ["AAAAAAGTCTTCAAAAAAAAAAAAAA"])
  let f7 = br_catalyse(br_enzyme("EcoRI"), "CCCCCCCC", true)
  assert_eq(f7, ["CCCCCCCC"])
  let f8 = br_catalyse(br_enzyme("BbsI"), "GAAGACAAGG", false)
  assert_eq(f8, ["GGGAAGACAA"])
  assert_eq(
    @src.bio_restriction_search(br_enzyme("EcoRI"), "GAATTC") catch {
      _ => abort("x")
    },
    [2],
  )
  assert_eq(
    @src.bio_restriction_catalyse(br_enzyme("SmaI"), "CCCGGG") catch {
      _ => abort("x")
    },
    ["CCC", "GGG"],
  )
  let raised = try {
    ignore(br_enzyme("Aba13301I").catalyse("GCAAAC"))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(raised)
  let raised = try {
    ignore(br_enzyme("HpyUM037X").catalyse("TNGGNAGGTGGNAG"))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(raised)
  assert_eq(br_search(br_enzyme("Aba13301I"), "GCAAAC", true), [1])
}

///|
test "bio_restriction_attributes" {
  assert_eq(br_enzyme("EcoRI").site, "GAATTC")
  assert_eq(br_enzyme("EcoRI").size, 6)
  assert_true(br_enzyme("EcoRI").is_palindromic())
  assert_true(!br_enzyme("BbsI").is_palindromic())
  assert_true((br_enzyme("EcoRI").frequency() - 4096.0).abs() < 1.0e-9)
  assert_true(br_enzyme("EcoRI").is_5overhang())
  assert_true(!br_enzyme("EcoRI").is_3overhang())
  assert_true(br_enzyme("KpnI").is_3overhang())
  assert_true(br_enzyme("SmaI").is_blunt())
  assert_true(br_enzyme("Aba13301I").is_unknown_end())
  assert_true(br_enzyme("EcoRI").is_defined())
  assert_true(br_enzyme("BbsI").is_ambiguous())
  assert_true(br_enzyme("Aba13301I").is_not_defined())
  assert_eq(br_enzyme("EcoRI").overhang_type(), "5' overhang")
  assert_eq(br_enzyme("SmaI").overhang_type(), "blunt")
  assert_eq(br_enzyme("KpnI").overhang_type(), "3' overhang")
  assert_eq(br_enzyme("Aba13301I").overhang_type(), "unknown")
  let ch = br_enzyme("EcoRI").characteristic()
  assert_eq(ch.0, Some(1))
  assert_eq(ch.1, Some(-1))
  assert_true(ch.2 is None)
  assert_true(ch.3 is None)
  assert_eq(ch.4, "GAATTC")
  assert_eq(br_enzyme("EcoRI").supplier_list(), [
    "B", "C", "I", "J", "K", "M", "N", "O", "Q", "R", "S", "V", "X", "Y",
  ])
  assert_true(br_enzyme("EcoRI").is_commercial())
  assert_true(!br_enzyme("Aba13301I").is_commercial())
  assert_true(
    (br_enzyme("AspNIH4III").frequency() - 5461.333333333333).abs() < 1.0e-9,
  )
  assert_eq(@src.bio_restriction_enzymes_by_supplier("N").length(), 238)
}

///|
test "bio_restriction_elucidate" {
  assert_eq(br_enzyme("EcoRI").elucidate(), "G^AATT_C")
  assert_eq(br_enzyme("KpnI").elucidate(), "G_GTAC^C")
  assert_eq(br_enzyme("SmaI").elucidate(), "CCC^_GGG")
  assert_eq(br_enzyme("BstYI").elucidate(), "R^GATC_Y")
  assert_eq(
    br_enzyme("BcgI").elucidate(),
    "cut twice, not yet implemented sorry.",
  )
  assert_eq(
    br_enzyme("BaeI").elucidate(),
    "cut twice, not yet implemented sorry.",
  )
  assert_eq(br_enzyme("Aba13301I").elucidate(), "? GCAAAC ?")
}

///|
test "bio_restriction_compatible_overhang" {
  assert_eq(br_enzyme("XhoI").compatible_overhang(br_enzyme("SalI")), true) // XhoI % SalI
  assert_eq(br_enzyme("SmaI").compatible_overhang(br_enzyme("EcoRV")), true) // SmaI % EcoRV
  assert_eq(br_enzyme("EcoRI").compatible_overhang(br_enzyme("KpnI")), false) // EcoRI % KpnI
  assert_eq(br_enzyme("BamHI").compatible_overhang(br_enzyme("BglII")), true) // BamHI % BglII
  assert_eq(br_enzyme("SmaI").compatible_overhang(br_enzyme("EcoRI")), false) // SmaI % EcoRI
  assert_eq(
    br_enzyme("Aba13301I").compatible_overhang(br_enzyme("EcoRI")),
    false,
  ) // Aba13301I % EcoRI
  assert_eq(br_enzyme("BspMI").compatible_overhang(br_enzyme("BspLU11I")), true) // BspMI % BspLU11I
  assert_eq(br_enzyme("BspLU11I").compatible_overhang(br_enzyme("BspMI")), true) // BspLU11I % BspMI
  assert_eq(br_enzyme("EcoRI").compatible_ends().length(), 51)
  assert_true(br_enzyme("EcoRI").compatible_ends().contains(br_enzyme("MfeI")))
  assert_true(br_enzyme("EcoRI").compatible_ends().contains(br_enzyme("BbsI")))
  assert_true(!br_enzyme("EcoRI").compatible_ends().contains(br_enzyme("XhoI")))
}

///|
test "bio_restriction_isoschizomers" {
  assert_true(br_enzyme("SacI").is_equischizomer(br_enzyme("SstI")))
  assert_true(!br_enzyme("SacI").is_neoschizomer(br_enzyme("SstI")))
  assert_true(br_enzyme("SacI").is_isoschizomer(br_enzyme("SstI")))
  assert_true(!br_enzyme("SmaI").is_equischizomer(br_enzyme("XmaI")))
  assert_true(br_enzyme("SmaI").is_neoschizomer(br_enzyme("XmaI")))
  assert_true(br_enzyme("SmaI").is_isoschizomer(br_enzyme("XmaI")))
  assert_true(!br_enzyme("EcoRI").is_isoschizomer(br_enzyme("BamHI")))
  assert_true(br_enzyme("SacI").equischizomers().contains(br_enzyme("SstI")))
  assert_true(br_enzyme("SmaI").neoschizomers().contains(br_enzyme("XmaI")))
  assert_true(br_enzyme("SmaI").isoschizomers().contains(br_enzyme("XmaI")))
}

///|
test "bio_restriction_batch" {
  let rb = @src.RestrictionBatch::new(names=[
    "EcoRI", "SmaI", "KpnI", "BamHI", "BbsI", "BstYI", "AvaI", "Aba13301I",
  ]) catch {
    _ => abort("batch")
  }
  assert_eq(rb.length(), 8)
  assert_true(rb.contains("EcoRI"))
  assert_true(!rb.contains("XhoI"))
  rb.add("XhoI") catch {
    _ => abort("add")
  }
  assert_eq(rb.length(), 9)
  let raised = try {
    ignore(rb.add("Nope123"))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(raised)
  rb.remove("XhoI")
  assert_eq(rb.length(), 8)
  assert_eq(rb.elements(), [
    "Aba13301I", "AvaI", "BamHI", "BbsI", "BstYI", "EcoRI", "KpnI", "SmaI",
  ])
  rb.add_supplier("N")
  assert_eq(rb.length(), 239)
  let bseq = "AAGAATTCAACCCGGGAAGGATCCAAGAAGACGGGGGGGGCC"
  let bs = rb.search(bseq) catch { _ => abort("search") }
  assert_eq(br_sites(bs, "EcoRI"), [4])
  assert_eq(br_sites(bs, "SmaI"), [14])
  assert_eq(br_sites(bs, "KpnI"), [])
  assert_eq(br_sites(bs, "BamHI"), [20])
  assert_eq(br_sites(bs, "BbsI"), [35])
  assert_eq(br_sites(bs, "BstYI"), [20])
  assert_eq(br_sites(bs, "AvaI"), [12])
  assert_eq(br_sites(bs, "Aba13301I"), [])
}

///|
test "bio_restriction_analysis" {
  let rb3 = @src.RestrictionBatch::new(names=[
    "EcoRI", "SmaI", "KpnI", "BamHI", "BbsI", "BstYI", "AvaI", "Aba13301I",
  ]) catch {
    _ => abort("batch")
  }
  let an = @src.RestrictionAnalysis::new(
    batch=rb3,
    seq="AAGAATTCAACCCGGGAAGGATCCAAGAAGACGGGGGGGGCC",
  ) catch {
    _ => abort("analysis")
  }
  assert_eq(an.sequence_length, 42)
  assert_eq(an.full().length(), 8)
  assert_eq(an.with_sites().length(), 6)
  assert_eq(an.without_site().length(), 2)
  assert_eq(an.with_n_sites(1).length(), 6)
  assert_eq(an.blunt().length(), 1)
  assert_eq(an.overhang5().length(), 5)
  assert_eq(an.overhang3().length(), 1)
  assert_eq(an.defined().length(), 5)
  assert_eq(an.ambiguous().length(), 2)
  assert_true(an.blunt().contains("SmaI"))
  assert_true(an.ambiguous().contains("AvaI"))
  assert_true(an.ambiguous().contains("BbsI"))
  assert_eq(br_sites(an.full(), "EcoRI"), [4])
  assert_eq(br_sites(an.full(), "BbsI"), [35])
  // range 1..25
  assert_eq(an.only_between(1, 25).length(), 5)
  assert_eq(an.between(1, 25).length(), 5)
  assert_eq(an.show_only_between(1, 25).length(), 5)
  assert_eq(an.only_outside(1, 25).length(), 1)
  assert_eq(an.outside(1, 25).length(), 1)
  assert_eq(an.do_not_cut(1, 25).length(), 3)
  assert_eq(br_sites(an.show_only_between(1, 25), "EcoRI"), [4])
  assert_eq(br_sites(an.show_only_between(1, 25), "SmaI"), [14])
  assert_eq(br_sites(an.show_only_between(1, 25), "BstYI"), [20])
  assert_eq(br_sites(an.show_only_between(1, 25), "AvaI"), [12])
  assert_eq(br_sites(an.show_only_between(1, 25), "BamHI"), [20])
  // range 20..40
  assert_eq(an.only_between(20, 40).length(), 3)
  assert_eq(an.between(20, 40).length(), 3)
  assert_eq(an.show_only_between(20, 40).length(), 3)
  assert_eq(an.only_outside(20, 40).length(), 3)
  assert_eq(an.outside(20, 40).length(), 3)
  assert_eq(an.do_not_cut(20, 40).length(), 5)
  assert_eq(br_sites(an.show_only_between(20, 40), "BstYI"), [20])
  assert_eq(br_sites(an.show_only_between(20, 40), "BbsI"), [35])
  assert_eq(br_sites(an.show_only_between(20, 40), "BamHI"), [20])
  // range 1..15
  assert_eq(an.only_between(1, 15).length(), 3)
  assert_eq(an.between(1, 15).length(), 3)
  assert_eq(an.show_only_between(1, 15).length(), 3)
  assert_eq(an.only_outside(1, 15).length(), 3)
  assert_eq(an.outside(1, 15).length(), 3)
  assert_eq(an.do_not_cut(1, 15).length(), 5)
  assert_eq(br_sites(an.show_only_between(1, 15), "EcoRI"), [4])
  assert_eq(br_sites(an.show_only_between(1, 15), "SmaI"), [14])
  assert_eq(br_sites(an.show_only_between(1, 15), "AvaI"), [12])
  // range -10..-1
  assert_eq(an.only_between(-10, -1).length(), 1)
  assert_eq(an.between(-10, -1).length(), 1)
  assert_eq(an.show_only_between(-10, -1).length(), 1)
  assert_eq(an.only_outside(-10, -1).length(), 5)
  assert_eq(an.outside(-10, -1).length(), 5)
  assert_eq(an.do_not_cut(-10, -1).length(), 7)
  assert_eq(br_sites(an.show_only_between(-10, -1), "BbsI"), [])
  // range 40..20
  assert_eq(an.only_between(40, 20).length(), 3)
  assert_eq(an.between(40, 20).length(), 3)
  assert_eq(an.show_only_between(40, 20).length(), 3)
  assert_eq(an.only_outside(40, 20).length(), 3)
  assert_eq(an.outside(40, 20).length(), 3)
  assert_eq(an.do_not_cut(40, 20).length(), 5)
  assert_eq(br_sites(an.show_only_between(40, 20), "BstYI"), [20])
  assert_eq(br_sites(an.show_only_between(40, 20), "BbsI"), [])
  assert_eq(br_sites(an.show_only_between(40, 20), "BamHI"), [20])
}
