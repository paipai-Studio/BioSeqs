///|
/// Tests for mirna_targets module (miRNA target prediction).

test "mir_extract_seed_basic" {
  // miRNA: "ATGCATGC" (8 nt)
  // seed(2-8): positions 1-7 = "TGCATGC"
  let seed = @src.mir_extract_seed("ATGCATGC")
  assert_eq(seed, "TGCATGC")
  // seed6(2-7): positions 1-6 = "TGCATG"
  let seed6 = @src.mir_extract_seed6("ATGCATGC")
  assert_eq(seed6, "TGCATG")
}

test "mir_extract_seed_short" {
  let seed = @src.mir_extract_seed("ACGT")
  assert_eq(seed, "CGT")
}

test "mir_find_seed_matches_8mer" {
  // miRNA: "ATGCATGC" -> seed = "TGCATGC", rc = "GCATGCA"
  // UTR: "AAGCATGCA" -> rc at pos 1, preceded by A -> 8mer
  let mirna = "ATGCATGC"
  let utr = "AAGCATGCA"
  let sites = @src.mir_find_seed_matches(mirna, utr, min_seed_type=2)
  assert_eq(sites.length(), 1)
  assert_eq(@src.mir_seed_type_label(sites[0].seed_type), "8mer")
}

test "mir_find_seed_matches_7mer_m8" {
  // miRNA: "ATGCATGC" -> seed = "TGCATGC", rc = "GCATGCA"
  // UTR: "TTGCATGCA" -> rc at pos 1, preceded by T -> 7mer-m8

  let mirna = "ATGCATGC"
  let utr = "TTGCATGCA"
  let sites = @src.mir_find_seed_matches(mirna, utr, min_seed_type=2)
  assert_eq(sites.length(), 1)
  assert_eq(@src.mir_seed_type_label(sites[0].seed_type), "7mer-m8")
}

test "mir_find_seed_matches_no_match" {
  let sites = @src.mir_find_seed_matches("ATGCATGC", "AAAAAAAAAA")
  assert_eq(sites.length(), 0)
}

test "mir_predict_targets_basic" {
  let mirnas = [("miR-A", "ATGCATGC")]
  let utrs = [("GENE1", "ATGCATGCA")]
  let results = @src.mir_predict_targets(mirnas, utrs)
  assert_true(results.length() > 0)
  assert_eq(results[0].mir_id, "miR-A")
}

test "mir_export_sites" {
  let mirnas = [("miR-B", "ATGCATGC")]
  let utrs = [("TEST1", "ATGCATGCA")]
  let results = @src.mir_predict_targets(mirnas, utrs)
  assert_true(results.length() > 0)
  let tsv = @src.mir_export_sites(results)
  assert_true(tsv.contains("mir_id"))
  assert_true(tsv.contains("miR-B"))
}
