///|
test "genome_info_db_basic_create" {
  let db = @src.GenomeInfoDb::new(genome_build="hg38", organism="Homo sapiens")
  let db2 = db.add_chrom("chr1", 248956422).add_chrom("chr2", 242193529)
  assert_eq(@src.GenomeInfoDb::num_chromosomes(db2), 2)
  assert_eq(@src.GenomeInfoDb::genome_length(db2), 491149951)
}

///|
test "genome_info_db_prebuilt_hg38" {
  let db = @src.genome_hg38()
  // hg38 should have 25 chromosomes (1-22, X, Y, M)
  assert_eq(@src.GenomeInfoDb::num_chromosomes(db), 25)
  // chr1 length is 248956422 in hg38
  match @src.GenomeInfoDb::seqLength(db, "chr1") {
    Some(len) => assert_eq(len, 248956422)
    None => assert_true(false)
  }
  // chrM is circular
  match @src.GenomeInfoDb::is_circular(db, "chrM") {
    Some(c) => assert_true(c)
    None => assert_true(false)
  }
  // chr1 should not be circular
  match @src.GenomeInfoDb::is_circular(db, "chr1") {
    Some(c) => assert_false(c)
    None => assert_true(false)
  }
}

///|
test "genome_info_db_prebuilt_hg19" {
  let db = @src.genome_hg19()
  assert_eq(@src.GenomeInfoDb::num_chromosomes(db), 25)
  // chr1 length is 249250621 in hg19
  match @src.GenomeInfoDb::seqLength(db, "chr1") {
    Some(len) => assert_eq(len, 249250621)
    None => assert_true(false)
  }
}

///|
test "genome_info_db_prebuilt_mm10" {
  let db = @src.genome_mm10()
  assert_eq(@src.GenomeInfoDb::num_chromosomes(db), 22)
  // chr1 length is 195471971 in mm10
  match @src.GenomeInfoDb::seqLength(db, "chr1") {
    Some(len) => assert_eq(len, 195471971)
    None => assert_true(false)
  }
}

///|
test "genome_info_db_centromere" {
  let db = @src.genome_hg38()
  // chr1 centromere in hg38
  match @src.GenomeInfoDb::centromere(db, "chr1") {
    Some(c) => {
      assert_eq(c.start, 121535434)
      assert_eq(c.end, 124595842)
    }
    None => assert_true(false)
  }
  // hg19 has no centromere data
  let db19 = @src.genome_hg19()
  match @src.GenomeInfoDb::centromere(db19, "chr1") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "genome_info_db_get_arm" {
  let db = @src.genome_hg38()
  // chr1 position before centromere is p arm
  assert_eq(@src.GenomeInfoDb::get_arm(db, "chr1", 1000000), "p")
  // chr1 position after centromere is q arm
  assert_eq(@src.GenomeInfoDb::get_arm(db, "chr1", 200000000), "q")
  // chr1 position within centromere
  assert_eq(@src.GenomeInfoDb::get_arm(db, "chr1", 123000000), "centromere")
}

///|
test "genome_info_db_seqnames" {
  let db = @src.genome_hg38()
  let names = @src.GenomeInfoDb::seqnames(db)
  assert_eq(names.length(), 25)
  assert_eq(names[0], "chr1")
  assert_eq(names[24], "chrM")
}

///|
test "genome_info_db_standard_chromosomes" {
  let db = @src.GenomeInfoDb::new(genome_build="hg38", organism="Homo sapiens")
  let db2 = db
    .add_chrom("chr1", 1000)
    .add_chrom("chr2", 2000)
    .add_chrom("chr1_random", 500)
    .add_chrom("chrUn_unplaced", 800)
    .add_chrom("chr1_KI270706v1_alt", 300)
  let std = @src.GenomeInfoDb::standard_chromosomes(db2)
  let non_std = @src.GenomeInfoDb::non_standard_chromosomes(db2)
  assert_eq(std.length(), 2)
  assert_eq(non_std.length(), 3)
}

///|
test "genome_info_db_filter_chromosomes" {
  let db = @src.genome_hg38()
  let chr1_only = @src.GenomeInfoDb::filter_chromosomes(db, "chr1")
  // Should match chr1 only (not chr10, chr11, etc. since contains "chr1")
  assert_true(chr1_only.length() >= 1)
  assert_true(chr1_only.contains("chr1"))
}

///|
test "genome_info_db_organism_from_build" {
  assert_eq(@src.organism_from_build("hg38"), "Homo sapiens")
  assert_eq(@src.organism_from_build("hg19"), "Homo sapiens")
  assert_eq(@src.organism_from_build("mm10"), "Mus musculus")
  assert_eq(@src.organism_from_build("rn6"), "Rattus norvegicus")
  assert_eq(@src.organism_from_build("danRer11"), "Danio rerio")
  assert_eq(@src.organism_from_build("dm6"), "Drosophila melanogaster")
  assert_eq(@src.organism_from_build("ce11"), "Caenorhabditis elegans")
  assert_eq(@src.organism_from_build("unknown_build"), "Unknown")
}

///|
test "genome_info_db_supported_builds" {
  let builds = @src.supported_genome_builds()
  assert_true(builds.contains("hg38"))
  assert_true(builds.contains("hg19"))
  assert_true(builds.contains("mm10"))
}

///|
test "genome_info_db_compare_builds" {
  assert_true(@src.compare_genome_builds("hg38", "hg38"))
  assert_false(@src.compare_genome_builds("hg38", "hg19"))
}

///|
test "genome_info_db_summary" {
  let db = @src.genome_hg38()
  let s = @src.GenomeInfoDb::summary(db)
  assert_true(s.contains("hg38"))
  assert_true(s.contains("Homo sapiens"))
  assert_true(s.contains("25 chromosomes"))
}

///|
test "genome_info_db_genome_length" {
  let db = @src.GenomeInfoDb::new(genome_build="test", organism="Test")
  let db2 = db.add_chrom("chr1", 1000).add_chrom("chr2", 2000)
  assert_eq(@src.GenomeInfoDb::genome_length(db2), 3000)
}
