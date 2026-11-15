///|
/// Tests for bioframe_assemblies.mbt — curated genome metadata.

// ---- bf_assemblies_available ----------------------------------------------

test "bf_assemblies_available: 9 assemblies" {
  let all = @src.bf_assemblies_available()
  assert_eq(all.length(), 9)
  assert_eq(all[0].name, "hg19")
  assert_eq(all[1].name, "hg38")
  assert_eq(all[2].name, "hs1")
}

// ---- bf_lookup_seqinfo ---------------------------------------------------

test "bf_lookup_seqinfo: hg38 succeeds with canonical seqinfo" {
  match @src.bf_lookup_seqinfo("hg38") {
    Some(infos) => {
      assert_true(infos.length() >= 25)
      // chr1
      assert_eq(infos[0].name, "chr1")
      assert_eq(infos[0].length, 248956422L)
    }
    None => assert_true(false)
  }
}

test "bf_lookup_seqinfo: unknown assembly → None" {
  match @src.bf_lookup_seqinfo("notAgenome") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

test "bf_lookup_seqinfo: hg19 chr1 length correct" {
  match @src.bf_lookup_seqinfo("hg19") {
    Some(infos) => {
      assert_eq(infos[0].name, "chr1")
      assert_eq(infos[0].length, 249250621L)
    }
    None => assert_true(false)
  }
}

// ---- bf_assembly_chromsizes -----------------------------------------------

test "bf_assembly_chromsizes: hg38 contains canonical chroms" {
  let cs = @src.bf_assembly_chromsizes("hg38")
  assert_true(cs.length() >= 25)
  let names = ["chr1", "chr22", "chrX", "chrY", "chrM"]
  for n in names {
    let mut found = false
    for pair in cs {
      if pair.0 == n { found = true; break }
    }
    assert_true(found)
  }
}

test "bf_assembly_chromsizes: unknown assembly → empty" {
  let cs = @src.bf_assembly_chromsizes("fake")
  assert_eq(cs.length(), 0)
}

// ---- bf_assembly_centromeres ----------------------------------------------

test "bf_assembly_centromeres: hg19 → 24 rows (1-22 + X + Y)" {
  let cents = @src.bf_assembly_centromeres("hg19")
  assert_eq(cents.length(), 24)
  assert_eq(cents[0].0, "chr1")
  assert_eq(cents[1].0, "chr2")
}

test "bf_assembly_centromeres: unknown assembly → empty" {
  let cents = @src.bf_assembly_centromeres("fake")
  assert_eq(cents.length(), 0)
}

test "bf_assembly_centromeres: mm10 → 21 rows" {
  let cents = @src.bf_assembly_centromeres("mm10")
  assert_eq(cents.length(), 21)
}

// ---- bf_assembly_cytobands ------------------------------------------------

test "bf_assembly_cytobands: hg38 → cytobands head non-empty" {
  let bands = @src.bf_assembly_cytobands("hg38")
  assert_true(bands.length() > 0)
  assert_eq(bands[0].chrom, "chr1")
  assert_eq(bands[0].name, "p36.33")
}

test "bf_assembly_cytobands: unknown assembly → empty" {
  let bands = @src.bf_assembly_cytobands("fake")
  assert_eq(bands.length(), 0)
}

// ---- bf_resolve_alias -----------------------------------------------------

test "bf_resolve_alias: MT → chrM (hg38)" {
  match @src.bf_resolve_alias("hg38", "MT") {
    Some(c) => assert_eq(c, "chrM")
    None => assert_true(false)
  }
}

test "bf_resolve_alias: 1 → chr1 (hg38)" {
  match @src.bf_resolve_alias("hg38", "1") {
    Some(c) => assert_eq(c, "chr1")
    None => assert_true(false)
  }
}

test "bf_resolve_alias: X → chrX (hg38)" {
  match @src.bf_resolve_alias("hg38", "X") {
    Some(c) => assert_eq(c, "chrX")
    None => assert_true(false)
  }
}

test "bf_resolve_alias: already canonical → returns as-is" {
  match @src.bf_resolve_alias("hg38", "chr1") {
    Some(c) => assert_eq(c, "chr1")
    None => assert_true(false)
  }
}

test "bf_resolve_alias: no assembly → None" {
  match @src.bf_resolve_alias("nonexistent", "MT") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

test "bf_resolve_alias: name not found → None" {
  match @src.bf_resolve_alias("hg38", "NOTCHROM") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

// ---- bf_assembly_genome_info ----------------------------------------------

test "bf_assembly_genome_info: hg38 → Some with chroms" {
  match @src.bf_assembly_genome_info("hg38") {
    Some(db) => {
      assert_true(db.seqnames().length() > 0)
    }
    None => assert_true(false)
  }
}

test "bf_assembly_genome_info: unknown → None" {
  match @src.bf_assembly_genome_info("fake") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}
