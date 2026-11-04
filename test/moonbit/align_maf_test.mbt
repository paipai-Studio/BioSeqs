///|
fn align_maf_test_parse(text : String) -> @src.AlignMafDocument {
  @src.align_maf_parse(text) catch {
    AlignMafError(message) => abort("valid MAF input failed: " + message)
  }
}

///|
fn align_maf_test_document() -> @src.AlignMafDocument {
  align_maf_test_parse(@src.align_maf_example_text())
}

///|
fn align_maf_test_index() -> @src.AlignMafIndex {
  @src.AlignMafIndex::create(align_maf_test_document(), "hg38.chr7") catch {
    AlignMafError(message) => abort("valid MAF index failed: " + message)
  }
}

///|
fn align_maf_test_rejects(text : String) -> Bool {
  try {
    ignore(@src.align_maf_parse(text))
    false
  } catch {
    AlignMafError(_) => true
  }
}

///|
fn align_maf_test_single(lines : String) -> String {
  "##maf version=1\n\n" + lines
}

///|
test "Bio.Align.maf parses the MAF version" {
  assert_eq(align_maf_test_document().version, "1")
}

///|
test "Bio.Align.maf parses scoring metadata" {
  assert_eq(align_maf_test_document().scoring, Some("blastz"))
}

///|
test "Bio.Align.maf parses program metadata" {
  assert_eq(align_maf_test_document().program, Some("multiz"))
}

///|
test "Bio.Align.maf parses header comments" {
  let document = align_maf_test_document()
  assert_eq(document.comments, ["generated fixture"])
}

///|
test "Bio.Align.maf parses three consecutive blocks" {
  assert_eq(align_maf_test_document().blocks.length(), 3)
}

///|
test "Bio.Align.maf parses a track line" {
  assert_true(align_maf_test_document().track is Some(_))
}

///|
test "Bio.Align.maf preserves a quoted track name" {
  let track = align_maf_test_document().track.unwrap()
  assert_eq(track.name, Some("Multiz demo"))
}

///|
test "Bio.Align.maf preserves a quoted track description" {
  let track = align_maf_test_document().track.unwrap()
  assert_eq(track.description, Some("MAF coordinate example"))
}

///|
test "Bio.Align.maf preserves mafDot metadata" {
  let track = align_maf_test_document().track.unwrap()
  assert_eq(track.maf_dot, Some("on"))
}

///|
test "Bio.Align.maf preserves visibility metadata" {
  let track = align_maf_test_document().track.unwrap()
  assert_eq(track.visibility, Some("pack"))
}

///|
test "Bio.Align.maf splits speciesOrder metadata" {
  let track = align_maf_test_document().track.unwrap()
  assert_eq(track.species_order, ["hg38", "mm39", "rn7"])
}

///|
test "Bio.Align.maf accepts a header-only document" {
  let document = align_maf_test_parse("##maf version=1\n")
  assert_eq(document.blocks.length(), 0)
}

///|
test "Bio.Align.maf accepts leading blank lines" {
  let document = align_maf_test_parse("\n\n##maf version=1\n")
  assert_eq(document.version, "1")
}

///|
test "Bio.Align.maf accepts CRLF input" {
  let text = "##maf version=1\r\n\r\na score=1\r\ns ref 0 2 + 10 AC\r\n"
  assert_eq(align_maf_test_parse(text).blocks.length(), 1)
}

///|
test "Bio.Align.maf accepts tab-separated component fields" {
  let text = align_maf_test_single("a\ns\tref\t0\t2\t+\t10\tAC\n")
  assert_eq(align_maf_test_parse(text).blocks[0].components.length(), 1)
}

///|
test "Bio.Align.maf parses a floating point score" {
  assert_eq(align_maf_test_document().blocks[0].score, Some(23262.0))
}

///|
test "Bio.Align.maf parses a positive pass annotation" {
  assert_eq(align_maf_test_document().blocks[0].pass_number, Some(1))
}

///|
test "Bio.Align.maf parses scientific score notation" {
  let text = align_maf_test_single("a score=1.25e2\ns ref 0 2 + 10 AC\n")
  assert_eq(align_maf_test_parse(text).blocks[0].score, Some(125.0))
}

///|
test "Bio.Align.maf parses an explicit positive score" {
  let text = align_maf_test_single("a score=+2.5\ns ref 0 2 + 10 AC\n")
  assert_eq(align_maf_test_parse(text).blocks[0].score, Some(2.5))
}

///|
test "Bio.Align.maf parses a negative score" {
  let text = align_maf_test_single("a score=-2.5\ns ref 0 2 + 10 AC\n")
  assert_eq(align_maf_test_parse(text).blocks[0].score, Some(-2.5))
}

///|
test "Bio.Align.maf allows an alignment without score" {
  let text = align_maf_test_single("a\ns ref 0 2 + 10 AC\n")
  assert_eq(align_maf_test_parse(text).blocks[0].score, None)
}

///|
test "Bio.Align.maf preserves aligned component fields" {
  let component = align_maf_test_document().blocks[0].components[0]
  assert_eq(component.source, "hg38.chr7")
  assert_eq(component.start, 100)
  assert_eq(component.size, 8)
  assert_eq(component.strand, "+")
  assert_eq(component.source_size, 1000)
  assert_eq(component.text, "AC-TGCAAT")
}

///|
test "Bio.Align.maf preserves a reverse component" {
  let component = align_maf_test_document().blocks[0].components[1]
  assert_eq(component.source, "mm39.chr5")
  assert_eq(component.strand, "-")
  assert_eq(component.forward_interval(), (442, 450))
}

///|
test "Bio.Align.maf attaches insertion context to the preceding component" {
  let insertion = align_maf_test_document().blocks[0].components[0].insertion
  assert_true(insertion is Some(_))
  assert_eq(insertion.unwrap().left_status, "C")
  assert_eq(insertion.unwrap().right_status, "I")
}

///|
test "Bio.Align.maf preserves insertion counts" {
  let insertion = align_maf_test_document().blocks[0].components[0].insertion.unwrap()
  assert_eq(insertion.left_count, 0)
  assert_eq(insertion.right_count, 2)
}

///|
test "Bio.Align.maf attaches aligned quality" {
  let quality = align_maf_test_document().blocks[0].components[0].quality
  assert_eq(quality, Some("98-765432"))
}

///|
test "Bio.Align.maf parses empty components" {
  let empty = align_maf_test_document().blocks[0].empty_components[0]
  assert_eq(empty.source, "rn7.chr1")
  assert_eq(empty.status, "I")
}

///|
test "Bio.Align.maf computes empty reverse intervals" {
  let text = align_maf_test_single("a\ns ref 0 2 + 10 AC\ne other 3 2 - 20 C\n")
  let empty = align_maf_test_parse(text).blocks[0].empty_components[0]
  assert_eq(empty.forward_interval(), (15, 17))
}

///|
test "Bio.Align.maf normalizes period gap characters" {
  let text = align_maf_test_single("a\ns ref 0 2 + 10 A.C\n")
  assert_eq(align_maf_test_parse(text).blocks[0].components[0].text, "A-C")
}

///|
test "Bio.Align.maf normalizes equals gap characters" {
  let text = align_maf_test_single("a\ns ref 0 2 + 10 A=C\n")
  assert_eq(align_maf_test_parse(text).blocks[0].components[0].text, "A-C")
}

///|
test "Bio.Align.maf normalizes underscore gap characters" {
  let text = align_maf_test_single("a\ns ref 0 2 + 10 A_C\n")
  assert_eq(align_maf_test_parse(text).blocks[0].components[0].text, "A-C")
}

///|
test "Bio.Align.maf ignores comments inside blocks" {
  let text = align_maf_test_single("a\n# block note\ns ref 0 2 + 10 AC\n")
  assert_eq(align_maf_test_parse(text).blocks.length(), 1)
}

///|
test "Bio.Align.maf preserves equal component widths" {
  assert_eq(align_maf_test_document().blocks[0].aligned_columns(), 9)
}

///|
test "Bio.Align.maf removes sequence gaps" {
  let component = align_maf_test_document().blocks[0].components[0]
  assert_eq(component.sequence(), "ACTGCAAT")
}

///|
test "Bio.Align.maf reports block sources in file order" {
  assert_eq(align_maf_test_document().blocks[0].sources(), [
    "hg38.chr7", "mm39.chr5",
  ])
}

///|
test "Bio.Align.maf reports document sources in first-seen order" {
  assert_eq(align_maf_test_document().sources(), [
    "hg38.chr7", "mm39.chr5", "rn7.chr1",
  ])
}

///|
test "Bio.Align.maf summary counts blocks" {
  assert_eq(align_maf_test_document().summary().block_count, 3)
}

///|
test "Bio.Align.maf summary counts aligned components" {
  assert_eq(align_maf_test_document().summary().component_count, 6)
}

///|
test "Bio.Align.maf summary counts empty components" {
  assert_eq(align_maf_test_document().summary().empty_component_count, 1)
}

///|
test "Bio.Align.maf summary counts alignment columns" {
  assert_eq(align_maf_test_document().summary().aligned_columns, 17)
}

///|
test "Bio.Align.maf summary counts reference bases" {
  assert_eq(align_maf_test_document().summary().reference_bases, 16)
}

///|
test "Bio.Align.maf summary counts distinct sources" {
  assert_eq(align_maf_test_document().summary().source_count, 3)
}

///|
test "Bio.Align.maf builds compressed coordinate breakpoints" {
  let path = @src.align_maf_coordinate_path(align_maf_test_document().blocks[0])
  assert_eq(path.length(), 6)
  assert_eq(path[0].column, 0)
  assert_eq(path[5].column, 9)
}

///|
test "Bio.Align.maf coordinate path uses absolute starts" {
  let path = @src.align_maf_coordinate_path(align_maf_test_document().blocks[0])
  assert_eq(path[0].positions, [100, 450])
}

///|
test "Bio.Align.maf coordinate path uses decreasing reverse coordinates" {
  let path = @src.align_maf_coordinate_path(align_maf_test_document().blocks[0])
  assert_eq(path[path.length() - 1].positions, [108, 442])
}

///|
test "Bio.Align.maf coordinate path marks first gap-pattern change" {
  let path = @src.align_maf_coordinate_path(align_maf_test_document().blocks[0])
  assert_eq(path[1].column, 2)
  assert_eq(path[1].positions, [102, 448])
}

///|
test "Bio.Align.maf coordinate path marks reference insertion boundary" {
  let path = @src.align_maf_coordinate_path(align_maf_test_document().blocks[0])
  assert_eq(path[2].column, 3)
  assert_eq(path[2].positions, [102, 447])
}

///|
test "Bio.Align.maf maps a reference base to reverse source coordinates" {
  let block = align_maf_test_document().blocks[0]
  assert_eq(
    @src.align_maf_map_position(block, "hg38.chr7", 100, "mm39.chr5"),
    Some(449),
  )
}

///|
test "Bio.Align.maf maps after a reference gap" {
  let block = align_maf_test_document().blocks[0]
  assert_eq(
    @src.align_maf_map_position(block, "hg38.chr7", 102, "mm39.chr5"),
    Some(446),
  )
}

///|
test "Bio.Align.maf maps reverse source coordinates back to reference" {
  let block = align_maf_test_document().blocks[0]
  assert_eq(
    @src.align_maf_map_position(block, "mm39.chr5", 449, "hg38.chr7"),
    Some(100),
  )
}

///|
test "Bio.Align.maf returns none when destination has a gap" {
  let block = align_maf_test_document().blocks[0]
  assert_eq(
    @src.align_maf_map_position(block, "mm39.chr5", 447, "hg38.chr7"),
    None,
  )
  assert_eq(
    @src.align_maf_map_position(block, "hg38.chr7", 104, "mm39.chr5"),
    None,
  )
}

///|
test "Bio.Align.maf returns none when source is absent" {
  let block = align_maf_test_document().blocks[0]
  assert_eq(@src.align_maf_map_position(block, "missing", 1, "hg38.chr7"), None)
}

///|
test "Bio.Align.maf returns none when destination is absent" {
  let block = align_maf_test_document().blocks[0]
  assert_eq(
    @src.align_maf_map_position(block, "hg38.chr7", 100, "missing"),
    None,
  )
}

///|
test "Bio.Align.maf returns none outside the source interval" {
  let block = align_maf_test_document().blocks[0]
  assert_eq(
    @src.align_maf_map_position(block, "hg38.chr7", 99, "mm39.chr5"),
    None,
  )
}

///|
test "Bio.Align.maf writes canonical header metadata" {
  let output = @src.align_maf_write(align_maf_test_document()) catch {
    AlignMafError(message) => abort("valid MAF write failed: " + message)
  }
  assert_true(output.contains("##maf version=1 scoring=blastz program=multiz"))
}

///|
test "Bio.Align.maf writes canonical track metadata" {
  let output = @src.align_maf_write(align_maf_test_document()) catch {
    AlignMafError(message) => abort("valid MAF write failed: " + message)
  }
  assert_true(output.has_prefix("track name=\"Multiz demo\""))
}

///|
test "Bio.Align.maf writes header comments" {
  let output = @src.align_maf_write(align_maf_test_document()) catch {
    AlignMafError(message) => abort("valid MAF write failed: " + message)
  }
  assert_true(output.contains("# generated fixture"))
}

///|
test "Bio.Align.maf writes insertion lines" {
  let output = @src.align_maf_write(align_maf_test_document()) catch {
    AlignMafError(message) => abort("valid MAF write failed: " + message)
  }
  assert_true(output.contains("i hg38.chr7 C 0 I 2"))
}

///|
test "Bio.Align.maf writes quality lines" {
  let output = @src.align_maf_write(align_maf_test_document()) catch {
    AlignMafError(message) => abort("valid MAF write failed: " + message)
  }
  assert_true(output.contains("q hg38.chr7 98-765432"))
}

///|
test "Bio.Align.maf writes empty component lines" {
  let output = @src.align_maf_write(align_maf_test_document()) catch {
    AlignMafError(message) => abort("valid MAF write failed: " + message)
  }
  assert_true(output.contains("e rn7.chr1 80 6 + 700 I"))
}

///|
test "Bio.Align.maf canonical output round-trips" {
  let document = align_maf_test_document()
  let output = @src.align_maf_write(document) catch {
    AlignMafError(message) => abort("valid MAF write failed: " + message)
  }
  assert_eq(align_maf_test_parse(output), document)
}

///|
test "Bio.Align.maf header-only output round-trips" {
  let document = align_maf_test_parse("##maf version=1\n")
  let output = @src.align_maf_write(document) catch {
    AlignMafError(message) => abort("valid MAF write failed: " + message)
  }
  assert_eq(align_maf_test_parse(output), document)
}

///|
test "Bio.Align.maf track writer escapes quotes" {
  let track = @src.AlignMafTrack::create(description=Some("quoted \"value\"")) catch {
    AlignMafError(message) => abort("valid track failed: " + message)
  }
  let document = @src.AlignMafDocument::create([], track=Some(track)) catch {
    AlignMafError(message) => abort("valid document failed: " + message)
  }
  let output = @src.align_maf_write(document) catch {
    AlignMafError(message) => abort("valid MAF write failed: " + message)
  }
  assert_eq(align_maf_test_parse(output), document)
}

///|
test "Bio.Align.maf track constructor copies speciesOrder" {
  let species = ["hg38"]
  let track = @src.AlignMafTrack::create(species_order=species) catch {
    AlignMafError(message) => abort("valid track failed: " + message)
  }
  species.push("mm39")
  assert_eq(track.species_order, ["hg38"])
}

///|
test "Bio.Align.maf document constructor copies blocks" {
  let blocks = align_maf_test_document().blocks.copy()
  let document = @src.AlignMafDocument::create(blocks) catch {
    AlignMafError(message) => abort("valid document failed: " + message)
  }
  ignore(blocks.pop())
  assert_eq(document.blocks.length(), 3)
}

///|
test "Bio.Align.maf index stores one entry per block" {
  assert_eq(align_maf_test_index().entry_count(), 3)
}

///|
test "Bio.Align.maf index preserves reference metadata" {
  let index = align_maf_test_index()
  assert_eq(index.reference, "hg38.chr7")
  assert_eq(index.source_size, 1000)
  assert_eq(index.reference_strand, "+")
}

///|
test "Bio.Align.maf index sorts reference intervals" {
  let index = align_maf_test_index()
  assert_eq(index.entries[0].start, 100)
  assert_eq(index.entries[1].start, 110)
  assert_eq(index.entries[2].start, 120)
}

///|
test "Bio.Align.maf index search finds an overlapping block" {
  let hits = align_maf_test_index().search(103, 104) catch {
    AlignMafError(message) => abort("valid search failed: " + message)
  }
  assert_eq(hits.length(), 1)
  assert_eq(hits[0].components[0].start, 100)
}

///|
test "Bio.Align.maf index search uses half-open left boundary" {
  let hits = align_maf_test_index().search(108, 110) catch {
    AlignMafError(message) => abort("valid search failed: " + message)
  }
  assert_eq(hits.length(), 0)
}

///|
test "Bio.Align.maf index search uses half-open right boundary" {
  let hits = align_maf_test_index().search(107, 108) catch {
    AlignMafError(message) => abort("valid search failed: " + message)
  }
  assert_eq(hits.length(), 1)
}

///|
test "Bio.Align.maf index search spans multiple blocks" {
  let hits = align_maf_test_index().search(107, 121) catch {
    AlignMafError(message) => abort("valid search failed: " + message)
  }
  assert_eq(hits.length(), 3)
}

///|
test "Bio.Align.maf index range search deduplicates blocks" {
  let hits = align_maf_test_index().search_ranges([100, 102], [104, 106]) catch {
    AlignMafError(message) => abort("valid range search failed: " + message)
  }
  assert_eq(hits.length(), 1)
}

///|
test "Bio.Align.maf index range search returns file order" {
  let hits = align_maf_test_index().search_ranges([120, 100], [121, 101]) catch {
    AlignMafError(message) => abort("valid range search failed: " + message)
  }
  assert_eq(hits.length(), 2)
  assert_eq(hits[0].components[0].start, 100)
  assert_eq(hits[1].components[0].start, 120)
}

///|
test "Bio.Align.maf splices two reference exons" {
  let spliced = align_maf_test_index().get_spliced([102, 110], [106, 113]) catch {
    AlignMafError(message) => abort("valid splice failed: " + message)
  }
  assert_eq(spliced.reference, "hg38.chr7")
  assert_eq(spliced.columns, 7)
  assert_eq(spliced.sequence("hg38.chr7"), Some("TGCAGGT"))
}

///|
test "Bio.Align.maf splicing preserves target gaps" {
  let spliced = align_maf_test_index().get_spliced([102, 110], [106, 113]) catch {
    AlignMafError(message) => abort("valid splice failed: " + message)
  }
  assert_eq(spliced.sequence("mm39.chr5"), Some("TG-AG-T"))
}

///|
test "Bio.Align.maf splicing fills missing species with gaps" {
  let spliced = align_maf_test_index().get_spliced([120], [124]) catch {
    AlignMafError(message) => abort("valid splice failed: " + message)
  }
  assert_eq(spliced.sequence("mm39.chr5"), Some("----"))
  assert_eq(spliced.sequence("rn7.chr1"), Some("AATA"))
}

///|
test "Bio.Align.maf splicing fills uncovered reference with N" {
  let spliced = align_maf_test_index().get_spliced([106], [112]) catch {
    AlignMafError(message) => abort("valid splice failed: " + message)
  }
  assert_eq(spliced.sequence("hg38.chr7"), Some("ATNNGG"))
}

///|
test "Bio.Align.maf splicing fills uncovered targets with gaps" {
  let spliced = align_maf_test_index().get_spliced([106], [112]) catch {
    AlignMafError(message) => abort("valid splice failed: " + message)
  }
  assert_eq(spliced.sequence("mm39.chr5"), Some("AT--G-"))
}

///|
test "Bio.Align.maf reverse-strand splicing reverse-complements output" {
  let spliced = align_maf_test_index().get_spliced(
    [102, 110],
    [106, 113],
    strand="-",
  ) catch {
    AlignMafError(message) => abort("valid splice failed: " + message)
  }
  assert_eq(spliced.sequence("hg38.chr7"), Some("ACCTGCA"))
  assert_eq(spliced.strand, "-")
}

///|
test "Bio.Align.maf spliced output records exon coordinates" {
  let spliced = align_maf_test_index().get_spliced([102, 110], [106, 113]) catch {
    AlignMafError(message) => abort("valid splice failed: " + message)
  }
  assert_eq(spliced.starts, [102, 110])
  assert_eq(spliced.ends, [106, 113])
}

///|
test "Bio.Align.maf spliced sequence lookup returns none for missing source" {
  let spliced = align_maf_test_index().get_spliced([102], [106]) catch {
    AlignMafError(message) => abort("valid splice failed: " + message)
  }
  assert_eq(spliced.sequence("missing"), None)
}

///|
test "Bio.Align.maf reports spliced ungapped length" {
  let spliced = align_maf_test_index().get_spliced([102, 110], [106, 113]) catch {
    AlignMafError(message) => abort("valid splice failed: " + message)
  }
  assert_eq(spliced.ungapped_length("mm39.chr5"), Some(5))
}

///|
test "Bio.Align.maf excludes uncovered N from ungapped length" {
  let spliced = align_maf_test_index().get_spliced([106], [112]) catch {
    AlignMafError(message) => abort("valid splice failed: " + message)
  }
  assert_eq(spliced.ungapped_length("hg38.chr7"), Some(4))
}

///|
test "Bio.Align.maf rejects empty input" {
  assert_true(align_maf_test_rejects(""))
}

///|
test "Bio.Align.maf rejects input without a header" {
  assert_true(align_maf_test_rejects("a\ns ref 0 2 + 10 AC\n"))
}

///|
test "Bio.Align.maf rejects a header without version" {
  assert_true(align_maf_test_rejects("##maf scoring=x\n"))
}

///|
test "Bio.Align.maf rejects unsupported versions" {
  assert_true(align_maf_test_rejects("##maf version=2\n"))
}

///|
test "Bio.Align.maf rejects duplicate version metadata" {
  assert_true(align_maf_test_rejects("##maf version=1 version=1\n"))
}

///|
test "Bio.Align.maf rejects duplicate scoring metadata" {
  assert_true(align_maf_test_rejects("##maf version=1 scoring=x scoring=y\n"))
}

///|
test "Bio.Align.maf rejects duplicate program metadata" {
  assert_true(align_maf_test_rejects("##maf version=1 program=x program=y\n"))
}

///|
test "Bio.Align.maf rejects unknown header metadata" {
  assert_true(align_maf_test_rejects("##maf version=1 unknown=x\n"))
}

///|
test "Bio.Align.maf rejects malformed header metadata" {
  assert_true(align_maf_test_rejects("##maf version\n"))
}

///|
test "Bio.Align.maf rejects a data line before an alignment line" {
  assert_true(align_maf_test_rejects("##maf version=1\ns ref 0 2 + 10 AC\n"))
}

///|
test "Bio.Align.maf rejects an empty track line" {
  assert_true(align_maf_test_rejects("track\n##maf version=1\n"))
}

///|
test "Bio.Align.maf rejects unterminated track quotes" {
  assert_true(align_maf_test_rejects("track name=\"broken\n##maf version=1\n"))
}

///|
test "Bio.Align.maf rejects unknown track variables" {
  assert_true(align_maf_test_rejects("track unknown=x\n##maf version=1\n"))
}

///|
test "Bio.Align.maf rejects duplicate track names" {
  assert_true(align_maf_test_rejects("track name=x name=y\n##maf version=1\n"))
}

///|
test "Bio.Align.maf rejects invalid mafDot values" {
  assert_true(align_maf_test_rejects("track mafDot=yes\n##maf version=1\n"))
}

///|
test "Bio.Align.maf rejects invalid visibility values" {
  assert_true(
    align_maf_test_rejects("track visibility=hide\n##maf version=1\n"),
  )
}

///|
test "Bio.Align.maf rejects empty speciesOrder" {
  assert_true(
    align_maf_test_rejects("track speciesOrder=\"\"\n##maf version=1\n"),
  )
}

///|
test "Bio.Align.maf rejects duplicate speciesOrder" {
  assert_true(
    align_maf_test_rejects(
      "track speciesOrder=x speciesOrder=y\n##maf version=1\n",
    ),
  )
}

///|
test "Bio.Align.maf rejects malformed alignment annotations" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a score\ns ref 0 2 + 10 AC\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects unknown alignment annotations" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a label=x\ns ref 0 2 + 10 AC\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects duplicate score annotations" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a score=1 score=2\ns ref 0 2 + 10 AC\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects duplicate pass annotations" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a pass=1 pass=2\ns ref 0 2 + 10 AC\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects zero pass annotations" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a pass=0\ns ref 0 2 + 10 AC\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects nonnumeric scores" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a score=abc\ns ref 0 2 + 10 AC\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects non-finite score magnitudes" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a score=1.0e301\ns ref 0 2 + 10 AC\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects malformed sequence rows" {
  assert_true(
    align_maf_test_rejects(align_maf_test_single("a\ns ref 0 2 + 10\n")),
  )
}

///|
test "Bio.Align.maf rejects negative starts" {
  assert_true(
    align_maf_test_rejects(align_maf_test_single("a\ns ref -1 2 + 10 AC\n")),
  )
}

///|
test "Bio.Align.maf rejects zero component sizes" {
  assert_true(
    align_maf_test_rejects(align_maf_test_single("a\ns ref 0 0 + 10 A\n")),
  )
}

///|
test "Bio.Align.maf rejects zero source sizes" {
  assert_true(
    align_maf_test_rejects(align_maf_test_single("a\ns ref 0 1 + 0 A\n")),
  )
}

///|
test "Bio.Align.maf rejects coordinates beyond source size" {
  assert_true(
    align_maf_test_rejects(align_maf_test_single("a\ns ref 9 2 + 10 AC\n")),
  )
}

///|
test "Bio.Align.maf rejects integer overflow" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\ns ref 2147483648 1 + 2147483647 A\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects additive coordinate overflow" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\ns ref 2147483647 1 + 2147483647 A\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects invalid component strands" {
  assert_true(
    align_maf_test_rejects(align_maf_test_single("a\ns ref 0 2 ? 10 AC\n")),
  )
}

///|
test "Bio.Align.maf rejects size and sequence disagreement" {
  assert_true(
    align_maf_test_rejects(align_maf_test_single("a\ns ref 0 3 + 10 AC\n")),
  )
}

///|
test "Bio.Align.maf rejects unequal component widths" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\ns ref 0 2 + 10 AC\ns other 0 3 + 10 ACG\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects blocks without aligned components" {
  assert_true(
    align_maf_test_rejects(align_maf_test_single("a\ne other 0 2 + 10 C\n")),
  )
}

///|
test "Bio.Align.maf rejects insertion rows before sequence rows" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\ni ref C 0 C 0\ns ref 0 2 + 10 AC\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects malformed insertion rows" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\ns ref 0 2 + 10 AC\ni ref C 0 C\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects insertion source mismatches" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\ns ref 0 2 + 10 AC\ni other C 0 C 0\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects unknown insertion status" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\ns ref 0 2 + 10 AC\ni ref X 0 C 0\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects duplicate insertion rows" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single(
        "a\ns ref 0 2 + 10 AC\ni ref C 0 C 0\ni ref C 0 C 0\n",
      ),
    ),
  )
}

///|
test "Bio.Align.maf rejects quality rows before sequence rows" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\nq ref 98\ns ref 0 2 + 10 AC\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects malformed quality rows" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\ns ref 0 2 + 10 AC\nq ref\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects quality source mismatches" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\ns ref 0 2 + 10 AC\nq other 98\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects quality width mismatches" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\ns ref 0 2 + 10 AC\nq ref 9\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects quality gap mismatches" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\ns ref 0 2 + 10 A-C\nq ref 999\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects duplicate quality rows" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\ns ref 0 2 + 10 AC\nq ref 98\nq ref 98\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects malformed empty rows" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\ns ref 0 2 + 10 AC\ne other 0 2 + 10\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects unknown empty status" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\ns ref 0 2 + 10 AC\ne other 0 2 + 10 X\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects empty coordinates beyond source size" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\ns ref 0 2 + 10 AC\ne other 9 2 + 10 C\n"),
    ),
  )
}

///|
test "Bio.Align.maf rejects unexpected line types" {
  assert_true(
    align_maf_test_rejects(
      align_maf_test_single("a\ns ref 0 2 + 10 AC\nx invalid\n"),
    ),
  )
}

///|
test "Bio.Align.maf track constructor rejects invalid mafDot" {
  let rejected = try {
    ignore(@src.AlignMafTrack::create(maf_dot=Some("maybe")))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf track constructor rejects invalid visibility" {
  let rejected = try {
    ignore(@src.AlignMafTrack::create(visibility=Some("hidden")))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf document constructor rejects unsupported version" {
  let rejected = try {
    ignore(@src.AlignMafDocument::create([], version="2"))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf index rejects empty documents" {
  let document = @src.AlignMafDocument::create([]) catch {
    AlignMafError(message) =>
      abort("empty document should be valid: " + message)
  }
  let rejected = try {
    ignore(@src.AlignMafIndex::create(document, "ref"))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf index rejects missing references" {
  let rejected = try {
    ignore(@src.AlignMafIndex::create(align_maf_test_document(), "missing"))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf index rejects inconsistent source sizes" {
  let text = align_maf_test_single(
    "a\ns ref 0 2 + 10 AC\n\na\ns ref 2 2 + 11 GT\n",
  )
  let rejected = try {
    ignore(@src.AlignMafIndex::create(align_maf_test_parse(text), "ref"))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf index rejects inconsistent reference strands" {
  let text = align_maf_test_single(
    "a\ns ref 0 2 + 10 AC\n\na\ns ref 2 2 - 10 GT\n",
  )
  let rejected = try {
    ignore(@src.AlignMafIndex::create(align_maf_test_parse(text), "ref"))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf index rejects duplicate references in one block" {
  let text = align_maf_test_single("a\ns ref 0 2 + 10 AC\ns ref 2 2 + 10 GT\n")
  let rejected = try {
    ignore(@src.AlignMafIndex::create(align_maf_test_parse(text), "ref"))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf search rejects empty ranges" {
  let rejected = try {
    ignore(align_maf_test_index().search(100, 100))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf search rejects negative ranges" {
  let rejected = try {
    ignore(align_maf_test_index().search(-1, 1))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf search rejects ranges beyond source size" {
  let rejected = try {
    ignore(align_maf_test_index().search(999, 1001))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf range search rejects empty arrays" {
  let rejected = try {
    ignore(align_maf_test_index().search_ranges([], []))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf range search rejects mismatched arrays" {
  let rejected = try {
    ignore(align_maf_test_index().search_ranges([1], [2, 3]))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf splicing rejects empty exon arrays" {
  let rejected = try {
    ignore(align_maf_test_index().get_spliced([], []))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf splicing rejects mismatched exon arrays" {
  let rejected = try {
    ignore(align_maf_test_index().get_spliced([100], [102, 104]))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf splicing rejects invalid strands" {
  let rejected = try {
    ignore(align_maf_test_index().get_spliced([100], [102], strand="?"))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf splicing rejects overlapping exons" {
  let rejected = try {
    ignore(align_maf_test_index().get_spliced([100, 103], [105, 106]))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf splicing rejects unsorted exons" {
  let rejected = try {
    ignore(align_maf_test_index().get_spliced([110, 100], [112, 102]))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf splicing rejects overlapping indexed blocks" {
  let text = align_maf_test_single(
    "a\ns ref 0 4 + 20 ACGT\n\na\ns ref 2 4 + 20 GTAA\n",
  )
  let index = @src.AlignMafIndex::create(align_maf_test_parse(text), "ref") catch {
    AlignMafError(message) => abort("overlap index should build: " + message)
  }
  let rejected = try {
    ignore(index.get_spliced([0], [6]))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.maf splicing rejects reverse reference indexes" {
  let text = align_maf_test_single("a\ns ref 0 4 - 20 ACGT\n")
  let index = @src.AlignMafIndex::create(align_maf_test_parse(text), "ref") catch {
    AlignMafError(message) => abort("reverse index should build: " + message)
  }
  let rejected = try {
    ignore(index.get_spliced([16], [20]))
    false
  } catch {
    AlignMafError(_) => true
  }
  assert_true(rejected)
}
