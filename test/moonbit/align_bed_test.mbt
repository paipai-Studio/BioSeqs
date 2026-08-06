///|
fn align_bed_test_parse(text : String) -> @src.AlignBedDocument {
  @src.align_bed_parse(text) catch {
    AlignBedError(message) => abort("valid BED input failed: " + message)
  }
}

///|
fn align_bed_test_document() -> @src.AlignBedDocument {
  align_bed_test_parse(@src.align_bed_example_text())
}

///|
fn align_bed_test_plus() -> @src.AlignBedAlignment {
  align_bed_test_document().alignments[0]
}

///|
fn align_bed_test_minus() -> @src.AlignBedAlignment {
  align_bed_test_document().alignments[1]
}

///|
fn align_bed_test_simple() -> @src.AlignBedAlignment {
  align_bed_test_document().alignments[2]
}

///|
fn align_bed_test_rejects(text : String) -> Bool {
  try {
    ignore(@src.align_bed_parse(text))
    false
  } catch {
    AlignBedError(_) => true
  }
}

///|
fn align_bed_test_coordinate(
  target : Int,
  query : Int,
) -> @src.AlignBedCoordinate {
  @src.AlignBedCoordinate::create(target, query) catch {
    AlignBedError(message) => abort("valid BED coordinate failed: " + message)
  }
}

///|
fn align_bed_test_construct(
  coordinates : Array[@src.AlignBedCoordinate],
) -> @src.AlignBedAlignment {
  @src.AlignBedAlignment::create("target", Some("query"), coordinates) catch {
    AlignBedError(message) => abort("valid BED alignment failed: " + message)
  }
}

///|
fn align_bed_test_constructor_rejects(
  coordinates : Array[@src.AlignBedCoordinate],
) -> Bool {
  try {
    ignore(@src.AlignBedAlignment::create("target", Some("query"), coordinates))
    false
  } catch {
    AlignBedError(_) => true
  }
}

///|
test "Bio.Align.bed parses three records" {
  assert_eq(align_bed_test_document().alignments.length(), 3)
}

///|
test "Bio.Align.bed parses an empty document" {
  assert_eq(align_bed_test_parse("").alignments.length(), 0)
}

///|
test "Bio.Align.bed parses a final newline" {
  assert_eq(align_bed_test_parse("chr1\t0\t10\n").alignments.length(), 1)
}

///|
test "Bio.Align.bed parses input without a final newline" {
  assert_eq(align_bed_test_parse("chr1\t0\t10").alignments.length(), 1)
}

///|
test "Bio.Align.bed accepts CRLF records" {
  assert_eq(
    align_bed_test_parse("chr1\t0\t10\r\nchr2\t5\t15\r\n").alignments.length(),
    2,
  )
}

///|
test "Bio.Align.bed accepts whitespace-delimited fields" {
  let alignment = align_bed_test_parse("chr1  10  20  tx1\n").alignments[0]
  assert_eq(alignment.query_id, Some("tx1"))
}

///|
test "Bio.Align.bed preserves target identifier" {
  assert_eq(align_bed_test_plus().target_id, "chr22")
}

///|
test "Bio.Align.bed preserves plus query identifier" {
  assert_eq(align_bed_test_plus().query_id, Some("mRNA1"))
}

///|
test "Bio.Align.bed preserves minus query identifier" {
  assert_eq(align_bed_test_minus().query_id, Some("mRNA2"))
}

///|
test "Bio.Align.bed BED3 query identifier is absent" {
  assert_eq(align_bed_test_simple().query_id, None)
}

///|
test "Bio.Align.bed records source column count" {
  assert_eq(align_bed_test_plus().source_columns, 12)
  assert_eq(align_bed_test_simple().source_columns, 3)
}

///|
test "Bio.Align.bed parses numeric score" {
  assert_eq(align_bed_test_plus().score.unwrap().number(), Some(960.0))
}

///|
test "Bio.Align.bed parses second numeric score" {
  assert_eq(align_bed_test_minus().score.unwrap().number(), Some(900.0))
}

///|
test "Bio.Align.bed BED3 score is absent" {
  assert_eq(align_bed_test_simple().score, None)
}

///|
test "Bio.Align.bed parses floating score" {
  let alignment = align_bed_test_parse("chr1\t0\t10\ttx\t2.5\n").alignments[0]
  assert_eq(alignment.score.unwrap().number(), Some(2.5))
}

///|
test "Bio.Align.bed parses scientific score" {
  let alignment = align_bed_test_parse("chr1\t0\t10\ttx\t1.25e2\n").alignments[0]
  assert_eq(alignment.score.unwrap().number(), Some(125.0))
}

///|
test "Bio.Align.bed parses explicit positive score" {
  let alignment = align_bed_test_parse("chr1\t0\t10\ttx\t+2.5\n").alignments[0]
  assert_eq(alignment.score.unwrap().number(), Some(2.5))
}

///|
test "Bio.Align.bed preserves textual score" {
  let alignment = align_bed_test_parse("chr1\t0\t10\ttx\tpass\n").alignments[0]
  assert_eq(alignment.score.unwrap().text_value(), Some("pass"))
}

///|
test "Bio.Align.bed parses thickStart" {
  assert_eq(align_bed_test_plus().thick_start, Some(1200))
}

///|
test "Bio.Align.bed parses thickEnd" {
  assert_eq(align_bed_test_plus().thick_end, Some(4900))
}

///|
test "Bio.Align.bed parses itemRgb" {
  assert_eq(align_bed_test_plus().item_rgb, Some("255,0,0"))
}

///|
test "Bio.Align.bed identifies plus strand" {
  assert_eq(align_bed_test_plus().strand(), "+")
}

///|
test "Bio.Align.bed identifies minus strand" {
  assert_eq(align_bed_test_minus().strand(), "-")
}

///|
test "Bio.Align.bed normalizes dot strand to plus" {
  let alignment = align_bed_test_parse("chr1\t0\t10\ttx\t0\t.\n").alignments[0]
  assert_eq(alignment.strand(), "+")
}

///|
test "Bio.Align.bed plus coordinate count" {
  assert_eq(align_bed_test_plus().coordinates.length(), 4)
}

///|
test "Bio.Align.bed plus first coordinate" {
  assert_eq(
    align_bed_test_plus().coordinates[0],
    align_bed_test_coordinate(1000, 0),
  )
}

///|
test "Bio.Align.bed plus first block end" {
  assert_eq(
    align_bed_test_plus().coordinates[1],
    align_bed_test_coordinate(1567, 567),
  )
}

///|
test "Bio.Align.bed plus second block start" {
  assert_eq(
    align_bed_test_plus().coordinates[2],
    align_bed_test_coordinate(4512, 567),
  )
}

///|
test "Bio.Align.bed plus final coordinate" {
  assert_eq(
    align_bed_test_plus().coordinates[3],
    align_bed_test_coordinate(5000, 1055),
  )
}

///|
test "Bio.Align.bed minus coordinate count" {
  assert_eq(align_bed_test_minus().coordinates.length(), 4)
}

///|
test "Bio.Align.bed minus first coordinate" {
  assert_eq(
    align_bed_test_minus().coordinates[0],
    align_bed_test_coordinate(2000, 832),
  )
}

///|
test "Bio.Align.bed minus first block end" {
  assert_eq(
    align_bed_test_minus().coordinates[1],
    align_bed_test_coordinate(2433, 399),
  )
}

///|
test "Bio.Align.bed minus second block start" {
  assert_eq(
    align_bed_test_minus().coordinates[2],
    align_bed_test_coordinate(5601, 399),
  )
}

///|
test "Bio.Align.bed minus final coordinate" {
  assert_eq(
    align_bed_test_minus().coordinates[3],
    align_bed_test_coordinate(6000, 0),
  )
}

///|
test "Bio.Align.bed BED3 has one continuous coordinate step" {
  assert_eq(align_bed_test_simple().coordinates.length(), 2)
  assert_eq(
    align_bed_test_simple().coordinates[0],
    align_bed_test_coordinate(100, 0),
  )
  assert_eq(
    align_bed_test_simple().coordinates[1],
    align_bed_test_coordinate(180, 80),
  )
}

///|
test "Bio.Align.bed plus block count" {
  assert_eq(align_bed_test_plus().blocks().length(), 2)
}

///|
test "Bio.Align.bed plus first block geometry" {
  let block = align_bed_test_plus().blocks()[0]
  assert_eq(block.target_start, 1000)
  assert_eq(block.target_end, 1567)
  assert_eq(block.query_start, 0)
  assert_eq(block.query_end, 567)
}

///|
test "Bio.Align.bed plus second block geometry" {
  let block = align_bed_test_plus().blocks()[1]
  assert_eq(block.target_start, 4512)
  assert_eq(block.target_end, 5000)
  assert_eq(block.query_start, 567)
  assert_eq(block.query_end, 1055)
}

///|
test "Bio.Align.bed minus first block query interval" {
  let block = align_bed_test_minus().blocks()[0]
  assert_eq(block.query_start, 399)
  assert_eq(block.query_end, 832)
}

///|
test "Bio.Align.bed target interval" {
  assert_eq(align_bed_test_plus().target_start(), 1000)
  assert_eq(align_bed_test_plus().target_end(), 5000)
}

///|
test "Bio.Align.bed plus query size" {
  assert_eq(align_bed_test_plus().query_size(), 1055)
}

///|
test "Bio.Align.bed minus query size" {
  assert_eq(align_bed_test_minus().query_size(), 832)
}

///|
test "Bio.Align.bed target span" {
  assert_eq(align_bed_test_plus().target_span(), 4000)
}

///|
test "Bio.Align.bed plus aligned bases" {
  assert_eq(align_bed_test_plus().aligned_bases(), 1055)
}

///|
test "Bio.Align.bed plus target skip count" {
  let counts = align_bed_test_plus().counts()
  assert_eq(counts.target_skip_bases, 2945)
  assert_eq(counts.target_skip_opens, 1)
}

///|
test "Bio.Align.bed minus target skip count" {
  let counts = align_bed_test_minus().counts()
  assert_eq(counts.target_skip_bases, 3168)
  assert_eq(counts.target_skip_opens, 1)
}

///|
test "Bio.Align.bed parsed records have no query skips" {
  let counts = align_bed_test_plus().counts()
  assert_eq(counts.query_skip_bases, 0)
  assert_eq(counts.query_skip_opens, 0)
}

///|
test "Bio.Align.bed alignment columns include target skips" {
  assert_eq(align_bed_test_plus().counts().columns, 4000)
}

///|
test "Bio.Align.bed maps plus target start" {
  assert_eq(align_bed_test_plus().map_target_position(1000), Some(0))
}

///|
test "Bio.Align.bed maps plus first block end residue" {
  assert_eq(align_bed_test_plus().map_target_position(1566), Some(566))
}

///|
test "Bio.Align.bed returns none inside target intron" {
  assert_eq(align_bed_test_plus().map_target_position(2000), None)
}

///|
test "Bio.Align.bed maps plus second block start" {
  assert_eq(align_bed_test_plus().map_target_position(4512), Some(567))
}

///|
test "Bio.Align.bed maps plus final target residue" {
  assert_eq(align_bed_test_plus().map_target_position(4999), Some(1054))
}

///|
test "Bio.Align.bed maps plus query start" {
  assert_eq(align_bed_test_plus().map_query_position(0), Some(1000))
}

///|
test "Bio.Align.bed maps plus query final residue" {
  assert_eq(align_bed_test_plus().map_query_position(1054), Some(4999))
}

///|
test "Bio.Align.bed maps minus first target residue" {
  assert_eq(align_bed_test_minus().map_target_position(2000), Some(831))
}

///|
test "Bio.Align.bed maps minus first exon final residue" {
  assert_eq(align_bed_test_minus().map_target_position(2432), Some(399))
}

///|
test "Bio.Align.bed maps minus second exon start" {
  assert_eq(align_bed_test_minus().map_target_position(5601), Some(398))
}

///|
test "Bio.Align.bed maps minus final target residue" {
  assert_eq(align_bed_test_minus().map_target_position(5999), Some(0))
}

///|
test "Bio.Align.bed maps reverse query zero" {
  assert_eq(align_bed_test_minus().map_query_position(0), Some(5999))
}

///|
test "Bio.Align.bed maps reverse query final residue" {
  assert_eq(align_bed_test_minus().map_query_position(831), Some(2000))
}

///|
test "Bio.Align.bed mapping rejects negative positions" {
  assert_eq(align_bed_test_plus().map_target_position(-1), None)
  assert_eq(align_bed_test_plus().map_query_position(-1), None)
}

///|
test "Bio.Align.bed overlap uses half-open intervals" {
  let alignment = align_bed_test_plus()
  assert_true(alignment.overlaps("chr22", 999, 1001))
  assert_false(alignment.overlaps("chr22", 5000, 5001))
}

///|
test "Bio.Align.bed overlap requires matching target" {
  assert_false(align_bed_test_plus().overlaps("chr7", 1000, 1001))
}

///|
test "Bio.Align.bed searches target intervals" {
  let hits = align_bed_test_document().search("chr22", 1500, 2500) catch {
    AlignBedError(message) => abort("valid BED search failed: " + message)
  }
  assert_eq(hits.length(), 2)
}

///|
test "Bio.Align.bed search respects half-open end" {
  let hits = align_bed_test_document().search("chr22", 6000, 6100) catch {
    AlignBedError(message) => abort("valid BED search failed: " + message)
  }
  assert_eq(hits.length(), 0)
}

///|
test "Bio.Align.bed search rejects empty interval" {
  let document = align_bed_test_document()
  assert_true(
    try {
      ignore(document.search("chr22", 10, 10))
      false
    } catch {
      AlignBedError(_) => true
    },
  )
}

///|
test "Bio.Align.bed preserves distinct target order" {
  assert_eq(align_bed_test_document().targets(), ["chr22", "chr7"])
}

///|
test "Bio.Align.bed preserves distinct query order" {
  assert_eq(align_bed_test_document().queries(), ["mRNA1", "mRNA2"])
}

///|
test "Bio.Align.bed summary alignment count" {
  assert_eq(align_bed_test_document().summary().alignment_count, 3)
}

///|
test "Bio.Align.bed summary target and query counts" {
  let summary = align_bed_test_document().summary()
  assert_eq(summary.target_count, 2)
  assert_eq(summary.query_count, 2)
}

///|
test "Bio.Align.bed summary strand counts" {
  let summary = align_bed_test_document().summary()
  assert_eq(summary.plus_strand_count, 2)
  assert_eq(summary.minus_strand_count, 1)
}

///|
test "Bio.Align.bed summary aligned bases" {
  assert_eq(align_bed_test_document().summary().aligned_bases, 1967)
}

///|
test "Bio.Align.bed summary target skips" {
  assert_eq(align_bed_test_document().summary().target_skip_bases, 6113)
}

///|
test "Bio.Align.bed writes canonical BED12 plus record" {
  let line = @src.align_bed_format(align_bed_test_plus()) catch {
    AlignBedError(message) => abort("valid BED writing failed: " + message)
  }
  assert_eq(
    line, "chr22\t1000\t5000\tmRNA1\t960\t+\t1200\t4900\t255,0,0\t2\t567,488,\t0,3512,\n",
  )
}

///|
test "Bio.Align.bed writes canonical BED12 minus record" {
  let line = @src.align_bed_format(align_bed_test_minus()) catch {
    AlignBedError(message) => abort("valid BED writing failed: " + message)
  }
  assert_eq(
    line, "chr22\t2000\t6000\tmRNA2\t900\t-\t2300\t5960\t0,255,0\t2\t433,399,\t0,3601,\n",
  )
}

///|
test "Bio.Align.bed writes BED3" {
  let line = @src.align_bed_format(align_bed_test_plus(), bed_columns=3) catch {
    AlignBedError(message) => abort("valid BED3 writing failed: " + message)
  }
  assert_eq(line, "chr22\t1000\t5000\n")
}

///|
test "Bio.Align.bed writes BED4" {
  let line = @src.align_bed_format(align_bed_test_plus(), bed_columns=4) catch {
    AlignBedError(message) => abort("valid BED4 writing failed: " + message)
  }
  assert_eq(line, "chr22\t1000\t5000\tmRNA1\n")
}

///|
test "Bio.Align.bed writes BED5" {
  let line = @src.align_bed_format(align_bed_test_plus(), bed_columns=5) catch {
    AlignBedError(message) => abort("valid BED5 writing failed: " + message)
  }
  assert_eq(line, "chr22\t1000\t5000\tmRNA1\t960\n")
}

///|
test "Bio.Align.bed writes BED6" {
  let line = @src.align_bed_format(align_bed_test_plus(), bed_columns=6) catch {
    AlignBedError(message) => abort("valid BED6 writing failed: " + message)
  }
  assert_eq(line, "chr22\t1000\t5000\tmRNA1\t960\t+\n")
}

///|
test "Bio.Align.bed writes BED7" {
  let line = @src.align_bed_format(align_bed_test_plus(), bed_columns=7) catch {
    AlignBedError(message) => abort("valid BED7 writing failed: " + message)
  }
  assert_eq(line, "chr22\t1000\t5000\tmRNA1\t960\t+\t1200\n")
}

///|
test "Bio.Align.bed writes BED8" {
  let line = @src.align_bed_format(align_bed_test_plus(), bed_columns=8) catch {
    AlignBedError(message) => abort("valid BED8 writing failed: " + message)
  }
  assert_eq(line, "chr22\t1000\t5000\tmRNA1\t960\t+\t1200\t4900\n")
}

///|
test "Bio.Align.bed writes BED9" {
  let line = @src.align_bed_format(align_bed_test_plus(), bed_columns=9) catch {
    AlignBedError(message) => abort("valid BED9 writing failed: " + message)
  }
  assert_eq(line, "chr22\t1000\t5000\tmRNA1\t960\t+\t1200\t4900\t255,0,0\n")
}

///|
test "Bio.Align.bed writes BED10" {
  let line = @src.align_bed_format(align_bed_test_plus(), bed_columns=10) catch {
    AlignBedError(message) => abort("valid BED10 writing failed: " + message)
  }
  assert_eq(line, "chr22\t1000\t5000\tmRNA1\t960\t+\t1200\t4900\t255,0,0\t2\n")
}

///|
test "Bio.Align.bed writes BED11" {
  let line = @src.align_bed_format(align_bed_test_plus(), bed_columns=11) catch {
    AlignBedError(message) => abort("valid BED11 writing failed: " + message)
  }
  assert_eq(
    line, "chr22\t1000\t5000\tmRNA1\t960\t+\t1200\t4900\t255,0,0\t2\t567,488,\n",
  )
}

///|
test "Bio.Align.bed BED12 round trip is stable" {
  let line = @src.align_bed_format(align_bed_test_plus()) catch {
    AlignBedError(message) => abort("valid BED writing failed: " + message)
  }
  let parsed = align_bed_test_parse(line).alignments[0]
  let second = @src.align_bed_format(parsed) catch {
    AlignBedError(message) => abort("valid BED writing failed: " + message)
  }
  assert_eq(second, line)
}

///|
test "Bio.Align.bed writes complete document" {
  let document = align_bed_test_document()
  let output = @src.align_bed_write(document, bed_columns=3) catch {
    AlignBedError(message) =>
      abort("valid BED document writing failed: " + message)
  }
  assert_eq(output.split("\n").length(), 4)
}

///|
test "Bio.Align.bed textual score writes unchanged" {
  let alignment = align_bed_test_parse("chr1\t0\t10\ttx\tpass\n").alignments[0]
  let output = @src.align_bed_format(alignment, bed_columns=5) catch {
    AlignBedError(message) =>
      abort("valid BED score writing failed: " + message)
  }
  assert_eq(output, "chr1\t0\t10\ttx\tpass\n")
}

///|
test "Bio.Align.bed writer supplies default query name" {
  let output = @src.align_bed_format(align_bed_test_simple(), bed_columns=4) catch {
    AlignBedError(message) => abort("valid BED writing failed: " + message)
  }
  assert_eq(output, "chr7\t100\t180\tquery\n")
}

///|
test "Bio.Align.bed writer supplies default score" {
  let output = @src.align_bed_format(align_bed_test_simple(), bed_columns=5) catch {
    AlignBedError(message) => abort("valid BED writing failed: " + message)
  }
  assert_eq(output, "chr7\t100\t180\tquery\t0\n")
}

///|
test "Bio.Align.bed writer supplies default thick interval and color" {
  let output = @src.align_bed_format(align_bed_test_simple(), bed_columns=9) catch {
    AlignBedError(message) => abort("valid BED writing failed: " + message)
  }
  assert_eq(output, "chr7\t100\t180\tquery\t0\t+\t100\t180\t0\n")
}

///|
test "Bio.Align.bed reads BED4" {
  let alignment = align_bed_test_parse("chr1\t10\t20\ttx\n").alignments[0]
  assert_eq(alignment.source_columns, 4)
  assert_eq(alignment.query_id, Some("tx"))
}

///|
test "Bio.Align.bed reads BED5" {
  let alignment = align_bed_test_parse("chr1\t10\t20\ttx\t7\n").alignments[0]
  assert_eq(alignment.source_columns, 5)
  assert_eq(alignment.score.unwrap().number(), Some(7.0))
}

///|
test "Bio.Align.bed reads BED6 reverse strand" {
  let alignment = align_bed_test_parse("chr1\t10\t20\ttx\t7\t-\n").alignments[0]
  assert_eq(alignment.coordinates[0].query, 10)
  assert_eq(alignment.coordinates[1].query, 0)
}

///|
test "Bio.Align.bed reads BED7" {
  let alignment = align_bed_test_parse("chr1\t10\t20\ttx\t7\t+\t12\n").alignments[0]
  assert_eq(alignment.thick_start, Some(12))
  assert_eq(alignment.thick_end, None)
}

///|
test "Bio.Align.bed reads BED8" {
  let alignment = align_bed_test_parse("chr1\t10\t20\ttx\t7\t+\t12\t18\n").alignments[0]
  assert_eq(alignment.thick_end, Some(18))
}

///|
test "Bio.Align.bed reads BED9" {
  let alignment = align_bed_test_parse(
      "chr1\t10\t20\ttx\t7\t+\t12\t18\t10,20,30\n",
    ).alignments[0]
  assert_eq(alignment.item_rgb, Some("10,20,30"))
}

///|
test "Bio.Align.bed reads lossless single-block BED10" {
  let alignment = align_bed_test_parse("chr1\t10\t20\ttx\t7\t+\t12\t18\t0\t1\n").alignments[0]
  assert_eq(alignment.blocks().length(), 1)
  assert_eq(alignment.source_columns, 10)
}

///|
test "Bio.Align.bed reads inferable contiguous BED11" {
  let alignment = align_bed_test_parse(
      "chr1\t10\t20\ttx\t7\t+\t12\t18\t0\t2\t4,6,\n",
    ).alignments[0]
  assert_eq(alignment.blocks().length(), 2)
  assert_eq(alignment.blocks()[1].target_start, 14)
}

///|
test "Bio.Align.bed accepts BED12 lists without trailing commas" {
  let alignment = align_bed_test_parse(
      "chr1\t10\t20\ttx\t7\t+\t12\t18\t0\t2\t4,3\t0,7\n",
    ).alignments[0]
  assert_eq(alignment.query_size(), 7)
}

///|
test "Bio.Align.bed rejects fewer than three columns" {
  assert_true(align_bed_test_rejects("chr1\t0\n"))
}

///|
test "Bio.Align.bed rejects more than twelve columns" {
  assert_true(
    align_bed_test_rejects(
      "chr1\t0\t10\ttx\t0\t+\t0\t10\t0\t1\t10,\t0,\textra\n",
    ),
  )
}

///|
test "Bio.Align.bed rejects interior blank records" {
  assert_true(align_bed_test_rejects("chr1\t0\t10\n\nchr2\t0\t10\n"))
}

///|
test "Bio.Align.bed rejects comment records" {
  assert_true(align_bed_test_rejects("# comment\nchr1\t0\t10\n"))
}

///|
test "Bio.Align.bed rejects track records" {
  assert_true(align_bed_test_rejects("track name=test\nchr1\t0\t10\n"))
}

///|
test "Bio.Align.bed rejects negative chromStart" {
  assert_true(align_bed_test_rejects("chr1\t-1\t10\n"))
}

///|
test "Bio.Align.bed rejects reversed target interval" {
  assert_true(align_bed_test_rejects("chr1\t20\t10\n"))
}

///|
test "Bio.Align.bed rejects empty target interval" {
  assert_true(align_bed_test_rejects("chr1\t10\t10\n"))
}

///|
test "Bio.Align.bed rejects integer overflow" {
  assert_true(align_bed_test_rejects("chr1\t0\t2147483648\n"))
}

///|
test "Bio.Align.bed rejects invalid strand" {
  assert_true(align_bed_test_rejects("chr1\t0\t10\ttx\t0\tx\n"))
}

///|
test "Bio.Align.bed rejects thickStart before interval" {
  assert_true(align_bed_test_rejects("chr1\t10\t20\ttx\t0\t+\t9\n"))
}

///|
test "Bio.Align.bed rejects thickEnd after interval" {
  assert_true(align_bed_test_rejects("chr1\t10\t20\ttx\t0\t+\t10\t21\n"))
}

///|
test "Bio.Align.bed rejects reversed thick interval" {
  assert_true(align_bed_test_rejects("chr1\t10\t20\ttx\t0\t+\t18\t12\n"))
}

///|
test "Bio.Align.bed rejects multi-block BED10" {
  assert_true(align_bed_test_rejects("chr1\t0\t20\ttx\t0\t+\t0\t20\t0\t2\n"))
}

///|
test "Bio.Align.bed rejects BED11 block count mismatch" {
  assert_true(
    align_bed_test_rejects("chr1\t0\t20\ttx\t0\t+\t0\t20\t0\t2\t20,\n"),
  )
}

///|
test "Bio.Align.bed rejects ambiguous gapped BED11" {
  assert_true(
    align_bed_test_rejects("chr1\t0\t20\ttx\t0\t+\t0\t20\t0\t2\t5,5,\n"),
  )
}

///|
test "Bio.Align.bed rejects zero block count" {
  assert_true(
    align_bed_test_rejects("chr1\t0\t20\ttx\t0\t+\t0\t20\t0\t0\t\t\n"),
  )
}

///|
test "Bio.Align.bed rejects blockSizes count mismatch" {
  assert_true(
    align_bed_test_rejects("chr1\t0\t20\ttx\t0\t+\t0\t20\t0\t2\t5,\t0,15,\n"),
  )
}

///|
test "Bio.Align.bed rejects blockStarts count mismatch" {
  assert_true(
    align_bed_test_rejects("chr1\t0\t20\ttx\t0\t+\t0\t20\t0\t2\t5,5,\t0,\n"),
  )
}

///|
test "Bio.Align.bed rejects empty block list item" {
  assert_true(
    align_bed_test_rejects("chr1\t0\t20\ttx\t0\t+\t0\t20\t0\t2\t5,,5,\t0,15,\n"),
  )
}

///|
test "Bio.Align.bed rejects noninteger block size" {
  assert_true(
    align_bed_test_rejects("chr1\t0\t20\ttx\t0\t+\t0\t20\t0\t2\t5,x,\t0,15,\n"),
  )
}

///|
test "Bio.Align.bed rejects zero block size" {
  assert_true(
    align_bed_test_rejects("chr1\t0\t20\ttx\t0\t+\t0\t20\t0\t2\t0,5,\t0,15,\n"),
  )
}

///|
test "Bio.Align.bed rejects first block after chromStart" {
  assert_true(
    align_bed_test_rejects("chr1\t0\t20\ttx\t0\t+\t0\t20\t0\t2\t5,5,\t1,15,\n"),
  )
}

///|
test "Bio.Align.bed rejects overlapping blocks" {
  assert_true(
    align_bed_test_rejects("chr1\t0\t20\ttx\t0\t+\t0\t20\t0\t2\t10,10,\t0,5,\n"),
  )
}

///|
test "Bio.Align.bed rejects unsorted blocks" {
  assert_true(
    align_bed_test_rejects("chr1\t0\t20\ttx\t0\t+\t0\t20\t0\t2\t5,5,\t0,0,\n"),
  )
}

///|
test "Bio.Align.bed rejects block beyond chromEnd" {
  assert_true(
    align_bed_test_rejects("chr1\t0\t20\ttx\t0\t+\t0\t20\t0\t2\t5,10,\t0,15,\n"),
  )
}

///|
test "Bio.Align.bed rejects final block before chromEnd" {
  assert_true(
    align_bed_test_rejects("chr1\t0\t20\ttx\t0\t+\t0\t20\t0\t2\t5,4,\t0,15,\n"),
  )
}

///|
test "Bio.Align.bed rejects block coordinate overflow" {
  assert_true(
    align_bed_test_rejects(
      "chr1\t0\t2147483647\ttx\t0\t+\t0\t2147483647\t0\t1\t2147483647,\t1,\n",
    ),
  )
}

///|
test "Bio.Align.bed rejects invalid numeric score magnitude" {
  assert_true(
    try {
      ignore(@src.AlignBedScore::numeric(1.0e301))
      false
    } catch {
      AlignBedError(_) => true
    },
  )
}

///|
test "Bio.Align.bed rejects empty textual score" {
  assert_true(
    try {
      ignore(@src.AlignBedScore::text(""))
      false
    } catch {
      AlignBedError(_) => true
    },
  )
}

///|
test "Bio.Align.bed rejects whitespace in textual score" {
  assert_true(
    try {
      ignore(@src.AlignBedScore::text("not valid"))
      false
    } catch {
      AlignBedError(_) => true
    },
  )
}

///|
test "Bio.Align.bed coordinate constructor rejects negative target" {
  assert_true(
    try {
      ignore(@src.AlignBedCoordinate::create(-1, 0))
      false
    } catch {
      AlignBedError(_) => true
    },
  )
}

///|
test "Bio.Align.bed coordinate constructor rejects negative query" {
  assert_true(
    try {
      ignore(@src.AlignBedCoordinate::create(0, -1))
      false
    } catch {
      AlignBedError(_) => true
    },
  )
}

///|
test "Bio.Align.bed constructor rejects one coordinate" {
  assert_true(
    align_bed_test_constructor_rejects([align_bed_test_coordinate(0, 0)]),
  )
}

///|
test "Bio.Align.bed constructor rejects zero-length segment" {
  assert_true(
    align_bed_test_constructor_rejects([
      align_bed_test_coordinate(0, 0),
      align_bed_test_coordinate(0, 0),
    ]),
  )
}

///|
test "Bio.Align.bed constructor rejects nonmonotonic target" {
  assert_true(
    align_bed_test_constructor_rejects([
      align_bed_test_coordinate(0, 0),
      align_bed_test_coordinate(5, 5),
      align_bed_test_coordinate(3, 7),
    ]),
  )
}

///|
test "Bio.Align.bed constructor rejects nonmonotonic query" {
  assert_true(
    align_bed_test_constructor_rejects([
      align_bed_test_coordinate(0, 0),
      align_bed_test_coordinate(5, 5),
      align_bed_test_coordinate(10, 3),
    ]),
  )
}

///|
test "Bio.Align.bed constructor rejects unequal aligned steps" {
  assert_true(
    align_bed_test_constructor_rejects([
      align_bed_test_coordinate(0, 0),
      align_bed_test_coordinate(5, 4),
    ]),
  )
}

///|
test "Bio.Align.bed constructor rejects gaps-only path" {
  assert_true(
    align_bed_test_constructor_rejects([
      align_bed_test_coordinate(0, 0),
      align_bed_test_coordinate(5, 0),
    ]),
  )
}

///|
test "Bio.Align.bed constructor rejects empty target identifier" {
  assert_true(
    try {
      ignore(
        @src.AlignBedAlignment::create("", Some("query"), [
          align_bed_test_coordinate(0, 0),
          align_bed_test_coordinate(5, 5),
        ]),
      )
      false
    } catch {
      AlignBedError(_) => true
    },
  )
}

///|
test "Bio.Align.bed constructor allows absent query identifier" {
  let alignment = @src.AlignBedAlignment::create("target", None, [
    align_bed_test_coordinate(0, 0),
    align_bed_test_coordinate(5, 5),
  ]) catch {
    AlignBedError(message) => abort("valid BED alignment failed: " + message)
  }
  assert_eq(alignment.query_id, None)
}

///|
test "Bio.Align.bed constructor rejects invalid source column count" {
  assert_true(
    try {
      ignore(
        @src.AlignBedAlignment::create(
          "target",
          Some("query"),
          [align_bed_test_coordinate(0, 0), align_bed_test_coordinate(5, 5)],
          source_columns=13,
        ),
      )
      false
    } catch {
      AlignBedError(_) => true
    },
  )
}

///|
test "Bio.Align.bed constructor rejects thickStart outside interval" {
  assert_true(
    try {
      ignore(
        @src.AlignBedAlignment::create(
          "target",
          Some("query"),
          [align_bed_test_coordinate(10, 0), align_bed_test_coordinate(20, 10)],
          thick_start=Some(9),
        ),
      )
      false
    } catch {
      AlignBedError(_) => true
    },
  )
}

///|
test "Bio.Align.bed constructor rejects itemRgb whitespace" {
  assert_true(
    try {
      ignore(
        @src.AlignBedAlignment::create(
          "target",
          Some("query"),
          [align_bed_test_coordinate(10, 0), align_bed_test_coordinate(20, 10)],
          item_rgb=Some("0 0 0"),
        ),
      )
      false
    } catch {
      AlignBedError(_) => true
    },
  )
}

///|
test "Bio.Align.bed alignment constructor copies coordinate array" {
  let coordinates = [
    align_bed_test_coordinate(0, 0),
    align_bed_test_coordinate(5, 5),
  ]
  let alignment = align_bed_test_construct(coordinates)
  ignore(coordinates.pop())
  assert_eq(alignment.coordinates.length(), 2)
}

///|
test "Bio.Align.bed document constructor copies alignment array" {
  let alignments = [align_bed_test_plus()]
  let document = @src.AlignBedDocument::create(alignments) catch {
    AlignBedError(message) => abort("valid BED document failed: " + message)
  }
  ignore(alignments.pop())
  assert_eq(document.alignments.length(), 1)
}

///|
test "Bio.Align.bed writer rejects invalid bedN" {
  assert_true(
    try {
      ignore(@src.align_bed_format(align_bed_test_plus(), bed_columns=2))
      false
    } catch {
      AlignBedError(_) => true
    },
  )
}

///|
test "Bio.Align.bed normalizes reversed target path when writing" {
  let alignment = align_bed_test_construct([
    align_bed_test_coordinate(20, 10),
    align_bed_test_coordinate(10, 0),
  ])
  let output = @src.align_bed_format(alignment, bed_columns=6) catch {
    AlignBedError(message) => abort("valid BED writing failed: " + message)
  }
  assert_eq(output, "target\t10\t20\tquery\t0\t+\n")
}

///|
test "Bio.Align.bed projects Biopython strand fixture into four blocks" {
  let alignment = align_bed_test_construct([
    align_bed_test_coordinate(0, 11),
    align_bed_test_coordinate(2, 11),
    align_bed_test_coordinate(2, 9),
    align_bed_test_coordinate(3, 8),
    align_bed_test_coordinate(4, 8),
    align_bed_test_coordinate(6, 6),
    align_bed_test_coordinate(6, 5),
    align_bed_test_coordinate(9, 2),
    align_bed_test_coordinate(10, 2),
    align_bed_test_coordinate(12, 0),
    align_bed_test_coordinate(15, 0),
  ])
  assert_eq(alignment.blocks().length(), 4)
}

///|
test "Bio.Align.bed writes Biopython strand fixture" {
  let score = @src.AlignBedScore::numeric(8.0) catch {
    AlignBedError(message) => abort("valid BED score failed: " + message)
  }
  let alignment = @src.AlignBedAlignment::create(
    "target",
    Some("query"),
    [
      align_bed_test_coordinate(0, 11),
      align_bed_test_coordinate(2, 11),
      align_bed_test_coordinate(2, 9),
      align_bed_test_coordinate(3, 8),
      align_bed_test_coordinate(4, 8),
      align_bed_test_coordinate(6, 6),
      align_bed_test_coordinate(6, 5),
      align_bed_test_coordinate(9, 2),
      align_bed_test_coordinate(10, 2),
      align_bed_test_coordinate(12, 0),
      align_bed_test_coordinate(15, 0),
    ],
    score=Some(score),
  ) catch {
    AlignBedError(message) => abort("valid BED fixture failed: " + message)
  }
  let output = @src.align_bed_format(alignment) catch {
    AlignBedError(message) => abort("valid BED writing failed: " + message)
  }
  assert_eq(
    output, "target\t2\t12\tquery\t8\t-\t2\t12\t0\t4\t1,2,3,2,\t0,2,4,8,\n",
  )
}

///|
test "Bio.Align.bed counts query-only skips in constructed path" {
  let alignment = align_bed_test_construct([
    align_bed_test_coordinate(0, 0),
    align_bed_test_coordinate(2, 2),
    align_bed_test_coordinate(2, 4),
    align_bed_test_coordinate(5, 7),
  ])
  let counts = alignment.counts()
  assert_eq(counts.query_skip_bases, 2)
  assert_eq(counts.query_skip_opens, 1)
}

///|
test "Bio.Align.bed stable document round trip by canonical text" {
  let first = @src.align_bed_write(align_bed_test_document()) catch {
    AlignBedError(message) => abort("valid BED writing failed: " + message)
  }
  let second = @src.align_bed_write(align_bed_test_parse(first)) catch {
    AlignBedError(message) => abort("valid BED writing failed: " + message)
  }
  assert_eq(second, first)
}
