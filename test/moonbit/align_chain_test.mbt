///|
fn align_chain_test_parse(text : String) -> @src.AlignChainAlignment {
  @src.align_chain_parse(text) catch {
    AlignChainError(message) => abort("valid chain input failed: " + message)
  }
}

///|
fn align_chain_test_parse_all(text : String) -> Array[@src.AlignChainAlignment] {
  @src.align_chain_parse_all(text) catch {
    AlignChainError(message) => abort("valid chain input failed: " + message)
  }
}

///|
fn align_chain_test_raises(text : String) -> Bool {
  try {
    ignore(@src.align_chain_parse(text))
    false
  } catch {
    AlignChainError(_) => true
  }
}

///|
fn align_chain_test_all_raises(text : String) -> Bool {
  try {
    ignore(@src.align_chain_parse_all(text))
    false
  } catch {
    AlignChainError(_) => true
  }
}

///|
fn align_chain_test_example() -> @src.AlignChainAlignment {
  let alignments = align_chain_test_parse_all(@src.align_chain_example_text())
  alignments[0]
}

///|
fn align_chain_test_reverse() -> @src.AlignChainAlignment {
  let alignments = align_chain_test_parse_all(@src.align_chain_example_text())
  alignments[1]
}

///|
fn align_chain_test_coordinate(
  target : Int,
  query : Int,
) -> @src.AlignChainCoordinate {
  @src.AlignChainCoordinate::create(target, query) catch {
    AlignChainError(message) => abort("valid coordinate failed: " + message)
  }
}

///|
fn align_chain_test_block(
  size : Int,
  target_gap : Int,
  query_gap : Int,
) -> @src.AlignChainBlock {
  @src.AlignChainBlock::create(size, target_gap~, query_gap~) catch {
    AlignChainError(message) => abort("valid block failed: " + message)
  }
}

///|
fn align_chain_test_constructor_raises(
  target_size : Int,
  query_size : Int,
  coordinates : Array[@src.AlignChainCoordinate],
) -> Bool {
  try {
    ignore(
      @src.AlignChainAlignment::create(
        "target", target_size, "query", query_size, coordinates,
      ),
    )
    false
  } catch {
    AlignChainError(_) => true
  }
}

///|
fn align_chain_test_from_blocks_raises(
  target_start : Int,
  target_end : Int,
  query_start : Int,
  query_end : Int,
  blocks : Array[@src.AlignChainBlock],
) -> Bool {
  try {
    ignore(
      @src.align_chain_from_blocks(
        "target", 100, "+", target_start, target_end, "query", 100, "+", query_start,
        query_end, blocks,
      ),
    )
    false
  } catch {
    AlignChainError(_) => true
  }
}

///|
test "align_chain example contains two records" {
  let alignments = align_chain_test_parse_all(@src.align_chain_example_text())
  assert_eq(alignments.length(), 2)
}

///|
test "align_chain official target identifier" {
  assert_eq(align_chain_test_example().target_id, "chr3")
}

///|
test "align_chain official target size" {
  assert_eq(align_chain_test_example().target_size, 198295559)
}

///|
test "align_chain official query identifier" {
  assert_eq(align_chain_test_example().query_id, "NR_046654.1")
}

///|
test "align_chain official query size" {
  assert_eq(align_chain_test_example().query_size, 181)
}

///|
test "align_chain official score" {
  assert_eq(align_chain_test_example().score, 176.0)
}

///|
test "align_chain official chain identifier" {
  assert_eq(align_chain_test_example().chain_id, "1")
}

///|
test "align_chain official target strand" {
  assert_eq(align_chain_test_example().target_strand(), "+")
}

///|
test "align_chain official query strand" {
  assert_eq(align_chain_test_example().query_strand(), "-")
}

///|
test "align_chain official coordinate count" {
  assert_eq(align_chain_test_example().coordinates.length(), 6)
}

///|
test "align_chain official first coordinate" {
  let coordinate = align_chain_test_example().coordinates[0]
  assert_eq(coordinate.target, 42530895)
  assert_eq(coordinate.query, 181)
}

///|
test "align_chain official first exon end coordinate" {
  let coordinate = align_chain_test_example().coordinates[1]
  assert_eq(coordinate.target, 42530958)
  assert_eq(coordinate.query, 118)
}

///|
test "align_chain official second exon start coordinate" {
  let coordinate = align_chain_test_example().coordinates[2]
  assert_eq(coordinate.target, 42532020)
  assert_eq(coordinate.query, 118)
}

///|
test "align_chain official second exon end coordinate" {
  let coordinate = align_chain_test_example().coordinates[3]
  assert_eq(coordinate.target, 42532095)
  assert_eq(coordinate.query, 43)
}

///|
test "align_chain official third exon start coordinate" {
  let coordinate = align_chain_test_example().coordinates[4]
  assert_eq(coordinate.target, 42532563)
  assert_eq(coordinate.query, 43)
}

///|
test "align_chain official final coordinate" {
  let coordinate = align_chain_test_example().coordinates[5]
  assert_eq(coordinate.target, 42532606)
  assert_eq(coordinate.query, 0)
}

///|
test "align_chain official normalized target span" {
  let alignment = align_chain_test_example()
  assert_eq(alignment.target_start(), 42530895)
  assert_eq(alignment.target_end(), 42532606)
}

///|
test "align_chain official normalized query span" {
  let alignment = align_chain_test_example()
  assert_eq(alignment.query_start(), 0)
  assert_eq(alignment.query_end(), 181)
}

///|
test "align_chain official operation path" {
  assert_eq(align_chain_test_example().operation_path(), "MDMDM")
}

///|
test "align_chain official alignment length" {
  assert_eq(align_chain_test_example().alignment_length(), 1711)
}

///|
test "align_chain official aligned bases" {
  assert_eq(align_chain_test_example().aligned_bases(), 181)
}

///|
test "align_chain official block count" {
  assert_eq(align_chain_test_example().blocks().length(), 3)
}

///|
test "align_chain official first block" {
  let block = align_chain_test_example().blocks()[0]
  assert_eq(block.size, 63)
  assert_eq(block.target_gap, 1062)
  assert_eq(block.query_gap, 0)
}

///|
test "align_chain official second block" {
  let block = align_chain_test_example().blocks()[1]
  assert_eq(block.size, 75)
  assert_eq(block.target_gap, 468)
  assert_eq(block.query_gap, 0)
}

///|
test "align_chain official final block" {
  let block = align_chain_test_example().blocks()[2]
  assert_eq(block.size, 43)
  assert_eq(block.target_gap, 0)
  assert_eq(block.query_gap, 0)
}

///|
test "align_chain official count aligned" {
  assert_eq(align_chain_test_example().counts().aligned, 181)
}

///|
test "align_chain official count target gaps" {
  assert_eq(align_chain_test_example().counts().target_gap_bases, 1530)
}

///|
test "align_chain official count query gaps" {
  assert_eq(align_chain_test_example().counts().query_gap_bases, 0)
}

///|
test "align_chain official count gap columns" {
  assert_eq(align_chain_test_example().counts().gap_columns, 1530)
}

///|
test "align_chain official count gap opens" {
  let counts = align_chain_test_example().counts()
  assert_eq(counts.target_gap_opens, 2)
  assert_eq(counts.query_gap_opens, 0)
}

///|
test "align_chain official count aligned blocks" {
  assert_eq(align_chain_test_example().counts().aligned_blocks, 3)
}

///|
test "align_chain official count path segments" {
  assert_eq(align_chain_test_example().counts().path_segments, 5)
}

///|
test "align_chain official map first target residue" {
  assert_eq(align_chain_test_example().map_target_position(42530895), Some(180))
}

///|
test "align_chain official map first exon end residue" {
  assert_eq(align_chain_test_example().map_target_position(42530957), Some(118))
}

///|
test "align_chain official target intron is unmapped" {
  assert_eq(align_chain_test_example().map_target_position(42530958), None)
}

///|
test "align_chain official map second exon start" {
  assert_eq(align_chain_test_example().map_target_position(42532020), Some(117))
}

///|
test "align_chain official map second exon end" {
  assert_eq(align_chain_test_example().map_target_position(42532094), Some(43))
}

///|
test "align_chain official second intron is unmapped" {
  assert_eq(align_chain_test_example().map_target_position(42532095), None)
}

///|
test "align_chain official map final exon start" {
  assert_eq(align_chain_test_example().map_target_position(42532563), Some(42))
}

///|
test "align_chain official map final target residue" {
  assert_eq(align_chain_test_example().map_target_position(42532605), Some(0))
}

///|
test "align_chain official map first query residue" {
  assert_eq(align_chain_test_example().map_query_position(180), Some(42530895))
}

///|
test "align_chain official map last query residue" {
  assert_eq(align_chain_test_example().map_query_position(0), Some(42532605))
}

///|
test "align_chain official outside target returns none" {
  let alignment = align_chain_test_example()
  assert_eq(alignment.map_target_position(-1), None)
  assert_eq(alignment.map_target_position(alignment.target_size), None)
}

///|
test "align_chain official outside query returns none" {
  let alignment = align_chain_test_example()
  assert_eq(alignment.map_query_position(-1), None)
  assert_eq(alignment.map_query_position(181), None)
}

///|
test "align_chain official target range first exon" {
  let mappings = align_chain_test_example().map_target_range(42530890, 42530960)
  assert_eq(mappings.length(), 1)
  assert_eq(mappings[0].target_start, 42530895)
  assert_eq(mappings[0].target_end, 42530958)
  assert_eq(mappings[0].query_start, 118)
  assert_eq(mappings[0].query_end, 181)
}

///|
test "align_chain official target range all exons" {
  let alignment = align_chain_test_example()
  let mappings = alignment.map_target_range(
    alignment.target_start(),
    alignment.target_end(),
  )
  assert_eq(mappings.length(), 3)
}

///|
test "align_chain official query range all exons" {
  let mappings = align_chain_test_example().map_query_range(0, 181)
  assert_eq(mappings.length(), 3)
}

///|
test "align_chain official aligned pairs endpoints" {
  let pairs = align_chain_test_example().aligned_pairs()
  assert_eq(pairs.length(), 181)
  assert_eq(pairs[0].target, 42530895)
  assert_eq(pairs[0].query, 180)
  assert_eq(pairs[180].target, 42532605)
  assert_eq(pairs[180].query, 0)
}

///|
test "align_chain aligned pair limit" {
  let raised = try {
    ignore(align_chain_test_example().aligned_pairs(limit=180))
    false
  } catch {
    AlignChainError(_) => true
  }
  assert_true(raised)
}

///|
test "align_chain official summary" {
  let summary = align_chain_test_example().summary()
  assert_true(summary.contains("chr3"))
  assert_true(summary.contains("aligned=181"))
  assert_true(summary.contains("gaps=1530"))
}

///|
test "align_chain reverse target strand" {
  assert_eq(align_chain_test_reverse().target_strand(), "-")
}

///|
test "align_chain reverse query strand" {
  assert_eq(align_chain_test_reverse().query_strand(), "+")
}

///|
test "align_chain reverse absolute target span" {
  let alignment = align_chain_test_reverse()
  assert_eq(alignment.target_start(), 800)
  assert_eq(alignment.target_end(), 900)
}

///|
test "align_chain reverse query span" {
  let alignment = align_chain_test_reverse()
  assert_eq(alignment.query_start(), 20)
  assert_eq(alignment.query_end(), 130)
}

///|
test "align_chain reverse coordinates" {
  let alignment = align_chain_test_reverse()
  assert_eq(alignment.coordinates.length(), 5)
  assert_eq(alignment.coordinates[0].target, 900)
  assert_eq(alignment.coordinates[4].target, 800)
  assert_eq(alignment.coordinates[0].query, 20)
  assert_eq(alignment.coordinates[4].query, 130)
}

///|
test "align_chain reverse operation path" {
  assert_eq(align_chain_test_reverse().operation_path(), "MDIM")
}

///|
test "align_chain reverse alignment length" {
  assert_eq(align_chain_test_reverse().alignment_length(), 120)
}

///|
test "align_chain reverse aligned bases" {
  assert_eq(align_chain_test_reverse().aligned_bases(), 90)
}

///|
test "align_chain reverse target gap bases" {
  assert_eq(align_chain_test_reverse().counts().target_gap_bases, 10)
}

///|
test "align_chain reverse query gap bases" {
  assert_eq(align_chain_test_reverse().counts().query_gap_bases, 20)
}

///|
test "align_chain reverse gap opens" {
  let counts = align_chain_test_reverse().counts()
  assert_eq(counts.target_gap_opens, 1)
  assert_eq(counts.query_gap_opens, 1)
}

///|
test "align_chain reverse first target mapping" {
  assert_eq(align_chain_test_reverse().map_target_position(899), Some(20))
}

///|
test "align_chain reverse first block end mapping" {
  assert_eq(align_chain_test_reverse().map_target_position(860), Some(59))
}

///|
test "align_chain reverse target gap mapping" {
  assert_eq(align_chain_test_reverse().map_target_position(859), None)
}

///|
test "align_chain reverse second block mapping" {
  assert_eq(align_chain_test_reverse().map_target_position(849), Some(80))
}

///|
test "align_chain reverse final mapping" {
  assert_eq(align_chain_test_reverse().map_target_position(800), Some(129))
}

///|
test "align_chain reverse inserted query is unmapped" {
  assert_eq(align_chain_test_reverse().map_query_position(60), None)
  assert_eq(align_chain_test_reverse().map_query_position(79), None)
}

///|
test "align_chain reverse query mapping" {
  assert_eq(align_chain_test_reverse().map_query_position(80), Some(849))
}

///|
test "align_chain reverse range strand metadata" {
  let mappings = align_chain_test_reverse().map_target_range(800, 900)
  assert_eq(mappings.length(), 2)
  assert_eq(mappings[0].target_strand, "-")
  assert_eq(mappings[0].query_strand, "+")
}

///|
test "align_chain invert identifiers" {
  let inverted = align_chain_test_reverse().invert()
  assert_eq(inverted.target_id, "queryA")
  assert_eq(inverted.query_id, "chr1")
}

///|
test "align_chain invert coordinates" {
  let original = align_chain_test_reverse()
  let inverted = original.invert()
  assert_eq(inverted.coordinates[0].target, original.coordinates[0].query)
  assert_eq(inverted.coordinates[0].query, original.coordinates[0].target)
}

///|
test "align_chain invert twice roundtrip" {
  let original = align_chain_test_reverse()
  assert_eq(original.invert().invert(), original)
}

///|
test "align_chain official canonical writer" {
  let output = @src.align_chain_write(align_chain_test_example())
  assert_eq(
    output,
    "chain 176 chr3 198295559 + 42530895 42532606 NR_046654.1 181 - 0 181 1\n" +
    "63\t1062\t0\n" +
    "75\t468\t0\n" +
    "43\n\n",
  )
}

///|
test "align_chain reverse canonical writer" {
  let output = @src.align_chain_write(align_chain_test_reverse())
  assert_eq(
    output,
    "chain 500 chr1 1000 - 100 200 queryA 500 + 20 130 reverse-target\n" +
    "40\t10\t20\n" +
    "50\n\n",
  )
}

///|
test "align_chain canonical single roundtrip" {
  let original = align_chain_test_example()
  let reparsed = align_chain_test_parse(@src.align_chain_write(original))
  assert_eq(reparsed, original)
}

///|
test "align_chain canonical multiple roundtrip" {
  let original = align_chain_test_parse_all(@src.align_chain_example_text())
  let reparsed = align_chain_test_parse_all(
    @src.align_chain_write_all(original),
  )
  assert_eq(reparsed, original)
}

///|
test "align_chain no identifier header" {
  let alignment = align_chain_test_parse(
    "chain 10 target 100 + 10 20 query 100 + 30 40\n10\n",
  )
  assert_eq(alignment.chain_id, "")
  assert_true(@src.align_chain_write(alignment).contains("query 100 + 30 40\n"))
}

///|
test "align_chain decimal score" {
  let alignment = align_chain_test_parse(
    "chain 10.5 target 100 + 10 20 query 100 + 30 40 id\n10\n",
  )
  assert_eq(alignment.score, 10.5)
}

///|
test "align_chain exponent score" {
  let alignment = align_chain_test_parse(
    "chain 1e2 target 100 + 10 20 query 100 + 30 40 id\n10\n",
  )
  assert_eq(alignment.score, 100.0)
}

///|
test "align_chain explicit positive score" {
  let alignment = align_chain_test_parse(
    "chain +10 target 100 + 10 20 query 100 + 30 40 id\n10\n",
  )
  assert_eq(alignment.score, 10.0)
}

///|
test "align_chain negative score preserved" {
  let alignment = align_chain_test_parse(
    "chain -10 target 100 + 10 20 query 100 + 30 40 id\n10\n",
  )
  assert_eq(alignment.score, -10.0)
}

///|
test "align_chain CRLF input" {
  let alignment = align_chain_test_parse(
    "chain 10 target 100 + 10 20 query 100 + 30 40 id\r\n10\r\n",
  )
  assert_eq(alignment.aligned_bases(), 10)
}

///|
test "align_chain tab separated input" {
  let alignment = align_chain_test_parse(
    "chain\t10\ttarget\t100\t+\t10\t20\tquery\t100\t+\t30\t40\tid\n10\n",
  )
  assert_eq(alignment.chain_id, "id")
}

///|
test "align_chain leading and trailing blank lines" {
  let alignments = align_chain_test_parse_all(
    "\n\nchain 10 target 100 + 10 20 query 100 + 30 40 id\n10\n\n",
  )
  assert_eq(alignments.length(), 1)
}

///|
test "align_chain records need no blank separator" {
  let alignments = align_chain_test_parse_all(
    "chain 10 target 100 + 10 20 query 100 + 30 40 one\n10\n" +
    "chain 20 target 100 + 20 30 query 100 + 40 50 two\n10\n",
  )
  assert_eq(alignments.length(), 2)
}

///|
test "align_chain size zero before initial target gap" {
  let alignment = align_chain_test_parse(
    "chain 10 target 100 + 10 25 query 100 + 30 40 id\n" + "0\t5\t0\n" + "10\n",
  )
  assert_eq(alignment.operation_path(), "DM")
  assert_eq(alignment.blocks()[0].size, 0)
}

///|
test "align_chain trailing target gap and zero final block" {
  let alignment = align_chain_test_parse(
    "chain 10 target 100 + 10 25 query 100 + 30 40 id\n" + "10\t5\t0\n" + "0\n",
  )
  assert_eq(alignment.operation_path(), "MD")
  assert_eq(alignment.blocks()[1].size, 0)
}

///|
test "align_chain block with both gap directions" {
  let alignment = align_chain_test_parse(
    "chain 20 target 100 + 10 35 query 100 + 20 47 id\n" + "10\t5\t7\n" + "10\n",
  )
  assert_eq(alignment.operation_path(), "MDIM")
}

///|
test "align_chain all strand combinations plus plus" {
  let alignment = align_chain_test_parse(
    "chain 10 target 100 + 10 20 query 100 + 30 40 pp\n10\n",
  )
  assert_eq(alignment.coordinates[0].target, 10)
  assert_eq(alignment.coordinates[0].query, 30)
}

///|
test "align_chain all strand combinations plus minus" {
  let alignment = align_chain_test_parse(
    "chain 10 target 100 + 10 20 query 100 - 30 40 pm\n10\n",
  )
  assert_eq(alignment.coordinates[0].target, 10)
  assert_eq(alignment.coordinates[0].query, 70)
  assert_eq(alignment.coordinates[1].query, 60)
}

///|
test "align_chain all strand combinations minus plus" {
  let alignment = align_chain_test_parse(
    "chain 10 target 100 - 10 20 query 100 + 30 40 mp\n10\n",
  )
  assert_eq(alignment.coordinates[0].target, 90)
  assert_eq(alignment.coordinates[1].target, 80)
  assert_eq(alignment.coordinates[0].query, 30)
}

///|
test "align_chain all strand combinations minus minus" {
  let alignment = align_chain_test_parse(
    "chain 10 target 100 - 10 20 query 100 - 30 40 mm\n10\n",
  )
  assert_eq(alignment.coordinates[0].target, 90)
  assert_eq(alignment.coordinates[0].query, 70)
  assert_eq(alignment.map_target_position(89), Some(69))
}

///|
test "align_chain constructor copies coordinates" {
  let coordinates : Array[@src.AlignChainCoordinate] = [
    align_chain_test_coordinate(10, 20),
    align_chain_test_coordinate(20, 30),
  ]
  let alignment = @src.AlignChainAlignment::create(
    "target", 100, "query", 100, coordinates,
  )
  coordinates.push(align_chain_test_coordinate(30, 40))
  assert_eq(alignment.coordinates.length(), 2)
}

///|
test "align_chain coordinate constructor" {
  let coordinate = @src.AlignChainCoordinate::create(12, 34)
  assert_eq(coordinate.target, 12)
  assert_eq(coordinate.query, 34)
}

///|
test "align_chain block constructor defaults" {
  let block = @src.AlignChainBlock::create(10)
  assert_eq(block.size, 10)
  assert_eq(block.target_gap, 0)
  assert_eq(block.query_gap, 0)
}

///|
test "align_chain from blocks constructor" {
  let alignment = @src.align_chain_from_blocks(
    "target",
    100,
    "+",
    10,
    35,
    "query",
    100,
    "+",
    20,
    47,
    [align_chain_test_block(10, 5, 7), align_chain_test_block(10, 0, 0)],
    score=42.0,
    chain_id="made",
  )
  assert_eq(alignment.operation_path(), "MDIM")
  assert_eq(alignment.score, 42.0)
  assert_eq(alignment.chain_id, "made")
}

///|
test "align_chain adjacent match segments canonicalize" {
  let alignment = @src.AlignChainAlignment::create(
    "target",
    100,
    "query",
    100,
    [
      align_chain_test_coordinate(10, 20),
      align_chain_test_coordinate(15, 25),
      align_chain_test_coordinate(20, 30),
    ],
  )
  let blocks = alignment.blocks()
  assert_eq(blocks.length(), 1)
  assert_eq(blocks[0].size, 10)
}

///|
test "align_chain find by identifier" {
  let alignments = align_chain_test_parse_all(@src.align_chain_example_text())
  let found = @src.align_chain_find_by_id(alignments, "reverse-target")
  assert_true(found is Some(_))
  match found {
    Some(alignment) => assert_eq(alignment.target_id, "chr1")
    None => ()
  }
}

///|
test "align_chain missing identifier lookup" {
  let alignments = align_chain_test_parse_all(@src.align_chain_example_text())
  assert_eq(@src.align_chain_find_by_id(alignments, "missing"), None)
}

///|
test "align_chain target overlap query" {
  let alignments = align_chain_test_parse_all(@src.align_chain_example_text())
  let matches = @src.align_chain_target_overlaps(
    alignments, "chr3", 42532000, 42532100,
  )
  assert_eq(matches.length(), 1)
}

///|
test "align_chain target overlap excludes touching endpoint" {
  let alignments = align_chain_test_parse_all(@src.align_chain_example_text())
  let matches = @src.align_chain_target_overlaps(
    alignments, "chr3", 42532606, 42532607,
  )
  assert_eq(matches.length(), 0)
}

///|
test "align_chain empty parse all" {
  assert_eq(align_chain_test_parse_all("").length(), 0)
}

///|
test "align_chain whitespace parse all" {
  assert_eq(align_chain_test_parse_all(" \n\t\n").length(), 0)
}

///|
test "align_chain parse requires one record" {
  assert_true(align_chain_test_raises(""))
}

///|
test "align_chain parse rejects multiple records" {
  assert_true(align_chain_test_raises(@src.align_chain_example_text()))
}

///|
test "align_chain rejects content before header" {
  assert_true(
    align_chain_test_all_raises(
      "comment\nchain 10 target 100 + 10 20 query 100 + 30 40 id\n10\n",
    ),
  )
}

///|
test "align_chain rejects short header" {
  assert_true(
    align_chain_test_raises("chain 10 target 100 + 10 20 query 100 + 30\n10\n"),
  )
}

///|
test "align_chain rejects long header" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 20 query 100 + 30 40 id extra\n10\n",
    ),
  )
}

///|
test "align_chain rejects wrong keyword" {
  assert_true(
    align_chain_test_raises(
      "chains 10 target 100 + 10 20 query 100 + 30 40 id\n10\n",
    ),
  )
}

///|
test "align_chain rejects invalid score" {
  assert_true(
    align_chain_test_raises(
      "chain score target 100 + 10 20 query 100 + 30 40 id\n10\n",
    ),
  )
}

///|
test "align_chain rejects incomplete exponent score" {
  assert_true(
    align_chain_test_raises(
      "chain 1e target 100 + 10 20 query 100 + 30 40 id\n10\n",
    ),
  )
}

///|
test "align_chain rejects non-finite score" {
  assert_true(
    align_chain_test_raises(
      "chain 1e999 target 100 + 10 20 query 100 + 30 40 id\n10\n",
    ),
  )
}

///|
test "align_chain rejects invalid target size" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target size + 10 20 query 100 + 30 40 id\n10\n",
    ),
  )
}

///|
test "align_chain rejects integer overflow" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 999999999999 + 10 20 query 100 + 30 40 id\n10\n",
    ),
  )
}

///|
test "align_chain rejects negative coordinate" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + -1 20 query 100 + 30 40 id\n10\n",
    ),
  )
}

///|
test "align_chain rejects target start after end" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 20 10 query 100 + 30 40 id\n10\n",
    ),
  )
}

///|
test "align_chain rejects target end beyond size" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 90 110 query 100 + 30 50 id\n20\n",
    ),
  )
}

///|
test "align_chain rejects query start after end" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 20 query 100 + 40 30 id\n10\n",
    ),
  )
}

///|
test "align_chain rejects query end beyond size" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 30 query 100 + 90 110 id\n20\n",
    ),
  )
}

///|
test "align_chain rejects target strand" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 ? 10 20 query 100 + 30 40 id\n10\n",
    ),
  )
}

///|
test "align_chain rejects query strand" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 20 query 100 ? 30 40 id\n10\n",
    ),
  )
}

///|
test "align_chain rejects missing block" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 20 query 100 + 30 40 id\n",
    ),
  )
}

///|
test "align_chain rejects blank before final block" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 20 query 100 + 30 40 id\n\n10\n",
    ),
  )
}

///|
test "align_chain rejects next header before final block" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 20 query 100 + 30 40 one\n" +
      "chain 10 target 100 + 20 30 query 100 + 40 50 two\n10\n",
    ),
  )
}

///|
test "align_chain rejects two-field block" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 20 query 100 + 30 40 id\n5 5\n5\n",
    ),
  )
}

///|
test "align_chain rejects four-field block" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 20 query 100 + 30 40 id\n5 1 1 1\n5\n",
    ),
  )
}

///|
test "align_chain rejects invalid block integer" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 20 query 100 + 30 40 id\nfive\n",
    ),
  )
}

///|
test "align_chain rejects negative block integer" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 20 query 100 + 30 40 id\n-10\n",
    ),
  )
}

///|
test "align_chain rejects zero gaps on non-final block" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 20 query 100 + 30 40 id\n5 0 0\n5\n",
    ),
  )
}

///|
test "align_chain rejects target span mismatch" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 25 query 100 + 30 40 id\n10\n",
    ),
  )
}

///|
test "align_chain rejects query span mismatch" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 20 query 100 + 30 45 id\n10\n",
    ),
  )
}

///|
test "align_chain rejects blocks exceeding target span" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 20 query 100 + 30 40 id\n11\n",
    ),
  )
}

///|
test "align_chain rejects blocks exceeding query span" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 20 query 100 + 30 39 id\n10\n",
    ),
  )
}

///|
test "align_chain rejects all-gap record" {
  assert_true(
    align_chain_test_raises(
      "chain 10 target 100 + 10 15 query 100 + 30 30 id\n0 5 0\n0\n",
    ),
  )
}

///|
test "align_chain coordinate constructor rejects negative target" {
  let raised = try {
    ignore(@src.AlignChainCoordinate::create(-1, 0))
    false
  } catch {
    AlignChainError(_) => true
  }
  assert_true(raised)
}

///|
test "align_chain coordinate constructor rejects negative query" {
  let raised = try {
    ignore(@src.AlignChainCoordinate::create(0, -1))
    false
  } catch {
    AlignChainError(_) => true
  }
  assert_true(raised)
}

///|
test "align_chain block constructor rejects negative size" {
  let raised = try {
    ignore(@src.AlignChainBlock::create(-1))
    false
  } catch {
    AlignChainError(_) => true
  }
  assert_true(raised)
}

///|
test "align_chain block constructor rejects negative target gap" {
  let raised = try {
    ignore(@src.AlignChainBlock::create(1, target_gap=-1))
    false
  } catch {
    AlignChainError(_) => true
  }
  assert_true(raised)
}

///|
test "align_chain constructor rejects empty target identifier" {
  let raised = try {
    ignore(
      @src.AlignChainAlignment::create("", 100, "query", 100, [
        align_chain_test_coordinate(10, 20),
        align_chain_test_coordinate(20, 30),
      ]),
    )
    false
  } catch {
    AlignChainError(_) => true
  }
  assert_true(raised)
}

///|
test "align_chain constructor rejects identifier whitespace" {
  let raised = try {
    ignore(
      @src.AlignChainAlignment::create("bad target", 100, "query", 100, [
        align_chain_test_coordinate(10, 20),
        align_chain_test_coordinate(20, 30),
      ]),
    )
    false
  } catch {
    AlignChainError(_) => true
  }
  assert_true(raised)
}

///|
test "align_chain constructor rejects negative target size" {
  assert_true(
    align_chain_test_constructor_raises(-1, 100, [
      align_chain_test_coordinate(0, 0),
      align_chain_test_coordinate(10, 10),
    ]),
  )
}

///|
test "align_chain constructor rejects coordinate outside target" {
  assert_true(
    align_chain_test_constructor_raises(10, 100, [
      align_chain_test_coordinate(0, 0),
      align_chain_test_coordinate(11, 11),
    ]),
  )
}

///|
test "align_chain constructor rejects coordinate outside query" {
  assert_true(
    align_chain_test_constructor_raises(100, 10, [
      align_chain_test_coordinate(0, 0),
      align_chain_test_coordinate(11, 11),
    ]),
  )
}

///|
test "align_chain constructor rejects one point" {
  assert_true(
    align_chain_test_constructor_raises(100, 100, [
      align_chain_test_coordinate(0, 0),
    ]),
  )
}

///|
test "align_chain constructor rejects duplicate points" {
  assert_true(
    align_chain_test_constructor_raises(100, 100, [
      align_chain_test_coordinate(10, 10),
      align_chain_test_coordinate(10, 10),
    ]),
  )
}

///|
test "align_chain constructor rejects unequal diagonal" {
  assert_true(
    align_chain_test_constructor_raises(100, 100, [
      align_chain_test_coordinate(10, 10),
      align_chain_test_coordinate(20, 25),
    ]),
  )
}

///|
test "align_chain constructor rejects target direction change" {
  assert_true(
    align_chain_test_constructor_raises(100, 100, [
      align_chain_test_coordinate(10, 10),
      align_chain_test_coordinate(20, 20),
      align_chain_test_coordinate(15, 20),
    ]),
  )
}

///|
test "align_chain constructor rejects query direction change" {
  assert_true(
    align_chain_test_constructor_raises(100, 100, [
      align_chain_test_coordinate(10, 10),
      align_chain_test_coordinate(20, 20),
      align_chain_test_coordinate(20, 15),
    ]),
  )
}

///|
test "align_chain constructor rejects gaps-only path" {
  assert_true(
    align_chain_test_constructor_raises(100, 100, [
      align_chain_test_coordinate(10, 10),
      align_chain_test_coordinate(20, 10),
      align_chain_test_coordinate(20, 30),
    ]),
  )
}

///|
test "align_chain from blocks rejects empty blocks" {
  assert_true(align_chain_test_from_blocks_raises(10, 20, 30, 40, []))
}

///|
test "align_chain from blocks rejects non-final zero gaps" {
  assert_true(
    align_chain_test_from_blocks_raises(10, 20, 30, 40, [
      align_chain_test_block(5, 0, 0),
      align_chain_test_block(5, 0, 0),
    ]),
  )
}

///|
test "align_chain from blocks rejects final gaps" {
  assert_true(
    align_chain_test_from_blocks_raises(10, 20, 30, 40, [
      align_chain_test_block(5, 5, 0),
    ]),
  )
}

///|
test "align_chain invalid target range" {
  let raised = try {
    ignore(align_chain_test_example().map_target_range(10, 10))
    false
  } catch {
    AlignChainError(_) => true
  }
  assert_true(raised)
}

///|
test "align_chain invalid query range" {
  let raised = try {
    ignore(align_chain_test_example().map_query_range(-1, 10))
    false
  } catch {
    AlignChainError(_) => true
  }
  assert_true(raised)
}

///|
test "align_chain invalid overlap query" {
  let raised = try {
    ignore(@src.align_chain_target_overlaps([], "chr1", 10, 10))
    false
  } catch {
    AlignChainError(_) => true
  }
  assert_true(raised)
}

///|
test "align_chain negative aligned pair limit" {
  let raised = try {
    ignore(align_chain_test_example().aligned_pairs(limit=-1))
    false
  } catch {
    AlignChainError(_) => true
  }
  assert_true(raised)
}
