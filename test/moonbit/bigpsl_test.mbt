///|
fn bigpsl_test_data() -> (Array[@src.BigBedTarget], Array[@src.BigPslAlignment]) {
  @src.bigpsl_example_data() catch {
    _ => abort("BigPsl example data should be valid")
  }
}

///|
fn bigpsl_test_file(
  compress? : Bool = true,
  items_per_slot? : Int = 2,
) -> @src.BigPslFile {
  let (targets, alignments) = bigpsl_test_data()
  let config = @src.BigPslWriteConfig::create(
    compress~,
    block_size=3,
    items_per_slot~,
  ) catch {
    _ => abort("BigPsl test configuration should be valid")
  }
  let bytes = @src.bigpsl_write(targets, alignments, config~) catch {
    _ => abort("BigPsl example should serialize")
  }
  @src.bigpsl_parse(bytes) catch {
    _ => abort("BigPsl example should parse")
  }
}

///|
fn bigpsl_test_alignment(
  target_coordinates? : Array[Int] = [10, 20],
  query_coordinates? : Array[Int] = [0, 10],
  sequence_type? : @src.BigPslSequenceType = @src.BigPslSequenceType::BigPslNucleotide,
) -> @src.BigPslAlignment {
  @src.BigPslAlignment::create(
    "chr1",
    100,
    "query",
    20,
    target_coordinates,
    query_coordinates,
    sequence_type~,
  ) catch {
    _ => abort("BigPsl test alignment should be valid")
  }
}

///|
fn bigpsl_test_custom_fields() -> Array[@src.BigBedField] {
  let schema = @src.bigpsl_schema()
  let fields : Array[@src.BigBedField] = []
  for index in 12..<schema.fields.length() {
    fields.push(schema.fields[index])
  }
  fields
}

///|
fn bigpsl_test_raw_bytes(
  extra_fields : Array[String],
  strand? : String = "+",
  block_sizes? : Array[Int] = [10],
  block_starts? : Array[Int] = [0],
  chrom_start? : Int = 10,
  chrom_end? : Int = 20,
) -> Array[Int] {
  let target = @src.BigBedTarget::create("chr1", 100) catch {
    _ => abort("raw BigPsl target should be valid")
  }
  let record = @src.BigBedRecord::create(
    "chr1",
    chrom_start,
    chrom_end,
    name="query",
    score=500,
    strand~,
    block_count=block_sizes.length(),
    block_sizes~,
    block_starts~,
    extra_fields~,
  ) catch {
    _ => abort("raw BigPsl record should be valid BigBed")
  }
  let config = @src.BigBedWriteConfig::create(
    bed_columns=12,
    schema_name="bigPsl",
    schema_comment="bigPsl pairwise alignment",
    custom_fields=bigpsl_test_custom_fields(),
  ) catch {
    _ => abort("raw BigPsl schema should be valid")
  }
  @src.bigbed_write([target], [record], config~) catch {
    _ => abort("raw BigPsl BigBed should serialize")
  }
}

///|
fn bigpsl_test_standard_extra() -> Array[String] {
  ["0", "10", "+", "20", "0", "", "", "100", "10", "0", "0", "0", "1"]
}

///|
test "bigpsl: default configuration matches Biopython defaults" {
  let config = @src.BigPslWriteConfig::default()
  assert_true(config.compress)
  assert_eq(config.block_size, 256)
  assert_eq(config.items_per_slot, 512)
  assert_true(!config.store_query_sequence)
  assert_true(!config.store_cds)
}

///|
test "bigpsl: custom configuration is preserved" {
  let config = @src.BigPslWriteConfig::create(
    compress=false,
    block_size=8,
    items_per_slot=4,
    store_query_sequence=true,
    store_cds=true,
  ) catch {
    _ => abort("custom configuration should be valid")
  }
  assert_true(!config.compress)
  assert_eq(config.block_size, 8)
  assert_eq(config.items_per_slot, 4)
  assert_true(config.store_query_sequence)
  assert_true(config.store_cds)
}

///|
test "bigpsl: configuration rejects blockSize below two" {
  let rejected = try {
    ignore(@src.BigPslWriteConfig::create(block_size=1))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: configuration rejects zero itemsPerSlot" {
  let rejected = try {
    ignore(@src.BigPslWriteConfig::create(items_per_slot=0))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: standard schema is bed12 plus thirteen fields" {
  let schema = @src.bigpsl_schema()
  assert_eq(schema.name, "bigPsl")
  assert_eq(schema.comment, "bigPsl pairwise alignment")
  assert_eq(schema.fields.length(), 25)
  assert_eq(schema.fields[11].name, "chromStarts")
  assert_eq(schema.fields[12].name, "oChromStart")
  assert_eq(schema.fields[24].name, "seqType")
}

///|
test "bigpsl: standard schema preserves array and long-string types" {
  let schema = @src.bigpsl_schema()
  assert_eq(schema.fields[16].as_type, "int[blockCount]")
  assert_eq(schema.fields[17].as_type, "lstring")
  assert_eq(schema.fields[18].as_type, "string")
}

///|
test "bigpsl: AutoSQL text contains the UCSC declaration" {
  let text = @src.bigpsl_schema().to_auto_sql()
  assert_true(text.contains("table bigPsl"))
  assert_true(text.contains("int[blockCount] oChromStarts"))
  assert_true(text.contains("uint seqType"))
}

///|
test "bigpsl: constructor preserves metadata and defaults counts" {
  let alignment = @src.BigPslAlignment::create(
    "chr1",
    100,
    "queryA",
    20,
    [10, 15, 20, 25],
    [0, 5, 5, 10],
    sequence_type=@src.BigPslSequenceType::BigPslNucleotide,
    score=700,
    thick_start=12,
    thick_end=24,
    item_rgb="10,20,30",
  ) catch {
    _ => abort("alignment should be valid")
  }
  assert_eq(alignment.target_name, "chr1")
  assert_eq(alignment.query_name, "queryA")
  assert_eq(alignment.score, 700)
  assert_eq(alignment.thick_start, 12)
  assert_eq(alignment.thick_end, 24)
  assert_eq(alignment.item_rgb, "10,20,30")
  assert_eq(alignment.matches, 10)
  assert_eq(alignment.mismatches, 0)
}

///|
test "bigpsl: constructor accepts explicit match categories" {
  let alignment = @src.BigPslAlignment::create(
    "chr1",
    100,
    "query",
    20,
    [10, 20],
    [0, 10],
    matches=6,
    mismatches=2,
    repeat_matches=1,
    n_count=1,
  ) catch {
    _ => abort("match categories should be valid")
  }
  assert_eq(alignment.matches, 6)
  assert_eq(alignment.mismatches, 2)
  assert_eq(alignment.repeat_matches, 1)
  assert_eq(alignment.n_count, 1)
}

///|
test "bigpsl: constructor rejects an empty target name" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create("", 100, "query", 20, [10, 20], [0, 10]),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: constructor rejects a zero target size" {
  let rejected = try {
    ignore(@src.BigPslAlignment::create("chr1", 0, "query", 20, [0, 1], [0, 1]))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: constructor rejects coordinate length mismatch" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create("chr1", 100, "query", 20, [10, 20], [
        0, 5, 10,
      ]),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: constructor rejects a one-point path" {
  let rejected = try {
    ignore(@src.BigPslAlignment::create("chr1", 100, "query", 20, [10], [0]))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: constructor rejects out-of-bounds target coordinates" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create("chr1", 100, "query", 20, [95, 105], [0, 10]),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: constructor rejects out-of-bounds query coordinates" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create("chr1", 100, "query", 20, [10, 20], [15, 25]),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: constructor rejects repeated coordinate points" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create("chr1", 100, "query", 20, [10, 10], [0, 0]),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: constructor rejects coordinate direction changes" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create("chr1", 100, "query", 20, [10, 20, 15], [
        0, 10, 15,
      ]),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: nucleotide alignment rejects reverse target coordinates" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create("chr1", 100, "query", 20, [20, 10], [0, 10]),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: nucleotide alignment requires equal diagonal steps" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create("chr1", 100, "query", 20, [10, 20], [0, 5]),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: protein alignment accepts three-to-one blocks" {
  let alignment = @src.BigPslAlignment::create(
    "chr1",
    100,
    "protein",
    20,
    [10, 40],
    [0, 10],
    sequence_type=@src.BigPslSequenceType::BigPslProtein,
  ) catch {
    _ => abort("translated alignment should be valid")
  }
  assert_eq(alignment.blocks()[0].target_size, 30)
  assert_eq(alignment.blocks()[0].query_size, 10)
}

///|
test "bigpsl: protein alignment rejects non-triplet blocks" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create(
        "chr1",
        100,
        "protein",
        20,
        [10, 39],
        [0, 10],
        sequence_type=@src.BigPslSequenceType::BigPslProtein,
      ),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: protein alignment rejects reverse query coordinates" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create(
        "chr1",
        100,
        "protein",
        20,
        [10, 40],
        [10, 0],
        sequence_type=@src.BigPslSequenceType::BigPslProtein,
      ),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: protein alignment accepts reverse target coordinates" {
  let alignment = @src.BigPslAlignment::create(
    "chr1",
    100,
    "protein",
    20,
    [40, 10],
    [0, 10],
    sequence_type=@src.BigPslSequenceType::BigPslProtein,
  ) catch {
    _ => abort("reverse translated alignment should be valid")
  }
  assert_true(alignment.is_reverse())
}

///|
test "bigpsl: constructor rejects a path without aligned blocks" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create("chr1", 100, "query", 20, [10, 20], [0, 0]),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: constructor rejects score above one thousand" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create(
        "chr1",
        100,
        "query",
        20,
        [10, 20],
        [0, 10],
        score=1001,
      ),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: constructor rejects an external thick interval" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create(
        "chr1",
        100,
        "query",
        20,
        [10, 20],
        [0, 10],
        thick_start=5,
        thick_end=20,
      ),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: constructor rejects partially supplied match categories" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create(
        "chr1",
        100,
        "query",
        20,
        [10, 20],
        [0, 10],
        matches=10,
      ),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: constructor rejects a match total inconsistent with blocks" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create(
        "chr1",
        100,
        "query",
        20,
        [10, 20],
        [0, 10],
        matches=9,
        mismatches=0,
        repeat_matches=0,
        n_count=0,
      ),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: constructor rejects a query sequence length mismatch" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create(
        "chr1",
        100,
        "query",
        20,
        [10, 20],
        [0, 10],
        query_sequence="ACGT",
      ),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: constructor rejects query sequence gaps" {
  let rejected = try {
    ignore(
      @src.BigPslAlignment::create(
        "chr1",
        100,
        "query",
        5,
        [10, 15],
        [0, 5],
        query_sequence="AC-GT",
      ),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: block extraction preserves oriented coordinates" {
  let alignment = bigpsl_test_alignment(target_coordinates=[10, 15, 20, 25], query_coordinates=[
    0, 5, 5, 10,
  ])
  let blocks = alignment.blocks()
  assert_eq(blocks.length(), 2)
  assert_eq(blocks[0].target_start, 10)
  assert_eq(blocks[0].target_end, 15)
  assert_eq(blocks[1].query_start, 5)
  assert_eq(blocks[1].query_end, 10)
}

///|
test "bigpsl: intervals include all blocks and gaps" {
  let alignment = bigpsl_test_alignment(target_coordinates=[10, 15, 20, 25], query_coordinates=[
    2, 7, 7, 12,
  ])
  assert_eq(alignment.target_interval(), (10, 25))
  assert_eq(alignment.query_interval(), (2, 12))
}

///|
test "bigpsl: forward coordinate mapping is bidirectional" {
  let alignment = bigpsl_test_alignment()
  assert_eq(alignment.target_to_query(12), Some(2))
  assert_eq(alignment.query_to_target_interval(2), Some((12, 13)))
  assert_eq(alignment.target_to_query(25), None)
}

///|
test "bigpsl: reverse query mapping uses reverse-complement positions" {
  let alignment = bigpsl_test_alignment(target_coordinates=[10, 20], query_coordinates=[
    15, 5,
  ])
  assert_eq(alignment.target_to_query(10), Some(14))
  assert_eq(alignment.target_to_query(19), Some(5))
  assert_eq(alignment.query_to_target_interval(14), Some((10, 11)))
}

///|
test "bigpsl: protein mapping returns codon intervals" {
  let alignment = @src.BigPslAlignment::create(
    "chr1",
    100,
    "protein",
    10,
    [10, 25],
    [2, 7],
    sequence_type=@src.BigPslSequenceType::BigPslProtein,
  ) catch {
    _ => abort("protein alignment should be valid")
  }
  assert_eq(alignment.target_to_query(10), Some(2))
  assert_eq(alignment.target_to_query(12), Some(2))
  assert_eq(alignment.target_to_query(13), Some(3))
  assert_eq(alignment.query_to_target_interval(3), Some((13, 16)))
}

///|
test "bigpsl: reverse protein mapping returns forward codon intervals" {
  let alignment = @src.BigPslAlignment::create(
    "chr1",
    100,
    "protein",
    10,
    [25, 10],
    [2, 7],
    sequence_type=@src.BigPslSequenceType::BigPslProtein,
  ) catch {
    _ => abort("reverse protein alignment should be valid")
  }
  assert_eq(alignment.target_to_query(24), Some(2))
  assert_eq(alignment.target_to_query(22), Some(2))
  assert_eq(alignment.query_to_target_interval(2), Some((22, 25)))
  assert_eq(alignment.query_to_target_interval(6), Some((10, 13)))
}

///|
test "bigpsl: gap counts distinguish query and target insertions" {
  let alignment = @src.BigPslAlignment::create(
    "chr1",
    100,
    "query",
    30,
    [10, 15, 15, 20, 23, 28],
    [0, 5, 7, 12, 12, 17],
  ) catch {
    _ => abort("gapped alignment should be valid")
  }
  let counts = alignment.counts()
  assert_eq(counts.aligned_query_units, 15)
  assert_eq(counts.query_insert_events, 1)
  assert_eq(counts.query_insert_units, 2)
  assert_eq(counts.target_insert_events, 1)
  assert_eq(counts.target_insert_bases, 3)
  assert_eq(counts.block_count, 3)
}

///|
test "bigpsl: nucleotide recount classifies matches mismatches and wildcard" {
  let alignment = @src.BigPslAlignment::create("chr1", 6, "query", 6, [0, 6], [
    0, 6,
  ]) catch {
    _ => abort("recount alignment should be valid")
  }
  let recounted = alignment.recount("AaCGNT", "AACGAC") catch {
    _ => abort("nucleotide recount should succeed")
  }
  assert_eq(recounted.matches, 4)
  assert_eq(recounted.mismatches, 1)
  assert_eq(recounted.repeat_matches, 0)
  assert_eq(recounted.n_count, 1)
  assert_eq(recounted.query_sequence, "AACGAC")
}

///|
test "bigpsl: lower-case mask reports repeat matches" {
  let alignment = @src.BigPslAlignment::create("chr1", 6, "query", 6, [0, 6], [
    0, 6,
  ]) catch {
    _ => abort("recount alignment should be valid")
  }
  let recounted = alignment.recount(
    "AaCGNT",
    "AACGAT",
    mask=@src.BigPslMaskMode::BigPslMaskLower,
  ) catch {
    _ => abort("lower-mask recount should succeed")
  }
  assert_eq(recounted.matches, 4)
  assert_eq(recounted.repeat_matches, 1)
  assert_eq(recounted.n_count, 1)
}

///|
test "bigpsl: upper-case mask reports repeat matches" {
  let alignment = @src.BigPslAlignment::create("chr1", 6, "query", 6, [0, 6], [
    0, 6,
  ]) catch {
    _ => abort("recount alignment should be valid")
  }
  let recounted = alignment.recount(
    "AaCGNT",
    "AACGAT",
    mask=@src.BigPslMaskMode::BigPslMaskUpper,
  ) catch {
    _ => abort("upper-mask recount should succeed")
  }
  assert_eq(recounted.matches, 1)
  assert_eq(recounted.repeat_matches, 4)
  assert_eq(recounted.n_count, 1)
}

///|
test "bigpsl: reverse nucleotide recount complements the query" {
  let alignment = @src.BigPslAlignment::create("chr1", 4, "query", 4, [0, 4], [
    4, 0,
  ]) catch {
    _ => abort("reverse recount alignment should be valid")
  }
  let recounted = alignment.recount("AACG", "CGTT") catch {
    _ => abort("reverse recount should succeed")
  }
  assert_eq(recounted.matches, 4)
  assert_eq(recounted.mismatches, 0)
}

///|
test "bigpsl: forward translated recount uses the standard genetic code" {
  let alignment = @src.BigPslAlignment::create(
    "chr1",
    9,
    "protein",
    3,
    [0, 9],
    [0, 3],
    sequence_type=@src.BigPslSequenceType::BigPslProtein,
  ) catch {
    _ => abort("protein recount alignment should be valid")
  }
  let recounted = alignment.recount("ATGAAATTT", "MKF") catch {
    _ => abort("protein recount should succeed")
  }
  assert_eq(recounted.matches, 3)
  assert_eq(recounted.mismatches, 0)
}

///|
test "bigpsl: reverse translated recount reverse-complements target codons" {
  let alignment = @src.BigPslAlignment::create(
    "chr1",
    9,
    "protein",
    3,
    [9, 0],
    [0, 3],
    sequence_type=@src.BigPslSequenceType::BigPslProtein,
  ) catch {
    _ => abort("reverse protein recount alignment should be valid")
  }
  let recounted = alignment.recount("AAATTTCAT", "MKF") catch {
    _ => abort("reverse protein recount should succeed")
  }
  assert_eq(recounted.matches, 3)
}

///|
test "bigpsl: translated recount supports X wildcard" {
  let alignment = @src.BigPslAlignment::create(
    "chr1",
    9,
    "protein",
    3,
    [0, 9],
    [0, 3],
    sequence_type=@src.BigPslSequenceType::BigPslProtein,
  ) catch {
    _ => abort("protein recount alignment should be valid")
  }
  let recounted = alignment.recount("ATGNNNTTT", "MXF", wildcard='X') catch {
    _ => abort("ambiguous protein recount should succeed")
  }
  assert_eq(recounted.matches, 2)
  assert_eq(recounted.n_count, 1)
}

///|
test "bigpsl: protein recount rejects repeat masking" {
  let alignment = @src.BigPslAlignment::create(
    "chr1",
    9,
    "protein",
    3,
    [0, 9],
    [0, 3],
    sequence_type=@src.BigPslSequenceType::BigPslProtein,
  ) catch {
    _ => abort("protein recount alignment should be valid")
  }
  let rejected = try {
    ignore(
      alignment.recount(
        "ATGAAATTT",
        "MKF",
        mask=@src.BigPslMaskMode::BigPslMaskLower,
      ),
    )
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: compressed binary round-trip preserves all example paths" {
  let (_, expected) = bigpsl_test_data()
  let observed = bigpsl_test_file().alignments() catch {
    _ => abort("BigPsl alignments should decode")
  }
  assert_eq(observed.length(), expected.length())
  for index in 0..<expected.length() {
    assert_eq(observed[index], expected[index])
  }
}

///|
test "bigpsl: uncompressed binary round-trip is supported" {
  let file = bigpsl_test_file(compress=false)
  assert_true(!file.is_compressed())
  let alignments = file.alignments() catch {
    _ => abort("uncompressed BigPsl should decode")
  }
  assert_eq(alignments.length(), 3)
}

///|
test "bigpsl: multiple data slots preserve target order" {
  let file = bigpsl_test_file(items_per_slot=1)
  assert_true(file.bed.blocks.length() >= 3)
  let alignments = file.alignments() catch {
    _ => abort("multi-block BigPsl should decode")
  }
  assert_eq(alignments[0].target_name, "chr1")
  assert_eq(alignments[2].target_name, "chr2")
}

///|
test "bigpsl: forward nucleotide path survives binary storage" {
  let alignment = (bigpsl_test_file().alignments() catch {
      _ => abort("alignments should decode")
    })[0]
  assert_eq(alignment.target_coordinates, [100, 110, 120, 130])
  assert_eq(alignment.query_coordinates, [0, 10, 10, 20])
}

///|
test "bigpsl: reverse nucleotide path survives binary storage" {
  let alignment = (bigpsl_test_file().alignments() catch {
      _ => abort("alignments should decode")
    })[1]
  assert_eq(alignment.target_coordinates, [200, 208, 208, 215])
  assert_eq(alignment.query_coordinates, [30, 22, 20, 13])
  assert_true(alignment.is_reverse())
}

///|
test "bigpsl: reverse translated path survives binary storage" {
  let alignment = (bigpsl_test_file().alignments() catch {
      _ => abort("alignments should decode")
    })[2]
  assert_eq(alignment.target_coordinates, [500, 488, 470, 470, 455])
  assert_eq(alignment.query_coordinates, [0, 4, 4, 6, 11])
  assert_eq(alignment.sequence_type, @src.BigPslSequenceType::BigPslProtein)
}

///|
test "bigpsl: query sequence and CDS round-trip when requested" {
  let target = @src.BigBedTarget::create("chr1", 100) catch {
    _ => abort("target should be valid")
  }
  let alignment = @src.BigPslAlignment::create(
    "chr1",
    100,
    "transcript",
    10,
    [10, 20],
    [0, 10],
    sequence_type=@src.BigPslSequenceType::BigPslNucleotide,
    query_sequence="AACCGGTTAA",
    cds="2..8",
  ) catch {
    _ => abort("annotated alignment should be valid")
  }
  let config = @src.BigPslWriteConfig::create(
    store_query_sequence=true,
    store_cds=true,
  ) catch {
    _ => abort("storage configuration should be valid")
  }
  let bytes = @src.bigpsl_write([target], [alignment], config~) catch {
    _ => abort("annotated BigPsl should serialize")
  }
  let parsed = @src.bigpsl_parse(bytes) catch {
    _ => abort("annotated BigPsl should parse")
  }
  let observed = (parsed.alignments() catch {
      _ => abort("annotated alignment should decode")
    })[0]
  assert_eq(observed.query_sequence, "AACCGGTTAA")
  assert_eq(observed.cds, "2..8")
}

///|
test "bigpsl: writer omits query sequence and CDS by default" {
  let target = @src.BigBedTarget::create("chr1", 100) catch {
    _ => abort("target should be valid")
  }
  let alignment = @src.BigPslAlignment::create(
    "chr1",
    100,
    "transcript",
    10,
    [10, 20],
    [0, 10],
    query_sequence="AACCGGTTAA",
    cds="2..8",
  ) catch {
    _ => abort("annotated alignment should be valid")
  }
  let bytes = @src.bigpsl_write([target], [alignment]) catch {
    _ => abort("BigPsl should serialize")
  }
  let parsed = @src.bigpsl_parse(bytes) catch {
    _ => abort("BigPsl should parse")
  }
  let observed = (parsed.alignments() catch {
      _ => abort("alignment should decode")
    })[0]
  assert_eq(observed.query_sequence, "")
  assert_eq(observed.cds, "")
}

///|
test "bigpsl: writer requires sequences when sequence storage is enabled" {
  let (targets, alignments) = bigpsl_test_data()
  let config = @src.BigPslWriteConfig::create(store_query_sequence=true) catch {
    _ => abort("storage configuration should be valid")
  }
  let rejected = try {
    ignore(@src.bigpsl_write(targets, alignments, config~))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: writer rejects an unknown target" {
  let target = @src.BigBedTarget::create("chr2", 100) catch {
    _ => abort("target should be valid")
  }
  let rejected = try {
    ignore(@src.bigpsl_write([target], [bigpsl_test_alignment()]))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: writer rejects a target length mismatch" {
  let target = @src.BigBedTarget::create("chr1", 101) catch {
    _ => abort("target should be valid")
  }
  let rejected = try {
    ignore(@src.bigpsl_write([target], [bigpsl_test_alignment()]))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: interval search uses the BigBed R-tree" {
  let file = bigpsl_test_file(items_per_slot=1)
  let first = file.search("chr1", start=105, end=106) catch {
    _ => abort("interval search should succeed")
  }
  assert_eq(first.length(), 1)
  assert_eq(first[0].query_name, "rna-forward")
  let second = file.search("chr1", start=210, end=211) catch {
    _ => abort("interval search should succeed")
  }
  assert_eq(second.length(), 1)
  assert_eq(second[0].query_name, "rna-reverse")
}

///|
test "bigpsl: interval search rejects an unknown target" {
  let rejected = try {
    ignore(bigpsl_test_file().search("chr9", start=1, end=2))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: query-name lookup returns matching alignments" {
  let results = bigpsl_test_file().find_by_query("protein-reverse") catch {
    _ => abort("query lookup should succeed")
  }
  assert_eq(results.length(), 1)
  assert_eq(results[0].target_name, "chr2")
}

///|
test "bigpsl: target lookup exposes chromosome sizes" {
  match bigpsl_test_file().target("chr2") {
    Some(target) => assert_eq(target.length, 800)
    None => abort("chr2 target should exist")
  }
  assert_eq(bigpsl_test_file().target("missing"), None)
}

///|
test "bigpsl: parser rejects a non-BigPsl schema" {
  let target = @src.BigBedTarget::create("chr1", 100) catch {
    _ => abort("target should be valid")
  }
  let record = @src.BigBedRecord::create("chr1", 10, 20) catch {
    _ => abort("record should be valid")
  }
  let config = @src.BigBedWriteConfig::create(bed_columns=3) catch {
    _ => abort("plain BigBed configuration should be valid")
  }
  let bytes = @src.bigbed_write([target], [record], config~) catch {
    _ => abort("plain BigBed should serialize")
  }
  let rejected = try {
    ignore(@src.bigpsl_parse(bytes))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: parser rejects an unknown seqType" {
  let extra = bigpsl_test_standard_extra()
  extra[12] = "9"
  let rejected = try {
    ignore(@src.bigpsl_parse(bigpsl_test_raw_bytes(extra)))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: parser rejects mismatched query block counts" {
  let extra = bigpsl_test_standard_extra()
  extra[4] = "0,5"
  let rejected = try {
    ignore(@src.bigpsl_parse(bigpsl_test_raw_bytes(extra)))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: parser rejects a query sequence length mismatch" {
  let extra = bigpsl_test_standard_extra()
  extra[5] = "ACGT"
  let rejected = try {
    ignore(@src.bigpsl_parse(bigpsl_test_raw_bytes(extra)))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: parser rejects a chromosome size mismatch" {
  let extra = bigpsl_test_standard_extra()
  extra[7] = "99"
  let rejected = try {
    ignore(@src.bigpsl_parse(bigpsl_test_raw_bytes(extra)))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: parser rejects inconsistent match categories" {
  let extra = bigpsl_test_standard_extra()
  extra[8] = "9"
  let rejected = try {
    ignore(@src.bigpsl_parse(bigpsl_test_raw_bytes(extra)))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: parser rejects an inconsistent other interval" {
  let extra = bigpsl_test_standard_extra()
  extra[0] = "1"
  let rejected = try {
    ignore(@src.bigpsl_parse(bigpsl_test_raw_bytes(extra)))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: parser rejects oStrand minus with BED strand plus" {
  let extra = bigpsl_test_standard_extra()
  extra[2] = "-"
  let rejected = try {
    ignore(@src.bigpsl_parse(bigpsl_test_raw_bytes(extra)))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: parser rejects a non-triplet protein block" {
  let extra = bigpsl_test_standard_extra()
  extra[8] = "0"
  extra[9] = "0"
  extra[10] = "0"
  extra[11] = "0"
  extra[12] = "2"
  let rejected = try {
    ignore(@src.bigpsl_parse(bigpsl_test_raw_bytes(extra)))
    false
  } catch {
    BigPslError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigpsl: nucleotide PSL export contains exact block and gap fields" {
  let alignment = bigpsl_test_data().1[0]
  let psl = alignment.to_psl()
  assert_true(psl.has_prefix("20\t0\t0\t0\t0\t0\t1\t10\t+"))
  assert_true(
    psl.contains("\trna-forward\t30\t0\t20\tchr1\t1000\t100\t130\t2\t"),
  )
  assert_true(psl.contains("10,10,\t0,10,\t100,120,"))
}

///|
test "bigpsl: reverse translated PSL export uses two-character strand" {
  let alignment = bigpsl_test_data().1[2]
  let psl = alignment.to_psl()
  assert_true(psl.contains("\t+-\tprotein-reverse\t"))
  assert_true(psl.contains("\t2\t4,5,\t0,6,"))
}

///|
test "bigpsl: file PSL export includes every query" {
  let text = bigpsl_test_file().to_psl() catch {
    _ => abort("PSL export should succeed")
  }
  assert_eq(text.split("\n").length(), 4)
  assert_true(text.contains("rna-forward"))
  assert_true(text.contains("rna-reverse"))
  assert_true(text.contains("protein-reverse"))
}

///|
test "bigpsl: alignment summary identifies type and block count" {
  let summary = bigpsl_test_data().1[2].summary()
  assert_true(summary.contains("chr2:455-500"))
  assert_true(summary.contains("protein-reverse"))
  assert_true(summary.contains("blocks=2"))
  assert_true(summary.contains("type=protein"))
}

///|
test "bigpsl: file summary aggregates sequence types and coverage" {
  let summary = bigpsl_test_file().summary() catch {
    _ => abort("file summary should succeed")
  }
  assert_eq(summary.alignment_count, 3)
  assert_eq(summary.target_count, 2)
  assert_eq(summary.nucleotide_count, 2)
  assert_eq(summary.protein_count, 1)
  assert_eq(summary.unknown_count, 0)
  assert_eq(summary.aligned_query_units, 44)
  assert_eq(summary.covered_target_bases, 90)
  assert_true(summary.compressed)
}
