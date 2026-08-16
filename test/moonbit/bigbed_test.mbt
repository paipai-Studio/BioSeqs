///|
fn bigbed_test_config(
  compress? : Bool = true,
  items_per_slot? : Int = 2,
) -> @src.BigBedWriteConfig {
  @src.BigBedWriteConfig::create(compress~, items_per_slot~, block_size=4) catch {
    _ => abort("BigBed test configuration should be valid")
  }
}

///|
fn bigbed_test_file(
  compress? : Bool = true,
  items_per_slot? : Int = 2,
) -> @src.BigBedFile {
  let (targets, records) = @src.bigbed_example_data()
  let config = bigbed_test_config(compress~, items_per_slot~)
  let bytes = @src.bigbed_write(targets, records, config~) catch {
    _ => abort("BigBed example should serialize")
  }
  @src.bigbed_parse(bytes) catch {
    _ => abort("BigBed example should parse")
  }
}

///|
fn bigbed_simple_record(
  block_count? : Int = 0,
  block_sizes? : Array[Int] = [],
  block_starts? : Array[Int] = [],
  extra_fields? : Array[String] = [],
) -> @src.BigBedRecord {
  @src.BigBedRecord::create(
    "chr1",
    10,
    40,
    name="item",
    score=500,
    strand="+",
    block_count~,
    block_sizes~,
    block_starts~,
    extra_fields~,
  ) catch {
    _ => abort("simple BigBed record should be valid")
  }
}

///|
fn bigbed_test_patch_u32(
  data : Array[Int],
  position : Int,
  value : Int,
) -> Unit {
  data[position] = value & 0xFF
  data[position + 1] = (value >> 8) & 0xFF
  data[position + 2] = (value >> 16) & 0xFF
  data[position + 3] = (value >> 24) & 0xFF
}

///|
fn bigbed_test_patch_u64(
  data : Array[Int],
  position : Int,
  value : Int,
) -> Unit {
  bigbed_test_patch_u32(data, position, value)
  bigbed_test_patch_u32(data, position + 4, 0)
}

///|
fn bigbed_test_replace_single_block(
  source : Array[Int],
  payload : Array[Int],
  uncompressed_size : Int,
) -> Array[Int] {
  let file = @src.bigbed_parse(source) catch {
    _ => abort("source BigBed should parse")
  }
  if file.blocks.length() != 1 {
    abort("fixture BigBed must contain one data block")
  }
  let block = file.blocks[0]
  let output : Array[Int] = []
  for i in 0..<block.data_offset {
    output.push(source[i])
  }
  for byte in payload {
    output.push(byte)
  }
  for i in file.header.full_index_offset..<source.length() {
    output.push(source[i])
  }
  let new_index_offset = file.header.full_index_offset +
    payload.length() -
    block.data_size
  bigbed_test_patch_u64(output, 24, new_index_offset)
  bigbed_test_patch_u32(output, 52, uncompressed_size)
  bigbed_test_patch_u64(output, new_index_offset + 32, new_index_offset)
  bigbed_test_patch_u64(output, new_index_offset + 76, payload.length())
  output
}

///|
test "bigbed: field constructor preserves metadata" {
  let field = @src.BigBedField::create(
    "float",
    "expression",
    comment="Mean expression",
  ) catch {
    _ => abort("field should be valid")
  }
  assert_eq(field.as_type, "float")
  assert_eq(field.name, "expression")
  assert_eq(field.comment, "Mean expression")
}

///|
test "bigbed: field rejects unsupported AutoSQL type" {
  let rejected = try {
    ignore(@src.BigBedField::create("matrix", "value"))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: default schema has twelve standard fields" {
  let schema = @src.BigBedSchema::default() catch {
    _ => abort("default schema should build")
  }
  assert_eq(schema.name, "bed")
  assert_eq(schema.fields.length(), 12)
  assert_eq(schema.fields[0].name, "chrom")
  assert_eq(schema.fields[11].name, "chromStarts")
}

///|
test "bigbed: reduced default schema follows BED column count" {
  let schema = @src.BigBedSchema::default(bed_columns=3) catch {
    _ => abort("BED3 schema should build")
  }
  assert_eq(schema.fields.length(), 3)
  assert_eq(schema.fields[2].name, "chromEnd")
}

///|
test "bigbed: schema rejects fewer than three fields" {
  let field = @src.BigBedField::create("string", "chrom") catch {
    _ => abort("field should build")
  }
  let rejected = try {
    ignore(@src.BigBedSchema::create("bad", "bad", [field]))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: AutoSQL standard schema round trip" {
  let schema = @src.BigBedSchema::default() catch {
    _ => abort("schema should build")
  }
  let parsed = @src.bigbed_parse_auto_sql(schema.to_auto_sql()) catch {
    _ => abort("generated AutoSQL should parse")
  }
  assert_eq(parsed, schema)
}

///|
test "bigbed: AutoSQL custom arrays parse" {
  let schema = @src.bigbed_parse_auto_sql(
    "table custom\n\"Custom BED\"\n(\n" +
    "string chrom; \"Chromosome\"\n" +
    "uint chromStart; \"Start\"\n" +
    "uint chromEnd; \"End\"\n" +
    "float[3] values; \"Values\"\n" +
    ")\n",
  ) catch {
    _ => abort("custom AutoSQL should parse")
  }
  assert_eq(schema.fields.length(), 4)
  assert_eq(schema.fields[3].as_type, "float[3]")
}

///|
test "bigbed: malformed AutoSQL is rejected" {
  let rejected = try {
    ignore(@src.bigbed_parse_auto_sql("record bed\n"))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: target constructor validates length" {
  let target = @src.BigBedTarget::create("chrX", 1000) catch {
    _ => abort("target should build")
  }
  assert_eq(target.name, "chrX")
  assert_eq(target.length, 1000)
  let rejected = try {
    ignore(@src.BigBedTarget::create("chrX", 0))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: record constructor fills BED defaults" {
  let record = @src.BigBedRecord::create("chr1", 10, 40) catch {
    _ => abort("record should build")
  }
  assert_eq(record.name, ".")
  assert_eq(record.score, 0)
  assert_eq(record.strand, ".")
  assert_eq(record.thick_start, 10)
  assert_eq(record.thick_end, 40)
}

///|
test "bigbed: record span and query size use blocks" {
  let record = bigbed_simple_record(block_count=2, block_sizes=[10, 15], block_starts=[
    0, 15,
  ])
  assert_eq(record.span(), 30)
  assert_eq(record.query_size(), 25)
  assert_eq(record.block_count(), 2)
}

///|
test "bigbed: record rejects invalid coordinates" {
  let rejected = try {
    ignore(@src.BigBedRecord::create("chr1", 20, 20))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: record rejects score outside BED range" {
  let rejected = try {
    ignore(@src.BigBedRecord::create("chr1", 0, 10, score=1001))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: record rejects invalid strand" {
  let rejected = try {
    ignore(@src.BigBedRecord::create("chr1", 0, 10, strand="?"))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: record rejects thick interval outside feature" {
  let rejected = try {
    ignore(
      @src.BigBedRecord::create("chr1", 10, 30, thick_start=5, thick_end=20),
    )
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: record rejects mismatched block arrays" {
  let rejected = try {
    ignore(
      @src.BigBedRecord::create("chr1", 0, 20, block_sizes=[10, 5], block_starts=[
        0,
      ]),
    )
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: record rejects a first block after chromStart" {
  let rejected = try {
    ignore(
      @src.BigBedRecord::create("chr1", 0, 20, block_sizes=[10], block_starts=[
        10,
      ]),
    )
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: overlap uses half-open boundaries" {
  let record = bigbed_simple_record()
  assert_true(record.overlaps("chr1", 0, 11))
  assert_true(record.overlaps("chr1", 39, 50))
  assert_true(!record.overlaps("chr1", 0, 10))
  assert_true(!record.overlaps("chr1", 40, 50))
}

///|
test "bigbed: forward BED12 coordinates include intron gaps" {
  let (_, records) = @src.bigbed_example_data()
  let coordinates = records[0].coordinates()
  assert_eq(coordinates[0], [100, 150, 180, 220, 250, 280])
  assert_eq(coordinates[1], [0, 50, 50, 90, 90, 120])
}

///|
test "bigbed: reverse BED12 coordinates reverse query direction" {
  let (_, records) = @src.bigbed_example_data()
  let coordinates = records[1].coordinates()
  assert_eq(coordinates[0], [250, 310, 360, 410])
  assert_eq(coordinates[1], [110, 50, 50, 0])
}

///|
test "bigbed: BED11 partial blocks use span coordinates" {
  let record = bigbed_simple_record(block_count=2, block_sizes=[10, 20])
  assert_eq(record.block_count(), 2)
  assert_eq(record.coordinates(), [[10, 40], [0, 30]])
}

///|
test "bigbed: default write configuration matches Biopython defaults" {
  let config = @src.BigBedWriteConfig::default()
  assert_eq(config.bed_columns, 12)
  assert_eq(config.items_per_slot, 512)
  assert_eq(config.block_size, 256)
  assert_true(config.compress)
}

///|
test "bigbed: configuration accepts custom fields" {
  let field = @src.BigBedField::create("float", "signal") catch {
    _ => abort("custom field should build")
  }
  let config = @src.BigBedWriteConfig::create(bed_columns=6, custom_fields=[
    field,
  ]) catch {
    _ => abort("custom config should build")
  }
  assert_eq(config.bed_columns, 6)
  assert_eq(config.custom_fields[0].name, "signal")
}

///|
test "bigbed: configuration rejects invalid BED columns" {
  let low = try {
    ignore(@src.BigBedWriteConfig::create(bed_columns=2))
    false
  } catch {
    BigBedError(_) => true
  }
  let high = try {
    ignore(@src.BigBedWriteConfig::create(bed_columns=13))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(low)
  assert_true(high)
}

///|
test "bigbed: configuration rejects invalid tree controls" {
  let slot = try {
    ignore(@src.BigBedWriteConfig::create(items_per_slot=0))
    false
  } catch {
    BigBedError(_) => true
  }
  let block = try {
    ignore(@src.BigBedWriteConfig::create(block_size=0))
    false
  } catch {
    BigBedError(_) => true
  }
  let non_branching = try {
    ignore(@src.BigBedWriteConfig::create(block_size=1))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(slot)
  assert_true(block)
  assert_true(non_branching)
}

///|
test "bigbed: compressed version 4 round trip" {
  let file = bigbed_test_file()
  assert_eq(file.header.version, 4)
  assert_true(file.header.little_endian)
  assert_true(file.is_compressed())
  assert_eq(file.record_count, 3)
}

///|
test "bigbed: compressed round trip preserves all BED12 records" {
  let file = bigbed_test_file()
  let records = file.records() catch { _ => abort("records should decode") }
  let (_, expected) = @src.bigbed_example_data()
  assert_eq(records, expected)
}

///|
test "bigbed: itemsPerSlot creates multiple indexed data blocks" {
  let file = bigbed_test_file(items_per_slot=1)
  assert_eq(file.blocks.length(), 3)
  assert_eq(file.summary().data_block_count, 3)
}

///|
test "bigbed: chromosome B+ tree preserves targets" {
  let file = bigbed_test_file()
  assert_eq(file.targets.length(), 2)
  assert_eq(file.targets[0].name, "chr1")
  assert_eq(file.targets[0].length, 1000000)
  assert_eq(file.targets[1].name, "chr2")
}

///|
test "bigbed: target lookup handles present and absent names" {
  let file = bigbed_test_file()
  match file.target("chr2") {
    Some(target) => assert_eq(target.length, 500000)
    None => abort("chr2 should exist")
  }
  assert_eq(file.target("chrM"), None)
}

///|
test "bigbed: search returns both overlapping chr1 transcripts" {
  let hits = bigbed_test_file().search("chr1", start=260, end=270) catch {
    _ => abort("region search should succeed")
  }
  assert_eq(hits.length(), 2)
  assert_eq(hits[0].name, "transcript-A")
  assert_eq(hits[1].name, "transcript-B")
}

///|
test "bigbed: search excludes record ending at query start" {
  let hits = bigbed_test_file().search("chr1", start=280, end=300) catch {
    _ => abort("boundary search should succeed")
  }
  assert_eq(hits.length(), 1)
  assert_eq(hits[0].name, "transcript-B")
}

///|
test "bigbed: search excludes record starting at query end" {
  let hits = bigbed_test_file().search("chr1", start=0, end=100) catch {
    _ => abort("boundary search should succeed")
  }
  assert_eq(hits.length(), 0)
}

///|
test "bigbed: search prunes other chromosomes" {
  let hits = bigbed_test_file(items_per_slot=1).search(
    "chr2",
    start=1000,
    end=1001,
  ) catch {
    _ => abort("chr2 search should succeed")
  }
  assert_eq(hits.length(), 1)
  assert_eq(hits[0].name, "transcript-C")
}

///|
test "bigbed: search defaults to the complete chromosome" {
  let hits = bigbed_test_file().search("chr1") catch {
    _ => abort("whole chromosome search should succeed")
  }
  assert_eq(hits.length(), 2)
}

///|
test "bigbed: search rejects unknown chromosome" {
  let rejected = try {
    ignore(bigbed_test_file().search("chrM"))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: search rejects invalid interval" {
  let rejected = try {
    ignore(bigbed_test_file().search("chr1", start=20, end=20))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: name lookup scans indexed records" {
  let hits = bigbed_test_file().find_by_name("transcript-C") catch {
    _ => abort("name lookup should succeed")
  }
  assert_eq(hits.length(), 1)
  assert_eq(hits[0].chrom, "chr2")
}

///|
test "bigbed: summary reports spans and compression" {
  let summary = bigbed_test_file().summary() catch {
    _ => abort("summary should succeed")
  }
  assert_eq(summary.record_count, 3)
  assert_eq(summary.target_count, 2)
  assert_eq(summary.covered_bases, 460)
  assert_eq(summary.minimum_span, 120)
  assert_eq(summary.maximum_span, 180)
  assert_true(summary.compressed)
}

///|
test "bigbed: BED text export preserves trailing block commas" {
  let text = bigbed_test_file().to_bed() catch {
    _ => abort("BED export should succeed")
  }
  assert_true(
    text.contains(
      "chr1\t100\t280\ttranscript-A\t960\t+\t120\t260\t40,120,220\t3\t50,40,30,\t0,80,150,\n",
    ),
  )
}

///|
test "bigbed: uncompressed round trip" {
  let file = bigbed_test_file(compress=false)
  assert_true(!file.is_compressed())
  let records = file.records() catch {
    _ => abort("uncompressed records should decode")
  }
  assert_eq(records.length(), 3)
}

///|
test "bigbed: parser decodes standard fixed Huffman zlib blocks" {
  let (targets, records) = @src.bigbed_example_data()
  let config = @src.BigBedWriteConfig::create(
    items_per_slot=3,
    block_size=4,
    compress=false,
  ) catch {
    _ => abort("fixture configuration should build")
  }
  let source = @src.bigbed_write(targets, records, config~) catch {
    _ => abort("fixture source should serialize")
  }
  let fixed = [
    120, 1, 99, 96, 96, 96, 72, 1, 98, 9, 70, 6, 134, 146, 162, 196, 188, 226, 228,
    162, 204, 130, 18, 93, 71, 78, 75, 51, 3, 78, 109, 78, 67, 35, 3, 78, 35, 32,
    203, 196, 64, 7, 200, 212, 49, 2, 114, 141, 57, 77, 13, 116, 128, 124, 99, 3,
    29, 78, 3, 29, 11, 160, 4, 144, 207, 0, 2, 191, 128, 120, 22, 170, 57, 78, 156,
    22, 230, 6, 156, 186, 156, 70, 64, 210, 216, 18, 104, 150, 17, 88, 139, 5, 144,
    197, 105, 102, 160, 99, 10, 54, 195, 208, 16, 104, 0, 80, 31, 195, 11, 102, 6,
    134, 4, 22, 20, 3, 156, 57, 205, 13, 192, 14, 49, 0, 82, 134, 32, 231, 128, 180,
    154, 152, 66, 181, 154, 3, 117, 2, 0, 27, 137, 36, 189,
  ]
  let bytes = bigbed_test_replace_single_block(source, fixed, 195)
  let file = @src.bigbed_parse(bytes) catch {
    _ => abort("fixed-Huffman BigBed should parse")
  }
  let decoded = file.records() catch {
    _ => abort("fixed-Huffman data should decode")
  }
  assert_eq(decoded, records)
}

///|
test "bigbed: parser decodes standard dynamic Huffman zlib blocks" {
  let (targets, records) = @src.bigbed_example_data()
  let config = @src.BigBedWriteConfig::create(
    items_per_slot=3,
    block_size=4,
    compress=false,
  ) catch {
    _ => abort("fixture configuration should build")
  }
  let source = @src.bigbed_write(targets, records, config~) catch {
    _ => abort("fixture source should serialize")
  }
  let dynamic = [
    120, 218, 85, 142, 49, 14, 194, 48, 16, 4, 15, 65, 179, 63, 224, 15, 92, 164,
    243, 57, 142, 157, 18, 248, 8, 8, 26, 26, 132, 66, 94, 147, 23, 242, 6, 42, 214,
    22, 77, 44, 157, 189, 43, 107, 70, 39, 34, 114, 231, 236, 55, 34, 243, 116, 125,
    190, 111, 211, 227, 53, 119, 71, 140, 131, 225, 128, 224, 6, 103, 234, 77, 25,
    213, 89, 35, 146, 41, 123, 52, 133, 105, 225, 7, 187, 212, 243, 229, 44, 107,
    207, 9, 37, 27, 58, 56, 239, 56, 210, 229, 13, 41, 76, 24, 76, 83, 115, 132,
    64, 1, 57, 249, 108, 69, 46, 187, 149, 224, 140, 108, 109, 17, 227, 19, 234,
    58, 21, 237, 211, 31, 205, 36, 127, 27, 137, 36, 189,
  ]
  let bytes = bigbed_test_replace_single_block(source, dynamic, 195)
  let file = @src.bigbed_parse(bytes) catch {
    _ => abort("dynamic-Huffman BigBed should parse")
  }
  let decoded = file.records() catch {
    _ => abort("dynamic-Huffman data should decode")
  }
  assert_eq(decoded, records)
}

///|
test "bigbed: multilevel chromosome and interval trees round trip" {
  let targets : Array[@src.BigBedTarget] = []
  let records : Array[@src.BigBedRecord] = []
  for i in 0..<10 {
    let name = "chr" + (9 - i).to_string()
    targets.push(
      @src.BigBedTarget::create(name, 1000) catch {
        _ => abort("multilevel target should build")
      },
    )
    records.push(
      @src.BigBedRecord::create(name, 10 + i, 30 + i) catch {
        _ => abort("multilevel record should build")
      },
    )
  }
  let config = @src.BigBedWriteConfig::create(
    bed_columns=3,
    items_per_slot=1,
    block_size=2,
    compress=false,
  ) catch {
    _ => abort("multilevel configuration should build")
  }
  let bytes = @src.bigbed_write(targets, records, config~) catch {
    _ => abort("multilevel BigBed should serialize")
  }
  let file = @src.bigbed_parse(bytes) catch {
    _ => abort("multilevel BigBed should parse")
  }
  assert_eq(file.targets, targets)
  assert_eq(file.blocks.length(), 10)
  assert_eq(file.records(), records)
  assert_eq(file.search("chr5", start=0, end=100).length(), 1)
}

///|
test "bigbed: R-tree bounds retain earlier long intervals" {
  let targets = [
    @src.BigBedTarget::create("chr1", 1000) catch {
      _ => abort("target should build")
    },
  ]
  let records = [
    @src.BigBedRecord::create("chr1", 0, 900) catch {
      _ => abort("long record should build")
    },
    @src.BigBedRecord::create("chr1", 100, 200) catch {
      _ => abort("short record should build")
    },
  ]
  let config = @src.BigBedWriteConfig::create(
    bed_columns=3,
    items_per_slot=2,
    block_size=2,
    compress=false,
  ) catch {
    _ => abort("configuration should build")
  }
  let bytes = @src.bigbed_write(targets, records, config~) catch {
    _ => abort("BigBed should serialize")
  }
  let file = @src.bigbed_parse(bytes) catch {
    _ => abort("BigBed should parse")
  }
  let hits = file.search("chr1", start=800, end=850) catch {
    _ => abort("long interval query should succeed")
  }
  assert_eq(hits.length(), 1)
  assert_eq(hits[0].chrom_end, 900)
}

///|
test "bigbed: BED3 round trip drops optional columns" {
  let targets = [
    @src.BigBedTarget::create("chr1", 100) catch {
      _ => abort("target should build")
    },
  ]
  let records = [bigbed_simple_record()]
  let config = @src.BigBedWriteConfig::create(bed_columns=3) catch {
    _ => abort("BED3 config should build")
  }
  let bytes = @src.bigbed_write(targets, records, config~) catch {
    _ => abort("BED3 should serialize")
  }
  let parsed = @src.bigbed_parse(bytes) catch {
    _ => abort("BED3 should parse")
  }
  let record = parsed.records()[0]
  assert_eq(record.name, ".")
  assert_eq(record.score, 0)
  assert_eq(parsed.header.defined_field_count, 3)
}

///|
test "bigbed: BED6 round trip preserves name score and strand" {
  let targets = [
    @src.BigBedTarget::create("chr1", 100) catch {
      _ => abort("target should build")
    },
  ]
  let record = @src.BigBedRecord::create(
    "chr1",
    10,
    40,
    name="feature",
    score=700,
    strand="-",
  ) catch {
    _ => abort("BED6 record should build")
  }
  let config = @src.BigBedWriteConfig::create(bed_columns=6) catch {
    _ => abort("BED6 config should build")
  }
  let parsed = @src.bigbed_parse(@src.bigbed_write(targets, [record], config~)) catch {
    _ => abort("BED6 round trip should succeed")
  }
  let result = parsed.records()[0]
  assert_eq(result.name, "feature")
  assert_eq(result.score, 700)
  assert_eq(result.strand, "-")
}

///|
test "bigbed: BED10 preserves declared block count" {
  let targets = [
    @src.BigBedTarget::create("chr1", 100) catch {
      _ => abort("target should build")
    },
  ]
  let config = @src.BigBedWriteConfig::create(bed_columns=10) catch {
    _ => abort("BED10 config should build")
  }
  let record = bigbed_simple_record(block_count=2)
  let parsed = @src.bigbed_parse(@src.bigbed_write(targets, [record], config~)) catch {
    _ => abort("BED10 round trip should succeed")
  }
  assert_eq(parsed.records()[0].block_count(), 2)
}

///|
test "bigbed: BED11 preserves block sizes without starts" {
  let targets = [
    @src.BigBedTarget::create("chr1", 100) catch {
      _ => abort("target should build")
    },
  ]
  let config = @src.BigBedWriteConfig::create(bed_columns=11) catch {
    _ => abort("BED11 config should build")
  }
  let record = bigbed_simple_record(block_count=2, block_sizes=[10, 20])
  let parsed = @src.bigbed_parse(@src.bigbed_write(targets, [record], config~)) catch {
    _ => abort("BED11 round trip should succeed")
  }
  let result = parsed.records()[0]
  assert_eq(result.block_count(), 2)
  assert_eq(result.block_sizes, [10, 20])
  assert_eq(result.block_starts, [])
}

///|
test "bigbed: custom AutoSQL annotation round trip" {
  let field = @src.BigBedField::create("float", "signal") catch {
    _ => abort("field should build")
  }
  let config = @src.BigBedWriteConfig::create(bed_columns=6, custom_fields=[
    field,
  ]) catch {
    _ => abort("config should build")
  }
  let targets = [
    @src.BigBedTarget::create("chr1", 100) catch {
      _ => abort("target should build")
    },
  ]
  let record = @src.BigBedRecord::create(
    "chr1",
    10,
    40,
    name="feature",
    score=700,
    strand="+",
    extra_fields=["2.75"],
  ) catch {
    _ => abort("custom record should build")
  }
  let parsed = @src.bigbed_parse(@src.bigbed_write(targets, [record], config~)) catch {
    _ => abort("custom BigBed should round trip")
  }
  let result = parsed.records()[0]
  assert_eq(result.annotation(parsed.schema, "signal", 6), Some("2.75"))
}

///|
test "bigbed: writer rejects unsorted records" {
  let (targets, records) = @src.bigbed_example_data()
  let rejected = try {
    ignore(@src.bigbed_write(targets, [records[1], records[0], records[2]]))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: writer rejects unknown target" {
  let targets = [
    @src.BigBedTarget::create("chr1", 100) catch {
      _ => abort("target should build")
    },
  ]
  let record = @src.BigBedRecord::create("chr2", 10, 20) catch {
    _ => abort("record should build")
  }
  let rejected = try {
    ignore(@src.bigbed_write(targets, [record]))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: writer rejects interval beyond target" {
  let targets = [
    @src.BigBedTarget::create("chr1", 20) catch {
      _ => abort("target should build")
    },
  ]
  let record = @src.BigBedRecord::create("chr1", 10, 30) catch {
    _ => abort("record should build")
  }
  let rejected = try {
    ignore(@src.bigbed_write(targets, [record]))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: writer rejects duplicate target names" {
  let targets = [
    @src.BigBedTarget::create("chr1", 100) catch {
      _ => abort("target should build")
    },
    @src.BigBedTarget::create("chr1", 200) catch {
      _ => abort("target should build")
    },
  ]
  let rejected = try {
    ignore(@src.bigbed_write(targets, []))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: writer rejects custom field count mismatch" {
  let field = @src.BigBedField::create("string", "category") catch {
    _ => abort("field should build")
  }
  let config = @src.BigBedWriteConfig::create(bed_columns=3, custom_fields=[
    field,
  ]) catch {
    _ => abort("config should build")
  }
  let targets = [
    @src.BigBedTarget::create("chr1", 100) catch {
      _ => abort("target should build")
    },
  ]
  let rejected = try {
    ignore(@src.bigbed_write(targets, [bigbed_simple_record()], config~))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: deterministic writer produces identical bytes" {
  let (targets, records) = @src.bigbed_example_data()
  let config = bigbed_test_config()
  let first = @src.bigbed_write(targets, records, config~) catch {
    _ => abort("first serialization should succeed")
  }
  let second = @src.bigbed_write(targets, records, config~) catch {
    _ => abort("second serialization should succeed")
  }
  assert_eq(first, second)
}

///|
test "bigbed: parser rejects invalid leading magic" {
  let (targets, records) = @src.bigbed_example_data()
  let bytes = @src.bigbed_write(targets, records) catch {
    _ => abort("serialization should succeed")
  }
  bytes[0] = 0
  let rejected = try {
    ignore(@src.bigbed_parse(bytes))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: parser rejects invalid trailing magic" {
  let (targets, records) = @src.bigbed_example_data()
  let bytes = @src.bigbed_write(targets, records) catch {
    _ => abort("serialization should succeed")
  }
  bytes[bytes.length() - 1] = 0
  let rejected = try {
    ignore(@src.bigbed_parse(bytes))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: parser rejects unsupported version" {
  let (targets, records) = @src.bigbed_example_data()
  let bytes = @src.bigbed_write(targets, records) catch {
    _ => abort("serialization should succeed")
  }
  bytes[4] = 3
  let rejected = try {
    ignore(@src.bigbed_parse(bytes))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: parser rejects truncated header" {
  let rejected = try {
    ignore(@src.bigbed_parse([0xEB, 0xF2, 0x89, 0x87]))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: parser rejects damaged R-tree magic" {
  let (targets, records) = @src.bigbed_example_data()
  let bytes = @src.bigbed_write(targets, records) catch {
    _ => abort("serialization should succeed")
  }
  let file = @src.bigbed_parse(bytes) catch {
    _ => abort("valid file should parse")
  }
  bytes[file.header.full_index_offset] = 0
  let rejected = try {
    ignore(@src.bigbed_parse(bytes))
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: decoder rejects damaged zlib checksum" {
  let (targets, records) = @src.bigbed_example_data()
  let config = bigbed_test_config()
  let bytes = @src.bigbed_write(targets, records, config~) catch {
    _ => abort("serialization should succeed")
  }
  let file = @src.bigbed_parse(bytes) catch {
    _ => abort("valid file should parse")
  }
  let block = file.blocks[0]
  let checksum_byte = block.data_offset + block.data_size - 1
  bytes[checksum_byte] = bytes[checksum_byte] ^ 1
  let damaged = @src.bigbed_parse(bytes) catch {
    _ => abort("index should still parse before data decoding")
  }
  let rejected = try {
    ignore(damaged.records())
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: decoder rejects malformed Huffman DEFLATE" {
  let (targets, records) = @src.bigbed_example_data()
  let config = bigbed_test_config()
  let bytes = @src.bigbed_write(targets, records, config~) catch {
    _ => abort("serialization should succeed")
  }
  let file = @src.bigbed_parse(bytes) catch {
    _ => abort("valid file should parse")
  }
  bytes[file.blocks[0].data_offset + 2] = 3
  let changed = @src.bigbed_parse(bytes) catch {
    _ => abort("index should still parse")
  }
  let rejected = try {
    ignore(changed.records())
    false
  } catch {
    BigBedError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigbed: empty file retains targets and zero records" {
  let targets = [
    @src.BigBedTarget::create("chr1", 100) catch {
      _ => abort("target should build")
    },
  ]
  let bytes = @src.bigbed_write(targets, []) catch {
    _ => abort("empty BigBed should serialize")
  }
  let file = @src.bigbed_parse(bytes) catch {
    _ => abort("empty BigBed should parse")
  }
  assert_eq(file.record_count, 0)
  assert_eq(file.records(), [])
  assert_eq(file.search("chr1"), [])
}
