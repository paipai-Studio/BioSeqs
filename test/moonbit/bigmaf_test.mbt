///|
fn bigmaf_test_data() -> (
  String,
  Array[@src.BigBedTarget],
  Array[@src.BigMafBlock],
) {
  @src.bigmaf_example_data() catch {
    _ => abort("BigMaf example data should be valid")
  }
}

///|
fn bigmaf_test_file(
  compress? : Bool = true,
  items_per_slot? : Int = 2,
) -> @src.BigMafFile {
  let (reference, targets, blocks) = bigmaf_test_data()
  let config = @src.BigMafWriteConfig::create(
    compress~,
    block_size=3,
    items_per_slot~,
  ) catch {
    _ => abort("BigMaf test configuration should be valid")
  }
  let bytes = @src.bigmaf_write(reference, targets, blocks, config~) catch {
    _ => abort("BigMaf example should serialize")
  }
  @src.bigmaf_parse(bytes) catch {
    _ => abort("BigMaf example should parse")
  }
}

///|
fn bigmaf_test_component(
  source? : String = "hg38.chr1",
  start? : Int = 10,
  size? : Int = 4,
  strand? : String = "+",
  source_size? : Int = 100,
  text? : String = "AC-GT",
) -> @src.BigMafComponent {
  @src.BigMafComponent::create(source, start, size, strand, source_size, text) catch {
    _ => abort("test component should be valid")
  }
}

///|
fn bigmaf_test_block(
  component? : @src.BigMafComponent = bigmaf_test_component(),
) -> @src.BigMafBlock {
  @src.BigMafBlock::create([component]) catch {
    _ => abort("test block should be valid")
  }
}

///|
fn bigmaf_test_close(
  actual : Double,
  expected : Double,
  tolerance? : Double = 1.0e-12,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "BigMaf value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn bigmaf_test_custom_bytes(
  payloads : Array[String],
  intervals : Array[(Int, Int)],
) -> Array[Int] {
  let target = @src.BigBedTarget::create("chr1", 100) catch {
    _ => abort("custom target should be valid")
  }
  let field = @src.BigBedField::create("lstring", "mafBlock") catch {
    _ => abort("custom field should be valid")
  }
  let config = @src.BigBedWriteConfig::create(
    bed_columns=3,
    schema_name="bedMaf",
    schema_comment="Bed3 with MAF block",
    custom_fields=[field],
  ) catch {
    _ => abort("custom config should be valid")
  }
  let records : Array[@src.BigBedRecord] = []
  for index in 0..<payloads.length() {
    records.push(
      @src.BigBedRecord::create(
        "chr1",
        intervals[index].0,
        intervals[index].1,
        extra_fields=[payloads[index]],
      ) catch {
        _ => abort("custom record should be valid")
      },
    )
  }
  @src.bigbed_write([target], records, config~) catch {
    _ => abort("custom BigBed should serialize")
  }
}

///|
test "bigmaf: default configuration matches Biopython defaults" {
  let config = @src.BigMafWriteConfig::default()
  assert_true(config.compress)
  assert_eq(config.block_size, 256)
  assert_eq(config.items_per_slot, 512)
}

///|
test "bigmaf: custom configuration is preserved" {
  let config = @src.BigMafWriteConfig::create(
    compress=false,
    block_size=8,
    items_per_slot=4,
  ) catch {
    _ => abort("configuration should be valid")
  }
  assert_true(!config.compress)
  assert_eq(config.block_size, 8)
  assert_eq(config.items_per_slot, 4)
}

///|
test "bigmaf: configuration rejects blockSize below two" {
  let rejected = try {
    ignore(@src.BigMafWriteConfig::create(block_size=1))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: configuration rejects zero itemsPerSlot" {
  let rejected = try {
    ignore(@src.BigMafWriteConfig::create(items_per_slot=0))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: standard AutoSQL declaration is bed3 plus mafBlock" {
  let schema = @src.bigmaf_schema()
  assert_eq(schema.name, "bedMaf")
  assert_eq(schema.comment, "Bed3 with MAF block")
  assert_eq(schema.fields.length(), 4)
  assert_eq(schema.fields[3].as_type, "lstring")
  assert_eq(schema.fields[3].name, "mafBlock")
}

///|
test "bigmaf: AutoSQL text contains the UCSC table declaration" {
  let text = @src.bigmaf_schema().to_auto_sql()
  assert_true(text.contains("table bedMaf"))
  assert_true(text.contains("lstring mafBlock"))
}

///|
test "bigmaf: insertion constructor preserves context" {
  let insertion = @src.BigMafInsertion::create("C", 0, "I", 12) catch {
    _ => abort("insertion context should be valid")
  }
  assert_eq(insertion.left_status, "C")
  assert_eq(insertion.left_count, 0)
  assert_eq(insertion.right_status, "I")
  assert_eq(insertion.right_count, 12)
}

///|
test "bigmaf: insertion rejects unknown status" {
  let rejected = try {
    ignore(@src.BigMafInsertion::create("X", 0, "C", 0))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: insertion rejects negative count" {
  let rejected = try {
    ignore(@src.BigMafInsertion::create("C", -1, "C", 0))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: component constructor preserves MAF fields" {
  let component = bigmaf_test_component()
  assert_eq(component.source, "hg38.chr1")
  assert_eq(component.start, 10)
  assert_eq(component.size, 4)
  assert_eq(component.strand, "+")
  assert_eq(component.source_size, 100)
  assert_eq(component.text, "AC-GT")
}

///|
test "bigmaf: component sequence removes alignment gaps" {
  assert_eq(bigmaf_test_component().sequence(), "ACGT")
}

///|
test "bigmaf: plus component reports forward interval" {
  assert_eq(bigmaf_test_component().forward_interval(), (10, 14))
}

///|
test "bigmaf: minus component converts MAF reverse coordinates" {
  let component = bigmaf_test_component(
    source="mm39.chr1",
    start=10,
    strand="-",
  )
  assert_eq(component.forward_interval(), (86, 90))
}

///|
test "bigmaf: minus component returns forward-oriented sequence" {
  let component = bigmaf_test_component(
    source="mm39.chr1",
    start=10,
    size=3,
    strand="-",
    text="ACG",
  )
  assert_eq(component.forward_sequence(), "CGT")
}

///|
test "bigmaf: plus column maps to source coordinate" {
  let component = bigmaf_test_component()
  assert_eq(component.column_to_source(0), Some(10))
  assert_eq(component.column_to_source(4), Some(13))
}

///|
test "bigmaf: alignment gap maps to no source coordinate" {
  assert_eq(bigmaf_test_component().column_to_source(2), None)
}

///|
test "bigmaf: plus source coordinate maps to alignment column" {
  let component = bigmaf_test_component()
  assert_eq(component.source_to_column(10), Some(0))
  assert_eq(component.source_to_column(12), Some(3))
  assert_eq(component.source_to_column(14), None)
}

///|
test "bigmaf: minus columns descend in forward coordinates" {
  let component = bigmaf_test_component(
    source="mm39.chr1",
    start=10,
    strand="-",
  )
  assert_eq(component.column_to_source(0), Some(89))
  assert_eq(component.column_to_source(4), Some(86))
}

///|
test "bigmaf: minus source coordinate maps back to column" {
  let component = bigmaf_test_component(
    source="mm39.chr1",
    start=10,
    strand="-",
  )
  assert_eq(component.source_to_column(89), Some(0))
  assert_eq(component.source_to_column(87), Some(3))
}

///|
test "bigmaf: component accepts aligned quality string" {
  let component = @src.BigMafComponent::create(
    "hg38.chr1",
    10,
    4,
    "+",
    100,
    "AC-GT",
    quality=Some("9F-73"),
  ) catch {
    _ => abort("quality should be valid")
  }
  assert_eq(component.quality, Some("9F-73"))
}

///|
test "bigmaf: component rejects quality width mismatch" {
  let rejected = try {
    ignore(
      @src.BigMafComponent::create(
        "hg38.chr1",
        10,
        4,
        "+",
        100,
        "AC-GT",
        quality=Some("9999"),
      ),
    )
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: component rejects quality gap mismatch" {
  let rejected = try {
    ignore(
      @src.BigMafComponent::create(
        "hg38.chr1",
        10,
        4,
        "+",
        100,
        "AC-GT",
        quality=Some("99F99"),
      ),
    )
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: component rejects non-gap size mismatch" {
  let rejected = try {
    ignore(@src.BigMafComponent::create("hg38.chr1", 10, 5, "+", 100, "AC-GT"))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: component rejects coordinate overflow" {
  let rejected = try {
    ignore(@src.BigMafComponent::create("hg38.chr1", 98, 4, "+", 100, "ACGT"))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: component rejects invalid strand" {
  let rejected = try {
    ignore(@src.BigMafComponent::create("hg38.chr1", 10, 4, ".", 100, "ACGT"))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: empty component reports plus interval" {
  let empty = @src.BigMafEmptyComponent::create(
    "dog.chr1", 20, 5, "+", 100, "M",
  ) catch {
    _ => abort("empty component should be valid")
  }
  assert_eq(empty.forward_interval(), (20, 25))
}

///|
test "bigmaf: empty component reports minus interval" {
  let empty = @src.BigMafEmptyComponent::create(
    "dog.chr1", 20, 5, "-", 100, "n",
  ) catch {
    _ => abort("empty component should be valid")
  }
  assert_eq(empty.forward_interval(), (75, 80))
}

///|
test "bigmaf: empty component rejects invalid status" {
  let rejected = try {
    ignore(@src.BigMafEmptyComponent::create("dog.chr1", 20, 5, "+", 100, "T"))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: block preserves score pass comments and empty components" {
  let empty = @src.BigMafEmptyComponent::create(
    "dog.chr1", 20, 5, "+", 100, "M",
  ) catch {
    _ => abort("empty component should be valid")
  }
  let block = @src.BigMafBlock::create(
    [bigmaf_test_component()],
    score=Some(42.5),
    pass_number=Some(2),
    empty_components=[empty],
    comments=["test"],
  ) catch {
    _ => abort("block should be valid")
  }
  assert_eq(block.score, Some(42.5))
  assert_eq(block.pass_number, Some(2))
  assert_eq(block.empty_components.length(), 1)
  assert_eq(block.comments, ["test"])
}

///|
test "bigmaf: block rejects unequal alignment widths" {
  let rejected = try {
    ignore(
      @src.BigMafBlock::create([
        bigmaf_test_component(),
        bigmaf_test_component(source="mouse.chr1", size=3, text="ACG"),
      ]),
    )
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: block rejects non-positive pass annotation" {
  let rejected = try {
    ignore(
      @src.BigMafBlock::create([bigmaf_test_component()], pass_number=Some(0)),
    )
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: block rejects empty component array" {
  let rejected = try {
    ignore(@src.BigMafBlock::create([]))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: embedded formatter uses semicolon-delimited MAF" {
  let (_, _, blocks) = bigmaf_test_data()
  let text = blocks[0].to_embedded_maf()
  assert_true(text.contains("a score=23262"))
  assert_true(text.contains(";s hg38.chr7 100 8 + 1000 ACGT-ACGT;"))
  assert_true(text.contains(";i mm39.chr5 C 0 I 12;"))
  assert_true(text.contains(";q hg38.chr7 9999-FFFF;"))
  assert_true(text.contains(";e canFam6.chr6 300 20 + 850 M;"))
}

///|
test "bigmaf: embedded block round trip preserves all annotations" {
  let (_, _, blocks) = bigmaf_test_data()
  let parsed = @src.bigmaf_parse_block(blocks[0].to_embedded_maf()) catch {
    _ => abort("embedded block should parse")
  }
  assert_eq(parsed, blocks[0])
}

///|
test "bigmaf: parser accepts comments before alignment line" {
  let block = @src.bigmaf_parse_block(
    "# first comment;#second;a score=1;s hg38.chr1 0 4 + 100 ACGT;",
  ) catch {
    _ => abort("commented block should parse")
  }
  assert_eq(block.comments, ["first comment", "second"])
}

///|
test "bigmaf: parser reads score and pass" {
  let block = @src.bigmaf_parse_block(
    "a score=1.25 pass=3;s hg38.chr1 0 4 + 100 ACGT;",
  ) catch {
    _ => abort("annotated block should parse")
  }
  assert_eq(block.score, Some(1.25))
  assert_eq(block.pass_number, Some(3))
}

///|
test "bigmaf: parser reads insertion quality and empty lines" {
  let block = @src.bigmaf_parse_block(
    "a;s hg38.chr1 0 4 + 100 AC-GT;i hg38.chr1 C 0 I 2;q hg38.chr1 99-FF;e dog.chr1 10 5 + 100 M;",
  ) catch {
    _ => abort("full MAF block should parse")
  }
  assert_true(block.components[0].insertion is Some(_))
  assert_eq(block.components[0].quality, Some("99-FF"))
  assert_eq(block.empty_components.length(), 1)
}

///|
test "bigmaf: parser rejects unknown alignment annotation" {
  let rejected = try {
    ignore(
      @src.bigmaf_parse_block("a score=1 custom=x;s hg38.chr1 0 4 + 100 ACGT;"),
    )
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: parser rejects malformed s line" {
  let rejected = try {
    ignore(@src.bigmaf_parse_block("a;s hg38.chr1 0 4 + 100;"))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: parser rejects insertion before component" {
  let rejected = try {
    ignore(
      @src.bigmaf_parse_block(
        "a;i hg38.chr1 C 0 C 0;s hg38.chr1 0 4 + 100 ACGT;",
      ),
    )
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: parser rejects mismatched insertion source" {
  let rejected = try {
    ignore(
      @src.bigmaf_parse_block(
        "a;s hg38.chr1 0 4 + 100 ACGT;i mouse.chr1 C 0 C 0;",
      ),
    )
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: parser rejects duplicate quality lines" {
  let rejected = try {
    ignore(
      @src.bigmaf_parse_block(
        "a;s hg38.chr1 0 4 + 100 ACGT;q hg38.chr1 9999;q hg38.chr1 9999;",
      ),
    )
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: parser rejects unknown line type" {
  let rejected = try {
    ignore(@src.bigmaf_parse_block("a;s hg38.chr1 0 4 + 100 ACGT;x unknown;"))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: parser rejects multiple alignment lines" {
  let rejected = try {
    ignore(
      @src.bigmaf_parse_block(
        "a;s hg38.chr1 0 4 + 100 ACGT;a;s hg38.chr1 4 4 + 100 ACGT;",
      ),
    )
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: parser rejects payload without alignment line" {
  let rejected = try {
    ignore(@src.bigmaf_parse_block("# comment;s hg38.chr1 0 4 + 100 ACGT;"))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: reference mapping follows shared alignment column" {
  let (_, _, blocks) = bigmaf_test_data()
  assert_eq(blocks[0].map_reference_to(100, "panTro6.chr6"), Some(200))
  assert_eq(blocks[0].map_reference_to(104, "panTro6.chr6"), Some(204))
}

///|
test "bigmaf: reference mapping returns none for target gap" {
  let (_, _, blocks) = bigmaf_test_data()
  assert_eq(blocks[1].map_reference_to(302, "mm39.chr5"), None)
}

///|
test "bigmaf: reference mapping handles reverse target component" {
  let (_, _, blocks) = bigmaf_test_data()
  assert_eq(blocks[0].map_reference_to(100, "mm39.chr5"), Some(649))
  assert_eq(blocks[0].map_reference_to(102, "mm39.chr5"), None)
}

///|
test "bigmaf: pairwise identity counts non-gap comparisons" {
  let block = @src.BigMafBlock::create([
    bigmaf_test_component(source="ref.chr1", size=4, text="AC-GT"),
    bigmaf_test_component(source="a.chr1", size=4, text="AC-AT"),
  ]) catch {
    _ => abort("identity block should be valid")
  }
  bigmaf_test_close(block.pairwise_identity(), 0.75)
}

///|
test "bigmaf: component lookup and source list preserve order" {
  let (_, _, blocks) = bigmaf_test_data()
  assert_eq(blocks[0].sources(), ["hg38.chr7", "panTro6.chr6", "mm39.chr5"])
  match blocks[0].component("mm39.chr5") {
    Some(component) => assert_eq(component.strand, "-")
    None => abort("mouse component should exist")
  }
  assert_eq(blocks[0].component("missing"), None)
}

///|
test "bigmaf: compressed binary round trip preserves blocks" {
  let file = bigmaf_test_file()
  let (_, _, expected) = bigmaf_test_data()
  assert_eq(file.blocks(), expected)
  assert_true(file.is_compressed())
}

///|
test "bigmaf: uncompressed binary round trip preserves blocks" {
  let file = bigmaf_test_file(compress=false)
  let (_, _, expected) = bigmaf_test_data()
  assert_eq(file.blocks(), expected)
  assert_true(!file.is_compressed())
}

///|
test "bigmaf: parsed binary exposes standard schema" {
  let file = bigmaf_test_file()
  assert_eq(file.bed.header.defined_field_count, 3)
  assert_eq(file.bed.header.field_count, 4)
  assert_eq(file.bed.schema, @src.bigmaf_schema())
}

///|
test "bigmaf: parser recovers reference prefix" {
  assert_eq(bigmaf_test_file().reference, "hg38")
}

///|
test "bigmaf: parser preserves chromosome targets" {
  let file = bigmaf_test_file()
  assert_eq(file.bed.targets.length(), 2)
  assert_eq(file.bed.targets[0].name, "chr7")
  assert_eq(file.bed.targets[1].name, "chr12")
}

///|
test "bigmaf: writer sorts blocks by target and start" {
  let (reference, targets, blocks) = bigmaf_test_data()
  let bytes = @src.bigmaf_write(reference, targets, [
    blocks[1],
    blocks[2],
    blocks[0],
  ]) catch {
    _ => abort("unsorted blocks should serialize")
  }
  let parsed = @src.bigmaf_parse(bytes) catch {
    _ => abort("sorted BigMaf should parse")
  }
  let ordered = parsed.blocks()
  assert_eq(ordered[0], blocks[0])
  assert_eq(ordered[1], blocks[1])
  assert_eq(ordered[2], blocks[2])
}

///|
test "bigmaf: interval search returns overlapping block" {
  let hits = bigmaf_test_file().search("chr7", start=104, end=105) catch {
    _ => abort("interval search should succeed")
  }
  assert_eq(hits.length(), 1)
  assert_eq(hits[0].score, Some(23262.0))
}

///|
test "bigmaf: interval search accepts reference-prefixed chromosome" {
  let hits = bigmaf_test_file().search("hg38.chr7", start=300, end=301) catch {
    _ => abort("prefixed interval search should succeed")
  }
  assert_eq(hits.length(), 1)
  assert_eq(hits[0].score, Some(5062.0))
}

///|
test "bigmaf: interval search uses half-open boundaries" {
  let hits = bigmaf_test_file().search("chr7", start=108, end=109) catch {
    _ => abort("boundary search should succeed")
  }
  assert_eq(hits.length(), 0)
}

///|
test "bigmaf: interval search rejects unknown chromosome" {
  let rejected = try {
    ignore(bigmaf_test_file().search("chrM"))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: all blocks retain scores" {
  let blocks = bigmaf_test_file().blocks() catch {
    _ => abort("all blocks should decode")
  }
  assert_eq(blocks.length(), 3)
  assert_eq(blocks[0].score, Some(23262.0))
  assert_eq(blocks[2].score, Some(6636.0))
}

///|
test "bigmaf: summary reports alignment contents" {
  let summary = bigmaf_test_file().summary() catch {
    _ => abort("summary should succeed")
  }
  assert_eq(summary.alignment_count, 3)
  assert_eq(summary.target_count, 2)
  assert_eq(summary.component_count, 8)
  assert_eq(summary.empty_component_count, 1)
  assert_eq(summary.aligned_columns, 23)
  assert_eq(summary.covered_reference_bases, 21)
  assert_true(summary.compressed)
}

///|
test "bigmaf: MAF export restores newline-delimited blocks" {
  let text = bigmaf_test_file().to_maf() catch {
    _ => abort("MAF export should succeed")
  }
  assert_true(text.has_prefix("##maf version=1\n\n"))
  assert_true(text.contains("a score=23262"))
  assert_true(text.contains("s hg38.chr12 40 7 + 800 GCA-GCTG"))
  assert_true(!text.contains(";s "))
}

///|
test "bigmaf: writer is deterministic" {
  let (reference, targets, blocks) = bigmaf_test_data()
  let first = @src.bigmaf_write(reference, targets, blocks) catch {
    _ => abort("first write should succeed")
  }
  let second = @src.bigmaf_write(reference, targets, blocks) catch {
    _ => abort("second write should succeed")
  }
  assert_eq(first, second)
}

///|
test "bigmaf: itemsPerSlot creates multiple indexed blocks" {
  let file = bigmaf_test_file(items_per_slot=1)
  assert_eq(file.bed.blocks.length(), 3)
  assert_eq(file.search("chr12").length(), 1)
}

///|
test "bigmaf: writer rejects inconsistent reference prefix" {
  let (_, targets, blocks) = bigmaf_test_data()
  let rejected = try {
    ignore(@src.bigmaf_write("hg19", targets, blocks))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: writer rejects dotted reference name" {
  let (_, targets, blocks) = bigmaf_test_data()
  let rejected = try {
    ignore(@src.bigmaf_write("hg38.patch", targets, blocks))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: writer rejects reference source size mismatch" {
  let target = @src.BigBedTarget::create("chr1", 101) catch {
    _ => abort("target should be valid")
  }
  let rejected = try {
    ignore(@src.bigmaf_write("hg38", [target], [bigmaf_test_block()]))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: writer rejects reverse reference component" {
  let target = @src.BigBedTarget::create("chr1", 100) catch {
    _ => abort("target should be valid")
  }
  let component = bigmaf_test_component(strand="-")
  let rejected = try {
    ignore(@src.bigmaf_write("hg38", [target], [bigmaf_test_block(component~)]))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: writer rejects empty targets" {
  let rejected = try {
    ignore(@src.bigmaf_write("hg38", [], []))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: parser rejects ordinary BigBed schema" {
  let target = @src.BigBedTarget::create("chr1", 100) catch {
    _ => abort("target should be valid")
  }
  let record = @src.BigBedRecord::create("chr1", 0, 4) catch {
    _ => abort("record should be valid")
  }
  let config = @src.BigBedWriteConfig::create(bed_columns=3) catch {
    _ => abort("config should be valid")
  }
  let bytes = @src.bigbed_write([target], [record], config~) catch {
    _ => abort("BigBed should serialize")
  }
  let rejected = try {
    ignore(@src.bigmaf_parse(bytes))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: parser rejects BED and MAF coordinate disagreement" {
  let bytes = bigmaf_test_custom_bytes(["a;s hg38.chr1 1 4 + 100 ACGT;"], [
    (0, 4),
  ])
  let rejected = try {
    ignore(@src.bigmaf_parse(bytes))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: parser rejects inconsistent reference prefixes" {
  let bytes = bigmaf_test_custom_bytes(
    ["a;s hg38.chr1 0 4 + 100 ACGT;", "a;s hg19.chr1 10 4 + 100 ACGT;"],
    [(0, 4), (10, 14)],
  )
  let rejected = try {
    ignore(@src.bigmaf_parse(bytes))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: parser rejects malformed embedded block" {
  let bytes = bigmaf_test_custom_bytes(["not-maf;"], [(0, 4)])
  let rejected = try {
    ignore(@src.bigmaf_parse(bytes))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: empty binary file retains targets" {
  let target = @src.BigBedTarget::create("chr1", 100) catch {
    _ => abort("target should be valid")
  }
  let bytes = @src.bigmaf_write("hg38", [target], []) catch {
    _ => abort("empty BigMaf should serialize")
  }
  let file = @src.bigmaf_parse(bytes) catch {
    _ => abort("empty BigMaf should parse")
  }
  assert_eq(file.reference, "")
  assert_eq(file.blocks(), [])
  assert_eq(file.target("chr1"), Some(target))
}

///|
test "bigmaf: target lookup accepts reference-prefixed name" {
  let file = bigmaf_test_file()
  assert_eq(file.target("hg38.chr7"), file.target("chr7"))
  assert_eq(file.target("hg38."), None)
}

///|
test "bigmaf: source coordinate outside component returns none" {
  assert_eq(bigmaf_test_component().source_to_column(9), None)
  assert_eq(bigmaf_test_component().source_to_column(14), None)
}

///|
test "bigmaf: column mapping rejects out-of-range column" {
  let rejected = try {
    ignore(bigmaf_test_component().column_to_source(5))
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}

///|
test "bigmaf: block rejects semicolon in comment" {
  let rejected = try {
    ignore(
      @src.BigMafBlock::create([bigmaf_test_component()], comments=[
        "invalid;comment",
      ]),
    )
    false
  } catch {
    BigMafError(_) => true
  }
  assert_true(rejected)
}
