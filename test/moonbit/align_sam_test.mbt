///|
fn align_sam_test_parse(text : String) -> @src.AlignSamDocument {
  @src.align_sam_parse(text) catch {
    AlignSamError(message) => abort("valid SAM failed: " + message)
  }
}

///|
fn align_sam_test_parse_raises(text : String) -> Bool {
  try {
    ignore(@src.align_sam_parse(text))
    false
  } catch {
    AlignSamError(_) => true
  }
}

///|
fn align_sam_test_cigar_raises(text : String) -> Bool {
  try {
    ignore(@src.align_sam_parse_cigar(text))
    false
  } catch {
    AlignSamError(_) => true
  }
}

///|
fn align_sam_test_record(
  flag : String,
  reference : String,
  position : String,
  mapq : String,
  cigar : String,
  mate_reference : String,
  mate_position : String,
  template_length : String,
  sequence : String,
  quality : String,
  tags : String,
) -> String {
  let suffix = if tags == "" { "" } else { "\t" + tags }
  "read\t" +
  flag +
  "\t" +
  reference +
  "\t" +
  position +
  "\t" +
  mapq +
  "\t" +
  cigar +
  "\t" +
  mate_reference +
  "\t" +
  mate_position +
  "\t" +
  template_length +
  "\t" +
  sequence +
  "\t" +
  quality +
  suffix +
  "\n"
}

///|
fn align_sam_test_mapped(
  cigar : String,
  sequence : String,
  quality : String,
  tags : String,
) -> String {
  "@SQ\tSN:chr1\tLN:100\n" +
  align_sam_test_record(
    "0", "chr1", "1", "60", cigar, "*", "0", "0", sequence, quality, tags,
  )
}

///|
fn align_sam_test_forward() -> @src.AlignSamAlignment {
  align_sam_test_parse(@src.align_sam_example_data()).alignments[0]
}

///|
fn align_sam_test_reverse() -> @src.AlignSamAlignment {
  align_sam_test_parse(@src.align_sam_example_data()).alignments[1]
}

///|
fn align_sam_test_reference() -> String {
  "N".repeat(4) + "ACGTTCGGAA" + "N".repeat(26)
}

///|
test "Bio.Align.sam parses header, references, and alignments" {
  let document = align_sam_test_parse(@src.align_sam_example_data())
  assert_eq(document.headers.length(), 3)
  assert_eq(document.references.length(), 1)
  assert_eq(document.alignments.length(), 3)
}

///|
test "Bio.Align.sam preserves HD metadata" {
  let document = align_sam_test_parse(@src.align_sam_example_data())
  let header = document.headers[0]
  assert_eq(header.record_type, "HD")
  assert_eq(header.fields[0].key, "VN")
  assert_eq(header.fields[0].value, "1.6")
  assert_eq(header.fields[1].key, "SO")
  assert_eq(header.fields[1].value, "coordinate")
}

///|
test "Bio.Align.sam exposes SQ metadata and lookup" {
  let document = align_sam_test_parse(@src.align_sam_example_data())
  match document.reference("chr1") {
    Some(reference) => {
      assert_eq(reference.name, "chr1")
      assert_eq(reference.length, 40)
      assert_eq(reference.annotations.length(), 1)
      assert_eq(reference.annotations[0].key, "AS")
      assert_eq(reference.annotations[0].value, "demo")
    }
    None => abort("chr1 reference was not indexed")
  }
  assert_true(document.reference("missing") is None)
}

///|
test "Bio.Align.sam parses comments and extension headers" {
  let document = align_sam_test_parse(
    "@HD\tVN:1.6\n@CO\tfirst\tsecond\n@ZZ\tID:extension\n",
  )
  assert_eq(document.headers.length(), 3)
  assert_eq(document.headers[1].record_type, "CO")
  assert_eq(document.headers[1].comment, "first\tsecond")
  assert_eq(document.headers[2].record_type, "ZZ")
  assert_eq(document.headers[2].fields[0].value, "extension")
}

///|
test "Bio.Align.sam accepts CRLF and blank lines" {
  let document = align_sam_test_parse(
    "\r\n@SQ\tSN:chr1\tLN:10\r\n\r\n" +
    "read\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\r\n",
  )
  assert_eq(document.alignments.length(), 1)
  assert_eq(document.alignments[0].query_name, "read")
}

///|
test "Bio.Align.sam parses an empty document" {
  let document = align_sam_test_parse("\n\r\n")
  assert_eq(document.headers.length(), 0)
  assert_eq(document.references.length(), 0)
  assert_eq(document.alignments.length(), 0)
}

///|
test "Bio.Align.sam forward record fields" {
  let alignment = align_sam_test_forward()
  assert_eq(alignment.query_name, "read_forward")
  assert_eq(alignment.reference_name, "chr1")
  assert_eq(alignment.reference_length, Some(40))
  assert_eq(alignment.reference_start(), Some(4))
  assert_eq(alignment.reference_end(), Some(14))
  assert_eq(alignment.mapq, Some(60))
  assert_false(alignment.is_reverse())
}

///|
test "Bio.Align.sam forward coordinate path" {
  let alignment = align_sam_test_forward()
  assert_eq(alignment.target_coordinates, [4, 8, 8, 10, 11, 14])
  assert_eq(alignment.query_coordinates, [3, 7, 8, 10, 10, 13])
  assert_eq(alignment.operations, [
    @src.SamAligned,
    @src.SamInsertion,
    @src.SamAligned,
    @src.SamDeletion,
    @src.SamAligned,
  ])
}

///|
test "Bio.Align.sam reverse sequence and qualities use biological orientation" {
  let alignment = align_sam_test_reverse()
  assert_true(alignment.is_reverse())
  assert_eq(alignment.query_sequence, "TTACGTACGTA")
  assert_eq(alignment.qualities.length(), 11)
  assert_eq(alignment.qualities[0], 41)
  assert_eq(alignment.qualities[10], 41)
  assert_eq(alignment.hard_clip_left, Some(1))
  assert_eq(alignment.hard_clip_right, Some(2))
}

///|
test "Bio.Align.sam reverse coordinate path decreases query coordinates" {
  let alignment = align_sam_test_reverse()
  assert_eq(alignment.target_coordinates, [19, 24, 27, 30])
  assert_eq(alignment.query_coordinates, [9, 4, 4, 1])
  assert_eq(alignment.operations, [
    @src.SamAligned,
    @src.SamSkipped,
    @src.SamAligned,
  ])
}

///|
test "Bio.Align.sam unmapped record has no alignment path" {
  let alignment = align_sam_test_parse(@src.align_sam_example_data()).alignments[2]
  assert_true(alignment.is_unmapped())
  assert_eq(alignment.reference_name, "*")
  assert_eq(alignment.target_coordinates.length(), 0)
  assert_eq(alignment.query_coordinates.length(), 0)
  assert_eq(alignment.operations.length(), 0)
  assert_eq(alignment.query_sequence, "ACGT")
}

///|
test "Bio.Align.sam flag helpers" {
  let document = align_sam_test_parse(
    align_sam_test_record(
      (0x1 | 0x4 | 0x100 | 0x800).to_string(),
      "*",
      "0",
      "0",
      "*",
      "*",
      "0",
      "0",
      "*",
      "*",
      "",
    ),
  )
  let alignment = document.alignments[0]
  assert_true(alignment.is_paired())
  assert_true(alignment.is_unmapped())
  assert_true(alignment.is_secondary())
  assert_true(alignment.is_supplementary())
  assert_false(alignment.is_reverse())
}

///|
test "Bio.Align.sam converts mate positions to zero based" {
  let alignment = align_sam_test_parse(
      "@SQ\tSN:chr1\tLN:20\n@SQ\tSN:chr2\tLN:20\n" +
      align_sam_test_record(
        "1", "chr1", "2", "30", "4M", "chr2", "8", "-5", "ACGT", "!!!!", "",
      ),
    ).alignments[0]
  assert_eq(alignment.mate_reference, Some("chr2"))
  assert_eq(alignment.mate_position, Some(7))
  assert_eq(alignment.template_length, Some(-5))
}

///|
test "Bio.Align.sam expands equal mate reference" {
  let alignment = align_sam_test_parse(
      "@SQ\tSN:chr1\tLN:20\n" +
      align_sam_test_record(
        "1", "chr1", "2", "30", "4M", "=", "8", "5", "ACGT", "!!!!", "",
      ),
    ).alignments[0]
  assert_eq(alignment.mate_reference, Some("chr1"))
  let formatted = @src.align_sam_format(alignment) catch {
    AlignSamError(message) => abort(message)
  }
  assert_eq(
    formatted,
    align_sam_test_record(
      "1", "chr1", "2", "30", "4M", "=", "8", "5", "ACGT", "!!!!", "",
    ),
  )
}

///|
test "Bio.Align.sam MAPQ 255 is represented as unknown" {
  let alignment = align_sam_test_parse(
      align_sam_test_record(
        "0", "contig", "1", "255", "4M", "*", "0", "0", "ACGT", "*", "",
      ),
    ).alignments[0]
  assert_true(alignment.mapq is None)
  assert_true(alignment.reference_length is None)
}

///|
test "Bio.Align.sam parses reference-less mapped SAM" {
  let alignment = align_sam_test_parse(
      align_sam_test_record(
        "0", "contig", "3", "20", "2M1D2M", "*", "0", "0", "ACGT", "!!!!", "",
      ),
    ).alignments[0]
  assert_eq(alignment.reference_start(), Some(2))
  assert_eq(alignment.reference_end(), Some(7))
  assert_true(alignment.reference_length is None)
}

///|
test "Bio.Align.sam CIGAR supports all alignment operations" {
  let cigar = @src.align_sam_parse_cigar("1H2S3M4I5D6N7=8X9S10H") catch {
    AlignSamError(message) => abort(message)
  }
  assert_eq(cigar.length(), 10)
  assert_eq(cigar[0].operation, "H")
  assert_eq(cigar[1].operation, "S")
  assert_eq(cigar[2].operation, "M")
  assert_eq(cigar[3].operation, "I")
  assert_eq(cigar[4].operation, "D")
  assert_eq(cigar[5].operation, "N")
  assert_eq(cigar[6].operation, "=")
  assert_eq(cigar[7].operation, "X")
  assert_eq(cigar[8].length, 9)
  assert_eq(cigar[9].length, 10)
}

///|
test "Bio.Align.sam CIGAR star has no elements" {
  let cigar = @src.align_sam_parse_cigar("*") catch {
    AlignSamError(message) => abort(message)
  }
  assert_eq(cigar.length(), 0)
  assert_eq(@src.align_sam_cigar_string(cigar), "*")
}

///|
test "Bio.Align.sam CIGAR round trip" {
  let cigar = @src.align_sam_parse_cigar("2H3S4M1I2=1X1D5N2M1S4H") catch {
    AlignSamError(message) => abort(message)
  }
  assert_eq(@src.align_sam_cigar_string(cigar), "2H3S4M1I2=1X1D5N2M1S4H")
}

///|
test "Bio.Align.sam reports CIGAR statistics" {
  let stats = align_sam_test_forward().stats()
  assert_eq(stats.reference_consumed, 10)
  assert_eq(stats.query_consumed, 15)
  assert_eq(stats.aligned, 9)
  assert_eq(stats.inserted, 1)
  assert_eq(stats.deleted, 1)
  assert_eq(stats.skipped, 0)
  assert_eq(stats.soft_clipped, 5)
  assert_eq(stats.hard_clipped, 0)
}

///|
test "Bio.Align.sam reports spliced and hard clipping statistics" {
  let stats = align_sam_test_reverse().stats()
  assert_eq(stats.reference_consumed, 11)
  assert_eq(stats.query_consumed, 11)
  assert_eq(stats.aligned, 8)
  assert_eq(stats.skipped, 3)
  assert_eq(stats.soft_clipped, 3)
  assert_eq(stats.hard_clipped, 3)
}

///|
test "Bio.Align.sam maps target bases to forward query bases" {
  let alignment = align_sam_test_forward()
  assert_eq(alignment.target_to_query(4), Some(3))
  assert_eq(alignment.target_to_query(7), Some(6))
  assert_eq(alignment.target_to_query(8), Some(8))
  assert_eq(alignment.target_to_query(10), None)
  assert_eq(alignment.target_to_query(13), Some(12))
}

///|
test "Bio.Align.sam maps forward query bases to target bases" {
  let alignment = align_sam_test_forward()
  assert_eq(alignment.query_to_target(3), Some(4))
  assert_eq(alignment.query_to_target(7), None)
  assert_eq(alignment.query_to_target(8), Some(8))
  assert_eq(alignment.query_to_target(12), Some(13))
  assert_eq(alignment.query_to_target(14), None)
}

///|
test "Bio.Align.sam maps target bases to reverse query bases" {
  let alignment = align_sam_test_reverse()
  assert_eq(alignment.target_to_query(19), Some(8))
  assert_eq(alignment.target_to_query(23), Some(4))
  assert_eq(alignment.target_to_query(24), None)
  assert_eq(alignment.target_to_query(27), Some(3))
  assert_eq(alignment.target_to_query(29), Some(1))
}

///|
test "Bio.Align.sam maps reverse query bases to target bases" {
  let alignment = align_sam_test_reverse()
  assert_eq(alignment.query_to_target(8), Some(19))
  assert_eq(alignment.query_to_target(4), Some(23))
  assert_eq(alignment.query_to_target(3), Some(27))
  assert_eq(alignment.query_to_target(1), Some(29))
  assert_eq(alignment.query_to_target(0), None)
}

///|
test "Bio.Align.sam parses scalar typed tags" {
  let alignment = align_sam_test_parse(
      align_sam_test_mapped(
        "4M", "ACGT", "!!!!", "AA:A:x\tXI:i:-2147483648\tFF:f:1.5\tZZ:Z:hello\tHX:H:0aFF",
      ),
    ).alignments[0]
  assert_eq(alignment.string_tag("AA"), Some("x"))
  assert_eq(alignment.integer_tag("XI"), Some(-2147483647 - 1))
  assert_eq(alignment.string_tag("ZZ"), Some("hello"))
  assert_eq(alignment.string_tag("HX"), Some("0aFF"))
  match alignment.tag("FF") {
    Some(tag) =>
      match tag.value {
        @src.SamTagFloat(value) => assert_eq(value, 1.5)
        _ => abort("FF did not retain float type")
      }
    None => abort("FF tag is missing")
  }
}

///|
test "Bio.Align.sam parses typed integer arrays" {
  let alignment = align_sam_test_parse(
      align_sam_test_mapped(
        "4M", "ACGT", "!!!!", "BC:B:c,-1,2\tBU:B:C,0,255\tBS:B:s,-2,3\tBT:B:S,0,65535\tBI:B:i,-4,5\tBJ:B:I,0,7",
      ),
    ).alignments[0]
  match alignment.tag("BC") {
    Some(tag) =>
      match tag.value {
        @src.SamTagIntegerArray(kind, values) => {
          assert_eq(kind, "c")
          assert_eq(values, [-1, 2])
        }
        _ => abort("BC did not retain integer array type")
      }
    None => abort("BC tag is missing")
  }
  match alignment.tag("BT") {
    Some(tag) =>
      match tag.value {
        @src.SamTagIntegerArray(kind, values) => {
          assert_eq(kind, "S")
          assert_eq(values, [0, 65535])
        }
        _ => abort("BT did not retain integer array type")
      }
    None => abort("BT tag is missing")
  }
  assert_eq(alignment.tags.length(), 6)
}

///|
test "Bio.Align.sam parses typed float arrays" {
  let alignment = align_sam_test_parse(
      align_sam_test_mapped("4M", "ACGT", "!!!!", "BF:B:f,1.5,-2.25"),
    ).alignments[0]
  match alignment.tag("BF") {
    Some(tag) =>
      match tag.value {
        @src.SamTagFloatArray(values) => assert_eq(values, [1.5, -2.25])
        _ => abort("BF did not retain float array type")
      }
    None => abort("BF tag is missing")
  }
}

///|
test "Bio.Align.sam tag lookup returns defensive array copies" {
  let alignment = align_sam_test_parse(
      align_sam_test_mapped("4M", "ACGT", "!!!!", "BI:B:i,1,2"),
    ).alignments[0]
  match alignment.tag("BI") {
    Some(tag) =>
      match tag.value {
        @src.SamTagIntegerArray(_, values) => values[0] = 99
        _ => abort("BI did not retain integer array type")
      }
    None => abort("BI tag is missing")
  }
  match alignment.tag("BI") {
    Some(tag) =>
      match tag.value {
        @src.SamTagIntegerArray(_, values) => assert_eq(values, [1, 2])
        _ => abort("BI did not retain integer array type")
      }
    None => abort("BI tag is missing")
  }
}

///|
test "Bio.Align.sam typed tag constructors" {
  let tags = [
    @src.align_sam_integer_tag("NM", 2),
    @src.align_sam_float_tag("ZF", 1.25),
    @src.align_sam_string_tag("RG", "group1"),
    @src.align_sam_character_tag("XC", "Y"),
    @src.align_sam_hex_tag("XH", "CAFE"),
    @src.align_sam_integer_array_tag("BI", "s", [-2, 3]),
    @src.align_sam_float_array_tag("BF", [0.5, 1.5]),
  ] catch {
    AlignSamError(message) => abort(message)
  }
  let alignment = @src.align_sam_create(
    "created",
    "chr1",
    20,
    2,
    "ACGT",
    "4M",
    tags~,
  ) catch {
    AlignSamError(message) => abort(message)
  }
  assert_eq(alignment.tags.length(), 7)
  assert_eq(alignment.integer_tag("NM"), Some(2))
  assert_eq(alignment.string_tag("RG"), Some("group1"))
}

///|
test "Bio.Align.sam reconstructs aligned rows from MD" {
  let rows = align_sam_test_forward().aligned_rows() catch {
    AlignSamError(message) => abort(message)
  }
  assert_eq(rows.0, "ACGT-TCGGAA")
  assert_eq(rows.1, "ACGTATC-GTA")
}

///|
test "Bio.Align.sam reconstructs skipped rows without a reference" {
  let rows = align_sam_test_reverse().aligned_rows() catch {
    AlignSamError(message) => abort(message)
  }
  assert_eq(rows.0, "CGTAC???GTA")
  assert_eq(rows.1, "CGTAC---GTA")
}

///|
test "Bio.Align.sam reconstructs rows from a full reference" {
  let rows = align_sam_test_forward().aligned_rows(
    reference_sequence=align_sam_test_reference(),
  ) catch {
    AlignSamError(message) => abort(message)
  }
  assert_eq(rows.0, "ACGT-TCGGAA")
  assert_eq(rows.1, "ACGTATC-GTA")
}

///|
test "Bio.Align.sam computes NM from MD, insertions, and deletions" {
  assert_eq(align_sam_test_forward().computed_nm(), Some(3))
  assert_eq(align_sam_test_reverse().computed_nm(), Some(0))
}

///|
test "Bio.Align.sam returns unknown NM for ambiguous M without MD" {
  let alignment = align_sam_test_parse(
      align_sam_test_mapped("4M", "ACGT", "!!!!", ""),
    ).alignments[0]
  assert_true(alignment.computed_nm() is None)
}

///|
test "Bio.Align.sam computes NM from explicit equal and mismatch CIGAR" {
  let alignment = align_sam_test_parse(
      align_sam_test_mapped("2=1X1I1D1=", "ACGTA", "!!!!!", ""),
    ).alignments[0]
  assert_eq(alignment.computed_nm(), Some(3))
}

///|
test "Bio.Align.sam calculates an MD tag from reference sequence" {
  let md = @src.align_sam_calculate_md(
    align_sam_test_forward(),
    align_sam_test_reference(),
  ) catch {
    AlignSamError(message) => abort(message)
  }
  assert_eq(md, "6^G1A1")
}

///|
test "Bio.Align.sam calculated MD replaces an existing tag on format" {
  let formatted = @src.align_sam_format(
    align_sam_test_forward(),
    calculate_md=true,
    reference_sequence=align_sam_test_reference(),
  ) catch {
    AlignSamError(message) => abort(message)
  }
  let alignment = align_sam_test_parse("@SQ\tSN:chr1\tLN:40\n" + formatted).alignments[0]
  assert_eq(alignment.string_tag("MD"), Some("6^G1A1"))
  assert_eq(alignment.tags.length(), 4)
}

///|
test "Bio.Align.sam record formatting round trips reverse orientation" {
  let original = align_sam_test_reverse()
  let formatted = @src.align_sam_format(original) catch {
    AlignSamError(message) => abort(message)
  }
  let parsed = align_sam_test_parse("@SQ\tSN:chr1\tLN:40\n" + formatted).alignments[0]
  assert_eq(parsed.query_sequence, original.query_sequence)
  assert_eq(parsed.qualities, original.qualities)
  assert_eq(parsed.query_coordinates, original.query_coordinates)
  assert_eq(parsed.hard_clip_left, original.hard_clip_left)
  assert_eq(parsed.hard_clip_right, original.hard_clip_right)
}

///|
test "Bio.Align.sam document write is canonical and lossless" {
  let text = @src.align_sam_example_data()
  let written = @src.align_sam_write(align_sam_test_parse(text)) catch {
    AlignSamError(message) => abort(message)
  }
  assert_eq(written, text)
}

///|
test "Bio.Align.sam creates forward mapped alignments" {
  let alignment = @src.align_sam_create(
    "created",
    "chr1",
    30,
    5,
    "AACGTT",
    "1S4M1S",
    mapq=42,
    qualities=[10, 11, 12, 13, 14, 15],
  ) catch {
    AlignSamError(message) => abort(message)
  }
  assert_eq(alignment.reference_start(), Some(5))
  assert_eq(alignment.reference_end(), Some(9))
  assert_eq(alignment.query_coordinates, [1, 5])
  assert_eq(alignment.query_sequence, "AACGTT")
  assert_eq(alignment.qualities, [10, 11, 12, 13, 14, 15])
  assert_eq(alignment.mapq, Some(42))
}

///|
test "Bio.Align.sam creates reverse mapped alignments" {
  let alignment = @src.align_sam_create(
    "reverse",
    "chr1",
    30,
    5,
    "AACGTT",
    "1H1S4M1S2H",
    reverse=true,
    qualities=[10, 11, 12, 13, 14, 15],
  ) catch {
    AlignSamError(message) => abort(message)
  }
  assert_true(alignment.is_reverse())
  assert_eq(alignment.query_sequence, "AACGTT")
  assert_eq(alignment.qualities, [10, 11, 12, 13, 14, 15])
  assert_eq(alignment.query_coordinates, [5, 1])
  assert_eq(alignment.hard_clip_left, Some(2))
  assert_eq(alignment.hard_clip_right, Some(1))
}

///|
test "Bio.Align.sam summary reports path statistics" {
  assert_eq(
    align_sam_test_forward().summary(),
    "SAM alignment read_forward -> chr1:4-14, aligned=9, inserted=1, deleted=1, skipped=0",
  )
  let unmapped = align_sam_test_parse(
      align_sam_test_record(
        "4", "*", "0", "0", "*", "*", "0", "0", "*", "*", "",
      ),
    ).alignments[0]
  assert_eq(unmapped.summary(), "SAM alignment read: unmapped")
}

///|
test "Bio.Align.sam preserves unknown sequence and qualities" {
  let alignment = align_sam_test_parse(
      align_sam_test_record(
        "0", "contig", "1", "10", "4M", "*", "0", "0", "*", "*", "",
      ),
    ).alignments[0]
  assert_false(alignment.sequence_known)
  assert_false(alignment.qualities_known)
  assert_eq(alignment.query_length(), 4)
  assert_eq(alignment.query_sequence, "")
  assert_eq(alignment.qualities.length(), 0)
}

///|
test "Bio.Align.sam rejects empty CIGAR text" {
  assert_true(align_sam_test_cigar_raises(""))
}

///|
test "Bio.Align.sam rejects CIGAR without a leading length" {
  assert_true(align_sam_test_cigar_raises("M4"))
}

///|
test "Bio.Align.sam rejects CIGAR without a trailing operation" {
  assert_true(align_sam_test_cigar_raises("4M2"))
}

///|
test "Bio.Align.sam rejects zero-length CIGAR operations" {
  assert_true(align_sam_test_cigar_raises("0M"))
}

///|
test "Bio.Align.sam rejects unknown CIGAR operations" {
  assert_true(align_sam_test_cigar_raises("4Q"))
}

///|
test "Bio.Align.sam rejects padding CIGAR operations" {
  assert_true(align_sam_test_cigar_raises("2M1P2M"))
}

///|
test "Bio.Align.sam rejects interior clipping" {
  assert_true(align_sam_test_cigar_raises("2M1S2M"))
  assert_true(align_sam_test_cigar_raises("2M1H2M"))
}

///|
test "Bio.Align.sam rejects invalid left clipping order" {
  assert_true(align_sam_test_cigar_raises("1S2H4M"))
}

///|
test "Bio.Align.sam rejects invalid right clipping order" {
  assert_true(align_sam_test_cigar_raises("4M2H1S"))
}

///|
test "Bio.Align.sam rejects clip-only CIGAR" {
  assert_true(align_sam_test_cigar_raises("2H4S"))
}

///|
test "Bio.Align.sam rejects headers after alignments" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_record(
        "4", "*", "0", "0", "*", "*", "0", "0", "*", "*", "",
      ) +
      "@HD\tVN:1.6\n",
    ),
  )
}

///|
test "Bio.Align.sam rejects HD after another header" {
  assert_true(align_sam_test_parse_raises("@SQ\tSN:chr1\tLN:10\n@HD\tVN:1.6\n"))
}

///|
test "Bio.Align.sam rejects duplicate HD" {
  assert_true(align_sam_test_parse_raises("@HD\tVN:1.6\n@HD\tVN:1.6\n"))
}

///|
test "Bio.Align.sam rejects HD without version" {
  assert_true(align_sam_test_parse_raises("@HD\tSO:coordinate\n"))
}

///|
test "Bio.Align.sam rejects malformed SQ metadata" {
  assert_true(align_sam_test_parse_raises("@SQ\tLN:10\n"))
  assert_true(align_sam_test_parse_raises("@SQ\tSN:chr1\n"))
  assert_true(align_sam_test_parse_raises("@SQ\tSN:chr1\tLN:0\n"))
  assert_true(align_sam_test_parse_raises("@SQ\tSN:chr1\tLN:10\tTP:branched\n"))
}

///|
test "Bio.Align.sam rejects duplicate references" {
  assert_true(
    align_sam_test_parse_raises("@SQ\tSN:chr1\tLN:10\n@SQ\tSN:chr1\tLN:20\n"),
  )
}

///|
test "Bio.Align.sam rejects duplicate header fields" {
  assert_true(align_sam_test_parse_raises("@HD\tVN:1.6\tVN:1.5\n"))
}

///|
test "Bio.Align.sam rejects alignment with too few fields" {
  assert_true(align_sam_test_parse_raises("read\t0\tchr1\n"))
}

///|
test "Bio.Align.sam validates FLAG range" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_record(
        "-1", "*", "0", "0", "*", "*", "0", "0", "*", "*", "",
      ),
    ),
  )
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_record(
        "65536", "*", "0", "0", "*", "*", "0", "0", "*", "*", "",
      ),
    ),
  )
}

///|
test "Bio.Align.sam validates MAPQ range" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_record(
        "0", "chr1", "1", "-1", "1M", "*", "0", "0", "A", "!", "",
      ),
    ),
  )
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_record(
        "0", "chr1", "1", "256", "1M", "*", "0", "0", "A", "!", "",
      ),
    ),
  )
}

///|
test "Bio.Align.sam validates mapped record fields" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_record(
        "0", "*", "1", "60", "1M", "*", "0", "0", "A", "!", "",
      ),
    ),
  )
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_record(
        "0", "chr1", "0", "60", "1M", "*", "0", "0", "A", "!", "",
      ),
    ),
  )
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_record(
        "0", "chr1", "1", "60", "*", "*", "0", "0", "A", "!", "",
      ),
    ),
  )
}

///|
test "Bio.Align.sam validates unmapped record fields" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_record(
        "4", "chr1", "1", "0", "*", "*", "0", "0", "A", "!", "",
      ),
    ),
  )
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_record(
        "4", "*", "0", "0", "1M", "*", "0", "0", "A", "!", "",
      ),
    ),
  )
}

///|
test "Bio.Align.sam rejects unknown references when SQ is present" {
  assert_true(
    align_sam_test_parse_raises(
      "@SQ\tSN:chr1\tLN:10\n" +
      align_sam_test_record(
        "0", "chr2", "1", "60", "1M", "*", "0", "0", "A", "!", "",
      ),
    ),
  )
}

///|
test "Bio.Align.sam rejects alignments beyond reference length" {
  assert_true(
    align_sam_test_parse_raises(
      "@SQ\tSN:chr1\tLN:10\n" +
      align_sam_test_record(
        "0", "chr1", "9", "60", "3M", "*", "0", "0", "ACG", "!!!", "",
      ),
    ),
  )
}

///|
test "Bio.Align.sam validates query length against CIGAR" {
  assert_true(
    align_sam_test_parse_raises(align_sam_test_mapped("4M", "ACG", "!!!", "")),
  )
}

///|
test "Bio.Align.sam rejects invalid sequence characters" {
  assert_true(
    align_sam_test_parse_raises(align_sam_test_mapped("4M", "AC*G", "!!!!", "")),
  )
}

///|
test "Bio.Align.sam requires QUAL and SEQ lengths to match" {
  assert_true(
    align_sam_test_parse_raises(align_sam_test_mapped("4M", "ACGT", "!!!", "")),
  )
  assert_true(
    align_sam_test_parse_raises(align_sam_test_mapped("4M", "*", "!!!!", "")),
  )
}

///|
test "Bio.Align.sam validates PHRED text range" {
  assert_true(
    align_sam_test_parse_raises(align_sam_test_mapped("1M", "A", " ", "")),
  )
}

///|
test "Bio.Align.sam rejects negative mate positions" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_record(
        "0", "chr1", "1", "60", "1M", "*", "-1", "0", "A", "!", "",
      ),
    ),
  )
}

///|
test "Bio.Align.sam rejects equal mate reference for unmapped record" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_record(
        "4", "*", "0", "0", "*", "=", "0", "0", "*", "*", "",
      ),
    ),
  )
}

///|
test "Bio.Align.sam rejects duplicate optional tags" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_mapped("1M", "A", "!", "NM:i:0\tNM:i:1"),
    ),
  )
}

///|
test "Bio.Align.sam validates optional tag names" {
  assert_true(
    align_sam_test_parse_raises(align_sam_test_mapped("1M", "A", "!", "A_:i:1")),
  )
  assert_true(
    align_sam_test_parse_raises(align_sam_test_mapped("1M", "A", "!", "1A:i:1")),
  )
}

///|
test "Bio.Align.sam rejects unknown optional tag types" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_mapped("1M", "A", "!", "ZZ:q:value"),
    ),
  )
}

///|
test "Bio.Align.sam validates character tags" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_mapped("1M", "A", "!", "AA:A:XY"),
    ),
  )
}

///|
test "Bio.Align.sam validates hexadecimal tags" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_mapped("1M", "A", "!", "HX:H:ABC"),
    ),
  )
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_mapped("1M", "A", "!", "HX:H:AZ"),
    ),
  )
}

///|
test "Bio.Align.sam validates array tag subtypes" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_mapped("1M", "A", "!", "BA:B:q,1"),
    ),
  )
  assert_true(
    align_sam_test_parse_raises(align_sam_test_mapped("1M", "A", "!", "BA:B:i")),
  )
}

///|
test "Bio.Align.sam validates array tag integer ranges" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_mapped("1M", "A", "!", "BA:B:c,128"),
    ),
  )
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_mapped("1M", "A", "!", "BA:B:C,-1"),
    ),
  )
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_mapped("1M", "A", "!", "BA:B:s,32768"),
    ),
  )
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_mapped("1M", "A", "!", "BA:B:S,65536"),
    ),
  )
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_mapped("1M", "A", "!", "BA:B:I,-1"),
    ),
  )
}

///|
test "Bio.Align.sam requires MD to use string datatype" {
  assert_true(
    align_sam_test_parse_raises(align_sam_test_mapped("1M", "A", "!", "MD:i:1")),
  )
}

///|
test "Bio.Align.sam requires AS to use integer datatype" {
  assert_true(
    align_sam_test_parse_raises(align_sam_test_mapped("1M", "A", "!", "AS:Z:1")),
  )
}

///|
test "Bio.Align.sam validates MD total reference consumption" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_mapped("4M", "ACGT", "!!!!", "MD:Z:3"),
    ),
  )
}

///|
test "Bio.Align.sam validates MD token grammar" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_mapped("1M", "A", "!", "MD:Z:A0"),
    ),
  )
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_mapped("1M", "A", "!", "MD:Z:0A"),
    ),
  )
}

///|
test "Bio.Align.sam aligns MD deletions with CIGAR D" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_mapped("1M", "A", "!", "MD:Z:0^A0"),
    ),
  )
  assert_true(
    align_sam_test_parse_raises(align_sam_test_mapped("1D", "*", "*", "MD:Z:1")),
  )
}

///|
test "Bio.Align.sam aligns MD mismatches with equal and mismatch CIGAR" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_mapped("1=", "A", "!", "MD:Z:0C0"),
    ),
  )
  assert_true(
    align_sam_test_parse_raises(align_sam_test_mapped("1X", "A", "!", "MD:Z:1")),
  )
  let alignment = align_sam_test_parse(
      align_sam_test_mapped("1X", "A", "!", "MD:Z:0C0"),
    ).alignments[0]
  assert_eq(alignment.computed_nm(), Some(1))
}

///|
test "Bio.Align.sam rejects MD on unmapped records" {
  assert_true(
    align_sam_test_parse_raises(
      align_sam_test_record(
        "4", "*", "0", "0", "*", "*", "0", "0", "A", "!", "MD:Z:0",
      ),
    ),
  )
}

///|
test "Bio.Align.sam aligned rows reject unmapped records" {
  let alignment = align_sam_test_parse(
      align_sam_test_record(
        "4", "*", "0", "0", "*", "*", "0", "0", "A", "!", "",
      ),
    ).alignments[0]
  let failed = try {
    ignore(alignment.aligned_rows())
    false
  } catch {
    AlignSamError(_) => true
  }
  assert_true(failed)
}

///|
test "Bio.Align.sam validates reference length for row reconstruction" {
  let failed = try {
    ignore(align_sam_test_forward().aligned_rows(reference_sequence="ACGT"))
    false
  } catch {
    AlignSamError(_) => true
  }
  assert_true(failed)
}

///|
test "Bio.Align.sam calculate MD requires mapped sequence and reference" {
  let unknown = align_sam_test_parse(
      align_sam_test_record(
        "0", "contig", "1", "60", "1M", "*", "0", "0", "*", "*", "",
      ),
    ).alignments[0]
  let failed_unknown = try {
    ignore(@src.align_sam_calculate_md(unknown, "A"))
    false
  } catch {
    AlignSamError(_) => true
  }
  assert_true(failed_unknown)
  let failed_length = try {
    ignore(
      @src.align_sam_calculate_md(align_sam_test_forward(), "A".repeat(39)),
    )
    false
  } catch {
    AlignSamError(_) => true
  }
  assert_true(failed_length)
}

///|
test "Bio.Align.sam format requires reference when calculating MD" {
  let failed = try {
    ignore(@src.align_sam_format(align_sam_test_forward(), calculate_md=true))
    false
  } catch {
    AlignSamError(_) => true
  }
  assert_true(failed)
}

///|
test "Bio.Align.sam constructor validates reference geometry" {
  let bad_length = try {
    ignore(@src.align_sam_create("r", "chr1", 0, 0, "A", "1M"))
    false
  } catch {
    AlignSamError(_) => true
  }
  let bad_start = try {
    ignore(@src.align_sam_create("r", "chr1", 10, -1, "A", "1M"))
    false
  } catch {
    AlignSamError(_) => true
  }
  let overrun = try {
    ignore(@src.align_sam_create("r", "chr1", 10, 9, "AA", "2M"))
    false
  } catch {
    AlignSamError(_) => true
  }
  assert_true(bad_length)
  assert_true(bad_start)
  assert_true(overrun)
}

///|
test "Bio.Align.sam constructor validates qualities and flags" {
  let bad_quality = try {
    ignore(
      @src.align_sam_create("r", "chr1", 10, 0, "AA", "2M", qualities=[20]),
    )
    false
  } catch {
    AlignSamError(_) => true
  }
  let bad_flag = try {
    ignore(@src.align_sam_create("r", "chr1", 10, 0, "A", "1M", flag=0x4))
    false
  } catch {
    AlignSamError(_) => true
  }
  let bad_mapq = try {
    ignore(@src.align_sam_create("r", "chr1", 10, 0, "A", "1M", mapq=256))
    false
  } catch {
    AlignSamError(_) => true
  }
  assert_true(bad_quality)
  assert_true(bad_flag)
  assert_true(bad_mapq)
}
