///|
/// CRAM format parsing example.
///
/// Demonstrates CRAM format parsing, checking magic numbers,
/// and available API functions.
///
/// Run: moon run examples/cram_demo/main.mbt
fn main {
  println("=== CRAM Format Parsing Demo ===")
  println("")

  println("1. Checking CRAM magic number...")
  let valid_cram = "CRAM"
  let invalid = "BAM\u{0001}"
  println("   Is CRAM magic ('CRAM'): \{@bio.is_cram_magic(valid_cram)}")
  println("   Is CRAM magic ('BAM...'): \{@bio.is_cram_magic(invalid)}")
  println("")

  println("2. CRAM Data Structures Overview...")
  println("")

  println("   CramVersion:")
  println("     - major: Int")
  println("     - minor: Int")
  println("")

  println("   CramCompression:")
  println("     - None: No compression")
  println("     - Gzip: Gzip compression")
  println("     - Bzip2: Bzip2 compression")
  println("     - Lzma: Lzma compression")
  println("     - Unknown(code): Unknown compression type")
  println("")

  println("   CramReference:")
  println("     - name: String")
  println("     - length: Int")
  println("     - md5: String")
  println("     - uri: String")
  println("")

  println("   CramHeader:")
  println("     - version: CramVersion")
  println("     - reference_count: Int")
  println("     - references: Array[CramReference]")
  println("     - compression: CramCompression")
  println("     - stream_count: Int")
  println("     - stream_info: Array[String]")
  println("     - sam_header: String")
  println("")

  println("   CramReadFeatures:")
  println("     - read_group: Int")
  println("     - flags: Int")
  println("     - mapping_quality: Int")
  println("     - read_name: String")
  println("     - cigar: Array[BamCigarElement]")
  println("     - sequence: String")
  println("     - quality: String")
  println("     - tags: Array[String]")
  println("")

  println("   CramRecord:")
  println("     - ref_id: Int")
  println("     - pos: Int")
  println("     - read_features: CramReadFeatures")
  println("")

  println("   CramFile:")
  println("     - header: CramHeader")
  println("     - records: Array[CramRecord]")
  println("")

  println("3. Available CRAM Functions...")
  println("")

  println("   Magic Number Detection:")
  println("     is_cram_magic(data): Bool")
  println("")

  println("   Parsing Functions:")
  println("     parse_cram(data): CramFile raise CramError")
  println("     parse_cram_base64(base64): CramFile raise CramError")
  println("     try_parse_cram_base64(base64): CramFile?")
  println("")

  println("   CramFile Methods:")
  println("     .len(): Int")
  println("     .num_references(): Int")
  println("     .get_reference_names(): Array[String]")
  println("     .get_reference_length(name): Int")
  println("     .to_bam(): BamFile")
  println("")

  println("   CramRecord Methods:")
  println("     .is_paired(): Bool")
  println("     .is_proper_pair(): Bool")
  println("     .is_unmapped(): Bool")
  println("     .is_reverse(): Bool")
  println("     .is_read1(): Bool")
  println("     .is_read2(): Bool")
  println("     .get_reference_name(header): String")
  println("")

  println("4. CRAM to BAM Conversion...")
  println("")
  println("   CRAM files can be converted to BAM for interoperability:")
  println("     let cram = parse_cram(data)")
  println("     let bam = cram.to_bam()")
  println("     let sam = bam.to_sam()")
  println("")

  println("5. CRAM vs BAM Comparison...")
  println("")
  println("   Feature          | CRAM         | BAM")
  println("   ---------------- |-------------|-------")
  println("   Magic Number     | CRAM\\0       | BAM\\1")
  println("   Compression      | Built-in    | BGZF")
  println("   Reference        | Embedded    | External")
  println("   Compression Ratio| Higher      | Lower")
  println("   Read Speed       | Slower      | Faster")
  println("   File Size        | Smaller     | Larger")
  println("")

  println("6. Usage Example...")
  println("")
  println("   let cram_data = ...")
  println("   if is_cram_magic(cram_data) {")
  println("     let cram = parse_cram(cram_data)")
  println("     println('References: ' ++ cram.num_references().show())")
  println("     println('Records: ' ++ cram.len().show())")
  println("     let bam = cram.to_bam()")
  println("     // Process BAM records...")
  println("   }")
  println("")

  println("=== Done ===")
}
