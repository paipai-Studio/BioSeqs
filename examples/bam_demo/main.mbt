///|
/// BAM format parsing example.
///
/// Demonstrates parsing raw BAM binary data, BGZF decompression, and accessing record information.
///
/// Run: moon run examples/bam_demo/main.mbt
fn main {
  println("=== BAM Format Parsing Demo ===")
  println("")

  let test_bam_base64 = "QkFNAS0AAABASEQJVk46MS42CVNPOmNvb3JkaW5hdGUKQFNRCVNOOmNocjEJTE46MTAwMAoBAAAABQAAAGNocjEA6AMAADkAAAAAAAAAZAAAAAAAAAA8AAAGAAEACAAAAP////8AAAAAAAAAAHJlYWQxAIAAAAAAAAAAAAAAAAAAAAA="

  println("1. Parsing base64-encoded BAM data...")

  let bam_opt = @bio.try_parse_bam_base64(test_bam_base64)
  if bam_opt is None {
    println("Error parsing BAM data")
    return
  }
  let bam = bam_opt.unwrap()

  println("Parsing successful!")
  println("")

  println("=== Header Information ===")
  println("References: \{bam.num_references()}")
  let ref_names = bam.get_reference_names()
  for name in ref_names {
    let len = bam.get_reference_length(name)
    println("  - \{name}: length = \{len}")
  }

  println("")
  println("=== Records ===")
  println("Total records: \{bam.len()}")
  println("")

  for i = 0; i < bam.len(); i = i + 1 {
    let record = bam.records[i]
    println("Record \{i + 1}:")
    println("  QNAME: \{record.qname}")
    println("  FLAG: \{record.flag}")
    println("  RNAME: \{record.get_reference_name(bam.header)}")
    println("  POS: \{record.pos + 1} (1-based)")
    println("  MAPQ: \{record.mapq}")
    println("  CIGAR: \{record.cigar.length()} operations")
    for j = 0; j < record.cigar.length(); j = j + 1 {
      let op = record.cigar[j]
      let op_str = match op.op {
        @bio.BamCigarOp::Match => "M"
        @bio.BamCigarOp::Insert => "I"
        @bio.BamCigarOp::Delete => "D"
        @bio.BamCigarOp::Skip => "N"
        @bio.BamCigarOp::SoftClip => "S"
        @bio.BamCigarOp::HardClip => "H"
        @bio.BamCigarOp::Pad => "P"
        @bio.BamCigarOp::Equal => "="
        @bio.BamCigarOp::Diff => "X"
      }
      println("    \{op.length}\{op_str}")
    }
    println("  SEQ: \{record.seq}")
    println("  QUAL: \{record.qual}")
    println("  TLEN: \{record.tlen}")
    println("  END: \{record.get_end() + 1} (1-based)")
    println("  Is paired: \{record.is_paired()}")
    println("  Is unmapped: \{record.is_unmapped()}")
    println("")
  }

  println("=== SAM Output ===")
  let sam_str = bam.to_sam()
  println(sam_str)

  println("")
  println("=== BGZF Decompression Demo ===")
  println("")

  let bgzf_magic = "\u001f\u008b\u0008\u0004"
  println("BGZF magic bytes detection:")
  println("  Is BGZF (magic): \{@bio.is_bgzf(bgzf_magic)}")
  let bam_magic = "BAM\u0001"
  println("  Is BGZF (BAM): \{@bio.is_bgzf(bam_magic)}")
  println("")

  println("BGZF functions available:")
  println("  - is_bgzf(data): Check if data is BGZF compressed")
  println("  - bgzf_decompress(data): Decompress BGZF data")
  println("  - parse_bam_from_bgzf(data): Parse BAM from BGZF compressed data")
  println("")

  println("=== Done ===")
}
