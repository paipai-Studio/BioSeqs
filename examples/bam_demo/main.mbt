///|
/// BAM format parsing example.
///
/// Demonstrates parsing raw BAM binary data and accessing record information.
///
/// Run: moon run examples/bam_demo/main.mbt
fn main {
  println("=== BAM Format Parsing Demo ===")
  println("")
  
  let test_bam_base64 = "QkFNAS0AAABASEQJVk46MS42CVNPOmNvb3JkaW5hdGUKQFNRCVNOOmNocjEJTE46MTAwMAoBAAAABQAAAGNocjEA6AMAADkAAAAAAAAAZAAAAAAAAAA8AAAGAAEACAAAAP////8AAAAAAAAAAHJlYWQxAIAAAAAAAAAAAAAAAAAAAAA="
  
  println("Parsing base64-encoded BAM data...")
  
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
  
  println("=== Done ===")
}
