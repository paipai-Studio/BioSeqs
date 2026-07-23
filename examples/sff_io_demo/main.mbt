///|
/// Bio.SeqIO.SffIO Demo - Standard Flowgram Format binary parsing.
/// Inspired by Biopython's Bio.SeqIO.SffIO module.

fn main {
  println("=== Bio.SeqIO.SffIO Demo ===")

  // 1. Create an example SFF file
  println("")
  println("1. Create an example SFF file:")
  let sff = @bio.sff_create_example_file()
  let header = sff.sff_header()
  println("   Flow chars: " + header.sff_flow_chars())
  println("   Key sequence: " + header.sff_key_sequence())
  println("   Number of reads: " + sff.sff_n_reads().to_string())
  println("   Total bases: " + sff.sff_total_bases().to_string())

  // 2. Inspect individual reads
  println("")
  println("2. Inspect reads:")
  let mut i = 0
  while i < sff.sff_n_reads() {
    let read = sff.sff_get_read(i)
    let mean_q = read.sff_mean_quality()
    println("   " + read.sff_name() +
      ": bases=" + read.sff_n_bases().to_string() +
      ", sequence=" + read.sff_bases() +
      ", mean_quality=" + mean_q.to_string())
    i = i + 1
  }

  // 3. Quality clipping
  println("")
  println("3. Quality clipping:")
  let read0 = sff.sff_get_read(0)
  let clipped = read0.sff_clipped_bases()
  let clipped_quals = read0.sff_clipped_qualities()
  println("   Original bases: " + read0.sff_bases())
  println("   Clipped bases: " + clipped)
  println("   Clipped qualities count: " + clipped_quals.length().to_string())

  // 4. Binary encode/decode round-trip
  println("")
  println("4. Binary encode/decode round-trip:")
  let bytes = @bio.sff_encode(sff)
  println("   Encoded " + bytes.length().to_string() + " bytes")
  let parsed_opt = @bio.sff_parse(bytes)
  if parsed_opt is Some(_) {
    let parsed = parsed_opt.unwrap()
    println("   Parsed " + parsed.sff_n_reads().to_string() + " reads back")
    let parsed_header = parsed.sff_header()
    println("   Flow chars match: " + (parsed_header.sff_flow_chars() == header.sff_flow_chars()).to_string())
    let r0 = parsed.sff_get_read(0)
    println("   Read 0 name: " + r0.sff_name())
    println("   Read 0 bases: " + r0.sff_bases())
  } else {
    println("   Parse failed!")
  }

  // 5. Header validation
  println("")
  println("5. Header validation:")
  let errors = header.sff_validate()
  println("   Validation errors: " + errors.length().to_string())

  // 6. Find read by name
  println("")
  println("6. Find read by name:")
  let found = sff.sff_find_read("READ0002")
  if found is Some(_) {
    println("   Found READ0002: " + found.unwrap().sff_bases())
  } else {
    println("   READ0002 not found")
  }

  println("")
  println("=== Bio.SeqIO.SffIO Demo Completed ===")
}
