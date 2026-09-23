///|
fn integer_array_text(values : Array[Int]) -> String {
  let output = StringBuilder()
  output.write_char('[')
  for index = 0; index < values.length(); index = index + 1 {
    if index > 0 {
      output.write_string(", ")
    }
    output.write_string(values[index].to_string())
  }
  output.write_char(']')
  output.to_string()
}

///|
fn main {
  println("=== Biopython Bio.Align.msf Demo ===")

  let alignment = @src.msf_parse(@src.msf_example_text()) catch {
    MsfError(message) => abort("MSF parsing failed: " + message)
  }

  println("\n1. Parse interleaved MSF metadata and rows")
  println("  " + alignment.summary())
  println("  title: " + alignment.metadata.title)
  for sequence in alignment.sequences {
    println(
      "  " +
      sequence.id +
      ": residues=" +
      sequence.length.to_string() +
      ", weight=" +
      sequence.weight.to_string() +
      ", checksum=" +
      sequence.checksum.to_string(),
    )
  }

  println("\n2. Build a Biopython-style coordinate path")
  let path = alignment.coordinate_path()
  for row = 0; row < alignment.num_sequences(); row = row + 1 {
    println(
      "  " + alignment.sequences[row].id + ": " + integer_array_text(path[row]),
    )
  }
  println(
    "  reference residue 4 -> query_one residue " +
    alignment.map_position(0, 1, 4).unwrap_or(-1).to_string(),
  )
  println(
    "  reference residue 3 -> query_one residue " +
    alignment.map_position(0, 1, 3).unwrap_or(-1).to_string() +
    " (-1 denotes a gap)",
  )

  println("\n3. Calculate pair statistics and consensus")
  let counts = alignment.pair_counts(0, 2) catch {
    MsfError(message) => abort("pair counting failed: " + message)
  }
  println(
    "  aligned=" +
    counts.aligned.to_string() +
    ", identities=" +
    counts.identities.to_string() +
    ", mismatches=" +
    counts.mismatches.to_string() +
    ", gap opens=" +
    counts.gap_opens.to_string(),
  )
  let consensus = alignment.consensus(minimum_fraction=0.67) catch {
    MsfError(message) => abort("consensus failed: " + message)
  }
  println("  consensus: " + consensus)
  println("  occupancy at column 1: " + alignment.occupancy()[1].to_string())

  println("\n4. Write canonical MSF and parse it again")
  let canonical = @src.msf_write(
    alignment,
    block_width=12,
    group_width=4,
    gap_character="~",
  ) catch {
    MsfError(message) => abort("MSF serialization failed: " + message)
  }
  let round_trip = @src.msf_parse(canonical) catch {
    MsfError(message) => abort("round-trip parsing failed: " + message)
  }
  println("  output bytes: " + canonical.length().to_string())
  println("  checksums valid: " + round_trip.checksums_valid().to_string())
  println(
    "  aligned rows preserved: " +
    (round_trip.sequences[1].aligned_sequence ==
    alignment.sequences[1].aligned_sequence).to_string(),
  )

  println("\n5. Preserve third-party declared-width mismatches")
  let mismatch_text = "!!AA_MULTIPLE_ALIGNMENT\n" +
    "MSF: 2 Type: P Check: 0 ..\n" +
    "Name: full Len: 4 Check: 0 Weight: 1.0\n" +
    "Name: short Len: 2 Check: 0 Weight: 1.0\n" +
    "//\n\nfull ACDE\nshort AC\n"
  let mismatch = @src.msf_parse(mismatch_text) catch {
    MsfError(message) => abort("length-mismatch parsing failed: " + message)
  }
  println(
    "  declared=" +
    mismatch.metadata.declared_length.to_string() +
    ", actual=" +
    mismatch.alignment_length().to_string() +
    ", mismatch=" +
    mismatch.length_mismatch().to_string(),
  )
}
