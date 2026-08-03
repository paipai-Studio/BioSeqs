///|
/// PHD module demo - Phred base-calling quality file parser.
///
/// Run: moon run examples/phd_demo
fn main {
  println("=== PHD Module Demo ===")
  println("Phred base-calling output file parser")
  println("Reference: Ewing B et al. (1998) Genome Res 8(3):175-85")
  println("")

  // ------------------------------------------------------------------
  // 1. Generate and parse sample PHD file
  // ------------------------------------------------------------------
  println("=== 1. Generate and parse sample PHD file ===")
  let text = @bio.phd_sample_text()
  let file = @bio.phd_parse(text)
  println("Parsed " + file.num_reads().to_string() + " read(s)")
  println("")

  // ------------------------------------------------------------------
  // 2. Display read information
  // ------------------------------------------------------------------
  println("=== 2. Read information ===")
  let read = file.get(0).unwrap()
  println("Read name: " + read.name())
  println("Length: " + read.length().to_string() + " bases")
  println("Sequence: " + read.sequence())
  println("")

  // ------------------------------------------------------------------
  // 3. Display quality scores
  // ------------------------------------------------------------------
  println("=== 3. Quality scores ===")
  let quals = read.quality()
  for i in 0..<quals.length() {
    let marker = if quals[i] >= 40 { " (high)" } else if quals[i] >= 20 { " (medium)" } else { " (low)" }
    println("  base " + (i + 1).to_string() + ": Q=" + quals[i].to_string() + marker)
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Display peak positions
  // ------------------------------------------------------------------
  println("=== 4. Peak positions (trace indices) ===")
  let peaks = read.peak_positions()
  for i in 0..<peaks.length() {
    println("  base " + (i + 1).to_string() + ": peak=" + peaks[i].to_string())
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Display comment metadata
  // ------------------------------------------------------------------
  println("=== 5. Comment metadata ===")
  let c = read.comment()
  println("Chromatogram file: " + c.chromat_file())
  println("Phred version: " + c.phred_version())
  println("Trace index range: " + c.trace_index_min().to_string() + "-" + c.trace_index_max().to_string())
  println("Trim: " + c.trim_start().to_string() + "-" + c.trim_end().to_string() + " (prob=" + c.trim_prob().to_string() + ")")
  println("Chemistry: " + c.chemistry())
  println("Dye: " + c.dye())
  println("")

  // ------------------------------------------------------------------
  // 6. Parse minimal PHD file
  // ------------------------------------------------------------------
  println("=== 6. Parse minimal PHD file ===")
  let minimal = "BEGIN_SEQUENCE r1\nBEGIN_DNA\nA 10 1\nT 20 2\nG 30 3\nC 40 4\nEND_DNA\nEND_SEQUENCE\n"
  let min_file = @bio.phd_parse(minimal)
  let min_read = min_file.get(0).unwrap()
  println("Read: " + min_read.name())
  println("Sequence: " + min_read.sequence())
  println("Quality: " + min_read.quality().to_string())
  println("")

  // ------------------------------------------------------------------
  // 7. Full read summary
  // ------------------------------------------------------------------
  println("=== 7. Full read summary ===")
  println(read.to_string())

  println("=== Demo complete ===")
}
