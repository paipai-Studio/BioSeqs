///|
/// abi module demo - ABI Sanger sequencing trace file parsing.
///
/// Demonstrates parsing of Applied Biosystems (ABIF) binary trace files
/// containing called bases, quality values, peak positions, and raw
/// chromatogram traces from Sanger sequencing.
///
/// Reference: Biopython Bio.SeqIO.AbiIO
/// Run: moon run examples/abi_demo
fn main {
  println("=== abi Module Demo ===")
  println("ABI Sanger sequencing trace file parsing")
  println("Reference: Biopython Bio.SeqIO.AbiIO")
  println("")

  // ------------------------------------------------------------------
  // 1. Sample trace creation and basic fields
  // ------------------------------------------------------------------
  println("=== 1. Sample Trace Creation ===")
  let trace = @bio.abi_sample_trace()
  println("Sample name: " + trace.sample_name)
  println("Capillary ID: " + trace.capillary_id)
  println("Run name: " + trace.run_name)
  println("Base order (FWO_): " + trace.base_order)
  println("Called bases: " + trace.bases)
  println("Number of bases: " + trace.num_bases.to_string())
  println("Number of trace points: " + trace.num_points.to_string())
  println("")

  // ------------------------------------------------------------------
  // 2. Run date and time
  // ------------------------------------------------------------------
  println("=== 2. Run Date and Time ===")
  println("Run date: " + trace.run_date.to_string())
  println("Run time: " + trace.run_time.to_string())
  println("")

  // ------------------------------------------------------------------
  // 3. Quality values
  // ------------------------------------------------------------------
  println("=== 3. Quality Values (Phred) ===")
  let mut i = 0
  while i < trace.qualities.length() {
    let q = trace.quality_at(i)
    let base = trace.bases[i:i + 1].to_owned()
    let label = if q >= 40 {
      "high"
    } else if q >= 20 {
      "medium"
    } else {
      "low"
    }
    println("  pos " + i.to_string() + ": base=" + base + " Q=" + q.to_string() + " (" + label + ")")
    i = i + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Peak positions
  // ------------------------------------------------------------------
  println("=== 4. Peak Positions ===")
  let mut p = 0
  let parts = StringBuilder::new()
  while p < trace.base_positions.length() {
    if p > 0 {
      parts.write_string(", ")
    }
    parts.write_string(trace.base_positions[p].to_string())
    p = p + 1
  }
  println("Base positions: [" + parts.to_string() + "]")
  println("")

  // ------------------------------------------------------------------
  // 5. Trace channels
  // ------------------------------------------------------------------
  println("=== 5. Trace Channels ===")
  let mut ch = 0
  while ch < trace.traces.length() {
    let channel = trace.traces[ch]
    let name = if ch < trace.trace_names.length() { trace.trace_names[ch] } else { "ch" + ch.to_string() }
    // Find max value in channel
    let mut max_val = 0
    let mut k = 0
    while k < channel.length() {
      if channel[k] > max_val {
        max_val = channel[k]
      }
      k = k + 1
    }
    println("  " + name + ": " + channel.length().to_string() + " points, max=" + max_val.to_string())
    ch = ch + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Trace for specific base
  // ------------------------------------------------------------------
  println("=== 6. Trace for Specific Base ===")
  let g_trace = trace.trace_for_base("G")
  println("G trace length: " + g_trace.length().to_string())
  // Find peak positions in G trace
  println("G trace peaks (>500):")
  let mut peak_idx = 0
  while peak_idx < g_trace.length() {
    if g_trace[peak_idx] > 500 {
      println("  point " + peak_idx.to_string() + ": " + g_trace[peak_idx].to_string())
    }
    peak_idx = peak_idx + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Binary encoding/decoding round trip
  // ------------------------------------------------------------------
  println("=== 7. Binary Encode/Decode ===")
  let bytes = @bio.abi_encode(trace)
  println("Encoded bytes: " + bytes.length().to_string() + " bytes")
  println("Magic (first 4 bytes): " + bytes[0].to_string() + " " + bytes[1].to_string() + " " + bytes[2].to_string() + " " + bytes[3].to_string() + " (ABIF)")
  let parsed = @bio.abi_parse(@bio.abi_sample_bytes())
  println("Parsed from sample bytes: bases=" + parsed.bases + ", order=" + parsed.base_order)
  println("")

  // ------------------------------------------------------------------
  // 8. Sequence length and statistics
  // ------------------------------------------------------------------
  println("=== 8. Sequence Statistics ===")
  println("Sequence length: " + trace.seq_length().to_string())
  let mut q_sum = 0
  for q in trace.qualities {
    q_sum = q_sum + q
  }
  let avg_q = if trace.qualities.length() > 0 {
    q_sum / trace.qualities.length()
  } else {
    0
  }
  println("Average quality: " + avg_q.to_string())
  println("")

  println("=== Demo complete ===")
}
