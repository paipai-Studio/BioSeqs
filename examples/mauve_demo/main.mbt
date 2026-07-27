///|
/// Mauve Genome Alignment Demo.
/// Demonstrates parsing Mauve alignments, detecting inversions
/// and breakpoints, computing coverage, and exporting results.

fn main {
  println("=== Mauve Genome Alignment Demo ===\n")

  println("1. Parsing Mauve alignment data...")
  let mauve_content = "##mauve version=1.0\n\n" +
    "a score=100\n" +
    "s\thuman.chr1\t0\t500\t+\t250000000\n" +
    "s\tmouse.chr1\t0\t480\t+\t240000000\n" +
    "s\trat.chr1\t0\t450\t+\t230000000\n\n" +
    "a score=80\n" +
    "s\thuman.chr1\t1000\t300\t-\t250000000\n" +
    "s\tmouse.chr1\t1200\t280\t-\t240000000\n" +
    "s\trat.chr1\t1100\t260\t-\t230000000\n\n" +
    "a score=60\n" +
    "s\thuman.chr1\t5000\t400\t+\t250000000\n" +
    "s\tmouse.chr1\t5200\t380\t+\t240000000\n" +
    "s\trat.chr1\t5100\t360\t+\t230000000\n\n"
  let ali = @src.parse_mauve(mauve_content)
  println("   Version: \{ali.version}")
  println("   Number of LCBs: \{@src.mauve_num_lcbs(ali)}")
  let all_names = @src.mauve_all_seq_names(ali)
  println("   Genomes: \{all_names}")

  println("\n2. Detecting inversions...")
  let inv_count = ali.inversions.length()
  println("   Detected inversions: \{inv_count}")
  let mut i = 0
  while i < ali.inversions.length() {
    let inv = ali.inversions[i]
    println("   Inversion \{i}: \{inv.seq_name} \{inv.start}-\{inv.end} (\{inv.length} bp)")
    println("     Affected LCBs: \{inv.affected_lcbs}")
    i = i + 1
  }

  println("\n3. Detecting breakpoints...")
  let bp_count = ali.breakpoints.length()
  println("   Detected breakpoints: \{bp_count}")
  let mut j = 0
  while j < ali.breakpoints.length() {
    let bp = ali.breakpoints[j]
    println("   Breakpoint \{j}: \{bp.lcb1_id} -> \{bp.lcb2_id} at \{bp.seq_name}:\{bp.position} (\{bp.breakpoint_type})")
    j = j + 1
  }

  println("\n4. Computing genome coverage...")
  let coverage = @src.mauve_genome_coverage(ali)
  let cov_keys = coverage.keys().collect()
  let mut k = 0
  while k < cov_keys.length() {
    let key = cov_keys[k]
    println("   \{key}: \{coverage[key]}% covered")
    k = k + 1
  }

  println("\n5. Generating summary report...")
  let summary = @src.mauve_summary(ali)
  println("\{summary}")

  println("6. BED format export...")
  let bed = @src.mauve_to_bed(ali)
  let bed_lines = bed.split("\n").to_array()
  println("   BED output lines: \{bed_lines.length()}")
  println("   First line: \{bed_lines[0]}")
  let inv_bed = @src.mauve_inversions_to_bed(ali)
  let inv_bed_lines = inv_bed.split("\n").to_array()
  println("   Inversions BED lines: \{inv_bed_lines.length()}")

  println("\n7. Filtering LCBs by score...")
  let filtered = @src.mauve_filter_lcbs(ali, 70.0)
  println("   Original LCBs: \{@src.mauve_num_lcbs(ali)}")
  println("   Filtered LCBs (score >= 70): \{@src.mauve_num_lcbs(filtered)}")

  println("\n8. Additional analyses...")
  let segments = @src.mauve_conserved_segments(ali)
  let seg_keys = segments.keys().collect()
  let mut m = 0
  while m < seg_keys.length() {
    let key = seg_keys[m]
    println("   Conserved segments for \{key}: \{segments[key]}")
    m = m + 1
  }
  let rate = @src.mauve_rearrangement_rate(ali)
  println("   Rearrangement rate: \{rate} breakpoints/Mb")

  println("\n=== Mauve Demo Complete ===")
}