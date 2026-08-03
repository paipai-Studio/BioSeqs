///|
/// GFA module demo - Graph Fragment Assembly format parser and writer.
///
/// Run: moon run examples/gfa_demo
fn main {
  println("=== GFA Module Demo ===")
  println("Graph Fragment Assembly format for genome assembly graphs")
  println("Reference: Li H et al. (2016) Bioinformatics 32:2203-2205")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample GFA graph
  // ------------------------------------------------------------------
  println("=== 1. Load sample GFA graph ===")
  let graph = @bio.gfa_sample_graph()
  println("Loaded sample GFA graph")
  println("")

  // ------------------------------------------------------------------
  // 2. Print graph summary (n_segments, n_links, n_paths)
  // ------------------------------------------------------------------
  println("=== 2. Graph summary ===")
  let n_segments = @bio.gfa_graph_n_segments(graph)
  let n_links = @bio.gfa_graph_n_links(graph)
  let n_paths = @bio.gfa_graph_n_paths(graph)
  println("n_segments: \{n_segments}")
  println("n_links: \{n_links}")
  println("n_paths: \{n_paths}")
  println("")
  println("Full summary: \{@bio.gfa_graph_summary(graph)}")
  println("")

  // ------------------------------------------------------------------
  // 3. Show all segments (name, sequence length, tags)
  // ------------------------------------------------------------------
  println("=== 3. Segments ===")
  let segments = graph.segments()
  for i in 0..<segments.length() {
    let seg = segments[i]
    let seq = seg.sequence()
    let seq_len = if seq == "*" { 0 } else { seq.length() }
    println("Segment \{i}:")
    println("  name: \{seg.name()}")
    println("  sequence: \{seq}")
    println("  sequence length: \{seq_len}")
    let tags = seg.tags()
    if tags.length() == 0 {
      println("  tags: (none)")
    } else {
      for t in 0..<tags.length() {
        println("  tag[\{t}]: \{tags[t].to_string()}")
      }
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Show all links (from, from_orient, to, to_orient, overlap)
  // ------------------------------------------------------------------
  println("=== 4. Links ===")
  let links = graph.links()
  if links.length() == 0 {
    println("(no links)")
  } else {
    for i in 0..<links.length() {
      let link = links[i]
      println("Link \{i}:")
      println("  from: \{link.from_segment()} (\{link.from_orient()})")
      println("  to:   \{link.to_segment()} (\{link.to_orient()})")
      println("  overlap: \{link.overlap()}")
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Show paths (name, segment list)
  // ------------------------------------------------------------------
  println("=== 5. Paths ===")
  let paths = graph.paths()
  if paths.length() == 0 {
    println("(no paths)")
  } else {
    for i in 0..<paths.length() {
      let path = paths[i]
      println("Path \{i}:")
      println("  name: \{path.name()}")
      let seg_list = path.segment_list()
      let overlaps = path.overlaps()
      let seg_sb = StringBuilder::new()
      for j in 0..<seg_list.length() {
        if j > 0 {
          seg_sb.write_string(", ")
        }
        seg_sb.write_string(seg_list[j])
      }
      println("  segments: \{seg_sb.to_string()}")
      let ov_sb = StringBuilder::new()
      for j in 0..<overlaps.length() {
        if j > 0 {
          ov_sb.write_string(", ")
        }
        ov_sb.write_string(overlaps[j])
      }
      println("  overlaps: \{ov_sb.to_string()}")
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Convert segments to SeqRecord array and show
  // ------------------------------------------------------------------
  println("=== 6. Convert segments to SeqRecord array ===")
  let records = @bio.gfa_get_segments_as_records(graph)
  println("Converted \{records.length()} segment(s) to SeqRecord(s)")
  for i in 0..<records.length() {
    let rec = records[i]
    println("  Record \{i}: id=\{rec.id}, name=\{rec.name}, length=\{rec.seq.to_string().length()}")
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Serialize graph to GFA text and print first few lines
  // ------------------------------------------------------------------
  println("=== 7. Serialize graph to GFA text ===")
  let gfa_text = @bio.gfa_to_string(graph)
  let lines = gfa_text.split("\n").to_array()
  let n_lines = lines.length()
  // Print at most the first 5 non-empty lines.
  let mut shown = 0
  let mut idx = 0
  while idx < n_lines && shown < 5 {
    let line = lines[idx]
    if line.length() > 0 {
      println("  |\{line}")
      shown = shown + 1
    }
    idx = idx + 1
  }
  println("  ... (total \{n_lines} line(s) in serialized output)")
  println("")

  // ------------------------------------------------------------------
  // 8. Parse a GFA string and verify round-trip
  // ------------------------------------------------------------------
  println("=== 8. Parse a GFA string and verify round-trip ===")
  let gfa_string = "H\tVN:Z:1.0\nS\tquery1\tACGTACGT\tLN:i:8\nS\tquery2\tTTTTGGGG\tLN:i:8\nL\tquery1\t+\tquery2\t-\t4M\nP\tpathA\tquery1+,query2-\t4M,4M\n"
  println("Input GFA string:")
  println(gfa_string)
  let parsed = @bio.gfa_parse(gfa_string)
  println("Parsed graph:")
  println("  segments: \{@bio.gfa_graph_n_segments(parsed)}")
  println("  links: \{@bio.gfa_graph_n_links(parsed)}")
  println("  paths: \{@bio.gfa_graph_n_paths(parsed)}")
  // Round-trip: serialize the parsed graph and re-parse, then compare counts.
  let re_text = @bio.gfa_to_string(parsed)
  let re_parsed = @bio.gfa_parse(re_text)
  let ok = @bio.gfa_graph_n_segments(re_parsed) ==
    @bio.gfa_graph_n_segments(parsed) &&
    @bio.gfa_graph_n_links(re_parsed) ==
    @bio.gfa_graph_n_links(parsed) &&
    @bio.gfa_graph_n_paths(re_parsed) ==
    @bio.gfa_graph_n_paths(parsed)
  println("Round-trip verify: \{ok}")
  println("")

  println("=== Demo complete ===")
}
