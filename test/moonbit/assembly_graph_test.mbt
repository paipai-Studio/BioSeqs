///|
/// Tests for assembly graph simplification and racon-lite polishing
/// (src/assembly_graph.mbt).

// ============================================================================
// Helpers
// ============================================================================

fn test_repeat(base : Char, len : Int32) -> String {
  let sb = StringBuilder()
  for _i in 0..<len {
    sb.write_char(base)
  }
  sb.to_string()
}

///|
fn test_ref() -> String {
  "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATT"
}

///|
fn test_mutate_base(seq : String, pos : Int32) -> String {
  let sb = StringBuilder()
  let mut i = 0
  while i < seq.length() {
    let c = seq.unsafe_get(i).unsafe_to_char()
    if i == pos {
      let replacement = if c == 'A' {
        'C'
      } else if c == 'C' {
        'G'
      } else if c == 'G' {
        'T'
      } else {
        'A'
      }
      sb.write_char(replacement)
    } else {
      sb.write_char(c)
    }
    i = i + 1
  }
  sb.to_string()
}

///|
fn test_revcomp(seq : String) -> String {
  let sb = StringBuilder()
  let mut i = seq.length() - 1
  while i >= 0 {
    let c = seq.unsafe_get(i).unsafe_to_char()
    let comp = if c == 'A' {
      'T'
    } else if c == 'C' {
      'G'
    } else if c == 'G' {
      'C'
    } else {
      'A'
    }
    sb.write_char(comp)
    i = i - 1
  }
  sb.to_string()
}

///|
fn test_segment(name : String, len : Int32, depth : Int32) -> @src.GfaSegment {
  let tags : Array[@src.GfaTag] = if depth >= 0 {
    [@src.GfaTag::new("DP", "i", depth.to_string())]
  } else {
    []
  }
  @src.GfaSegment::new(name, test_repeat('A', len), tags)
}

///|
fn test_link(from : String, to : String) -> @src.GfaLink {
  @src.GfaLink::new(from, "+", to, "+", "0M", [])
}

///|
fn test_tip_graph(tip_len : Int32) -> @src.GfaGraph {
  let g = @src.GfaGraph::new()
  g.add_segment(test_segment("a", 200, -1))
  g.add_segment(test_segment("b", 200, -1))
  g.add_segment(test_segment("c", 200, -1))
  g.add_segment(test_segment("tip", tip_len, -1))
  g.add_link(test_link("a", "b"))
  g.add_link(test_link("b", "c"))
  g.add_link(test_link("b", "tip"))
  g
}

///|
fn test_bubble_graph() -> @src.GfaGraph {
  let g = @src.GfaGraph::new()
  g.add_segment(test_segment("entry", 50, -1))
  g.add_segment(test_segment("a1", 100, 10))
  g.add_segment(test_segment("a2", 100, 10))
  g.add_segment(test_segment("b1", 100, 2))
  g.add_segment(test_segment("exit", 50, -1))
  g.add_link(test_link("entry", "a1"))
  g.add_link(test_link("a1", "a2"))
  g.add_link(test_link("a2", "exit"))
  g.add_link(test_link("entry", "b1"))
  g.add_link(test_link("b1", "exit"))
  g
}

///|
fn test_equal_support_bubble_graph() -> @src.GfaGraph {
  let g = @src.GfaGraph::new()
  g.add_segment(test_segment("entry", 50, 1))
  g.add_segment(test_segment("a1", 100, 1))
  g.add_segment(test_segment("b1", 150, 1))
  g.add_segment(test_segment("exit", 50, 1))
  g.add_link(test_link("entry", "a1"))
  g.add_link(test_link("a1", "exit"))
  g.add_link(test_link("entry", "b1"))
  g.add_link(test_link("b1", "exit"))
  g
}

///|
fn test_simplify_with_tip_len(
  graph : @src.GfaGraph,
  min_tip_len : Int32,
) -> @src.AssemblySimplification? {
  try {
    let cfg = @src.assembly_simplify_config(min_tip_len~)
    Some(@src.assembly_simplify(graph, config=cfg))
  } catch {
    _ => None
  }
}

///|
fn test_polish_default(
  reference : String,
  reads : Array[String],
) -> @src.PolishResult? {
  try {
    let cfg = @src.polish_config(k=15, min_coverage=3, min_fraction=0.6)
    Some(@src.racon_polish(reference, reads, config=cfg))
  } catch {
    _ => None
  }
}

///|
fn test_assert_segment_present(
  graph : @src.GfaGraph,
  name : String,
) -> Unit raise {
  match @src.gfa_get_segment(graph, name) {
    Some(seg) => assert_eq(seg.name(), name)
    None => assert_true(false)
  }
}

///|
fn test_assert_segment_absent(
  graph : @src.GfaGraph,
  name : String,
) -> Unit raise {
  match @src.gfa_get_segment(graph, name) {
    Some(_) => assert_true(false)
    None => ()
  }
}

///|
fn test_reads(seq : String, n : Int32) -> Array[String] {
  let reads : Array[String] = Array(capacity=n)
  for _i in 0..<n {
    reads.push(seq)
  }
  reads
}

// ============================================================================
// Segment tag helpers
// ============================================================================

///|
test "assembly_segment_length prefers sequence and falls back to LN tag" {
  let with_seq = @src.GfaSegment::new("s1", "ACGTACGT", [])
  assert_eq(@src.assembly_segment_length(with_seq), 8)
  let tags : Array[@src.GfaTag] = [@src.GfaTag::new("LN", "i", "42")]
  let no_seq = @src.GfaSegment::new("s2", "*", tags)
  assert_eq(@src.assembly_segment_length(no_seq), 42)
}

///|
test "assembly_segment_depth parses DP then KC tags" {
  let dp : Array[@src.GfaTag] = [@src.GfaTag::new("DP", "i", "12")]
  assert_eq(
    @src.assembly_segment_depth(@src.GfaSegment::new("s", "A", dp)),
    12.0,
  )
  let kc : Array[@src.GfaTag] = [@src.GfaTag::new("KC", "i", "3")]
  assert_eq(
    @src.assembly_segment_depth(@src.GfaSegment::new("s", "A", kc)),
    3.0,
  )
  assert_eq(
    @src.assembly_segment_depth(@src.GfaSegment::new("s", "A", [])),
    1.0,
  )
}

// ============================================================================
// Tip removal
// ============================================================================

///|
test "short dead-end tip is removed" {
  let g = test_tip_graph(30)
  let result = @src.assembly_simplify(g)
  assert_eq(result.stats().tips_removed(), 1)
  assert_eq(result.stats().segments_before(), 4)
  assert_eq(result.stats().segments_after(), 3)
  assert_eq(@src.gfa_graph_n_segments(result.graph()), 3)
  assert_eq(@src.gfa_graph_n_links(result.graph()), 2)
  test_assert_segment_absent(result.graph(), "tip")
}

///|
test "long dead-end is retained above min_tip_len" {
  let g = test_tip_graph(150)
  match test_simplify_with_tip_len(g, 100) {
    Some(result) => {
      assert_eq(result.stats().tips_removed(), 0)
      assert_eq(result.stats().segments_after(), 4)
    }
    None => assert_true(false)
  }
}

///|
test "tip shorter than threshold is trimmed" {
  let g = test_tip_graph(50)
  match test_simplify_with_tip_len(g, 100) {
    Some(result) => assert_eq(result.stats().tips_removed(), 1)
    None => assert_true(false)
  }
}

///|
test "isolated short segment is removed" {
  let g = @src.GfaGraph::new()
  g.add_segment(test_segment("main", 300, -1))
  g.add_segment(test_segment("orphan", 20, -1))
  let result = @src.assembly_simplify(g)
  assert_eq(result.stats().tips_removed(), 1)
  assert_eq(@src.gfa_graph_n_segments(result.graph()), 1)
}

// ============================================================================
// Containment collapse
// ============================================================================

///|
test "contained segment is collapsed" {
  let g = @src.GfaGraph::new()
  g.add_segment(test_segment("chr", 300, -1))
  g.add_segment(test_segment("dup", 100, -1))
  g.add_containment(
    @src.GfaContainment::new("chr", "+", "dup", "+", 0, "100M", []),
  )
  let result = @src.assembly_simplify(g)
  assert_eq(result.stats().contained_removed(), 1)
  assert_eq(result.stats().segments_after(), 1)
  test_assert_segment_absent(result.graph(), "dup")
  assert_eq(@src.gfa_graph_n_links(result.graph()), 0)
}

///|
test "containment collapse keeps at least one segment" {
  let g = @src.GfaGraph::new()
  g.add_segment(test_segment("only", 100, -1))
  g.add_containment(
    @src.GfaContainment::new("only", "+", "only", "+", 0, "100M", []),
  )
  let result = @src.assembly_simplify(g)
  assert_eq(result.stats().contained_removed(), 0)
  assert_eq(result.stats().segments_after(), 1)
}

// ============================================================================
// Bubble detection and resolution
// ============================================================================

///|
test "simple bubble is detected" {
  let g = test_bubble_graph()
  let bubbles = @src.assembly_detect_bubbles(g)
  assert_eq(bubbles.length(), 1)
  assert_eq(bubbles[0].entry(), "entry")
  assert_eq(bubbles[0].exit(), "exit")
  assert_eq(bubbles[0].branch_a(), ["a1", "a2", "exit"])
  assert_eq(bubbles[0].branch_b(), ["b1", "exit"])
  assert_eq(bubbles[0].support_a(), 20.0)
  assert_eq(bubbles[0].support_b(), 2.0)
}

///|
test "bubble resolution keeps the higher-support branch" {
  let g = test_bubble_graph()
  let result = @src.assembly_simplify(g)
  assert_eq(result.stats().bubbles_detected(), 1)
  assert_eq(result.stats().bubbles_resolved(), 1)
  test_assert_segment_present(result.graph(), "a1")
  test_assert_segment_absent(result.graph(), "b1")
  assert_eq(result.stats().segments_after(), 4)
}

///|
test "equal-support bubble keeps the longer branch" {
  let g = test_equal_support_bubble_graph()
  let result = @src.assembly_simplify(g)
  assert_eq(result.stats().bubbles_detected(), 1)
  assert_eq(result.stats().bubbles_resolved(), 1)
  test_assert_segment_absent(result.graph(), "a1")
  test_assert_segment_present(result.graph(), "b1")
}

///|
test "bubble longer than max_bubble_len is ignored" {
  let g = test_bubble_graph()
  let bubbles = @src.assembly_detect_bubbles(g, max_bubble_len=150)
  assert_eq(bubbles.length(), 0)
}

///|
test "cross-linked branch is not a simple bubble" {
  let g = @src.GfaGraph::new()
  g.add_segment(test_segment("entry", 50, 1))
  g.add_segment(test_segment("s1", 100, 1))
  g.add_segment(test_segment("s2", 100, 1))
  g.add_segment(test_segment("exit", 50, 1))
  g.add_link(test_link("entry", "s1"))
  g.add_link(test_link("entry", "s2"))
  // s2 同时被 entry 与 s1 指向, 分支起点入度为 2 -> 非简单 bubble.
  g.add_link(test_link("s1", "s2"))
  g.add_link(test_link("s2", "exit"))
  let bubbles = @src.assembly_detect_bubbles(g)
  assert_eq(bubbles.length(), 0)
}

// ============================================================================
// Paths and stats output
// ============================================================================

///|
test "path referencing removed segment is dropped" {
  let g = test_tip_graph(30)
  g.add_path(@src.GfaPath::new("p_tip", ["a+", "b+", "tip+"], ["0M", "0M"], []))
  g.add_path(@src.GfaPath::new("p_main", ["a+", "b+", "c+"], ["0M", "0M"], []))
  let result = @src.assembly_simplify(g)
  assert_eq(result.stats().paths_dropped(), 1)
  assert_eq(@src.gfa_graph_n_paths(result.graph()), 1)
}

///|
test "simplify stats JSON contains all counters" {
  let g = test_tip_graph(30)
  let result = @src.assembly_simplify(g)
  let json = @src.assembly_simplify_stats_to_json(result.stats())
  assert_true(json.contains("\"tips_removed\":1"))
  assert_true(json.contains("\"segments_before\":4"))
  assert_true(json.contains("\"segments_after\":3"))
}

// ============================================================================
// Racon-lite polishing
// ============================================================================

///|
test "polishing corrects a supported substitution" {
  let reference = test_ref()
  let read = test_mutate_base(reference[0:80].to_owned(), 20)
  let reads = test_reads(read, 10)
  match test_polish_default(reference, reads) {
    Some(result) => {
      assert_eq(result.reads_used(), 10)
      assert_eq(result.reads_total(), 10)
      assert_eq(result.edits().length(), 1)
      assert_eq(result.edits()[0].kind(), "substitution")
      assert_eq(result.edits()[0].position(), 20)
      assert_eq(result.consensus()[20:21].to_owned(), "C")
      assert_true(result.mean_coverage() > 4.0)
    }
    None => assert_true(false)
  }
}

///|
test "low coverage does not change the reference" {
  let reference = test_ref()
  let read = test_mutate_base(reference[0:80].to_owned(), 20)
  let reads = test_reads(read, 2)
  match test_polish_default(reference, reads) {
    Some(result) => {
      assert_eq(result.reads_used(), 2)
      assert_eq(result.edits().length(), 0)
      assert_eq(result.consensus(), reference)
    }
    None => assert_true(false)
  }
}

///|
test "uncovered region keeps the reference sequence" {
  let reference = test_ref()
  let read = test_mutate_base(reference[0:80].to_owned(), 20)
  let reads = test_reads(read, 10)
  match test_polish_default(reference, reads) {
    Some(result) => {
      assert_eq(result.consensus().length(), reference.length())
      assert_eq(
        result.consensus()[80:160].to_owned(),
        reference[80:160].to_owned(),
      )
    }
    None => assert_true(false)
  }
}

///|
test "polishing inserts a consistently supported insertion" {
  let reference = test_ref()
  let read = reference[0:60].to_owned() + "GG" + reference[60:160].to_owned()
  let reads = test_reads(read, 10)
  match test_polish_default(reference, reads) {
    Some(result) => {
      assert_eq(result.consensus().length(), 162)
      assert_true(result.consensus().contains("GG"))
      let mut insertion_found = false
      for edit in result.edits() {
        if edit.kind() == "insertion" && edit.new_base() == "GG" {
          insertion_found = true
        }
      }
      assert_true(insertion_found)
    }
    None => assert_true(false)
  }
}

///|
test "reverse-complement reads are anchored and leave no edits" {
  let reference = test_ref()
  let read = test_revcomp(reference[10:90].to_owned())
  let reads = test_reads(read, 10)
  match test_polish_default(reference, reads) {
    Some(result) => {
      assert_eq(result.reads_used(), 10)
      assert_eq(result.edits().length(), 0)
      assert_eq(result.consensus(), reference)
    }
    None => assert_true(false)
  }
}

///|
test "noisy reads keep the consensus when errors are unsupported" {
  let reference = test_ref()
  // 少数读段携带随机错误, 多数读段保持正确 -> 不修改.
  let correct = reference[0:100].to_owned()
  let noisy = test_mutate_base(correct, 50)
  let reads = test_reads(correct, 10)
  reads.push(noisy)
  match test_polish_default(reference, reads) {
    Some(result) => {
      assert_eq(result.edits().length(), 0)
      assert_eq(result.consensus(), reference)
    }
    None => assert_true(false)
  }
}

///|
test "short reference shrinks k automatically" {
  let reference = test_ref()[0:20].to_owned()
  let reads = test_reads(reference, 5)
  match test_polish_default(reference, reads) {
    Some(result) => {
      assert_eq(result.consensus(), reference)
      assert_eq(result.reads_used(), 5)
    }
    None => assert_true(false)
  }
}

///|
test "polish config validates k" {
  let mut raised = false
  ignore(@src.polish_config(k=2)) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "polish rejects empty reference" {
  let mut raised = false
  ignore(@src.racon_polish("", [])) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "polish result JSON reports edits" {
  let reference = test_ref()
  let read = test_mutate_base(reference[0:80].to_owned(), 20)
  let reads = test_reads(read, 10)
  match test_polish_default(reference, reads) {
    Some(result) => {
      let json = @src.polish_result_to_json(result)
      assert_true(json.contains("\"kind\":\"substitution\""))
      assert_true(json.contains("\"reads_used\":10"))
    }
    None => assert_true(false)
  }
}
