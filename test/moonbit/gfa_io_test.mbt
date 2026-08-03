///|
/// Tests for the GFA I/O module (Bio.SeqIO.GfaIO equivalent).
///
/// These tests exercise the public API exposed in `src/gfa_io.mbt`:
/// the tag/header/segment/link/containment/path/graph structs, the
/// parser entry points (`gfa_parse`, `gfa_parse_*`), the serializers
/// (`gfa_*_to_line`, `gfa_to_string`), and the utility helpers
/// (`gfa_get_segment`, `gfa_get_segments_as_records`, `gfa_reverse_link`,
/// `gfa_graph_summary`, `gfa_sample_graph`).

// ============================================================================
// 1. GfaTag - creation, accessors, to_string
// ============================================================================

test "GfaTag::new and accessors" {
  let tag = @src.GfaTag::new("LN", "i", "100")
  assert_eq(tag.name(), "LN")
  assert_eq(tag.tag_type(), "i")
  assert_eq(tag.value(), "100")
}

test "GfaTag::to_string integer" {
  let tag = @src.GfaTag::new("LN", "i", "100")
  assert_eq(tag.to_string(), "LN:i:100")
}

test "GfaTag::to_string string type" {
  let tag = @src.GfaTag::new("VN", "Z", "1.0")
  assert_eq(tag.to_string(), "VN:Z:1.0")
}

test "GfaTag::to_string float type" {
  let tag = @src.GfaTag::new("GC", "f", "0.42")
  assert_eq(tag.to_string(), "GC:f:0.42")
}

test "gfa_tag_to_string matches method" {
  let tag = @src.GfaTag::new("OR", "Z", "sample")
  assert_eq(@src.gfa_tag_to_string(tag), "OR:Z:sample")
  assert_eq(@src.gfa_tag_to_string(tag), tag.to_string())
}

// ============================================================================
// 2. GfaHeader - creation and version
// ============================================================================

test "GfaHeader::new with version and tags" {
  let tags = [@src.GfaTag::new("VN", "Z", "1.0")]
  let header = @src.GfaHeader::new("1.0", tags)
  assert_eq(header.version(), "1.0")
  assert_eq(header.tags().length(), 1)
  assert_eq(header.tags()[0].name(), "VN")
  assert_eq(header.tags()[0].value(), "1.0")
}

test "GfaHeader::new empty" {
  let header = @src.GfaHeader::new("", [])
  assert_eq(header.version(), "")
  assert_eq(header.tags().length(), 0)
}

test "GfaHeader::new with multiple tags" {
  let tags = [
    @src.GfaTag::new("VN", "Z", "1.1"),
    @src.GfaTag::new("OR", "Z", "sample"),
  ]
  let header = @src.GfaHeader::new("1.1", tags)
  assert_eq(header.version(), "1.1")
  assert_eq(header.tags().length(), 2)
}

// ============================================================================
// 3. GfaSegment - creation, accessors, to_string
// ============================================================================

test "GfaSegment::new and accessors" {
  let seg = @src.GfaSegment::new(
    "seq1",
    "ACGT",
    [@src.GfaTag::new("LN", "i", "4")],
  )
  assert_eq(seg.name(), "seq1")
  assert_eq(seg.sequence(), "ACGT")
  assert_eq(seg.tags().length(), 1)
  assert_eq(seg.tags()[0].name(), "LN")
  assert_eq(seg.tags()[0].tag_type(), "i")
  assert_eq(seg.tags()[0].value(), "4")
}

test "GfaSegment::new with no tags" {
  let seg = @src.GfaSegment::new("seq2", "*", [])
  assert_eq(seg.name(), "seq2")
  assert_eq(seg.sequence(), "*")
  assert_eq(seg.tags().length(), 0)
}

test "GfaSegment::to_string with tag" {
  let seg = @src.GfaSegment::new(
    "seq1",
    "ACGT",
    [@src.GfaTag::new("LN", "i", "4")],
  )
  let s = seg.to_string()
  assert_true(s.contains("GfaSegment(name=seq1"))
  assert_true(s.contains("sequence=ACGT"))
  assert_true(s.contains("LN:i:4"))
  assert_true(s.contains("tags=["))
}

test "GfaSegment::to_string empty tags" {
  let seg = @src.GfaSegment::new("only", "ACGT", [])
  let s = seg.to_string()
  assert_true(s.contains("GfaSegment(name=only"))
  assert_true(s.contains("tags=[]"))
}

// ============================================================================
// 4. GfaLink - creation and accessors
// ============================================================================

test "GfaLink::new and accessors" {
  let link = @src.GfaLink::new(
    "seq1",
    "+",
    "seq2",
    "-",
    "8M",
    [@src.GfaTag::new("MQ", "i", "60")],
  )
  assert_eq(link.from_segment(), "seq1")
  assert_eq(link.from_orient(), "+")
  assert_eq(link.to_segment(), "seq2")
  assert_eq(link.to_orient(), "-")
  assert_eq(link.overlap(), "8M")
  assert_eq(link.tags().length(), 1)
  assert_eq(link.tags()[0].name(), "MQ")
}

test "GfaLink::new with no tags" {
  let link = @src.GfaLink::new("a", "+", "b", "+", "0M", [])
  assert_eq(link.tags().length(), 0)
  assert_eq(link.overlap(), "0M")
}

// ============================================================================
// 4b. GfaContainment - creation and accessors
// ============================================================================

test "GfaContainment::new and accessors" {
  let cont = @src.GfaContainment::new(
    "seq1",
    "+",
    "seq2",
    "+",
    5,
    "8M",
    [@src.GfaTag::new("MQ", "i", "60")],
  )
  assert_eq(cont.container(), "seq1")
  assert_eq(cont.container_orient(), "+")
  assert_eq(cont.contained(), "seq2")
  assert_eq(cont.contained_orient(), "+")
  assert_eq(cont.position(), 5)
  assert_eq(cont.overlap(), "8M")
  assert_eq(cont.tags().length(), 1)
}

// ============================================================================
// 5. GfaPath - creation, accessors
// ============================================================================

test "GfaPath::new and accessors" {
  let path = @src.GfaPath::new(
    "path1",
    ["seq1+", "seq2-"],
    ["8M", "4M"],
    [@src.GfaTag::new("NO", "i", "1")],
  )
  assert_eq(path.name(), "path1")
  assert_eq(path.segment_list().length(), 2)
  assert_eq(path.segment_list()[0], "seq1+")
  assert_eq(path.segment_list()[1], "seq2-")
  assert_eq(path.overlaps().length(), 2)
  assert_eq(path.overlaps()[0], "8M")
  assert_eq(path.overlaps()[1], "4M")
  assert_eq(path.tags().length(), 1)
}

test "GfaPath::new empty" {
  let path = @src.GfaPath::new("empty", [], [], [])
  assert_eq(path.name(), "empty")
  assert_eq(path.segment_list().length(), 0)
  assert_eq(path.overlaps().length(), 0)
  assert_eq(path.tags().length(), 0)
}

// ============================================================================
// 6. GfaGraph - creation and add operations
// ============================================================================

test "GfaGraph::new is empty" {
  let graph = @src.GfaGraph::new()
  assert_eq(graph.headers().length(), 0)
  assert_eq(graph.segments().length(), 0)
  assert_eq(graph.links().length(), 0)
  assert_eq(graph.containments().length(), 0)
  assert_eq(graph.paths().length(), 0)
}

test "GfaGraph add operations" {
  let graph = @src.GfaGraph::new()
  graph.add_header(
    @src.GfaHeader::new("1.0", [@src.GfaTag::new("VN", "Z", "1.0")]),
  )
  graph.add_segment(@src.GfaSegment::new("seq1", "ACGT", []))
  graph.add_link(@src.GfaLink::new("seq1", "+", "seq2", "+", "4M", []))
  graph.add_containment(
    @src.GfaContainment::new("seq1", "+", "seq2", "+", 5, "4M", []),
  )
  graph.add_path(@src.GfaPath::new("p1", ["seq1+"], ["4M"], []))
  assert_eq(graph.headers().length(), 1)
  assert_eq(graph.segments().length(), 1)
  assert_eq(graph.links().length(), 1)
  assert_eq(graph.containments().length(), 1)
  assert_eq(graph.paths().length(), 1)
}

// ============================================================================
// 7. gfa_parse_tag - valid and invalid inputs
// ============================================================================

test "gfa_parse_tag valid LN:i:100" {
  match @src.gfa_parse_tag("LN:i:100") {
    Some(t) => {
      assert_eq(t.name(), "LN")
      assert_eq(t.tag_type(), "i")
      assert_eq(t.value(), "100")
    }
    None => assert_true(false)
  }
}

test "gfa_parse_tag valid VN:Z:1.0" {
  match @src.gfa_parse_tag("VN:Z:1.0") {
    Some(t) => {
      assert_eq(t.name(), "VN")
      assert_eq(t.tag_type(), "Z")
      assert_eq(t.value(), "1.0")
    }
    None => assert_true(false)
  }
}

test "gfa_parse_tag invalid bad (no colons)" {
  match @src.gfa_parse_tag("bad") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

test "gfa_parse_tag missing type (one colon)" {
  match @src.gfa_parse_tag("LN:100") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

test "gfa_parse_tag empty string" {
  match @src.gfa_parse_tag("") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

test "gfa_parse_tag with value containing colon" {
  // The value itself may legally contain ':' (e.g. JSON-style). The parser
  // splits on the first two colons only.
  match @src.gfa_parse_tag("J:Z:{\"a\":1}") {
    Some(t) => {
      assert_eq(t.name(), "J")
      assert_eq(t.tag_type(), "Z")
      assert_eq(t.value(), "{\"a\":1}")
    }
    None => assert_true(false)
  }
}

// ============================================================================
// 8. gfa_parse - full GFA document with H, S, L, P lines
// ============================================================================

test "gfa_parse full document" {
  let content = "H\tVN:Z:1.0\nS\tseq1\tACGTACGT\tLN:i:8\nS\tseq2\tTTTTGGGG\tLN:i:8\nL\tseq1\t+\tseq2\t-\t8M\nP\tpath1\tseq1+,seq2-\t8M,4M\n"
  let graph = @src.gfa_parse(content)
  assert_eq(graph.headers().length(), 1)
  assert_eq(graph.headers()[0].version(), "1.0")
  assert_eq(graph.segments().length(), 2)
  assert_eq(graph.segments()[0].name(), "seq1")
  assert_eq(graph.segments()[0].sequence(), "ACGTACGT")
  assert_eq(graph.segments()[0].tags().length(), 1)
  assert_eq(graph.segments()[0].tags()[0].name(), "LN")
  assert_eq(graph.segments()[1].name(), "seq2")
  assert_eq(graph.segments()[1].sequence(), "TTTTGGGG")
  assert_eq(graph.links().length(), 1)
  assert_eq(graph.links()[0].from_segment(), "seq1")
  assert_eq(graph.links()[0].from_orient(), "+")
  assert_eq(graph.links()[0].to_segment(), "seq2")
  assert_eq(graph.links()[0].to_orient(), "-")
  assert_eq(graph.links()[0].overlap(), "8M")
  assert_eq(graph.paths().length(), 1)
  assert_eq(graph.paths()[0].name(), "path1")
  assert_eq(graph.paths()[0].segment_list().length(), 2)
  assert_eq(graph.paths()[0].segment_list()[0], "seq1+")
  assert_eq(graph.paths()[0].segment_list()[1], "seq2-")
  assert_eq(graph.paths()[0].overlaps().length(), 2)
}

test "gfa_parse handles CRLF line endings" {
  let content = "H\tVN:Z:1.0\r\nS\tseq1\tACGT\r\n"
  let graph = @src.gfa_parse(content)
  assert_eq(graph.headers().length(), 1)
  assert_eq(graph.headers()[0].version(), "1.0")
  assert_eq(graph.segments().length(), 1)
  assert_eq(graph.segments()[0].name(), "seq1")
  assert_eq(graph.segments()[0].sequence(), "ACGT")
}

test "gfa_parse skips blank and unknown lines" {
  let content = "\nH\tVN:Z:1.0\nX\tunknown\n\nS\tseq1\tACGT\n"
  let graph = @src.gfa_parse(content)
  assert_eq(graph.headers().length(), 1)
  assert_eq(graph.segments().length(), 1)
  assert_eq(graph.segments()[0].name(), "seq1")
}

test "gfa_parse_containment" {
  let fields = ["C", "seq1", "+", "seq2", "+", "5", "8M"]
  let cont = @src.gfa_parse_containment(fields)
  assert_eq(cont.container(), "seq1")
  assert_eq(cont.contained(), "seq2")
  assert_eq(cont.position(), 5)
  assert_eq(cont.overlap(), "8M")
}

// ============================================================================
// 9. gfa_to_string round-trip
// ============================================================================

test "gfa_to_string round-trip preserves data" {
  let content = "H\tVN:Z:1.0\nS\tseq1\tACGTACGT\tLN:i:8\nS\tseq2\tTTTTGGGG\tLN:i:8\nL\tseq1\t+\tseq2\t-\t8M\nP\tpath1\tseq1+,seq2-\t8M,4M\n"
  let graph1 = @src.gfa_parse(content)
  let serialized = @src.gfa_to_string(graph1)
  let graph2 = @src.gfa_parse(serialized)
  assert_eq(graph2.headers().length(), graph1.headers().length())
  assert_eq(graph2.segments().length(), graph1.segments().length())
  assert_eq(graph2.links().length(), graph1.links().length())
  assert_eq(graph2.paths().length(), graph1.paths().length())
  assert_eq(graph2.headers()[0].version(), "1.0")
  assert_eq(graph2.segments()[0].name(), "seq1")
  assert_eq(graph2.segments()[0].sequence(), "ACGTACGT")
  assert_eq(graph2.segments()[0].tags()[0].to_string(), "LN:i:8")
  assert_eq(graph2.links()[0].from_segment(), "seq1")
  assert_eq(graph2.links()[0].to_segment(), "seq2")
  assert_eq(graph2.links()[0].to_orient(), "-")
  assert_eq(graph2.paths()[0].segment_list()[0], "seq1+")
  assert_eq(graph2.paths()[0].overlaps()[1], "4M")
}

test "gfa_to_string matches expected lines" {
  let content = "H\tVN:Z:1.0\nS\tseq1\tACGT\tLN:i:4\nL\tseq1\t+\tseq2\t-\t8M\nP\tp1\tseq1+,seq2-\t8M,4M\n"
  let graph = @src.gfa_parse(content)
  let out = @src.gfa_to_string(graph)
  assert_true(out.contains("H\tVN:Z:1.0\n"))
  assert_true(out.contains("S\tseq1\tACGT\tLN:i:4\n"))
  assert_true(out.contains("L\tseq1\t+\tseq2\t-\t8M\n"))
  assert_true(out.contains("P\tp1\tseq1+,seq2-\t8M,4M\n"))
}

// ============================================================================
// Per-record serializers
// ============================================================================

test "gfa_header_to_line" {
  let header = @src.GfaHeader::new(
    "1.0",
    [@src.GfaTag::new("VN", "Z", "1.0")],
  )
  assert_eq(@src.gfa_header_to_line(header), "H\tVN:Z:1.0")
}

test "gfa_header_to_line multiple tags" {
  let header = @src.GfaHeader::new(
    "1.0",
    [
      @src.GfaTag::new("VN", "Z", "1.0"),
      @src.GfaTag::new("OR", "Z", "sample"),
    ],
  )
  assert_eq(@src.gfa_header_to_line(header), "H\tVN:Z:1.0\tOR:Z:sample")
}

test "gfa_segment_to_line" {
  let seg = @src.GfaSegment::new(
    "seq1",
    "ACGT",
    [@src.GfaTag::new("LN", "i", "4")],
  )
  assert_eq(@src.gfa_segment_to_line(seg), "S\tseq1\tACGT\tLN:i:4")
}

test "gfa_segment_to_line no tags" {
  let seg = @src.GfaSegment::new("seq2", "TTTT", [])
  assert_eq(@src.gfa_segment_to_line(seg), "S\tseq2\tTTTT")
}

test "gfa_link_to_line" {
  let link = @src.GfaLink::new("seq1", "+", "seq2", "-", "8M", [])
  assert_eq(@src.gfa_link_to_line(link), "L\tseq1\t+\tseq2\t-\t8M")
}

test "gfa_link_to_line with tag" {
  let link = @src.GfaLink::new(
    "a",
    "+",
    "b",
    "+",
    "4M",
    [@src.GfaTag::new("MQ", "i", "60")],
  )
  assert_eq(@src.gfa_link_to_line(link), "L\ta\t+\tb\t+\t4M\tMQ:i:60")
}

test "gfa_path_to_line" {
  let path = @src.GfaPath::new("p1", ["seq1+", "seq2-"], ["8M", "4M"], [])
  assert_eq(@src.gfa_path_to_line(path), "P\tp1\tseq1+,seq2-\t8M,4M")
}

test "gfa_containment_to_line" {
  let cont = @src.GfaContainment::new(
    "seq1",
    "+",
    "seq2",
    "+",
    5,
    "8M",
    [@src.GfaTag::new("MQ", "i", "60")],
  )
  assert_eq(
    @src.gfa_containment_to_line(cont),
    "C\tseq1\t+\tseq2\t+\t5\t8M\tMQ:i:60",
  )
}

// ============================================================================
// Per-record parsers
// ============================================================================

test "gfa_parse_header extracts version" {
  let fields = ["H", "VN:Z:1.0"]
  let header = @src.gfa_parse_header(fields)
  assert_eq(header.version(), "1.0")
  assert_eq(header.tags().length(), 1)
}

test "gfa_parse_segment fields" {
  let fields = ["S", "seq1", "ACGTACGT", "LN:i:8"]
  let seg = @src.gfa_parse_segment(fields)
  assert_eq(seg.name(), "seq1")
  assert_eq(seg.sequence(), "ACGTACGT")
  assert_eq(seg.tags().length(), 1)
  assert_eq(seg.tags()[0].name(), "LN")
  assert_eq(seg.tags()[0].value(), "8")
}

test "gfa_parse_segment defaults sequence to *" {
  let fields = ["S", "seq1"]
  let seg = @src.gfa_parse_segment(fields)
  assert_eq(seg.name(), "seq1")
  assert_eq(seg.sequence(), "*")
  assert_eq(seg.tags().length(), 0)
}

test "gfa_parse_link fields" {
  let fields = ["L", "seq1", "+", "seq2", "-", "8M"]
  let link = @src.gfa_parse_link(fields)
  assert_eq(link.from_segment(), "seq1")
  assert_eq(link.from_orient(), "+")
  assert_eq(link.to_segment(), "seq2")
  assert_eq(link.to_orient(), "-")
  assert_eq(link.overlap(), "8M")
}

test "gfa_parse_path fields" {
  let fields = ["P", "path1", "seq1+,seq2-", "8M,4M"]
  let path = @src.gfa_parse_path(fields)
  assert_eq(path.name(), "path1")
  assert_eq(path.segment_list().length(), 2)
  assert_eq(path.segment_list()[0], "seq1+")
  assert_eq(path.segment_list()[1], "seq2-")
  assert_eq(path.overlaps().length(), 2)
  assert_eq(path.overlaps()[0], "8M")
  assert_eq(path.overlaps()[1], "4M")
}

// ============================================================================
// 10. gfa_graph_n_segments / n_links / n_paths
// ============================================================================

test "gfa_graph_n_segments/links/paths on sample" {
  let graph = @src.gfa_sample_graph()
  assert_eq(@src.gfa_graph_n_segments(graph), 3)
  assert_eq(@src.gfa_graph_n_links(graph), 2)
  assert_eq(@src.gfa_graph_n_paths(graph), 1)
}

test "gfa_graph_n_segments/links/paths on empty" {
  let graph = @src.GfaGraph::new()
  assert_eq(@src.gfa_graph_n_segments(graph), 0)
  assert_eq(@src.gfa_graph_n_links(graph), 0)
  assert_eq(@src.gfa_graph_n_paths(graph), 0)
}

// ============================================================================
// 11. gfa_get_segment - found and not found
// ============================================================================

test "gfa_get_segment found" {
  let graph = @src.gfa_sample_graph()
  match @src.gfa_get_segment(graph, "seq2") {
    Some(seg) => {
      assert_eq(seg.name(), "seq2")
      assert_eq(seg.sequence(), "TTTTGGGG")
    }
    None => assert_true(false)
  }
}

test "gfa_get_segment found first" {
  let graph = @src.gfa_sample_graph()
  match @src.gfa_get_segment(graph, "seq1") {
    Some(seg) => assert_eq(seg.sequence(), "ACGTACGT")
    None => assert_true(false)
  }
}

test "gfa_get_segment not found" {
  let graph = @src.gfa_sample_graph()
  match @src.gfa_get_segment(graph, "missing") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

test "gfa_get_segment on empty graph" {
  let graph = @src.GfaGraph::new()
  match @src.gfa_get_segment(graph, "anything") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

// ============================================================================
// 12. gfa_get_segments_as_records - convert to SeqRecord array
// ============================================================================

test "gfa_get_segments_as_records sample graph" {
  let graph = @src.gfa_sample_graph()
  let records = @src.gfa_get_segments_as_records(graph)
  assert_eq(records.length(), 3)
  assert_eq(records[0].id, "seq1")
  assert_eq(records[0].name, "seq1")
  assert_eq(records[0].description, "seq1")
  assert_eq(records[0].seq.to_string(), "ACGTACGT")
  assert_eq(records[1].id, "seq2")
  assert_eq(records[1].seq.to_string(), "TTTTGGGG")
  assert_eq(records[2].id, "seq3")
  assert_eq(records[2].seq.to_string(), "CCCCAAAA")
}

test "gfa_get_segments_as_records empty graph" {
  let graph = @src.GfaGraph::new()
  let records = @src.gfa_get_segments_as_records(graph)
  assert_eq(records.length(), 0)
}

test "gfa_get_segments_as_records star sequence becomes empty" {
  let graph = @src.GfaGraph::new()
  graph.add_segment(@src.GfaSegment::new("masked", "*", []))
  let records = @src.gfa_get_segments_as_records(graph)
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "masked")
  assert_eq(records[0].seq.to_string(), "")
}

// ============================================================================
// 13. gfa_reverse_link - from/to swapped, orientations swapped
// ============================================================================

test "gfa_reverse_link swaps segments and orientations" {
  let link = @src.GfaLink::new(
    "seq1",
    "+",
    "seq2",
    "-",
    "8M",
    [@src.GfaTag::new("MQ", "i", "60")],
  )
  let rev = @src.gfa_reverse_link(link)
  assert_eq(rev.from_segment(), "seq2")
  assert_eq(rev.from_orient(), "-")
  assert_eq(rev.to_segment(), "seq1")
  assert_eq(rev.to_orient(), "+")
  assert_eq(rev.overlap(), "8M")
  assert_eq(rev.tags().length(), 1)
  assert_eq(rev.tags()[0].name(), "MQ")
}

test "gfa_reverse_link double-reverse is identity" {
  let link = @src.GfaLink::new("a", "+", "b", "-", "4M", [])
  let rev2 = @src.gfa_reverse_link(@src.gfa_reverse_link(link))
  assert_eq(rev2.from_segment(), "a")
  assert_eq(rev2.from_orient(), "+")
  assert_eq(rev2.to_segment(), "b")
  assert_eq(rev2.to_orient(), "-")
  assert_eq(rev2.overlap(), "4M")
}

test "gfa_reverse_link preserves tags" {
  let link = @src.GfaLink::new(
    "x",
    "+",
    "y",
    "+",
    "0M",
    [
      @src.GfaTag::new("MQ", "i", "40"),
      @src.GfaTag::new("NM", "i", "0"),
    ],
  )
  let rev = @src.gfa_reverse_link(link)
  assert_eq(rev.tags().length(), 2)
  assert_eq(rev.tags()[0].name(), "MQ")
  assert_eq(rev.tags()[1].name(), "NM")
}

// ============================================================================
// 14. gfa_graph_summary - non-empty string
// ============================================================================

test "gfa_graph_summary non-empty on sample" {
  let graph = @src.gfa_sample_graph()
  let summary = @src.gfa_graph_summary(graph)
  assert_true(summary.length() > 0)
  assert_true(summary.contains("GfaGraph"))
  assert_true(summary.contains("segments=3"))
  assert_true(summary.contains("links=2"))
  assert_true(summary.contains("paths=1"))
  assert_true(summary.contains("total_sequence_length=24"))
}

test "gfa_graph_summary on empty graph" {
  let graph = @src.GfaGraph::new()
  let summary = @src.gfa_graph_summary(graph)
  assert_true(summary.contains("headers=0"))
  assert_true(summary.contains("segments=0"))
  assert_true(summary.contains("links=0"))
  assert_true(summary.contains("paths=0"))
  assert_true(summary.contains("total_sequence_length=0"))
}

test "gfa_graph_summary ignores star sequences" {
  let graph = @src.GfaGraph::new()
  graph.add_segment(@src.GfaSegment::new("a", "ACGT", []))
  graph.add_segment(@src.GfaSegment::new("b", "*", []))
  let summary = @src.gfa_graph_summary(graph)
  assert_true(summary.contains("segments=2"))
  assert_true(summary.contains("total_sequence_length=4"))
}

// ============================================================================
// 15. Sample graph validation
// ============================================================================

test "gfa_sample_graph structure" {
  let graph = @src.gfa_sample_graph()
  assert_eq(graph.headers().length(), 1)
  assert_eq(graph.headers()[0].version(), "1.0")
  assert_eq(graph.segments().length(), 3)
  assert_eq(graph.links().length(), 2)
  assert_eq(graph.paths().length(), 1)
  assert_eq(graph.containments().length(), 0)
}

test "gfa_sample_graph segments" {
  let graph = @src.gfa_sample_graph()
  assert_eq(graph.segments()[0].name(), "seq1")
  assert_eq(graph.segments()[0].sequence(), "ACGTACGT")
  assert_eq(graph.segments()[1].name(), "seq2")
  assert_eq(graph.segments()[1].sequence(), "TTTTGGGG")
  assert_eq(graph.segments()[2].name(), "seq3")
  assert_eq(graph.segments()[2].sequence(), "CCCCAAAA")
  // Each segment carries an LN tag matching its sequence length.
  assert_eq(graph.segments()[0].tags()[0].name(), "LN")
  assert_eq(graph.segments()[0].tags()[0].value(), "8")
}

test "gfa_sample_graph links" {
  let graph = @src.gfa_sample_graph()
  let l0 = graph.links()[0]
  assert_eq(l0.from_segment(), "seq1")
  assert_eq(l0.to_segment(), "seq2")
  assert_eq(l0.overlap(), "4M")
  let l1 = graph.links()[1]
  assert_eq(l1.from_segment(), "seq2")
  assert_eq(l1.to_segment(), "seq3")
  assert_eq(l1.overlap(), "4M")
}

test "gfa_sample_graph path" {
  let graph = @src.gfa_sample_graph()
  let p = graph.paths()[0]
  assert_eq(p.name(), "path1")
  assert_eq(p.segment_list().length(), 3)
  assert_eq(p.segment_list()[0], "seq1+")
  assert_eq(p.segment_list()[1], "seq2-")
  assert_eq(p.segment_list()[2], "seq3+")
  assert_eq(p.overlaps().length(), 2)
}

test "gfa_sample_graph round-trips through gfa_to_string" {
  let graph = @src.gfa_sample_graph()
  let text = @src.gfa_to_string(graph)
  let reparsed = @src.gfa_parse(text)
  assert_eq(@src.gfa_graph_n_segments(reparsed), 3)
  assert_eq(@src.gfa_graph_n_links(reparsed), 2)
  assert_eq(@src.gfa_graph_n_paths(reparsed), 1)
  assert_eq(reparsed.headers()[0].version(), "1.0")
}

// ============================================================================
// 16. Edge cases - empty graph, single segment, no links
// ============================================================================

test "edge case: empty graph serializes to empty string" {
  let graph = @src.GfaGraph::new()
  let s = @src.gfa_to_string(graph)
  assert_eq(s, "")
}

test "edge case: parse empty content" {
  let graph = @src.gfa_parse("")
  assert_eq(@src.gfa_graph_n_segments(graph), 0)
  assert_eq(@src.gfa_graph_n_links(graph), 0)
  assert_eq(@src.gfa_graph_n_paths(graph), 0)
}

test "edge case: single segment, no links" {
  let graph = @src.GfaGraph::new()
  graph.add_segment(@src.GfaSegment::new("only", "ACGT", []))
  assert_eq(@src.gfa_graph_n_segments(graph), 1)
  assert_eq(@src.gfa_graph_n_links(graph), 0)
  assert_eq(@src.gfa_graph_n_paths(graph), 0)
  let out = @src.gfa_to_string(graph)
  assert_eq(out, "S\tonly\tACGT\n")
}

test "edge case: single segment round-trip" {
  let graph = @src.GfaGraph::new()
  graph.add_segment(@src.GfaSegment::new("only", "ACGT", []))
  let text = @src.gfa_to_string(graph)
  let reparsed = @src.gfa_parse(text)
  assert_eq(@src.gfa_graph_n_segments(reparsed), 1)
  assert_eq(reparsed.segments()[0].name(), "only")
  assert_eq(reparsed.segments()[0].sequence(), "ACGT")
}

test "edge case: segment with star sequence" {
  let content = "S\tseq1\t*\tLN:i:100\n"
  let graph = @src.gfa_parse(content)
  assert_eq(@src.gfa_graph_n_segments(graph), 1)
  assert_eq(graph.segments()[0].sequence(), "*")
  assert_eq(graph.segments()[0].tags()[0].name(), "LN")
}

test "edge case: only headers" {
  let content = "H\tVN:Z:1.0\nH\tOR:Z:test\n"
  let graph = @src.gfa_parse(content)
  assert_eq(graph.headers().length(), 2)
  assert_eq(@src.gfa_graph_n_segments(graph), 0)
}

test "edge case: unknown record type ignored" {
  let content = "H\tVN:Z:1.0\nZ\tunknown\tfield\nS\tseq1\tACGT\n"
  let graph = @src.gfa_parse(content)
  assert_eq(graph.headers().length(), 1)
  assert_eq(@src.gfa_graph_n_segments(graph), 1)
}
