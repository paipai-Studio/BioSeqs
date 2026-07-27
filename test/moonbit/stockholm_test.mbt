///|
/// Tests for Stockholm format parser module.

test "stockholm_parse_header" {
  let content = "# STOCKHOLM 1.0\nseq1/1-10 ACGTACGTAC\n//\n"
  let ali = @src.stockholm_parse(content)
  assert_eq(ali.version, "1.0")
  assert_eq(ali.blocks.length(), 1)
}

test "stockholm_parse_empty" {
  let ali = @src.stockholm_parse("")
  assert_eq(ali.blocks.length(), 0)
  assert_eq(ali.gf_annotations.keys().length(), 0)
}

test "stockholm_parse_header_only" {
  let content = "# STOCKHOLM 1.0\n//\n"
  let ali = @src.stockholm_parse(content)
  assert_eq(ali.blocks.length(), 0)
}

test "stockholm_parse_single_sequence" {
  let content = "# STOCKHOLM 1.0\nseq1/1-10 ACGTACGTAC\n//\n"
  let ali = @src.stockholm_parse(content)
  assert_eq(ali.blocks.length(), 1)
  assert_eq(ali.blocks[0].sequences.length(), 1)
  assert_eq(ali.blocks[0].sequences[0].name, "seq1")
  assert_eq(ali.blocks[0].sequences[0].aligned_seq, "ACGTACGTAC")
  assert_eq(ali.blocks[0].sequences[0].start, 1)
  assert_eq(ali.blocks[0].sequences[0].end, 10)
}

test "stockholm_parse_multiple_sequences" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "seq2/1-10 ACGTACGTAC\n" +
    "seq3/1-10 ACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  assert_eq(ali.blocks.length(), 1)
  assert_eq(ali.blocks[0].sequences.length(), 3)
}

test "stockholm_parse_gf_annotations" {
  let content = "# STOCKHOLM 1.0\n" +
    "#=GF ID   TestFamily\n" +
    "#=GF AC   PF00001.20\n" +
    "#=GF DE   A test protein family\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  assert_eq(ali.gf_annotations.keys().length(), 3)
  assert_eq(ali.gf_annotations["ID"], "TestFamily")
  assert_eq(ali.gf_annotations["AC"], "PF00001.20")
  assert_eq(ali.gf_annotations["DE"], "A test protein family")
}

test "stockholm_parse_gc_ss_cons" {
  let content = "# STOCKHOLM 1.0\n" +
    "#=GC SS_cons    ........((((...)))).....\n" +
    "seq1/1-22 ACGTACGTACGTACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  assert_eq(ali.blocks[0].secondary_structure, "........((((...)))).....")
}

test "stockholm_parse_gr_annotation" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "#=GR seq1    ........\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  assert_true(ali.blocks[0].gr_annotations.keys().length() > 0)
}

test "stockholm_parse_gs_annotation" {
  let content = "# STOCKHOLM 1.0\n" +
    "#=GS seq1 some annotation value\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  assert_true(ali.gs_annotations.keys().length() > 0)
}

test "stockholm_parse_multiple_blocks" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "seq2/1-10 ACGTACGTAC\n" +
    "\n" +
    "seq1/11-20 ACGTACGTAC\n" +
    "seq2/11-20 ACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  assert_eq(ali.blocks.length(), 2)
  assert_eq(ali.blocks[0].sequences.length(), 2)
  assert_eq(ali.blocks[1].sequences.length(), 2)
}

test "stockholm_sequence_struct" {
  let seq = @src.StockholmSequence::new(
    "test",
    "PF00001",
    "Test description",
    "ACGTAC",
    1,
    6,
  )
  assert_eq(seq.name, "test")
  assert_eq(seq.accession, "PF00001")
  assert_eq(seq.description, "Test description")
  assert_eq(seq.aligned_seq, "ACGTAC")
  assert_eq(seq.start, 1)
  assert_eq(seq.end, 6)
}

test "stockholm_block_struct" {
  let seqs : Array[@src.StockholmSequence] = Array::new()
  let gc : Map[String, String] = Map([], capacity=0)
  let gr : Map[String, String] = Map([], capacity=0)
  let ann : Array[String] = Array::new()
  let block = @src.StockholmBlock::new(seqs, "ss", "cons", ann, gc, gr)
  assert_eq(block.secondary_structure, "ss")
  assert_eq(block.consensus, "cons")
}

test "stockholm_alignment_struct" {
  let blocks : Array[@src.StockholmBlock] = Array::new()
  let gf : Map[String, String] = Map([], capacity=0)
  let gs : Map[String, String] = Map([], capacity=0)
  let markup : Array[String] = Array::new()
  let ali = @src.StockholmAlignment::new("1.0", blocks, gf, gs, markup)
  assert_eq(ali.version, "1.0")
  assert_eq(ali.blocks.length(), 0)
}

test "stockholm_percent_identity_identical" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "seq2/1-10 ACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  let pid = @src.stockholm_percent_identity(ali.blocks[0])
  assert_true(pid > 0.99)
}

test "stockholm_percent_identity_different" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "seq2/1-10 TTTTTTTTTT\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  let pid = @src.stockholm_percent_identity(ali.blocks[0])
  assert_true(pid < 0.3)
}

test "stockholm_percent_identity_with_gaps" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "seq2/1-10 ACGT--GTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  let pid = @src.stockholm_percent_identity(ali.blocks[0])
  assert_true(pid >= 0.0)
  assert_true(pid <= 1.0)
}

test "stockholm_conservation_single_column" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1/1-5 ACGT.\n" +
    "seq2/1-5 ACGT.\n" +
    "seq3/1-5 ACGT-\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  let cons = @src.stockholm_conservation(ali.blocks[0])
  assert_eq(cons.length(), 5)
  assert_true(cons[0] > 0.99)
  assert_true(cons[4] >= 0.0)
}

test "stockholm_conservation_all_gap_column" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1/1-5 ACGT-\n" +
    "seq2/1-5 ACGT-\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  let cons = @src.stockholm_conservation(ali.blocks[0])
  assert_eq(cons[4], 0.0)
}

test "stockholm_to_fasta_basic" {
  let content = "# STOCKHOLM 1.0\n" +
    "#=GF AC   PF00001\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "seq2/1-10 ACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  let fasta = @src.stockholm_to_fasta(ali)
  assert_true(fasta.contains(">seq1"))
  assert_true(fasta.contains(">seq2"))
  assert_true(fasta.contains("ACGTACGTAC"))
}

test "stockholm_to_fasta_no_duplicates" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "\n" +
    "seq1/11-20 ACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  let fasta = @src.stockholm_to_fasta(ali)
  assert_true(fasta.contains(">seq1"))
}

test "stockholm_merge_blocks_basic" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "\n" +
    "seq1/11-20 ACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  assert_eq(ali.blocks.length(), 2)
  let merged = @src.stockholm_merge_blocks(ali)
  assert_eq(merged.blocks.length(), 1)
  assert_eq(merged.blocks[0].sequences.length(), 1)
}

test "stockholm_merge_blocks_sequence_concat" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1/1-5 ACGTG\n" +
    "\n" +
    "seq1/6-10 TACGT\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  let merged = @src.stockholm_merge_blocks(ali)
  let seq = merged.blocks[0].sequences[0]
  assert_eq(seq.aligned_seq, "ACGTGTACGT")
}

test "stockholm_write_roundtrip" {
  let content = "# STOCKHOLM 1.0\n" +
    "#=GF ID   TestFamily\n" +
    "#=GF AC   PF00001.20\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "seq2/1-10 ACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  let output = @src.stockholm_write(ali)
  assert_true(output.contains("# STOCKHOLM"))
  assert_true(output.contains("seq1"))
  assert_true(output.contains("seq2"))
  assert_true(output.contains("//"))
}

test "stockholm_write_contains_gf" {
  let content = "# STOCKHOLM 1.0\n" +
    "#=GF ID   TestFamily\n" +
    "#=GF AC   PF00001\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  let output = @src.stockholm_write(ali)
  assert_true(output.contains("#=GF ID"))
  assert_true(output.contains("#=GF AC"))
}

test "stockholm_sequence_without_range" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1 ACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  assert_eq(ali.blocks[0].sequences.length(), 1)
  assert_eq(ali.blocks[0].sequences[0].name, "seq1")
  assert_eq(ali.blocks[0].sequences[0].start, 1)
}

test "stockholm_gap_handling" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1/1-10 ACGT--GTAC\n" +
    "seq2/1-10 ACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  assert_eq(ali.blocks[0].sequences[0].aligned_seq, "ACGT--GTAC")
  assert_eq(ali.blocks[0].sequences[1].aligned_seq, "ACGTACGTAC")
}

test "stockholm_sample_content" {
  let content = @src.sample_stockholm_content()
  let ali = @src.stockholm_parse(content)
  assert_eq(ali.blocks.length(), 1)
  assert_eq(ali.gf_annotations.keys().length(), 3)
  assert_eq(ali.blocks[0].sequences.length(), 3)
}

test "stockholm_sample_alignment" {
  let ali = @src.sample_stockholm_alignment()
  assert_eq(ali.version, "1.0")
  assert_eq(ali.blocks.length(), 1)
  assert_true(ali.blocks[0].secondary_structure.length() > 0)
}

test "stockholm_multiline_gf" {
  let content = "# STOCKHOLM 1.0\n" +
    "#=GF DE   This is a long description\n" +
    "#=GF DE   that spans multiple lines\n" +
    "seq1/1-5 ACGTG\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  assert_true(ali.gf_annotations.keys().length() > 0)
}

test "stockholm_conservation_with_mixed_chars" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1/1-5 ACGTN\n" +
    "seq2/1-5 ACGTN\n" +
    "seq3/1-5 ACGTA\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  let cons = @src.stockholm_conservation(ali.blocks[0])
  assert_eq(cons.length(), 5)
  assert_true(cons[3] > 0.99)
  assert_true(cons[4] >= 0.66)
}

test "stockholm_percent_identity_single_sequence" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  let pid = @src.stockholm_percent_identity(ali.blocks[0])
  assert_eq(pid, 0.0)
}

test "stockholm_to_fasta_empty_alignment" {
  let content = "# STOCKHOLM 1.0\n//\n"
  let ali = @src.stockholm_parse(content)
  let fasta = @src.stockholm_to_fasta(ali)
  assert_eq(fasta.length(), 0)
}

test "stockholm_merge_single_block" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  let merged = @src.stockholm_merge_blocks(ali)
  assert_eq(merged.blocks.length(), 1)
  assert_eq(merged.blocks[0].sequences.length(), 1)
}

test "stockholm_parse_with_comment_lines" {
  let content = "# STOCKHOLM 1.0\n" +
    "# Some comment\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  assert_eq(ali.blocks[0].sequences.length(), 1)
  assert_true(ali.markup.length() > 0)
}

test "stockholm_ss_cons_roundtrip" {
  let content = "# STOCKHOLM 1.0\n" +
    "#=GC SS_cons    ........((((...)))).....\n" +
    "seq1/1-22 ACGTACGTACGTACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  let output = @src.stockholm_write(ali)
  assert_true(output.contains("SS_cons"))
  assert_true(output.contains("........((((...))))....."))
}

test "stockholm_parse_range_with_dash" {
  let content = "# STOCKHOLM 1.0\n" +
    "seq1/1-100 ACGTACGT\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  assert_eq(ali.blocks[0].sequences[0].start, 1)
  assert_eq(ali.blocks[0].sequences[0].end, 100)
}

test "stockholm_parse_accession_from_gf" {
  let content = "# STOCKHOLM 1.0\n" +
    "#=GF AC   PF12345.6\n" +
    "#=GF DE   Test family\n" +
    "seq1/1-10 ACGTACGTAC\n" +
    "//\n"
  let ali = @src.stockholm_parse(content)
  assert_eq(ali.gf_annotations["AC"], "PF12345.6")
  assert_eq(ali.gf_annotations["DE"], "Test family")
}