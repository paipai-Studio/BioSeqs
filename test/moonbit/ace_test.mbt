///|
/// Tests for ACE contig assembly format parser module.

test "ace_read_new" {
  let read = @src.AceRead::new(
    "read1",
    "ATCGATCG",
    [30, 30, 30, 30, 30, 30, 30, 30],
    1,
    8,
    "+",
    "",
    "chem1",
    "dye1",
  )
  assert_eq(read.read_id, "read1")
  assert_eq(read.sequence, "ATCGATCG")
  assert_eq(read.quality.length(), 8)
  assert_eq(read.clip_start, 1)
  assert_eq(read.clip_end, 8)
  assert_eq(read.sequence_length, 8)
  assert_eq(read.strand, "+")
  assert_eq(read.chemistry, "chem1")
  assert_eq(read.dye, "dye1")
}

test "ace_contig_new" {
  let read = @src.AceRead::new(
    "read1",
    "ATCGATCG",
    [30, 30, 30, 30, 30, 30, 30, 30],
    1,
    8,
    "+",
    "",
    "",
    "",
  )
  let contig = @src.AceContig::new(
    "contig1",
    "ATCGATCG",
    [read],
    [40, 40, 40, 40, 40, 40, 40, 40],
    "",
    false,
  )
  assert_eq(contig.contig_name, "contig1")
  assert_eq(contig.sequence, "ATCGATCG")
  assert_eq(contig.length, 8)
  assert_eq(contig.reads.length(), 1)
  assert_eq(contig.base_qualities.length(), 8)
  assert_false(contig.reverse_complemented)
}

test "ace_alignment_new" {
  let qual : Array[Int] = [30, 30, 30]
  let ali = @src.AceAlignment::new("contig1", 10, "ATC", qual)
  assert_eq(ali.contig_id, "contig1")
  assert_eq(ali.start_position, 10)
  assert_eq(ali.sequence, "ATC")
  assert_eq(ali.quality.length(), 3)
}

test "ace_data_new" {
  let contigs : Array[@src.AceContig] = []
  let data = @src.AceData::new(contigs, 0)
  assert_eq(data.total_contigs, 0)
  assert_eq(data.total_reads, 0)
  assert_eq(data.base_qual_threshold, 0)
}

test "ace_parse_empty" {
  let data = @src.ace_parse("")
  assert_eq(data.contigs.length(), 0)
  assert_eq(data.total_contigs, 0)
}

test "ace_parse_only_header" {
  let data = @src.ace_parse("AF contig1 100\n")
  assert_eq(data.contigs.length(), 0)
  assert_eq(data.total_contigs, 0)
}

test "ace_parse_single_read" {
  let content = "AF contig1 10\n" +
    "RD read1 10 chem1 dye1\n" +
    "ATCGATCGAT\n" +
    "q 30 30 30 30 30 30 30 30 30 30\n" +
    "CT contig1 10 0 10\n" +
    "ATCGATCGAT\n" +
    "q 40 40 40 40 40 40 40 40 40 40\n" +
    "AF read1 1 10 1 10 1 10 +\n"
  let data = @src.ace_parse(content)
  assert_eq(data.contigs.length(), 1)
  assert_eq(data.contigs[0].contig_name, "contig1")
  assert_eq(data.contigs[0].length, 10)
  assert_eq(data.contigs[0].sequence, "ATCGATCGAT")
  assert_eq(data.contigs[0].reads.length(), 1)
  assert_eq(data.contigs[0].reads[0].read_id, "read1")
}

test "ace_parse_multiple_reads" {
  let content = "AF contig1 20\n" +
    "RD read1 10 chem1 dye1\n" +
    "ATCGATCGAT\n" +
    "q 30 30 30 30 30 30 30 30 30 30\n" +
    "RD read2 10 chem2 dye2\n" +
    "GCTAGCTAGC\n" +
    "q 25 25 25 25 25 25 25 25 25 25\n" +
    "CT contig1 20 0 20\n" +
    "ATCGATCGATGCTAGCTAGC\n" +
    "q 40 40 40 40 40 40 40 40 40 40 35 35 35 35 35 35 35 35 35 35\n" +
    "AF read1 1 10 1 10 1 10 +\n" +
    "AF read2 11 20 1 10 1 10 -\n"
  let data = @src.ace_parse(content)
  assert_eq(data.contigs.length(), 1)
  assert_eq(data.contigs[0].reads.length(), 2)
}

test "ace_parse_multiple_contigs" {
  let content = "AF contig1 10\n" +
    "RD read1 10 chem1 dye1\n" +
    "AAAAAAAAAA\n" +
    "q 30 30 30 30 30 30 30 30 30 30\n" +
    "CT contig1 10 0 10\n" +
    "AAAAAAAAAA\n" +
    "q 40 40 40 40 40 40 40 40 40 40\n" +
    "AF read1 1 10 1 10 1 10 +\n" +
    "AF contig2 10\n" +
    "RD read2 10 chem2 dye2\n" +
    "CCCCCCCCCC\n" +
    "q 25 25 25 25 25 25 25 25 25 25\n" +
    "CT contig2 10 0 10\n" +
    "CCCCCCCCCC\n" +
    "q 35 35 35 35 35 35 35 35 35 35\n" +
    "AF read2 1 10 1 10 1 10 +\n"
  let data = @src.ace_parse(content)
  assert_eq(data.contigs.length(), 2)
  assert_eq(data.total_contigs, 2)
  assert_eq(data.total_reads, 2)
}

test "ace_parse_strands" {
  let content = "AF contig1 10\n" +
    "RD read1 10 chem1 dye1\n" +
    "AAAAAAAAAA\n" +
    "q 30 30 30 30 30 30 30 30 30 30\n" +
    "RD read2 10 chem2 dye2\n" +
    "CCCCCCCCCC\n" +
    "q 25 25 25 25 25 25 25 25 25 25\n" +
    "CT contig1 10 0 10\n" +
    "AAAAAAAAAA\n" +
    "q 40 40 40 40 40 40 40 40 40 40\n" +
    "AF read1 1 10 1 10 1 10 +\n" +
    "AF read2 1 10 1 10 1 10 -\n"
  let data = @src.ace_parse(content)
  assert_eq(data.contigs[0].reads.length(), 2)
  assert_eq(data.contigs[0].reads[0].strand, "+")
  assert_eq(data.contigs[0].reads[1].strand, "-")
}

test "ace_parse_quality_scores" {
  let content = "AF contig1 5\n" +
    "RD read1 5 chem1 dye1\n" +
    "ACGTG\n" +
    "q 10 20 30 40 50\n" +
    "CT contig1 5 0 5\n" +
    "ACGTG\n" +
    "q 15 25 35 45 55\n" +
    "AF read1 1 5 1 5 1 5 +\n"
  let data = @src.ace_parse(content)
  assert_eq(data.contigs[0].reads[0].quality.length(), 5)
  assert_eq(data.contigs[0].reads[0].quality[0], 10)
  assert_eq(data.contigs[0].reads[0].quality[4], 50)
  assert_eq(data.contigs[0].base_qualities.length(), 5)
  assert_eq(data.contigs[0].base_qualities[0], 15)
  assert_eq(data.contigs[0].base_qualities[4], 55)
}

test "ace_parse_reads_only" {
  let content = "RD read1 5 chem1 dye1\n" +
    "ACGTG\n" +
    "q 10 20 30 40 50\n" +
    "RD read2 5 chem2 dye2\n" +
    "TGCAC\n" +
    "q 15 25 35 45 55\n"
  let reads = @src.ace_parse_reads(content)
  assert_eq(reads.length(), 2)
  assert_eq(reads[0].read_id, "read1")
  assert_eq(reads[0].sequence, "ACGTG")
  assert_eq(reads[1].read_id, "read2")
}

test "ace_parse_contigs_only" {
  let content = "RD read1 5 chem1 dye1\n" +
    "ACGTG\n" +
    "q 10 20 30 40 50\n" +
    "CT contig1 5 0 5\n" +
    "ACGTG\n" +
    "q 15 25 35 45 55\n" +
    "AF read1 1 5 1 5 1 5 +\n"
  let contigs = @src.ace_parse_contigs(content)
  assert_eq(contigs.length(), 1)
  assert_eq(contigs[0].contig_name, "contig1")
  assert_eq(contigs[0].reads.length(), 1)
}

test "ace_get_contig" {
  let content = "AF contig1 10\n" +
    "RD read1 10 chem1 dye1\n" +
    "AAAAAAAAAA\n" +
    "q 30 30 30 30 30 30 30 30 30 30\n" +
    "CT contig1 10 0 10\n" +
    "AAAAAAAAAA\n" +
    "q 40 40 40 40 40 40 40 40 40 40\n" +
    "AF read1 1 10 1 10 1 10 +\n"
  let data = @src.ace_parse(content)
  let contig = @src.ace_get_contig(data, "contig1")
  match contig {
    Some(c) => assert_eq(c.contig_name, "contig1")
    None => assert_true(false)
  }
  let not_found = @src.ace_get_contig(data, "nonexistent")
  match not_found {
    None => assert_true(true)
    Some(_) => assert_true(false)
  }
}

test "ace_contig_length" {
  let content = "AF contig1 10\n" +
    "RD read1 10 chem1 dye1\n" +
    "AAAAAAAAAA\n" +
    "q 30 30 30 30 30 30 30 30 30 30\n" +
    "CT contig1 10 0 10\n" +
    "AAAAAAAAAA\n" +
    "q 40 40 40 40 40 40 40 40 40 40\n" +
    "AF read1 1 10 1 10 1 10 +\n"
  let data = @src.ace_parse(content)
  assert_eq(@src.ace_contig_length(data, "contig1"), 10)
  assert_eq(@src.ace_contig_length(data, "nonexistent"), 0)
}

test "ace_contig_reads" {
  let content = "AF contig1 10\n" +
    "RD read1 10 chem1 dye1\n" +
    "AAAAAAAAAA\n" +
    "q 30 30 30 30 30 30 30 30 30 30\n" +
    "CT contig1 10 0 10\n" +
    "AAAAAAAAAA\n" +
    "q 40 40 40 40 40 40 40 40 40 40\n" +
    "AF read1 1 10 1 10 1 10 +\n"
  let data = @src.ace_parse(content)
  let reads = @src.ace_contig_reads(data, "contig1")
  assert_eq(reads.length(), 1)
  assert_eq(reads[0].read_id, "read1")
  let empty_reads = @src.ace_contig_reads(data, "nonexistent")
  assert_eq(empty_reads.length(), 0)
}

test "ace_read_coverage" {
  let read1 = @src.AceRead::new(
    "read1",
    "AAAAAAAAAA",
    [30, 30, 30, 30, 30, 30, 30, 30, 30, 30],
    1,
    10,
    "+",
    "",
    "",
    "",
  )
  let read2 = @src.AceRead::new(
    "read2",
    "CCCCCCCCCC",
    [25, 25, 25, 25, 25, 25, 25, 25, 25, 25],
    1,
    10,
    "-",
    "",
    "",
    "",
  )
  let contig = @src.AceContig::new(
    "contig1",
    "ACGTACGTAC",
    [read1, read2],
    [40, 40, 40, 40, 40, 40, 40, 40, 40, 40],
    "",
    false,
  )
  let coverage = @src.ace_read_coverage(contig)
  assert_eq(coverage.get(1), Some(2))
  assert_eq(coverage.get(5), Some(2))
  assert_eq(coverage.get(10), Some(2))
}

test "ace_read_coverage_partial" {
  let read1 = @src.AceRead::new(
    "read1",
    "AAAAAA",
    [30, 30, 30, 30, 30, 30],
    1,
    6,
    "+",
    "",
    "",
    "",
  )
  let read2 = @src.AceRead::new(
    "read2",
    "CCCCCC",
    [25, 25, 25, 25, 25, 25],
    4,
    9,
    "-",
    "",
    "",
    "",
  )
  let contig = @src.AceContig::new(
    "contig1",
    "ACGTACGTACGTAC",
    [read1, read2],
    [40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40],
    "",
    false,
  )
  let coverage = @src.ace_read_coverage(contig)
  assert_eq(coverage.get(1), Some(1))
  assert_eq(coverage.get(2), Some(1))
  assert_eq(coverage.get(3), Some(1))
  assert_eq(coverage.get(4), Some(2))
  assert_eq(coverage.get(5), Some(2))
  assert_eq(coverage.get(6), Some(2))
  assert_eq(coverage.get(7), Some(1))
  assert_eq(coverage.get(8), Some(1))
  assert_eq(coverage.get(9), Some(1))
  assert_eq(coverage.get(10), Some(0))
}

test "ace_consensus_sequence" {
  let read = @src.AceRead::new(
    "read1",
    "ATCG",
    [30, 30, 30, 30],
    1,
    4,
    "+",
    "",
    "",
    "",
  )
  let contig = @src.AceContig::new(
    "contig1",
    "GCTA",
    [read],
    [40, 40, 40, 40],
    "",
    false,
  )
  assert_eq(@src.ace_consensus_sequence(contig), "GCTA")
}

test "ace_contig_gc_content" {
  let read = @src.AceRead::new(
    "read1",
    "AATT",
    [30, 30, 30, 30],
    1,
    4,
    "+",
    "",
    "",
    "",
  )
  let contig = @src.AceContig::new(
    "contig1",
    "GCGC",
    [read],
    [40, 40, 40, 40],
    "",
    false,
  )
  let gc = @src.ace_contig_gc_content(contig)
  assert_eq(gc, 1.0)
}

test "ace_contig_gc_content_mixed" {
  let read = @src.AceRead::new(
    "read1",
    "ACGT",
    [30, 30, 30, 30],
    1,
    4,
    "+",
    "",
    "",
    "",
  )
  let contig = @src.AceContig::new(
    "contig1",
    "ACGT",
    [read],
    [40, 40, 40, 40],
    "",
    false,
  )
  let gc = @src.ace_contig_gc_content(contig)
  assert_eq(gc, 0.5)
}

test "ace_contig_gc_content_empty" {
  let contig = @src.AceContig::new(
    "contig1",
    "",
    [],
    [],
    "",
    false,
  )
  let gc = @src.ace_contig_gc_content(contig)
  assert_eq(gc, 0.0)
}

test "ace_to_string_basic" {
  let read = @src.AceRead::new(
    "read1",
    "ATCG",
    [30, 30, 30, 30],
    1,
    4,
    "+",
    "",
    "chem1",
    "dye1",
  )
  let contig = @src.AceContig::new(
    "contig1",
    "ATCG",
    [read],
    [40, 40, 40, 40],
    "",
    false,
  )
  let data = @src.AceData::new([contig], 0)
  let output = @src.ace_to_string(data)
  assert_true(output.contains("AF contig1"))
  assert_true(output.contains("RD read1"))
  assert_true(output.contains("CT contig1"))
  assert_true(output.contains("AF read1"))
}

test "ace_to_string_roundtrip" {
  let content = "AF contig1 10\n" +
    "RD read1 10 chem1 dye1\n" +
    "AAAAAAAAAA\n" +
    "q 30 30 30 30 30 30 30 30 30 30\n" +
    "CT contig1 10 0 10\n" +
    "AAAAAAAAAA\n" +
    "q 40 40 40 40 40 40 40 40 40 40\n" +
    "AF read1 1 10 1 10 1 10 +\n"
  let data = @src.ace_parse(content)
  let output = @src.ace_to_string(data)
  let data2 = @src.ace_parse(output)
  assert_eq(data2.contigs.length(), 1)
  assert_eq(data2.contigs[0].contig_name, "contig1")
  assert_eq(data2.contigs[0].length, 10)
  assert_eq(data2.contigs[0].reads.length(), 1)
  assert_eq(data2.contigs[0].reads[0].read_id, "read1")
}

test "ace_sample_data" {
  let data = @src.sample_ace_data()
  assert_eq(data.contigs.length(), 1)
  assert_eq(data.total_reads, 2)
  assert_eq(data.total_contigs, 1)
  assert_eq(data.contigs[0].contig_name, "contig1")
  assert_eq(data.contigs[0].length, 60)
  assert_eq(data.contigs[0].reads.length(), 2)
}

test "ace_multiline_sequence" {
  let mut seq = ""
  let mut i = 0
  while i < 120 {
    seq = seq + "A"
    i = i + 1
  }
  let mut qual_str = ""
  let mut qi = 0
  while qi < 120 {
    if qi > 0 {
      qual_str = qual_str + " "
    }
    qual_str = qual_str + "30"
    qi = qi + 1
  }
  let content = "AF contig1 120\n" +
    "RD read1 120 chem1 dye1\n" +
    seq + "\n" +
    "q " + qual_str + "\n" +
    "CT contig1 120 0 120\n" +
    seq + "\n" +
    "q " + qual_str + "\n" +
    "AF read1 1 120 1 120 1 120 +\n"
  let data = @src.ace_parse(content)
  assert_eq(data.contigs.length(), 1)
  assert_eq(data.contigs[0].length, 120)
  assert_eq(data.contigs[0].sequence.length(), 120)
  assert_eq(data.contigs[0].reads[0].sequence.length(), 120)
}

test "ace_total_reads_count" {
  let content = "AF contig1 10\n" +
    "RD read1 10 chem1 dye1\n" +
    "AAAAAAAAAA\n" +
    "q 30 30 30 30 30 30 30 30 30 30\n" +
    "CT contig1 10 0 10\n" +
    "AAAAAAAAAA\n" +
    "q 40 40 40 40 40 40 40 40 40 40\n" +
    "AF read1 1 10 1 10 1 10 +\n" +
    "AF contig2 10\n" +
    "RD read2 10 chem2 dye2\n" +
    "CCCCCCCCCC\n" +
    "q 25 25 25 25 25 25 25 25 25 25\n" +
    "CT contig2 10 0 10\n" +
    "CCCCCCCCCC\n" +
    "q 35 35 35 35 35 35 35 35 35 35\n" +
    "AF read2 1 10 1 10 1 10 +\n"
  let data = @src.ace_parse(content)
  assert_eq(data.total_reads, 2)
}

test "ace_base_qualities_preserved" {
  let content = "AF contig1 3\n" +
    "RD read1 3 chem1 dye1\n" +
    "ACG\n" +
    "q 10 20 30\n" +
    "CT contig1 3 0 3\n" +
    "ACG\n" +
    "q 50 60 70\n" +
    "AF read1 1 3 1 3 1 3 +\n"
  let data = @src.ace_parse(content)
  assert_eq(data.contigs[0].base_qualities.length(), 3)
  assert_eq(data.contigs[0].base_qualities[0], 50)
  assert_eq(data.contigs[0].base_qualities[1], 60)
  assert_eq(data.contigs[0].base_qualities[2], 70)
}

test "ace_read_quality_preserved" {
  let content = "AF contig1 3\n" +
    "RD read1 3 chem1 dye1\n" +
    "ACG\n" +
    "q 10 20 30\n" +
    "CT contig1 3 0 3\n" +
    "ACG\n" +
    "q 50 60 70\n" +
    "AF read1 1 3 1 3 1 3 +\n"
  let data = @src.ace_parse(content)
  assert_eq(data.contigs[0].reads[0].quality.length(), 3)
  assert_eq(data.contigs[0].reads[0].quality[0], 10)
  assert_eq(data.contigs[0].reads[0].quality[1], 20)
  assert_eq(data.contigs[0].reads[0].quality[2], 30)
}

test "ace_coverage_empty" {
  let contig = @src.AceContig::new(
    "contig1",
    "ACGT",
    [],
    [40, 40, 40, 40],
    "",
    false,
  )
  let coverage = @src.ace_read_coverage(contig)
  assert_eq(coverage.get(1), Some(0))
  assert_eq(coverage.get(2), Some(0))
  assert_eq(coverage.get(3), Some(0))
  assert_eq(coverage.get(4), Some(0))
}

test "ace_contig_reads_strand_info" {
  let read1 = @src.AceRead::new(
    "read1",
    "AAAAAA",
    [30, 30, 30, 30, 30, 30],
    1,
    6,
    "+",
    "",
    "chem1",
    "dye1",
  )
  let read2 = @src.AceRead::new(
    "read2",
    "CCCCCC",
    [25, 25, 25, 25, 25, 25],
    1,
    6,
    "-",
    "",
    "chem2",
    "dye2",
  )
  let contig = @src.AceContig::new(
    "contig1",
    "ACGTACGT",
    [read1, read2],
    [40, 40, 40, 40, 40, 40, 40, 40],
    "",
    false,
  )
  assert_eq(contig.reads[0].chemistry, "chem1")
  assert_eq(contig.reads[0].dye, "dye1")
  assert_eq(contig.reads[1].chemistry, "chem2")
  assert_eq(contig.reads[1].dye, "dye2")
}

test "ace_roundtrip_coverage" {
  let read1 = @src.AceRead::new(
    "read1",
    "AAAAAAAAAA",
    [30, 30, 30, 30, 30, 30, 30, 30, 30, 30],
    1,
    10,
    "+",
    "",
    "",
    "",
  )
  let contig = @src.AceContig::new(
    "contig1",
    "AAAAAAAAAA",
    [read1],
    [40, 40, 40, 40, 40, 40, 40, 40, 40, 40],
    "",
    false,
  )
  let data = @src.AceData::new([contig], 0)
  let output = @src.ace_to_string(data)
  let parsed = @src.ace_parse(output)
  assert_eq(parsed.contigs.length(), 1)
  assert_eq(parsed.contigs[0].reads.length(), 1)
  assert_eq(parsed.contigs[0].reads[0].strand, "+")
}

test "ace_quality_to_string_format" {
  let qual : Array[Int] = [10, 20, 30]
  let read = @src.AceRead::new(
    "read1",
    "ACG",
    qual,
    1,
    3,
    "+",
    "",
    "",
    "",
  )
  let contig = @src.AceContig::new(
    "contig1",
    "ACG",
    [read],
    [40, 40, 40],
    "",
    false,
  )
  let data = @src.AceData::new([contig], 0)
  let output = @src.ace_to_string(data)
  assert_true(output.contains("10 20 30"))
}

test "ace_consensus_from_contig" {
  let read = @src.AceRead::new(
    "read1",
    "AAAA",
    [30, 30, 30, 30],
    1,
    4,
    "+",
    "",
    "",
    "",
  )
  let contig = @src.AceContig::new(
    "contig1",
    "ACGT",
    [read],
    [40, 40, 40, 40],
    "",
    false,
  )
  assert_eq(@src.ace_consensus_sequence(contig), "ACGT")
}

test "ace_coverage_values" {
  let read1 = @src.AceRead::new(
    "read1",
    "AAAA",
    [30, 30, 30, 30],
    2,
    3,
    "+",
    "",
    "",
    "",
  )
  let contig = @src.AceContig::new(
    "contig1",
    "ACGTACGT",
    [read1],
    [40, 40, 40, 40, 40, 40, 40, 40],
    "",
    false,
  )
  let cov = @src.ace_read_coverage(contig)
  assert_eq(cov.get(1), Some(0))
  assert_eq(cov.get(2), Some(1))
  assert_eq(cov.get(3), Some(1))
  assert_eq(cov.get(4), Some(0))
}

test "ace_reverse_complemented_default" {
  let contig = @src.AceContig::new(
    "contig1",
    "ACGT",
    [],
    [],
    "",
    false,
  )
  assert_false(contig.reverse_complemented)
}

test "ace_data_threshold" {
  let contig = @src.AceContig::new(
    "contig1",
    "ACGT",
    [],
    [],
    "",
    false,
  )
  let data = @src.AceData::new([contig], 30)
  assert_eq(data.base_qual_threshold, 30)
}