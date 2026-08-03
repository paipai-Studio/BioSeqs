///|
/// Tests for Bio.Wise module - Wise2 (GeneWise/GenomeWise) output parser.

// ---------------------------------------------------------------------------
// WiseBlockType enum
// ---------------------------------------------------------------------------

test "wise_block_type_to_string_exon" {
  let bt = @src.wise_make_block_type("exon")
  assert_eq(bt.to_string(), "exon")
}

test "wise_block_type_to_string_intron" {
  let bt = @src.wise_make_block_type("intron")
  assert_eq(bt.to_string(), "intron")
}

test "wise_block_type_to_string_match" {
  let bt = @src.wise_make_block_type("match")
  assert_eq(bt.to_string(), "match")
}

test "wise_block_type_to_string_mismatch" {
  let bt = @src.wise_make_block_type("mismatch")
  assert_eq(bt.to_string(), "mismatch")
}

test "wise_block_type_to_string_insertion" {
  let bt = @src.wise_make_block_type("insertion")
  assert_eq(bt.to_string(), "insertion")
}

test "wise_block_type_to_string_deletion" {
  let bt = @src.wise_make_block_type("deletion")
  assert_eq(bt.to_string(), "deletion")
}

test "wise_block_type_eq_via_is_pattern" {
  // Enum variants can be compared using `is` pattern matching.
  let bt = @src.wise_make_block_type("exon")
  assert_true(bt is @src.WiseBlockType::WiseExon)
  assert_false(bt is @src.WiseBlockType::WiseIntron)
}

// ---------------------------------------------------------------------------
// WiseExon
// ---------------------------------------------------------------------------

test "wise_exon_new_defaults" {
  let e = @src.WiseExon::new()
  assert_eq(e.start, 0)
  assert_eq(e.end, 0)
  assert_eq(e.phase, 0)
  assert_eq(e.end_phase, 0)
  assert_eq(e.score, 0.0)
  assert_eq(e.protein_start, 0)
  assert_eq(e.protein_end, 0)
}

test "wise_exon_new_with_params" {
  let e = @src.WiseExon::new(
    start=100,
    end=199,
    phase=0,
    end_phase=2,
    score=45.5,
    protein_start=1,
    protein_end=33,
  )
  assert_eq(e.start, 100)
  assert_eq(e.end, 199)
  assert_eq(e.phase, 0)
  assert_eq(e.end_phase, 2)
  assert_eq(e.score, 45.5)
  assert_eq(e.protein_start, 1)
  assert_eq(e.protein_end, 33)
}

test "wise_exon_length" {
  let e = @src.WiseExon::new(start=0, end=99)
  assert_eq(e.length(), 100)
}

test "wise_exon_length_single_base" {
  let e = @src.WiseExon::new(start=50, end=50)
  assert_eq(e.length(), 1)
}

// ---------------------------------------------------------------------------
// WiseIntron
// ---------------------------------------------------------------------------

test "wise_intron_new_defaults" {
  let i = @src.WiseIntron::new()
  assert_eq(i.start, 0)
  assert_eq(i.end, 0)
  assert_eq(i.donor_score, 0.0)
  assert_eq(i.acceptor_score, 0.0)
  // length auto-calculated: end - start + 1 = 1
  assert_eq(i.length, 1)
}

test "wise_intron_new_with_params" {
  let i = @src.WiseIntron::new(
    start=100,
    end=199,
    donor_score=12.5,
    acceptor_score=8.3,
  )
  assert_eq(i.start, 100)
  assert_eq(i.end, 199)
  assert_eq(i.donor_score, 12.5)
  assert_eq(i.acceptor_score, 8.3)
  // length auto-calculated: 199 - 100 + 1 = 100
  assert_eq(i.length, 100)
}

test "wise_intron_length_auto_calc_normal" {
  let i = @src.WiseIntron::new(start=200, end=299)
  assert_eq(i.length, 100)
}

test "wise_intron_length_auto_calc_end_before_start" {
  // When end < start, length should be 0
  let i = @src.WiseIntron::new(start=300, end=100)
  assert_eq(i.length, 0)
}

// ---------------------------------------------------------------------------
// WiseAlignmentColumn
// ---------------------------------------------------------------------------

test "wise_alignment_column_new" {
  let c = @src.WiseAlignmentColumn::new(
    protein_char='M',
    gene_codon="ATG",
    match_type=" ",
    gene_position=0,
    protein_position=0,
  )
  assert_eq(c.protein_char, 'M')
  assert_eq(c.gene_codon, "ATG")
  assert_eq(c.match_type, " ")
  assert_eq(c.gene_position, 0)
  assert_eq(c.protein_position, 0)
}

test "wise_alignment_column_new_default_match_type" {
  let c = @src.WiseAlignmentColumn::new(
    protein_char='X',
    gene_codon="TAG",
  )
  assert_eq(c.protein_char, 'X')
  assert_eq(c.gene_codon, "TAG")
  // Default match_type is " "
  assert_eq(c.match_type, " ")
  assert_eq(c.gene_position, 0)
  assert_eq(c.protein_position, 0)
}

// ---------------------------------------------------------------------------
// WiseResult construction and basic accessors
// ---------------------------------------------------------------------------

test "wise_result_new_empty" {
  let r = @src.WiseResult::new()
  assert_eq(r.protein_id, "")
  assert_eq(r.gene_id, "")
  assert_eq(r.score, 0.0)
  assert_eq(r.bits_score, 0.0)
  assert_eq(r.get_num_exons(), 0)
  assert_eq(r.get_num_introns(), 0)
  assert_eq(r.get_total_exon_length(), 0)
}

test "wise_result_add_and_get_exons" {
  let r = @src.WiseResult::new()
  let e1 = @src.WiseExon::new(start=0, end=99)
  let e2 = @src.WiseExon::new(start=200, end=299)
  r.add_exon(e1)
  r.add_exon(e2)
  let exons = r.get_exons()
  assert_eq(exons.length(), 2)
  assert_eq(exons[0].start, 0)
  assert_eq(exons[1].start, 200)
}

test "wise_result_add_and_get_introns" {
  let r = @src.WiseResult::new()
  let i1 = @src.WiseIntron::new(start=100, end=199)
  r.add_intron(i1)
  let introns = r.get_introns()
  assert_eq(introns.length(), 1)
  assert_eq(introns[0].start, 100)
  assert_eq(introns[0].end, 199)
}

test "wise_result_get_num_exons" {
  let r = @src.WiseResult::new()
  assert_eq(r.get_num_exons(), 0)
  r.add_exon(@src.WiseExon::new(start=0, end=50))
  assert_eq(r.get_num_exons(), 1)
  r.add_exon(@src.WiseExon::new(start=60, end=110))
  r.add_exon(@src.WiseExon::new(start=120, end=170))
  assert_eq(r.get_num_exons(), 3)
}

test "wise_result_get_num_introns" {
  let r = @src.WiseResult::new()
  assert_eq(r.get_num_introns(), 0)
  r.add_intron(@src.WiseIntron::new(start=51, end=59))
  assert_eq(r.get_num_introns(), 1)
  r.add_intron(@src.WiseIntron::new(start=111, end=119))
  assert_eq(r.get_num_introns(), 2)
}

test "wise_result_get_total_exon_length" {
  let r = @src.WiseResult::new()
  // Exon 1: 0-99 -> length 100
  r.add_exon(@src.WiseExon::new(start=0, end=99))
  // Exon 2: 200-299 -> length 100
  r.add_exon(@src.WiseExon::new(start=200, end=299))
  // Exon 3: 400-479 -> length 80
  r.add_exon(@src.WiseExon::new(start=400, end=479))
  assert_eq(r.get_total_exon_length(), 280)
}

test "wise_result_get_total_exon_length_empty" {
  let r = @src.WiseResult::new()
  assert_eq(r.get_total_exon_length(), 0)
}

test "wise_result_set_protein_id" {
  let r = @src.WiseResult::new()
  r.set_protein_id("P12345")
  assert_eq(r.protein_id, "P12345")
}

test "wise_result_set_gene_id" {
  let r = @src.WiseResult::new()
  r.set_gene_id("G67890")
  assert_eq(r.gene_id, "G67890")
}

test "wise_result_set_score" {
  let r = @src.WiseResult::new()
  r.set_score(145.32)
  assert_eq(r.score, 145.32)
}

test "wise_result_set_bits_score" {
  let r = @src.WiseResult::new()
  r.set_bits_score(52.18)
  assert_eq(r.bits_score, 52.18)
}

test "wise_result_add_alignment_column" {
  let r = @src.WiseResult::new()
  let c = @src.WiseAlignmentColumn::new(
    protein_char='M',
    gene_codon="ATG",
  )
  r.add_alignment_column(c)
  let aln = r.get_alignment()
  assert_eq(aln.length(), 1)
  assert_eq(aln[0].protein_char, 'M')
}

// ---------------------------------------------------------------------------
// wise_parse
// ---------------------------------------------------------------------------

test "wise_parse_protein_id" {
  let r = @src.wise_parse("Protein: P12345\n")
  assert_eq(r.protein_id, "P12345")
}

test "wise_parse_gene_id" {
  let r = @src.wise_parse("Gene: G67890\n")
  assert_eq(r.gene_id, "G67890")
}

test "wise_parse_score" {
  let r = @src.wise_parse("Score: 145.32\n")
  // Use tolerance for floating-point comparison (parse_double precision)
  assert_true(r.score > 145.31 && r.score < 145.33)
}

test "wise_parse_bits_score" {
  let r = @src.wise_parse("Bits: 52.18\n")
  assert_eq(r.bits_score, 52.18)
}

test "wise_parse_exons" {
  let r = @src.wise_parse(
    "Exon 1: 0-99 (phase 0) score 45.5 protein 1-33\nExon 2: 200-299 (phase 0) score 52.3 protein 34-66\n",
  )
  assert_eq(r.get_num_exons(), 2)
  let exons = r.get_exons()
  assert_eq(exons[0].start, 0)
  assert_eq(exons[0].end, 99)
  assert_eq(exons[0].phase, 0)
  assert_eq(exons[0].score, 45.5)
  assert_eq(exons[0].protein_start, 1)
  assert_eq(exons[0].protein_end, 33)
  assert_eq(exons[1].start, 200)
  assert_eq(exons[1].end, 299)
  assert_eq(exons[1].score, 52.3)
}

test "wise_parse_introns" {
  let r = @src.wise_parse(
    "Intron 1: 100-199 donor 12.5 acceptor 8.3\nIntron 2: 300-399 donor 10.2 acceptor 9.1\n",
  )
  assert_eq(r.get_num_introns(), 2)
  let introns = r.get_introns()
  assert_eq(introns[0].start, 100)
  assert_eq(introns[0].end, 199)
  assert_eq(introns[0].donor_score, 12.5)
  assert_eq(introns[0].acceptor_score, 8.3)
  assert_eq(introns[0].length, 100)
  assert_eq(introns[1].start, 300)
  assert_eq(introns[1].end, 399)
  assert_eq(introns[1].donor_score, 10.2)
  assert_eq(introns[1].acceptor_score, 9.1)
}

test "wise_parse_empty_input" {
  let r = @src.wise_parse("")
  assert_eq(r.protein_id, "")
  assert_eq(r.gene_id, "")
  assert_eq(r.score, 0.0)
  assert_eq(r.bits_score, 0.0)
  assert_eq(r.get_num_exons(), 0)
  assert_eq(r.get_num_introns(), 0)
}

test "wise_parse_full_sample" {
  let r = @src.wise_parse(@src.wise_sample_output())
  assert_eq(r.protein_id, "P12345")
  assert_eq(r.gene_id, "G67890")
  // Use tolerance for floating-point comparison (parse_double precision)
  assert_true(r.score > 145.31 && r.score < 145.33)
  assert_eq(r.bits_score, 52.18)
  assert_eq(r.get_num_exons(), 3)
  assert_eq(r.get_num_introns(), 2)
  assert_eq(r.get_total_exon_length(), 280)
}

// ---------------------------------------------------------------------------
// Sample data
// ---------------------------------------------------------------------------

test "wise_sample_output_content" {
  let s = @src.wise_sample_output()
  assert_true(s.contains("Protein: P12345"))
  assert_true(s.contains("Gene: G67890"))
  assert_true(s.contains("Score: 145.32"))
  assert_true(s.contains("Bits: 52.18"))
  assert_true(s.contains("Exon 1: 0-99"))
  assert_true(s.contains("Intron 1: 100-199"))
}

test "wise_sample_result" {
  let r = @src.wise_sample()
  assert_eq(r.protein_id, "P12345")
  assert_eq(r.gene_id, "G67890")
  assert_eq(r.get_num_exons(), 3)
  assert_eq(r.get_num_introns(), 2)
}

// ---------------------------------------------------------------------------
// wise_write
// ---------------------------------------------------------------------------

test "wise_write_basic" {
  let r = @src.WiseResult::new()
  r.set_protein_id("P12345")
  r.set_gene_id("G67890")
  r.set_score(145.32)
  r.set_bits_score(52.18)
  let s = @src.wise_write(r)
  assert_true(s.has_prefix("GeneWise Result"))
  assert_true(s.contains("Protein: P12345"))
  assert_true(s.contains("Gene: G67890"))
  assert_true(s.contains("Score: 145.32"))
  assert_true(s.contains("Bits: 52.18"))
}

test "wise_write_with_exons_and_introns" {
  let r = @src.WiseResult::new()
  r.set_protein_id("P1")
  r.set_gene_id("G1")
  r.add_exon(@src.WiseExon::new(start=0, end=99, phase=0, score=45.5, protein_start=1, protein_end=33))
  r.add_intron(@src.WiseIntron::new(start=100, end=199, donor_score=12.5, acceptor_score=8.3))
  let s = @src.wise_write(r)
  assert_true(s.contains("Exon 1: 0-99"))
  assert_true(s.contains("(phase 0)"))
  assert_true(s.contains("score 45.5"))
  assert_true(s.contains("protein 1-33"))
  assert_true(s.contains("Intron 1: 100-199"))
  assert_true(s.contains("donor 12.5"))
  assert_true(s.contains("acceptor 8.3"))
}

test "wise_write_roundtrip" {
  // Parse sample, write it back, verify key fields are present
  let r = @src.wise_sample()
  let s = @src.wise_write(r)
  assert_true(s.contains("Protein: P12345"))
  assert_true(s.contains("Gene: G67890"))
  assert_true(s.contains("Exon 1: 0-99"))
  assert_true(s.contains("Exon 2: 200-299"))
  assert_true(s.contains("Exon 3: 400-479"))
  assert_true(s.contains("Intron 1: 100-199"))
  assert_true(s.contains("Intron 2: 300-399"))
}

// ---------------------------------------------------------------------------
// wise_gene_structure
// ---------------------------------------------------------------------------

test "wise_gene_structure_basic" {
  let r = @src.WiseResult::new()
  r.add_exon(@src.WiseExon::new(start=0, end=99, phase=0))
  r.add_exon(@src.WiseExon::new(start=200, end=299, phase=0))
  r.add_intron(@src.WiseIntron::new(start=100, end=199))
  let s = @src.wise_gene_structure(r)
  assert_true(s.contains("Exon 1: 0-99"))
  assert_true(s.contains("length 100 nt"))
  assert_true(s.contains("phase 0"))
  assert_true(s.contains("Exon 2: 200-299"))
  assert_true(s.contains("Intron 1: 100-199"))
}

test "wise_gene_structure_empty" {
  let r = @src.WiseResult::new()
  let s = @src.wise_gene_structure(r)
  // Empty result should produce an empty structure string
  assert_eq(s.length(), 0)
}

test "wise_gene_structure_sample" {
  let r = @src.wise_sample()
  let s = @src.wise_gene_structure(r)
  // Three exons with lengths 100, 100, 80
  assert_true(s.contains("Exon 1: 0-99 (length 100 nt, phase 0)"))
  assert_true(s.contains("Exon 2: 200-299 (length 100 nt, phase 0)"))
  assert_true(s.contains("Exon 3: 400-479 (length 80 nt, phase 0)"))
  assert_true(s.contains("Intron 1: 100-199 (length 100 nt)"))
  assert_true(s.contains("Intron 2: 300-399 (length 100 nt)"))
}

// ---------------------------------------------------------------------------
// wise_translate_gene
// ---------------------------------------------------------------------------

test "wise_translate_gene_single_exon" {
  // Exon 0-8 (9 nt) -> ATG GCC GGT -> M A G
  let r = @src.WiseResult::new()
  r.add_exon(@src.WiseExon::new(start=0, end=8))
  let gene = "ATGGCCGGT"
  let protein = @src.wise_translate_gene(r, gene)
  assert_eq(protein, "MAG")
}

test "wise_translate_gene_multiple_codons" {
  // Exon 0-11 (12 nt) -> ATG GCC GGT AAA -> M A G K
  let r = @src.WiseResult::new()
  r.add_exon(@src.WiseExon::new(start=0, end=11))
  let gene = "ATGGCCGGTAAA"
  let protein = @src.wise_translate_gene(r, gene)
  assert_eq(protein, "MAGK")
}

test "wise_translate_gene_stop_codon" {
  // Exon 0-8 -> ATG TAA TAG -> M * *
  let r = @src.WiseResult::new()
  r.add_exon(@src.WiseExon::new(start=0, end=8))
  let gene = "ATGTAATAG"
  let protein = @src.wise_translate_gene(r, gene)
  assert_eq(protein.length(), 3)
  assert_eq(protein.unsafe_get(0), 'M')
  assert_eq(protein.unsafe_get(1), '*')
  assert_eq(protein.unsafe_get(2), '*')
}

// ---------------------------------------------------------------------------
// wise_percent_identity
// ---------------------------------------------------------------------------

test "wise_percent_identity_full" {
  // Translation matches protein exactly -> 100.0%
  let r = @src.WiseResult::new()
  r.add_exon(@src.WiseExon::new(start=0, end=8))
  let gene = "ATGGCCGGT"
  let protein = "MAG"
  let pct = @src.wise_percent_identity(r, protein, gene)
  assert_eq(pct, 100.0)
}

test "wise_percent_identity_none" {
  // Translation does not match protein -> 0.0%
  let r = @src.WiseResult::new()
  r.add_exon(@src.WiseExon::new(start=0, end=8))
  let gene = "ATGGCCGGT" // translates to "MAG"
  let protein = "XXX"
  let pct = @src.wise_percent_identity(r, protein, gene)
  assert_eq(pct, 0.0)
}

test "wise_percent_identity_partial" {
  // Translation: MAG, protein: MAP -> 2/3 match -> ~66.67%
  let r = @src.WiseResult::new()
  r.add_exon(@src.WiseExon::new(start=0, end=8))
  let gene = "ATGGCCGGT" // translates to "MAG"
  let protein = "MAP"
  let pct = @src.wise_percent_identity(r, protein, gene)
  assert_true(pct > 66.0 && pct < 67.0)
}

test "wise_percent_identity_empty_protein" {
  let r = @src.WiseResult::new()
  r.add_exon(@src.WiseExon::new(start=0, end=8))
  let gene = "ATGGCCGGT"
  let pct = @src.wise_percent_identity(r, "", gene)
  assert_eq(pct, 0.0)
}

// ---------------------------------------------------------------------------
// wise_summary
// ---------------------------------------------------------------------------

test "wise_summary_basic" {
  let r = @src.wise_sample()
  let s = @src.wise_summary(r)
  assert_true(s.contains("GeneWise Result Summary:"))
  assert_true(s.contains("Protein ID: P12345"))
  assert_true(s.contains("Gene ID: G67890"))
  assert_true(s.contains("Score: 145.32"))
  assert_true(s.contains("Bits: 52.18"))
  assert_true(s.contains("Exons: 3"))
  assert_true(s.contains("Introns: 2"))
  assert_true(s.contains("Total exon length: 280 nt"))
}

test "wise_summary_empty" {
  let r = @src.WiseResult::new()
  let s = @src.wise_summary(r)
  assert_true(s.contains("GeneWise Result Summary:"))
  assert_true(s.contains("Exons: 0"))
  assert_true(s.contains("Introns: 0"))
  assert_true(s.contains("Total exon length: 0 nt"))
}

// ---------------------------------------------------------------------------
// GenomeWiseSegment
// ---------------------------------------------------------------------------

test "genome_wise_segment_new_defaults" {
  let s = @src.GenomeWiseSegment::new()
  assert_eq(s.segment_id, "")
  assert_eq(s.start, 0)
  assert_eq(s.end, 0)
  assert_eq(s.score, 0.0)
  // wise_result should be an empty WiseResult
  assert_eq(s.wise_result.get_num_exons(), 0)
}

test "genome_wise_segment_new_with_params" {
  let s = @src.GenomeWiseSegment::new(
    segment_id="seg1",
    start=1000,
    end=2000,
    score=85.5,
  )
  assert_eq(s.segment_id, "seg1")
  assert_eq(s.start, 1000)
  assert_eq(s.end, 2000)
  assert_eq(s.score, 85.5)
}

// ---------------------------------------------------------------------------
// GenomeWiseResult
// ---------------------------------------------------------------------------

test "genome_wise_result_new_empty" {
  let r = @src.GenomeWiseResult::new()
  assert_eq(r.gene_id, "")
  assert_eq(r.total_score, 0.0)
  assert_eq(r.get_num_segments(), 0)
}

test "genome_wise_result_add_and_get_segments" {
  let r = @src.GenomeWiseResult::new()
  let s1 = @src.GenomeWiseSegment::new(segment_id="seg1", start=100, end=200)
  let s2 = @src.GenomeWiseSegment::new(segment_id="seg2", start=300, end=400)
  r.add_segment(s1)
  r.add_segment(s2)
  let segs = r.get_segments()
  assert_eq(segs.length(), 2)
  assert_eq(segs[0].segment_id, "seg1")
  assert_eq(segs[1].segment_id, "seg2")
}

test "genome_wise_result_get_num_segments" {
  let r = @src.GenomeWiseResult::new()
  assert_eq(r.get_num_segments(), 0)
  r.add_segment(@src.GenomeWiseSegment::new(segment_id="s1"))
  assert_eq(r.get_num_segments(), 1)
  r.add_segment(@src.GenomeWiseSegment::new(segment_id="s2"))
  r.add_segment(@src.GenomeWiseSegment::new(segment_id="s3"))
  assert_eq(r.get_num_segments(), 3)
}
