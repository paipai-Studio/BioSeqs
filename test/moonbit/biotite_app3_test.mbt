///|
/// Tests for biotite_app3.mbt — request construction and response parsing.
test "biotite3_alphafold_request_and_parse" {
  assert_true(
    @src.biotite3_alphafold_request_url().has_prefix("https://api.esmatlas.com"),
  )
  assert_eq(@src.biotite3_alphafold_request_body("ACDE"), "ACDE")
  let pdb = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C\nEND\n"
  let result = @src.biotite3_alphafold_parse_pdb(pdb)
  assert_eq(result.status, "parsed")
  assert_eq(result.confidence_scores.length(), 1)
  assert_eq(result.sequence, "A")
  assert_true(@moonbit.assert_num_close(result.mean_confidence, 90.0, 1.0e-9))
  let normalized = @src.biotite3_alphafold_result_to_pdb(result)
  assert_true(normalized.has_prefix("ATOM"))
}

///|
test "biotite3_alphafold_bfactor_clamps_plddt" {
  let scores = @src.biotite3_alphafold_confidence_to_bfactor([120.0, -5.0, 70.0])
  assert_eq(scores.length(), 3)
  assert_true(@moonbit.assert_num_close(scores[0], 100.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(scores[1], 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(scores[2], 70.0, 1.0e-9))
}

///|
test "biotite3_sra_url_and_tsv_parse" {
  let url = @src.biotite3_sra_search_url("SRR000001")
  assert_true(url.has_prefix("https://www.ebi.ac.uk/ena/portal/api/filereport"))
  assert_true(url.contains("SRR000001"))
  let fastq = @src.biotite3_sra_fastq_url("SRR000001", 10)
  assert_true(fastq.contains("limit=10"))
  let tsv = "run_accession\tstudy_accession\tsample_accession\texperiment_accession\tlibrary_strategy\tlibrary_source\tlibrary_selection\tread_count\tread_length\tscientific_name\n" +
    "SRR1\tSRP1\tSRS1\tSRX1\tWGS\tGENOMIC\tRANDOM\t100\t150\tE. coli\n"
  let records = @src.biotite3_sra_parse_tsv(tsv)
  assert_eq(records.length(), 1)
  assert_eq(records[0].run_accession, "SRR1")
  assert_eq(records[0].n_reads, 100)
  assert_eq(records[0].read_length, 150)
  assert_eq(records[0].organism, "E. coli")
  let xml = @src.biotite3_sra_metadata_xml(records[0])
  assert_true(xml.contains("SRR1"))
  assert_true(xml.contains("read_count=\"100\""))
}

///|
test "biotite3_clustal_parse_and_local_alignment" {
  let aln = "CLUSTAL O(1.2.4) multiple sequence alignment\n\nseq1  ACDE\nseq2  AC-E\n"
  let parsed = @src.biotite3_clustal_omega_parse_aln(aln)
  assert_eq(parsed.status, "parsed")
  assert_eq(parsed.names.length(), 2)
  assert_eq(parsed.aligned_sequences[0], "ACDE")
  assert_eq(parsed.aligned_sequences[1], "AC-E")
  assert_true(parsed.alignment_score > 0.0)
  let local_result = @src.biotite3_clustal_omega_local(["ACGT", "ACGT"], [
    "a", "b",
  ])
  assert_eq(local_result.status, "computed-local")
  assert_eq(local_result.aligned_sequences.length(), 2)
  assert_true(@src.biotite3_clustal_omega_to_fasta(local_result).contains(">a"))
}
