///|
/// Tests for scikit-bio alignment module.

///|
test "dna_creation" {
  let dna = @src.DNA::new("ACGT")
  assert_eq(dna.to_string(), "ACGT")
  assert_eq(dna.length(), 4)
}

///|
test "rna_creation" {
  let rna = @src.RNA::new("ACGU")
  assert_eq(rna.to_string(), "ACGU")
  assert_eq(rna.length(), 4)
}

///|
test "protein_creation" {
  let prot = @src.Protein::new("ACDEFGHIKLMNPQRSTVWY")
  assert_eq(prot.to_string(), "ACDEFGHIKLMNPQRSTVWY")
  assert_eq(prot.length(), 20)
}

///|
test "alignment_tabular_msa_creation" {
  let dna1 = @src.DNA::new("ACGT")
  let dna2 = @src.DNA::new("ACCT")
  let seqs : Array[@src.DNA] = Array::new()
  seqs.push(dna1)
  seqs.push(dna2)
  let msa = @src.TabularMSA::from_dna(seqs)

  assert_eq(msa.num_records(), 2)
  assert_eq(msa.get_alignment_length(), 4)
}

///|
test "alignment_tabular_msa_column" {
  let dna1 = @src.DNA::new("ACGT")
  let dna2 = @src.DNA::new("ACCT")
  let seqs : Array[@src.DNA] = Array::new()
  seqs.push(dna1)
  seqs.push(dna2)
  let msa = @src.TabularMSA::from_dna(seqs)

  assert_eq(msa.column(0).unwrap(), "AA")
  assert_eq(msa.column(1).unwrap(), "CC")
  assert_eq(msa.column(2).unwrap(), "GC")
  assert_eq(msa.column(3).unwrap(), "TT")
}

///|
test "alignment_tabular_msa_consensus" {
  let dna1 = @src.DNA::new("ACGT")
  let dna2 = @src.DNA::new("ACCT")
  let dna3 = @src.DNA::new("AGGT")
  let seqs : Array[@src.DNA] = Array::new()
  seqs.push(dna1)
  seqs.push(dna2)
  seqs.push(dna3)
  let msa = @src.TabularMSA::from_dna(seqs)

  let consensus = msa.consensus()
  assert_eq(consensus.to_string(), "ACGT")
}

///|
test "alignment_tabular_msa_conservation" {
  let dna1 = @src.DNA::new("ACGT")
  let dna2 = @src.DNA::new("ACCT")
  let dna3 = @src.DNA::new("AGGT")
  let seqs : Array[@src.DNA] = Array::new()
  seqs.push(dna1)
  seqs.push(dna2)
  seqs.push(dna3)
  let msa = @src.TabularMSA::from_dna(seqs)

  let col2 = msa.column(2).unwrap()
  assert_eq(col2, "GCG")

  let conservation = msa.conservation()
  assert_eq(conservation.length(), 4)
  assert_eq(conservation[0], 1.0)
  assert_eq(conservation[1], 0.6666666666666666)
  assert_eq(conservation[2], 0.6666666666666666)
  assert_eq(conservation[3], 1.0)
}

///|
test "alignment_global_pairwise_nucleotide" {
  let seq1 = @src.DNA::new("ACGT")
  let seq2 = @src.DNA::new("AGT")

  let (msa, score, pos1, pos2) = @src.global_pairwise_align_nucleotide(
    seq1, seq2,
  )

  assert_eq(msa.num_records(), 2)
  assert_eq(pos1.0, 0)
  assert_eq(pos1.1, 4)
  assert_eq(pos2.0, 0)
  assert_eq(pos2.1, 3)
}

///|
test "alignment_local_pairwise_nucleotide" {
  let seq1 = @src.DNA::new("ACGTACGT")
  let seq2 = @src.DNA::new("CGT")

  let (msa, score, pos1, pos2) = @src.local_pairwise_align_nucleotide(
    seq1, seq2,
  )

  assert_eq(msa.num_records(), 2)

  let col0 = msa.column(0).unwrap()
  assert_eq(col0, "CC")
  let col1 = msa.column(1).unwrap()
  assert_eq(col1, "GG")
  let col2 = msa.column(2).unwrap()
  assert_eq(col2, "TT")
}

///|
test "alignment_global_pairwise_protein" {
  let seq1 = @src.Protein::new("ACDE")
  let seq2 = @src.Protein::new("ACE")

  let (msa, score, pos1, pos2) = @src.global_pairwise_align_protein(seq1, seq2)

  assert_eq(msa.num_records(), 2)
}

///|
test "alignment_local_pairwise_protein" {
  let seq1 = @src.Protein::new("ACDEFGHI")
  let seq2 = @src.Protein::new("EFG")

  let (msa, score, pos1, pos2) = @src.local_pairwise_align_protein(seq1, seq2)

  assert_eq(msa.num_records(), 2)
}

///|
test "alignment_score" {
  let dna1 = @src.DNA::new("ACGT")
  let dna2 = @src.DNA::new("ACGT")
  let seqs : Array[@src.DNA] = Array::new()
  seqs.push(dna1)
  seqs.push(dna2)
  let msa = @src.TabularMSA::from_dna(seqs)

  let score = @src.align_score(msa)
  assert_eq(score, 6.0)
}

///|
test "alignment_global_score" {
  let seq1 = @src.DNA::new("ATCG")
  let seq2 = @src.DNA::new("ATAG")

  let (msa, _, _, _) = @src.global_pairwise_align_nucleotide(seq1, seq2)

  assert_eq(msa.get_alignment_length(), 4)
}
