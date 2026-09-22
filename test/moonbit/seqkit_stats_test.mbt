///|
/// Tests for the SeqKit `stats` port. Values are hand-computed per seqkit's
/// `stat.go` semantics (GC% = G/C count ÷ total length, quartiles truncated,
/// AvgQual = -10·log10(Σ10^(-q/10)/len)).

///|
test "seqkit_stats_fasta" {
  let records = @src.parse_fasta(">s1\nACGTACGTAC\n>s2\nACGT\n>s3\nACGTACGT\n")
  // lengths [10,4,8]; GC counts [5,2,4]; no gaps/N
  let st = @src.seqkit_stats(records, true, false)
  assert_eq(st.num_seqs, 3)
  assert_eq(st.sum_len, 22)
  assert_eq(st.min_len, 4)
  assert_eq(st.max_len, 10)
  assert_true(assert_close(st.avg_len, 22.0 / 3.0, 0.05))
  assert_eq(st.q1, 4)
  assert_eq(st.q2, 8)
  assert_eq(st.q3, 10)
  assert_eq(st.sum_gap, 0)
  assert_eq(st.n50, 8)
  assert_eq(st.l50, 2)
  assert_true(assert_close(st.gc_pct, 50.0, 0.01))
  assert_eq(st.sum_n, 0)
  assert_true(assert_close(st.avg_qual, 0.0, 1.0e-12))
}

///|
test "seqkit_stats_fastq_qual" {
  let records = @src.parse_fastq("@s1\nACGT\n+\n????\n")
  // all qualities phred 30 ('?' = 33+30)
  let st = @src.seqkit_stats(records, true, true)
  assert_eq(st.num_seqs, 1)
  assert_eq(st.sum_len, 4)
  assert_true(assert_close(st.gc_pct, 50.0, 0.01))
  assert_true(assert_close(st.q20_pct, 100.0, 0.01))
  assert_true(assert_close(st.q30_pct, 100.0, 0.01))
  assert_true(assert_close(st.avg_qual, 30.0, 0.05))
}

///|
test "seqkit_stats_gap_n" {
  let records = @src.parse_fasta(">s1\nAC-GT.NACNNT\n")
  // "AC-GT.NACNNT": len 12, gaps '-','.' =2, N =3, GC: A0 C2 G1 T1 = C+G=3
  let st = @src.seqkit_stats(records, true, false)
  assert_eq(st.sum_len, 12)
  assert_eq(st.sum_gap, 2)
  assert_eq(st.sum_n, 3)
}