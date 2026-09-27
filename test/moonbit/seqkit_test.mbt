///|
/// Tests for the SeqKit `seq`/`fq2fa`/`fx2tab`/`tab2fx` ports.

///|
test "seqkit_seq_transforms" {
  assert_eq(@src.seqkit_reverse("ACGTAA"), "AATGCA")
  assert_eq(@src.seqkit_complement("ACGTAA"), "TGCATT")
  assert_eq(@src.seqkit_reverse_complement("ACGTAA"), "TTACGT")
  assert_eq(@src.seqkit_remove_gaps("AC-GT.A"), "ACGTA")
  // case-preserving + pass-through of gaps/N
  assert_eq(@src.seqkit_complement("acgt-n"), "tgca-n")
  assert_eq(@src.seqkit_complement("ACRY"), "TGYR")
}

///|
test "seqkit_fx2tab_tab2fx_fasta" {
  let records = @src.parse_fasta(">s1\nACGT\n>s2\nTT\n")
  let tab = @src.seqkit_fx2tab(records, false)
  assert_eq(tab, "s1\tACGT\ns2\tTT")
  let fx = @src.seqkit_tab2fx(tab)
  assert_eq(fx, ">s1\nACGT\n>s2\nTT\n")
  // round-trip: parse back gives same ids + seqs
  let back = @src.parse_fasta(fx)
  assert_eq(back.length(), 2)
  assert_eq(back[0].id, "s1")
  assert_eq(back[0].seq.to_string(), "ACGT")
}

///|
test "seqkit_fx2tab_tab2fx_fastq" {
  let records = @src.parse_fastq("@s1\nACGT\n+\n????\n")
  let tab = @src.seqkit_fx2tab(records, true)
  assert_eq(tab, "s1\tACGT\t????")
  let fx = @src.seqkit_tab2fx(tab)
  assert_eq(fx, "@s1\nACGT\n+\n????\n")
  // fq2fa drops quality
  let fa = @src.seqkit_fq2fa(records)
  assert_eq(fa, ">s1\nACGT\n")
}

// ============================================================
// seqkit amplicon
// ============================================================

test "seqkit_amplicon_basic" {
  // seq = AAACCCGGGTTT；fwd=AAA 在 0-2，rev=AAAC 的 rc=GTTT 在 8-11
  // amplicon = seq[3:8] = CCCGG
  match @src.seqkit_amplicon("AAACCCGGGTTT", "AAA", "AAAC") {
    Some(a) => assert_eq(a, "CCCGG")
    None => fail("should find amplicon")
  }
}

test "seqkit_amplicon_missing_primer" {
  assert_true(@src.seqkit_amplicon("ACGT", "TTT", "AAA") is None)
  assert_true(@src.seqkit_amplicon("ACGT", "ACG", "CCC") is None)
}

test "seqkit_amplicon_inverted_empty" {
  // rev primer 的反向互补出现在 fwd 之前 → 区域为空/inverted → None
  match @src.seqkit_amplicon("AAAACGTTT", "ACG", "TTT") {
    None => ()
    Some(_) => fail("should be None for inverted region")
  }
}
