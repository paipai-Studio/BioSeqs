///|
/// Tests for biotite_application_mafft.mbt.
test "mafft command default" {
  let cmd = @src.biotite_mafft_command("input.fasta")
  assert_eq(cmd, "mafft --quiet --anysymbol input.fasta")
}

///|
test "mafft command global pair" {
  let cmd = @src.biotite_mafft_command(
    "input.fasta",
    bin_path="/usr/bin/mafft",
    global_pair=true,
  )
  assert_eq(
    cmd, "/usr/bin/mafft --quiet --anysymbol --globalpair --maxiterate 1000 input.fasta",
  )
}

///|
test "mafft command explicit max iterate" {
  let cmd = @src.biotite_mafft_command(
    "input.fasta",
    quiet=false,
    any_symbol=false,
    max_iterate=100,
  )
  assert_eq(cmd, "mafft --maxiterate 100 input.fasta")
}

///|
test "mafft parse fasta" {
  let (names, seqs) = @src.biotite_mafft_parse_fasta(
    ">alpha desc\nAC\nGU\n>beta\nAG\nGU\n",
  )
  assert_eq(names.length(), 2)
  assert_eq(seqs.length(), 2)
  assert_eq(names[0], "alpha")
  assert_eq(names[1], "beta")
  assert_eq(seqs[0], "ACGU")
  assert_eq(seqs[1], "AGGU")
}

///|
test "mafft local alignment" {
  let result = @src.biotite_mafft_local(["ACGU", "ACGC"], names=["s1", "s2"])
  assert_eq(result.status, "computed-local")
  assert_eq(result.aligned_sequences.length(), 2)
  assert_eq(result.names.length(), 2)
  assert_eq(
    result.aligned_sequences[0].length(),
    result.aligned_sequences[1].length(),
  )
  assert_true(result.command.contains("mafft"))
}

///|
test "mafft to fasta roundtrip" {
  let result = @src.biotite_mafft_local(["ACGU", "ACGC"], names=["s1", "s2"])
  let fasta = @src.biotite_mafft_to_fasta(result)
  let (names, seqs) = @src.biotite_mafft_parse_fasta(fasta)
  assert_eq(names[0], "s1")
  assert_eq(names[1], "s2")
  assert_eq(seqs[0], result.aligned_sequences[0])
  assert_eq(seqs[1], result.aligned_sequences[1])
}

///|
test "mafft parse empty" {
  let (names, seqs) = @src.biotite_mafft_parse_fasta("")
  assert_eq(names.length(), 0)
  assert_eq(seqs.length(), 0)
}
