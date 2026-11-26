///|
/// Tests for biotite_application_muscle.mbt.
test "muscle5 command" {
  let cmd = @src.biotite_muscle5_command("in.fasta", "out.fasta")
  assert_eq(cmd, "muscle -align in.fasta -output out.fasta")
}

///|
test "muscle5 super5 command" {
  let cmd = @src.biotite_muscle5_command(
    "in.fasta",
    "out.fasta",
    bin_path="/usr/bin/muscle",
    super5=true,
  )
  assert_eq(cmd, "/usr/bin/muscle -super5 in.fasta -output out.fasta")
}

///|
test "muscle3 command" {
  let cmd = @src.biotite_muscle3_command("in.fasta", "out.fasta")
  assert_eq(cmd, "muscle -in in.fasta -out out.fasta")
}

///|
test "muscle parse fasta" {
  let (names, seqs) = @src.biotite_muscle_parse_fasta(
    ">seq_a\nMKTA\n>seq_b\nMKAA\n",
  )
  assert_eq(names.length(), 2)
  assert_eq(names[0], "seq_a")
  assert_eq(names[1], "seq_b")
  assert_eq(seqs[0], "MKTA")
  assert_eq(seqs[1], "MKAA")
}

///|
test "muscle local alignment" {
  let result = @src.biotite_muscle_local(["MKTA", "MKAA"], names=["a", "b"])
  assert_eq(result.status, "computed-local")
  assert_eq(result.aligned_sequences.length(), 2)
  assert_eq(
    result.aligned_sequences[0].length(),
    result.aligned_sequences[1].length(),
  )
  assert_true(result.command.contains("muscle"))
}

///|
test "muscle to fasta roundtrip" {
  let result = @src.biotite_muscle_local(["MKTA", "MKAA"], names=["a", "b"])
  let fasta = @src.biotite_muscle_to_fasta(result)
  let (names, seqs) = @src.biotite_muscle_parse_fasta(fasta)
  assert_eq(names[0], "a")
  assert_eq(names[1], "b")
  assert_eq(seqs[0], result.aligned_sequences[0])
  assert_eq(seqs[1], result.aligned_sequences[1])
}

///|
test "muscle identical sequences align unchanged" {
  let result = @src.biotite_muscle_local(["AAAA", "AAAA"])
  assert_eq(result.aligned_sequences[0], "AAAA")
  assert_eq(result.aligned_sequences[1], "AAAA")
}
