///|
/// Tests for the auto-generated built-in substitution matrix name lists.

///|
test "substitution_matrices_model_names" {
  @src.subs_ensure_registry()
  let dna = @src.bio_dna_model_names()
  let protein = @src.bio_protein_model_names()
  let dist = @src.bio_distance_model_names()
  // 非空
  assert_true(dna.length() > 0)
  assert_true(protein.length() > 0)
  // 已知模型名
  assert_true(dna.contains("blastn"))
  assert_true(protein.contains("blosum62"))
  assert_true(dna.contains("benner22"))
  // distance = ["identity"] + dna + protein
  assert_eq(dist[0], "identity")
  assert_eq(dist.length(), 1 + dna.length() + protein.length())
  // 顺序一致
  assert_eq(dist[1], dna[0])
  assert_eq(dist[1 + dna.length()], protein[0])
}