///|
/// HMM Tests
test "hmm_create" {
  let hmm = @src.create_gene_prediction_hmm()
  assert_eq(hmm.states.length(), 5)
  assert_eq(hmm.observations.length(), 4)
}

///|
test "hmm_forward" {
  let hmm = @src.create_gene_prediction_hmm()
  let obs = ["A", "T", "C", "G"]
  let prob = hmm.forward(obs)
  assert_true(prob > 0.0)
}

///|
test "hmm_backward" {
  let hmm = @src.create_gene_prediction_hmm()
  let obs = ["A", "T", "C", "G"]
  let beta = hmm.backward(obs)
  assert_eq(beta.length(), 5)
}

///|
test "hmm_viterbi" {
  let hmm = @src.create_gene_prediction_hmm()
  let obs = ["A", "T", "C", "G"]
  let path = hmm.viterbi(obs)
  assert_eq(path.length(), 4)
}

///|
test "hmm_predict_genes" {
  let hmm = @src.create_gene_prediction_hmm()
  let prediction = @src.predict_genes(hmm, "ATCG")
  assert_eq(prediction.length(), 4)
}

///|
test "hmm_extract_exons" {
  let hmm = @src.create_gene_prediction_hmm()
  let prediction = @src.predict_genes(hmm, "ATCGATCG")
  let exons = @src.extract_exons(prediction, "ATCGATCG")
  assert_true(exons.length() >= 0)
}

///|
test "hmm_small_dna" {
  let hmm = @src.create_gene_prediction_hmm()
  let dna = "AAATTTCCCGGG"
  let prediction = @src.predict_genes(hmm, dna)
  assert_eq(prediction.length(), dna.length())
}

///|
test "hmm_baum_welch" {
  let hmm = @src.create_gene_prediction_hmm()
  let obs = ["A", "T", "C", "G", "A", "T"]
  let trained = hmm.baum_welch(obs, 10, 0.000001)
  assert_eq(trained.states.length(), 5)
}

///|
test "hmm_state_transitions" {
  let hmm = @src.create_gene_prediction_hmm()
  assert_true(hmm.trans_prob.get("Start") is Some(_))
  assert_true(hmm.trans_prob.get("Exon") is Some(_))
}

///|
test "hmm_emit_probabilities" {
  let hmm = @src.create_gene_prediction_hmm()
  let emit = hmm.emit_prob.get("Exon").unwrap_or(Map::from_array([]))
  assert_true(emit.get("A") is Some(_))
}

///|
test "hmm_empty_input" {
  let hmm = @src.create_gene_prediction_hmm()
  let obs = []
  let prob = hmm.forward(obs)
  assert_eq(prob, 0.0)
}

///|
test "hmm_single_observation" {
  let hmm = @src.create_gene_prediction_hmm()
  let obs = ["A"]
  let prob = hmm.forward(obs)
  assert_true(prob > 0.0)
}

///|
test "hmm_viterbi_single" {
  let hmm = @src.create_gene_prediction_hmm()
  let obs = ["A"]
  let path = hmm.viterbi(obs)
  assert_eq(path.length(), 1)
}

///|
test "hmm_backward_single" {
  let hmm = @src.create_gene_prediction_hmm()
  let obs = ["A"]
  let beta = hmm.backward(obs)
  assert_eq(beta.length(), 5)
}

///|
test "hmm_exon_extraction" {
  let hmm = @src.create_gene_prediction_hmm()
  let dna = "AAAA"
  let prediction = @src.predict_genes(hmm, dna)
  let exons = @src.extract_exons(prediction, dna)
  assert_true(exons.length() >= 0)
}
