///|
/// Hidden Markov Model Demo - Gene Prediction
///
/// This demo shows how to use HMM for gene prediction and sequence analysis.
fn main {
  println("=== Hidden Markov Model (HMM) Demo ===")
  println("")

  println("1. Create gene prediction HMM:")
  let hmm = @src.create_gene_prediction_hmm()
  println("   States: \{hmm.states}")
  println("   Observations: \{hmm.observations}")
  println("")

  println("2. Forward algorithm - compute probability:")
  let obs = ["A", "T", "C", "G"]
  let prob = hmm.forward(obs)
  println("   Observation: \{obs}")
  println("   Probability: \{prob}")
  println("")

  println("3. Backward algorithm:")
  let beta = hmm.backward(obs)
  println("   Number of states: \{beta.length()}")
  println("")

  println("4. Viterbi algorithm - most probable state sequence:")
  let path = hmm.viterbi(obs)
  println("   Observation: \{obs}")
  println("   State path: \{path}")
  println("")

  println("5. Gene prediction on DNA sequence:")
  let dna = "AAATTTCCCGGGAAATTT"
  println("   DNA: \{dna}")
  let prediction = @src.predict_genes(hmm, dna)
  println("   Prediction: \{prediction}")
  println("")

  println("6. Extract exons from prediction:")
  let exons = @src.extract_exons(prediction, dna)
  println("   Exons found: \{exons.length()}")
  for i = 0; i < exons.length(); i = i + 1 {
    println("     Exon \{i}: \{exons[i]}")
  }
  println("")

  println("7. Baum-Welch training:")
  let training_obs = ["A", "T", "C", "G", "A", "T", "C", "G"]
  let trained_hmm = hmm.baum_welch(training_obs, 10, 0.000001)
  println("   Training observations: \{training_obs}")
  println("   Trained states: \{trained_hmm.states.length()}")
  println("")

  println("8. Predict on trained model:")
  let trained_prediction = @src.predict_genes(trained_hmm, "ATCG")
  println("   DNA: ATCG")
  println("   Prediction: \{trained_prediction}")
  println("")

  println("=== Demo Complete ===")
}
