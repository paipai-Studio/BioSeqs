///|
/// Bio.HMM (classic Markov model) Demo.
///
/// Demonstrates the port of Biopython's classic Bio.HMM package:
/// MarkovModelBuilder, Viterbi decoding, the scaled forward/backward
/// algorithms, and the Baum-Welch / known-state trainers.

///|
fn main raise {
  // A two-state CpG-style model: state "0" = neutral, state "1" = CpG island.
  let builder = @src.MarkovModelBuilder::new(
    state_alphabet=["0", "1"],
    emission_alphabet=["A", "C", "G", "T"],
  )
  builder.allow_all_transitions()
  let initial : Map[String, Double] = Map([], capacity=2)
  initial.set("0", 1.0)
  builder.set_initial_probabilities(initial)
  builder.set_transition_score(from_state="0", to_state="0", probability=0.6)
  builder.set_transition_score(from_state="0", to_state="1", probability=0.4)
  builder.set_transition_score(from_state="1", to_state="0", probability=0.1)
  builder.set_transition_score(from_state="1", to_state="1", probability=0.9)
  for s in ["A", "C", "G", "T"] {
    builder.set_emission_score(seq_state="0", emission_state=s, probability=0.25)
  }
  builder.set_emission_score(seq_state="1", emission_state="A", probability=0.125)
  builder.set_emission_score(seq_state="1", emission_state="C", probability=0.375)
  builder.set_emission_score(seq_state="1", emission_state="G", probability=0.375)
  builder.set_emission_score(seq_state="1", emission_state="T", probability=0.125)
  let model = builder.get_markov_model()

  println("=== 1. MarkovModelBuilder ===")
  println("states:    \{model.get_state_alphabet()}")
  println("emissions: \{model.get_emission_alphabet()}")
  println(
    "transition 0 -> 1: \{model.get_transition_probability("0", "1")}",
  )
  println(
    "emission 1 -> G:   \{model.get_emission_probability("1", "G")}",
  )

  println("\n=== 2. Viterbi decoding ===")
  let sequence = ["A", "C", "G", "C", "G", "A", "T", "C"]
  let decoded = model.viterbi(sequence)
  println("sequence:   \{sequence}")
  println("state path: \{decoded.state_path()}")
  println("log prob:   \{decoded.prob()}")

  println("\n=== 3. Scaled forward algorithm ===")
  let training_seq = @src.TrainingSequence::new(
    emissions=sequence,
    state_path=[],
  )
  let dp = @src.ScaledDp::new(markov_model=model, sequence=training_seq)
  let forward = dp.forward_algorithm()
  println("P(sequence) = \{forward.seq_prob()}")
  let f0 = forward.forward_vars().get("0", 0)
  let f1 = forward.forward_vars().get("1", 0)
  println("scaled forward variables at position 0: [\{f0}, \{f1}]")

  println("\n=== 4. Scaled backward algorithm ===")
  let backward = dp.backward_algorithm()
  let b_last0 = backward.get("0", 7)
  let b_last1 = backward.get("1", 7)
  println("scaled backward variables at position 7: [\{b_last0}, \{b_last1}]")

  println("\n=== 5. Baum-Welch training (unknown state path) ===")
  let t1 = @src.TrainingSequence::new(
    emissions=["A", "C", "G", "C", "G", "A", "T", "C"],
    state_path=[],
  )
  let t2 = @src.TrainingSequence::new(
    emissions=["G", "G", "C", "T", "A", "T", "A"],
    state_path=[],
  )
  let trainer = @src.BaumWelchTrainer::new(markov_model=model)
  let trained = trainer.train([t1, t2], fn(change : Double, iterations : Int) {
    change < 1.0e-12 || iterations >= 100
  })
  println(
    "trained transition 0 -> 1: \{trained.get_transition_probability("0", "1")}",
  )
  println(
    "trained transition 1 -> 1: \{trained.get_transition_probability("1", "1")}",
  )
  println(
    "trained emission 1 -> G:   \{trained.get_emission_probability("1", "G")}",
  )

  println("\n=== 6. Known-state training ===")
  let known = @src.TrainingSequence::new(
    emissions=["x", "y", "x"],
    state_path=["A", "B", "B"],
  )
  // a small two-letter model to show state-path supervised training
  let known_builder = @src.MarkovModelBuilder::new(
    state_alphabet=["A", "B"],
    emission_alphabet=["x", "y"],
  )
  known_builder.allow_transition(from_state="A", to_state="B")
  known_builder.allow_transition(from_state="B", to_state="B")
  let known_initial : Map[String, Double] = Map([], capacity=2)
  known_initial.set("A", 1.0)
  known_builder.set_initial_probabilities(known_initial)
  let known_model = known_builder.get_markov_model()
  let known_trainer = @src.KnownStateTrainer::new(markov_model=known_model)
  let known_trained = known_trainer.train([known])
  println(
    "trained transition A -> B: \{known_trained.get_transition_probability("A", "B")}",
  )
  println(
    "trained emission A -> x:   \{known_trained.get_emission_probability("A", "x")}",
  )

  println("\nDone.")
}
