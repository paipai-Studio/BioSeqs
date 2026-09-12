///|
/// Tests for the classic Bio.HMM package port (MarkovModel,
/// DynamicProgramming and Trainer). Expected values were generated with
/// the original Biopython 1.84 implementation.

///|
fn hmm_approx(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() < tol
}

///|
fn hmm_approx_vec(a : Array[Double], b : Array[Double], tol : Double) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i = 0; i < a.length(); i = i + 1 {
    if (a[i] - b[i]).abs() >= tol {
      return false
    }
  }
  true
}

///|
/// Build the CpG-style two-state model used by the Biopython tutorial:
/// states ('0', '1'), emissions ACGT, all transitions allowed.
fn hmm_build_cpg() -> @src.HiddenMarkovModel raise {
  let b = @src.MarkovModelBuilder::new(state_alphabet=["0", "1"], emission_alphabet=[
    "A", "C", "G", "T",
  ])
  b.allow_all_transitions()
  let initial : Map[String, Double] = Map([], capacity=2)
  initial.set("0", 1.0)
  b.set_initial_probabilities(initial)
  b.set_transition_score(from_state="0", to_state="0", probability=0.6)
  b.set_transition_score(from_state="0", to_state="1", probability=0.4)
  b.set_transition_score(from_state="1", to_state="0", probability=0.1)
  b.set_transition_score(from_state="1", to_state="1", probability=0.9)
  for s in ["A", "C", "G", "T"] {
    b.set_emission_score(seq_state="0", emission_state=s, probability=0.25)
  }
  b.set_emission_score(seq_state="1", emission_state="A", probability=0.125)
  b.set_emission_score(seq_state="1", emission_state="C", probability=0.375)
  b.set_emission_score(seq_state="1", emission_state="G", probability=0.375)
  b.set_emission_score(seq_state="1", emission_state="T", probability=0.125)
  b.get_markov_model()
}

///|
test "hmm_builder_defaults" {
  let b = @src.MarkovModelBuilder::new(state_alphabet=["0", "1"], emission_alphabet=[
    "A", "C", "G", "T",
  ])
  // by default: no transitions, all emissions allowed with zero
  // probability and a pseudo count of 1
  assert_eq(b.emission_pair_count(), 8)
  for s in ["0", "1"] {
    for e in ["A", "C", "G", "T"] {
      assert_true(hmm_approx(b.emission_score(s, e), 0.0, 1.0e-12))
      assert_true(hmm_approx(b.emission_pseudocount(s, e), 1.0, 1.0e-12))
    }
  }
  assert_eq(b.transition_count(), 0)
}

///|
test "hmm_set_initial_probabilities_residual" {
  let b = @src.MarkovModelBuilder::new(state_alphabet=["A", "B"], emission_alphabet=[
    "x", "y",
  ])
  // fully specified
  let p1 : Map[String, Double] = Map([], capacity=2)
  p1.set("A", 1.0)
  b.set_initial_probabilities(p1)
  assert_true(hmm_approx(b.initial_probability("A"), 1.0, 1.0e-12))
  assert_true(hmm_approx(b.initial_probability("B"), 0.0, 1.0e-12))
  // residual is distributed evenly over the unset states
  let p2 : Map[String, Double] = Map([], capacity=2)
  p2.set("A", 0.5)
  b.set_initial_probabilities(p2)
  assert_true(hmm_approx(b.initial_probability("B"), 0.5, 1.0e-12))
  // empty distribution gives equal probabilities
  let p3 : Map[String, Double] = Map([], capacity=2)
  b.set_initial_probabilities(p3)
  assert_true(hmm_approx(b.initial_probability("A"), 0.5, 1.0e-12))
  assert_true(hmm_approx(b.initial_probability("B"), 0.5, 1.0e-12))
}

///|
test "hmm_set_initial_probabilities_errors" {
  let b = @src.MarkovModelBuilder::new(state_alphabet=["A", "B"], emission_alphabet=[
    "x", "y",
  ])
  // unknown state
  let bad : Map[String, Double] = Map([], capacity=2)
  bad.set("Z", 1.0)
  let r1 = try {
    b.set_initial_probabilities(bad)
    false
  } catch {
    HmmError(_) => true
  }
  assert_true(r1)
  // total probability above 1
  let over : Map[String, Double] = Map([], capacity=2)
  over.set("A", 0.7)
  over.set("B", 0.7)
  let r2 = try {
    b.set_initial_probabilities(over)
    false
  } catch {
    HmmError(_) => true
  }
  assert_true(r2)
}

///|
test "hmm_get_markov_model_requires_initial_probabilities" {
  let b = @src.MarkovModelBuilder::new(state_alphabet=["A", "B"], emission_alphabet=[
    "x", "y",
  ])
  let r = try {
    ignore(b.get_markov_model())
    false
  } catch {
    HmmError(_) => true
  }
  assert_true(r)
}

///|
test "hmm_transitions_management" {
  let b = @src.MarkovModelBuilder::new(state_alphabet=["A", "B"], emission_alphabet=[
    "x", "y",
  ])
  b.allow_transition(from_state="A", to_state="B")
  assert_eq(b.transition_count(), 1)
  // duplicate allow raises
  let r1 = try {
    b.allow_transition(from_state="A", to_state="B")
    false
  } catch {
    HmmError(_) => true
  }
  assert_true(r1)
  // default probability is 0 and pseudo count is 1
  assert_true(hmm_approx(b.transition_score("A", "B"), 0.0, 1.0e-12))
  assert_true(hmm_approx(b.transition_pseudocount("A", "B"), 1.0, 1.0e-12))
  b.set_transition_score(from_state="A", to_state="B", probability=0.75)
  assert_true(hmm_approx(b.transition_score("A", "B"), 0.75, 1.0e-12))
  b.set_transition_pseudocount(from_state="A", to_state="B", count=2.0)
  assert_true(hmm_approx(b.transition_pseudocount("A", "B"), 2.0, 1.0e-12))
  b.destroy_transition(from_state="A", to_state="B")
  assert_eq(b.transition_count(), 0)
  // destroying a disallowed transition raises
  let r2 = try {
    b.destroy_transition(from_state="A", to_state="B")
    false
  } catch {
    HmmError(_) => true
  }
  assert_true(r2)
  // unknown states raise
  let r3 = try {
    b.allow_transition(from_state="A", to_state="Z")
    false
  } catch {
    HmmError(_) => true
  }
  assert_true(r3)
}

///|
test "hmm_allow_all_transitions_preserves_scores" {
  let b = @src.MarkovModelBuilder::new(state_alphabet=["A", "B"], emission_alphabet=[
    "x", "y",
  ])
  b.allow_transition(from_state="A", to_state="B", probability=0.5)
  b.allow_all_transitions()
  assert_eq(b.transition_count(), 4)
  assert_true(hmm_approx(b.transition_score("A", "B"), 0.5, 1.0e-12))
  assert_true(hmm_approx(b.transition_score("B", "A"), 0.0, 1.0e-12))
}

///|
test "hmm_set_equal_probabilities" {
  // reference: Biopython 1.84 behaviour, including the original quirk
  // that initial probabilities are divided by the number of transitions
  let b = @src.MarkovModelBuilder::new(state_alphabet=["A", "B"], emission_alphabet=[
    "x", "y",
  ])
  b.allow_transition(from_state="A", to_state="B")
  b.allow_transition(from_state="B", to_state="A")
  b.set_equal_probabilities()
  assert_true(hmm_approx(b.initial_probability("A"), 0.5, 1.0e-12))
  assert_true(hmm_approx(b.initial_probability("B"), 0.5, 1.0e-12))
  assert_true(hmm_approx(b.transition_score("A", "B"), 0.5, 1.0e-12))
  assert_true(hmm_approx(b.transition_score("B", "A"), 0.5, 1.0e-12))
  assert_true(hmm_approx(b.emission_score("A", "x"), 0.25, 1.0e-12))
  assert_true(hmm_approx(b.emission_score("B", "y"), 0.25, 1.0e-12))
}

///|
test "hmm_set_random_probabilities_deterministic" {
  let b = @src.MarkovModelBuilder::new(state_alphabet=["A", "B"], emission_alphabet=[
    "x", "y",
  ])
  b.allow_all_transitions()
  @src.hmm_seed(42)
  b.set_random_probabilities()
  let init_sum = b.initial_probability("A") + b.initial_probability("B")
  assert_true(hmm_approx(init_sum, 1.0, 1.0e-9))
  let trans_sum = b.transition_score("A", "A") + b.transition_score("A", "B")
  assert_true(hmm_approx(trans_sum, 1.0, 1.0e-9))
  let emis_sum = b.emission_score("A", "x") + b.emission_score("A", "y")
  assert_true(hmm_approx(emis_sum, 1.0, 1.0e-9))
  // reproducible with the same seed
  let b2 = @src.MarkovModelBuilder::new(state_alphabet=["A", "B"], emission_alphabet=[
    "x", "y",
  ])
  b2.allow_all_transitions()
  @src.hmm_seed(42)
  b2.set_random_probabilities()
  assert_true(
    hmm_approx(
      b2.transition_score("A", "A"),
      b.transition_score("A", "A"),
      1.0e-30,
    ),
  )
  assert_true(
    hmm_approx(b2.emission_score("B", "y"), b.emission_score("B", "y"), 1.0e-30),
  )
  // random transitions without allowed transitions raise
  let b3 = @src.MarkovModelBuilder::new(state_alphabet=["A", "B"], emission_alphabet=[
    "x", "y",
  ])
  let r = try {
    let _ = b3.set_random_transition_probabilities()
    false
  } catch {
    HmmError(_) => true
  }
  assert_true(r)
}

///|
test "hmm_viterbi_cpg" {
  let mm = hmm_build_cpg()
  let result = mm.viterbi(["A", "C", "G", "C", "G", "A", "T", "C"])
  assert_eq(result.state_path(), ["0", "1", "1", "1", "1", "1", "1", "1"])
  assert_true(hmm_approx(result.prob(), -11.997777535359305, 1.0e-9))
}

///|
test "hmm_viterbi_zero_probability_transition" {
  // a transition with score 0 is avoided by the log transform
  let b = @src.MarkovModelBuilder::new(state_alphabet=["A", "B"], emission_alphabet=[
    "x", "y",
  ])
  b.allow_transition(from_state="A", to_state="A")
  b.allow_transition(from_state="A", to_state="B")
  b.allow_transition(from_state="B", to_state="A")
  b.allow_transition(from_state="B", to_state="B")
  b.set_transition_score(from_state="A", to_state="A", probability=0.5)
  b.set_transition_score(from_state="A", to_state="B", probability=0.5)
  b.set_transition_score(from_state="B", to_state="A", probability=1.0)
  b.set_transition_score(from_state="B", to_state="B", probability=0.0)
  for s in ["x", "y"] {
    b.set_emission_score(seq_state="A", emission_state=s, probability=0.5)
    b.set_emission_score(seq_state="B", emission_state=s, probability=0.5)
  }
  let initial : Map[String, Double] = Map([], capacity=2)
  initial.set("A", 1.0)
  b.set_initial_probabilities(initial)
  let mm = b.get_markov_model()
  let result = mm.viterbi(["x", "x", "y"])
  assert_eq(result.state_path(), ["A", "B", "A"])
  assert_true(hmm_approx(result.prob(), -2.772588722239781, 1.0e-9))
}

///|
test "hmm_viterbi_no_incoming_transitions_raises" {
  // state B has no incoming transitions: the original raises ValueError
  let b = @src.MarkovModelBuilder::new(state_alphabet=["A", "B"], emission_alphabet=[
    "x", "y",
  ])
  b.allow_transition(from_state="A", to_state="A")
  let initial : Map[String, Double] = Map([], capacity=2)
  initial.set("A", 1.0)
  b.set_initial_probabilities(initial)
  let mm = b.get_markov_model()
  let r = try {
    ignore(mm.viterbi(["x", "y"]))
    false
  } catch {
    HmmError(_) => true
  }
  assert_true(r)
}

///|
test "hmm_training_sequence_validation" {
  let r = try {
    ignore(@src.TrainingSequence::new(emissions=["x", "y"], state_path=["A"]))
    false
  } catch {
    HmmError(_) => true
  }
  assert_true(r)
  // empty state path is allowed for unlabelled training
  let ts = @src.TrainingSequence::new(emissions=["x", "y"], state_path=[])
  assert_eq(ts.emissions().length(), 2)
  assert_eq(ts.states().length(), 0)
}

///|
test "hmm_forward_algorithm_cpg" {
  let mm = hmm_build_cpg()
  let ts = @src.TrainingSequence::new(
    emissions=["A", "C", "G", "C", "G", "A", "T", "C"],
    state_path=[],
  )
  let dp = @src.ScaledDp::new(markov_model=mm, sequence=ts)
  let fwd = dp.forward_algorithm()
  assert_true(hmm_approx(fwd.seq_prob(), 0.20820867050392472, 1.0e-12))
  assert_true(
    hmm_approx_vec(
      [fwd.forward_vars().get("0", 0), fwd.forward_vars().get("1", 0)],
      [0.75, 0.25],
      1.0e-12,
    ),
  )
  assert_true(
    hmm_approx_vec(
      [fwd.forward_vars().get("0", 7), fwd.forward_vars().get("1", 7)],
      [0.21641734100784948, 0.7835826589921504],
      1.0e-12,
    ),
  )
}

///|
test "hmm_backward_algorithm_cpg" {
  let mm = hmm_build_cpg()
  let ts = @src.TrainingSequence::new(
    emissions=["A", "C", "G", "C", "G", "A", "T", "C"],
    state_path=[],
  )
  let dp = @src.ScaledDp::new(markov_model=mm, sequence=ts)
  let _ = dp.forward_algorithm()
  let bwd = dp.backward_algorithm()
  assert_true(
    hmm_approx_vec(
      [bwd.get("0", 0), bwd.get("1", 0)],
      [0.23374849786929278, 0.446281131809817],
      1.0e-12,
    ),
  )
  assert_true(
    hmm_approx_vec(
      [bwd.get("0", 6), bwd.get("1", 6)],
      [0.6456571220550653, 0.3631821311559743],
      1.0e-12,
    ),
  )
  assert_true(
    hmm_approx_vec([bwd.get("0", 7), bwd.get("1", 7)], [0.6, 0.1], 1.0e-12),
  )
}

///|
test "hmm_baum_welch_cpg" {
  let mm = hmm_build_cpg()
  let t1 = @src.TrainingSequence::new(
    emissions=["A", "C", "G", "C", "G", "A", "T", "C"],
    state_path=[],
  )
  let t2 = @src.TrainingSequence::new(
    emissions=["G", "G", "C", "T", "A", "T", "A"],
    state_path=[],
  )
  let trainer = @src.BaumWelchTrainer::new(markov_model=mm)
  let trained = trainer.train([t1, t2], fn(change : Double, iterations : Int) {
    change < 1.0e-12 || iterations >= 100
  })
  // reference values from the original Biopython implementation
  assert_true(
    hmm_approx(
      trained.get_transition_probability("0", "0"),
      0.7269202953382092,
      1.0e-9,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_transition_probability("0", "1"),
      0.2730797046617907,
      1.0e-9,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_transition_probability("1", "0"),
      0.4627381457672509,
      1.0e-9,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_transition_probability("1", "1"),
      0.5372618542327491,
      1.0e-9,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("0", "A"),
      0.3919273904378798,
      1.0e-9,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("0", "C"),
      0.16357350192453832,
      1.0e-9,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("0", "G"),
      0.08791659993166538,
      1.0e-9,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("0", "T"),
      0.3565825077059167,
      1.0e-9,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("1", "A"),
      0.05842511383179566,
      1.0e-9,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("1", "C"),
      0.37346131792199144,
      1.0e-9,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("1", "G"),
      0.5205622956341125,
      1.0e-9,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("1", "T"),
      0.04755127261210037,
      1.0e-9,
    ),
  )
}

///|
test "hmm_known_state_trainer_cpg" {
  let mm = hmm_build_cpg()
  let ts = @src.TrainingSequence::new(
    emissions=["A", "C", "G", "C", "G", "A", "T", "C"],
    state_path=["0", "0", "1", "1", "1", "1", "1", "0"],
  )
  let trainer = @src.KnownStateTrainer::new(markov_model=mm)
  let trained = trainer.train([ts])
  assert_true(
    hmm_approx(trained.get_transition_probability("0", "0"), 0.5, 1.0e-12),
  )
  assert_true(
    hmm_approx(trained.get_transition_probability("0", "1"), 0.5, 1.0e-12),
  )
  assert_true(
    hmm_approx(
      trained.get_transition_probability("1", "0"),
      0.2857142857142857,
      1.0e-12,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_transition_probability("1", "1"),
      0.7142857142857143,
      1.0e-12,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("0", "A"),
      0.2857142857142857,
      1.0e-12,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("0", "C"),
      0.42857142857142855,
      1.0e-12,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("0", "G"),
      0.14285714285714285,
      1.0e-12,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("0", "T"),
      0.14285714285714285,
      1.0e-12,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("1", "A"),
      0.2222222222222222,
      1.0e-12,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("1", "C"),
      0.2222222222222222,
      1.0e-12,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("1", "G"),
      0.3333333333333333,
      1.0e-12,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("1", "T"),
      0.2222222222222222,
      1.0e-12,
    ),
  )
}

///|
test "hmm_known_state_trainer_unexpected_transition" {
  let b = @src.MarkovModelBuilder::new(state_alphabet=["A", "B"], emission_alphabet=[
    "x", "y",
  ])
  b.allow_transition(from_state="A", to_state="B")
  let initial : Map[String, Double] = Map([], capacity=2)
  initial.set("A", 1.0)
  b.set_initial_probabilities(initial)
  let mm = b.get_markov_model()
  let trainer = @src.KnownStateTrainer::new(markov_model=mm)
  // the path contains the disallowed A -> A transition
  let ts = @src.TrainingSequence::new(emissions=["x", "y", "x", "y"], state_path=[
    "A", "A", "B", "B",
  ])
  let r = try {
    ignore(trainer.train([ts]))
    false
  } catch {
    HmmError(_) => true
  }
  assert_true(r)
}

///|
test "hmm_known_state_trainer_pseudocount_smoothing" {
  // emissions not observed in training are smoothed by the pseudo counts
  let b = @src.MarkovModelBuilder::new(state_alphabet=["A", "B"], emission_alphabet=[
    "x", "y",
  ])
  b.allow_transition(from_state="A", to_state="B")
  b.allow_transition(from_state="B", to_state="B")
  let initial : Map[String, Double] = Map([], capacity=2)
  initial.set("A", 1.0)
  b.set_initial_probabilities(initial)
  let mm = b.get_markov_model()
  let trainer = @src.KnownStateTrainer::new(markov_model=mm)
  let ts = @src.TrainingSequence::new(emissions=["x", "y", "x"], state_path=[
    "A", "B", "B",
  ])
  let trained = trainer.train([ts])
  assert_true(
    hmm_approx(trained.get_transition_probability("A", "B"), 1.0, 1.0e-12),
  )
  assert_true(
    hmm_approx(trained.get_transition_probability("B", "B"), 1.0, 1.0e-12),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("A", "x"),
      0.6666666666666666,
      1.0e-12,
    ),
  )
  assert_true(
    hmm_approx(
      trained.get_emission_probability("A", "y"),
      0.3333333333333333,
      1.0e-12,
    ),
  )
  assert_true(
    hmm_approx(trained.get_emission_probability("B", "x"), 0.5, 1.0e-12),
  )
  assert_true(
    hmm_approx(trained.get_emission_probability("B", "y"), 0.5, 1.0e-12),
  )
}

///|
test "hmm_transitions_from_to" {
  let b = @src.MarkovModelBuilder::new(state_alphabet=["A", "B"], emission_alphabet=[
    "x", "y",
  ])
  b.allow_transition(from_state="B", to_state="A")
  b.allow_transition(from_state="A", to_state="B")
  let initial : Map[String, Double] = Map([], capacity=2)
  initial.set("A", 1.0)
  b.set_initial_probabilities(initial)
  let mm = b.get_markov_model()
  // insertion order is preserved like the Python dict
  assert_eq(mm.transitions_from("B"), ["A"])
  assert_eq(mm.transitions_from("A"), ["B"])
  assert_eq(mm.transitions_to("A"), ["B"])
  assert_eq(mm.transitions_to("B"), ["A"])
  assert_eq(mm.transitions_from("Z"), [])
  assert_eq(mm.transitions_to("Z"), [])
}
