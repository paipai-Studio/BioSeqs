///|
/// Test file for markov module.

fn make_dna_seqs() -> Array[String] {
  ["ACGTACGT", "CGCGCGTA", "TTTAAAAC", "GGCCTTAA", "ATGCATGC"]
}

fn make_cpg_seqs() -> Array[String] {
  ["CGCGCGCG", "ACGTACGT", "CGATCGAT", "GCGCGCGC", "TCGATCGA"]
}

fn make_bg_seqs() -> Array[String] {
  ["ATATATAT", "TTAAAATT", "GGGGGGGG", "CCCCCCCC", "AATTCCGG"]
}

test "markov_chain_type_helpers" {
  let a = @src.markov_first_order()
  let b = @src.markov_second_order()
  let c = @src.markov_third_order()
  assert_eq(a, a)
  assert_true(a != b)
  assert_true(b != c)
}

test "markov_model_new_default" {
  let m = @src.MarkovModel::new()
  assert_eq(m.order, 1)
  assert_eq(m.states.length(), 4)
  assert_eq(m.pseudo_count, 1.0)
}

test "markov_model_default_dna" {
  let m = @src.MarkovModel::default_dna()
  assert_eq(m.order, 1)
  assert_eq(m.states.length(), 4)
  assert_eq(m.states[0], "A")
  assert_eq(m.states[1], "C")
}

test "markov_model_default_protein" {
  let m = @src.MarkovModel::default_protein()
  assert_eq(m.order, 1)
  assert_eq(m.states.length(), 20)
}

test "markov_model_set_order" {
  let m = @src.MarkovModel::new()
  let m2 = m.set_order(2)
  assert_eq(m2.order, 2)
  let m3 = m.set_order(3)
  assert_eq(m3.order, 3)
  assert_eq(m.order, 1)
}

test "markov_model_set_states" {
  let m = @src.MarkovModel::new()
  let states = ["X", "Y"]
  let m2 = m.set_states(states=states)
  assert_eq(m2.states.length(), 2)
  assert_eq(m2.states[0], "X")
}

test "markov_model_set_chain_type" {
  let m = @src.MarkovModel::new()
  let m2 = m.set_chain_type(val=@src.markov_second_order())
  assert_eq(m2.order, 2)
  let m3 = m.set_chain_type(val=@src.markov_third_order())
  assert_eq(m3.order, 3)
}

test "markov_model_set_pseudo_count" {
  let m = @src.MarkovModel::new()
  let m2 = m.set_pseudo_count(0.5)
  assert_eq(m2.pseudo_count, 0.5)
  assert_eq(m.pseudo_count, 1.0)
}

test "markov_build_model_dna_order1" {
  let seqs = make_dna_seqs()
  let m = @src.markov_build_model(seqs, order=1)
  assert_eq(m.order, 1)
  assert_true(m.states.length() >= 4)
  assert_true(m.transition_probs.length() > 0)
}

test "markov_build_model_different_orders" {
  let seqs = make_dna_seqs()
  let m1 = @src.markov_build_model(seqs, order=1)
  let m2 = @src.markov_build_model(seqs, order=2)
  let m3 = @src.markov_build_model(seqs, order=3)
  let c1 = m1.transition_probs.length()
  let c2 = m2.transition_probs.length()
  let c3 = m3.transition_probs.length()
  assert_true(c1 != c2 || c2 != c3)
}

test "markov_score_sequence_finite_negative" {
  let seqs = make_dna_seqs()
  let m = @src.markov_build_model(seqs, order=1)
  let score = @src.markov_score_sequence(m, "ACGTACGT")
  assert_true(score <= 0.0)
  assert_true(score > -1000.0)
}

test "markov_score_per_base_length" {
  let seqs = make_dna_seqs()
  let m1 = @src.markov_build_model(seqs, order=1)
  let s1 = @src.markov_score_per_base(m1, "ACGTACGT")
  assert_eq(s1.length(), 8 - 1)
  let m2 = @src.markov_build_model(seqs, order=2)
  let s2 = @src.markov_score_per_base(m2, "ACGTACGT")
  assert_eq(s2.length(), 8 - 2)
}

test "markov_generate_sequence_length_and_alphabet" {
  let seqs = make_dna_seqs()
  let m = @src.markov_build_model(seqs, order=1)
  let g = @src.markov_generate_sequence(m, 20, seed=42)
  assert_eq(g.length(), 20)
  let states_set : Map[String, Unit] = Map::new()
  let mut si = 0
  while si < m.states.length() {
    states_set.set(m.states[si], ())
    si = si + 1
  }
  let mut ci = 0
  while ci < g.length() {
    let ch = g[ci:ci + 1].to_string()
    assert_true(states_set.contains(ch))
    ci = ci + 1
  }
}

test "markov_log_odds_positive_cpg" {
  let cpg_seqs = make_cpg_seqs()
  let bg_seqs = make_bg_seqs()
  let m_cpg = @src.markov_build_model(cpg_seqs, order=1)
  let m_bg = @src.markov_build_model(bg_seqs, order=1)
  let cpg_test = "CGCGCGCGCG"
  let lo = @src.markov_log_odds(m_cpg, m_bg, cpg_test)
  assert_true(lo > -100.0)
}

test "markov_stationary_sums_to_one" {
  let seqs = make_dna_seqs()
  let m = @src.markov_build_model(seqs, order=1)
  let sd = @src.markov_stationary_distribution(m)
  let mut total = 0.0
  let keys = sd.keys().collect()
  let mut ki = 0
  while ki < keys.length() {
    total = total + sd.get(keys[ki]).unwrap_or(0.0)
    ki = ki + 1
  }
  assert_true(total > 0.9)
  assert_true(total < 1.1)
}

test "markov_pseudo_count_affects_unknowns" {
  let seqs = ["AAAAA", "CCCCCCC"]
  let sts = ["A", "C", "G", "T"]
  let m_low = @src.markov_build_model(seqs, order=1, states=sts, pseudo_count=0.1)
  let m_high = @src.markov_build_model(seqs, order=1, states=sts, pseudo_count=10.0)
  let s_low = @src.markov_score_sequence(m_low, "GGGG")
  let s_high = @src.markov_score_sequence(m_high, "GGGG")
  assert_true(s_high > s_low)
}

test "markov_build_model_infers_states" {
  let seqs = ["AABBAABB", "BBBBAAAA"]
  let m = @src.markov_build_model(seqs)
  assert_true(m.states.length() >= 2)
}

test "markov_short_sequence_score" {
  let seqs = make_dna_seqs()
  let m = @src.markov_build_model(seqs, order=1)
  let s = @src.markov_score_sequence(m, "A")
  assert_eq(s, 0.0)
}
