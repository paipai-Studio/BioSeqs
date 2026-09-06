// Bio.HMM - Hidden Markov Model code.
//
// Port of the classic Biopython Bio.HMM package (MarkovModel,
// DynamicProgramming and Trainer modules, last shipped with real code in
// Biopython 1.84 before the module was retired). It provides:
//   - MarkovModelBuilder: build a HiddenMarkovModel with explicit states,
//     emissions, transitions and pseudo counts.
//   - HiddenMarkovModel::viterbi: most probable state path decoding.
//   - ScaledDp: scaled forward/backward algorithms (Durbin et al. p78).
//   - BaumWelchTrainer: EM parameter estimation from unlabelled sequences.
//   - KnownStateTrainer: maximum likelihood estimation from labelled paths.
//
// Reference: Durbin et al., "Biological Sequence Analysis" (1998).

///|
/// Error raised by the hidden Markov model implementation when a state,
/// transition or emission is not allowed by the model, or when the model
/// has not been fully initialized.
pub suberror HmmError {
  HmmError(String)
}

///|
fn[T] hmm_fail(msg : String) -> T raise HmmError {
  raise HmmError(msg)
}

// ---------------------------------------------------------------------------
// Random helpers
// ---------------------------------------------------------------------------

///|
let hmm_rng_state : Ref[Int] = Ref(123456789)

///|
/// Seed the internal random number generator used by the builder's
/// set_random_* methods. The generator is deterministic (Park-Miller).
pub fn hmm_seed(seed : Int) -> Unit {
  hmm_rng_state.val = if seed <= 0 { 123456789 } else { seed % 2147483647 }
}

///|
/// Uniform random number in (0, 1).
fn hmm_uniform() -> Double {
  let a = 48271
  let m = 2147483647
  let q = 44488
  let r = 3399
  let s = hmm_rng_state.val
  let hi = s / q
  let lo = s % q
  let mut t = a * lo - r * hi
  if t <= 0 {
    t = t + m
  }
  hmm_rng_state.val = t
  t.to_double() / m.to_double()
}

///|
/// Return an array of n random numbers summing to 1.0.
fn hmm_gen_random_array(n : Int) -> Array[Double] {
  let rand_array : Array[Double] = Array::new(capacity=n)
  let mut total = 0.0
  for i = 0; i < n; i = i + 1 {
    let v = hmm_uniform()
    rand_array.push(v)
    total = total + v
  }
  for i = 0; i < n; i = i + 1 {
    rand_array[i] = rand_array[i] / total
  }
  rand_array
}

// ---------------------------------------------------------------------------
// Two-level probability map helpers
// ---------------------------------------------------------------------------

///|
fn hmm_map2() -> Map[String, Map[String, Double]] {
  Map([], capacity=16)
}

///|
/// Get an entry from a two-level map, raising when the key is absent
/// (mirrors the Python KeyError).
fn hmm_get2(
  m : Map[String, Map[String, Double]],
  k1 : String,
  k2 : String,
) -> Double raise HmmError {
  match m.get(k1) {
    Some(inner) =>
      match inner.get(k2) {
        Some(v) => v
        None => hmm_fail("Missing probability key (\{k1}, \{k2})")
      }
    None => hmm_fail("Missing probability key (\{k1}, \{k2})")
  }
}

///|
fn hmm_contains2(
  m : Map[String, Map[String, Double]],
  k1 : String,
  k2 : String,
) -> Bool {
  match m.get(k1) {
    Some(inner) => inner.contains(k2)
    None => false
  }
}

///|
fn hmm_set2(
  m : Map[String, Map[String, Double]],
  k1 : String,
  k2 : String,
  v : Double,
) -> Unit {
  match m.get(k1) {
    Some(inner) => inner.set(k2, v)
    None => {
      let inner : Map[String, Double] = Map([], capacity=8)
      inner.set(k2, v)
      m.set(k1, inner)
    }
  }
}

///|
fn hmm_remove2(
  m : Map[String, Map[String, Double]],
  k1 : String,
  k2 : String,
) -> Unit {
  match m.get(k1) {
    Some(inner) => inner.remove(k2)
    None => ()
  }
}

///|
fn hmm_copy2(
  m : Map[String, Map[String, Double]],
) -> Map[String, Map[String, Double]] {
  let out = hmm_map2()
  for k1, inner in m {
    let copy : Map[String, Double] = Map([], capacity=8)
    for k2, v in inner {
      copy.set(k2, v)
    }
    out.set(k1, copy)
  }
  out
}

///|
fn hmm_map2_size(m : Map[String, Map[String, Double]]) -> Int {
  let mut n = 0
  for _k1, inner in m {
    n = n + inner.length()
  }
  n
}

// ---------------------------------------------------------------------------
// MarkovModelBuilder
// ---------------------------------------------------------------------------

///|
/// Interface to build up a Markov Model.
///
/// This class is designed to try to separate the task of specifying the
/// Markov Model from the actual model itself, so this builder class should
/// be used to create Markov models instead of trying to initialize a
/// HiddenMarkovModel directly.
pub struct MarkovModelBuilder {
  priv state_alphabet : Array[String]
  priv emission_alphabet : Array[String]
  priv initial_prob : Map[String, Double]
  priv mut transition_prob : Map[String, Map[String, Double]]
  priv mut transition_pseudo : Map[String, Map[String, Double]]
  // Allowed transitions in insertion order (Python dicts keep insertion
  // order; an explicit list reproduces that behaviour deterministically).
  priv mut transition_order : Array[(String, String)]
  priv emission_prob : Map[String, Map[String, Double]]
  priv emission_pseudo : Map[String, Map[String, Double]]
}

///|
/// The default pseudo counts to use.
let hmm_default_pseudo : Double = 1.0

///|
/// Initialize a builder to create Markov Models.
///
/// - state_alphabet: all of the letters that can appear in the states.
/// - emission_alphabet: all of the letters that can be emitted by the HMM.
pub fn MarkovModelBuilder::new(
  state_alphabet~ : Array[String],
  emission_alphabet~ : Array[String],
) -> MarkovModelBuilder {
  let states : Array[String] = Array::new(capacity=state_alphabet.length())
  for s in state_alphabet {
    states.push(s)
  }
  let emissions : Array[String] = Array::new(capacity=emission_alphabet.length())
  for s in emission_alphabet {
    emissions.push(s)
  }
  let emission_prob = hmm_map2()
  let emission_pseudo = hmm_map2()
  for s in states {
    for e in emissions {
      hmm_set2(emission_prob, s, e, 0.0)
      hmm_set2(emission_pseudo, s, e, hmm_default_pseudo)
    }
  }
  {
    state_alphabet: states,
    emission_alphabet: emissions,
    initial_prob: Map([], capacity=8),
    transition_prob: hmm_map2(),
    transition_pseudo: hmm_map2(),
    transition_order: [],
    emission_prob,
    emission_pseudo,
  }
}

///|
/// Return the markov model corresponding with the current parameters.
///
/// Each markov model returned by a call to this function is unique
/// (ie. they don't influence each other).
pub fn MarkovModelBuilder::get_markov_model(
  self : MarkovModelBuilder,
) -> HiddenMarkovModel raise HmmError {
  // user must set initial probabilities
  if self.initial_prob.length() == 0 {
    hmm_fail(
      "set_initial_probabilities must be called to fully initialize the Markov model",
    )
  }
  let initial = Map([], capacity=8)
  for k, v in self.initial_prob {
    initial.set(k, v)
  }
  let order : Array[(String, String)] = Array::new(capacity=self.transition_order.length())
  for key in self.transition_order {
    order.push(key)
  }
  HiddenMarkovModel::build(
    self.state_alphabet,
    self.emission_alphabet,
    initial,
    hmm_copy2(self.transition_prob),
    hmm_copy2(self.emission_prob),
    hmm_copy2(self.transition_pseudo),
    hmm_copy2(self.emission_pseudo),
    order,
  )
}

///|
/// Set initial state probabilities.
///
/// initial_prob is a map mapping states to probabilities. If initial
/// probabilities are set for all states then they should add up to 1.
/// Otherwise the sum should be <= 1; the residual probability is divided
/// up evenly between all the states for which the initial probability has
/// not been set. For example, with state alphabet ('A', 'B') calling
/// set_initial_prob({}) results in P('A') = 0.5 and P('B') = 0.5.
pub fn MarkovModelBuilder::set_initial_probabilities(
  self : MarkovModelBuilder,
  initial_prob : Map[String, Double],
) -> Unit raise HmmError {
  // ensure that all referenced states are valid
  for state, _v in initial_prob {
    if !(self.state_alphabet.contains(state)) {
      hmm_fail("State \{state} was not found in the sequence alphabet")
    }
  }
  // replace any previous probabilities (original: self.initial_prob = copy)
  self.initial_prob.clear()
  for state, v in initial_prob {
    self.initial_prob.set(state, v)
  }
  // distribute the residual probability, if any
  let num_states_not_set = self.state_alphabet.length() - self.initial_prob.length()
  if num_states_not_set < 0 {
    hmm_fail("Initial probabilities can't exceed # of states")
  }
  let mut prob_sum = 0.0
  for _state, v in self.initial_prob {
    prob_sum = prob_sum + v
  }
  if prob_sum > 1.0 {
    hmm_fail("Total initial probability cannot exceed 1.0")
  }
  if num_states_not_set > 0 {
    let prob = (1.0 - prob_sum) / num_states_not_set.to_double()
    for state in self.state_alphabet {
      if !self.initial_prob.contains(state) {
        self.initial_prob.set(state, prob)
      }
    }
  }
}

///|
/// Reset all probabilities to be an average value.
///
/// Resets the values of all initial probabilities and all allowed
/// transitions and all allowed emissions to be equal to 1 divided by the
/// number of possible elements. Note that, exactly like the original
/// implementation, the initial probability denominator is the number of
/// allowed transitions.
pub fn MarkovModelBuilder::set_equal_probabilities(
  self : MarkovModelBuilder,
) -> Unit raise HmmError {
  let n_trans = self.transition_order.length()
  if n_trans == 0 {
    hmm_fail("No transitions have been allowed yet; cannot set equal probabilities")
  }
  // set initial state probabilities (original quirk: divided by the number
  // of transitions, not states)
  let new_initial_prob = 1.0 / n_trans.to_double()
  for state in self.state_alphabet {
    self.initial_prob.set(state, new_initial_prob)
  }
  // set the transitions
  let new_trans_prob = 1.0 / n_trans.to_double()
  for key in self.transition_order {
    hmm_set2(self.transition_prob, key.0, key.1, new_trans_prob)
  }
  // set the emissions
  let new_emission_prob = 1.0 / hmm_map2_size(self.emission_prob).to_double()
  for state in self.state_alphabet {
    for emission in self.emission_alphabet {
      hmm_set2(self.emission_prob, state, emission, new_emission_prob)
    }
  }
}

///|
/// Set all initial state probabilities to a randomly generated
/// distribution. Returns the map containing the initial probabilities.
pub fn MarkovModelBuilder::set_random_initial_probabilities(
  self : MarkovModelBuilder,
) -> Map[String, Double] {
  let initial_freqs = hmm_gen_random_array(self.state_alphabet.length())
  let mut idx = self.state_alphabet.length() - 1
  for state in self.state_alphabet {
    self.initial_prob.set(state, initial_freqs[idx])
    idx = idx - 1
  }
  self.initial_prob
}

///|
/// Set all allowed transition probabilities to a randomly generated
/// distribution. Returns the map containing the transition probabilities.
pub fn MarkovModelBuilder::set_random_transition_probabilities(
  self : MarkovModelBuilder,
) -> Map[String, Map[String, Double]] raise HmmError {
  if self.transition_order.length() == 0 {
    hmm_fail(
      "No transitions have been allowed yet. Allow some or all transitions by calling allow_transition or allow_all_transitions first.",
    )
  }
  // group destination states by source state, keeping first-occurrence
  // order of the source states and insertion order of the destinations
  let groups : Map[String, Array[String]] = Map([], capacity=8)
  let from_order : Array[String] = []
  for key in self.transition_order {
    match groups.get(key.0) {
      Some(dests) => dests.push(key.1)
      None => {
        let dests : Array[String] = [key.1]
        groups.set(key.0, dests)
        from_order.push(key.0)
      }
    }
  }
  for from_state in from_order {
    let dests = groups[from_state]
    let freqs = hmm_gen_random_array(dests.length())
    let mut idx = dests.length() - 1
    for to_state in dests {
      hmm_set2(self.transition_prob, from_state, to_state, freqs[idx])
      idx = idx - 1
    }
  }
  self.transition_prob
}

///|
/// Set all allowed emission probabilities to a randomly generated
/// distribution. Returns the map containing the emission probabilities.
pub fn MarkovModelBuilder::set_random_emission_probabilities(
  self : MarkovModelBuilder,
) -> Map[String, Map[String, Double]] raise HmmError {
  if hmm_map2_size(self.emission_prob) == 0 {
    hmm_fail("No emissions have been allowed yet. Allow some or all emissions.")
  }
  for state in self.state_alphabet {
    let freqs = hmm_gen_random_array(self.emission_alphabet.length())
    let mut idx = self.emission_alphabet.length() - 1
    for symbol in self.emission_alphabet {
      hmm_set2(self.emission_prob, state, symbol, freqs[idx])
      idx = idx - 1
    }
  }
  self.emission_prob
}

///|
/// Set all probabilities to randomly generated numbers.
///
/// Resets probabilities of all initial states, transitions, and emissions
/// to random values.
pub fn MarkovModelBuilder::set_random_probabilities(
  self : MarkovModelBuilder,
) -> Unit raise HmmError {
  let _ = self.set_random_initial_probabilities()
  let _ = self.set_random_transition_probabilities()
  let _ = self.set_random_emission_probabilities()
}

///|
/// Create transitions between all states.
///
/// By default all transitions within the alphabet are disallowed; this is
/// a convenience function to change this to allow all possible transitions.
pub fn MarkovModelBuilder::allow_all_transitions(self : MarkovModelBuilder) -> Unit {
  // first get all probabilities and pseudo counts set to the default values
  let all_probs = hmm_map2()
  let all_pseudo = hmm_map2()
  for s1 in self.state_alphabet {
    for s2 in self.state_alphabet {
      hmm_set2(all_probs, s1, s2, 0.0)
      hmm_set2(all_pseudo, s1, s2, hmm_default_pseudo)
    }
  }
  // now set any probabilities and pseudo counts that were previously set
  for key in self.transition_order {
    let p = match self.transition_prob.get(key.0) {
      Some(inner) => inner.get(key.1).unwrap_or(0.0)
      None => 0.0
    }
    hmm_set2(all_probs, key.0, key.1, p)
    let q = match self.transition_pseudo.get(key.0) {
      Some(inner) => inner.get(key.1).unwrap_or(hmm_default_pseudo)
      None => hmm_default_pseudo
    }
    hmm_set2(all_pseudo, key.0, key.1, q)
  }
  // finally reinitialize the transition probs and pseudo counts
  self.transition_prob = all_probs
  self.transition_pseudo = all_pseudo
  let order : Array[(String, String)] = []
  for s1 in self.state_alphabet {
    for s2 in self.state_alphabet {
      order.push((s1, s2))
    }
  }
  self.transition_order = order
}

///|
/// Set a transition as being possible between the two states.
///
/// probability and pseudocount are optional arguments specifying the
/// probabilities and pseudo counts for the transition; when not supplied
/// they default to 0 and 1 respectively.
pub fn MarkovModelBuilder::allow_transition(
  self : MarkovModelBuilder,
  from_state~ : String,
  to_state~ : String,
  probability? : Double,
  pseudocount? : Double,
) -> Unit raise HmmError {
  // check the sanity of adding these states
  for state in [from_state, to_state] {
    if !(self.state_alphabet.contains(state)) {
      hmm_fail("State \{state} was not found in the sequence alphabet")
    }
  }
  // ensure that the states are not already set
  if !(self.transition_order.contains((from_state, to_state))) {
    // set the initial probability
    let prob = probability.unwrap_or(0.0)
    hmm_set2(self.transition_prob, from_state, to_state, prob)
    // set the initial pseudocounts
    let pseudo = pseudocount.unwrap_or(hmm_default_pseudo)
    hmm_set2(self.transition_pseudo, from_state, to_state, pseudo)
    self.transition_order.push((from_state, to_state))
  } else {
    hmm_fail("Transition from \{from_state} to \{to_state} is already allowed.")
  }
}

///|
/// Restrict transitions between the two states.
pub fn MarkovModelBuilder::destroy_transition(
  self : MarkovModelBuilder,
  from_state~ : String,
  to_state~ : String,
) -> Unit raise HmmError {
  if self.transition_order.contains((from_state, to_state)) {
    hmm_remove2(self.transition_prob, from_state, to_state)
    hmm_remove2(self.transition_pseudo, from_state, to_state)
    let mut idx = -1
    for i = 0; i < self.transition_order.length(); i = i + 1 {
      if self.transition_order[i] == (from_state, to_state) {
        idx = i
        break
      }
    }
    if idx >= 0 {
      let _ = self.transition_order.remove(idx)
    }
  } else {
    hmm_fail("Transition from \{from_state} to \{to_state} is already disallowed.")
  }
}

///|
/// Set the probability of a transition between two states.
pub fn MarkovModelBuilder::set_transition_score(
  self : MarkovModelBuilder,
  from_state~ : String,
  to_state~ : String,
  probability~ : Double,
) -> Unit raise HmmError {
  if self.transition_order.contains((from_state, to_state)) {
    hmm_set2(self.transition_prob, from_state, to_state, probability)
  } else {
    hmm_fail("Transition from \{from_state} to \{to_state} is not allowed.")
  }
}

///|
/// Set the default pseudocount for a transition.
pub fn MarkovModelBuilder::set_transition_pseudocount(
  self : MarkovModelBuilder,
  from_state~ : String,
  to_state~ : String,
  count~ : Double,
) -> Unit raise HmmError {
  if self.transition_order.contains((from_state, to_state)) {
    hmm_set2(self.transition_pseudo, from_state, to_state, count)
  } else {
    hmm_fail("Transition from \{from_state} to \{to_state} is not allowed.")
  }
}

///|
/// Set the probability of an emission from a particular state.
pub fn MarkovModelBuilder::set_emission_score(
  self : MarkovModelBuilder,
  seq_state~ : String,
  emission_state~ : String,
  probability~ : Double,
) -> Unit raise HmmError {
  if hmm_contains2(self.emission_prob, seq_state, emission_state) {
    hmm_set2(self.emission_prob, seq_state, emission_state, probability)
  } else {
    hmm_fail("Emission of \{emission_state} from \{seq_state} is not allowed.")
  }
}

///|
/// Set the default pseudocount for an emission.
pub fn MarkovModelBuilder::set_emission_pseudocount(
  self : MarkovModelBuilder,
  seq_state~ : String,
  emission_state~ : String,
  count~ : Double,
) -> Unit raise HmmError {
  if hmm_contains2(self.emission_pseudo, seq_state, emission_state) {
    hmm_set2(self.emission_pseudo, seq_state, emission_state, count)
  } else {
    hmm_fail("Emission of \{emission_state} from \{seq_state} is not allowed.")
  }
}

// --- inspection helpers (not in the original Python API) ---

///|
/// Number of currently allowed transitions.
pub fn MarkovModelBuilder::transition_count(self : MarkovModelBuilder) -> Int {
  self.transition_order.length()
}

///|
/// Number of allowed (state, emission) pairs.
pub fn MarkovModelBuilder::emission_pair_count(
  self : MarkovModelBuilder,
) -> Int {
  hmm_map2_size(self.emission_prob)
}

///|
/// Get the current transition probability (0 when not allowed).
pub fn MarkovModelBuilder::transition_score(
  self : MarkovModelBuilder,
  from_state : String,
  to_state : String,
) -> Double {
  match self.transition_prob.get(from_state) {
    Some(inner) => inner.get(to_state).unwrap_or(0.0)
    None => 0.0
  }
}

///|
/// Get the current transition pseudo count (0 when not allowed).
pub fn MarkovModelBuilder::transition_pseudocount(
  self : MarkovModelBuilder,
  from_state : String,
  to_state : String,
) -> Double {
  match self.transition_pseudo.get(from_state) {
    Some(inner) => inner.get(to_state).unwrap_or(0.0)
    None => 0.0
  }
}

///|
/// Get the current emission probability.
pub fn MarkovModelBuilder::emission_score(
  self : MarkovModelBuilder,
  seq_state : String,
  emission_state : String,
) -> Double {
  match self.emission_prob.get(seq_state) {
    Some(inner) => inner.get(emission_state).unwrap_or(0.0)
    None => 0.0
  }
}

///|
/// Get the current emission pseudo count.
pub fn MarkovModelBuilder::emission_pseudocount(
  self : MarkovModelBuilder,
  seq_state : String,
  emission_state : String,
) -> Double {
  match self.emission_pseudo.get(seq_state) {
    Some(inner) => inner.get(emission_state).unwrap_or(0.0)
    None => 0.0
  }
}

///|
/// Get the current initial probability of a state (raises when unset).
pub fn MarkovModelBuilder::initial_probability(
  self : MarkovModelBuilder,
  state : String,
) -> Double raise HmmError {
  match self.initial_prob.get(state) {
    Some(v) => v
    None => hmm_fail("Missing initial probability for state \{state}")
  }
}

// ---------------------------------------------------------------------------
// HiddenMarkovModel
// ---------------------------------------------------------------------------

///|
/// Represent a hidden markov model that can be used for state estimation.
pub struct HiddenMarkovModel {
  priv state_alphabet : Array[String]
  priv emission_alphabet : Array[String]
  priv initial_prob : Map[String, Double]
  priv mut transition_prob : Map[String, Map[String, Double]]
  priv mut emission_prob : Map[String, Map[String, Double]]
  priv transition_pseudo : Map[String, Map[String, Double]]
  priv emission_pseudo : Map[String, Map[String, Double]]
  priv transitions_from : Map[String, Array[String]]
  priv transitions_to : Map[String, Array[String]]
}

///|
fn HiddenMarkovModel::build(
  state_alphabet : Array[String],
  emission_alphabet : Array[String],
  initial_prob : Map[String, Double],
  transition_prob : Map[String, Map[String, Double]],
  emission_prob : Map[String, Map[String, Double]],
  transition_pseudo : Map[String, Map[String, Double]],
  emission_pseudo : Map[String, Map[String, Double]],
  transition_order : Array[(String, String)],
) -> HiddenMarkovModel {
  let transitions_from : Map[String, Array[String]] = Map([], capacity=8)
  let transitions_to : Map[String, Array[String]] = Map([], capacity=8)
  for key in transition_order {
    match transitions_from.get(key.0) {
      Some(dests) => dests.push(key.1)
      None => {
        transitions_from.set(key.0, [key.1])
      }
    }
    match transitions_to.get(key.1) {
      Some(srcs) => srcs.push(key.0)
      None => {
        transitions_to.set(key.1, [key.0])
      }
    }
  }
  {
    state_alphabet,
    emission_alphabet,
    initial_prob,
    transition_prob,
    emission_prob,
    transition_pseudo,
    emission_pseudo,
    transitions_from,
    transitions_to,
  }
}

///|
/// Get the state alphabet of the model.
pub fn HiddenMarkovModel::get_state_alphabet(
  self : HiddenMarkovModel,
) -> Array[String] {
  self.state_alphabet
}

///|
/// Get the emission alphabet of the model.
pub fn HiddenMarkovModel::get_emission_alphabet(
  self : HiddenMarkovModel,
) -> Array[String] {
  self.emission_alphabet
}

///|
/// Get the initial probability of a state.
pub fn HiddenMarkovModel::get_initial_probability(
  self : HiddenMarkovModel,
  state : String,
) -> Double raise HmmError {
  match self.initial_prob.get(state) {
    Some(v) => v
    None => hmm_fail("Missing initial probability for state \{state}")
  }
}

///|
/// Get the transition probability from one state to another.
pub fn HiddenMarkovModel::get_transition_probability(
  self : HiddenMarkovModel,
  from_state : String,
  to_state : String,
) -> Double raise HmmError {
  hmm_get2(self.transition_prob, from_state, to_state)
}

///|
/// Get the probability of a state emitting a symbol.
pub fn HiddenMarkovModel::get_emission_probability(
  self : HiddenMarkovModel,
  state : String,
  symbol : String,
) -> Double raise HmmError {
  hmm_get2(self.emission_prob, state, symbol)
}

///|
/// Get the default (pseudo count) transitions of the model. Note that, as
/// in the original implementation, the returned map is the internal pseudo
/// count table used as the starting point of each training iteration.
pub fn HiddenMarkovModel::get_blank_transitions(
  self : HiddenMarkovModel,
) -> Map[String, Map[String, Double]] {
  self.transition_pseudo
}

///|
/// Get the starting default (pseudo count) emissions for each sequence.
pub fn HiddenMarkovModel::get_blank_emissions(
  self : HiddenMarkovModel,
) -> Map[String, Map[String, Double]] {
  self.emission_pseudo
}

///|
/// Get all destination states which can transition from source state.
///
/// An empty list is returned if the state has no outgoing transitions.
pub fn HiddenMarkovModel::transitions_from(
  self : HiddenMarkovModel,
  state_letter : String,
) -> Array[String] {
  match self.transitions_from.get(state_letter) {
    Some(dests) => dests
    None => []
  }
}

///|
/// Get all source states which can transition to destination state.
///
/// An empty list is returned if the state is unreachable.
pub fn HiddenMarkovModel::transitions_to(
  self : HiddenMarkovModel,
  state_letter : String,
) -> Array[String] {
  match self.transitions_to.get(state_letter) {
    Some(srcs) => srcs
    None => []
  }
}

///|
/// The result of the Viterbi decoding: the most probable state path and
/// its (log-space) probability.
pub struct HmmDecodeResult {
  priv state_path : Array[String]
  priv prob : Double
}

///|
/// Get the decoded state path.
pub fn HmmDecodeResult::state_path(self : HmmDecodeResult) -> Array[String] {
  self.state_path
}

///|
/// Get the log-space probability of the decoded path.
pub fn HmmDecodeResult::prob(self : HmmDecodeResult) -> Double {
  self.prob
}

///|
/// Log transform of a probability map: values <= 0 become -inf.
fn hmm_log_transform_2(
  probability : Map[String, Map[String, Double]],
) -> Map[String, Map[String, Double]] {
  let out = hmm_map2()
  for k1, inner in probability {
    for k2, v in inner {
      let lv = if v > 0.0 { @math.ln(v) } else { @double.neg_infinity }
      hmm_set2(out, k1, k2, lv)
    }
  }
  out
}

///|
fn hmm_log_transform_1(
  probability : Map[String, Double],
) -> Map[String, Double] {
  let out : Map[String, Double] = Map([], capacity=8)
  for k, v in probability {
    let lv = if v > 0.0 { @math.ln(v) } else { @double.neg_infinity }
    out.set(k, lv)
  }
  out
}

///|
/// Maximum likelihood variable container keyed by (state, position).
pub struct HmmDpVars {
  priv values : Map[String, Map[Int, Double]]
}

///|
/// Get the value for (state, position), raising when it was never set
/// (mirrors the Python KeyError).
pub fn HmmDpVars::get(
  self : HmmDpVars,
  state : String,
  pos : Int,
) -> Double raise HmmError {
  match self.values.get(state) {
    Some(inner) =>
      match inner.get(pos) {
        Some(v) => v
        None =>
          hmm_fail("No dynamic programming value for state \{state} at \{pos}")
      }
    None => hmm_fail("No dynamic programming value for state \{state} at \{pos}")
  }
}

///|
fn HmmDpVars::set(
  self : HmmDpVars,
  state : String,
  pos : Int,
  value : Double,
) -> Unit {
  match self.values.get(state) {
    Some(inner) => inner.set(pos, value)
    None => {
      let inner : Map[Int, Double] = Map([], capacity=16)
      inner.set(pos, value)
      self.values.set(state, inner)
    }
  }
}

///|
/// The result of the scaled forward algorithm.
pub struct HmmForwardResult {
  priv forward_vars : HmmDpVars
  priv seq_prob : Double
}

///|
/// Get the forward variables keyed by (state, position).
pub fn HmmForwardResult::forward_vars(self : HmmForwardResult) -> HmmDpVars {
  self.forward_vars
}

///|
/// Get the calculated probability of the sequence.
pub fn HmmForwardResult::seq_prob(self : HmmForwardResult) -> Double {
  self.seq_prob
}

///|
/// Calculate the most probable state path using the Viterbi algorithm.
///
/// This implements the Viterbi algorithm (see pgs 55-57 in Durbin et al)
/// to allow decoding of the state path, given a sequence of emissions.
/// The returned probability is in log space.
pub fn HiddenMarkovModel::viterbi(
  self : HiddenMarkovModel,
  sequence : Array[String],
  state_alphabet? : Array[String],
) -> HmmDecodeResult raise HmmError {
  let letters = match state_alphabet {
    Some(letters) => letters
    None => self.state_alphabet
  }
  // calculate logarithms of the initial, transition, and emission probs
  let log_initial = hmm_log_transform_1(self.initial_prob)
  let log_trans = hmm_log_transform_2(self.transition_prob)
  let log_emission = hmm_log_transform_2(self.emission_prob)
  let viterbi_probs : HmmDpVars = { values: Map([], capacity=16) }
  let pred_state_seq : Map[String, String] = Map([], capacity=16)
  let seq_len = sequence.length()
  // --- recursion
  for i = 0; i < seq_len; i = i + 1 {
    // loop over all of the possible i-th states in the state path
    for cur_state in letters {
      // e_{l}(x_{i})
      let emission_part = match log_emission.get(cur_state) {
        Some(inner) =>
          match inner.get(sequence[i]) {
            Some(v) => v
            None =>
              hmm_fail("Emission of \{sequence[i]} from \{cur_state} is not allowed.")
          }
        None =>
          hmm_fail("Emission of \{sequence[i]} from \{cur_state} is not allowed.")
      }
      let max_prob : Double = if i == 0 {
        // for the first state, use the initial probability rather than
        // looking back to previous states
        match log_initial.get(cur_state) {
          Some(v) => v
          None => hmm_fail("Missing initial probability for state \{cur_state}")
        }
      } else {
        // loop over all possible (i-1)-th previous states
        let mut best = -1.7976931348623157e308
        let mut found = false
        for prev_state in self.transitions_to(cur_state) {
          let trans_part = match log_trans.get(prev_state) {
            Some(inner) =>
              match inner.get(cur_state) {
                Some(v) => v
                None =>
                  hmm_fail(
                    "Transition from \{prev_state} to \{cur_state} is not allowed.",
                  )
              }
            None =>
              hmm_fail("Transition from \{prev_state} to \{cur_state} is not allowed.")
          }
          let viterbi_part = viterbi_probs.get(prev_state, i - 1)
          let cur_prob = viterbi_part + trans_part
          if !found || cur_prob > best {
            best = cur_prob
            found = true
          }
        }
        if !found {
          hmm_fail("No allowed transitions into state \{cur_state}")
        }
        best
      }
      // v_{k}(i)
      viterbi_probs.set(cur_state, i, emission_part + max_prob)
      if i > 0 {
        // get the most likely prev_state leading to cur_state
        for prev_state in self.transitions_to(cur_state) {
          let trans_part = match log_trans.get(prev_state) {
            Some(inner) =>
              match inner.get(cur_state) {
                Some(v) => v
                None =>
                  hmm_fail(
                    "Transition from \{prev_state} to \{cur_state} is not allowed.",
                  )
              }
            None =>
              hmm_fail("Transition from \{prev_state} to \{cur_state} is not allowed.")
          }
          let cur_prob = viterbi_probs.get(prev_state, i - 1) + trans_part
          if cur_prob == max_prob {
            pred_state_seq.set("\{i - 1}:\{cur_state}", prev_state)
            break
          }
        }
      }
    }
  }
  // --- termination
  // calculate the probability of the state path
  let mut state_path_prob = -1.7976931348623157e308
  let mut last_state = ""
  for state in letters {
    let p = viterbi_probs.get(state, seq_len - 1)
    if p > state_path_prob {
      state_path_prob = p
      last_state = state
    } else if p == state_path_prob {
      // like the original dict iteration, later equal states win
      last_state = state
    }
  }
  // --- traceback
  let traceback_seq : Array[String] = Array::new(capacity=seq_len)
  let mut state = last_state
  traceback_seq.push(state)
  for i = seq_len - 1; i >= 1; i = i - 1 {
    match pred_state_seq.get("\{i - 1}:\{state}") {
      Some(prev) => state = prev
      None => hmm_fail("Didn't find the state to trace from at position \{i - 1}")
    }
    traceback_seq.push(state)
  }
  // put the traceback sequence in the proper orientation
  let traceback_seq = traceback_seq.rev()
  { state_path: traceback_seq, prob: state_path_prob }
}

// ---------------------------------------------------------------------------
// Scaled dynamic programming (forward/backward algorithms)
// ---------------------------------------------------------------------------

///|
/// Implement forward and backward algorithms using a rescaling approach.
///
/// This scales the f and b variables so that they remain within a
/// manageable numerical interval during calculations (Durbin et al. p78).
pub struct ScaledDp {
  priv mm : HiddenMarkovModel
  priv seq : TrainingSequence
  priv s_values : Map[Int, Double]
}

///|
/// Initialize the scaled approach to calculating probabilities.
///
/// - markov_model: the current Markov model we are working with.
/// - sequence: a TrainingSequence that must have a set of emissions.
pub fn ScaledDp::new(
  markov_model~ : HiddenMarkovModel,
  sequence~ : TrainingSequence,
) -> ScaledDp {
  {
    mm: markov_model,
    seq: sequence,
    s_values: Map([], capacity=16),
  }
}

///|
/// Calculate the next scaling variable for a sequence position.
fn ScaledDp::calculate_s_value(
  self : ScaledDp,
  seq_pos : Int,
  previous_vars : HmmDpVars,
) -> Double raise HmmError {
  let state_letters = self.mm.state_alphabet
  let mut s_value = 0.0
  for main_state in state_letters {
    let emission = self.mm.get_emission_probability(
      main_state,
      self.seq.emissions[seq_pos],
    )
    // now sum over all of the previous vars and transitions
    let mut trans_and_var_sum = 0.0
    for second_state in self.mm.transitions_from(main_state) {
      // the value of the previous f or b value
      let var_value = previous_vars.get(second_state, seq_pos - 1)
      // the transition probability
      let trans_value = self.mm.get_transition_probability(second_state, main_state)
      trans_and_var_sum = trans_and_var_sum + var_value * trans_value
    }
    s_value = s_value + emission * trans_and_var_sum
  }
  s_value
}

///|
/// Calculate the value of the forward recursion. Returns None when the
/// state has no outgoing transitions.
fn ScaledDp::forward_recursion(
  self : ScaledDp,
  cur_state : String,
  sequence_pos : Int,
  forward_vars : HmmDpVars,
) -> Double? raise HmmError {
  // calculate the s value, if we haven't done so already
  if !self.s_values.contains(sequence_pos) {
    self.s_values.set(sequence_pos, self.calculate_s_value(sequence_pos, forward_vars))
  }
  // e_{l}(x_{i})
  let seq_letter = self.seq.emissions[sequence_pos]
  let cur_emission_prob = self.mm.get_emission_probability(cur_state, seq_letter)
  // divide by the scaling value
  let scale_emission_prob = cur_emission_prob / self.s_values[sequence_pos]
  // loop over all of the possible states at the position
  let mut state_pos_sum = 0.0
  let mut have_transition = false
  for second_state in self.mm.transitions_from(cur_state) {
    have_transition = true
    // get the previous forward_var values f_{k}(i - 1)
    let prev_forward = forward_vars.get(second_state, sequence_pos - 1)
    // a_{kl}
    let cur_trans_prob = self.mm.get_transition_probability(second_state, cur_state)
    state_pos_sum = state_pos_sum + prev_forward * cur_trans_prob
  }
  // if we have the possibility of having a transition return the value
  if have_transition {
    Some(scale_emission_prob * state_pos_sum)
  } else {
    None
  }
}

///|
/// Calculate the value of the backward recursion. Returns None when the
/// state has no outgoing transitions.
fn ScaledDp::backward_recursion(
  self : ScaledDp,
  cur_state : String,
  sequence_pos : Int,
  backward_vars : HmmDpVars,
) -> Double? raise HmmError {
  // calculate the s value, if we haven't done so already
  if !self.s_values.contains(sequence_pos) {
    self.s_values.set(sequence_pos, self.calculate_s_value(sequence_pos, backward_vars))
  }
  let mut state_pos_sum = 0.0
  let mut have_transition = false
  for second_state in self.mm.transitions_from(cur_state) {
    have_transition = true
    // e_{l}(x_{i + 1})
    let seq_letter = self.seq.emissions[sequence_pos + 1]
    let cur_emission_prob = self.mm.get_emission_probability(cur_state, seq_letter)
    // get the previous backward_var value b_{l}(i + 1)
    let prev_backward = backward_vars.get(second_state, sequence_pos + 1)
    // the transition probability a_{kl}
    let cur_transition_prob = self.mm.get_transition_probability(cur_state, second_state)
    state_pos_sum = state_pos_sum + cur_emission_prob * prev_backward * cur_transition_prob
  }
  if have_transition {
    Some(state_pos_sum / self.s_values[sequence_pos])
  } else {
    None
  }
}

///|
/// Calculate sequence probability using the forward algorithm.
///
/// This implements the forward algorithm, as described on p57-58 of
/// Durbin et al. Returns the forward variables and the probability of
/// the sequence.
pub fn ScaledDp::forward_algorithm(self : ScaledDp) -> HmmForwardResult raise HmmError {
  // all of the different letters that the state path can be in
  let state_letters = self.mm.state_alphabet
  // -- initialize the algorithm: f_{0}(0) = 1, f_{k}(0) = 0 for k > 0
  let forward_vars : HmmDpVars = { values: Map([], capacity=16) }
  forward_vars.set(state_letters[0], -1, 1.0)
  for k = 1; k < state_letters.length(); k = k + 1 {
    forward_vars.set(state_letters[k], -1, 0.0)
  }
  // -- now do the recursion step
  let emissions_len = self.seq.emissions.length()
  for i = 0; i < emissions_len; i = i + 1 {
    for main_state in state_letters {
      let forward_value = self.forward_recursion(main_state, i, forward_vars)
      match forward_value {
        Some(v) => forward_vars.set(main_state, i, v)
        None => ()
      }
    }
  }
  // -- termination step - calculate the probability of the sequence
  let first_state = state_letters[0]
  let mut seq_prob = 0.0
  for state_item in state_letters {
    // f_{k}(L)
    let forward_value = forward_vars.get(state_item, emissions_len - 1)
    // a_{k0}
    let transition_value = self.mm.get_transition_probability(state_item, first_state)
    seq_prob = seq_prob + forward_value * transition_value
  }
  { forward_vars, seq_prob }
}

///|
/// Calculate sequence probability using the backward algorithm.
///
/// This implements the backward algorithm, as described on p58-59 of
/// Durbin et al. The termination step is skipped to avoid recalculations;
/// sequence probabilities should be obtained from the forward algorithm.
pub fn ScaledDp::backward_algorithm(self : ScaledDp) -> HmmDpVars raise HmmError {
  let state_letters = self.mm.state_alphabet
  // -- initialize: b_{k}(L) = a_{k0} for all k
  let backward_vars : HmmDpVars = { values: Map([], capacity=16) }
  let last_pos = self.seq.emissions.length() - 1
  for state in state_letters {
    backward_vars.set(state, last_pos, self.mm.get_transition_probability(state, state_letters[0]))
  }
  // -- recursion: (i = L - 1 ... 1) in 0-based positions
  let emissions_len = self.seq.emissions.length()
  for i = emissions_len - 2; i >= 0; i = i - 1 {
    for main_state in state_letters {
      let backward_value = self.backward_recursion(main_state, i, backward_vars)
      match backward_value {
        Some(v) => backward_vars.set(main_state, i, v)
        None => ()
      }
    }
  }
  backward_vars
}

// ---------------------------------------------------------------------------
// Trainers
// ---------------------------------------------------------------------------

///|
/// Hold a training sequence with emissions and optionally, a state path.
pub struct TrainingSequence {
  priv emissions : Array[String]
  priv states : Array[String]
}

///|
/// Initialize a training sequence.
///
/// - emissions: the sequence of emissions in the training sequence.
/// - state_path: the sequence of states. If there is no known state path,
///   then the sequence of states should be empty.
pub fn TrainingSequence::new(
  emissions~ : Array[String],
  state_path~ : Array[String],
) -> TrainingSequence raise HmmError {
  if state_path.length() > 0 && emissions.length() != state_path.length() {
    hmm_fail("State path does not match associated emissions.")
  }
  let emissions_copy : Array[String] = Array::new(capacity=emissions.length())
  for e in emissions {
    emissions_copy.push(e)
  }
  let states_copy : Array[String] = Array::new(capacity=state_path.length())
  for s in state_path {
    states_copy.push(s)
  }
  { emissions: emissions_copy, states: states_copy }
}

///|
/// Get the emission sequence of the training sequence.
pub fn TrainingSequence::emissions(self : TrainingSequence) -> Array[String] {
  self.emissions
}

///|
/// Get the state path of the training sequence (empty when unknown).
pub fn TrainingSequence::states(self : TrainingSequence) -> Array[String] {
  self.states
}

///|
/// Calculate the log likelihood of the training sequences.
///
/// - probabilities: the probabilities of each training sequence under the
///   current parameters, calculated using the forward algorithm.
fn hmm_log_likelihood(probabilities : Array[Double]) -> Double {
  let mut total_likelihood = 0.0
  for probability in probabilities {
    total_likelihood = total_likelihood + @math.ln(probability)
  }
  total_likelihood
}

///|
/// Safe two-level map getter for keys that are guaranteed to exist
/// (used when iterating the map's own keys).
fn hmm_get2_unchecked(
  m : Map[String, Map[String, Double]],
  k1 : String,
  k2 : String,
) -> Double {
  match m.get(k1) {
    Some(inner) => inner.get(k2).unwrap_or(0.0)
    None => 0.0
  }
}

///|
/// Calculate the maximum likelihood estimator for transitions or emissions.
///
/// The counts are grouped by the first letter of the key and normalized
/// with formula 3.18 of Durbin et al:
///   a_{kl} = A_{kl} / sum(A_{kl'}), e_{k}(b) = E_{k}(b) / sum(E_{k}(b'))
fn hmm_ml_estimator(
  counts : Map[String, Map[String, Double]],
) -> Map[String, Map[String, Double]] {
  // get an ordered list of all items (lexicographic like Python's sorted)
  let all_ordered : Array[(String, String)] = []
  for k1, inner in counts {
    for k2, _v in inner {
      all_ordered.push((k1, k2))
    }
  }
  all_ordered.sort_by(fn(a : (String, String), b : (String, String)) -> Int {
    let c = a.0.compare(b.0)
    if c != 0 {
      c
    } else {
      a.1.compare(b.1)
    }
  })
  let ml_estimation = hmm_map2()
  let mut cur_letter : String? = None
  let mut cur_letter_counts = 0.0
  for pos = 0; pos < all_ordered.length(); pos = pos + 1 {
    let cur_item = all_ordered[pos]
    // if we are on a new letter (ie. the first letter of the tuple)
    let new_letter = match cur_letter {
      Some(letter) => letter != cur_item.0
      None => true
    }
    if new_letter {
      cur_letter = Some(cur_item.0)
      // count up the total counts for this letter
      cur_letter_counts = hmm_get2_unchecked(counts, cur_item.0, cur_item.1)
      // add counts for all other items with the same first letter
      let mut cur_position = pos + 1
      while cur_position < all_ordered.length() &&
        all_ordered[cur_position].0 == cur_item.0 {
        cur_letter_counts = cur_letter_counts +
          hmm_get2_unchecked(
            counts,
            all_ordered[cur_position].0,
            all_ordered[cur_position].1,
          )
        cur_position = cur_position + 1
      }
    }
    // now calculate the ml and add it to the estimation
    let cur_ml = hmm_get2_unchecked(counts, cur_item.0, cur_item.1) / cur_letter_counts
    hmm_set2(ml_estimation, cur_item.0, cur_item.1, cur_ml)
  }
  ml_estimation
}

///|
/// Get a maximum likelihood estimation of transition and emission.
fn hmm_estimate_params(
  transition_counts : Map[String, Map[String, Double]],
  emission_counts : Map[String, Map[String, Double]],
) -> (Map[String, Map[String, Double]], Map[String, Map[String, Double]]) {
  let ml_transitions = hmm_ml_estimator(transition_counts)
  let ml_emissions = hmm_ml_estimator(emission_counts)
  (ml_transitions, ml_emissions)
}

///|
/// Trainer that uses the Baum-Welch algorithm to estimate parameters.
///
/// These should be used when a training sequence for an HMM has unknown
/// paths for the actual states, and you need to make an estimation of the
/// model parameters from the observed emissions.
///
/// This algorithm is guaranteed to converge to a local maximum, but not
/// necessarily to the global maxima, so use with care!
pub struct BaumWelchTrainer {
  priv markov_model : HiddenMarkovModel
}

///|
/// Initialize the trainer with a model that has some initial estimates
/// we can build from.
pub fn BaumWelchTrainer::new(markov_model~ : HiddenMarkovModel) -> BaumWelchTrainer {
  { markov_model }
}

///|
/// Estimate the parameters using training sequences.
///
/// The algorithm for this is taken from Durbin et al. p64.
///
/// - training_seqs: training sequences used for estimating the parameters.
/// - stopping_criteria: a function that when passed the change in log
///   likelihood and the iteration number will indicate if we should stop
///   the estimation iterations.
pub fn BaumWelchTrainer::train(
  self : BaumWelchTrainer,
  training_seqs : Array[TrainingSequence],
  stopping_criteria : (Double, Int) -> Bool,
) -> HiddenMarkovModel raise HmmError {
  let mut prev_log_likelihood : Double? = None
  let mut num_iterations = 1
  while true {
    let transition_count = self.markov_model.get_blank_transitions()
    let emission_count = self.markov_model.get_blank_emissions()
    // remember all of the sequence probabilities
    let all_probabilities : Array[Double] = []
    for training_seq in training_seqs {
      // calculate the forward and backward variables
      let dp = ScaledDp::new(markov_model=self.markov_model, sequence=training_seq)
      let forward_result = dp.forward_algorithm()
      let backward_var = dp.backward_algorithm()
      all_probabilities.push(forward_result.seq_prob())
      // update the counts for transitions and emissions
      let _ = self.update_transitions(
        transition_count,
        training_seq,
        forward_result.forward_vars(),
        backward_var,
        forward_result.seq_prob(),
      )
      let _ = self.update_emissions(
        emission_count,
        training_seq,
        forward_result.forward_vars(),
        backward_var,
        forward_result.seq_prob(),
      )
    }
    // update the markov model with the new probabilities
    let (ml_transitions, ml_emissions) = hmm_estimate_params(transition_count, emission_count)
    self.markov_model.transition_prob = ml_transitions
    self.markov_model.emission_prob = ml_emissions
    let cur_log_likelihood = hmm_log_likelihood(all_probabilities)
    // if we have previously calculated the log likelihood (ie. not the
    // first round), see if we can finish
    match prev_log_likelihood {
      Some(prev) => {
        let log_likelihood_change = (cur_log_likelihood.abs() - prev.abs()).abs()
        if stopping_criteria(log_likelihood_change, num_iterations) {
          break
        }
      }
      None => ()
    }
    // set up for another round of iterations
    prev_log_likelihood = Some(cur_log_likelihood)
    num_iterations = num_iterations + 1
  }
  self.markov_model
}

///|
/// Add the contribution of a new training sequence to the transitions.
///
/// This calculates A_{kl} (the estimated transition counts from state k
/// to state l) using formula 3.20 in Durbin et al.
pub fn BaumWelchTrainer::update_transitions(
  self : BaumWelchTrainer,
  transition_counts : Map[String, Map[String, Double]],
  training_seq : TrainingSequence,
  forward_vars : HmmDpVars,
  backward_vars : HmmDpVars,
  training_seq_prob : Double,
) -> Map[String, Map[String, Double]] raise HmmError {
  // set up the transition and emission probabilities we are using
  let transitions = self.markov_model.transition_prob
  let emissions = self.markov_model.emission_prob
  // loop over the possible combinations of state path letters
  for k in self.markov_model.state_alphabet {
    for l in self.markov_model.transitions_from(k) {
      let mut estimated_counts = 0.0
      // now loop over the entire training sequence
      for i = 0; i < training_seq.emissions.length() - 1; i = i + 1 {
        // the forward value of k at the current position
        let forward_value = forward_vars.get(k, i)
        // the backward value of l in the next position
        let backward_value = backward_vars.get(l, i + 1)
        // the probability of a transition from k to l
        let trans_value = hmm_get2(transitions, k, l)
        // the probability of getting the emission at the next pos
        let emm_value = hmm_get2(emissions, l, training_seq.emissions[i + 1])
        estimated_counts = estimated_counts +
          forward_value * trans_value * emm_value * backward_value
      }
      // update the transition approximation
      let prev = hmm_get2(transition_counts, k, l)
      hmm_set2(transition_counts, k, l, prev + estimated_counts / training_seq_prob)
    }
  }
  transition_counts
}

///|
/// Add the contribution of a new training sequence to the emissions.
///
/// This calculates E_{k}(b) (the estimated emission counts of letter b
/// from state k) using formula 3.21 in Durbin et al.
pub fn BaumWelchTrainer::update_emissions(
  self : BaumWelchTrainer,
  emission_counts : Map[String, Map[String, Double]],
  training_seq : TrainingSequence,
  forward_vars : HmmDpVars,
  backward_vars : HmmDpVars,
  training_seq_prob : Double,
) -> Map[String, Map[String, Double]] raise HmmError {
  // loop over the possible combinations of state path letters
  for k in self.markov_model.state_alphabet {
    // now loop over all of the possible emissions
    for b in self.markov_model.emission_alphabet {
      let mut expected_times = 0.0
      // finally loop over the entire training sequence
      for i = 0; i < training_seq.emissions.length(); i = i + 1 {
        // only count the forward and backward probability if the
        // emission at the position is the same as b
        if training_seq.emissions[i] == b {
          // f_{k}(i) b_{k}(i)
          expected_times = expected_times +
            forward_vars.get(k, i) * backward_vars.get(k, i)
        }
      }
      // add to E_{k}(b)
      let prev = hmm_get2(emission_counts, k, b)
      hmm_set2(emission_counts, k, b, prev + expected_times / training_seq_prob)
    }
  }
  emission_counts
}

///|
/// Estimate probabilities with known state sequences.
///
/// This should be used for direct estimation of emission and transition
/// probabilities when both the state path and emission sequence are
/// known for the training examples.
pub struct KnownStateTrainer {
  priv markov_model : HiddenMarkovModel
}

///|
/// Initialize the class.
pub fn KnownStateTrainer::new(markov_model~ : HiddenMarkovModel) -> KnownStateTrainer {
  { markov_model }
}

///|
/// Estimate the Markov Model parameters with known state paths.
///
/// This trainer requires that both the state and the emissions are known
/// for all of the training sequences. Training counts all of the
/// transitions and emissions, and uses this to estimate the parameters
/// of the model.
pub fn KnownStateTrainer::train(
  self : KnownStateTrainer,
  training_seqs : Array[TrainingSequence],
) -> HiddenMarkovModel raise HmmError {
  // count up all of the transitions and emissions
  let transition_counts = self.markov_model.get_blank_transitions()
  let emission_counts = self.markov_model.get_blank_emissions()
  for training_seq in training_seqs {
    let _ = self.count_emissions(training_seq, emission_counts)
    let _ = self.count_transitions(training_seq.states, transition_counts)
  }
  // update the markov model from the counts
  let (ml_transitions, ml_emissions) = hmm_estimate_params(transition_counts, emission_counts)
  self.markov_model.transition_prob = ml_transitions
  self.markov_model.emission_prob = ml_emissions
  self.markov_model
}

///|
/// Add emissions from the training sequence to the current counts.
pub fn KnownStateTrainer::count_emissions(
  self : KnownStateTrainer,
  training_seq : TrainingSequence,
  emission_counts : Map[String, Map[String, Double]],
) -> Map[String, Map[String, Double]] raise HmmError {
  for index = 0; index < training_seq.emissions.length(); index = index + 1 {
    let cur_state = training_seq.states[index]
    let cur_emission = training_seq.emissions[index]
    if hmm_contains2(emission_counts, cur_state, cur_emission) {
      let prev = hmm_get2(emission_counts, cur_state, cur_emission)
      hmm_set2(emission_counts, cur_state, cur_emission, prev + 1.0)
    } else {
      hmm_fail("Unexpected emission (\{cur_state}, \{cur_emission})")
    }
  }
  emission_counts
}

///|
/// Add transitions from the training sequence to the current counts.
pub fn KnownStateTrainer::count_transitions(
  self : KnownStateTrainer,
  state_seq : Array[String],
  transition_counts : Map[String, Map[String, Double]],
) -> Map[String, Map[String, Double]] raise HmmError {
  for cur_pos = 0; cur_pos < state_seq.length() - 1; cur_pos = cur_pos + 1 {
    let cur_state = state_seq[cur_pos]
    let next_state = state_seq[cur_pos + 1]
    if hmm_contains2(transition_counts, cur_state, next_state) {
      let prev = hmm_get2(transition_counts, cur_state, next_state)
      hmm_set2(transition_counts, cur_state, next_state, prev + 1.0)
    } else {
      hmm_fail("Unexpected transition (\{cur_state}, \{next_state})")
    }
  }
  transition_counts
}
