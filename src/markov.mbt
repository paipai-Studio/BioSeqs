///|
/// Bio.Markov - Higher-order Markov chain modeling for biological sequences.
/// Inspired by Biopython Bio.Markov / MarkovModel, used for sequence probability
/// scoring, generation, CpG island detection etc.

///|
/// Markov chain order type.
pub enum ChainType {
  FirstOrder
  SecondOrder
  ThirdOrder
} derive(Eq, Debug)

///|
/// Create a FirstOrder chain type.
pub fn markov_first_order() -> ChainType {
  FirstOrder
}

///|
/// Create a SecondOrder chain type.
pub fn markov_second_order() -> ChainType {
  SecondOrder
}

///|
/// Create a ThirdOrder chain type.
pub fn markov_third_order() -> ChainType {
  ThirdOrder
}

///|
/// Higher-order Markov chain model.
pub struct MarkovModel {
  order : Int
  chain_type : ChainType
  states : Array[String]
  transition_counts : Map[String, Map[String, Double]]
  transition_probs : Map[String, Map[String, Double]]
  initial_probs : Map[String, Double]
  pseudo_count : Double
}

fn markov_empty_map_kv() -> Map[String, Map[String, Double]] {
  Map([], capacity=16)
}

fn markov_empty_map_d() -> Map[String, Double] {
  Map([], capacity=16)
}

fn markov_empty_map_u() -> Map[String, Unit] {
  Map([], capacity=16)
}

///|
/// Create a default MarkovModel: order=1, states=ACGT, pseudo=1.0.
pub fn MarkovModel::new() -> MarkovModel {
  MarkovModel::{
    order: 1,
    chain_type: FirstOrder,
    states: ["A", "C", "G", "T"],
    transition_counts: markov_empty_map_kv(),
    transition_probs: markov_empty_map_kv(),
    initial_probs: markov_empty_map_d(),
    pseudo_count: 1.0,
  }
}

///|
/// Create a default DNA 4-state MarkovModel.
pub fn MarkovModel::default_dna() -> MarkovModel {
  MarkovModel::{
    order: 1,
    chain_type: FirstOrder,
    states: ["A", "C", "G", "T"],
    transition_counts: markov_empty_map_kv(),
    transition_probs: markov_empty_map_kv(),
    initial_probs: markov_empty_map_d(),
    pseudo_count: 1.0,
  }
}

///|
/// Create a default protein 20-state MarkovModel.
pub fn MarkovModel::default_protein() -> MarkovModel {
  MarkovModel::{
    order: 1,
    chain_type: FirstOrder,
    states: [
      "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
      "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
    ],
    transition_counts: markov_empty_map_kv(),
    transition_probs: markov_empty_map_kv(),
    initial_probs: markov_empty_map_d(),
    pseudo_count: 1.0,
  }
}

///|
/// Set order, returning new model.
pub fn MarkovModel::set_order(self : MarkovModel, val : Int) -> MarkovModel {
  let ct = if val == 1 {
    FirstOrder
  } else if val == 2 {
    SecondOrder
  } else if val == 3 {
    ThirdOrder
  } else {
    FirstOrder
  }
  MarkovModel::{
    order: val,
    chain_type: ct,
    states: self.states.copy(),
    transition_counts: self.transition_counts.copy(),
    transition_probs: self.transition_probs.copy(),
    initial_probs: self.initial_probs.copy(),
    pseudo_count: self.pseudo_count,
  }
}

///|
/// Set states, returning new model.
pub fn MarkovModel::set_states(self : MarkovModel, states~ : Array[String]) -> MarkovModel {
  MarkovModel::{
    order: self.order,
    chain_type: self.chain_type,
    states: states.copy(),
    transition_counts: self.transition_counts.copy(),
    transition_probs: self.transition_probs.copy(),
    initial_probs: self.initial_probs.copy(),
    pseudo_count: self.pseudo_count,
  }
}

///|
/// Set chain_type, returning new model.
pub fn MarkovModel::set_chain_type(self : MarkovModel, val~ : ChainType) -> MarkovModel {
  let ord = match val {
    FirstOrder => 1
    SecondOrder => 2
    ThirdOrder => 3
  }
  MarkovModel::{
    order: ord,
    chain_type: val,
    states: self.states.copy(),
    transition_counts: self.transition_counts.copy(),
    transition_probs: self.transition_probs.copy(),
    initial_probs: self.initial_probs.copy(),
    pseudo_count: self.pseudo_count,
  }
}

///|
/// Set pseudo_count, returning new model.
pub fn MarkovModel::set_pseudo_count(self : MarkovModel, val : Double) -> MarkovModel {
  MarkovModel::{
    order: self.order,
    chain_type: self.chain_type,
    states: self.states.copy(),
    transition_counts: self.transition_counts.copy(),
    transition_probs: self.transition_probs.copy(),
    initial_probs: self.initial_probs.copy(),
    pseudo_count: val,
  }
}

///|
/// Infer states from sequences if not provided.
fn markov_infer_states(sequences : Array[String]) -> Array[String] {
  let set = markov_empty_map_u()
  let mut si = 0
  while si < sequences.length() {
    let s = sequences[si]
    let mut ci = 0
    while ci < s.length() {
      let ch = substring(s, ci, 1)
      if !set.contains(ch) {
        set.set(ch, ())
      }
      ci = ci + 1
    }
    si = si + 1
  }
  set.keys().collect()
}

///|
/// Build a Markov chain model from training sequences.
pub fn markov_build_model(
  sequences : Array[String],
  order? : Int = 1,
  states? : Array[String] = [],
  pseudo_count? : Double = 1.0,
) -> MarkovModel {
  let actual_order = if order < 1 { 1 } else if order > 3 { 3 } else { order }
  let actual_states = if states.length() > 0 {
    states.copy()
  } else {
    markov_infer_states(sequences)
  }
  let ct = if actual_order == 1 {
    FirstOrder
  } else if actual_order == 2 {
    SecondOrder
  } else {
    ThirdOrder
  }
  let pseudo = pseudo_count

  let transition_counts = markov_empty_map_kv()
  let initial_counts = markov_empty_map_d()

  let mut si = 0
  while si < sequences.length() {
    let seq = sequences[si]
    let n = seq.length()
    if n > actual_order {
      let start_kmer = substring(seq, 0, actual_order)
      if initial_counts.contains(start_kmer) {
        initial_counts.set(start_kmer, initial_counts.get(start_kmer).unwrap() + 1.0)
      } else {
        initial_counts.set(start_kmer, 1.0)
      }
    }
    let mut i = actual_order
    while i < n {
      let context = substring(seq, i - actual_order, actual_order)
      let next_ch = substring(seq, i, 1)
      if !transition_counts.contains(context) {
        transition_counts.set(context, markov_empty_map_d())
      }
      let row = transition_counts.get(context).unwrap()
      if row.contains(next_ch) {
        row.set(next_ch, row.get(next_ch).unwrap() + 1.0)
      } else {
        row.set(next_ch, 1.0)
      }
      i = i + 1
    }
    si = si + 1
  }

  let total_initial : Double = {
    let mut s = 0.0
    let keys = initial_counts.keys().collect()
    let mut ki = 0
    while ki < keys.length() {
      s = s + initial_counts.get(keys[ki]).unwrap_or(0.0)
      ki = ki + 1
    }
    s
  }

  let initial_probs = markov_empty_map_d()
  let uniform_init = 1.0 / actual_states.length().to_double()
  if total_initial > 0.0 {
    let keys = initial_counts.keys().collect()
    let mut ki = 0
    while ki < keys.length() {
      let k = keys[ki]
      let v = initial_counts.get(k).unwrap_or(0.0)
      initial_probs.set(k, v / total_initial)
      ki = ki + 1
    }
  } else {
    let mut si2 = 0
    while si2 < actual_states.length() {
      initial_probs.set(actual_states[si2], uniform_init)
      si2 = si2 + 1
    }
  }

  let transition_probs = markov_empty_map_kv()
  let uniform_prob = 1.0 / actual_states.length().to_double()
  let context_keys = transition_counts.keys().collect()
  let mut cki = 0
  while cki < context_keys.length() {
    let ctx = context_keys[cki]
    let row_counts = transition_counts.get(ctx).unwrap()
    let mut row_total = 0.0
    let rk = row_counts.keys().collect()
    let mut rki = 0
    while rki < rk.length() {
      row_total = row_total + row_counts.get(rk[rki]).unwrap_or(0.0)
      rki = rki + 1
    }
    let prob_row = markov_empty_map_d()
    if row_total > 0.0 {
      let mut si3 = 0
      while si3 < actual_states.length() {
        let st = actual_states[si3]
        let cnt = row_counts.get(st).unwrap_or(0.0) + pseudo
        let denom = row_total + pseudo * actual_states.length().to_double()
        prob_row.set(st, cnt / denom)
        si3 = si3 + 1
      }
    } else {
      let mut si3 = 0
      while si3 < actual_states.length() {
        prob_row.set(actual_states[si3], uniform_prob)
        si3 = si3 + 1
      }
    }
    transition_probs.set(ctx, prob_row)
    cki = cki + 1
  }

  MarkovModel::{
    order: actual_order,
    chain_type: ct,
    states: actual_states,
    transition_counts,
    transition_probs,
    initial_probs,
    pseudo_count: pseudo,
  }
}

///|
/// Get transition probability with pseudo-count smoothing for unknowns.
fn markov_get_trans_prob(
  model : MarkovModel,
  context : String,
  next : String,
) -> Double {
  let n = model.states.length().to_double()
  match model.transition_probs.get(context) {
    Some(row) =>
      match row.get(next) {
        Some(v) => v
        None => {
          // For a known context but unknown next state:
          // Smooth with pseudo_count. Compute total for this row by summing entries,
          // then P = pseudo / (total + pseudo * N). Higher pseudo -> closer to uniform.
          let mut row_total = 0.0
          let values = row.values().collect()
          let mut vi = 0
          while vi < values.length() {
            row_total = row_total + values[vi]
            vi = vi + 1
          }
          let denom = row_total + model.pseudo_count * n
          if denom > 0.0 { model.pseudo_count / denom } else { 1.0 / n }
        }
      }
    None => {
      // Completely unseen context. Use pseudo-weighted uniform:
      // P = pseudo / (1 + pseudo * N). Higher pseudo -> higher prob, closer to 1/N.
      let denom = 1.0 + model.pseudo_count * n
      if denom > 0.0 { model.pseudo_count / denom } else { 1.0 / n }
    }
  }
}

///|
/// Get initial probability with smoothing for unknowns.
fn markov_get_init_prob(model : MarkovModel, kmer : String) -> Double {
  match model.initial_probs.get(kmer) {
    Some(v) => v
    None => {
      let n = model.initial_probs.length()
      if n > 0 {
        let denom = model.pseudo_count * (n + 1).to_double()
        if denom > 0.0 { model.pseudo_count / denom } else { 1.0 / n.to_double() }
      } else {
        1.0 / model.states.length().to_double()
      }
    }
  }
}

///|
/// Compute log probability of a sequence under the model.
pub fn markov_score_sequence(model : MarkovModel, sequence : String) -> Double {
  let n = sequence.length()
  if n <= model.order {
    return 0.0
  }
  let start_kmer = substring(sequence, 0, model.order)
  let init_p = markov_get_init_prob(model, start_kmer)
  let tiny = model.pseudo_count * 0.0001 + 0.000001
  let mut log_sum = if init_p > 0.0 { @math.ln(init_p) } else { @math.ln(tiny) }
  let mut i = model.order
  while i < n {
    let context = substring(sequence, i - model.order, model.order)
    let next_ch = substring(sequence, i, 1)
    let p = markov_get_trans_prob(model, context, next_ch)
    let lp = if p > 0.0 { @math.ln(p) } else { @math.ln(tiny) }
    log_sum = log_sum + lp
    i = i + 1
  }
  log_sum
}

///|
/// Compute per-position log probabilities.
pub fn markov_score_per_base(model : MarkovModel, sequence : String) -> Array[Double] {
  let n = sequence.length()
  let result : Array[Double] = Array::new()
  if n <= model.order {
    return result
  }
  let tiny = model.pseudo_count * 0.0001 + 0.000001
  let mut i = model.order
  while i < n {
    let context = substring(sequence, i - model.order, model.order)
    let next_ch = substring(sequence, i, 1)
    let p = markov_get_trans_prob(model, context, next_ch)
    let lp = if p > 0.0 { @math.ln(p) } else { @math.ln(tiny) }
    result.push(lp)
    i = i + 1
  }
  result
}

///|
/// Weighted sample from a probability map.
fn markov_weighted_sample(
  probs : Map[String, Double],
  states : Array[String],
  r_val : Double,
) -> String {
  let mut cumul = 0.0
  let mut si = 0
  while si < states.length() {
    let st = states[si]
    let p = probs.get(st).unwrap_or(0.0)
    cumul = cumul + p
    if r_val <= cumul && cumul > 0.0 {
      return st
    }
    si = si + 1
  }
  if states.length() > 0 { states[states.length() - 1] } else { "" }
}

///|
/// Generate a random sequence of given length from the model.
pub fn markov_generate_sequence(
  model : MarkovModel,
  length : Int,
  seed? : Int = 0,
) -> String {
  if length <= 0 {
    return ""
  }

  let mask = 0x7fffffff
  let modulus = 2147483647.0
  let mut rng_state : Int = seed + length

  let n_states = model.states.length()
  let init_keys = model.initial_probs.keys().collect()
  let mut context = ""
  if length >= model.order && init_keys.length() > 0 {
    let total_init : Double = {
      let mut s = 0.0
      let mut ki = 0
      while ki < init_keys.length() {
        s = s + model.initial_probs.get(init_keys[ki]).unwrap_or(0.0)
        ki = ki + 1
      }
      s
    }
    rng_state = (rng_state * 1103515245 + 12345) & mask
    let r_val = rng_state.to_double() / modulus
    let mut cumul = 0.0
    let mut ki = 0
    let mut found = false
    while ki < init_keys.length() && !found {
      let k = init_keys[ki]
      let p = model.initial_probs.get(k).unwrap_or(0.0) / total_init
      cumul = cumul + p
      if r_val <= cumul && cumul > 0.0 {
        context = k
        found = true
      }
      ki = ki + 1
    }
    if !found && init_keys.length() > 0 {
      context = init_keys[0]
    }
  }
  if context.length() == 0 && n_states > 0 {
    let mut ci = 0
    while ci < model.order {
      context = context + model.states[0]
      ci = ci + 1
    }
  }

  let take = if model.order > length { length } else { model.order }
  let mut result = if context.length() >= take { substring(context, 0, take) } else { context }

  let mut step = result.length()
  while step < length {
    let ctx = if result.length() >= model.order {
      substring(result, result.length() - model.order, model.order)
    } else {
      result
    }
    rng_state = (rng_state * 1103515245 + 12345) & mask
    let r_next = rng_state.to_double() / modulus
    let next_char = match model.transition_probs.get(ctx) {
      Some(row) => markov_weighted_sample(row, model.states, r_next)
      None => {
        if n_states > 0 {
          let idx_f = r_next * n_states.to_double()
          let idx = idx_f.to_int()
          let abs_idx = if idx < 0 { -idx } else { idx }
          let bounded = abs_idx % n_states
          model.states[bounded]
        } else {
          ""
        }
      }
    }
    result = result + next_char
    step = step + 1
  }

  if result.length() > length {
    substring(result, 0, length)
  } else {
    result
  }
}

///|
/// Compute log odds ratio between CpG model and background model.
pub fn markov_log_odds(
  model_cpg : MarkovModel,
  model_bg : MarkovModel,
  sequence : String,
) -> Double {
  let s_cpg = markov_score_sequence(model_cpg, sequence)
  let s_bg = markov_score_sequence(model_bg, sequence)
  s_cpg - s_bg
}

///|
/// Compute stationary distribution via iterative multiplication.
pub fn markov_stationary_distribution(model : MarkovModel) -> Map[String, Double] {
  let states = model.states
  let n = states.length()
  let dist = markov_empty_map_d()
  if n == 0 {
    return dist
  }
  let uniform = 1.0 / n.to_double()
  let mut si = 0
  while si < n {
    dist.set(states[si], uniform)
    si = si + 1
  }

  let max_iter = 100
  let tol = 1.0e-8
  let mut iter = 0
  while iter < max_iter {
    let new_dist = markov_empty_map_d()
    let mut ti = 0
    while ti < n {
      let to_state = states[ti]
      let mut total = 0.0
      let mut fi = 0
      while fi < n {
        let from_state = states[fi]
        let d = dist.get(from_state).unwrap_or(0.0)
        let tp = markov_get_trans_prob(model, from_state, to_state)
        total = total + d * tp
        fi = fi + 1
      }
      new_dist.set(to_state, total)
      ti = ti + 1
    }
    let mut delta = 0.0
    let mut ki = 0
    while ki < n {
      let k = states[ki]
      let d1 = dist.get(k).unwrap_or(0.0)
      let d2 = new_dist.get(k).unwrap_or(0.0)
      let diff = d1 - d2
      delta = delta + (if diff < 0.0 { -diff } else { diff })
      ki = ki + 1
    }
    let mut si2 = 0
    while si2 < n {
      dist.set(states[si2], new_dist.get(states[si2]).unwrap_or(0.0))
      si2 = si2 + 1
    }
    if delta < tol {
      iter = max_iter
    }
    iter = iter + 1
  }
  dist
}
