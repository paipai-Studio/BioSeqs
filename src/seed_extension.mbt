///|
/// X-drop ungapped seed extension — extends a seed (diagonal hit) in
/// both directions without gaps, stopping when the running score drops
/// more than `x_drop` below the best score seen so far.
///
/// This is a core component of BLAST and MMseqs2 pipelines: after
/// finding a k-mer seed, the seed is extended ungapped to produce
/// an HSP (High-scoring Segment Pair), which is then passed to the
/// gapped alignment stage.
///
/// Reference: Altschul et al., "Gapped BLAST and PSI-BLAST: a new
/// generation of protein database search programs", NAR 1997.

///|
/// Result of X-drop seed extension.
pub struct SeedExtensionResult {
  query_start : Int32
  query_end : Int32
  target_start : Int32
  target_end : Int32
  score : Double
  aligned_query : String
  aligned_target : String
  identity : Double
}

///|
/// Construct a SeedExtensionResult.
pub fn SeedExtensionResult::new(
  query_start~ : Int32,
  query_end~ : Int32,
  target_start~ : Int32,
  target_end~ : Int32,
  score~ : Double,
  aligned_query~ : String,
  aligned_target~ : String,
  identity~ : Double,
) -> SeedExtensionResult {
  SeedExtensionResult::{
    query_start,
    query_end,
    target_start,
    target_end,
    score,
    aligned_query,
    aligned_target,
    identity,
  }
}

///|
/// Get the alignment score.
pub fn SeedExtensionResult::get_score(self : SeedExtensionResult) -> Double {
  self.score
}

///|
/// Get the percent identity.
pub fn SeedExtensionResult::get_identity(self : SeedExtensionResult) -> Double {
  self.identity
}

///|
/// Get the aligned query string.
pub fn SeedExtensionResult::get_aligned_query(
  self : SeedExtensionResult,
) -> String {
  self.aligned_query
}

///|
/// Get the aligned target string.
pub fn SeedExtensionResult::get_aligned_target(
  self : SeedExtensionResult,
) -> String {
  self.aligned_target
}

///|
/// Get the query range [start, end).
pub fn SeedExtensionResult::get_query_range(
  self : SeedExtensionResult,
) -> (Int32, Int32) {
  (self.query_start, self.query_end)
}

///|
/// Get the target range [start, end).
pub fn SeedExtensionResult::get_target_range(
  self : SeedExtensionResult,
) -> (Int32, Int32) {
  (self.target_start, self.target_end)
}

///|
/// A seed (diagonal hit) to be extended.
pub struct Seed {
  query_pos : Int32
  target_pos : Int32
  length : Int32
}

///|
/// Construct a Seed.
pub fn Seed::new(
  query_pos~ : Int32,
  target_pos~ : Int32,
  length~ : Int32 = 1,
) -> Seed {
  Seed::{ query_pos, target_pos, length, }
}

///|
/// Extend a seed to the right (forward) without gaps, using X-drop.
fn extend_right(
  query : String,
  target : String,
  q_start : Int32,
  t_start : Int32,
  x_drop : Double,
  match_reward : Double,
  mismatch_penalty : Double,
) -> (Int32, Double) {
  let mut q = q_start
  let mut t = t_start
  let mut best_score = 0.0
  let mut current_score = 0.0
  let mut best_len = 0
  while q < query.length() && t < target.length() {
    let qc = query.unsafe_get(q).unsafe_to_char()
    let tc = target.unsafe_get(t).unsafe_to_char()
    if qc == tc {
      current_score = current_score + match_reward
    } else {
      current_score = current_score + mismatch_penalty
    }
    if current_score > best_score {
      best_score = current_score
      best_len = q - q_start + 1
    }
    if best_score - current_score > x_drop {
      break
    }
    q = q + 1
    t = t + 1
  }
  (best_len, best_score)
}

///|
/// Extend a seed to the left (backward) without gaps, using X-drop.
fn extend_left(
  query : String,
  target : String,
  q_start : Int32,
  t_start : Int32,
  x_drop : Double,
  match_reward : Double,
  mismatch_penalty : Double,
) -> (Int32, Double) {
  let mut q = q_start - 1
  let mut t = t_start - 1
  let mut best_score = 0.0
  let mut current_score = 0.0
  let mut best_len = 0
  while q >= 0 && t >= 0 {
    let qc = query.unsafe_get(q).unsafe_to_char()
    let tc = target.unsafe_get(t).unsafe_to_char()
    if qc == tc {
      current_score = current_score + match_reward
    } else {
      current_score = current_score + mismatch_penalty
    }
    if current_score > best_score {
      best_score = current_score
      best_len = q_start - q
    }
    if best_score - current_score > x_drop {
      break
    }
    q = q - 1
    t = t - 1
  }
  (best_len, best_score)
}

///|
/// Extract a substring [start, end).
fn seed_substring(s : String, start : Int32, end : Int32) -> String {
  let mut result = ""
  let mut i = start
  while i < end && i < s.length() {
    result = result + s.unsafe_get(i).unsafe_to_char().to_string()
    i = i + 1
  }
  result
}

///|
/// Extend a seed in both directions using X-drop termination.
///
/// The seed is extended left and right without gaps. Extension stops
/// when the running score drops more than `x_drop` below the best
/// score seen in that direction.
pub fn extend_seed(
  query : String,
  target : String,
  seed : Seed,
  x_drop~ : Double = 20.0,
  match_reward~ : Double = 2.0,
  mismatch_penalty~ : Double = -1.0,
) -> SeedExtensionResult {
  let seed_end_q = seed.query_pos + seed.length
  let seed_end_t = seed.target_pos + seed.length
  let (right_len, right_score) = extend_right(
    query, target, seed_end_q, seed_end_t, x_drop, match_reward, mismatch_penalty,
  )
  let (left_len, left_score) = extend_left(
    query,
    target,
    seed.query_pos,
    seed.target_pos,
    x_drop,
    match_reward,
    mismatch_penalty,
  )
  let q_start = seed.query_pos - left_len
  let t_start = seed.target_pos - left_len
  let q_end = seed_end_q + right_len
  let t_end = seed_end_t + right_len
  let aligned_query = seed_substring(query, q_start, q_end)
  let aligned_target = seed_substring(target, t_start, t_end)
  let mut seed_score = 0.0
  for i in 0..<seed.length {
    let qc = query.unsafe_get(seed.query_pos + i).unsafe_to_char()
    let tc = target.unsafe_get(seed.target_pos + i).unsafe_to_char()
    if qc == tc {
      seed_score = seed_score + match_reward
    } else {
      seed_score = seed_score + mismatch_penalty
    }
  }
  let total_score = left_score + seed_score + right_score
  let aln_len = aligned_query.length()
  let mut matches = 0
  for i in 0..<aln_len {
    if aligned_query.unsafe_get(i) == aligned_target.unsafe_get(i) {
      matches = matches + 1
    }
  }
  let identity = if aln_len > 0 {
    matches.to_double() / aln_len.to_double() * 100.0
  } else {
    0.0
  }
  SeedExtensionResult::new(
    query_start=q_start,
    query_end=q_end,
    target_start=t_start,
    target_end=t_end,
    score=total_score,
    aligned_query~,
    aligned_target~,
    identity~,
  )
}
