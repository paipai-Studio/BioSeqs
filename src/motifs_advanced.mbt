///|
/// Bio.motifs advanced module - JASPAR format, TRANSFAC format, advanced motif comparison.
/// Includes JASPAR database format parsing, TRANSFAC format, motif alignment,
/// pseudocount optimization, and motif clustering.

///|
/// JASPAR motif record.
pub struct JasparMotif {
  matrix_id : String
  name : String
  collection : String
  tf_class : String
  tf_family : String
  species : String
  tax_group : String
  acc : Array[String]
  pwm : Array[Array[Double]]
  alphabet : Array[String]
} derive(Eq)

///|
/// TRANSFAC motif record.
pub struct TransfacMotif {
  accession : String
  id : String
  name : String
  description : String
  pwm : Array[Array[Double]]
  alphabet : Array[String]
  binding_sites : Array[String]
} derive(Eq)

///|
/// Motif alignment result.
pub struct MotifAlignmentResult {
  score : Double
  offset : Int
  strand : String
  aligned_length : Int
} derive(Eq)

///|
/// Clustered motif group.
pub struct MotifCluster {
  cluster_id : Int
  motifs : Array[String]
  representative : String
  avg_similarity : Double
} derive(Eq)

///|
let dna_alphabet_jaspar : Array[String] = ["A", "C", "G", "T"]

///|
pub fn JasparMotif::new(
  matrix_id : String,
  name : String,
  collection : String,
  tf_class : String,
  tf_family : String,
  species : String,
  tax_group : String,
  acc : Array[String],
  pwm : Array[Array[Double]],
  alphabet : Array[String],
) -> JasparMotif {
  JasparMotif::{
    matrix_id,
    name,
    collection,
    tf_class,
    tf_family,
    species,
    tax_group,
    acc,
    pwm,
    alphabet,
  }
}

///|
pub fn TransfacMotif::new(
  accession : String,
  id : String,
  name : String,
  description : String,
  pwm : Array[Array[Double]],
  alphabet : Array[String],
  binding_sites : Array[String],
) -> TransfacMotif {
  TransfacMotif::{
    accession,
    id,
    name,
    description,
    pwm,
    alphabet,
    binding_sites,
  }
}

///|
pub fn MotifAlignmentResult::new(
  score : Double,
  offset : Int,
  strand : String,
  aligned_length : Int,
) -> MotifAlignmentResult {
  MotifAlignmentResult::{ score, offset, strand, aligned_length }
}

///|
pub fn MotifCluster::new(
  cluster_id : Int,
  motifs : Array[String],
  representative : String,
  avg_similarity : Double,
) -> MotifCluster {
  MotifCluster::{ cluster_id, motifs, representative, avg_similarity }
}

// ============================================================================
// Helper functions
// ============================================================================

///|
fn adv_parse_float(s : String) -> Double? {
  if s == "." || s.length() == 0 {
    return None
  }
  let mut result = 0.0
  let mut sign = 1.0
  let mut i = 0
  let zero = '0'.to_int()
  
  if s.unsafe_get(0).to_int() == '-'.to_int() {
    sign = -1.0
    i = 1
  }
  
  let mut has_dot = false
  let mut decimal_place = 1.0
  
  while i < s.length() {
    let c = s.unsafe_get(i).to_int()
    if c == '.'.to_int() {
      has_dot = true
      i = i + 1
      continue
    }
    if c >= zero && c <= zero + 9 {
      let digit = (c - zero).to_double()
      result = result * 10.0 + digit
      if has_dot {
        decimal_place = decimal_place * 10.0
      }
    }
    i = i + 1
  }
  
  Some(result * sign / decimal_place)
}

///|
fn adv_log2(x : Double) -> Double {
  if x <= 0.0 { return 0.0 }
  let ln2 = 0.6931471805599453
  let term = (x - 1.0) / (x + 1.0)
  let term2 = term * term
  let mut sum = term
  let mut cur_term = term
  let mut i = 1
  while i <= 50 {
    cur_term = cur_term * term2
    sum = sum + cur_term / (2 * i + 1).to_double()
    i = i + 1
  }
  let ln_x = 2.0 * sum
  ln_x / ln2
}

///|
fn adv_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  let mut guess = x / 2.0
  let mut i = 0
  while i < 20 {
    guess = (guess + x / guess) / 2.0
    i = i + 1
  }
  guess
}

///|
fn adv_index_of(s : String, sub : String) -> Int {
  let s_len = s.length()
  let sub_len = sub.length()
  if sub_len == 0 || s_len < sub_len {
    return -1
  }
  let mut i = 0
  while i <= s_len - sub_len {
    let mut matched = true
    let mut j = 0
    while j < sub_len && matched {
      if s.unsafe_get(i + j) != sub.unsafe_get(j) {
        matched = false
      }
      j = j + 1
    }
    if matched {
      return i
    }
    i = i + 1
  }
  -1
}

///|
fn split_whitespace(s : String) -> Array[String] {
  let result : Array[String] = Array::new()
  let len = s.length()
  let mut i = 0
  while i < len {
    // Skip whitespace
    while i < len && (s.unsafe_get(i) == ' ' || s.unsafe_get(i) == '\t') {
      i = i + 1
    }
    if i >= len {
      break
    }
    // Find end of token
    let start = i
    while i < len && s.unsafe_get(i) != ' ' && s.unsafe_get(i) != '\t' {
      i = i + 1
    }
    let token = s[start:i].to_owned()
    result.push(token)
  }
  result
}

// ============================================================================
// JASPAR format parsing
// ============================================================================

///|
/// Parse a JASPAR PFM (Position Frequency Matrix) format string.
pub fn parse_jaspar_pfm(content : String) -> Array[JasparMotif] {
  let motifs : Array[JasparMotif] = Array::new()
  let lines_view = content.split("\n")
  let lines = lines_view.to_array()
  
  let mut i = 0
  while i < lines.length() {
    let line_view = lines[i].trim()
    let line = line_view.to_owned()
    if line.has_prefix(">") {
      let header = line[1:].to_owned()
      let parts_view = header.split(" ")
      let parts = parts_view.to_array()
      let mut matrix_id = ""
      let mut name = ""
      if parts.length() >= 1 {
        matrix_id = parts[0].trim().to_owned()
      }
      if parts.length() >= 2 {
        name = parts[1].trim().to_owned()
      }
      
      let mut a_counts : Array[Double] = Array::new()
      let mut c_counts : Array[Double] = Array::new()
      let mut g_counts : Array[Double] = Array::new()
      let mut t_counts : Array[Double] = Array::new()
      
      let mut j = i + 1
      while j < lines.length() && !lines[j].trim().has_prefix(">") {
        let row_view = lines[j].trim()
        let row_line = row_view.to_owned()
        if row_line.length() > 0 {
          let (base, counts) = parse_jaspar_row(row_line)
          if base == "A" {
            a_counts = counts
          } else if base == "C" {
            c_counts = counts
          } else if base == "G" {
            g_counts = counts
          } else if base == "T" {
            t_counts = counts
          }
        }
        j = j + 1
      }
      
      let width = a_counts.length()
      if width > 0 {
        let pwm = Array::make(width, Array::make(4, 0.0))
        for pos = 0; pos < width; pos = pos + 1 {
          pwm[pos][0] = a_counts[pos]
          pwm[pos][1] = c_counts[pos]
          pwm[pos][2] = g_counts[pos]
          pwm[pos][3] = t_counts[pos]
        }
        
        let motif = JasparMotif::new(
          matrix_id,
          name,
          "",
          "",
          "",
          "",
          "",
          [],
          pwm,
          dna_alphabet_jaspar,
        )
        motifs.push(motif)
      }
      
      i = j
    } else {
      i = i + 1
    }
  }
  
  motifs
}

///|
fn parse_jaspar_row(line : String) -> (String, Array[Double]) {
  let counts : Array[Double] = Array::new()
  let mut base = ""
  
  if line.length() >= 1 {
    base = line[0:1].to_owned().to_upper()
  }
  
  let start_bracket = adv_index_of(line, "[")
  let end_bracket = adv_index_of(line, "]")
  
  if start_bracket >= 0 && end_bracket > start_bracket {
    let numbers_str = line[start_bracket + 1:end_bracket].to_owned()
    let number_parts = split_whitespace(numbers_str)
    for k = 0; k < number_parts.length(); k = k + 1 {
      let num_str = number_parts[k]
      if num_str.length() > 0 {
        let parsed = adv_parse_float(num_str)
        match parsed {
          Some(v) => counts.push(v)
          None => counts.push(0.0)
        }
      }
    }
  } else {
    let parts = split_whitespace(line)
    if parts.length() > 1 {
      for k = 1; k < parts.length(); k = k + 1 {
        let num_str = parts[k]
        if num_str.length() > 0 {
          let parsed = adv_parse_float(num_str)
          match parsed {
            Some(v) => counts.push(v)
            None => counts.push(0.0)
          }
        }
      }
    }
  }
  
  (base, counts)
}

///|
/// Convert JASPAR PFM (counts) to PWM (probabilities) with pseudocount.
pub fn jaspar_pfm_to_pwm(
  pfm : Array[Array[Double]],
  pseudocount : Double,
) -> Array[Array[Double]] {
  if pfm.length() == 0 {
    return []
  }
  let width = pfm.length()
  let alphabet_size = pfm[0].length()
  let pwm = Array::make(width, Array::make(alphabet_size, 0.0))
  
  for pos = 0; pos < width; pos = pos + 1 {
    let row = pfm[pos]
    let mut total = 0.0
    for i = 0; i < row.length(); i = i + 1 {
      total = total + row[i]
    }
    let norm_total = total + pseudocount * alphabet_size.to_double()
    
    for alpha_idx = 0; alpha_idx < alphabet_size; alpha_idx = alpha_idx + 1 {
      pwm[pos][alpha_idx] = (row[alpha_idx] + pseudocount) / norm_total
    }
  }
  
  pwm
}

// ============================================================================
// TRANSFAC format parsing
// ============================================================================

///|
/// Parse TRANSFAC matrix format.
pub fn parse_transfac(content : String) -> Array[TransfacMotif] {
  let motifs : Array[TransfacMotif] = Array::new()
  let lines_view = content.split("\n")
  let lines = lines_view.to_array()
  
  let mut i = 0
  while i < lines.length() {
    let line_view = lines[i].trim()
    let line = line_view.to_owned()
    if line.has_prefix("AC ") {
      let accession = line[2:].trim().to_owned()
      let mut id = ""
      let mut name = ""
      let mut description = ""
      let matrix_rows : Array[Array[Double]] = Array::new()
      let binding_sites : Array[String] = Array::new()
      let mut in_matrix = false
      
      let mut j = i + 1
      while j < lines.length() && !lines[j].trim().has_prefix("//") {
        let row_view = lines[j].trim()
        let row_line = row_view.to_owned()
        
        if row_line.has_prefix("ID ") {
          id = row_line[2:].trim().to_owned()
        } else if row_line.has_prefix("NA ") {
          name = row_line[2:].trim().to_owned()
        } else if row_line.has_prefix("DE ") {
          description = row_line[2:].trim().to_owned()
        } else if row_line.has_prefix("P0 ") {
          in_matrix = true
        } else if in_matrix && row_line.length() > 0 {
          let parts = split_whitespace(row_line)
          if parts.length() >= 5 {
            let row : Array[Double] = Array::new()
            for k = 1; k < 5 && k < parts.length(); k = k + 1 {
              let parsed = adv_parse_float(parts[k])
              match parsed {
                Some(v) => row.push(v)
                None => row.push(0.0)
              }
            }
            matrix_rows.push(row)
          }
        }
        
        j = j + 1
      }
      
      if matrix_rows.length() > 0 {
        let width = matrix_rows.length()
        let pwm = Array::make(width, Array::make(4, 0.0))
        for pos = 0; pos < width; pos = pos + 1 {
          pwm[pos][0] = matrix_rows[pos][0]
          pwm[pos][1] = matrix_rows[pos][1]
          pwm[pos][2] = matrix_rows[pos][2]
          pwm[pos][3] = matrix_rows[pos][3]
        }
        
        let motif = TransfacMotif::new(
          accession,
          id,
          name,
          description,
          pwm,
          dna_alphabet_jaspar,
          binding_sites,
        )
        motifs.push(motif)
      }
      
      i = j + 1
    } else {
      i = i + 1
    }
  }
  
  motifs
}

// ============================================================================
// Advanced motif comparison
// ============================================================================

///|
/// Find the optimal alignment between two motifs using sliding window approach.
pub fn optimal_motif_alignment(
  motif1_pwm : Array[Array[Double]],
  motif2_pwm : Array[Array[Double]],
) -> MotifAlignmentResult {
  if motif1_pwm.length() == 0 || motif2_pwm.length() == 0 {
    return MotifAlignmentResult::new(0.0, 0, "+", 0)
  }
  
  let len1 = motif1_pwm.length()
  let len2 = motif2_pwm.length()
  let mut best_score = -1.0
  let mut best_offset = 0
  let mut best_strand = "+"
  let mut best_aligned_len = 0
  
  for offset = -(len2 - 1); offset < len1; offset = offset + 1 {
    let (score, aligned_len) = calculate_alignment_score(
      motif1_pwm,
      motif2_pwm,
      offset,
    )
    if score > best_score || (score == best_score && aligned_len > best_aligned_len) {
      best_score = score
      best_offset = offset
      best_strand = "+"
      best_aligned_len = aligned_len
    }
  }
  
  let rev_motif2 = reverse_complement_pwm(motif2_pwm)
  for offset = -(len2 - 1); offset < len1; offset = offset + 1 {
    let (score, aligned_len) = calculate_alignment_score(
      motif1_pwm,
      rev_motif2,
      offset,
    )
    if score > best_score || (score == best_score && aligned_len > best_aligned_len) {
      best_score = score
      best_offset = offset
      best_strand = "-"
      best_aligned_len = aligned_len
    }
  }
  
  MotifAlignmentResult::new(best_score, best_offset, best_strand, best_aligned_len)
}

///|
fn calculate_alignment_score(
  pwm1 : Array[Array[Double]],
  pwm2 : Array[Array[Double]],
  offset : Int,
) -> (Double, Int) {
  let len1 = pwm1.length()
  let len2 = pwm2.length()
  let mut score = 0.0
  let mut aligned_len = 0
  
  let start1 = if offset >= 0 { offset } else { 0 }
  let start2 = if offset >= 0 { 0 } else { -offset }
  
  let mut i = start1
  let mut j = start2
  while i < len1 && j < len2 {
    let row1 = pwm1[i]
    let row2 = pwm2[j]
    let mut pos_score = 0.0
    for k = 0; k < row1.length() && k < row2.length(); k = k + 1 {
      pos_score = pos_score + row1[k] * row2[k]
    }
    score = score + pos_score
    aligned_len = aligned_len + 1
    i = i + 1
    j = j + 1
  }
  
  if aligned_len > 0 {
    score = score / aligned_len.to_double()
  }
  
  (score, aligned_len)
}

///|
fn reverse_complement_pwm(pwm : Array[Array[Double]]) -> Array[Array[Double]] {
  let width = pwm.length()
  if width == 0 {
    return []
  }
  let alphabet_size = pwm[0].length()
  let rev_pwm = Array::make(width, Array::make(alphabet_size, 0.0))
  
  for pos = 0; pos < width; pos = pos + 1 {
    let orig_pos = width - 1 - pos
    if alphabet_size == 4 {
      rev_pwm[pos][0] = pwm[orig_pos][3]
      rev_pwm[pos][1] = pwm[orig_pos][2]
      rev_pwm[pos][2] = pwm[orig_pos][1]
      rev_pwm[pos][3] = pwm[orig_pos][0]
    } else {
      for k = 0; k < alphabet_size; k = k + 1 {
        rev_pwm[pos][k] = pwm[orig_pos][k]
      }
    }
  }
  
  rev_pwm
}

///|
/// Calculate the Kullback-Leibler divergence between two PWMs.
pub fn motif_kl_divergence(
  pwm1 : Array[Array[Double]],
  pwm2 : Array[Array[Double]],
) -> Double {
  if pwm1.length() == 0 || pwm2.length() == 0 {
    return 0.0
  }
  
  let min_len = if pwm1.length() < pwm2.length() {
    pwm1.length()
  } else {
    pwm2.length()
  }
  
  let mut total_kl = 0.0
  for pos = 0; pos < min_len; pos = pos + 1 {
    let row1 = pwm1[pos]
    let row2 = pwm2[pos]
    let mut pos_kl = 0.0
    for k = 0; k < row1.length() && k < row2.length(); k = k + 1 {
      let p = row1[k]
      let q = row2[k]
      if p > 0.0 && q > 0.0 {
        pos_kl = pos_kl + p * adv_log2(p / q)
      }
    }
    total_kl = total_kl + pos_kl
  }
  
  total_kl / min_len.to_double()
}

///|
/// Calculate Jensen-Shannon divergence between two PWMs.
pub fn motif_js_divergence(
  pwm1 : Array[Array[Double]],
  pwm2 : Array[Array[Double]],
) -> Double {
  if pwm1.length() == 0 || pwm2.length() == 0 {
    return 0.0
  }
  
  let min_len = if pwm1.length() < pwm2.length() {
    pwm1.length()
  } else {
    pwm2.length()
  }
  
  let mut total_js = 0.0
  for pos = 0; pos < min_len; pos = pos + 1 {
    let row1 = pwm1[pos]
    let row2 = pwm2[pos]
    let alphabet_size = row1.length()
    let m_row = Array::make(alphabet_size, 0.0)
    
    for k = 0; k < alphabet_size; k = k + 1 {
      m_row[k] = (row1[k] + row2[k]) / 2.0
    }
    
    let mut kl1 = 0.0
    let mut kl2 = 0.0
    for k = 0; k < alphabet_size; k = k + 1 {
      if row1[k] > 0.0 && m_row[k] > 0.0 {
        kl1 = kl1 + row1[k] * adv_log2(row1[k] / m_row[k])
      }
      if row2[k] > 0.0 && m_row[k] > 0.0 {
        kl2 = kl2 + row2[k] * adv_log2(row2[k] / m_row[k])
      }
    }
    
    total_js = total_js + (kl1 + kl2) / 2.0
  }
  
  total_js / min_len.to_double()
}

// ============================================================================
// Motif clustering
// ============================================================================

///|
/// Cluster motifs based on pairwise similarity using simple hierarchical clustering.
pub fn cluster_motifs(
  motif_names : Array[String],
  motif_pwms : Array[Array[Array[Double]]],
  threshold : Double,
) -> Array[MotifCluster] {
  let n = motif_names.length()
  if n == 0 {
    return []
  }
  
  let similarity_matrix = Array::make(n, Array::make(n, 0.0))
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      if i == j {
        similarity_matrix[i][j] = 1.0
      } else if j < i {
        similarity_matrix[i][j] = similarity_matrix[j][i]
      } else {
        let alignment = optimal_motif_alignment(motif_pwms[i], motif_pwms[j])
        similarity_matrix[i][j] = alignment.score
      }
    }
  }
  
  let clusters : Array[MotifCluster] = Array::new()
  let assigned = Array::make(n, false)
  
  for i = 0; i < n; i = i + 1 {
    if !assigned[i] {
      let cluster_members : Array[String] = Array::new()
      let cluster_indices : Array[Int] = Array::new()
      
      for j = i; j < n; j = j + 1 {
        if !assigned[j] {
          let mut all_similar = true
          for k = 0; k < cluster_indices.length() && all_similar; k = k + 1 {
            if similarity_matrix[j][cluster_indices[k]] < threshold {
              all_similar = false
            }
          }
          if all_similar {
            cluster_members.push(motif_names[j])
            cluster_indices.push(j)
            assigned[j] = true
          }
        }
      }
      
      let mut best_rep_idx = 0
      let mut best_avg_sim = 0.0
      for k = 0; k < cluster_indices.length(); k = k + 1 {
        let mut avg_sim = 0.0
        let idx = cluster_indices[k]
        for m = 0; m < cluster_indices.length(); m = m + 1 {
          avg_sim = avg_sim + similarity_matrix[idx][cluster_indices[m]]
        }
        avg_sim = avg_sim / cluster_indices.length().to_double()
        if avg_sim > best_avg_sim {
          best_avg_sim = avg_sim
          best_rep_idx = k
        }
      }
      
      let cluster = MotifCluster::new(
        clusters.length(),
        cluster_members,
        cluster_members[best_rep_idx],
        best_avg_sim,
      )
      clusters.push(cluster)
    }
  }
  
  clusters
}

// ============================================================================
// Pseudocount optimization
// ============================================================================

///|
/// Calculate optimal pseudocount based on sequence count and information content.
pub fn optimal_pseudocount(num_sequences : Int, width : Int) -> Double {
  let n = num_sequences.to_double()
  let w = width.to_double()
  let sqrt_n = adv_sqrt(n)
  sqrt_n / (n + w)
}

///|
/// Calculate motif GC content.
pub fn motif_gc_content_adv(pwm : Array[Array[Double]]) -> Double {
  if pwm.length() == 0 {
    return 0.0
  }
  let mut gc = 0.0
  for pos = 0; pos < pwm.length(); pos = pos + 1 {
    let row = pwm[pos]
    if row.length() >= 4 {
      gc = gc + row[1] + row[2]
    }
  }
  gc / pwm.length().to_double()
}

// ============================================================================
// Motif shuffling and randomization
// ============================================================================

///|
/// Shuffle the columns of a PWM to create a randomized control motif.
pub fn shuffle_motif_columns(pwm : Array[Array[Double]]) -> Array[Array[Double]] {
  let width = pwm.length()
  if width == 0 {
    return []
  }
  let alphabet_size = pwm[0].length()
  
  let indices = Array::make(width, 0)
  for i = 0; i < width; i = i + 1 {
    indices[i] = i
  }
  
  for i = width - 1; i > 0; i = i - 1 {
    let j = (i * 7 + 3) % (i + 1)
    let temp = indices[i]
    indices[i] = indices[j]
    indices[j] = temp
  }
  
  let shuffled = Array::make(width, Array::make(alphabet_size, 0.0))
  for i = 0; i < width; i = i + 1 {
    shuffled[i] = pwm[indices[i]].copy()
  }
  
  shuffled
}

///|
/// Create a motif with randomized nucleotide frequencies for each position.
pub fn randomize_motif(width : Int, alphabet_size : Int) -> Array[Array[Double]] {
  let pwm = Array::make(width, Array::make(alphabet_size, 0.0))
  
  for pos = 0; pos < width; pos = pos + 1 {
    let row : Array[Double] = Array::new()
    let mut total = 0.0
    
    for i = 0; i < alphabet_size; i = i + 1 {
      let val = adv_sin((pos * 17 + i * 31).to_double()).abs() + 0.1
      row.push(val)
      total = total + val
    }
    
    for i = 0; i < alphabet_size; i = i + 1 {
      pwm[pos][i] = row[i] / total
    }
  }
  
  pwm
}

///|
fn adv_sin(x : Double) -> Double {
  let mut result = 0.0
  let mut term = x
  let mut sign = 1.0
  let mut factorial = 1.0
  let mut power = x
  
  let n = 15
  let mut i = 0
  while i < n {
    result = result + sign * power / factorial
    sign = -sign
    power = power * x * x
    factorial = factorial * (2 * i + 2).to_double() * (2 * i + 3).to_double()
    i = i + 1
  }
  
  result
}

// ============================================================================
// Motif format conversion helpers
// ============================================================================

///|
/// Convert PWM to TRANSFAC matrix format string.
pub fn pwm_to_transfac(
  pwm : Array[Array[Double]],
  name : String,
  accession : String,
) -> String {
  if pwm.length() == 0 {
    return ""
  }
  
  let mut result = "AC " + accession + "\n"
  result = result + "ID " + name + "\n"
  result = result + "NA " + name + "\n"
  result = result + "DE " + name + " matrix\n"
  result = result + "P0      A      C      G      T\n"
  
  for pos = 0; pos < pwm.length(); pos = pos + 1 {
    let row = pwm[pos]
    let pos_str = if pos < 10 {
      "0" + pos.to_string()
    } else {
      pos.to_string()
    }
    
    let mut row_str = pos_str + "    "
    for k = 0; k < row.length(); k = k + 1 {
      let val_str = row[k].to_string()
      if val_str.length() < 6 {
        let pad = 6 - val_str.length()
        let mut padding = ""
        let mut p = 0
        while p < pad {
          padding = padding + " "
          p = p + 1
        }
        row_str = row_str + padding + val_str
      } else {
        row_str = row_str + " " + val_str[0:5].to_owned()
      }
    }
    result = result + row_str + "\n"
  }
  
  result = result + "XX\n"
  result = result + "//\n"
  
  result
}

///|
/// Convert PWM to JASPAR PFM format string.
pub fn pwm_to_jaspar(
  pwm : Array[Array[Double]],
  matrix_id : String,
  name : String,
) -> String {
  if pwm.length() == 0 {
    return ""
  }
  
  let width = pwm.length()
  let mut result = ">" + matrix_id + " " + name + "\n"
  
  let bases = ["A", "C", "G", "T"]
  for b = 0; b < 4; b = b + 1 {
    let mut values_str = bases[b] + " [ "
    for pos = 0; pos < width; pos = pos + 1 {
      let val = pwm[pos][b]
      let int_val = (val * 100.0).to_int()
      values_str = values_str + int_val.to_string() + " "
    }
    values_str = values_str + "]\n"
    result = result + values_str
  }
  
  result
}
