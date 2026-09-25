///|
/// motif_clusterbuster - ClusterBuster format parsing for clustered motifs.
/// Inspired by Biopython Bio.motifs (clusterbuster support).
/// Format:
///   >cluster_name score
///   motif_name
///   [p_A p_C p_G p_T]
///   [p_A p_C p_G p_T]
///   ...

///|
/// A motif defined in ClusterBuster format: name + PFM (probabilities).
pub struct ClusterBusterMotif {
  name : String
  pfm : Array[Array[Double]] // columns x 4 (A C G T)
}

///|
/// A cluster of motifs with an aggregate score.
pub struct ClusterBusterCluster {
  name : String
  score : Double
  motifs : Array[ClusterBusterMotif]
}

///|
/// Parse a string of doubles inside brackets "[a b c d]".
fn clusterbuster_parse_probs(line : String) -> Array[Double]? {
  let trimmed = line.trim(chars=" \t\r").to_owned()
  if trimmed.length() < 2 || trimmed[0] != 91 || trimmed[trimmed.length() - 1] != 93 {
    return None // '[' ... ']'
  }
  let inner = trimmed[1:trimmed.length() - 1].to_owned()
  let fields = inner.split(" ").to_array()
  let probs : Array[Double] = []
  let mut i = 0
  while i < fields.length() {
    let f = fields[i].trim(chars=" \t").to_owned()
    if f.length() > 0 {
      match clusterbuster_atof_opt(f) {
        Some(v) => probs.push(v)
        None => return None
      }
    }
    i = i + 1
  }
  if probs.length() == 4 {
    Some(probs)
  } else {
    None
  }
}

///|
/// Parse double or None.
fn clusterbuster_atof_opt(s : String) -> Double? {
  let t = s.trim(chars=" \t\r").to_owned()
  if t.length() == 0 {
    return None
  }
  let neg = t[0] == 45
  let start = if t[0] == 45 || t[0] == 43 { 1 } else { 0 }
  let mut int_part = 0.0
  let mut frac_part = 0.0
  let mut in_frac = false
  let mut frac_div = 10.0
  let mut any_digit = false
  let mut k = start
  while k < t.length() {
    let c = t[k]
    if c >= 48 && c <= 57 {
      any_digit = true
      let digit = (c - 48).to_int().to_double()
      if in_frac {
        frac_part = frac_part + digit / frac_div
        frac_div = frac_div * 10.0
      } else {
        int_part = int_part * 10.0 + digit
      }
    } else if c == 46 && !(in_frac) {
      in_frac = true
    } else if c == 101 || c == 69 {
      // scientific notation exponent: parse rest as int
      let exp_str = t[k + 1:t.length()].to_owned()
      match clusterbuster_atof_opt(exp_str) {
        Some(exp_v) => {
          if !(any_digit) {
            return None
          }
          let base = int_part + frac_part
          let value = base * @math.pow(10.0, exp_v)
          return Some(if neg { -value } else { value })
        }
        None => return None
      }
    } else if c == 45 && k == start {
      // handled above
      continue
    } else if c == 43 && k == start {
      continue
    } else {
      return None
    }
    k = k + 1
  }
  if !(any_digit) {
    return None
  }
  let value = int_part + frac_part
  Some(if neg { -value } else { value })
}

///|
/// Parse ClusterBuster format content into clusters.
pub fn clusterbuster_parse(content : String) -> Array[ClusterBusterCluster] {
  let clusters : Array[ClusterBusterCluster] = []
  let lines = content.split("\n").to_array()
  let mut current_cluster : ClusterBusterCluster? = None
  let mut current_motif : ClusterBusterMotif? = None
  let mut i = 0
  while i < lines.length() {
    let raw = lines[i]
    let line = raw.trim(chars=" \t\r").to_owned()
    if line.length() == 0 || line[0] == 35 {
      // skip empty and comment lines
      i = i + 1
      continue
    }
    if line[0] == 62 {
      // '>' cluster header
      // finalize previous motif
      match current_cluster {
        Some(cl) => {
          match current_motif {
            Some(m) => {
              cl.motifs.push(m)
              current_motif = None
            }
            None => continue
          }
          clusters.push(cl)
        }
        None => ()
      }
      // parse header: >name score
      let header = line[1:line.length()].to_owned()
      let fields = header.split(" ").to_array()
      let name = if fields.length() > 0 { fields[0].to_owned() } else { "" }
      let score = if fields.length() > 1 {
        match clusterbuster_atof_opt(fields[1].to_owned()) {
          Some(v) => v
          None => 0.0
        }
      } else {
        0.0
      }
      current_cluster = Some(ClusterBusterCluster::{ name, score, motifs: [] })
    } else {
      // Either motif name or PFM row
      match clusterbuster_parse_probs(line) {
        Some(probs) => {
          // PFM row belongs to current motif
          match current_motif {
            Some(m) => {
              m.pfm.push(probs)
              current_motif = Some(m)
            }
            None => ()
          }
        }
        None => {
          // Motif name line: finalize previous motif first
          match current_cluster {
            Some(cl) => {
              match current_motif {
                Some(m) => cl.motifs.push(m)
                None => ()
              }
            }
            None => ()
          }
          current_motif = Some(ClusterBusterMotif::{ name: line, pfm: [] })
        }
      }
    }
    i = i + 1
  }
  // Finalize last motif and cluster
  match current_cluster {
    Some(cl) => {
      match current_motif {
        Some(m) => cl.motifs.push(m)
        None => ()
      }
      clusters.push(cl)
    }
    None => ()
  }
  clusters
}

///|
/// Serialize clusters to ClusterBuster format.
pub fn clusterbuster_to_string(clusters : Array[ClusterBusterCluster]) -> String {
  let sb = StringBuilder()
  let mut c = 0
  while c < clusters.length() {
    let cl = clusters[c]
    sb.write_string(">")
    sb.write_string(cl.name)
    sb.write_string(" ")
    sb.write_string(cl.score.to_string())
    sb.write_string("\n")
    let mut m = 0
    while m < cl.motifs.length() {
      let motif = cl.motifs[m]
      sb.write_string(motif.name)
      sb.write_string("\n")
      let mut p = 0
      while p < motif.pfm.length() {
        sb.write_string("[")
        let mut a = 0
        while a < motif.pfm[p].length() {
          if a > 0 {
            sb.write_string(" ")
          }
          sb.write_string(motif.pfm[p][a].to_string())
          a = a + 1
        }
        sb.write_string("]\n")
        p = p + 1
      }
      m = m + 1
    }
    c = c + 1
  }
  sb.to_string()
}

///|
/// Convert a ClusterBuster motif PFM to MotifPFM (counts = probs * 100).
pub fn clusterbuster_to_motif_pfm(motif : ClusterBusterMotif) -> MotifPFM {
  let alphabet = ["A", "C", "G", "T"]
  let counts : Array[Array[Double]] = []
  let mut p = 0
  while p < motif.pfm.length() {
    let row : Array[Double] = Array::make(4, 0.0)
    let mut a = 0
    while a < 4 {
      row[a] = motif.pfm[p][a] * 100.0
      a = a + 1
    }
    counts.push(row)
    p = p + 1
  }
  MotifPFM::new(alphabet, counts)
}

///|
/// Convert a MotifPFM to a ClusterBusterMotif (probabilities = counts / total).
pub fn motif_pfm_to_clusterbuster(name : String, pfm : MotifPFM) -> ClusterBusterMotif {
  let pfm_matrix : Array[Array[Double]] = []
  let mut p = 0
  while p < pfm.counts.length() {
    let total = pfm.counts[p][0] + pfm.counts[p][1] + pfm.counts[p][2] + pfm.counts[p][3]
    let row : Array[Double] = Array::make(4, 0.0)
    let mut a = 0
    while a < 4 {
      row[a] = if total > 0.0 { pfm.counts[p][a] / total } else { 0.25 }
      a = a + 1
    }
    pfm_matrix.push(row)
    p = p + 1
  }
  ClusterBusterMotif::{ name, pfm: pfm_matrix }
}

///|
/// Compute cluster score: sum of max-probability per position across motifs.
pub fn clusterbuster_cluster_score(cluster : ClusterBusterCluster) -> Double {
  let mut total = 0.0
  let mut m = 0
  while m < cluster.motifs.length() {
    let motif = cluster.motifs[m]
    let mut p = 0
    while p < motif.pfm.length() {
      let mut max_v = 0.0
      let mut a = 0
      while a < motif.pfm[p].length() {
        if motif.pfm[p][a] > max_v {
          max_v = motif.pfm[p][a]
        }
        a = a + 1
      }
      total = total + max_v
      p = p + 1
    }
    m = m + 1
  }
  total
}

///|
/// Example ClusterBuster content.
pub fn clusterbuster_example() -> Array[ClusterBusterCluster] {
  clusterbuster_parse(">cluster1 2.5\nmotifA\n[0.7 0.1 0.1 0.1]\n[0.1 0.1 0.7 0.1]\n")
}