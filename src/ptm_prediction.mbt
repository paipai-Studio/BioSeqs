///|
/// Post-Translational Modification (PTM) site prediction from protein sequence.
///
/// Biopython equivalent: Not directly in Biopython, but commonly used alongside
/// Bio.ExPASy for protein analysis. This module implements sequence-based
/// PTM site prediction using position-specific scoring and rule-based approaches.
///
/// Covered PTM types:
///   - Phosphorylation (Ser/Thr/Tyr)
///   - N-glycosylation (Asn-X-Ser/Thr, X != Pro)
///   - O-glycosylation (Ser/Thr)
///   - Acetylation (Lys)
///   - Ubiquitination (Lys)
///   - Sumoylation (Lys) — psi-K-X-E motif
///   - Methylation (Lys/Arg)
///   - S-nitrosylation (Cys)
///   - Disulfide bond (Cys-Cys)
pub struct PtmSite {
  position : Int32
  residue : String
  ptm_type : String
  score : Double
  motif : String
}

///|
pub fn PtmSite::position(self : PtmSite) -> Int32 {
  self.position
}

///|
pub fn PtmSite::residue(self : PtmSite) -> String {
  self.residue
}

///|
pub fn PtmSite::ptm_type(self : PtmSite) -> String {
  self.ptm_type
}

///|
pub fn PtmSite::score(self : PtmSite) -> Double {
  self.score
}

///|
pub fn PtmSite::motif(self : PtmSite) -> String {
  self.motif
}

///|
pub struct PtmResult {
  sequence : String
  sites : Array[PtmSite]
  site_count : Int32
}

///|
pub fn PtmResult::sequence(self : PtmResult) -> String {
  self.sequence
}

///|
pub fn PtmResult::sites(self : PtmResult) -> Array[PtmSite] {
  self.sites
}

///|
pub fn PtmResult::site_count(self : PtmResult) -> Int32 {
  self.site_count
}

///|
fn get_window(seq : String, pos : Int32, radius : Int32) -> String {
  let start = if pos - radius < 0 { 0 } else { pos - radius }
  let end = if pos + radius + 1 > seq.length() {
    seq.length()
  } else {
    pos + radius + 1
  }
  seq[start:end].to_owned()
}

///|
fn predict_phosphorylation(seq : String) -> Array[PtmSite] {
  let sites : Array[PtmSite] = []
  let n = seq.length()
  let mut i = 0
  while i < n {
    let aa = seq[i:i + 1].to_owned()
    if aa == "S" || aa == "T" || aa == "Y" {
      let window = get_window(seq, i, 5)
      let mut score = 0.5
      if i + 1 < n && seq[i + 1:i + 2].to_owned() == "P" {
        score = score + 0.2
      }
      if i - 1 >= 0 && seq[i - 1:i].to_owned() == "P" {
        score = score + 0.15
      }
      let mut j = 1
      while j <= 3 && i + j < n {
        let next_aa = seq[i + j:i + j + 1].to_owned()
        if next_aa == "D" || next_aa == "E" {
          score = score + 0.1
        }
        j = j + 1
      }
      let mut k = 2
      while k <= 3 && i - k >= 0 {
        let prev_aa = seq[i - k:i - k + 1].to_owned()
        if prev_aa == "R" || prev_aa == "K" {
          score = score + 0.1
        }
        k = k + 1
      }
      if aa == "Y" {
        let mut m = -3
        while m <= 3 && i + m >= 0 && i + m < n {
          if m != 0 {
            let nearby = seq[i + m:i + m + 1].to_owned()
            if nearby == "E" || nearby == "D" {
              score = score + 0.05
            }
          }
          m = m + 1
        }
      }
      if score > 1.0 {
        score = 1.0
      }
      sites.push({
        position: i,
        residue: aa,
        ptm_type: "phosphorylation",
        score,
        motif: window,
      })
    }
    i = i + 1
  }
  sites
}

///|
fn predict_n_glycosylation(seq : String) -> Array[PtmSite] {
  let sites : Array[PtmSite] = []
  let n = seq.length()
  let mut i = 0
  while i + 2 < n {
    let aa0 = seq[i:i + 1].to_owned()
    let aa1 = seq[i + 1:i + 2].to_owned()
    let aa2 = seq[i + 2:i + 3].to_owned()
    if aa0 == "N" && aa1 != "P" && (aa2 == "S" || aa2 == "T") {
      let motif = seq[i:i + 3].to_owned()
      sites.push({
        position: i,
        residue: "N",
        ptm_type: "n_glycosylation",
        score: 0.9,
        motif,
      })
    }
    i = i + 1
  }
  sites
}

///|
fn predict_o_glycosylation(seq : String) -> Array[PtmSite] {
  let sites : Array[PtmSite] = []
  let n = seq.length()
  let mut i = 0
  while i < n {
    let aa = seq[i:i + 1].to_owned()
    if aa == "S" || aa == "T" {
      let window = get_window(seq, i, 3)
      let mut score = 0.3
      let mut pro_count = 0
      let mut j = -3
      while j <= 3 && i + j >= 0 && i + j < n {
        if j != 0 && seq[i + j:i + j + 1].to_owned() == "P" {
          pro_count = pro_count + 1
        }
        j = j + 1
      }
      score = score + pro_count.to_double() * 0.15
      if score > 1.0 {
        score = 1.0
      }
      if score > 0.45 {
        sites.push({
          position: i,
          residue: aa,
          ptm_type: "o_glycosylation",
          score,
          motif: window,
        })
      }
    }
    i = i + 1
  }
  sites
}

///|
fn predict_acetylation(seq : String) -> Array[PtmSite] {
  let sites : Array[PtmSite] = []
  let n = seq.length()
  let mut i = 0
  while i < n {
    let aa = seq[i:i + 1].to_owned()
    if aa == "K" {
      let window = get_window(seq, i, 4)
      let mut score = 0.4
      if i + 1 < n {
        let next = seq[i + 1:i + 2].to_owned()
        if next == "G" || next == "A" || next == "S" {
          score = score + 0.2
        }
      }
      if i - 1 >= 0 {
        let prev = seq[i - 1:i].to_owned()
        if prev == "G" || prev == "A" || prev == "S" {
          score = score + 0.15
        }
      }
      let mut j = 2
      while j <= 4 && i + j < n {
        let nearby = seq[i + j:i + j + 1].to_owned()
        if nearby == "E" || nearby == "D" {
          score = score + 0.05
        }
        j = j + 1
      }
      if score > 1.0 {
        score = 1.0
      }
      sites.push({
        position: i,
        residue: "K",
        ptm_type: "acetylation",
        score,
        motif: window,
      })
    }
    i = i + 1
  }
  sites
}

///|
fn predict_ubiquitination(seq : String) -> Array[PtmSite] {
  let sites : Array[PtmSite] = []
  let n = seq.length()
  let mut i = 0
  while i < n {
    let aa = seq[i:i + 1].to_owned()
    if aa == "K" {
      let window = get_window(seq, i, 3)
      let mut score = 0.35
      let mut j = -2
      while j <= 2 && i + j >= 0 && i + j < n {
        if j != 0 {
          let nearby = seq[i + j:i + j + 1].to_owned()
          if nearby == "L" || nearby == "I" || nearby == "V" || nearby == "F" {
            score = score + 0.1
          }
        }
        j = j + 1
      }
      if score > 1.0 {
        score = 1.0
      }
      sites.push({
        position: i,
        residue: "K",
        ptm_type: "ubiquitination",
        score,
        motif: window,
      })
    }
    i = i + 1
  }
  sites
}

///|
fn predict_sumoylation(seq : String) -> Array[PtmSite] {
  let sites : Array[PtmSite] = []
  let n = seq.length()
  let mut i = 0
  while i < n {
    let aa = seq[i:i + 1].to_owned()
    if aa == "K" {
      let window = get_window(seq, i, 3)
      let mut score = 0.2
      let has_hydrophobic = if i - 1 >= 0 {
        let prev = seq[i - 1:i].to_owned()
        prev == "L" ||
        prev == "I" ||
        prev == "V" ||
        prev == "F" ||
        prev == "M" ||
        prev == "A"
      } else {
        false
      }
      let has_acidic = if i + 2 < n {
        let next2 = seq[i + 2:i + 3].to_owned()
        next2 == "E" || next2 == "D"
      } else {
        false
      }
      if has_hydrophobic {
        score = score + 0.3
      }
      if has_acidic {
        score = score + 0.3
      }
      if score > 1.0 {
        score = 1.0
      }
      if score > 0.5 {
        sites.push({
          position: i,
          residue: "K",
          ptm_type: "sumoylation",
          score,
          motif: window,
        })
      }
    }
    i = i + 1
  }
  sites
}

///|
fn predict_methylation(seq : String) -> Array[PtmSite] {
  let sites : Array[PtmSite] = []
  let n = seq.length()
  let mut i = 0
  while i < n {
    let aa = seq[i:i + 1].to_owned()
    if aa == "K" || aa == "R" {
      let window = get_window(seq, i, 3)
      let mut score = 0.35
      if aa == "R" {
        if i + 1 < n && seq[i + 1:i + 2].to_owned() == "G" {
          score = score + 0.25
          if i + 2 < n && seq[i + 2:i + 3].to_owned() == "R" {
            score = score + 0.2
          }
        }
      }
      if aa == "K" {
        if i + 1 < n {
          let next = seq[i + 1:i + 2].to_owned()
          if next == "A" || next == "S" {
            score = score + 0.15
          }
        }
      }
      if score > 1.0 {
        score = 1.0
      }
      sites.push({
        position: i,
        residue: aa,
        ptm_type: "methylation",
        score,
        motif: window,
      })
    }
    i = i + 1
  }
  sites
}

///|
fn predict_s_nitrosylation(seq : String) -> Array[PtmSite] {
  let sites : Array[PtmSite] = []
  let n = seq.length()
  let mut i = 0
  while i < n {
    let aa = seq[i:i + 1].to_owned()
    if aa == "C" {
      let window = get_window(seq, i, 4)
      let mut score = 0.3
      let has_acidic = if i - 1 >= 0 {
        let prev = seq[i - 1:i].to_owned()
        prev == "D" || prev == "E"
      } else {
        false
      }
      let has_basic = if i + 1 < n {
        let next = seq[i + 1:i + 2].to_owned()
        next == "R" || next == "K" || next == "H"
      } else {
        false
      }
      if has_acidic {
        score = score + 0.25
      }
      if has_basic {
        score = score + 0.25
      }
      if score > 1.0 {
        score = 1.0
      }
      if score > 0.5 {
        sites.push({
          position: i,
          residue: "C",
          ptm_type: "s_nitrosylation",
          score,
          motif: window,
        })
      }
    }
    i = i + 1
  }
  sites
}

///|
fn predict_disulfide_bonds(seq : String) -> Array[PtmSite] {
  let sites : Array[PtmSite] = []
  let n = seq.length()
  let mut i = 0
  while i < n {
    let aa = seq[i:i + 1].to_owned()
    if aa == "C" {
      let window = get_window(seq, i, 3)
      sites.push({
        position: i,
        residue: "C",
        ptm_type: "disulfide_bond",
        score: 0.7,
        motif: window,
      })
    }
    i = i + 1
  }
  sites
}

///|
fn push_all(dst : Array[PtmSite], src : Array[PtmSite]) -> Unit {
  let mut i = 0
  while i < src.length() {
    dst.push(src[i])
    i = i + 1
  }
}

///|
/// Predict all PTM sites in a protein sequence.
pub fn predict_all_ptm_sites(seq : String) -> PtmResult {
  let upper_seq = seq.to_upper()
  let all_sites : Array[PtmSite] = []
  push_all(all_sites, predict_phosphorylation(upper_seq))
  push_all(all_sites, predict_n_glycosylation(upper_seq))
  push_all(all_sites, predict_o_glycosylation(upper_seq))
  push_all(all_sites, predict_acetylation(upper_seq))
  push_all(all_sites, predict_ubiquitination(upper_seq))
  push_all(all_sites, predict_sumoylation(upper_seq))
  push_all(all_sites, predict_methylation(upper_seq))
  push_all(all_sites, predict_s_nitrosylation(upper_seq))
  push_all(all_sites, predict_disulfide_bonds(upper_seq))
  { sequence: seq, sites: all_sites, site_count: all_sites.length(), }
}

///|
/// Predict specific PTM type only.
/// ptm_type: "phosphorylation", "n_glycosylation", "o_glycosylation",
///           "acetylation", "ubiquitination", "sumoylation",
///           "methylation", "s_nitrosylation", "disulfide_bond"
pub fn predict_ptm_sites(seq : String, ptm_type : String) -> Array[PtmSite] {
  let upper_seq = seq.to_upper()
  match ptm_type {
    "phosphorylation" => predict_phosphorylation(upper_seq)
    "n_glycosylation" => predict_n_glycosylation(upper_seq)
    "o_glycosylation" => predict_o_glycosylation(upper_seq)
    "acetylation" => predict_acetylation(upper_seq)
    "ubiquitination" => predict_ubiquitination(upper_seq)
    "sumoylation" => predict_sumoylation(upper_seq)
    "methylation" => predict_methylation(upper_seq)
    "s_nitrosylation" => predict_s_nitrosylation(upper_seq)
    "disulfide_bond" => predict_disulfide_bonds(upper_seq)
    _ => []
  }
}

///|
/// Count PTM sites by type.
pub fn count_ptm_by_type(sites : Array[PtmSite]) -> Map[String, Int32] {
  let counts : Map[String, Int32] = Map([])
  let mut i = 0
  while i < sites.length() {
    let t = sites[i].ptm_type
    match counts.get(t) {
      Some(c) => counts[t] = c + 1
      None => counts[t] = 1
    }
    i = i + 1
  }
  counts
}

///|
/// Filter PTM sites by minimum score threshold.
pub fn filter_ptm_by_score(
  sites : Array[PtmSite],
  min_score : Double,
) -> Array[PtmSite] {
  let result : Array[PtmSite] = []
  let mut i = 0
  while i < sites.length() {
    if sites[i].score >= min_score {
      result.push(sites[i])
    }
    i = i + 1
  }
  result
}
