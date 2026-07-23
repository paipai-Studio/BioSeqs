///|
/// Prosite - Protein motif database search and analysis.
/// Inspired by Biopython Bio.Prosite.
/// Provides Prosite pattern matching and motif analysis.

///|
pub struct PrositePattern {
  accession : String
  name : String
  pattern : String
  description : String
  profile : String
  comments : String
}

///|
pub fn PrositePattern::new(accession : String, name : String, pattern : String) -> PrositePattern {
  PrositePattern::{
    accession,
    name,
    pattern,
    description: "",
    profile: "",
    comments: ""
  }
}

///|
pub struct PrositeMatch {
  pattern_accession : String
  pattern_name : String
  start : Int
  end : Int
  matched_sequence : String
  score : Double
}

///|
pub fn PrositeMatch::new(pattern_accession : String, pattern_name : String, start : Int, end : Int, matched_sequence : String) -> PrositeMatch {
  PrositeMatch::{
    pattern_accession,
    pattern_name,
    start,
    end,
    matched_sequence,
    score: 0.0
  }
}

///|
pub struct PrositeEntry {
  accession : String
  name : String
  description : String
  pattern : String
  profile : String
  type_ : String
  taxonomic_range : String
  created : String
  modified : String
  entries_with_motif : Int
  literature_references : Array[String]
}

///|
pub fn PrositeEntry::new(accession : String, name : String) -> PrositeEntry {
  PrositeEntry::{
    accession,
    name,
    description: "",
    pattern: "",
    profile: "",
    type_: "",
    taxonomic_range: "",
    created: "",
    modified: "",
    entries_with_motif: 0,
    literature_references: Array::new()
  }
}

///|
pub fn prosite_pattern_to_regex(pattern : String) -> String {
  let mut regex = ""
  let mut i = 0
  
  while i < pattern.length() {
    let c = pattern.unsafe_get(i).to_int()
    
    if c == '['.to_int() {
      regex = regex + "["
      i = i + 1
      while i < pattern.length() && pattern.unsafe_get(i).to_int() != ']'.to_int() {
        let pc = pattern.unsafe_get(i).to_int()
        if pc == '-'.to_int() {
          regex = regex + "-"
        } else {
          let cc = char_at(pattern, i)
          regex = regex + cc
        }
        i = i + 1
      }
      if i < pattern.length() && pattern.unsafe_get(i).to_int() == ']'.to_int() {
        regex = regex + "]"
        i = i + 1
      }
    } else if c == '{'.to_int() {
      regex = regex + "[^"
      i = i + 1
      while i < pattern.length() && pattern.unsafe_get(i).to_int() != '}'.to_int() {
        let cc = char_at(pattern, i)
        regex = regex + cc
        i = i + 1
      }
      if i < pattern.length() && pattern.unsafe_get(i).to_int() == '}'.to_int() {
        regex = regex + "]"
        i = i + 1
      }
    } else if c == '('.to_int() {
      let mut count = 0
      i = i + 1
      while i < pattern.length() && pattern.unsafe_get(i).to_int() != ')'.to_int() {
        let pc = pattern.unsafe_get(i).to_int()
        if pc >= '0'.to_int() && pc <= '9'.to_int() {
          count = count * 10 + (pc - '0'.to_int())
        }
        i = i + 1
      }
      if i < pattern.length() && pattern.unsafe_get(i).to_int() == ')'.to_int() {
        i = i + 1
      }
      if count > 0 {
        regex = regex + "{" + count.to_string() + "}"
      }
    } else if c == '<'.to_int() {
      regex = regex + "^"
      i = i + 1
    } else if c == '>'.to_int() {
      regex = regex + "$"
      i = i + 1
    } else if c == '.'.to_int() {
      regex = regex + "."
      i = i + 1
    } else {
      let cc = char_at(pattern, i)
      regex = regex + cc
      i = i + 1
    }
  }
  
  regex
}

///|
pub fn prosite_search(pattern : String, sequence : String) -> Array[PrositeMatch] {
  let matches : Array[PrositeMatch] = Array::new()
  
  if pattern.length() == 0 || sequence.length() == 0 {
    return matches
  }
  
  let normalized_seq = to_lowercase(sequence)
  let normalized_pattern = to_lowercase(pattern)
  let mut start_pos = 0
  
  while start_pos <= normalized_seq.length() {
    let mut pattern_idx = 0
    let mut seq_idx = start_pos
    let mut matched = true
    
    while seq_idx < normalized_seq.length() && pattern_idx < pattern.length() && matched {
      let p_char = normalized_pattern.unsafe_get(pattern_idx).to_int()
      
      if p_char == '['.to_int() {
        let mut in_bracket = true
        let mut bracket_content = ""
        pattern_idx = pattern_idx + 1
        
        while pattern_idx < normalized_pattern.length() && in_bracket {
          let bc = normalized_pattern.unsafe_get(pattern_idx).to_int()
          if bc == ']'.to_int() {
            in_bracket = false
          } else {
            bracket_content = bracket_content + char_at(normalized_pattern, pattern_idx)
          }
          pattern_idx = pattern_idx + 1
        }
        
        let s_char = char_at(normalized_seq, seq_idx)
        if contains(bracket_content, s_char) {
          seq_idx = seq_idx + 1
        } else {
          matched = false
        }
      } else if p_char == '{'.to_int() {
        let mut in_brace = true
        let mut brace_content = ""
        pattern_idx = pattern_idx + 1
        
        while pattern_idx < normalized_pattern.length() && in_brace {
          let bc = normalized_pattern.unsafe_get(pattern_idx).to_int()
          if bc == '}'.to_int() {
            in_brace = false
          } else {
            brace_content = brace_content + char_at(normalized_pattern, pattern_idx)
          }
          pattern_idx = pattern_idx + 1
        }
        
        let s_char = char_at(normalized_seq, seq_idx)
        if !contains(brace_content, s_char) {
          seq_idx = seq_idx + 1
        } else {
          matched = false
        }
      } else if p_char == '('.to_int() {
        let mut count_str = ""
        pattern_idx = pattern_idx + 1
        
        while pattern_idx < normalized_pattern.length() && normalized_pattern.unsafe_get(pattern_idx).to_int() != ')'.to_int() {
          count_str = count_str + char_at(normalized_pattern, pattern_idx)
          pattern_idx = pattern_idx + 1
        }
        if pattern_idx < normalized_pattern.length() && normalized_pattern.unsafe_get(pattern_idx).to_int() == ')'.to_int() {
          pattern_idx = pattern_idx + 1
        }
        
        let count = parse_int(count_str)
        seq_idx = seq_idx + count
      } else if p_char == '.'.to_int() || p_char == 'x'.to_int() {
        seq_idx = seq_idx + 1
        pattern_idx = pattern_idx + 1
      } else if p_char == '-'.to_int() {
        pattern_idx = pattern_idx + 1
      } else {
        let pc = char_at(normalized_pattern, pattern_idx)
        let sc = char_at(normalized_seq, seq_idx)
        if pc == sc {
          seq_idx = seq_idx + 1
          pattern_idx = pattern_idx + 1
        } else {
          matched = false
        }
      }
    }
    
    if matched && pattern_idx == pattern.length() {
      let matched_seq = substring(normalized_seq, start_pos, seq_idx - start_pos)
      let ps_match = PrositeMatch::new("", "Custom pattern", start_pos + 1, seq_idx, matched_seq)
      matches.push(ps_match)
    }
    
    start_pos = start_pos + 1
  }
  
  matches
}

///|
pub fn prosite_scan(sequence : String, patterns : Array[PrositePattern]) -> Array[PrositeMatch] {
  let all_matches : Array[PrositeMatch] = Array::new()
  
  for pattern in patterns {
    let matches = prosite_search(pattern.pattern, sequence)
    let mut i = 0
    while i < matches.length() {
      let m = matches[i]
      all_matches.push(PrositeMatch::new(
        pattern.accession,
        pattern.name,
        m.start,
        m.end,
        m.matched_sequence
      ))
      i = i + 1
    }
  }
  
  all_matches
}

///|
pub fn prosite_get_pattern(accession : String) -> PrositePattern? {
  let patterns = Map([
    ("PS00001", ("ATP/GTP-binding site motif", "[AG]-X(4)-G-K-[ST]")),
    ("PS00002", ("cAMP- and cGMP-dependent protein kinase phosphorylation site", "[RK]-X-[ST]")),
    ("PS00003", ("Protein kinase C phosphorylation site", "[ST]-X-[RK]")),
    ("PS00004", ("Casein kinase II phosphorylation site", "[ST]-X(2)-[DE]")),
    ("PS00005", ("Tyrosine kinase phosphorylation site", "[RK]-X(2)-[DE]-X(3)-Y")),
    ("PS00006", ("N-myristoylation site", "G-X(2)-[AGCS]-X(2)-[STAGC]-[LIVMFYW]")),
    ("PS00007", ("CAMP-dependent protein kinase phosphorylation site", "[RK]-X(2)-[ST]")),
    ("PS00008", ("Amidation site", "[AG]-G-[RK]-R"))
  ], capacity=8)
  
  if patterns.contains(accession) {
    let info = patterns.get(accession).unwrap()
    Some(PrositePattern::new(accession, info.0, info.1))
  } else {
    None
  }
}

///|
pub fn prosite_calculate_score(match_ : PrositeMatch, pattern : PrositePattern) -> Double {
  let mut score = 100.0
  
  if match_.matched_sequence.length() < pattern.pattern.length() {
    score = score * 0.5
  }
  
  let mut mismatches = 0
  let mut i = 0
  while i < match_.matched_sequence.length() && i < pattern.pattern.length() {
    let pc = char_at(pattern.pattern, i)
    let mc = char_at(match_.matched_sequence, i)
    
    if pc != mc && pattern.pattern.unsafe_get(i).to_int() != '.'.to_int() {
      mismatches = mismatches + 1
    }
    i = i + 1
  }
  
  score = score - mismatches.to_double() * 10.0
  
  if score < 0.0 { 0.0 } else { score }
}

///|
pub fn prosite_create_example_patterns() -> Array[PrositePattern] {
  let patterns : Array[PrositePattern] = Array::new()
  
  let ps1 = PrositePattern::new("PS00001", "ATP/GTP-binding site motif", "[AG]-X(4)-G-K-[ST]")
  patterns.push(ps1)
  
  let ps2 = PrositePattern::new("PS00003", "Protein kinase C phosphorylation site", "[ST]-X-[RK]")
  patterns.push(ps2)
  
  let ps3 = PrositePattern::new("PS00006", "N-myristoylation site", "G-X(2)-[AGCS]-X(2)-[STAGC]-[LIVMFYW]")
  patterns.push(ps3)
  
  patterns
}

///|
pub fn prosite_create_example_entry() -> PrositeEntry {
  let entry = PrositeEntry::new("PS00001", "ATP/GTP-binding site motif")
  PrositeEntry::{
    accession: entry.accession,
    name: entry.name,
    description: "ATP/GTP-binding proteins typically possess a phosphate-binding loop (P-loop) which is involved in binding the phosphate groups of the nucleotide.",
    pattern: "[AG]-X(4)-G-K-[ST]",
    profile: "",
    type_: "PATTERN",
    taxonomic_range: "Eukaryota; Bacteria; Archaea",
    created: "15-NOV-1990",
    modified: "20-FEB-2002",
    entries_with_motif: 1234,
    literature_references: ["PMID:12345678", "PMID:87654321"]
  }
}