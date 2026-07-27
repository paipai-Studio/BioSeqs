///|
/// Bio.Align.Mauve - Mauve genome alignment format parser.
/// Mauve is a multiple genome alignment algorithm that identifies
/// conserved segments and inversions. This module supports:
/// - Parsing Mauve alignment output (XMFA/MAF format with Mauve extensions)
/// - Locally Collinear Blocks (LCBs)
/// - Inversion detection
/// - Progressive alignment support

// ============================================================================
// Data Structures
// ============================================================================

///|
/// Represents a Locally Collinear Block (LCB) in a Mauve alignment.
/// LCBs are conserved regions between genomes that appear
/// in the same order without rearrangements.
pub struct MauveLCB {
  /// Block identifier
  id : String
  /// Block score
  mut score : Double
  /// Sequence entries in this block
  sequences : Array[MauveSequence]
  /// Whether this block is collinear (no rearrangements)
  is_collinear : Bool
} derive(Eq)

///|
/// Represents a sequence entry in a Mauve alignment block.
pub struct MauveSequence {
  /// Sequence name (genome.chromosome)
  name : String
  /// Start coordinate
  start : Int
  /// End coordinate
  end : Int
  /// Length of the aligned region
  length : Int
  /// Strand/orientation
  strand : String
  /// Total sequence size
  seq_size : Int
  /// Whether this sequence segment is inverted
  is_inverted : Bool
  /// Group number (for progressive alignment)
  group : Int
  /// Parent block ID (for progressive alignment)
  parent_id : String
} derive(Eq)

///|
/// Represents a breakpoint between LCBs.
pub struct MauveBreakpoint {
  /// The first LCB ID
  lcb1_id : String
  /// The second LCB ID
  lcb2_id : String
  /// Sequence name where breakpoint occurs
  seq_name : String
  /// Position in the sequence
  position : Int
  /// Breakpoint type (inversion, translocation, etc.)
  breakpoint_type : String
} derive(Eq)

///|
/// Represents an inversion in a Mauve alignment.
pub struct MauveInversion {
  /// Sequence name
  seq_name : String
  /// Start position of inversion
  start : Int
  /// End position of inversion
  end : Int
  /// Length of inversion
  length : Int
  /// Affected LCB IDs
  affected_lcbs : Array[String]
} derive(Eq)

///|
/// Configuration for Mauve alignment parsing and analysis.
pub struct MauveConfig {
  /// Minimum LCB score to report
  min_lcb_score : Double
  /// Minimum block length
  min_block_length : Int
  /// Whether to detect inversions
  detect_inversions : Bool
  /// Whether to detect breakpoints
  detect_breakpoints : Bool
} derive(Eq)

///|
/// Default Mauve configuration.
pub fn MauveConfig::default() -> MauveConfig {
  MauveConfig::{
    min_lcb_score: 0.0,
    min_block_length: 50,
    detect_inversions: true,
    detect_breakpoints: true,
  }
}

///|
/// Represents a complete Mauve alignment.
pub struct MauveAlignment {
  /// Version string
  version : String
  /// All LCBs found in the alignment
  lcbs : Array[MauveLCB]
  /// All sequences across all LCBs
  sequences : Array[MauveSequence]
  /// Detected inversions
  inversions : Array[MauveInversion]
  /// Detected breakpoints
  breakpoints : Array[MauveBreakpoint]
  /// Progressive alignment groups
  groups : Map[Int, Array[String]]
} derive(Eq)

// ============================================================================
// MauveSequence Implementation
// ============================================================================

///|
/// Create a new MauveSequence.
pub fn new_mauve_sequence(
  name : String,
  start : Int,
  end : Int,
  strand : String,
  seq_size : Int,
) -> MauveSequence {
  let length = if end > start { end - start } else { start - end }
  MauveSequence::{
    name,
    start,
    end,
    length,
    strand,
    seq_size,
    is_inverted: strand == "-",
    group: 0,
    parent_id: "",
  }
}

///|
/// Get the coverage of this sequence segment.
pub fn mauve_seq_coverage(seq : MauveSequence) -> Double {
  if seq.seq_size == 0 {
    return 0.0
  }
  (seq.length.to_double() / seq.seq_size.to_double()) * 100.0
}

// ============================================================================
// MauveLCB Implementation
// ============================================================================

///|
/// Create a new empty MauveLCB.
pub fn new_mauve_lcb(id : String) -> MauveLCB {
  MauveLCB::{
    id,
    score: 0.0,
    sequences: Array::new(),
    is_collinear: true,
  }
}

///|
/// Set the score of an LCB.
pub fn mauve_lcb_set_score(lcb : MauveLCB, score : Double) -> Unit {
  lcb.score = score
}

///|
/// Add a sequence to the LCB.
pub fn mauve_lcb_add_seq(lcb : MauveLCB, seq : MauveSequence) -> Unit {
  lcb.sequences.push(seq)
}

///|
/// Get the number of sequences in the LCB.
pub fn mauve_lcb_num_seqs(lcb : MauveLCB) -> Int {
  lcb.sequences.length()
}

///|
/// Check if all sequences in the LCB are in the same orientation.
pub fn mauve_lcb_is_consistent(lcb : MauveLCB) -> Bool {
  if lcb.sequences.length() <= 1 {
    return true
  }
  let first_strand = lcb.sequences[0].strand
  for s = 0; s < lcb.sequences.length(); s = s + 1 {
    if lcb.sequences[s].strand != first_strand {
      return false
    }
  }
  true
}

///|
/// Get the sequence names in this LCB.
pub fn mauve_lcb_seq_names(lcb : MauveLCB) -> Array[String] {
  let names = Array::new()
  for s = 0; s < lcb.sequences.length(); s = s + 1 {
    names.push(lcb.sequences[s].name)
  }
  names
}

// ============================================================================
// MauveAlignment Implementation
// ============================================================================

///|
/// Create a new empty MauveAlignment.
pub fn new_mauve_alignment() -> MauveAlignment {
  MauveAlignment::{
    version: "1.0",
    lcbs: Array::new(),
    sequences: Array::new(),
    inversions: Array::new(),
    breakpoints: Array::new(),
    groups: Map([], capacity=10),
  }
}

///|
/// Add an LCB to the alignment.
pub fn mauve_add_lcb(ali : MauveAlignment, lcb : MauveLCB) -> Unit {
  ali.lcbs.push(lcb)
}

///|
/// Get the number of LCBs.
pub fn mauve_num_lcbs(ali : MauveAlignment) -> Int {
  ali.lcbs.length()
}

///|
/// Get all sequence names across all LCBs.
pub fn mauve_all_seq_names(ali : MauveAlignment) -> Array[String] {
  let names = Array::new()
  for l = 0; l < ali.lcbs.length(); l = l + 1 {
    let lcb = ali.lcbs[l]
    for s = 0; s < lcb.sequences.length(); s = s + 1 {
      let name = lcb.sequences[s].name
      if !names.contains(name) {
        names.push(name)
      }
    }
  }
  names
}

///|
/// Get sequences for a specific genome.
pub fn mauve_get_genome_seqs(ali : MauveAlignment, genome_name : String) -> Array[MauveSequence] {
  let seqs = Array::new()
  let prefix_len = genome_name.length()
  for l = 0; l < ali.lcbs.length(); l = l + 1 {
    let lcb = ali.lcbs[l]
    for s = 0; s < lcb.sequences.length(); s = s + 1 {
      let seq = lcb.sequences[s]
      let matches = if seq.name.length() >= prefix_len {
        starts_with(seq.name, 0, genome_name)
      } else {
        false
      }
      if matches || seq.name == genome_name {
        seqs.push(seq)
      }
    }
  }
  seqs
}

// ============================================================================
// Mauve Parser (XMFA/MAF with Mauve extensions)
// ============================================================================

///|
/// Parse a Mauve alignment from MAF/XMFA format.
pub fn parse_mauve(content : String) -> MauveAlignment {
  let mut result = new_mauve_alignment()
  let lines = content.split("\n").to_array()
  let mut current_lcb : MauveLCB? = None
  let mut lcb_counter = 0
  
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    let trimmed = trim(line)
    
    if trimmed.length() == 0 {
      match current_lcb {
        Some(lcb) => {
          result.lcbs.push(lcb)
          current_lcb = None
          ignore(())
        }
        None => ignore(())
      }
      continue
    }
    
    if starts_with(trimmed, 0, "##maf") || starts_with(trimmed, 0, "##mauve") {
      if contains(trimmed, "version=") {
        let parts = split_by_char(trimmed, ' '.to_int())
        for p = 0; p < parts.length(); p = p + 1 {
          let part = parts[p]
          if starts_with(part, 0, "version=") {
            result = MauveAlignment::{ ..result, version: part[8:part.length()].to_owned() }
          }
        }
      }
      continue
    }
    
    if starts_with(trimmed, 0, "a ") {
      lcb_counter = lcb_counter + 1
      let lcb_id = "" + "lcb_" + lcb_counter.to_string()
      current_lcb = Some(new_mauve_lcb(lcb_id))
      
      let parts = split_by_char(trimmed, ' '.to_int())
      for p = 0; p < parts.length(); p = p + 1 {
        let part = parts[p]
        if starts_with(part, 0, "score=") {
          let score_val = parse_double(part[6:part.length()].to_owned())
          match current_lcb {
            Some(lcb) => {
              let score = if score_val is Some(_) { score_val.unwrap() } else { 0.0 }
              let new_lcb = MauveLCB::{ ..lcb, score: score }
              current_lcb = Some(new_lcb)
              ignore(())
            }
            None => ignore(())
          }
        } else if starts_with(part, 0, "collinear=") {
          let val = part[10:part.length()].to_owned()
          match current_lcb {
            Some(lcb) => {
              let new_lcb = MauveLCB::{ ..lcb, is_collinear: val == "yes" || val == "1" }
              current_lcb = Some(new_lcb)
              ignore(())
            }
            None => ignore(())
          }
        }
      }
      continue
    }
    
    if starts_with(trimmed, 0, "s ") || starts_with(trimmed, 0, "s\t") {
      if current_lcb is None {
        lcb_counter = lcb_counter + 1
        current_lcb = Some(new_mauve_lcb("lcb_" + lcb_counter.to_string()))
      }
      
      let mut parts = split_by_tab(trimmed)
      if parts.length() < 6 {
        parts = split_by_char(trimmed, ' '.to_int())
      }
      if parts.length() >= 6 {
        let name = parts[1]
        let start = parse_int(parts[2])
        let length = parse_int(parts[3])
        let strand = parts[4]
        let seq_size = parse_int(parts[5])
        let end = start + length
        
        let seq = new_mauve_sequence(name, start, end, strand, seq_size)
        match current_lcb {
          Some(lcb) => {
            lcb.sequences.push(seq)
            ignore(())
          }
          None => ignore(())
        }
        result.sequences.push(seq)
      }
      continue
    }
    
    if starts_with(trimmed, 0, "g ") {
      let parts = split_by_char(trimmed, ' '.to_int())
      if parts.length() >= 3 {
        let group_num = parse_int(parts[1])
        let seq_name = parts[2]
        let existing = result.groups.get(group_num)
        let group_arr = if existing is Some(_) {
          existing.unwrap()
        } else {
          let arr = Array::new()
          result.groups[group_num] = arr
          arr
        }
        group_arr.push(seq_name)
      }
      continue
    }
  }
  
  match current_lcb {
    Some(lcb) => {
      result.lcbs.push(lcb)
      ignore(())
    }
    None => ignore(())
  }
  
  result = detect_mauve_inversions(result)
  result = detect_mauve_breakpoints(result)
  
  result
}

// ============================================================================
// Inversion Detection
// ============================================================================

///|
/// Detect inversions in the alignment.
pub fn detect_mauve_inversions(ali : MauveAlignment) -> MauveAlignment {
  let result = ali
  let seq_names = mauve_all_seq_names(result)
  
  for n = 0; n < seq_names.length(); n = n + 1 {
    let name = seq_names[n]
    let seqs = mauve_get_genome_seqs(result, name)
    if seqs.length() < 2 {
      continue
    }
    
    for i = 1; i < seqs.length(); i = i + 1 {
      let prev = seqs[i - 1]
      let curr = seqs[i]
      
      if prev.strand != curr.strand {
        let inv_start = if prev.start < curr.start { prev.start } else { curr.start }
        let inv_end = if prev.end > curr.end { prev.end } else { curr.end }
        let inv = MauveInversion::{
          seq_name: name,
          start: inv_start,
          end: inv_end,
          length: if curr.end - prev.start > 0 { curr.end - prev.start } else { prev.start - curr.end },
          affected_lcbs: Array::new(),
        }
        
        for l = 0; l < result.lcbs.length(); l = l + 1 {
          let lcb = result.lcbs[l]
          for s = 0; s < lcb.sequences.length(); s = s + 1 {
            let seq = lcb.sequences[s]
            if seq.name == name && 
               ((seq.start >= inv.start && seq.end <= inv.end) ||
                (seq.end >= inv.start && seq.start <= inv.end)) {
              inv.affected_lcbs.push(lcb.id)
            }
          }
        }
        
        result.inversions.push(inv)
      }
    }
  }
  result
}

///|
/// Detect breakpoints between LCBs.
pub fn detect_mauve_breakpoints(ali : MauveAlignment) -> MauveAlignment {
  let result = ali
  if result.lcbs.length() < 2 {
    return result
  }
  
  for i = 1; i < result.lcbs.length(); i = i + 1 {
    let prev_lcb = result.lcbs[i - 1]
    let curr_lcb = result.lcbs[i]
    
    for ps = 0; ps < prev_lcb.sequences.length(); ps = ps + 1 {
      let prev_seq = prev_lcb.sequences[ps]
      for cs = 0; cs < curr_lcb.sequences.length(); cs = cs + 1 {
        let curr_seq = curr_lcb.sequences[cs]
        if prev_seq.name == curr_seq.name {
          let gap_start = if prev_seq.end > curr_seq.start { prev_seq.end } else { curr_seq.start }
          let gap_end = if prev_seq.start < curr_seq.end { prev_seq.start } else { curr_seq.end }
          
          if gap_start > gap_end {
            let bp_type = if prev_seq.strand != curr_seq.strand { "inversion" } else { "translocation" }
            let bp = MauveBreakpoint::{
              lcb1_id: prev_lcb.id,
              lcb2_id: curr_lcb.id,
              seq_name: prev_seq.name,
              position: gap_start,
              breakpoint_type: bp_type,
            }
            result.breakpoints.push(bp)
          }
        }
      }
    }
  }
  result
}

// ============================================================================
// Mauve Analysis Functions
// ============================================================================

///|
/// Create a new inversion.
pub fn new_mauve_inversion(
  seq_name : String,
  start : Int,
  end : Int,
  length : Int,
) -> MauveInversion {
  MauveInversion::{
    seq_name,
    start,
    end,
    length,
    affected_lcbs: Array::new(),
  }
}

///|
/// Calculate genome coverage for each genome.
pub fn mauve_genome_coverage(ali : MauveAlignment) -> Map[String, Double] {
  let coverage = Map([], capacity=10)
  let seq_names = mauve_all_seq_names(ali)
  
  for n = 0; n < seq_names.length(); n = n + 1 {
    let name = seq_names[n]
    let seqs = mauve_get_genome_seqs(ali, name)
    let mut total_aligned = 0.0
    let mut max_end = 0
    let mut min_start = 99999999
    
    for s = 0; s < seqs.length(); s = s + 1 {
      let seq = seqs[s]
      total_aligned = total_aligned + seq.length.to_double()
      if seq.end > max_end { max_end = seq.end }
      if seq.start < min_start { min_start = seq.start }
    }
    
    let total_range = (max_end - min_start).to_double()
    if total_range > 0.0 {
      coverage[name] = (total_aligned / total_range) * 100.0
    } else {
      coverage[name] = 0.0
    }
  }
  
  coverage
}

///|
/// Count the number of conserved segments per genome.
pub fn mauve_conserved_segments(ali : MauveAlignment) -> Map[String, Int] {
  let segments = Map([], capacity=10)
  let seq_names = mauve_all_seq_names(ali)
  
  for n = 0; n < seq_names.length(); n = n + 1 {
    let name = seq_names[n]
    let count = mauve_get_genome_seqs(ali, name).length()
    segments[name] = count
  }
  
  segments
}

///|
/// Calculate the rearrangement rate (breakpoints per Mb).
pub fn mauve_rearrangement_rate(ali : MauveAlignment) -> Double {
  let total_size = ali.sequences.length()
  if total_size == 0 {
    return 0.0
  }
  (ali.breakpoints.length().to_double() / total_size.to_double()) * 1000000.0
}

///|
/// Get syntenic regions (consecutive collinear LCBs).
pub fn mauve_syntenic_regions(ali : MauveAlignment) -> Array[Array[String]] {
  let regions = Array::new()
  let mut current_region = Array::new()
  let mut prev_collinear = true
  
  for l = 0; l < ali.lcbs.length(); l = l + 1 {
    let lcb = ali.lcbs[l]
    if lcb.is_collinear && prev_collinear {
      current_region.push(lcb.id)
    } else {
      if current_region.length() > 0 {
        regions.push(current_region)
      }
      current_region = Array::new()
      if lcb.is_collinear {
        current_region.push(lcb.id)
      }
    }
    prev_collinear = lcb.is_collinear
  }
  
  if current_region.length() > 0 {
    regions.push(current_region)
  }
  
  regions
}

///|
/// Filter LCBs by score threshold.
pub fn mauve_filter_lcbs(ali : MauveAlignment, min_score : Double) -> MauveAlignment {
  let mut result = MauveAlignment::{ ..ali, lcbs: Array::new(), inversions: Array::new(), breakpoints: Array::new() }
  
  for l = 0; l < ali.lcbs.length(); l = l + 1 {
    let lcb = ali.lcbs[l]
    if lcb.score >= min_score {
      result.lcbs.push(lcb)
    }
  }
  
  result = detect_mauve_inversions(result)
  result = detect_mauve_breakpoints(result)
  
  result
}

// ============================================================================
// Mauve Progressive Alignment
// ============================================================================

///|
/// Initialize progressive alignment.
pub struct MauveProgressiveAligner {
  /// Pending sequences to add
  pending_seqs : Array[String]
  /// Current alignment
  alignment : MauveAlignment
  /// Guide tree (order of addition)
  mut guide_tree : Array[Array[String]]
}

///|
/// Create a new progressive aligner.
pub fn new_mauve_progressive_aligner() -> MauveProgressiveAligner {
  MauveProgressiveAligner::{
    pending_seqs: Array::new(),
    alignment: new_mauve_alignment(),
    guide_tree: Array::new(),
  }
}

///|
/// Add a sequence to the progressive alignment.
pub fn mauve_add_sequence(aligner : MauveProgressiveAligner, seq_name : String) -> Unit {
  aligner.pending_seqs.push(seq_name)
}

///|
/// Build the guide tree for progressive alignment.
pub fn mauve_build_guide_tree(aligner : MauveProgressiveAligner) -> Unit {
  let groups = Array::new()
  for s = 0; s < aligner.pending_seqs.length(); s = s + 1 {
    let group = Array::new()
    group.push(aligner.pending_seqs[s])
    groups.push(group)
  }
  aligner.guide_tree = groups
}

///|
/// Execute progressive alignment by adding sequences one at a time.
pub fn mauve_execute_progressive(aligner : MauveProgressiveAligner) -> MauveAlignment {
  mauve_build_guide_tree(aligner)
  
  let result = aligner.alignment
  
  for i = 0; i < aligner.guide_tree.length(); i = i + 1 {
    result.groups[i] = aligner.guide_tree[i]
  }
  
  result
}

// ============================================================================
// Statistics Output
// ============================================================================

///|
/// Generate a summary report of the Mauve alignment.
pub fn mauve_summary(ali : MauveAlignment) -> String {
  let mut result = ""
  result = result + "Mauve Alignment Summary\n"
  result = result + "====================\n\n"
  
  result = result + "Version: "
  result = result + ali.version
  result = result + "\n\n"
  
  result = result + "Locally Collinear Blocks (LCBs): "
  result = result + ali.lcbs.length().to_string()
  result = result + "\n\n"
  
  result = result + "Genomes aligned:\n"
  let seq_names = mauve_all_seq_names(ali)
  let coverage = mauve_genome_coverage(ali)
  for n = 0; n < seq_names.length(); n = n + 1 {
    let name = seq_names[n]
    result = result + "  - "
    result = result + name
    result = result + " (coverage: "
    result = result + coverage[name].to_string()
    result = result + "%)\n"
  }
  
  result = result + "\n"
  result = result + "Inversions detected: "
  result = result + ali.inversions.length().to_string()
  result = result + "\n"
  
  for i = 0; i < ali.inversions.length(); i = i + 1 {
    let inv = ali.inversions[i]
    result = result + "  - "
    result = result + inv.seq_name
    result = result + ": "
    result = result + inv.start.to_string()
    result = result + "-"
    result = result + inv.end.to_string()
    result = result + " ("
    result = result + inv.length.to_string()
    result = result + " bp)\n"
  }
  
  result = result + "\n"
  result = result + "Breakpoints detected: "
  result = result + ali.breakpoints.length().to_string()
  result = result + "\n"
  
  let bp_types = Map([], capacity=5)
  for b = 0; b < ali.breakpoints.length(); b = b + 1 {
    let bp = ali.breakpoints[b]
    let existing = bp_types.get(bp.breakpoint_type)
    if existing is Some(_) {
      bp_types[bp.breakpoint_type] = existing.unwrap() + 1
    } else {
      bp_types[bp.breakpoint_type] = 1
    }
  }
  
  let bp_keys = bp_types.keys().collect()
  for k = 0; k < bp_keys.length(); k = k + 1 {
    let key = bp_keys[k]
    result = result + "  "
    result = result + key
    result = result + ": "
    result = result + bp_types[key].to_string()
    result = result + "\n"
  }
  
  result = result + "\n"
  result = result + "Syntenic regions:\n"
  let syntenic = mauve_syntenic_regions(ali)
  result = result + "  Number: "
  result = result + syntenic.length().to_string()
  result = result + "\n"
  
  let segments = mauve_conserved_segments(ali)
  result = result + "\nConserved segments per genome:\n"
  let seg_keys = segments.keys().collect()
  for k = 0; k < seg_keys.length(); k = k + 1 {
    let key = seg_keys[k]
    result = result + "  "
    result = result + key
    result = result + ": "
    result = result + segments[key].to_string()
    result = result + " segments\n"
  }
  
  result = result + "\nRearrangement rate: "
  result = result + mauve_rearrangement_rate(ali).to_string()
  result = result + " breakpoints/Mb\n"
  
  result
}

///|
/// Export alignment to BED format for visualization.
pub fn mauve_to_bed(ali : MauveAlignment) -> String {
  let mut result = ""
  result = result + "track name=\"MauveAlignment\" description=\"Mauve Alignment Blocks\"\n"
  
  for l = 0; l < ali.lcbs.length(); l = l + 1 {
    let lcb = ali.lcbs[l]
    for s = 0; s < lcb.sequences.length(); s = s + 1 {
      let seq = lcb.sequences[s]
      result = result + seq.name
      result = result + "\t"
      result = result + seq.start.to_string()
      result = result + "\t"
      result = result + seq.end.to_string()
      result = result + "\t"
      result = result + lcb.id
      result = result + "\t"
      result = result + lcb.score.to_string()
      result = result + "\t"
      result = result + seq.strand
      result = result + "\n"
    }
  }
  
  result
}

///|
/// Export inversions to BED format.
pub fn mauve_inversions_to_bed(ali : MauveAlignment) -> String {
  let mut result = ""
  result = result + "track name=\"Inversions\" color=\"100,0,0\"\n"
  
  for i = 0; i < ali.inversions.length(); i = i + 1 {
    let inv = ali.inversions[i]
    result = result + inv.seq_name
    result = result + "\t"
    result = result + inv.start.to_string()
    result = result + "\t"
    result = result + inv.end.to_string()
    result = result + "\tInversion\t"
    result = result + inv.length.to_string()
    result = result + "\t-\n"
  }
  
  result
}