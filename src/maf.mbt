///|
/// Bio.Align.MAF - Multiple Alignment Format parser.
/// MAF is a format for storing multiple genome alignments,
/// originally developed with the Mauve aligner.
/// It stores alignment blocks with sequence coordinates,
/// strand information, and conservation tracking.

// ============================================================================
// Data Structures
// ============================================================================

///|
/// Represents a single sequence entry in an alignment block.
pub struct MAFSequence {
  /// Sequence name (e.g., "hg19.chr1")
  seq_name : String
  /// Start position in the sequence
  start : Int
  /// Length of the alignment region
  length : Int
  /// Strand ('+' or '-')
  strand : String
  /// Total size of the sequence
  seq_size : Int
  /// The aligned sequence string (may contain dashes for gaps)
  aligned_seq : String
}

///|
/// Represents an event line (e) in a MAF block.
/// Events describe features like inversions and breakpoints.
pub struct MAFEvent {
  /// Event type (e.g., "i" for inversion, "b" for breakpoint)
  event_type : String
  /// Sequence name
  seq_name : String
  /// Position
  position : Int
  /// Additional information
  description : String
}

///|
/// Represents a comment/annotation line in a MAF block.
pub struct MAFComment {
  /// Comment key (text after '#')
  key : String
}

///|
/// Represents a single alignment block in a MAF file.
pub struct MAFBlock {
  /// Block score (e.g., from 'a' line)
  mut score : Double?
  /// Other attributes from 'a' line
  mut attributes : Map[String, String]
  /// Sequence entries
  sequences : Array[MAFSequence]
  /// Comment lines within the block
  mut comments : Array[MAFComment]
  /// Event lines
  mut events : Array[MAFEvent]
  /// Whether this block has been parsed
  mut is_complete : Bool
}

///|
/// Represents a complete MAF file with header and blocks.
pub struct MAFAliignment {
  /// Version string (e.g., "1")
  mut version : String
  /// Alignment score
  mut score : Double?
  /// Free comments at the top level
  mut comments : Array[String]
  /// The alignment blocks
  blocks : Array[MAFBlock]
}

///|
/// Parser configuration for MAF files.
pub struct MAFConfig {
  /// Whether to store raw sequences
  store_sequences : Bool
  /// Minimum block score to keep
  min_score : Double
}

///|
/// Default MAF configuration.
pub fn MAFConfig::default() -> MAFConfig {
  MAFConfig::{
    store_sequences: true,
    min_score: 0.0,
  }
}

// ============================================================================
// MAFSequence Implementation
// ============================================================================

///|
/// Create a new MAFSequence.
pub fn new_maf_sequence(
  seq_name : String,
  start : Int,
  length : Int,
  strand : String,
  seq_size : Int,
  aligned_seq : String,
) -> MAFSequence {
  MAFSequence::{
    seq_name,
    start,
    length,
    strand,
    seq_size,
    aligned_seq,
  }
}

///|
/// Get the end position (start + length).
pub fn maf_seq_end(seq : MAFSequence) -> Int {
  seq.start + seq.length
}

///|
/// Get the sequence (without gaps).
pub fn maf_seq_sequence(seq : MAFSequence) -> String {
  let mut result = ""
  for i = 0; i < seq.aligned_seq.length(); i = i + 1 {
    let c = seq.aligned_seq.unsafe_get(i)
    if c.to_int() != '-'.to_int() && c.to_int() != '.'.to_int() {
      result = result + seq.aligned_seq[i:i+1].to_owned()
    }
  }
  result
}

///|
/// Calculate the fraction of the sequence that aligns (no gaps).
pub fn maf_seq_coverage(seq : MAFSequence) -> Double {
  let aligned_len = seq.aligned_seq.length()
  if aligned_len == 0 {
    return 0.0
  }
  let matches = maf_seq_sequence(seq).length()
  matches.to_double() / aligned_len.to_double()
}

// ============================================================================
// MAFBlock Implementation
// ============================================================================

///|
/// Create a new empty MAFBlock.
pub fn new_maf_block() -> MAFBlock {
  MAFBlock::{
    score: None,
    attributes: Map([], capacity=10),
    sequences: Array::new(),
    comments: Array::new(),
    events: Array::new(),
    is_complete: false,
  }
}

///|
/// Add a sequence to the block.
pub fn maf_block_add_seq(block : MAFBlock, seq : MAFSequence) -> MAFBlock {
  block.sequences.push(seq)
  block
}

///|
/// Add a comment to the block.
pub fn maf_block_add_comment(block : MAFBlock, comment : MAFComment) -> MAFBlock {
  block.comments.push(comment)
  block
}

///|
/// Add an event to the block.
pub fn maf_block_add_event(block : MAFBlock, event : MAFEvent) -> MAFBlock {
  block.events.push(event)
  block
}

///|
/// Set the score of a block.
pub fn maf_block_set_score(block : MAFBlock, score : Double?) -> MAFBlock {
  block.score = score
  block
}

///|
/// Get the number of sequences in the block.
pub fn maf_block_num_seqs(block : MAFBlock) -> Int {
  block.sequences.length()
}

///|
/// Get the alignment length (number of columns).
pub fn maf_block_alignment_length(block : MAFBlock) -> Int {
  if block.sequences.length() == 0 {
    return 0
  }
  block.sequences[0].aligned_seq.length()
}

///|
/// Calculate the percent identity of the alignment block.
/// Uses all pairwise comparisons to compute identity.
pub fn maf_block_percent_identity(block : MAFBlock) -> Double {
  let n_seqs = block.sequences.length()
  if n_seqs < 2 {
    return 100.0
  }
  let aln_len = maf_block_alignment_length(block)
  if aln_len == 0 {
    return 0.0
  }
  
  let mut total_matches = 0
  let mut total_pairs = 0
  
  for col = 0; col < aln_len; col = col + 1 {
    let bases = Array::new()
    for s = 0; s < block.sequences.length(); s = s + 1 {
      let seq = block.sequences[s]
      let c = seq.aligned_seq.unsafe_get(col)
      if c.to_int() != '-'.to_int() && c.to_int() != '.'.to_int() {
        bases.push(seq.aligned_seq[col:col+1].to_owned())
      }
    }
    let n_bases = bases.length()
    if n_bases >= 2 {
      total_pairs = total_pairs + n_bases * (n_bases - 1) / 2
      for i = 0; i < n_bases; i = i + 1 {
        for j = i + 1; j < n_bases; j = j + 1 {
          if bases[i] == bases[j] {
            total_matches = total_matches + 1
          }
        }
      }
    }
  }
  
  if total_pairs == 0 {
    return 0.0
  }
  (total_matches.to_double() / total_pairs.to_double()) * 100.0
}

///|
/// Calculate block coverage (fraction of positions aligned).
pub fn maf_block_coverage(block : MAFBlock) -> Double {
  let n_seqs = block.sequences.length()
  if n_seqs == 0 {
    return 0.0
  }
  let mut total_coverage = 0.0
  for s = 0; s < block.sequences.length(); s = s + 1 {
    total_coverage = total_coverage + maf_seq_coverage(block.sequences[s])
  }
  total_coverage / n_seqs.to_double()
}

// ============================================================================
// MAFAliignment Implementation
// ============================================================================

///|
/// Create a new empty MAFAliignment.
pub fn new_maf_alignment() -> MAFAliignment {
  MAFAliignment::{
    version: "1",
    score: None,
    comments: Array::new(),
    blocks: Array::new(),
  }
}

///|
/// Add a block to the alignment.
pub fn maf_add_block(ali : MAFAliignment, block : MAFBlock) -> MAFAliignment {
  ali.blocks.push(block)
  ali
}

///|
/// Get the number of blocks.
pub fn maf_num_blocks(ali : MAFAliignment) -> Int {
  ali.blocks.length()
}

///|
/// Get total alignment length across all blocks.
pub fn maf_total_length(ali : MAFAliignment) -> Int {
  let mut total = 0
  for b = 0; b < ali.blocks.length(); b = b + 1 {
    total = total + maf_block_alignment_length(ali.blocks[b])
  }
  total
}

///|
/// Get all sequence names across all blocks.
pub fn maf_all_seq_names(ali : MAFAliignment) -> Array[String] {
  let names = Array::new()
  for b = 0; b < ali.blocks.length(); b = b + 1 {
    let block = ali.blocks[b]
    for s = 0; s < block.sequences.length(); s = s + 1 {
      let seq = block.sequences[s]
      let found = names.contains(seq.seq_name)
      if !found {
        names.push(seq.seq_name)
      }
    }
  }
  names
}

///|
/// Get the average percent identity across all blocks.
pub fn maf_average_percent_identity(ali : MAFAliignment) -> Double {
  if ali.blocks.length() == 0 {
    return 0.0
  }
  let mut total = 0.0
  for b = 0; b < ali.blocks.length(); b = b + 1 {
    total = total + maf_block_percent_identity(ali.blocks[b])
  }
  total / ali.blocks.length().to_double()
}

// ============================================================================
// MAF Parser
// ============================================================================

///|
/// Parse a MAF file from a string.
/// Returns a MAFAliignment with all blocks.
pub fn parse_maf(content : String) -> MAFAliignment {
  let mut result = new_maf_alignment()
  let lines = content.split("\n").to_array()
  let mut current_block : MAFBlock? = None
  
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    let trimmed = trim(line)
    
    if trimmed.length() == 0 {
      // Empty line - marks end of block
      match current_block {
        Some(block) => {
          block.is_complete = true
          result = maf_add_block(result, block)
          current_block = None
        }
        None => ()
      }
      continue
    }
    
    if starts_with(trimmed, 0, "##maf") {
      // Header line
      if starts_with(trimmed, 0, "##maf version=") {
        let parts = split_by_char(trimmed, ' '.to_int())
        for p = 0; p < parts.length(); p = p + 1 {
          let part = parts[p]
          if starts_with(part, 0, "version=") {
            result.version = part[8:part.length()].to_owned()
          } else if starts_with(part, 0, "score=") {
            let score_val = parse_double(part[6:part.length()].to_owned())
            result.score = score_val
          }
        }
      }
      continue
    }
    
    if starts_with(trimmed, 0, "#") {
      // Comment line
      match current_block {
        Some(block) => {
          let comment = MAFComment::{ key: trimmed }
          current_block = Some(maf_block_add_comment(block, comment))
        }
        None => {
          result.comments.push(trimmed)
        }
      }
      continue
    }
    
    if starts_with(trimmed, 0, "a ") {
      // Block attributes line
      current_block = Some(new_maf_block())
      let parts = split_by_char(trimmed, ' '.to_int())
      for p = 0; p < parts.length(); p = p + 1 {
        let part = parts[p]
        if starts_with(part, 0, "score=") {
          let score_val = parse_double(part[6:part.length()].to_owned())
          match current_block {
            Some(block) => {
              block.score = score_val
              current_block = Some(block)
            }
            None => ()
          }
        } else if contains(part, "=") {
          let kv = split_by_char(part, '='.to_int())
          if kv.length() == 2 {
            match current_block {
              Some(block) => {
                block.attributes[kv[0]] = kv[1]
                current_block = Some(block)
              }
              None => ()
            }
          }
        }
      }
      continue
    }
    
    if starts_with(trimmed, 0, "s ") || starts_with(trimmed, 0, "s\t") {
      // Sequence line
      if current_block is None {
        current_block = Some(new_maf_block())
      }
      let mut parts = split_by_tab(trimmed)
      if parts.length() < 6 {
        parts = split_by_char(trimmed, ' '.to_int())
      }
      if parts.length() >= 6 {
        let seq_name = parts[1]
        let start = parse_int(parts[2])
        let length = parse_int(parts[3])
        let strand = parts[4]
        let seq_size = parse_int(parts[5])
        let aligned_seq = if parts.length() >= 7 { parts[6] } else { "" }
        let seq = new_maf_sequence(seq_name, start, length, strand, seq_size, aligned_seq)
        match current_block {
          Some(block) => {
            current_block = Some(maf_block_add_seq(block, seq))
          }
          None => ()
        }
      }
      continue
    }
    
    if starts_with(trimmed, 0, "e ") {
      // Event line
      if current_block is None {
        current_block = Some(new_maf_block())
      }
      let parts = split_by_char(trimmed, ' '.to_int())
      if parts.length() >= 3 {
        let event_type = parts[1]
        let seq_name = parts[2]
        let position = if parts.length() >= 4 { parse_int(parts[3]) } else { 0 }
        let description = if parts.length() >= 5 { parts[4] } else { "" }
        let event = MAFEvent::{
          event_type,
          seq_name,
          position,
          description,
        }
        match current_block {
          Some(block) => {
            current_block = Some(maf_block_add_event(block, event))
          }
          None => ()
        }
      }
      continue
    }
    
    if starts_with(trimmed, 0, "b ") {
      // Block info line (marks start of new block with boundary info)
      match current_block {
        Some(block) => {
          block.is_complete = true
          result = maf_add_block(result, block)
        }
        None => ()
      }
      current_block = Some(new_maf_block())
      continue
    }
  }
  
  // Add final block if exists
  match current_block {
    Some(block) => {
      block.is_complete = true
      result = maf_add_block(result, block)
    }
    None => ()
  }
  
  result
}

// ============================================================================
// MAF Writer
// ============================================================================

///|
/// Write a MAFAliignment to a MAF-formatted string.
pub fn write_maf(ali : MAFAliignment) -> String {
  let mut result = ""
  
  // Write header
  result = result + "##maf version="
  result = result + ali.version
  match ali.score {
    Some(score_val) => {
      result = result + " score="
      result = result + score_val.to_string()
    }
    None => ()
  }
  result = result + "\n\n"
  
  // Write top-level comments
  for c = 0; c < ali.comments.length(); c = c + 1 {
    result = result + ali.comments[c]
    result = result + "\n"
  }
  
  // Write blocks
  for b = 0; b < ali.blocks.length(); b = b + 1 {
    let block = ali.blocks[b]
    // Write block attributes
    result = result + "a"
    match block.score {
      Some(score_val) => {
        result = result + " score="
        result = result + score_val.to_string()
      }
      None => ()
    }
    let attr_keys = block.attributes.keys().collect()
    for k = 0; k < attr_keys.length(); k = k + 1 {
      let key = attr_keys[k]
      result = result + " "
      result = result + key
      result = result + "="
      result = result + block.attributes[key]
    }
    result = result + "\n"
    
    // Write sequences
    for s = 0; s < block.sequences.length(); s = s + 1 {
      let seq = block.sequences[s]
      result = result + "s "
      result = result + seq.seq_name
      result = result + " "
      result = result + seq.start.to_string()
      result = result + " "
      result = result + seq.length.to_string()
      result = result + " "
      result = result + seq.strand
      result = result + " "
      result = result + seq.seq_size.to_string()
      result = result + " "
      result = result + seq.aligned_seq
      result = result + "\n"
    }
    
    // Write comments
    for c = 0; c < block.comments.length(); c = c + 1 {
      result = result + block.comments[c].key
      result = result + "\n"
    }
    
    // Write events
    for e = 0; e < block.events.length(); e = e + 1 {
      let event = block.events[e]
      result = result + "e "
      result = result + event.event_type
      result = result + " "
      result = result + event.seq_name
      result = result + " "
      result = result + event.position.to_string()
      if event.description.length() > 0 {
        result = result + " "
        result = result + event.description
      }
      result = result + "\n"
    }
    
    result = result + "\n"
  }
  
  result
}

// ============================================================================
// Block Operations
// ============================================================================

///|
/// Select sequences by name and create a new alignment.
pub fn maf_select_seqs(ali : MAFAliignment, names : Array[String]) -> MAFAliignment {
  let mut result = new_maf_alignment()
  result.version = ali.version
  result.score = ali.score
  result.comments = ali.comments
  
  for b = 0; b < ali.blocks.length(); b = b + 1 {
    let block = ali.blocks[b]
    let mut new_block = new_maf_block()
    new_block.score = block.score
    new_block.attributes = block.attributes
    new_block.events = block.events
    new_block.comments = block.comments
    new_block.is_complete = block.is_complete
    
    for s = 0; s < block.sequences.length(); s = s + 1 {
      let seq = block.sequences[s]
      if names.contains(seq.seq_name) {
        new_block = maf_block_add_seq(new_block, seq)
      }
    }
    
    if new_block.sequences.length() > 0 {
      result = maf_add_block(result, new_block)
    }
  }
  
  result
}

///|
/// Filter blocks by minimum alignment length.
pub fn maf_filter_by_length(ali : MAFAliignment, min_length : Int) -> MAFAliignment {
  let mut result = new_maf_alignment()
  result.version = ali.version
  result.score = ali.score
  result.comments = ali.comments
  
  for b = 0; b < ali.blocks.length(); b = b + 1 {
    let block = ali.blocks[b]
    let aln_len = maf_block_alignment_length(block)
    if aln_len >= min_length {
      result = maf_add_block(result, block)
    }
  }
  
  result
}

///|
/// Filter blocks by minimum percent identity.
pub fn maf_filter_by_identity(ali : MAFAliignment, min_identity : Double) -> MAFAliignment {
  let mut result = new_maf_alignment()
  result.version = ali.version
  result.score = ali.score
  result.comments = ali.comments
  
  for b = 0; b < ali.blocks.length(); b = b + 1 {
    let block = ali.blocks[b]
    let identity = maf_block_percent_identity(block)
    if identity >= min_identity {
      result = maf_add_block(result, block)
    }
  }
  
  result
}

///|
/// Merge consecutive blocks with the same sequences.
pub fn maf_merge_blocks(ali : MAFAliignment) -> MAFAliignment {
  if ali.blocks.length() <= 1 {
    return ali
  }
  
  let mut result = new_maf_alignment()
  result.version = ali.version
  result.score = ali.score
  result.comments = ali.comments
  
  let mut current_block = ali.blocks[0]
  
  for i = 1; i < ali.blocks.length(); i = i + 1 {
    let block = ali.blocks[i]
    
    // Check if blocks have the same sequences
    let same_seqs = current_block.sequences.length() == block.sequences.length()
    let can_merge = same_seqs && 
      current_block.events.length() == 0 && 
      block.events.length() == 0
    
    if can_merge {
      // Merge: concatenate sequences
      let mut new_block = new_maf_block()
      new_block.score = current_block.score
      new_block.attributes = current_block.attributes
      
      for j = 0; j < current_block.sequences.length(); j = j + 1 {
        let seq1 = current_block.sequences[j]
        let seq2 = block.sequences[j]
        let merged_seq = new_maf_sequence(
          seq1.seq_name,
          seq1.start,
          seq1.length + seq2.length,
          seq1.strand,
          seq1.seq_size,
          seq1.aligned_seq + seq2.aligned_seq,
        )
        new_block = maf_block_add_seq(new_block, merged_seq)
      }
      
      current_block = new_block
    } else {
      result = maf_add_block(result, current_block)
      current_block = block
    }
  }
  
  result = maf_add_block(result, current_block)
  result
}

///|
/// Get coverage for each sequence name across all blocks.
pub fn maf_sequence_coverage(ali : MAFAliignment) -> Map[String, Double] {
  let coverage = Map([], capacity=10)
  let seq_names = maf_all_seq_names(ali)
  
  for n = 0; n < seq_names.length(); n = n + 1 {
    let name = seq_names[n]
    let mut total_alignment = 0.0
    let mut total_length = 0.0
    
    for b = 0; b < ali.blocks.length(); b = b + 1 {
      let block = ali.blocks[b]
      for s = 0; s < block.sequences.length(); s = s + 1 {
        let seq = block.sequences[s]
        if seq.seq_name == name {
          total_alignment = total_alignment + seq.length.to_double()
          total_length = total_length + seq.seq_size.to_double()
        }
      }
    }
    
    if total_length > 0.0 {
      coverage[name] = (total_alignment / total_length) * 100.0
    } else {
      coverage[name] = 0.0
    }
  }
  
  coverage
}

// ============================================================================
// Statistics
// ============================================================================

///|
/// Calculate alignment statistics for a MAFAliignment.
pub struct MAFStats {
  num_blocks : Int
  total_length : Int
  num_sequences : Int
  avg_percent_identity : Double
  avg_block_coverage : Double
  min_block_length : Int
  max_block_length : Int
  median_block_length : Double
}

///|
/// Compute alignment statistics.
pub fn maf_compute_stats(ali : MAFAliignment) -> MAFStats {
  let num_blocks = ali.blocks.length()
  let total_length = maf_total_length(ali)
  let num_sequences = maf_all_seq_names(ali).length()
  let avg_percent_identity = maf_average_percent_identity(ali)
  
  let mut total_coverage = 0.0
  let lengths = Array::new()
  let mut min_len = 99999999
  let mut max_len = 0
  
  for b = 0; b < ali.blocks.length(); b = b + 1 {
    let block = ali.blocks[b]
    total_coverage = total_coverage + maf_block_coverage(block)
    let blen = maf_block_alignment_length(block)
    lengths.push(blen)
    if blen < min_len { min_len = blen }
    if blen > max_len { max_len = blen }
  }
  
  let avg_block_coverage = if num_blocks > 0 { total_coverage / num_blocks.to_double() } else { 0.0 }
  
  // Calculate median
  lengths.sort()
  let median = if lengths.length() > 0 {
    let mid = lengths.length() / 2
    if lengths.length() % 2 == 0 {
      (lengths[mid - 1] + lengths[mid]).to_double() / 2.0
    } else {
      lengths[mid].to_double()
    }
  } else {
    0.0
  }
  
  MAFStats::{
    num_blocks,
    total_length,
    num_sequences,
    avg_percent_identity,
    avg_block_coverage,
    min_block_length: if min_len == 99999999 { 0 } else { min_len },
    max_block_length: max_len,
    median_block_length: median,
  }
}

///|
/// Format MAFStats as a display string.
pub fn maf_stats_to_string(stats : MAFStats) -> String {
  let mut result = ""
  result = result + "MAF Alignment Statistics:\n"
  result = result + "  Number of blocks: "
  result = result + stats.num_blocks.to_string()
  result = result + "\n"
  result = result + "  Total alignment length: "
  result = result + stats.total_length.to_string()
  result = result + " columns\n"
  result = result + "  Number of sequences: "
  result = result + stats.num_sequences.to_string()
  result = result + "\n"
  result = result + "  Average percent identity: "
  result = result + stats.avg_percent_identity.to_string()
  result = result + "%\n"
  result = result + "  Average block coverage: "
  result = result + stats.avg_block_coverage.to_string()
  result = result + "\n"
  result = result + "  Min block length: "
  result = result + stats.min_block_length.to_string()
  result = result + "\n"
  result = result + "  Max block length: "
  result = result + stats.max_block_length.to_string()
  result = result + "\n"
  result = result + "  Median block length: "
  result = result + stats.median_block_length.to_string()
  result = result + "\n"
  result
}
