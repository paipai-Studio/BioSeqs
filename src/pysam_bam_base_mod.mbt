/// pysam_bam_base_mod.mbt — Base modification (MM/ML tags) parsing
/// Items 9-10: MM/ML标签解析, 碱基修饰迭代

///|
/// A single base modification call.
pub struct PysamBaseMod {
  position : Int32 // 0-based position in read
  base : Int32 // canonical base (byte value)
  mod_type : String // "m" (5mC), "h" (5hmC), "o" (5fC), "a" (6mA), "n" (5caC), etc.
  mod_prob : Int32 // probability (0-255 from ML tag)
  strand : Bool // true = forward
}

///|
/// Parsed MM tag entry.
pub struct PysamMMEntry {
  canonical_base : String // "C", "A", "G", "T"
  mod_type : String // "m", "h", "o", "a", etc.
  strand : Bool // true = forward, false = reverse
  positions : Array[Int32] // delta-encoded positions
}

///|
/// Parse MM tag (modified bases) string.
/// Format: "C+m,5,10,20;A+a,3,7;C-h,1,4"
pub fn pysam_parse_mm_tag(mm_tag : String) -> Array[PysamMMEntry] {
  let entries : Array[PysamMMEntry] = Array(capacity=0)
  let parts = mm_tag.split(";").to_array()
  for i = 0; i < parts.length(); i = i + 1 {
    let part = parts[i].to_owned()
    if part.length() < 3 {
      continue
    }
    // Parse canonical_base + mod_type + strand
    let canonical_base = part[0:1].to_owned()
    let mod_char = part[1].to_int()
    let mod_type = part[2:3].to_owned()
    let strand = mod_char == 43 // '+' = 43
    // Parse positions (delta-encoded)
    let positions : Array[Int32] = Array(capacity=0)
    let comma_idx = part.find(",")
    match comma_idx {
      None => continue
      Some(ci) => {
        let pos_str = part[ci + 1:].to_owned()
        let pos_parts = pos_str.split(",").to_array()
        let mut prev = 0
        for p = 0; p < pos_parts.length(); p = p + 1 {
          let delta = parse_int(pos_parts[p].to_owned())
          prev = prev + delta
          positions.push(prev)
        }
      }
    }
    entries.push({ canonical_base, mod_type, strand, positions, })
  }
  entries
}

///|
/// Parse ML tag (modification probabilities) string.
/// Returns array of probability values (0-255).
pub fn pysam_parse_ml_tag(ml_tag : String) -> Array[Int32] {
  let result : Array[Int32] = Array(capacity=0)
  let parts = ml_tag.split(",").to_array()
  for i = 0; i < parts.length(); i = i + 1 {
    result.push(parse_int(parts[i].to_owned()))
  }
  result
}

///|
/// Iterate base modifications from MM and ML tags.
/// Returns all modification calls with positions and probabilities.
pub fn pysam_iter_base_mods(
  mm_tag : String,
  ml_tag : String,
  read_seq : String,
) -> Array[PysamBaseMod] {
  let mm_entries = pysam_parse_mm_tag(mm_tag)
  let ml_probs = pysam_parse_ml_tag(ml_tag)
  let result : Array[PysamBaseMod] = Array(capacity=0)
  let mut prob_idx = 0
  for e = 0; e < mm_entries.length(); e = e + 1 {
    let entry = mm_entries[e]
    let canonical_byte = entry.canonical_base[0].to_int()
    // Walk through read to find canonical base positions
    let mut canonical_pos = 0
    for r = 0; r < read_seq.length(); r = r + 1 {
      if read_seq[r].to_int() == canonical_byte {
        // Check if this position is in the modification list
        for p = 0; p < entry.positions.length(); p = p + 1 {
          if entry.positions[p] == canonical_pos {
            let prob = if prob_idx < ml_probs.length() {
              ml_probs[prob_idx]
            } else {
              0
            }
            prob_idx = prob_idx + 1
            result.push({
              position: r,
              base: canonical_byte,
              mod_type: entry.mod_type,
              mod_prob: prob,
              strand: entry.strand,
            })
          }
        }
        canonical_pos = canonical_pos + 1
      }
    }
  }
  result
}

///|
/// Count modifications by type.
pub fn pysam_count_base_mods(
  mods : Array[PysamBaseMod],
) -> Array[(String, Int32)] {
  let types : Array[String] = Array(capacity=0)
  let counts : Array[Int32] = Array(capacity=0)
  for i = 0; i < mods.length(); i = i + 1 {
    let mt = mods[i].mod_type
    let mut found = false
    for j = 0; j < types.length(); j = j + 1 {
      if types[j] == mt {
        counts[j] = counts[j] + 1
        found = true
        break
      }
    }
    if !found {
      types.push(mt)
      counts.push(1)
    }
  }
  let result : Array[(String, Int32)] = Array(capacity=0)
  for j = 0; j < types.length(); j = j + 1 {
    result.push((types[j], counts[j]))
  }
  result
}

///|
/// Filter modifications by probability threshold.
pub fn pysam_filter_base_mods(
  mods : Array[PysamBaseMod],
  min_prob~ : Int32 = 128,
) -> Array[PysamBaseMod] {
  let result : Array[PysamBaseMod] = Array(capacity=0)
  for i = 0; i < mods.length(); i = i + 1 {
    if mods[i].mod_prob >= min_prob {
      result.push(mods[i])
    }
  }
  result
}

///|
/// Generate modification summary string.
pub fn pysam_base_mod_summary(mods : Array[PysamBaseMod]) -> String {
  let counts = pysam_count_base_mods(mods)
  let sb = StringBuilder()
  let _ = sb.write("Total modifications: " + mods.length().to_string() + "\n")
  for i = 0; i < counts.length(); i = i + 1 {
    let (mt, count) = counts[i]
    let _ = sb.write("  " + mt + ": " + count.to_string() + "\n")
  }
  sb.to_string()
}
