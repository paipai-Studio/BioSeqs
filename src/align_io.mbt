///|
/// AlignIO - Alignment format parsing inspired by Biopython Bio.AlignIO.

///|
pub struct Alignment {
  sequences : Array[(String, String)]
  alphabet : String
  description : String
}

///|
pub fn alignment(
  sequences : Array[(String, String)],
  alphabet : String,
  description : String,
) -> Alignment {
  Alignment::{ sequences, alphabet, description }
}

///|
pub fn parse_clustalw(content : String) -> Alignment {
  let sequences : Array[(String, String)] = Array::new()
  let name_map : Map[String, String] = Map([], capacity=10)
  let n = content.length()

  let mut i = 0
  while i < n {
    let line_start = i
    while i < n && content.unsafe_get(i).to_int() != '\n'.to_int() {
      i = i + 1
    }
    let line = substring(content, line_start, i - line_start)
    i = i + 1

    if line.length() == 0 {
      continue
    }

    let trimmed = trim_string(line)
    if trimmed.length() == 0 {
      continue
    }

    let skip = starts_with(trimmed, 0, "//") ||
      starts_with(trimmed, 0, "CLUSTAL") ||
      starts_with(trimmed, 0, "MUSCLE")
    if skip {
      continue
    }

    let parts = split_by_whitespace(trimmed)
    if parts.length() >= 2 {
      let name = parts[0]
      let seq = parts[parts.length() - 1]

      match name_map.get(name) {
        Some(existing) => name_map[name] = existing + seq
        None => name_map[name] = seq
      }
    }
  }

  let keys = name_map.keys().collect()
  for i = 0; i < keys.length(); i = i + 1 {
    sequences.push((keys[i], name_map[keys[i]]))
  }

  Alignment::{ sequences, alphabet: "DNA", description: "ClustalW alignment" }
}

///|
pub fn parse_fasta_alignment(content : String) -> Alignment {
  let sequences : Array[(String, String)] = Array::new()
  let n = content.length()

  let mut current_name = ""
  let mut current_seq = ""
  let mut i = 0

  while i < n {
    if content.unsafe_get(i).to_int() == '>'.to_int() {
      if current_name.length() > 0 {
        sequences.push((current_name, current_seq))
      }
      i = i + 1
      let name_start = i
      while i < n && content.unsafe_get(i).to_int() != '\n'.to_int() {
        i = i + 1
      }
      current_name = substring(content, name_start, i - name_start)
      current_seq = ""
      i = i + 1
    } else if content.unsafe_get(i).to_int() != '\n'.to_int() &&
      content.unsafe_get(i).to_int() != '\r'.to_int() {
      let seq_start = i
      while i < n && content.unsafe_get(i).to_int() != '\n'.to_int() {
        i = i + 1
      }
      current_seq = current_seq + substring(content, seq_start, i - seq_start)
      i = i + 1
    } else {
      i = i + 1
    }
  }

  if current_name.length() > 0 {
    sequences.push((current_name, current_seq))
  }

  Alignment::{ sequences, alphabet: "DNA", description: "FASTA alignment" }
}

///|
pub fn parse_stockholm(content : String) -> Alignment {
  let sequences : Array[(String, String)] = Array::new()
  let name_map : Map[String, String] = Map([], capacity=10)
  let n = content.length()

  let mut i = 0
  while i < n {
    let line_start = i
    while i < n && content.unsafe_get(i).to_int() != '\n'.to_int() {
      i = i + 1
    }
    let line = substring(content, line_start, i - line_start)
    i = i + 1

    if line.length() == 0 {
      continue
    }

    let trimmed = trim_string(line)
    if trimmed.length() == 0 {
      continue
    }

    let skip = starts_with(trimmed, 0, "#") || starts_with(trimmed, 0, "//")
    if skip {
      continue
    }

    let parts = split_by_whitespace(trimmed)
    if parts.length() >= 2 {
      let name = parts[0]
      let seq = parts[parts.length() - 1]

      match name_map.get(name) {
        Some(existing) => name_map[name] = existing + seq
        None => name_map[name] = seq
      }
    }
  }

  let keys = name_map.keys().collect()
  for i = 0; i < keys.length(); i = i + 1 {
    sequences.push((keys[i], name_map[keys[i]]))
  }

  Alignment::{ sequences, alphabet: "DNA", description: "Stockholm alignment" }
}

///|
pub fn write_clustalw(align : Alignment) -> String {
  let mut result = "CLUSTAL W (1.83) multiple sequence alignment\n\n"

  let mut max_name_len = 0
  for i = 0; i < align.sequences.length(); i = i + 1 {
    if align.sequences[i].0.length() > max_name_len {
      max_name_len = align.sequences[i].0.length()
    }
  }

  let seq_len = if align.sequences.length() > 0 {
    align.sequences[0].1.length()
  } else {
    0
  }
  let block_size = 60

  for start = 0; start < seq_len; start = start + block_size {
    for i = 0; i < align.sequences.length(); i = i + 1 {
      let name = align.sequences[i].0
      let seq = align.sequences[i].1

      let padding = max_name_len - name.length() + 4
      let mut line = name

      for p = 0; p < padding; p = p + 1 {
        line = line + " "
      }

      let end = if start + block_size < seq_len {
        start + block_size
      } else {
        seq_len
      }
      line = line + substring(seq, start, end - start)

      if i < align.sequences.length() - 1 || start + block_size < seq_len {
        line = line + "\n"
      }

      result = result + line
    }

    if start + block_size < seq_len {
      result = result + "\n\n"
    }
  }

  result = result + "\n//\n"
  result
}

///|
pub fn write_fasta_alignment(align : Alignment) -> String {
  let mut result = ""

  for i = 0; i < align.sequences.length(); i = i + 1 {
    result = result + ">" + align.sequences[i].0 + "\n"

    let seq = align.sequences[i].1
    let line_len = 80

    for start = 0; start < seq.length(); start = start + line_len {
      let end = if start + line_len < seq.length() {
        start + line_len
      } else {
        seq.length()
      }
      result = result + substring(seq, start, end - start) + "\n"
    }
  }

  result
}

///|
pub fn alignment_length(align : Alignment) -> Int {
  if align.sequences.length() > 0 {
    align.sequences[0].1.length()
  } else {
    0
  }
}

///|
pub fn num_sequences(align : Alignment) -> Int {
  align.sequences.length()
}

///|
pub fn get_sequence(align : Alignment, name : String) -> String? {
  for i = 0; i < align.sequences.length(); i = i + 1 {
    if align.sequences[i].0 == name {
      return Some(align.sequences[i].1)
    }
  }
  None
}

///|
pub fn column_consensus(align : Alignment, column : Int) -> String {
  if align.sequences.length() == 0 || column >= align.sequences[0].1.length() {
    return "N"
  }

  let counts : Map[String, Int] = Map([], capacity=4)

  for i = 0; i < align.sequences.length(); i = i + 1 {
    let c = substring(align.sequences[i].1, column, 1)

    match counts.get(c) {
      Some(existing) => counts[c] = existing + 1
      None => counts[c] = 1
    }
  }

  let mut max_count = 0
  let mut consensus = "N"

  let keys = counts.keys().collect()
  for i = 0; i < keys.length(); i = i + 1 {
    if counts[keys[i]] > max_count {
      max_count = counts[keys[i]]
      consensus = keys[i]
    }
  }

  consensus
}

///|
pub fn consensus_sequence(align : Alignment) -> String {
  let len = alignment_length(align)
  let mut consensus = ""

  for i = 0; i < len; i = i + 1 {
    consensus = consensus + column_consensus(align, i)
  }

  consensus
}

///|
fn split_by_whitespace(s : String) -> Array[String] {
  let result : Array[String] = Array::new()
  let mut current_word = ""

  for i = 0; i < s.length(); i = i + 1 {
    let c = char_at(s, i)

    if c == " " || c == "\t" {
      if current_word.length() > 0 {
        result.push(current_word)
        current_word = ""
      }
    } else {
      current_word = current_word + c
    }
  }

  if current_word.length() > 0 {
    result.push(current_word)
  }

  result
}

///|
fn trim_string(s : String) -> String {
  let mut start = 0
  while start < s.length() {
    let c = char_at(s, start)
    if c != " " && c != "\t" && c != "\n" && c != "\r" {
      break
    }
    start = start + 1
  }

  let mut end = s.length()
  while end > start {
    let c = char_at(s, end - 1)
    if c != " " && c != "\t" && c != "\n" && c != "\r" {
      break
    }
    end = end - 1
  }

  substring(s, start, end - start)
}
