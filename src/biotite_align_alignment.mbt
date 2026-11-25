///|
/// biotite_align_alignment.mbt — a lightweight Alignment object over
/// already aligned strings.

///|
pub(all) struct BiotiteAlignment {
  names : Array[String]
  sequences : Array[String]
  score : Double
}

///|
/// Build an alignment from already aligned sequences.
/// The score is the fraction of fully conserved columns.
pub fn biotite_alignment_from_strings(
  sequences : Array[String],
  names : Array[String],
) -> BiotiteAlignment {
  let n = sequences.length()
  if n == 0 {
    return { names, sequences, score: 0.0, }
  }
  let len = sequences[0].length()
  let mut identical = 0
  for col in 0..<len {
    let mut same = true
    for i in 1..<n {
      if col >= sequences[i].length() ||
        sequences[i][col:col + 1].to_owned() !=
        sequences[0][col:col + 1].to_owned() {
        same = false
        break
      }
    }
    if same {
      identical = identical + 1
    }
  }
  let score = if len > 0 {
    identical.to_double() / len.to_double()
  } else {
    0.0
  }
  { names, sequences, score, }
}

///|
pub fn biotite_alignment_length(alignment : BiotiteAlignment) -> Int32 {
  if alignment.sequences.length() == 0 {
    0
  } else {
    alignment.sequences[0].length()
  }
}

///|
pub fn biotite_alignment_num_sequences(alignment : BiotiteAlignment) -> Int32 {
  alignment.sequences.length()
}

///|
pub fn biotite_alignment_get_sequence(
  alignment : BiotiteAlignment,
  index : Int32,
) -> String? {
  if index >= 0 && index < alignment.sequences.length() {
    Some(alignment.sequences[index])
  } else {
    None
  }
}

///|
pub fn biotite_alignment_get_name(
  alignment : BiotiteAlignment,
  index : Int32,
) -> String? {
  if index >= 0 && index < alignment.names.length() {
    Some(alignment.names[index])
  } else {
    None
  }
}

///|
/// Average pairwise identity over all sequence pairs (gaps compared literally).
pub fn biotite_alignment_identity(alignment : BiotiteAlignment) -> Double {
  let n = alignment.sequences.length()
  if n < 2 {
    return 1.0
  }
  let len = biotite_alignment_length(alignment)
  let mut total = 0.0
  let mut pairs = 0
  for i in 0..<n {
    for j in (i + 1)..<n {
      let mut matches = 0
      for col in 0..<len {
        if col < alignment.sequences[i].length() &&
          col < alignment.sequences[j].length() &&
          alignment.sequences[i][col:col + 1].to_owned() ==
          alignment.sequences[j][col:col + 1].to_owned() {
          matches = matches + 1
        }
      }
      total = total + matches.to_double() / len.to_double()
      pairs = pairs + 1
    }
  }
  if pairs == 0 {
    0.0
  } else {
    total / pairs.to_double()
  }
}

///|
pub fn biotite_alignment_is_gapped(alignment : BiotiteAlignment) -> Bool {
  for seq in alignment.sequences {
    for ch in seq {
      if ch == '-' {
        return true
      }
    }
  }
  false
}
