/// biotite_seq_core3.mbt — Biotite round 3: sequence core layer
/// 6 features: AlphabetMapper, common_alphabet, SequenceProfile,
/// PositionalSequence/PurePositionalSequence, Annotation/AnnotatedSequence,
/// Location with defect.

// #1 AlphabetMapper — symbol code conversion between alphabets
pub struct Biotite3AlphabetMapper {
  source_alphabet : Array[String]
  target_alphabet : Array[String]
  mapping : Array[Int]
}

pub fn biotite3_alphabet_mapper(source : Array[String], target : Array[String]) -> Biotite3AlphabetMapper {
  let mapping = Array(capacity=source.length())
  for i in 0..<source.length() {
    let mut target_idx = -1
    for j in 0..<target.length() {
      if target[j] == source[i] {
        target_idx = j
        break
      }
    }
    mapping.push(target_idx)
  }
  { source_alphabet: source, target_alphabet: target, mapping }
}

pub fn biotite3_alphabet_mapper_map(mapper : Biotite3AlphabetMapper, codes : Array[Int]) -> Array[Int] {
  let result = Array(capacity=codes.length())
  for i in 0..<codes.length() {
    if codes[i] >= 0 && codes[i] < mapper.mapping.length() {
      result.push(mapper.mapping[codes[i]])
    } else {
      result.push(-1)
    }
  }
  result
}

pub fn biotite3_alphabet_mapper_map_symbols(mapper : Biotite3AlphabetMapper, symbols : Array[String]) -> Array[String] {
  let result = Array(capacity=symbols.length())
  for i in 0..<symbols.length() {
    let mut src_idx = -1
    for j in 0..<mapper.source_alphabet.length() {
      if mapper.source_alphabet[j] == symbols[i] {
        src_idx = j
        break
      }
    }
    if src_idx >= 0 && mapper.mapping[src_idx] >= 0 && mapper.mapping[src_idx] < mapper.target_alphabet.length() {
      result.push(mapper.target_alphabet[mapper.mapping[src_idx]])
    } else {
      result.push("")
    }
  }
  result
}

// #2 common_alphabet — find common super alphabet
pub fn biotite3_common_alphabet(alphabets : Array[Array[String]]) -> Array[String] {
  let result : Array[String] = []
  let seen : Map[String, Bool] = Map([])
  for i in 0..<alphabets.length() {
    for j in 0..<alphabets[i].length() {
      let sym = alphabets[i][j]
      if !seen.contains(sym) {
        seen[sym] = true
        result.push(sym)
      }
    }
  }
  result
}

// #3 SequenceProfile — position-specific symbol frequencies
pub struct Biotite3SequenceProfile {
  alphabet : Array[String]
  frequencies : Array[Array[Double]]
  n_sequences : Int
}

pub fn biotite3_sequence_profile(alphabet : Array[String], sequences : Array[Array[Int]]) -> Biotite3SequenceProfile {
  let n_pos = {
    let mut max_len = 0
    for i in 0..<sequences.length() {
      if sequences[i].length() > max_len {
        max_len = sequences[i].length()
      }
    }
    max_len
  }
  let frequencies = Array(capacity=n_pos)
  for pos in 0..<n_pos {
    let col = Array(capacity=alphabet.length())
    for _s in 0..<alphabet.length() {
      col.push(0.0)
    }
    let mut count = 0
    for seq_idx in 0..<sequences.length() {
      if pos < sequences[seq_idx].length() {
        let code = sequences[seq_idx][pos]
        if code >= 0 && code < alphabet.length() {
          col[code] = col[code] + 1.0
        }
        count = count + 1
      }
    }
    if count > 0 {
      for s in 0..<alphabet.length() {
        col[s] = col[s] / count.to_double()
      }
    }
    frequencies.push(col)
  }
  { alphabet, frequencies, n_sequences: sequences.length() }
}

pub fn biotite3_profile_consensus(profile : Biotite3SequenceProfile) -> Array[Int] {
  let result = Array(capacity=profile.frequencies.length())
  for pos in 0..<profile.frequencies.length() {
    let mut best_idx = 0
    let mut best_freq = -1.0
    for s in 0..<profile.frequencies[pos].length() {
      if profile.frequencies[pos][s] > best_freq {
        best_freq = profile.frequencies[pos][s]
        best_idx = s
      }
    }
    result.push(best_idx)
  }
  result
}

pub fn biotite3_profile_entropy(profile : Biotite3SequenceProfile) -> Array[Double] {
  let result = Array(capacity=profile.frequencies.length())
  for pos in 0..<profile.frequencies.length() {
    let mut h = 0.0
    for s in 0..<profile.frequencies[pos].length() {
      let p = profile.frequencies[pos][s]
      if p > 0.0 {
        h = h - p * @math.ln(p)
      }
    }
    result.push(h)
  }
  result
}

// #4 PositionalSequence / PurePositionalSequence
pub struct Biotite3PositionalSequence {
  symbols : Array[String]
  positions : Array[Double]
}

pub fn biotite3_positional_sequence(symbols : Array[String], positions : Array[Double]) -> Biotite3PositionalSequence {
  { symbols, positions }
}

pub fn biotite3_pure_positional_sequence(n : Int) -> Biotite3PositionalSequence {
  let symbols = Array(capacity=n)
  let positions = Array(capacity=n)
  for i in 0..<n {
    symbols.push(i.to_string())
    positions.push(i.to_double())
  }
  { symbols, positions }
}

pub fn biotite3_positional_sequence_get(seq : Biotite3PositionalSequence, index : Int) -> (String, Double) {
  (seq.symbols[index], seq.positions[index])
}

pub fn biotite3_positional_sequence_slice(seq : Biotite3PositionalSequence, start : Int, stop : Int) -> Biotite3PositionalSequence {
  let symbols = Array(capacity=stop - start)
  let positions = Array(capacity=stop - start)
  for i in start..<stop {
    if i < seq.symbols.length() {
      symbols.push(seq.symbols[i])
      positions.push(seq.positions[i])
    }
  }
  { symbols, positions }
}

// #5 Annotation + AnnotatedSequence
pub struct Biotite3Location {
  first : Int
  last : Int
  strand : Int
  defect : Int
}

pub fn biotite3_location(first : Int, last : Int, strand : Int, defect : Int) -> Biotite3Location {
  { first, last, strand, defect }
}

pub struct Biotite3Feature {
  key : String
  quals : Array[(String, String)]
  loc : Biotite3Location
}

pub fn biotite3_feature(key : String, loc : Biotite3Location, quals : Array[(String, String)]) -> Biotite3Feature {
  { key, quals, loc }
}

pub struct Biotite3Annotation {
  features : Array[Biotite3Feature]
}

pub fn biotite3_annotation() -> Biotite3Annotation {
  { features: [] }
}

pub fn biotite3_annotation_add(annot : Biotite3Annotation, feature : Biotite3Feature) -> Biotite3Annotation {
  let new_features = Array(capacity=annot.features.length() + 1)
  for i in 0..<annot.features.length() {
    new_features.push(annot.features[i])
  }
  new_features.push(feature)
  { features: new_features }
}

pub fn biotite3_annotation_merge(annot : Biotite3Annotation, key : String) -> Biotite3Annotation {
  let merged = Array(capacity=0)
  let used : Array[Bool] = Array(capacity=annot.features.length())
  for _i in 0..<annot.features.length() {
    used.push(false)
  }
  for i in 0..<annot.features.length() {
    if used[i] || annot.features[i].key != key {
      continue
    }
    let mut min_first = annot.features[i].loc.first
    let mut max_last = annot.features[i].loc.last
    let strand = annot.features[i].loc.strand
    used[i] = true
    let _start_j = i + 1
    for j in _start_j..<annot.features.length() {
      if !used[j] && annot.features[j].key == key && annot.features[j].loc.strand == strand {
        if annot.features[j].loc.first < min_first {
          min_first = annot.features[j].loc.first
        }
        if annot.features[j].loc.last > max_last {
          max_last = annot.features[j].loc.last
        }
        used[j] = true
      }
    }
    merged.push(biotite3_feature(key, biotite3_location(min_first, max_last, strand, 0), []))
  }
  for i in 0..<annot.features.length() {
    if !used[i] {
      merged.push(annot.features[i])
    }
  }
  { features: merged }
}

pub fn biotite3_annotation_split(annot : Biotite3Annotation, pos : Int) -> (Biotite3Annotation, Biotite3Annotation) {
  let left = Array(capacity=0)
  let right = Array(capacity=0)
  for i in 0..<annot.features.length() {
    let f = annot.features[i]
    if f.loc.last < pos {
      left.push(f)
    } else if f.loc.first >= pos {
      right.push(f)
    } else {
      left.push(biotite3_feature(f.key, biotite3_location(f.loc.first, pos - 1, f.loc.strand, f.loc.defect), f.quals))
      right.push(biotite3_feature(f.key, biotite3_location(pos, f.loc.last, f.loc.strand, f.loc.defect), f.quals))
    }
  }
  ({ features: left }, { features: right })
}

pub struct Biotite3AnnotatedSequence {
  sequence : String
  annotation : Biotite3Annotation
}

pub fn biotite3_annotated_sequence(sequence : String, annotation : Biotite3Annotation) -> Biotite3AnnotatedSequence {
  { sequence, annotation }
}

pub fn biotite3_annotated_sequence_get_features(seq : Biotite3AnnotatedSequence, key : String) -> Array[Biotite3Feature] {
  let result = Array(capacity=0)
  for i in 0..<seq.annotation.features.length() {
    if seq.annotation.features[i].key == key {
      result.push(seq.annotation.features[i])
    }
  }
  result
}

pub fn biotite3_annotated_sequence_subsequence(seq : Biotite3AnnotatedSequence, start : Int, stop : Int) -> Biotite3AnnotatedSequence {
  let sub_seq = if stop <= seq.sequence.length() { seq.sequence[start:stop].to_owned() } else { seq.sequence[start:].to_owned() }
  let (left, _right) = biotite3_annotation_split(seq.annotation, start)
  biotite3_annotated_sequence(sub_seq, left)
}

// #6 Location defect constants and helpers
pub fn biotite3_location_defect_none() -> Int { 0 }
pub fn biotite3_location_defect_beyond_left() -> Int { 1 }
pub fn biotite3_location_defect_beyond_right() -> Int { 2 }
pub fn biotite3_location_defect_between() -> Int { 4 }

pub fn biotite3_location_overlap(loc1 : Biotite3Location, loc2 : Biotite3Location) -> Bool {
  !(loc1.last < loc2.first || loc2.last < loc1.first)
}

pub fn biotite3_location_contains(outer : Biotite3Location, inner : Biotite3Location) -> Bool {
  outer.first <= inner.first && inner.last <= outer.last
}

pub fn biotite3_location_to_string(loc : Biotite3Location) -> String {
  let strand_str = if loc.strand < 0 { "(-)" } else if loc.strand > 0 { "(+)" } else { "" }
  let defect_str = if loc.defect != 0 { "<" + loc.defect.to_string() + ">" } else { "" }
  defect_str + loc.first.to_string() + ".." + loc.last.to_string() + strand_str
}

pub fn _biotite_seq_core3_dummy() -> Unit {
  let _ = Biotite3AlphabetMapper::{ source_alphabet: [], target_alphabet: [], mapping: [] }
  let _ = Biotite3SequenceProfile::{ alphabet: [], frequencies: [], n_sequences: 0 }
  let _ = Biotite3PositionalSequence::{ symbols: [], positions: [] }
  let _ = Biotite3Location::{ first: 0, last: 0, strand: 0, defect: 0 }
  let _ = Biotite3Feature::{ key: "", quals: [], loc: biotite3_location(0, 0, 0, 0) }
  let _ = Biotite3Annotation::{ features: [] }
  let _ = Biotite3AnnotatedSequence::{ sequence: "", annotation: biotite3_annotation() }
}
