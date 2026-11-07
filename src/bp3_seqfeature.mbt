///|
/// bp3_seqfeature.mbt — Biopython round 3: CompoundLocation + InvertedRepeatDetection
/// #13 CompoundLocation — compound locations for multi-region features
/// #15 InvertedRepeatDetection — detect inverted repeats / hairpins

///|
/// #13 CompoundLocation — compound locations for features spanning multiple regions
pub struct Bp3SimpleLocation {
  start : Int32
  end : Int32
  strand : Int32
}

///|
pub fn Bp3SimpleLocation::new(
  start : Int32,
  end : Int32,
  strand~ : Int32 = 1,
) -> Bp3SimpleLocation {
  { start, end, strand, }
}

///|
pub fn Bp3SimpleLocation::length(self : Bp3SimpleLocation) -> Int32 {
  if self.end > self.start {
    self.end - self.start
  } else {
    0
  }
}

///|
pub struct Bp3CompoundLocation {
  parts : Array[Bp3SimpleLocation]
  operator : String
}

///|
pub fn Bp3CompoundLocation::new(
  parts : Array[Bp3SimpleLocation],
  operator~ : String = "join",
) -> Bp3CompoundLocation {
  { parts, operator, }
}

///|
pub fn Bp3CompoundLocation::length(self : Bp3CompoundLocation) -> Int32 {
  let mut total = 0
  for i in 0..<self.parts.length() {
    total = total + self.parts[i].length()
  }
  total
}

///|
pub fn Bp3CompoundLocation::start(self : Bp3CompoundLocation) -> Int32 {
  if self.parts.length() == 0 {
    return 0
  }
  let mut min_start = self.parts[0].start
  for i in 1..<self.parts.length() {
    if self.parts[i].start < min_start {
      min_start = self.parts[i].start
    }
  }
  min_start
}

///|
pub fn Bp3CompoundLocation::end(self : Bp3CompoundLocation) -> Int32 {
  if self.parts.length() == 0 {
    return 0
  }
  let mut max_end = self.parts[0].end
  for i in 1..<self.parts.length() {
    if self.parts[i].end > max_end {
      max_end = self.parts[i].end
    }
  }
  max_end
}

///|
pub fn Bp3CompoundLocation::is_contiguous(self : Bp3CompoundLocation) -> Bool {
  self.parts.length() <= 1
}

///|
pub fn Bp3CompoundLocation::extract(
  self : Bp3CompoundLocation,
  sequence : String,
) -> String {
  let sb = StringBuilder()
  for i in 0..<self.parts.length() {
    let p = self.parts[i]
    let s = if p.start >= 0 && p.end <= sequence.length() && p.start < p.end {
      sequence[p.start:p.end].to_owned()
    } else {
      ""
    }
    let _ = sb.write(s)
  }
  let result = sb.to_string()
  if self.parts.length() > 0 && self.parts[0].strand == -1 {
    bp3_reverse_complement(result)
  } else {
    result
  }
}

///|
pub fn Bp3CompoundLocation::to_genbank(self : Bp3CompoundLocation) -> String {
  let sb = StringBuilder()
  let _ = sb.write(self.operator + "(")
  for i in 0..<self.parts.length() {
    let p = self.parts[i]
    if i > 0 {
      let _ = sb.write(",")
    }
    let _ = sb.write(p.start.to_string() + ".." + p.end.to_string())
  }
  let _ = sb.write(")")
  if self.parts.length() > 0 && self.parts[0].strand == -1 {
    let _ = sb.write("::complement")
  }
  sb.to_string()
}

///|
fn bp3_reverse_complement(seq : String) -> String {
  let comp : Map[Int32, Int32] = Map([
    (65, 84),
    (84, 65),
    (71, 67),
    (67, 71),
    (97, 116),
    (116, 97),
    (103, 99),
    (99, 103),
  ])
  let sb = StringBuilder()
  for i in 0..<seq.length() {
    let ch = seq[seq.length() - 1 - i]
    let v = comp.get(ch.to_int())
    match v {
      Some(c) => {
        let _ = sb.write_char(Int32::unsafe_to_char(c))
      }
      None => {
        let _ = sb.write_char(ch.unsafe_to_char())
      }
    }
  }
  sb.to_string()
}

///|
/// #15 InvertedRepeatDetection — detect inverted repeats / hairpin structures
pub struct Bp3InvertedRepeat {
  start1 : Int32
  end1 : Int32
  start2 : Int32
  end2 : Int32
  stem_length : Int32
  loop_length : Int32
  sequence : String
  energy : Double
}

///|
pub fn bp3_find_inverted_repeats(
  sequence : String,
  min_stem~ : Int32 = 6,
  max_loop~ : Int32 = 50,
  min_loop~ : Int32 = 3,
) -> Array[Bp3InvertedRepeat] {
  let results = Array(capacity=0)
  let n = sequence.length()
  let comp : Map[Int32, Int32] = Map([
    (65, 84),
    (84, 65),
    (71, 67),
    (67, 71),
    (97, 116),
    (116, 97),
    (103, 99),
    (99, 103),
  ])
  for i in 0..<n {
    let max_stem = (n - i) / 2
    if max_stem < min_stem {
      continue
    }
    for stem_len in min_stem..<=max_stem {
      for loop_len in min_loop..<=max_loop {
        let j = i + stem_len + loop_len
        if j + stem_len > n {
          break
        }
        let mut match_count = 0
        for k in 0..<stem_len {
          let c1 = sequence[i + k].to_int()
          let c2 = sequence[j + stem_len - 1 - k].to_int()
          let expected = comp.get(c1)
          match expected {
            Some(e) => if c2 == e { match_count = match_count + 1 }
            None => ()
          }
        }
        if match_count == stem_len {
          let stem_seq = sequence[i:i + stem_len].to_owned()
          let energy = bp3_hairpin_energy(stem_seq, loop_len)
          results.push({
            start1: i,
            end1: i + stem_len,
            start2: j,
            end2: j + stem_len,
            stem_length: stem_len,
            loop_length: loop_len,
            sequence: stem_seq,
            energy,
          })
          break
        }
      }
    }
  }
  results
}

///|
fn bp3_hairpin_energy(stem : String, loop_len : Int32) -> Double {
  let mut energy = 0.0
  let n = stem.length()
  let _half_n = n / 2
  for i in 0..<_half_n {
    let c1 = stem[i:i + 1].to_owned()
    let c2 = stem[n - 1 - i:n - i].to_owned()
    let pair = if c1 < c2 { c1 + c2 } else { c2 + c1 }
    match pair {
      "AT" => energy = energy - 2.0
      "CG" => energy = energy - 3.0
      "GT" => energy = energy - 1.0
      _ => energy = energy + 0.5
    }
  }
  if loop_len < 3 {
    energy = energy + 3.0
  } else if loop_len < 5 {
    energy = energy + 1.0
  }
  energy
}

///|
pub fn bp3_find_palindromes(
  sequence : String,
  min_len~ : Int32 = 6,
  max_len~ : Int32 = 20,
) -> Array[(Int32, Int32, String)] {
  let results = Array(capacity=0)
  let n = sequence.length()
  let comp : Map[Int32, Int32] = Map([
    (65, 84),
    (84, 65),
    (71, 67),
    (67, 71),
    (97, 116),
    (116, 97),
    (103, 99),
    (99, 103),
  ])
  for i in 0..<n {
    let max_l = if i + max_len <= n { max_len } else { n - i }
    for l in min_len..<=max_l {
      let half = l / 2
      let mut is_pal = true
      for k in 0..<half {
        let c1 = sequence[i + k].to_int()
        let c2 = sequence[i + l - 1 - k].to_int()
        let expected = comp.get(c1)
        match expected {
          Some(e) => if c2 != e { is_pal = false }
          None => is_pal = false
        }
      }
      if is_pal {
        let s = sequence[i:i + l].to_owned()
        results.push((i, i + l, s))
      }
    }
  }
  results
}

///|
pub fn _bp3_seqfeature_dummy() -> Unit {
  let _ = Bp3SimpleLocation::{ start: 0, end: 0, strand: 0, }
  let _ = Bp3CompoundLocation::{ parts: [], operator: "", }
  let _ = Bp3InvertedRepeat::{
    start1: 0,
    end1: 0,
    start2: 0,
    end2: 0,
    stem_length: 0,
    loop_length: 0,
    sequence: "",
    energy: 0.0,
  }
}
