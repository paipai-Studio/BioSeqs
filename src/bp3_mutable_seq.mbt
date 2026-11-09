///|
/// bp3_mutable_seq.mbt — Biopython round 3: MutableSeq + UnknownSeq
/// #1 MutableSeq — in-place mutable sequence operations
/// #2 UnknownSeq — memory-efficient unknown-content sequence

///|
/// #1 MutableSeq — editable sequence with in-place operations
pub struct Bp3MutableSeq {
  data : Array[Int32]
  alphabet : String
}

///|
pub fn Bp3MutableSeq::new(
  seq : String,
  alphabet~ : String = "ACGT",
) -> Bp3MutableSeq {
  let data = Array(capacity=seq.length())
  for i in 0..<seq.length() {
    data.push(seq[i].to_int())
  }
  { data, alphabet, }
}

///|
pub fn Bp3MutableSeq::to_string(self : Bp3MutableSeq) -> String {
  let sb = StringBuilder()
  for i in 0..<self.data.length() {
    let _ = sb.write_char(self.data[i].unsafe_to_char())
  }
  sb.to_string()
}

///|
pub fn Bp3MutableSeq::length(self : Bp3MutableSeq) -> Int32 {
  self.data.length()
}

///|
pub fn Bp3MutableSeq::append(self : Bp3MutableSeq, ch : Int32) -> Unit {
  self.data.push(ch)
}

///|
pub fn Bp3MutableSeq::insert(
  self : Bp3MutableSeq,
  pos : Int32,
  ch : Int32,
) -> Unit {
  if pos >= 0 && pos <= self.data.length() {
    self.data.insert(pos, ch)
  }
}

///|
pub fn Bp3MutableSeq::pop(self : Bp3MutableSeq, pos~ : Int32 = -1) -> Int32 {
  let idx = if pos < 0 { self.data.length() - 1 } else { pos }
  if idx >= 0 && idx < self.data.length() {
    self.data.remove(idx)
  } else {
    0
  }
}

///|
pub fn Bp3MutableSeq::remove(self : Bp3MutableSeq, pos : Int32) -> Unit {
  if pos >= 0 && pos < self.data.length() {
    ignore(self.data.remove(pos))
  }
}

///|
pub fn Bp3MutableSeq::reverse(self : Bp3MutableSeq) -> Unit {
  let n = self.data.length()
  let _half_n = n / 2
  for i in 0..<_half_n {
    let tmp = self.data[i]
    self.data[i] = self.data[n - 1 - i]
    self.data[n - 1 - i] = tmp
  }
}

///|
pub fn Bp3MutableSeq::complement(self : Bp3MutableSeq) -> Unit {
  let comp_map : Map[Int32, Int32] = Map([
    (65, 84),
    (84, 65),
    (71, 67),
    (67, 71),
    (97, 116),
    (116, 97),
    (103, 99),
    (99, 103),
  ])
  for i in 0..<self.data.length() {
    let v = comp_map.get(self.data[i])
    match v {
      Some(c) => self.data[i] = c
      None => ()
    }
  }
}

///|
pub fn Bp3MutableSeq::extend_seq(
  self : Bp3MutableSeq,
  other : Bp3MutableSeq,
) -> Unit {
  for i in 0..<other.data.length() {
    self.data.push(other.data[i])
  }
}

///|
pub fn Bp3MutableSeq::set(
  self : Bp3MutableSeq,
  pos : Int32,
  ch : Int32,
) -> Unit {
  if pos >= 0 && pos < self.data.length() {
    self.data[pos] = ch
  }
}

///|
pub fn Bp3MutableSeq::get(self : Bp3MutableSeq, pos : Int32) -> Int32 {
  if pos >= 0 && pos < self.data.length() {
    self.data[pos]
  } else {
    0
  }
}

///|
pub fn Bp3MutableSeq::delete_range(
  self : Bp3MutableSeq,
  start : Int32,
  end : Int32,
) -> Unit {
  if start < 0 || end > self.data.length() || start >= end {
    return
  }
  let count = end - start
  for _i in 0..<count {
    ignore(self.data.remove(start))
  }
}

///|
/// #2 UnknownSeq — unknown content but known length
pub struct Bp3UnknownSeq {
  length_ : Int32
  alphabet : String
  description : String
}

///|
pub fn Bp3UnknownSeq::new(
  length_ : Int32,
  alphabet~ : String = "ACGT",
  description~ : String = "",
) -> Bp3UnknownSeq {
  { length_, alphabet, description, }
}

///|
pub fn Bp3UnknownSeq::length(self : Bp3UnknownSeq) -> Int32 {
  self.length_
}

///|
pub fn Bp3UnknownSeq::to_string(self : Bp3UnknownSeq) -> String {
  let sb = StringBuilder()
  let fill = if self.alphabet.length() > 0 { self.alphabet[0] } else { 78 }
  for _i in 0..<self.length_ {
    let _ = sb.write_char(fill.unsafe_to_char())
  }
  sb.to_string()
}

///|
pub fn Bp3UnknownSeq::is_defined(_self : Bp3UnknownSeq) -> Bool {
  false
}

///|
pub fn Bp3UnknownSeq::count(_self : Bp3UnknownSeq, _sub : String) -> Int32 {
  0
}

///|
pub fn Bp3UnknownSeq::find(_self : Bp3UnknownSeq, _sub : String) -> Int32 {
  -1
}

///|
pub fn Bp3UnknownSeq::slice(
  self : Bp3UnknownSeq,
  start : Int32,
  end : Int32,
) -> Bp3UnknownSeq {
  let s = if start < 0 { 0 } else { start }
  let e = if end > self.length_ { self.length_ } else { end }
  let new_len = if e > s { e - s } else { 0 }
  { length_: new_len, alphabet: self.alphabet, description: self.description, }
}

///|
pub fn Bp3UnknownSeq::concat(
  self : Bp3UnknownSeq,
  other : Bp3UnknownSeq,
) -> Bp3UnknownSeq {
  {
    length_: self.length_ + other.length_,
    alphabet: self.alphabet,
    description: self.description,
  }
}

///|
pub fn _bp3_mutable_seq_dummy() -> Unit {
  let _ = Bp3MutableSeq::{ data: [], alphabet: "", }
  let _ = Bp3UnknownSeq::{ length_: 0, alphabet: "", description: "", }
}
