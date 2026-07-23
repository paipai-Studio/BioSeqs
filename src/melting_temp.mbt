///|
/// Bio.SeqUtils.MeltingTemp - DNA melting temperature calculation.

///|
pub struct MeltingTempResult {
  tm : Double
  method : String
}

///|
pub fn MeltingTempResult::new(tm : Double, method : String) -> MeltingTempResult {
  MeltingTempResult::{ tm, method }
}

///|
pub fn mt_wallace(seq : String) -> MeltingTempResult {
  let seq_upper = seq.to_upper()
  
  let mut gc_count = 0
  let mut at_count = 0
  
  let mut idx = 0
  while idx < seq_upper.length() {
    let c = substring(seq_upper, idx, 1)
    if c == "G" || c == "C" {
      gc_count = gc_count + 1
    } else if c == "A" || c == "T" {
      at_count = at_count + 1
    }
    idx = idx + 1
  }
  
  let tm = 4.0 * gc_count.to_double() + 2.0 * at_count.to_double()
  
  MeltingTempResult::new(tm, "wallace")
}

///|
pub fn mt_gc_content(seq : String) -> MeltingTempResult {
  let seq_upper = seq.to_upper()
  let length = seq_upper.length().to_double()
  
  let mut gc_count = 0
  
  let mut idx = 0
  while idx < seq_upper.length() {
    let c = substring(seq_upper, idx, 1)
    if c == "G" || c == "C" {
      gc_count = gc_count + 1
    }
    idx = idx + 1
  }
  
  let tm = 64.9 + 41.0 * (gc_count.to_double() - 16.4) / length
  
  MeltingTempResult::new(tm, "gc_content")
}

///|
pub fn create_example_dna() -> String {
  "AGCTGATCGATCGATCGATCGATCGATCGATCGATCG"
}

///|
pub fn recommend_tm_method(seq : String) -> String {
  let length = seq.length()
  
  if length <= 14 {
    "wallace"
  } else {
    "gc_content"
  }
}