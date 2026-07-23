///|
/// Bio.Kmer - k-mer counting and analysis.

///|
pub struct BioKmerCounter {
  k : Int
  counts : Map[String, Int]
  sequence : String
}

///|
pub fn BioKmerCounter::new(k : Int) -> BioKmerCounter {
  BioKmerCounter::{
    k,
    counts: Map([], capacity=1000),
    sequence: ""
  }
}

///|
pub fn BioKmerCounter::count_kmers(self : BioKmerCounter, seq : String) -> BioKmerCounter {
  let n = seq.length() - self.k + 1
  let counts = self.counts.copy()
  
  let mut i = 0
  while i < n {
    let kmer = substring(seq, i, self.k)
    if counts.contains(kmer) {
      counts.set(kmer, counts.get(kmer).unwrap() + 1)
    } else {
      counts.set(kmer, 1)
    }
    i = i + 1
  }
  
  BioKmerCounter::{ k: self.k, counts: counts, sequence: seq }
}

///|
pub fn BioKmerCounter::count(self : BioKmerCounter, kmer : String) -> Int {
  if self.counts.contains(kmer) {
    self.counts.get(kmer).unwrap()
  } else {
    0
  }
}

///|
pub fn BioKmerCounter::total_kmers(self : BioKmerCounter) -> Int {
  let mut total = 0
  for count in self.counts.values() {
    total = total + count
  }
  total
}

///|
pub fn BioKmerCounter::unique_kmers(self : BioKmerCounter) -> Int {
  self.counts.keys().length()
}

///|
pub fn BioKmerCounter::frequency(self : BioKmerCounter, kmer : String) -> Double {
  let total = self.total_kmers()
  if total == 0 {
    0.0
  } else {
    self.count(kmer).to_double() / total.to_double()
  }
}

///|
pub fn BioKmerCounter::get_top_kmers(self : BioKmerCounter, n : Int) -> Array[(String, Int)] {
  let items : Array[(String, Int)] = Array::new()
  for key in self.counts.keys() {
    items.push((key, self.counts.get(key).unwrap()))
  }
  
  let mut swapped = true
  let mut j = 0
  while swapped {
    swapped = false
    j = j + 1
    let mut i = 0
    while i < items.length() - j {
      if items[i].1 < items[i + 1].1 {
        let temp = items[i]
        items[i] = items[i + 1]
        items[i + 1] = temp
        swapped = true
      }
      i = i + 1
    }
  }
  
  let result : Array[(String, Int)] = Array::new()
  let mut k = 0
  while k < n && k < items.length() {
    result.push(items[k])
    k = k + 1
  }
  
  result
}

///|
pub fn BioKmerCounter::jaccard(self : BioKmerCounter, other : BioKmerCounter) -> Double {
  let mut intersection = 0
  let mut union = 0
  
  for key in self.counts.keys() {
    if other.counts.contains(key) {
      intersection = intersection + 1
    }
    union = union + 1
  }
  
  for key in other.counts.keys() {
    if !self.counts.contains(key) {
      union = union + 1
    }
  }
  
  if union == 0 {
    0.0
  } else {
    intersection.to_double() / union.to_double()
  }
}

///|
pub fn BioKmerCounter::hamming_distance(self : BioKmerCounter, other : BioKmerCounter) -> Int {
  let mut distance = 0
  
  for key in self.counts.keys() {
    let diff = self.counts.get(key).unwrap() - other.count(key)
    let abs_diff = if diff > 0 { diff } else { -diff }
    distance = distance + abs_diff
  }
  
  for key in other.counts.keys() {
    if !self.counts.contains(key) {
      distance = distance + other.count(key)
    }
  }
  
  distance
}

///|
pub struct KmerFrequency {
  kmer : String
  count : Int
  frequency : Double
}

///|
pub fn KmerFrequency::new(kmer : String, count : Int, frequency : Double) -> KmerFrequency {
  KmerFrequency::{ kmer, count, frequency }
}

///|
pub fn bio_count_kmers(seq : String, k : Int) -> Map[String, Int] {
  let counter = BioKmerCounter::new(k)
  let counter = counter.count_kmers(seq)
  counter.counts
}

///|
pub fn kmer_spectrum(seq : String, min_k : Int, max_k : Int) -> Array[Int] {
  let spectrum : Array[Int] = Array::new()
  
  let mut k = min_k
  while k <= max_k {
    let counter = BioKmerCounter::new(k)
    let counter = counter.count_kmers(seq)
    spectrum.push(counter.unique_kmers())
    k = k + 1
  }
  
  spectrum
}

///|
pub fn create_example_kmer_counter() -> BioKmerCounter {
  let counter = BioKmerCounter::new(3)
  let counter = counter.count_kmers("ATCGATCGATCG")
  counter
}