///|
pub struct BSGenome {
  name : String
  organism : String
  provider : String
  release : String
  chromosomes : Array[(String, String)]
  chromosome_lengths : Array[(String, Int)]
}

///|
pub fn BSGenome::new(
  name : String,
  organism : String,
  provider : String,
  release : String,
) -> BSGenome {
  BSGenome::{
    name,
    organism,
    provider,
    release,
    chromosomes: Array::new(),
    chromosome_lengths: Array::new(),
  }
}

///|
pub fn BSGenome::add_chromosome(
  self : BSGenome,
  name : String,
  sequence : String,
) -> BSGenome {
  let new_chromosomes = self.chromosomes.copy()
  new_chromosomes.push((name, sequence))
  let new_lengths = self.chromosome_lengths.copy()
  new_lengths.push((name, sequence.length()))
  BSGenome::{
    name: self.name,
    organism: self.organism,
    provider: self.provider,
    release: self.release,
    chromosomes: new_chromosomes,
    chromosome_lengths: new_lengths,
  }
}

///|
pub fn BSGenome::chromosome_names(self : BSGenome) -> Array[String] {
  let names : Array[String] = Array::new()
  let mut i = 0
  while i < self.chromosomes.length() {
    names.push(self.chromosomes[i].0)
    i = i + 1
  }
  names
}

///|
pub fn BSGenome::chromosome_length(self : BSGenome, chr : String) -> Int? {
  let mut i = 0
  while i < self.chromosome_lengths.length() {
    if self.chromosome_lengths[i].0 == chr {
      return Some(self.chromosome_lengths[i].1)
    }
    i = i + 1
  }
  None
}

///|
pub fn BSGenome::chromosome_sequence(self : BSGenome, chr : String) -> String? {
  let mut i = 0
  while i < self.chromosomes.length() {
    if self.chromosomes[i].0 == chr {
      return Some(self.chromosomes[i].1)
    }
    i = i + 1
  }
  None
}

///|
pub fn BSGenome::get_sequence(
  self : BSGenome,
  chr : String,
  start : Int,
  end : Int,
) -> String? {
  let seq = self.chromosome_sequence(chr)
  if seq is None {
    return None
  }
  let s = seq.unwrap()
  let len = s.length()
  if start < 0 || end > len || start >= end {
    return None
  }
  Some(s[start:end].to_owned())
}

///|
pub fn BSGenome::get_subseq(
  self : BSGenome,
  chr : String,
  start : Int,
  length : Int,
) -> String? {
  let seq = self.chromosome_sequence(chr)
  if seq is None {
    return None
  }
  let s = seq.unwrap()
  let len = s.length()
  if start < 0 || start + length > len {
    return None
  }
  Some(s[start:start + length].to_owned())
}

///|
pub fn BSGenome::total_length(self : BSGenome) -> Int {
  let mut total = 0
  let mut i = 0
  while i < self.chromosome_lengths.length() {
    total = total + self.chromosome_lengths[i].1
    i = i + 1
  }
  total
}

///|
pub fn BSGenome::chromosome_count(self : BSGenome) -> Int {
  self.chromosomes.length()
}

///|
pub fn BSGenome::get_chromosome_names(self : BSGenome) -> Array[String] {
  let names : Array[String] = Array::new()
  let mut i = 0
  while i < self.chromosomes.length() {
    names.push(self.chromosomes[i].0)
    i = i + 1
  }
  names
}

///|
pub fn BSGenome::has_chromosome(self : BSGenome, chr : String) -> Bool {
  let mut i = 0
  while i < self.chromosomes.length() {
    if self.chromosomes[i].0 == chr {
      return true
    }
    i = i + 1
  }
  false
}

///|
pub fn BSGenome::extract_gene(
  self : BSGenome,
  chr : String,
  start : Int,
  end : Int,
  strand : String,
) -> String? {
  let seq = self.get_sequence(chr, start, end)
  if seq is None {
    return None
  }
  let s = seq.unwrap()
  if strand == "-" {
    Some(rev_comp(s))
  } else {
    Some(s)
  }
}

///|
fn rev_comp(s : String) -> String {
  let n = s.length()
  let result : Array[String] = Array::new()
  let mut i = n - 1
  while i >= 0 {
    let c = s.unsafe_get(i).unsafe_to_char().to_string().to_upper()
    if c == "A" {
      result.push("T")
    } else if c == "T" {
      result.push("A")
    } else if c == "C" {
      result.push("G")
    } else if c == "G" {
      result.push("C")
    } else {
      result.push(c)
    }
    i = i - 1
  }
  result.join("")
}

///|
pub fn create_demo_genome() -> BSGenome {
  let genome = BSGenome::new("Homo sapiens", "Homo sapiens", "UCSC", "hg38")
  let genome = genome.add_chromosome(
    "chr1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGT",
  )
  let genome = genome.add_chromosome(
    "chr2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA",
  )
  let genome = genome.add_chromosome(
    "chr3", "ATCGATCGATCGATCGATCGATCGATCGATCGATC",
  )
  let genome = genome.add_chromosome(
    "chrX", "TGCAATGCAATGCAATGCAATGCAATGCAATGCAA",
  )
  genome.add_chromosome("chrY", "CGATCGATCGATCGATCGATCGATCGATCGATCGA")
}
