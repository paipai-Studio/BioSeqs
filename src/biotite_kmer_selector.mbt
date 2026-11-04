/// biotite_kmer_selector.mbt — K-mer subset selection
/// Inspired by biotite.sequence.align MinimizerSelector, SyncmerSelector, MincodeSelector

///|
/// Minimizer selector: pick k-mers with minimum hash value in each window.
pub struct BiotiteMinimizerSelector {
  k : Int32 // k-mer size
  window_size : Int32 // window size (in k-mers)
}

///|
/// Create a minimizer selector.
pub fn BiotiteMinimizerSelector::new(
  k~ : Int32 = 15,
  window_size~ : Int32 = 10,
) -> BiotiteMinimizerSelector {
  { k, window_size, }
}

///|
/// Hash a k-mer string to an integer (for minimizer comparison).
fn biotite_kmer_hash_str(s : String) -> Int32 {
  let mut h = 0
  for i = 0; i < s.length(); i = i + 1 {
    h = h * 31 + s[i].to_int()
  }
  if h < 0 {
    -h
  } else {
    h
  }
}

///|
/// Select minimizers from a sequence.
/// Returns positions of selected k-mers.
pub fn BiotiteMinimizerSelector::select(
  self : BiotiteMinimizerSelector,
  sequence : String,
) -> Array[Int32] {
  let n = sequence.length()
  if n < self.k {
    return []
  }
  let n_kmers = n - self.k + 1
  if n_kmers < self.window_size {
    // Select the minimum from all k-mers
    let mut min_pos = 0
    let mut min_hash = biotite_kmer_hash_str(sequence[0:self.k].to_owned())
    for i = 1; i < n_kmers; i = i + 1 {
      let h = biotite_kmer_hash_str(sequence[i:i + self.k].to_owned())
      if h < min_hash {
        min_hash = h
        min_pos = i
      }
    }
    return [min_pos]
  }
  let selected : Array[Int32] = Array(capacity=0)
  let mut last_selected = -1
  for w_start = 0; w_start <= n_kmers - self.window_size; w_start = w_start + 1 {
    let w_end = w_start + self.window_size
    let mut min_pos = w_start
    let mut min_hash = biotite_kmer_hash_str(
      sequence[w_start:w_start + self.k].to_owned(),
    )
    for i = w_start + 1; i < w_end; i = i + 1 {
      let h = biotite_kmer_hash_str(sequence[i:i + self.k].to_owned())
      if h < min_hash {
        min_hash = h
        min_pos = i
      }
    }
    if min_pos != last_selected {
      selected.push(min_pos)
      last_selected = min_pos
    }
  }
  selected
}

///|
/// Syncmer selector: pick k-mers whose internal minimizer is at a specific position.
/// A syncmer is a k-mer where the smallest s-mer occurs at a designated position.
pub struct BiotiteSyncmerSelector {
  k : Int32 // k-mer size
  s : Int32 // s-mer size (internal minimizer)
  t : Int32 // target position of minimizer within k-mer (0-indexed)
}

///|
/// Create a syncmer selector.
pub fn BiotiteSyncmerSelector::new(
  k~ : Int32 = 15,
  s~ : Int32 = 10,
  t~ : Int32 = 0,
) -> BiotiteSyncmerSelector {
  { k, s, t, }
}

///|
/// Select syncmers from a sequence.
pub fn BiotiteSyncmerSelector::select(
  self : BiotiteSyncmerSelector,
  sequence : String,
) -> Array[Int32] {
  let n = sequence.length()
  if n < self.k {
    return []
  }
  let selected : Array[Int32] = Array(capacity=0)
  for pos = 0; pos <= n - self.k; pos = pos + 1 {
    // Find the position of the minimum s-mer within this k-mer
    let kmer = sequence[pos:pos + self.k].to_owned()
    if kmer.length() < self.s {
      continue
    }
    let mut min_s_pos = 0
    let mut min_s_hash = biotite_kmer_hash_str(kmer[0:self.s].to_owned())
    for j = 1; j <= self.k - self.s; j = j + 1 {
      let h = biotite_kmer_hash_str(kmer[j:j + self.s].to_owned())
      if h < min_s_hash {
        min_s_hash = h
        min_s_pos = j
      }
    }
    if min_s_pos == self.t {
      selected.push(pos)
    }
  }
  selected
}

///|
/// Mincode selector: select k-mers whose hash value is below a threshold.
/// This gives a fixed-rate sampling of k-mers.
pub struct BiotiteMincodeSelector {
  k : Int32 // k-mer size
  threshold : Int32 // hash threshold (0 to hash_max)
  hash_max : Int32 // maximum hash value
}

///|
/// Create a mincode selector with a target density.
pub fn BiotiteMincodeSelector::new(
  k~ : Int32 = 15,
  density~ : Double = 0.1,
) -> BiotiteMincodeSelector {
  let hash_max = 2147483647 // Int max
  let threshold = (hash_max.to_double() * density).to_int()
  { k, threshold, hash_max, }
}

///|
/// Select mincode k-mers from a sequence.
pub fn BiotiteMincodeSelector::select(
  self : BiotiteMincodeSelector,
  sequence : String,
) -> Array[Int32] {
  let n = sequence.length()
  if n < self.k {
    return []
  }
  let selected : Array[Int32] = Array(capacity=0)
  for pos = 0; pos <= n - self.k; pos = pos + 1 {
    let kmer = sequence[pos:pos + self.k].to_owned()
    let h = biotite_kmer_hash_str(kmer)
    if h <= self.threshold {
      selected.push(pos)
    }
  }
  selected
}

///|
/// Get the expected density of selected k-mers.
pub fn BiotiteMincodeSelector::density(self : BiotiteMincodeSelector) -> Double {
  self.threshold.to_double() / self.hash_max.to_double()
}
