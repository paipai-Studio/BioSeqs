///|
/// InteractionSet - Chromatin interaction data inspired by Bioconductor InteractionSet.
/// Manages Hi-C interactions, anchor pairs, and interaction matrices.

///|
/// Anchor for chromatin interaction (genomic region).
pub struct Anchor {
  chrom : String
  start : Int
  end : Int
}

///|
/// Construct Anchor.
pub fn Anchor::new(chrom~ : String, start~ : Int, end~ : Int) -> Anchor {
  Anchor::{ chrom, start, end }
}

///|
/// Get anchor width.
pub fn Anchor::width(self : Anchor) -> Int {
  self.end - self.start + 1
}

///|
/// Check if two anchors overlap.
pub fn anchors_overlap(a1 : Anchor, a2 : Anchor) -> Bool {
  if a1.chrom != a2.chrom {
    return false
  }
  a1.start <= a2.end && a2.start <= a1.end
}

///|
/// A genomic interaction (pair of anchors).
pub struct GInteraction {
  anchor1 : Anchor
  anchor2 : Anchor
  count : Int
}

///|
/// Construct GInteraction.
pub fn GInteraction::new(
  anchor1~ : Anchor,
  anchor2~ : Anchor,
  count? : Int = 0,
) -> GInteraction {
  GInteraction::{ anchor1, anchor2, count }
}

///|
/// Get interaction distance (bp between anchor midpoints on same chrom).
pub fn GInteraction::distance(self : GInteraction) -> Int? {
  if self.anchor1.chrom != self.anchor2.chrom {
    return None
  }
  let mid1 = (self.anchor1.start + self.anchor1.end) / 2
  let mid2 = (self.anchor2.start + self.anchor2.end) / 2
  if mid1 > mid2 {
    Some(mid1 - mid2)
  } else {
    Some(mid2 - mid1)
  }
}

///|
/// Check if interaction is intra-chromosomal.
pub fn GInteraction::is_intra_chromosomal(self : GInteraction) -> Bool {
  self.anchor1.chrom == self.anchor2.chrom
}

///|
/// Check if interaction is inter-chromosomal.
pub fn GInteraction::is_inter_chromosomal(self : GInteraction) -> Bool {
  self.anchor1.chrom != self.anchor2.chrom
}

///|
/// InteractionSet - a collection of genomic interactions.
pub struct InteractionSet {
  interactions : Array[GInteraction]
  bin_size : Int
}

///|
/// Construct InteractionSet.
pub fn InteractionSet::new(
  interactions~ : Array[GInteraction],
  bin_size? : Int = 0,
) -> InteractionSet {
  InteractionSet::{ interactions, bin_size }
}

///|
/// Get number of interactions.
pub fn InteractionSet::num_interactions(self : InteractionSet) -> Int {
  self.interactions.length()
}

///|
/// Get all interactions.
pub fn InteractionSet::interactions(
  self : InteractionSet,
) -> Array[GInteraction] {
  self.interactions
}

///|
/// Get total interaction count.
pub fn InteractionSet::total_count(self : InteractionSet) -> Int {
  let mut total = 0
  for i in self.interactions {
    total = total + i.count
  }
  total
}

///|
/// Filter interactions by minimum count.
pub fn InteractionSet::filter_by_count(
  self : InteractionSet,
  min_count : Int,
) -> InteractionSet {
  let filtered = self.interactions.filter(fn(i) { i.count >= min_count })
  InteractionSet::new(interactions=filtered, bin_size=self.bin_size)
}

///|
/// Get intra-chromosomal interactions.
pub fn InteractionSet::intra_chromosomal(
  self : InteractionSet,
) -> InteractionSet {
  let filtered = self.interactions.filter(fn(i) { i.is_intra_chromosomal() })
  InteractionSet::new(interactions=filtered, bin_size=self.bin_size)
}

///|
/// Get inter-chromosomal interactions.
pub fn InteractionSet::inter_chromosomal(
  self : InteractionSet,
) -> InteractionSet {
  let filtered = self.interactions.filter(fn(i) { i.is_inter_chromosomal() })
  InteractionSet::new(interactions=filtered, bin_size=self.bin_size)
}

///|
/// Get interactions on a specific chromosome.
pub fn InteractionSet::interactions_on_chrom(
  self : InteractionSet,
  chrom : String,
) -> InteractionSet {
  let filtered = self.interactions.filter(fn(i) {
    i.anchor1.chrom == chrom || i.anchor2.chrom == chrom
  })
  InteractionSet::new(interactions=filtered, bin_size=self.bin_size)
}

///|
/// Build an interaction matrix (chrom1 x chrom2).
pub fn InteractionSet::interaction_matrix(
  self : InteractionSet,
  chroms : Array[String],
) -> Array[Array[Int]] {
  let n = chroms.length()
  let matrix : Array[Array[Int]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Int] = Array::make(n, 0)
    matrix.push(row)
    i = i + 1
  }
  for interaction in self.interactions {
    let mut idx1 = -1
    let mut idx2 = -1
    let mut j = 0
    while j < n {
      if chroms[j] == interaction.anchor1.chrom {
        idx1 = j
      }
      if chroms[j] == interaction.anchor2.chrom {
        idx2 = j
      }
      j = j + 1
    }
    if idx1 >= 0 && idx2 >= 0 {
      matrix[idx1][idx2] = matrix[idx1][idx2] + interaction.count
      if idx1 != idx2 {
        matrix[idx2][idx1] = matrix[idx2][idx1] + interaction.count
      }
    }
  }
  matrix
}

///|
/// Get distance distribution (histogram of interaction distances).
pub fn InteractionSet::distance_distribution(
  self : InteractionSet,
  bin_width : Int,
) -> Array[(Int, Int)] {
  let result : Array[(Int, Int)] = Array::new()
  let distances : Array[Int] = Array::new()
  for i in self.interactions {
    match i.distance() {
      Some(d) => distances.push(d)
      None => ()
    }
  }
  if distances.length() == 0 {
    return result
  }
  let mut max_dist = 0
  for d in distances {
    if d > max_dist {
      max_dist = d
    }
  }
  let num_bins = max_dist / bin_width + 1
  let counts : Array[Int] = Array::make(num_bins, 0)
  for d in distances {
    let bin = d / bin_width
    counts[bin] = counts[bin] + 1
  }
  let mut b = 0
  while b < num_bins {
    if counts[b] > 0 {
      result.push((b * bin_width, counts[b]))
    }
    b = b + 1
  }
  result
}

///|
/// Get top N interactions by count.
pub fn InteractionSet::top_interactions(
  self : InteractionSet,
  n : Int,
) -> Array[GInteraction] {
  let sorted = self.interactions
  let len = sorted.length()
  let mut i = 0
  while i < len {
    let mut j = 0
    while j < len - i - 1 {
      if sorted[j].count < sorted[j + 1].count {
        let tmp = sorted[j]
        sorted[j] = sorted[j + 1]
        sorted[j + 1] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  let result : Array[GInteraction] = Array::new()
  let mut k = 0
  while k < n && k < sorted.length() {
    result.push(sorted[k])
    k = k + 1
  }
  result
}

///|
/// Calculate the fraction of inter-chromosomal interactions.
pub fn InteractionSet::inter_fraction(self : InteractionSet) -> Double {
  if self.interactions.length() == 0 {
    return 0.0
  }
  let inter = self.interactions
    .filter(fn(i) { i.is_inter_chromosomal() })
    .length()
  inter.to_double() / self.interactions.length().to_double()
}

///|
/// Get anchor regions from interactions.
pub fn InteractionSet::anchors(self : InteractionSet) -> Array[Anchor] {
  let result : Array[Anchor] = Array::new()
  for i in self.interactions {
    result.push(i.anchor1)
    result.push(i.anchor2)
  }
  result
}

///|
/// Create a simple Hi-C interaction set from bin data.
pub fn create_binned_interactions(
  chrom : String,
  bins : Array[(Int, Int)],
  counts : Array[Array[Int]],
) -> InteractionSet {
  let interactions : Array[GInteraction] = Array::new()
  let n = bins.length()
  let mut i = 0
  while i < n {
    let mut j = i
    while j < n {
      let c = counts[i][j]
      if c > 0 {
        let a1 = Anchor::new(chrom~, start=bins[i].0, end=bins[i].1)
        let a2 = Anchor::new(chrom~, start=bins[j].0, end=bins[j].1)
        interactions.push(GInteraction::new(anchor1=a1, anchor2=a2, count=c))
      }
      j = j + 1
    }
    i = i + 1
  }
  InteractionSet::new(interactions~, bin_size=0)
}

///|
/// Get interaction summary.
pub fn InteractionSet::summary(self : InteractionSet) -> String {
  let n = self.num_interactions()
  let total = self.total_count()
  let intra = self.intra_chromosomal().num_interactions()
  let inter = self.inter_chromosomal().num_interactions()
  "InteractionSet: " +
  n.to_string() +
  " interactions, " +
  total.to_string() +
  " total counts, " +
  intra.to_string() +
  " intra-chr, " +
  inter.to_string() +
  " inter-chr"
}
