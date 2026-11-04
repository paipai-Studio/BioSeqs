/// Bin genome into views with labels.
///
/// Bioframe equivalent: `bin_genome` with view support
///
/// Divide a genome (or view) into fixed-size bins, optionally labeled.

///|
/// A genome bin
pub struct BfBin {
  chrom : String
  start : Int32
  end : Int32
  bin_id : Int32
  label : String
}

///|
pub fn BfBin::chrom(self : BfBin) -> String {
  self.chrom
}

///|
pub fn BfBin::start(self : BfBin) -> Int32 {
  self.start
}

///|
pub fn BfBin::end(self : BfBin) -> Int32 {
  self.end
}

///|
pub fn BfBin::bin_id(self : BfBin) -> Int32 {
  self.bin_id
}

///|
pub fn BfBin::label(self : BfBin) -> String {
  self.label
}

///|
/// Bin a genome into fixed-size bins.
///
/// chromsizes: array of (chrom, length)
/// bin_size: width of each bin
/// label: optional label for all bins (default "bin")
///
/// The last bin per chromosome may be shorter than bin_size.
pub fn bf_bin_genome(
  chromsizes : Array[(String, Int32)],
  bin_size : Int32,
  label~ : String = "bin",
) -> Array[BfBin] {
  let result : Array[BfBin] = Array(capacity=0)
  let mut bin_id : Int32 = 0
  for cs in chromsizes {
    let chrom = cs.0
    let length = cs.1
    let mut pos : Int32 = 0
    while pos < length {
      let end = if pos + bin_size > length { length } else { pos + bin_size }
      result.push({ chrom, start: pos, end, bin_id, label, })
      bin_id = bin_id + 1
      pos = pos + bin_size
    }
  }
  result
}

///|
/// Bin a specific view region into fixed-size bins.
pub fn bf_bin_view(view : BfViewFrame, bin_size : Int32) -> Array[BfBin] {
  let result : Array[BfBin] = Array(capacity=0)
  let mut bin_id : Int32 = 0
  for vr in view.regions {
    let mut pos = vr.start
    while pos < vr.end {
      let end = if pos + bin_size > vr.end { vr.end } else { pos + bin_size }
      result.push({ chrom: vr.chrom, start: pos, end, bin_id, label: vr.name, })
      bin_id = bin_id + 1
      pos = pos + bin_size
    }
  }
  result
}

///|
/// Bin genome with step (sliding window bins).
/// Overlapping bins with given step size.
pub fn bf_bin_genome_stepped(
  chromsizes : Array[(String, Int32)],
  bin_size : Int32,
  step : Int32,
  label~ : String = "bin",
) -> Array[BfBin] {
  let result : Array[BfBin] = Array(capacity=0)
  let mut bin_id : Int32 = 0
  for cs in chromsizes {
    let chrom = cs.0
    let length = cs.1
    let mut pos : Int32 = 0
    while pos + bin_size <= length {
      result.push({ chrom, start: pos, end: pos + bin_size, bin_id, label, })
      bin_id = bin_id + 1
      pos = pos + step
    }
  }
  result
}

///|
/// Convert bins to BedInterval array (dropping metadata).
pub fn bf_bins_to_intervals(bins : Array[BfBin]) -> Array[BedInterval] {
  let result : Array[BedInterval] = Array(capacity=0)
  for bin in bins {
    result.push(BedInterval::{
      chrom: bin.chrom,
      start: bin.start,
      end: bin.end,
    })
  }
  result
}

///|
/// Select bins by label.
pub fn bf_select_bins_by_label(
  bins : Array[BfBin],
  label : String,
) -> Array[BfBin] {
  let result : Array[BfBin] = Array(capacity=0)
  for bin in bins {
    if bin.label == label {
      result.push(bin)
    }
  }
  result
}

///|
/// Count bins per label.
pub fn bf_count_bins_by_label(bins : Array[BfBin]) -> Map[String, Int32] {
  let counts : Map[String, Int32] = Map([])
  for bin in bins {
    match counts.get(bin.label) {
      Some(c) => counts[bin.label] = c + 1
      None => counts[bin.label] = 1
    }
  }
  counts
}
