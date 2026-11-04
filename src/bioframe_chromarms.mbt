/// Make chromosome arms from centromere positions.
///
/// Bioframe equivalent: `make_chromarms`
///
/// Given centromere positions for each chromosome, creates two arms per
/// chromosome: p-arm (0..centromere) and q-arm (centromere..end).

///|
/// A chromosome arm
pub struct BfChromArm {
  name : String
  chrom : String
  start : Int32
  end : Int32
  arm : Char
}

///|
pub fn BfChromArm::name(self : BfChromArm) -> String {
  self.name
}

///|
pub fn BfChromArm::chrom(self : BfChromArm) -> String {
  self.chrom
}

///|
pub fn BfChromArm::start(self : BfChromArm) -> Int32 {
  self.start
}

///|
pub fn BfChromArm::end(self : BfChromArm) -> Int32 {
  self.end
}

///|
pub fn BfChromArm::arm(self : BfChromArm) -> Char {
  self.arm
}

///|
/// Make chromosome arms from centromere positions.
///
/// chromsizes: array of (chrom, length)
/// centromeres: array of (chrom, position) where position is the centromere midpoint
///
/// Returns array of ChromArm, each with name "chrN_p" or "chrN_q".
pub fn bf_make_chromarms(
  chromsizes : Array[(String, Int32)],
  centromeres : Array[(String, Int32)],
) -> Array[BfChromArm] {
  let cent_map : Map[String, Int32] = Map([])
  for c in centromeres {
    cent_map[c.0] = c.1
  }
  let result : Array[BfChromArm] = Array(capacity=0)
  for cs in chromsizes {
    let chrom = cs.0
    let length = cs.1
    match cent_map.get(chrom) {
      Some(cent_pos) => {
        // p-arm: 0..centromere
        if cent_pos > 0 {
          result.push({
            name: chrom + "_p",
            chrom,
            start: 0,
            end: cent_pos,
            arm: 'p',
          })
        }
        // q-arm: centromere..end
        if cent_pos < length {
          result.push({
            name: chrom + "_q",
            chrom,
            start: cent_pos,
            end: length,
            arm: 'q',
          })
        }
      }
      None =>
        // No centromere: whole chromosome as one arm
        result.push({ name: chrom, chrom, start: 0, end: length, arm: 'n', })
    }
  }
  result
}

///|
/// Make chromarms and return as a viewframe.
pub fn bf_make_chromarms_view(
  chromsizes : Array[(String, Int32)],
  centromeres : Array[(String, Int32)],
) -> BfViewFrame {
  let arms = bf_make_chromarms(chromsizes, centromeres)
  let regions : Array[BfViewRegion] = Array(capacity=0)
  for arm in arms {
    regions.push({
      name: arm.name,
      chrom: arm.chrom,
      start: arm.start,
      end: arm.end,
    })
  }
  { regions, }
}

///|
/// Get only p-arms from chromarms.
pub fn bf_get_p_arms(arms : Array[BfChromArm]) -> Array[BfChromArm] {
  let result : Array[BfChromArm] = Array(capacity=0)
  for arm in arms {
    if arm.arm == 'p' {
      result.push(arm)
    }
  }
  result
}

///|
/// Get only q-arms from chromarms.
pub fn bf_get_q_arms(arms : Array[BfChromArm]) -> Array[BfChromArm] {
  let result : Array[BfChromArm] = Array(capacity=0)
  for arm in arms {
    if arm.arm == 'q' {
      result.push(arm)
    }
  }
  result
}
