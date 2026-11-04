/// Named constants for BAM FLAG field bits.
///
/// pysam equivalent: `pysam.FPAIRED`, `pysam.FPROPER_PAIR`, etc.

///|
/// BAM flag bit constants
pub const BAM_FPAIRED : Int32 = 0x1

///|
pub const BAM_FPROPER_PAIR : Int32 = 0x2

///|
pub const BAM_FUNMAP : Int32 = 0x4

///|
pub const BAM_FMUNMAP : Int32 = 0x8

///|
pub const BAM_FREVERSE : Int32 = 0x10

///|
pub const BAM_FMREVERSE : Int32 = 0x20

///|
pub const BAM_FREAD1 : Int32 = 0x40

///|
pub const BAM_FREAD2 : Int32 = 0x80

///|
pub const BAM_FSECONDARY : Int32 = 0x100

///|
pub const BAM_FQCFAIL : Int32 = 0x200

///|
pub const BAM_FDUP : Int32 = 0x400

///|
pub const BAM_FSUPPLEMENTARY : Int32 = 0x800

///|
/// Check if a specific flag bit is set
pub fn bam_flag_has(flag : Int32, bit : Int32) -> Bool {
  flag.land(bit) != 0
}

///|
/// Set a flag bit (returns new flag value)
pub fn bam_flag_set(flag : Int32, bit : Int32) -> Int32 {
  flag.lor(bit)
}

///|
/// Clear a flag bit (returns new flag value)
pub fn bam_flag_clear(flag : Int32, bit : Int32) -> Int32 {
  flag.land(bit.lnot())
}

///|
/// Toggle a flag bit (returns new flag value)
pub fn bam_flag_toggle(flag : Int32, bit : Int32) -> Int32 {
  flag.lxor(bit)
}

///|
/// Convert a flag integer to a human-readable description string.
pub fn bam_flag_to_string(flag : Int32) -> String {
  let parts : Array[String] = Array(capacity=0)
  if bam_flag_has(flag, BAM_FPAIRED) {
    parts.push("paired")
  }
  if bam_flag_has(flag, BAM_FPROPER_PAIR) {
    parts.push("proper_pair")
  }
  if bam_flag_has(flag, BAM_FUNMAP) {
    parts.push("unmapped")
  }
  if bam_flag_has(flag, BAM_FMUNMAP) {
    parts.push("mate_unmapped")
  }
  if bam_flag_has(flag, BAM_FREVERSE) {
    parts.push("reverse")
  }
  if bam_flag_has(flag, BAM_FMREVERSE) {
    parts.push("mate_reverse")
  }
  if bam_flag_has(flag, BAM_FREAD1) {
    parts.push("read1")
  }
  if bam_flag_has(flag, BAM_FREAD2) {
    parts.push("read2")
  }
  if bam_flag_has(flag, BAM_FSECONDARY) {
    parts.push("secondary")
  }
  if bam_flag_has(flag, BAM_FQCFAIL) {
    parts.push("qcfail")
  }
  if bam_flag_has(flag, BAM_FDUP) {
    parts.push("duplicate")
  }
  if bam_flag_has(flag, BAM_FSUPPLEMENTARY) {
    parts.push("supplementary")
  }
  parts.join(", ")
}

///|
/// Parse a flag string like "paired,reverse" into a flag integer.
pub fn bam_flag_from_string(s : String) -> Int32 {
  let mut flag : Int32 = 0
  let parts = s.split(",").to_array()
  for part in parts {
    let p = part.trim().to_lower()
    match p {
      "paired" => flag = bam_flag_set(flag, BAM_FPAIRED)
      "proper_pair" => flag = bam_flag_set(flag, BAM_FPROPER_PAIR)
      "unmapped" => flag = bam_flag_set(flag, BAM_FUNMAP)
      "mate_unmapped" => flag = bam_flag_set(flag, BAM_FMUNMAP)
      "reverse" => flag = bam_flag_set(flag, BAM_FREVERSE)
      "mate_reverse" => flag = bam_flag_set(flag, BAM_FMREVERSE)
      "read1" => flag = bam_flag_set(flag, BAM_FREAD1)
      "read2" => flag = bam_flag_set(flag, BAM_FREAD2)
      "secondary" => flag = bam_flag_set(flag, BAM_FSECONDARY)
      "qcfail" => flag = bam_flag_set(flag, BAM_FQCFAIL)
      "duplicate" => flag = bam_flag_set(flag, BAM_FDUP)
      "supplementary" => flag = bam_flag_set(flag, BAM_FSUPPLEMENTARY)
      _ => ()
    }
  }
  flag
}
