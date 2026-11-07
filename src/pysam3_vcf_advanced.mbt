// Pysam third-round: VCF/BCF advanced operations (items 10-15)
// VCF writing, normalization, Ti/Tv, genotype concordance, heterozygosity, BCF

///|
pub struct Pysam3VcfHeader {
  version : String
  contigs : Array[(String, Int32)] // (contig_name, length)
  samples : Array[String]
  info_fields : Array[String]
  format_fields : Array[String]
}

///|
pub struct Pysam3VcfRecord {
  chrom : String
  pos : Int32 // 1-based
  id_ : String
  ref_ : String
  alt : String
  qual : Double
  filter_ : String
  info : String
  format_ : String
  genotypes : Array[String]
}

///|
pub fn pysam3_vcf_header_to_text(header : Pysam3VcfHeader) -> String {
  let sb = StringBuilder()
  let _ = sb.write("##fileformat=VCFv4.2\n")
  for i in 0..<header.contigs.length() {
    let (name, length) = header.contigs[i]
    let _ = sb.write("##contig=<ID=")
    let _ = sb.write(name)
    let _ = sb.write(",length=")
    let _ = sb.write(length.to_string())
    let _ = sb.write(">\n")
  }
  for i in 0..<header.info_fields.length() {
    let _ = sb.write("##INFO=")
    let _ = sb.write(header.info_fields[i])
    let _ = sb.write("\n")
  }
  for i in 0..<header.format_fields.length() {
    let _ = sb.write("##FORMAT=")
    let _ = sb.write(header.format_fields[i])
    let _ = sb.write("\n")
  }
  let _ = sb.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
  if header.samples.length() > 0 {
    let _ = sb.write("\tFORMAT")
    for i in 0..<header.samples.length() {
      let _ = sb.write("\t")
      let _ = sb.write(header.samples[i])
    }
  }
  let _ = sb.write("\n")
  sb.to_string()
}

///|
pub fn pysam3_vcf_record_to_text(rec : Pysam3VcfRecord) -> String {
  let sb = StringBuilder()
  let _ = sb.write(rec.chrom)
  let _ = sb.write("\t")
  let _ = sb.write(rec.pos.to_string())
  let _ = sb.write("\t")
  let _ = sb.write(rec.id_)
  let _ = sb.write("\t")
  let _ = sb.write(rec.ref_)
  let _ = sb.write("\t")
  let _ = sb.write(rec.alt)
  let _ = sb.write("\t")
  if rec.qual < 0.0 {
    let _ = sb.write(".")
  } else {
    let _ = sb.write(rec.qual.to_string())
  }
  let _ = sb.write("\t")
  let _ = sb.write(rec.filter_)
  let _ = sb.write("\t")
  let _ = sb.write(rec.info)
  if rec.format_.length() > 0 {
    let _ = sb.write("\t")
    let _ = sb.write(rec.format_)
    for i in 0..<rec.genotypes.length() {
      let _ = sb.write("\t")
      let _ = sb.write(rec.genotypes[i])
    }
  }
  sb.to_string()
}

///|
pub fn pysam3_vcf_write(
  header : Pysam3VcfHeader,
  records : Array[Pysam3VcfRecord],
) -> String {
  let sb = StringBuilder()
  let _ = sb.write(pysam3_vcf_header_to_text(header))
  for i in 0..<records.length() {
    let _ = sb.write(pysam3_vcf_record_to_text(records[i]))
    let _ = sb.write("\n")
  }
  sb.to_string()
}

// --- VCF normalization ---

///|
pub fn pysam3_vcf_left_align(
  rec : Pysam3VcfRecord,
  ref_seq : String,
) -> Pysam3VcfRecord {
  // Left-align indel: shift position left while ref/alt prefix matches
  let ref_ = rec.ref_
  let alt = rec.alt
  if ref_.length() <= 1 || alt.length() <= 1 {
    return rec
  }
  let mut pos = rec.pos
  let mut r = ref_
  let mut a = alt
  let mut shifted = true
  while shifted && pos > 1 {
    shifted = false
    let ref_idx = pos - 2
    if ref_idx >= 0 && ref_idx < ref_seq.length() {
      let prev_base = ref_seq[ref_idx].to_int()
      let r_last = r[r.length() - 1].to_int()
      let a_last = a[a.length() - 1].to_int()
      if r_last == prev_base && a_last == prev_base {
        // Can shift left
        let _sb_r = StringBuilder()
        let _ = _sb_r.write_char(Int32::unsafe_to_char(prev_base))
        let _ = _sb_r.write(r[0:r.length() - 1].to_owned())
        let new_r = _sb_r.to_string()
        let _sb_a = StringBuilder()
        let _ = _sb_a.write_char(Int32::unsafe_to_char(prev_base))
        let _ = _sb_a.write(a[0:a.length() - 1].to_owned())
        let new_a = _sb_a.to_string()
        r = new_r
        a = new_a
        pos = pos - 1
        shifted = true
      }
    }
  }
  { ..rec, pos, ref_: r, alt: a, }
}

///|
pub fn pysam3_vcf_split_multiallelic(
  rec : Pysam3VcfRecord,
) -> Array[Pysam3VcfRecord] {
  let alts = rec.alt.split(",").to_array()
  if alts.length() <= 1 {
    return [rec]
  }
  let result : Array[Pysam3VcfRecord] = []
  for i in 0..<alts.length() {
    result.push({ ..rec, alt: alts[i].to_owned(), })
  }
  result
}

// --- Ti/Tv ratio ---

///|
pub fn pysam3_is_transition(ref_ : String, alt : String) -> Bool {
  if ref_.length() != 1 || alt.length() != 1 {
    return false
  }
  let r = ref_[0].to_int()
  let a = alt[0].to_int()
  let is_purine = fn(b : Int32) -> Bool {
    b == 65 || b == 71 || b == 97 || b == 103
  }
  let is_pyrimidine = fn(b : Int32) -> Bool {
    b == 67 || b == 84 || b == 99 || b == 116
  }
  if is_purine(r) && is_purine(a) {
    true
  } else if is_pyrimidine(r) && is_pyrimidine(a) {
    true
  } else {
    false
  }
}

///|
pub fn pysam3_ti_tv_count(records : Array[Pysam3VcfRecord]) -> (Int32, Int32) {
  let mut ti = 0
  let mut tv = 0
  for i in 0..<records.length() {
    let rec = records[i]
    if rec.ref_.length() == 1 && rec.alt.length() == 1 {
      if pysam3_is_transition(rec.ref_, rec.alt) {
        ti = ti + 1
      } else {
        tv = tv + 1
      }
    }
  }
  (ti, tv)
}

///|
pub fn pysam3_ti_tv_ratio(records : Array[Pysam3VcfRecord]) -> Double {
  let (ti, tv) = pysam3_ti_tv_count(records)
  if tv == 0 {
    return 0.0
  }
  ti.to_double() / tv.to_double()
}

// --- Genotype concordance ---

///|
pub fn pysam3_genotype_concordance(gt1 : String, gt2 : String) -> Bool {
  // Compare genotypes, handling phase
  let normalize = fn(gt : String) -> String {
    let parts = gt.split("/").to_array()
    if parts.length() == 1 {
      let parts2 = gt.split("|").to_array()
      if parts2.length() > 1 {
        let arr : Array[String] = []
        for j in 0..<parts2.length() {
          arr.push(parts2[j].to_owned())
        }
        arr.sort()
        let sb = StringBuilder()
        for j in 0..<arr.length() {
          if j > 0 {
            let _ = sb.write("/")
          }
          let _ = sb.write(arr[j])
        }
        sb.to_string()
      } else {
        gt
      }
    } else {
      let arr : Array[String] = []
      for j in 0..<parts.length() {
        arr.push(parts[j].to_owned())
      }
      arr.sort()
      let sb = StringBuilder()
      for j in 0..<arr.length() {
        if j > 0 {
          let _ = sb.write("/")
        }
        let _ = sb.write(arr[j])
      }
      sb.to_string()
    }
  }
  normalize(gt1) == normalize(gt2)
}

///|
pub fn pysam3_vcf_concordance(
  records1 : Array[Pysam3VcfRecord],
  records2 : Array[Pysam3VcfRecord],
) -> (Int32, Int32, Int32) {
  // returns (concordant, discordant, missing)
  let mut concordant = 0
  let mut discordant = 0
  let mut missing = 0
  for i in 0..<records1.length() {
    let r1 = records1[i]
    let mut found = false
    for j in 0..<records2.length() {
      let r2 = records2[j]
      if r1.chrom == r2.chrom && r1.pos == r2.pos {
        found = true
        if r1.genotypes.length() > 0 && r2.genotypes.length() > 0 {
          if pysam3_genotype_concordance(r1.genotypes[0], r2.genotypes[0]) {
            concordant = concordant + 1
          } else {
            discordant = discordant + 1
          }
        } else {
          missing = missing + 1
        }
        break
      }
    }
    if not(found) {
      missing = missing + 1
    }
  }
  (concordant, discordant, missing)
}

// --- Heterozygosity ---

///|
pub fn pysam3_is_heterozygous(gt : String) -> Bool {
  let parts = gt.split("/").to_array()
  if parts.length() == 2 {
    parts[0].to_owned() != parts[1].to_owned()
  } else {
    let parts2 = gt.split("|").to_array()
    if parts2.length() == 2 {
      parts2[0].to_owned() != parts2[1].to_owned()
    } else {
      false
    }
  }
}

///|
pub fn pysam3_heterozygosity_rate(
  records : Array[Pysam3VcfRecord],
  sample_idx : Int32,
) -> Double {
  let mut het = 0
  let mut total = 0
  for i in 0..<records.length() {
    let rec = records[i]
    if sample_idx < rec.genotypes.length() {
      total = total + 1
      if pysam3_is_heterozygous(rec.genotypes[sample_idx]) {
        het = het + 1
      }
    }
  }
  if total == 0 {
    return 0.0
  }
  het.to_double() / total.to_double()
}

// --- BCF binary support (header marker) ---

///|
pub fn pysam3_bcf_magic() -> Array[Int32] {
  // BCF magic bytes: BCF\2\2
  [66, 67, 70, 2, 2]
}

///|
pub fn pysam3_bcf_encode_record(rec : Pysam3VcfRecord) -> Array[Int32] {
  // Simplified BCF record encoding (text representation as byte array)
  let text = pysam3_vcf_record_to_text(rec)
  let bytes : Array[Int32] = []
  for i in 0..<text.length() {
    bytes.push(text[i].to_int())
  }
  bytes
}

///|
pub fn pysam3_bcf_decode_record(bytes : Array[Int32]) -> String {
  let sb = StringBuilder()
  for i in 0..<bytes.length() {
    let _ = sb.write_char(Int32::unsafe_to_char(bytes[i]))
  }
  sb.to_string()
}

///|
pub fn pysam3_vcf_header_dummy() -> Pysam3VcfHeader {
  {
    version: "VCFv4.2",
    contigs: [],
    samples: [],
    info_fields: [],
    format_fields: [],
  }
}

///|
pub fn pysam3_vcf_record_dummy() -> Pysam3VcfRecord {
  {
    chrom: "",
    pos: 0,
    id_: ".",
    ref_: "A",
    alt: "T",
    qual: 0.0,
    filter_: "PASS",
    info: ".",
    format_: "",
    genotypes: [],
  }
}
