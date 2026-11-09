// Pysam third-round: Variant set operations (items 29-31)
// VCF merge, set operations, variant density

///|
pub fn pysam3_vcf_sort_records(
  records : Array[Pysam3VcfRecord],
) -> Array[Pysam3VcfRecord] {
  // Sort by (chrom, pos) — simple chromosome string sort + position
  let arr = []
  for i in 0..<records.length() {
    arr.push(records[i])
  }
  arr.sort_by(fn(a : Pysam3VcfRecord, b : Pysam3VcfRecord) -> Int32 {
    if a.chrom < b.chrom {
      -1
    } else if a.chrom > b.chrom {
      1
    } else if a.pos < b.pos {
      -1
    } else if a.pos > b.pos {
      1
    } else {
      0
    }
  })
  arr
}

///|
pub fn pysam3_vcf_merge(
  files : Array[Array[Pysam3VcfRecord]],
) -> Array[Pysam3VcfRecord] {
  // Merge multiple VCF record arrays, sort, and deduplicate
  let all : Array[Pysam3VcfRecord] = []
  for i in 0..<files.length() {
    for j in 0..<files[i].length() {
      all.push(files[i][j])
    }
  }
  let sorted = pysam3_vcf_sort_records(all)
  // Deduplicate
  let result : Array[Pysam3VcfRecord] = []
  for i in 0..<sorted.length() {
    if i > 0 {
      let prev = sorted[i - 1]
      let cur = sorted[i]
      if prev.chrom == cur.chrom &&
        prev.pos == cur.pos &&
        prev.ref_ == cur.ref_ &&
        prev.alt == cur.alt {
        continue
      }
    }
    result.push(sorted[i])
  }
  result
}

// --- Set operations ---

///|
fn pysam3_vcf_record_key(rec : Pysam3VcfRecord) -> String {
  rec.chrom + ":" + rec.pos.to_string() + ":" + rec.ref_ + ":" + rec.alt
}

///|
pub fn pysam3_vcf_intersect(
  records1 : Array[Pysam3VcfRecord],
  records2 : Array[Pysam3VcfRecord],
) -> Array[Pysam3VcfRecord] {
  // Variants present in both files
  let result : Array[Pysam3VcfRecord] = []
  for i in 0..<records1.length() {
    let key1 = pysam3_vcf_record_key(records1[i])
    let mut found = false
    for j in 0..<records2.length() {
      let key2 = pysam3_vcf_record_key(records2[j])
      if key1 == key2 {
        found = true
        break
      }
    }
    if found {
      result.push(records1[i])
    }
  }
  result
}

///|
pub fn pysam3_vcf_union(
  records1 : Array[Pysam3VcfRecord],
  records2 : Array[Pysam3VcfRecord],
) -> Array[Pysam3VcfRecord] {
  // All unique variants from both files
  let all : Array[Pysam3VcfRecord] = []
  for i in 0..<records1.length() {
    all.push(records1[i])
  }
  for i in 0..<records2.length() {
    all.push(records2[i])
  }
  let sorted = pysam3_vcf_sort_records(all)
  let result : Array[Pysam3VcfRecord] = []
  for i in 0..<sorted.length() {
    if i > 0 {
      let prev = sorted[i - 1]
      let cur = sorted[i]
      if pysam3_vcf_record_key(prev) == pysam3_vcf_record_key(cur) {
        continue
      }
    }
    result.push(sorted[i])
  }
  result
}

///|
pub fn pysam3_vcf_difference(
  records1 : Array[Pysam3VcfRecord],
  records2 : Array[Pysam3VcfRecord],
) -> Array[Pysam3VcfRecord] {
  // Variants in records1 but NOT in records2
  let result : Array[Pysam3VcfRecord] = []
  for i in 0..<records1.length() {
    let key1 = pysam3_vcf_record_key(records1[i])
    let mut found = false
    for j in 0..<records2.length() {
      let key2 = pysam3_vcf_record_key(records2[j])
      if key1 == key2 {
        found = true
        break
      }
    }
    if not(found) {
      result.push(records1[i])
    }
  }
  result
}

///|
pub fn pysam3_vcf_symmetric_diff(
  records1 : Array[Pysam3VcfRecord],
  records2 : Array[Pysam3VcfRecord],
) -> Array[Pysam3VcfRecord] {
  // Variants in either but not both
  let diff1 = pysam3_vcf_difference(records1, records2)
  let diff2 = pysam3_vcf_difference(records2, records1)
  pysam3_vcf_union(diff1, diff2)
}

// --- Variant density ---

///|
pub struct Pysam3VariantDensity {
  chrom : String
  start : Int32
  end_ : Int32
  count : Int32
  density : Double // variants per bp
}

///|
pub fn pysam3_variant_density(
  records : Array[Pysam3VcfRecord],
  chrom : String,
  window_size : Int32,
) -> Array[Pysam3VariantDensity] {
  // Calculate variant density in sliding windows
  let chrom_records : Array[Pysam3VcfRecord] = []
  for i in 0..<records.length() {
    if records[i].chrom == chrom {
      chrom_records.push(records[i])
    }
  }
  if chrom_records.length() == 0 {
    return []
  }
  let sorted = pysam3_vcf_sort_records(chrom_records)
  let min_pos = sorted[0].pos
  let max_pos = sorted[sorted.length() - 1].pos
  let result : Array[Pysam3VariantDensity] = []
  let mut win_start = min_pos
  while win_start <= max_pos {
    let win_end = win_start + window_size
    let mut count = 0
    for i in 0..<sorted.length() {
      if sorted[i].pos >= win_start && sorted[i].pos < win_end {
        count = count + 1
      }
    }
    result.push({
      chrom,
      start: win_start,
      end_: win_end,
      count,
      density: count.to_double() / window_size.to_double(),
    })
    win_start = win_end
  }
  result
}

///|
pub fn pysam3_variant_density_summary(
  density : Array[Pysam3VariantDensity],
) -> (Double, Double, Double) {
  // (mean_density, max_density, total_variants)
  if density.length() == 0 {
    return (0.0, 0.0, 0.0)
  }
  let mut sum = 0.0
  let mut max_d = 0.0
  let mut total = 0.0
  for i in 0..<density.length() {
    let d = density[i]
    sum = sum + d.density
    if d.density > max_d {
      max_d = d.density
    }
    total = total + d.count.to_double()
  }
  (sum / density.length().to_double(), max_d, total)
}

///|
pub fn pysam3_snp_density_by_chrom(
  records : Array[Pysam3VcfRecord],
) -> Array[(String, Int32, Double)] {
  // (chrom, snp_count, snps_per_kb)
  let chrom_counts : Array[(String, Int32)] = []
  for i in 0..<records.length() {
    let rec = records[i]
    if rec.ref_.length() == 1 && rec.alt.length() == 1 {
      let mut found = false
      for j in 0..<chrom_counts.length() {
        let (c, cnt) = chrom_counts[j]
        if c == rec.chrom {
          chrom_counts[j] = (c, cnt + 1)
          found = true
          break
        }
      }
      if not(found) {
        chrom_counts.push((rec.chrom, 1))
      }
    }
  }
  // Calculate per-kb density (assuming we don't have chrom lengths, use position range)
  let result : Array[(String, Int32, Double)] = []
  for i in 0..<chrom_counts.length() {
    let (chrom, count) = chrom_counts[i]
    // Find max position for this chrom
    let mut max_pos = 0
    for j in 0..<records.length() {
      if records[j].chrom == chrom && records[j].pos > max_pos {
        max_pos = records[j].pos
      }
    }
    let per_kb = if max_pos > 0 {
      count.to_double() / (max_pos / 1000 + 1).to_double()
    } else {
      0.0
    }
    result.push((chrom, count, per_kb))
  }
  result
}
