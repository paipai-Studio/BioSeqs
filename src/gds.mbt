///|
/// SeqArray - GDS-style hierarchical container for sequence variant data.
/// Inspired by the Bioconductor SeqArray package (Zheng et al.). Stores
/// sample/variant/genotype data in a tree of typed nodes and provides
/// subset filters over samples and variants, mirroring seqSetFilter and
/// seqGetData from the R package.

///|
/// Value stored in a GDS node.
pub(all) enum GdsValue {
  GdsString(String)
  GdsInt(Int)
  GdsDouble(Double)
  GdsStringArray(Array[String])
  GdsIntArray(Array[Int])
  GdsDoubleArray(Array[Double])
  GdsIntMatrix(Array[Array[Int]])
  GdsInt3D(Array[Array[Array[Int]]])
}

///|
/// A node in the GDS tree. Folder nodes have children; leaf nodes carry
/// a value.

///|
pub struct GdsNode {
  name : String
  description : String
  mut value : GdsValue?
  children : Array[GdsNode]
}

///|
/// Create a leaf node holding a value.
pub fn GdsNode::leaf(
  name : String,
  description : String,
  value : GdsValue,
) -> GdsNode {
  GdsNode::{ name, description, value: Some(value), children: [], }
}

///|
/// Create a folder node.
pub fn GdsNode::folder(name : String, description : String) -> GdsNode {
  GdsNode::{ name, description, value: None, children: [], }
}

///|
/// A GDS file: rooted tree of nodes plus simple path addressing.
pub struct GdsFile {
  root : GdsNode
}

///|
/// Create an empty GDS file.
pub fn gds_file_new() -> GdsFile {
  GdsFile::{ root: GdsNode::folder("", "root"), }
}

///|
/// Add a child node to a folder addressed by path (e.g. "annotation/info").
/// Intermediate folders are created as needed. If a path component
/// currently holds a leaf value it is converted into a folder. Returns
/// the new node.
pub fn GdsFile::add_node(
  self : GdsFile,
  path : String,
  node : GdsNode,
) -> GdsNode {
  let mut current = self.root
  if path.length() > 0 {
    let parts = gds_split_path(path)
    for part in parts {
      let mut found : GdsNode? = None
      for child in current.children {
        if child.name == part {
          found = Some(child)
          break
        }
      }
      match found {
        Some(existing) => {
          if existing.value is Some(_) {
            existing.value = None
          }
          current = existing
        }
        None => {
          let folder = GdsNode::folder(part, "")
          current.children.push(folder)
          current = folder
        }
      }
    }
  }
  current.children.push(node)
  node
}

///|
/// Fetch the node at the given path; raises when absent.
pub fn GdsFile::get_node(
  self : GdsFile,
  path : String,
) -> GdsNode raise SeqIOError {
  let mut current = self.root
  for part in gds_split_path(path) {
    let mut found : GdsNode? = None
    for child in current.children {
      if child.name == part {
        found = Some(child)
        break
      }
    }
    match found {
      Some(next) => current = next
      None => raise SeqIOError("GDS node not found: " + path)
    }
  }
  current
}

///|
/// Fetch the node at the given path; returns None when absent.
pub fn GdsFile::try_get_node(self : GdsFile, path : String) -> GdsNode? {
  let mut current = self.root
  for part in gds_split_path(path) {
    let mut found : GdsNode? = None
    for child in current.children {
      if child.name == part {
        found = Some(child)
        break
      }
    }
    match found {
      Some(next) => current = next
      None => return None
    }
  }
  Some(current)
}

///|
/// Check whether a node exists.
pub fn GdsFile::has_node(self : GdsFile, path : String) -> Bool {
  self.try_get_node(path) is Some(_)
}

///|
fn gds_split_path(path : String) -> Array[String] {
  let parts : Array[String] = []
  let mut start = 0
  let mut i = 0
  while i <= path.length() {
    if i == path.length() || path.unsafe_get(i) == '/' {
      if i > start {
        parts.push(path.substring(start~, end=i))
      }
      start = i + 1
    }
    i = i + 1
  }
  parts
}

///|
/// Parse one GT field token ("0/1", "1|0", "0", "./.") into two allele
/// indices; missing alleles become -1. A haploid token ("0") duplicates
/// the single allele into the first slot and leaves the second missing.
fn gds_parse_gt_token(token : String) -> (Int, Int) {
  let n = token.length()
  let alleles : Array[Int] = []
  let mut current = 0
  let mut has_digits = false
  let mut i = 0
  while i < n {
    let ch = token.unsafe_get(i)
    if ch == '/' || ch == '|' {
      alleles.push(if has_digits { current } else { -1 })
      current = 0
      has_digits = false
    } else if ch >= '0' && ch <= '9' {
      current = current * 10 + (ch.to_int() - '0'.to_int())
      has_digits = true
    }
    i = i + 1
  }
  if has_digits || n == 0 {
    alleles.push(if has_digits { current } else { -1 })
  } else if alleles.length() > 0 {
    alleles.push(-1)
  }
  let a1 = if alleles.length() > 0 { alleles[0] } else { -1 }
  let a2 = if alleles.length() > 1 { alleles[1] } else { -1 }
  (a1, a2)
}

///|
/// Import a parsed VCF into a SeqArray-style GDS container, mirroring
/// seqVCF2GDS. Stores sample.id, variant.id, chromosome, position,
/// allele (REF,ALT,...) and genotypes (2 x sample x variant dosage-free
/// allele indices, -1 for missing).
pub fn seq_vcf2gds(vcf : VcfFile) -> SeqVarData {
  let file = gds_file_new()
  let sample_names = vcf.get_sample_names()
  let n_samples = sample_names.length()
  let n_variants = vcf.len()
  let _ = n_variants

  let _ = file.add_node(
    "",
    GdsNode::leaf(
      "sample.id",
      "sample identifiers",
      GdsValue::GdsStringArray(sample_names),
    ),
  )

  let variant_ids : Array[String] = Array::new(capacity=vcf.len())
  let chroms : Array[String] = Array::new(capacity=vcf.len())
  let positions : Array[Int] = Array::new(capacity=vcf.len())
  let alleles : Array[Array[String]] = Array::new(capacity=vcf.len())
  // genotypes[variant][sample] = [a1, a2]
  let genotypes : Array[Array[Array[Int]]] = Array::new(capacity=vcf.len())

  let mut v = 0
  while v < vcf.len() {
    let rec = vcf.records[v]
    let vid = if rec.id.length() > 0 && rec.id[0] != "." {
      rec.id[0]
    } else {
      rec.chrom + ":" + rec.pos.to_string() + ":" + rec.ref_seq
    }
    variant_ids.push(vid)
    chroms.push(rec.chrom)
    positions.push(rec.pos)
    let alts : Array[String] = Array::new(capacity=rec.alt.length() + 1)
    alts.push(rec.ref_seq)
    for a in rec.alt {
      alts.push(a)
    }
    alleles.push(alts)

    let per_sample : Array[Array[Int]] = Array::new(capacity=n_samples)
    for s in 0..<n_samples {
      let name = sample_names[s]
      let mut a1 = -1
      let mut a2 = -1
      match rec.samples.get(name) {
        Some(fields) =>
          if fields.length() > 0 {
            let (f1, f2) = gds_parse_gt_token(fields[0])
            a1 = f1
            a2 = f2
          }
        None => ()
      }
      per_sample.push([a1, a2])
    }
    genotypes.push(per_sample)
    v = v + 1
  }

  let _ = file.add_node(
    "",
    GdsNode::leaf(
      "variant.id",
      "variant identifiers",
      GdsValue::GdsStringArray(variant_ids),
    ),
  )
  let _ = file.add_node(
    "",
    GdsNode::leaf("chromosome", "chromosome", GdsValue::GdsStringArray(chroms)),
  )
  let _ = file.add_node(
    "",
    GdsNode::leaf("position", "position", GdsValue::GdsIntArray(positions)),
  )
  let _ = file.add_node(
    "",
    GdsNode::leaf(
      "allele",
      "alleles",
      GdsValue::GdsStringArray(gds_flatten_alleles(alleles)),
    ),
  )
  let _ = file.add_node(
    "",
    GdsNode::leaf(
      "genotype",
      "genotype data (2 x sample x variant)",
      GdsValue::GdsInt3D(genotypes),
    ),
  )
  SeqVarData::new(file)
}

///|
fn gds_flatten_alleles(alleles : Array[Array[String]]) -> Array[String] {
  let out : Array[String] = []
  for entry in alleles {
    let mut joined = ""
    for idx in 0..<entry.length() {
      if idx > 0 {
        joined = joined + ","
      }
      joined = joined + entry[idx]
    }
    out.push(joined)
  }
  out
}

///|
/// SeqArray data container with sample/variant filters, mirroring
/// SeqVarGDSClass plus seqSetFilter.
pub struct SeqVarData {
  file : GdsFile
  mut sample_filter : Array[Bool]
  mut variant_filter : Array[Bool]
}

///|
/// Build a SeqVarData wrapper around a GDS file with all-pass filters.
pub fn SeqVarData::new(file : GdsFile) -> SeqVarData {
  let n_samples = seq_gds_count(file, "sample.id")
  let n_variants = seq_gds_count(file, "variant.id")
  SeqVarData::{
    file,
    sample_filter: Array::make(n_samples, true),
    variant_filter: Array::make(n_variants, true),
  }
}

///|
fn seq_gds_count(file : GdsFile, node : String) -> Int {
  match file.try_get_node(node) {
    Some(node_ref) =>
      match node_ref.value {
        Some(GdsValue::GdsStringArray(items)) => items.length()
        Some(GdsValue::GdsIntArray(items)) => items.length()
        _ => 0
      }
    None => 0
  }
}

///|
/// Number of samples passing the current filter.
pub fn SeqVarData::n_samples(self : SeqVarData) -> Int {
  let mut n = 0
  for keep in self.sample_filter {
    if keep {
      n = n + 1
    }
  }
  n
}

///|
/// Total number of variants (unfiltered).
pub fn SeqVarData::n_variants_total(self : SeqVarData) -> Int {
  seq_gds_count(self.file, "variant.id")
}

///|
/// Number of variants passing the current filter.
pub fn SeqVarData::n_variants(self : SeqVarData) -> Int {
  let mut n = 0
  for keep in self.variant_filter {
    if keep {
      n = n + 1
    }
  }
  n
}

///|
/// Filtered sample ids in file order.
pub fn SeqVarData::sample_ids(self : SeqVarData) -> Array[String] {
  seq_gds_strings_filtered(self.file, "sample.id", self.sample_filter)
}

///|
/// Filtered variant ids in file order.
pub fn SeqVarData::variant_ids(self : SeqVarData) -> Array[String] {
  seq_gds_strings_filtered(self.file, "variant.id", self.variant_filter)
}

///|
/// Filtered chromosome names.
pub fn SeqVarData::chromosomes(self : SeqVarData) -> Array[String] {
  seq_gds_strings_filtered(self.file, "chromosome", self.variant_filter)
}

///|
/// Filtered positions.
pub fn SeqVarData::positions(self : SeqVarData) -> Array[Int] {
  let out : Array[Int] = []
  match self.file.try_get_node("position") {
    Some(node_ref) =>
      match node_ref.value {
        Some(GdsValue::GdsIntArray(items)) => {
          let mut i = 0
          while i < items.length() {
            if i < self.variant_filter.length() && self.variant_filter[i] {
              out.push(items[i])
            }
            i = i + 1
          }
        }
        _ => ()
      }
    None => ()
  }
  out
}

///|
/// Filtered allele strings ("REF,ALT[,...]").
pub fn SeqVarData::alleles(self : SeqVarData) -> Array[String] {
  seq_gds_strings_filtered(self.file, "allele", self.variant_filter)
}

///|
fn seq_gds_strings_filtered(
  file : GdsFile,
  node : String,
  filter : Array[Bool],
) -> Array[String] {
  let out : Array[String] = []
  match file.try_get_node(node) {
    Some(node_ref) =>
      match node_ref.value {
        Some(GdsValue::GdsStringArray(items)) => {
          let mut i = 0
          while i < items.length() {
            if i < filter.length() && filter[i] {
              out.push(items[i])
            }
            i = i + 1
          }
        }
        _ => ()
      }
    None => ()
  }
  out
}

///|
/// Reset both filters to include everything, mirroring seqResetFilter.
pub fn SeqVarData::reset_filter(self : SeqVarData) -> Unit {
  let mut i = 0
  while i < self.sample_filter.length() {
    self.sample_filter[i] = true
    i = i + 1
  }
  i = 0
  while i < self.variant_filter.length() {
    self.variant_filter[i] = true
    i = i + 1
  }
}

///|
/// Keep only the given variant indices, mirroring seqSetFilter(object,
/// variant.id=...).
pub fn SeqVarData::set_variant_filter(
  self : SeqVarData,
  indices : Array[Int],
) -> Unit raise SeqIOError {
  let total = self.n_variants_total()
  let keep = Array::make(total, false)
  for idx in indices {
    if idx < 0 || idx >= total {
      raise SeqIOError("variant index out of range: " + idx.to_string())
    }
    keep[idx] = true
  }
  self.variant_filter = keep
}

///|
/// Keep only variants on the given chromosomes.
pub fn SeqVarData::set_variant_chrom_filter(
  self : SeqVarData,
  chroms : Array[String],
) -> Unit {
  let total = self.n_variants_total()
  let keep = Array::make(total, false)
  match self.file.try_get_node("chromosome") {
    Some(node_ref) =>
      match node_ref.value {
        Some(GdsValue::GdsStringArray(items)) => {
          let mut i = 0
          while i < items.length() && i < total {
            for want in chroms {
              if items[i] == want {
                keep[i] = true
                break
              }
            }
            i = i + 1
          }
        }
        _ => ()
      }
    None => ()
  }
  self.variant_filter = keep
}

///|
/// Keep only variants inside the half-open position range [start, end).
pub fn SeqVarData::set_variant_range_filter(
  self : SeqVarData,
  start : Int,
  end : Int,
) -> Unit {
  let total = self.n_variants_total()
  let keep = Array::make(total, false)
  match self.file.try_get_node("position") {
    Some(node_ref) =>
      match node_ref.value {
        Some(GdsValue::GdsIntArray(items)) => {
          let mut i = 0
          while i < items.length() && i < total {
            if items[i] >= start && items[i] < end {
              keep[i] = true
            }
            i = i + 1
          }
        }
        _ => ()
      }
    None => ()
  }
  self.variant_filter = keep
}

///|
/// Keep only the given sample indices, mirroring seqSetFilter(object,
/// sample.id=...).
pub fn SeqVarData::set_sample_filter(
  self : SeqVarData,
  indices : Array[Int],
) -> Unit raise SeqIOError {
  let total = self.sample_filter.length()
  let keep = Array::make(total, false)
  for idx in indices {
    if idx < 0 || idx >= total {
      raise SeqIOError("sample index out of range: " + idx.to_string())
    }
    keep[idx] = true
  }
  self.sample_filter = keep
}

///|
/// Filtered genotype dosages: rows are variants, columns are filtered
/// samples; value is the number of non-reference alleles (0/1/2) and -1
/// marks missing calls. Haploid calls ("1") count the single allele.
pub fn SeqVarData::genotype_dosages(self : SeqVarData) -> Array[Array[Int]] {
  let dosages : Array[Array[Int]] = []
  match self.file.try_get_node("genotype") {
    Some(node_ref) =>
      match node_ref.value {
        Some(GdsValue::GdsInt3D(data)) => {
          let mut v = 0
          while v < data.length() {
            if v < self.variant_filter.length() && self.variant_filter[v] {
              let row : Array[Int] = []
              let mut s = 0
              while s < data[v].length() {
                if s < self.sample_filter.length() && self.sample_filter[s] {
                  let a1 = data[v][s][0]
                  let a2 = data[v][s][1]
                  let dosage = if a1 < 0 {
                    -1
                  } else if a2 < 0 {
                    // haploid call: count the single allele
                    if a1 > 0 {
                      1
                    } else {
                      0
                    }
                  } else {
                    (if a1 > 0 { 1 } else { 0 }) + (if a2 > 0 { 1 } else { 0 })
                  }
                  row.push(dosage)
                }
                s = s + 1
              }
              dosages.push(row)
            }
            v = v + 1
          }
        }
        _ => ()
      }
    None => ()
  }
  dosages
}

///|
/// Per-variant allele frequency of the alternate alleles over filtered
/// samples (missing calls excluded). Multi-allelic sites count every
/// non-reference allele; haploid calls contribute one allele.
pub fn SeqVarData::allele_frequencies(self : SeqVarData) -> Array[Double] {
  let out : Array[Double] = []
  match self.file.try_get_node("genotype") {
    Some(node_ref) =>
      match node_ref.value {
        Some(GdsValue::GdsInt3D(data)) => {
          let mut v = 0
          while v < data.length() {
            if v < self.variant_filter.length() && self.variant_filter[v] {
              let mut called = 0
              let mut alt = 0
              let mut s = 0
              while s < data[v].length() {
                if s < self.sample_filter.length() && self.sample_filter[s] {
                  for a in data[v][s] {
                    if a >= 0 {
                      called = called + 1
                      if a > 0 {
                        alt = alt + 1
                      }
                    }
                  }
                }
                s = s + 1
              }
              out.push(
                if called > 0 {
                  alt.to_double() / called.to_double()
                } else {
                  0.0
                },
              )
            }
            v = v + 1
          }
        }
        _ => ()
      }
    None => ()
  }
  out
}

///|
/// Per-variant missing call rate over filtered samples.
pub fn SeqVarData::missing_rates(self : SeqVarData) -> Array[Double] {
  let out : Array[Double] = []
  for row in self.genotype_dosages() {
    let mut missing = 0
    for dosage in row {
      if dosage < 0 {
        missing = missing + 1
      }
    }
    out.push(
      if row.length() > 0 {
        missing.to_double() / row.length().to_double()
      } else {
        0.0
      },
    )
  }
  out
}

///|
/// Summary of the filtered dataset, mirroring seqSummary.
pub struct SeqVarSummary {
  n_samples : Int
  n_variants : Int
  n_chromosomes : Int
  mean_allele_frequency : Double
  mean_missing_rate : Double
}

///|
/// Compute a summary over the current filter.
pub fn SeqVarData::summary(self : SeqVarData) -> SeqVarSummary {
  let freqs = self.allele_frequencies()
  let missing = self.missing_rates()
  let mut freq_sum = 0.0
  for f in freqs {
    freq_sum = freq_sum + f
  }
  let mut missing_sum = 0.0
  for m in missing {
    missing_sum = missing_sum + m
  }
  let chroms = self.chromosomes()
  let seen : Array[String] = []
  for c in chroms {
    let mut found = false
    for s in seen {
      if s == c {
        found = true
        break
      }
    }
    if !found {
      seen.push(c)
    }
  }
  SeqVarSummary::{
    n_samples: self.n_samples(),
    n_variants: self.n_variants(),
    n_chromosomes: seen.length(),
    mean_allele_frequency: if freqs.length() > 0 {
      freq_sum / freqs.length().to_double()
    } else {
      0.0
    },
    mean_missing_rate: if missing.length() > 0 {
      missing_sum / missing.length().to_double()
    } else {
      0.0
    },
  }
}

///|
/// Apply a function to each filtered variant and collect the results.
pub fn SeqVarData::apply_per_variant(
  self : SeqVarData,
  f : (Array[Int], Int) -> Double,
) -> Array[Double] {
  let dosages = self.genotype_dosages()
  let out : Array[Double] = []
  let mut v = 0
  while v < dosages.length() {
    out.push(f(dosages[v], v))
    v = v + 1
  }
  out
}
