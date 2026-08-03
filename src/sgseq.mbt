// SGSeq - Splice variant detection from RNA-seq data.
// Inspired by Bioconductor SGSeq package.
//
// Detects splice variants from RNA-seq junction reads by building
// splice graphs and identifying alternative splicing events:
//   - Skipped Exon (SE)
//   - Alternative 5' Splice Site (A5SS)
//   - Alternative 3' Splice Site (A3SS)
//   - Mutually Exclusive Exons (MXE)
//   - Retained Intron (RI)
//
// References:
//   - Goldstein LD et al. (2016) Prediction and quantification of splice
//     isoforms from RNA-seq data with SGSeq. Bioinformatics 32:1188-1195.

// ===========================================================================
// Data structures
// ===========================================================================

///|
/// Type of splice variant.
pub enum SGSeqVariantType {
  VariantSE // Skipped Exon
  VariantA5SS // Alternative 5' Splice Site
  VariantA3SS // Alternative 3' Splice Site
  VariantMXE // Mutually Exclusive Exons
  VariantRI // Retained Intron
}

///|
/// Convert variant type to string.
pub fn SGSeqVariantType::to_string(self : SGSeqVariantType) -> String {
  match self {
    VariantSE => "SE"
    VariantA5SS => "A5SS"
    VariantA3SS => "A3SS"
    VariantMXE => "MXE"
    VariantRI => "RI"
  }
}

///|
/// Factory function for cross-package SGSeqVariantType construction.
pub fn sgseq_make_variant_type(name : String) -> SGSeqVariantType {
  match name {
    "SE" => VariantSE
    "A5SS" => VariantA5SS
    "A3SS" => VariantA3SS
    "MXE" => VariantMXE
    "RI" => VariantRI
    _ => VariantSE
  }
}

///|
/// A splice junction from split-read alignments.
pub struct SGSeqJunction {
  pub chr : String
  pub start : Int // donor (5') site (1-based, inclusive)
  pub end : Int // acceptor (3') site (1-based, inclusive)
  pub strand : String // "+", "-", or "*"
  pub count : Int // number of supporting reads
}

///|
/// Construct a SGSeqJunction.
pub fn SGSeqJunction::new(
  chr~ : String,
  start~ : Int,
  end~ : Int,
  strand~ : String = "+",
  count~ : Int = 1,
) -> SGSeqJunction {
  SGSeqJunction::{ chr, start, end, strand, count }
}

///|
/// Get junction span (length of intron).
pub fn SGSeqJunction::intron_length(self : SGSeqJunction) -> Int {
  self.end - self.start - 1
}

///|
/// Check if two junctions share the same donor (5') site.
pub fn SGSeqJunction::same_donor(
  self : SGSeqJunction,
  other : SGSeqJunction,
) -> Bool {
  self.chr == other.chr && self.start == other.start
}

///|
/// Check if two junctions share the same acceptor (3') site.
pub fn SGSeqJunction::same_acceptor(
  self : SGSeqJunction,
  other : SGSeqJunction,
) -> Bool {
  self.chr == other.chr && self.end == other.end
}

///|
/// An exon inferred from junction coordinates.
pub struct SGSeqExon {
  pub chr : String
  pub start : Int
  pub end : Int
  pub strand : String
  pub count : Int
  pub mut is_cassette : Bool // potentially skipped exon
}

///|
/// Construct a SGSeqExon.
pub fn SGSeqExon::new(
  chr~ : String,
  start~ : Int,
  end~ : Int,
  strand~ : String = "+",
) -> SGSeqExon {
  SGSeqExon::{
    chr,
    start,
    end,
    strand,
    count: 0,
    is_cassette: false,
  }
}

///|
/// Get exon length.
pub fn SGSeqExon::length(self : SGSeqExon) -> Int {
  self.end - self.start + 1
}

///|
/// Check if this exon contains a position.
pub fn SGSeqExon::contains(self : SGSeqExon, pos : Int) -> Bool {
  pos >= self.start && pos <= self.end
}

///|
/// A splice graph feature (exon or junction).
pub struct SGSeqFeature {
  pub type_ : String // "E" for exon, "J" for junction
  pub chr : String
  pub start : Int
  pub end : Int
  pub strand : String
  pub count : Int
}

///|
/// Construct a SGSeqFeature.
pub fn SGSeqFeature::new(
  type_~ : String,
  chr~ : String,
  start~ : Int,
  end~ : Int,
  strand~ : String = "+",
  count~ : Int = 0,
) -> SGSeqFeature {
  SGSeqFeature::{ type_, chr, start, end, strand, count }
}

///|
/// A detected splice variant event.
pub struct SGSeqVariant {
  pub chr : String
  pub start : Int
  pub end : Int
  pub strand : String
  pub variant_type : SGSeqVariantType
  pub mut inclusion_count : Int // reads supporting inclusion
  pub mut exclusion_count : Int // reads supporting exclusion
  pub mut psi : Double // Percent Spliced In (0..1)
  pub mut details : String // event-specific coordinates
}

///|
/// Construct a SGSeqVariant.
pub fn SGSeqVariant::new(
  chr~ : String,
  start~ : Int,
  end~ : Int,
  strand~ : String = "+",
  variant_type~ : SGSeqVariantType,
) -> SGSeqVariant {
  SGSeqVariant::{
    chr,
    start,
    end,
    strand,
    variant_type,
    inclusion_count: 0,
    exclusion_count: 0,
    psi: 0.0,
    details: "",
  }
}

///|
/// Get the variant type as string.
pub fn SGSeqVariant::type_name(self : SGSeqVariant) -> String {
  self.variant_type.to_string()
}

///|
/// Set inclusion count.
pub fn SGSeqVariant::set_inclusion_count(self : SGSeqVariant, val : Int) -> Unit {
  self.inclusion_count = val
}

///|
/// Set exclusion count.
pub fn SGSeqVariant::set_exclusion_count(self : SGSeqVariant, val : Int) -> Unit {
  self.exclusion_count = val
}

///|
/// Set PSI value.
pub fn SGSeqVariant::set_psi(self : SGSeqVariant, val : Double) -> Unit {
  self.psi = val
}

///|
/// Set details string.
pub fn SGSeqVariant::set_details(self : SGSeqVariant, val : String) -> Unit {
  self.details = val
}

///|
/// Compute PSI from inclusion/exclusion counts.
pub fn SGSeqVariant::compute_psi(self : SGSeqVariant) -> Double {
  let total = self.inclusion_count + self.exclusion_count
  if total == 0 {
    0.0
  } else {
    self.inclusion_count.to_double() / total.to_double()
  }
}

///|
/// Check if variant is significant (PSI between 0.1 and 0.9).
pub fn SGSeqVariant::is_significant(self : SGSeqVariant) -> Bool {
  let p = self.compute_psi()
  p > 0.1 && p < 0.9 && (self.inclusion_count + self.exclusion_count) >= 10
}

// ===========================================================================
// Junction parsing
// ===========================================================================

///|
/// Parse junctions from a tab-delimited text (one junction per line).
/// Format: chr<TAB>start<TAB>end<TAB>strand<TAB>count
/// Or STAR SJ.out.tab format: chr<TAB>start<TAB>end<TAB>strand(0/1/2)<TAB>intron_motif<TAB>annotated<TAB>unique<TAB>multi<TAB>max_overhang
pub fn sgseq_parse_junctions(text : String) -> Array[SGSeqJunction] {
  let junctions = Array::new()
  let lines = text.split("\n").to_array()
  for line in lines {
    let trimmed = line.to_string()
    if trimmed.length() == 0 {
      continue
    }
    // Skip comment lines
    if trimmed.has_prefix("#") {
      continue
    }
    let parts = trimmed.split("\t").to_array()
    if parts.length() < 3 {
      continue
    }
    let chr = parts[0].to_owned()
    let start = parse_int(parts[1].to_owned()) + 1 // convert 0-based to 1-based
    let end = parse_int(parts[2].to_owned())
    let strand = if parts.length() >= 4 {
      match parts[3].to_owned() {
        "0" | "+" => "+"
        "1" | "-" => "-"
        _ => "*"
      }
    } else {
      "*"
    }
    let count = if parts.length() >= 5 {
      parse_int(parts[4].to_owned())
    } else {
      1
    }
    if count > 0 && end > start {
      junctions.push(SGSeqJunction::new(
        chr=chr,
        start=start,
        end=end,
        strand=strand,
        count=count,
      ))
    }
  }
  junctions
}

// ===========================================================================
// Splice graph construction
// ===========================================================================

///|
/// Build splice graph exons from junctions.
/// Exons are segments between junction donor/acceptor sites.
pub fn sgseq_build_splice_graph(
  junctions : Array[SGSeqJunction],
) -> Array[SGSeqExon] {
  if junctions.length() == 0 {
    return []
  }
  // Collect all splice sites
  let sites = Array::new()
  for j in junctions {
    sites.push(j.start)
    sites.push(j.end)
  }
  // Sort sites (bubble sort)
  let mut i = 0
  while i < sites.length() {
    let mut j = i + 1
    while j < sites.length() {
      if sites[j] < sites[i] {
        let tmp = sites[i]
        sites[i] = sites[j]
        sites[j] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  // Remove duplicates
  let unique_sites = Array::new()
  for s in sites {
    if unique_sites.length() == 0 || unique_sites[unique_sites.length() - 1] != s {
      unique_sites.push(s)
    }
  }
  // Exons are between consecutive splice sites
  // Donor site = end of upstream exon (inclusive)
  // Acceptor site = start of downstream exon (inclusive)
  // So exon k: [unique_sites[2k-1]+1, unique_sites[2k]-1] approximately
  // Actually: junction (start, end) means:
  //   - upstream exon ends at 'start' (last base of exon)
  //   - downstream exon starts at 'end' (first base of exon)
  // So exon between two junctions: [j1.end, j2.start]
  let exons = Array::new()
  let chr = junctions[0].chr
  let strand = junctions[0].strand
  // Group junctions by chr
  let by_chr : Map[String, Array[SGSeqJunction]] = Map([], capacity=10)
  for j in junctions {
    let arr = by_chr.get(j.chr).or(Array::new())
    arr.push(j)
    by_chr[j.chr] = arr
  }
  for chr_name in by_chr.keys() {
    let chr_junctions = by_chr[chr_name]
    // Sort junctions by start
    let mut i = 0
    while i < chr_junctions.length() {
      let mut j = i + 1
      while j < chr_junctions.length() {
        if chr_junctions[j].start < chr_junctions[i].start {
          let tmp = chr_junctions[i]
          chr_junctions[i] = chr_junctions[j]
          chr_junctions[j] = tmp
        }
        j = j + 1
      }
      i = i + 1
    }
    // Exon 1: starts before first junction's donor, ends at first junction's start
    // Exon between j1.end and j2.start
    // Last exon: starts at last junction's end
    if chr_junctions.length() == 0 {
      continue
    }
    // First exon: ends at first junction's donor (start)
    exons.push(SGSeqExon::new(
      chr=chr_name,
      start=1,
      end=chr_junctions[0].start,
      strand=chr_junctions[0].strand,
    ))
    // Middle exons: between j[k].end and j[k+1].start
    let mut k = 0
    while k < chr_junctions.length() - 1 {
      exons.push(SGSeqExon::new(
        chr=chr_name,
        start=chr_junctions[k].end,
        end=chr_junctions[k + 1].start,
        strand=chr_junctions[k].strand,
      ))
      k = k + 1
    }
    // Last exon: starts at last junction's acceptor (end)
    let last = chr_junctions[chr_junctions.length() - 1]
    exons.push(SGSeqExon::new(
      chr=chr_name,
      start=last.end,
      end=last.end + 1000, // arbitrary end
      strand=last.strand,
    ))
  }
  exons
}

///|
/// Convert junctions and exons to features.
pub fn sgseq_to_features(
  junctions : Array[SGSeqJunction],
  exons : Array[SGSeqExon],
) -> Array[SGSeqFeature] {
  let features = Array::new()
  for e in exons {
    features.push(SGSeqFeature::new(
      type_="E",
      chr=e.chr,
      start=e.start,
      end=e.end,
      strand=e.strand,
      count=e.count,
    ))
  }
  for j in junctions {
    features.push(SGSeqFeature::new(
      type_="J",
      chr=j.chr,
      start=j.start,
      end=j.end,
      strand=j.strand,
      count=j.count,
    ))
  }
  features
}

// ===========================================================================
// Variant detection
// ===========================================================================

///|
/// Detect Skipped Exon (SE) events.
/// A cassette exon is skipped when a junction spans over it.
/// SE: exon E with junction (upstream_end, E_start) and (E_end, downstream_start)
///     and also a skipping junction (upstream_end, downstream_start).
pub fn sgseq_detect_se(
  junctions : Array[SGSeqJunction],
  exons : Array[SGSeqExon],
) -> Array[SGSeqVariant] {
  let variants = Array::new()
  for exon in exons {
    if exon.length() < 3 {
      continue
    }
    // Find junctions ending at exon.start (upstream junction into this exon)
    let upstream = Array::new()
    let downstream = Array::new()
    let skipping = Array::new()
    for j in junctions {
      if j.chr != exon.chr {
        continue
      }
      // Junction ending at this exon's start (acceptor)
      if j.end == exon.start {
        upstream.push(j)
      }
      // Junction starting at this exon's end (donor)
      if j.start == exon.end {
        downstream.push(j)
      }
      // Junction spanning over this exon (skipping)
      if j.start < exon.start && j.end > exon.end {
        skipping.push(j)
      }
    }
    // If we have both upstream and downstream and at least one skipping
    if upstream.length() > 0 && downstream.length() > 0 && skipping.length() > 0 {
      let mut inc_count = 0
      for j in upstream {
        inc_count = inc_count + j.count
      }
      for j in downstream {
        inc_count = inc_count + j.count
      }
      let mut exc_count = 0
      for j in skipping {
        exc_count = exc_count + j.count
      }
      let v = SGSeqVariant::new(
        chr=exon.chr,
        start=exon.start,
        end=exon.end,
        strand=exon.strand,
        variant_type=VariantSE,
      )
      v.inclusion_count = inc_count / 2 // each inclusion read uses 2 junctions
      v.exclusion_count = exc_count
      v.psi = v.compute_psi()
      v.details = "exon " + exon.start.to_string() + "-" + exon.end.to_string()
      variants.push(v)
      exon.is_cassette = true
    }
  }
  variants
}

///|
/// Detect Alternative 5' Splice Site (A5SS) events.
/// Same acceptor (3') site, different donor (5') sites.
pub fn sgseq_detect_a5ss(
  junctions : Array[SGSeqJunction],
) -> Array[SGSeqVariant] {
  let variants = Array::new()
  // Group by (chr, end) = same acceptor
  let groups : Map[String, Array[SGSeqJunction]] = Map([], capacity=20)
  for j in junctions {
    let key = j.chr + ":" + j.end.to_string()
    let arr = groups.get(key).or(Array::new())
    arr.push(j)
    groups[key] = arr
  }
  for key in groups.keys() {
    let group = groups[key]
    if group.length() < 2 {
      continue
    }
    // Different donor sites (start) with same acceptor (end)
    let starts = Array::new()
    for j in group {
      let mut found = false
      for s in starts {
        if s == j.start {
          found = true
        }
      }
      if !found {
        starts.push(j.start)
      }
    }
    if starts.length() >= 2 {
      // Sort starts
      let mut i = 0
      while i < starts.length() {
        let mut j = i + 1
        while j < starts.length() {
          if starts[j] < starts[i] {
            let tmp = starts[i]
            starts[i] = starts[j]
            starts[j] = tmp
          }
          j = j + 1
        }
        i = i + 1
      }
      let v = SGSeqVariant::new(
        chr=group[0].chr,
        start=starts[0],
        end=group[0].end,
        strand=group[0].strand,
        variant_type=VariantA5SS,
      )
      let mut total = 0
      for j in group {
        total = total + j.count
      }
      v.inclusion_count = total
      v.psi = 0.5 // simplified
      v.details = "donors " +
        starts[0].to_string() +
        "," +
        starts[starts.length() - 1].to_string() +
        " acceptor " +
        group[0].end.to_string()
      variants.push(v)
    }
  }
  variants
}

///|
/// Detect Alternative 3' Splice Site (A3SS) events.
/// Same donor (5') site, different acceptor (3') sites.
pub fn sgseq_detect_a3ss(
  junctions : Array[SGSeqJunction],
) -> Array[SGSeqVariant] {
  let variants = Array::new()
  // Group by (chr, start) = same donor
  let groups : Map[String, Array[SGSeqJunction]] = Map([], capacity=20)
  for j in junctions {
    let key = j.chr + ":" + j.start.to_string()
    let arr = groups.get(key).or(Array::new())
    arr.push(j)
    groups[key] = arr
  }
  for key in groups.keys() {
    let group = groups[key]
    if group.length() < 2 {
      continue
    }
    // Different acceptor sites (end) with same donor (start)
    let ends = Array::new()
    for j in group {
      let mut found = false
      for e in ends {
        if e == j.end {
          found = true
        }
      }
      if !found {
        ends.push(j.end)
      }
    }
    if ends.length() >= 2 {
      // Sort ends
      let mut i = 0
      while i < ends.length() {
        let mut j = i + 1
        while j < ends.length() {
          if ends[j] < ends[i] {
            let tmp = ends[i]
            ends[i] = ends[j]
            ends[j] = tmp
          }
          j = j + 1
        }
        i = i + 1
      }
      let v = SGSeqVariant::new(
        chr=group[0].chr,
        start=group[0].start,
        end=ends[ends.length() - 1],
        strand=group[0].strand,
        variant_type=VariantA3SS,
      )
      let mut total = 0
      for j in group {
        total = total + j.count
      }
      v.inclusion_count = total
      v.psi = 0.5
      v.details = "donor " +
        group[0].start.to_string() +
        " acceptors " +
        ends[0].to_string() +
        "," +
        ends[ends.length() - 1].to_string()
      variants.push(v)
    }
  }
  variants
}

///|
/// Detect Mutually Exclusive Exons (MXE) events.
/// Two exons where junctions skip one or the other but never both.
pub fn sgseq_detect_mxe(
  junctions : Array[SGSeqJunction],
  exons : Array[SGSeqExon],
) -> Array[SGSeqVariant] {
  let variants = Array::new()
  // Look for pairs of exons that are mutually exclusive
  let mut i = 0
  while i < exons.length() {
    let mut j = i + 1
    while j < exons.length() {
      let e1 = exons[i]
      let e2 = exons[j]
      if e1.chr != e2.chr {
        j = j + 1
        continue
      }
      // Check if there's a junction from before e1 to after e2 (skipping both)
      // and junctions including e1 but not e2, and including e2 but not e1
      let mut skip_both = 0
      let mut include_e1 = 0
      let mut include_e2 = 0
      for jc in junctions {
        if jc.chr != e1.chr {
          continue
        }
        // Junction spanning both exons
        if jc.start < e1.start && jc.end > e2.end {
          skip_both = skip_both + jc.count
        }
        // Junction ending at e1.start (including e1, skipping e2)
        if jc.end == e1.start {
          include_e1 = include_e1 + jc.count
        }
        // Junction ending at e2.start (including e2, skipping e1)
        if jc.end == e2.start {
          include_e2 = include_e2 + jc.count
        }
      }
      // MXE: both exons have inclusion but neither is always included
      if include_e1 > 0 && include_e2 > 0 && skip_both == 0 {
        let v = SGSeqVariant::new(
          chr=e1.chr,
          start=e1.start,
          end=e2.end,
          strand=e1.strand,
          variant_type=VariantMXE,
        )
        v.inclusion_count = include_e1
        v.exclusion_count = include_e2
        v.psi = v.compute_psi()
        v.details = "exon1 " +
          e1.start.to_string() +
          "-" +
          e1.end.to_string() +
          " exon2 " +
          e2.start.to_string() +
          "-" +
          e2.end.to_string()
        variants.push(v)
      }
      j = j + 1
    }
    i = i + 1
  }
  variants
}

///|
/// Detect Retained Intron (RI) events.
/// An intron with very low junction support relative to flanking exon coverage.
pub fn sgseq_detect_ri(
  junctions : Array[SGSeqJunction],
  exons : Array[SGSeqExon],
) -> Array[SGSeqVariant] {
  let variants = Array::new()
  for j in junctions {
    // Check if this junction has very low count relative to neighbors
    // An RI event: the intron (j.start+1 to j.end-1) is sometimes retained
    // Look for evidence of retention: coverage in the intron region
    // Simplified: if junction count is low and flanking exons exist
    let intron_start = j.start + 1
    let intron_end = j.end - 1
    // Check if any exon overlaps the intron (retention evidence)
    let mut retention_exon = false
    for e in exons {
      if e.chr == j.chr && e.start <= intron_end && e.end >= intron_start {
        retention_exon = true
      }
    }
    if retention_exon && j.count < 10 {
      let v = SGSeqVariant::new(
        chr=j.chr,
        start=j.start,
        end=j.end,
        strand=j.strand,
        variant_type=VariantRI,
      )
      v.exclusion_count = j.count
      v.inclusion_count = 50 // simplified: assume moderate retention
      v.psi = v.compute_psi()
      v.details = "intron " +
        j.start.to_string() +
        "-" +
        j.end.to_string() +
        " count " +
        j.count.to_string()
      variants.push(v)
    }
  }
  variants
}

///|
/// Detect all types of splice variants.
pub fn sgseq_detect_variants(
  junctions : Array[SGSeqJunction],
  exons : Array[SGSeqExon],
) -> Array[SGSeqVariant] {
  let variants = Array::new()
  for v in sgseq_detect_se(junctions, exons) {
    variants.push(v)
  }
  for v in sgseq_detect_a5ss(junctions) {
    variants.push(v)
  }
  for v in sgseq_detect_a3ss(junctions) {
    variants.push(v)
  }
  for v in sgseq_detect_mxe(junctions, exons) {
    variants.push(v)
  }
  for v in sgseq_detect_ri(junctions, exons) {
    variants.push(v)
  }
  variants
}

///|
/// Quantify PSI for a variant (recompute from junctions).
pub fn sgseq_quantify_psi(
  variant : SGSeqVariant,
  junctions : Array[SGSeqJunction],
) -> Double {
  let mut inc = 0
  let mut exc = 0
  for j in junctions {
    if j.chr != variant.chr {
      continue
    }
    match variant.variant_type {
      VariantSE => {
        // Inclusion: junctions to/from the exon
        if j.end == variant.start || j.start == variant.end {
          inc = inc + j.count
        }
        // Exclusion: junction spanning the exon
        if j.start < variant.start && j.end > variant.end {
          exc = exc + j.count
        }
      }
      _ => {
        // Simplified: use counts near variant
        if j.start >= variant.start && j.end <= variant.end {
          inc = inc + j.count
        }
      }
    }
  }
  let total = inc + exc
  if total == 0 {
    0.0
  } else {
    inc.to_double() / total.to_double()
  }
}

///|
/// Filter variants by minimum read support.
pub fn sgseq_filter_variants(
  variants : Array[SGSeqVariant],
  min_reads~ : Int = 10,
) -> Array[SGSeqVariant] {
  let result = Array::new()
  for v in variants {
    if v.inclusion_count + v.exclusion_count >= min_reads {
      result.push(v)
    }
  }
  result
}

///|
/// Generate sample junctions for testing.
/// Simulates a gene with 5 exons and various splicing events.
pub fn sgseq_sample_junctions() -> Array[SGSeqJunction] {
  [
    // Normal junctions (E1-E2, E2-E3, E3-E4, E4-E5)
    SGSeqJunction::new(chr="chr1", start=1000, end=2000, strand="+", count=100),
    SGSeqJunction::new(chr="chr1", start=2100, end=3000, strand="+", count=95),
    SGSeqJunction::new(chr="chr1", start=3100, end=4000, strand="+", count=90),
    SGSeqJunction::new(chr="chr1", start=4100, end=5000, strand="+", count=85),
    // SE event: skipping E3 (junction from E2 donor to E4 acceptor)
    SGSeqJunction::new(chr="chr1", start=2100, end=4000, strand="+", count=20),
    // A5SS event: alternative 5' site for E2 acceptor
    SGSeqJunction::new(chr="chr1", start=1500, end=2000, strand="+", count=15),
    // A3SS event: alternative 3' site for E3 donor
    SGSeqJunction::new(chr="chr1", start=3100, end=3500, strand="+", count=10),
  ]
}

///|
/// Generate sample exons for testing.
pub fn sgseq_sample_exons() -> Array[SGSeqExon] {
  [
    SGSeqExon::new(chr="chr1", start=1, end=1000, strand="+"),
    SGSeqExon::new(chr="chr1", start=2000, end=2100, strand="+"),
    SGSeqExon::new(chr="chr1", start=3000, end=3100, strand="+"),
    SGSeqExon::new(chr="chr1", start=4000, end=4100, strand="+"),
    SGSeqExon::new(chr="chr1", start=5000, end=6000, strand="+"),
  ]
}

///|
/// Format a variant as a string.
pub fn sgseq_variant_to_string(v : SGSeqVariant) -> String {
  v.type_name() +
    " " +
    v.chr +
    ":" +
    v.start.to_string() +
    "-" +
    v.end.to_string() +
    " (" +
    v.strand +
    ") inc=" +
    v.inclusion_count.to_string() +
    " exc=" +
    v.exclusion_count.to_string() +
    " PSI=" +
    v.psi.to_string() +
    " [" +
    v.details +
    "]"
}

///|
/// Summarize variants by type.
pub fn sgseq_summary(
  variants : Array[SGSeqVariant],
) -> Map[String, Int] {
  let counts : Map[String, Int] = Map([], capacity=5)
  for v in variants {
    let name = v.type_name()
    let cur = counts.get(name).or(0)
    counts[name] = cur + 1
  }
  counts
}
