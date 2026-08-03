///|
/// Bio.PopGen.GenePop - GenePop file format parsing and analysis.
/// Mirrors Biopython's Bio.PopGen.GenePop module.
///
/// The GenePop file format stores multi-locus genotype data for one or more
/// populations. The structure is:
///
/// ```
/// <optional title line>
/// Pop
/// <individual_name>, <locus1> <locus2> ... <locusN>
/// <individual_name>, <locus1> <locus2> ... <locusN>
/// Pop
/// <individual_name>, <locus1> <locus2> ... <locusN>
/// ```
///
/// Each locus entry is either:
///   * diploid 2-digit: a 4-character string "AABB" (allele1=AA, allele2=BB)
///   * diploid 3-digit: a 6-character string "AABBBB" (allele1=AAA, allele2=BBB)
///   * haploid: a single allele (no second allele)
///   * missing: "0000" (2-digit) or "000000" (3-digit) or "0"
///
/// The locus digit size is auto-detected from the first valid entry.

// ============================================================================
// Data structures
// ============================================================================

///|
/// A single genotype at one locus for one individual.
pub struct GenePopGenotype {
  /// First allele code (0 = missing).
  allele1 : Int
  /// Second allele code; None for haploid loci.
  allele2 : Int?
  /// True if this genotype represents missing data.
  is_missing : Bool
}

///|
pub fn GenePopGenotype::new(allele1 : Int, allele2 : Int?) -> GenePopGenotype {
  let missing = allele1 == 0 && (allele2 == Some(0) || allele2 == None)
  GenePopGenotype::{ allele1, allele2, is_missing: missing }
}

///|
/// A diploid genotype (both alleles present).
pub fn GenePopGenotype::diploid(allele1 : Int, allele2 : Int) -> GenePopGenotype {
  let missing = allele1 == 0 && allele2 == 0
  GenePopGenotype::{ allele1, allele2: Some(allele2), is_missing: missing }
}

///|
/// A haploid genotype (only one allele).
pub fn GenePopGenotype::haploid(allele1 : Int) -> GenePopGenotype {
  GenePopGenotype::{ allele1, allele2: None, is_missing: allele1 == 0 }
}

///|
/// An individual with a name and a genotype per locus.
pub struct GenePopIndividual {
  name : String
  genotypes : Array[GenePopGenotype]
}

///|
pub fn GenePopIndividual::new(
  name : String,
  genotypes : Array[GenePopGenotype],
) -> GenePopIndividual {
  GenePopIndividual::{ name, genotypes }
}

///|
/// A population: a named group of individuals. The name is auto-assigned
/// (e.g. "Pop_1") when not explicitly given in the file.
pub struct GenePopPopulation {
  name : String
  individuals : Array[GenePopIndividual]
}

///|
pub fn GenePopPopulation::new(
  name : String,
  individuals : Array[GenePopIndividual],
) -> GenePopPopulation {
  GenePopPopulation::{ name, individuals }
}

///|
/// A complete GenePop record: title, locus names, and populations.
pub struct GenePopRecord {
  title : String
  loci_names : Array[String]
  populations : Array[GenePopPopulation]
  /// Number of digits per allele (2 or 3), 0 if unknown.
  locus_digit_size : Int
}

///|
pub fn GenePopRecord::new(
  title : String,
  loci_names : Array[String],
  populations : Array[GenePopPopulation],
  locus_digit_size : Int,
) -> GenePopRecord {
  GenePopRecord::{
    title,
    loci_names,
    populations,
    locus_digit_size,
  }
}

// ============================================================================
// Internal helpers
// ============================================================================

///|
fn gp_starts_with(s : String, prefix : String) -> Bool {
  let n = s.length()
  let p = prefix.length()
  if p > n {
    return false
  }
  let mut i = 0
  while i < p {
    if s.unsafe_get(i).to_int() != prefix.unsafe_get(i).to_int() {
      return false
    }
    i = i + 1
  }
  true
}

///|
/// Case-insensitive string equality for short ASCII strings.
fn gp_iequals(a : String, b : String) -> Bool {
  let na = a.length()
  let nb = b.length()
  if na != nb {
    return false
  }
  let mut i = 0
  while i < na {
    let ca = a.unsafe_get(i).to_int()
    let cb = b.unsafe_get(i).to_int()
    let la = if ca >= 'A'.to_int() && ca <= 'Z'.to_int() {
      ca + 32
    } else {
      ca
    }
    let lb = if cb >= 'A'.to_int() && cb <= 'Z'.to_int() {
      cb + 32
    } else {
      cb
    }
    if la != lb {
      return false
    }
    i = i + 1
  }
  true
}

///|
/// Trim leading and trailing ASCII whitespace.
fn gp_trim(s : String) -> String {
  let n = s.length()
  let mut start = 0
  let mut end = n
  while start < end {
    let c = s.unsafe_get(start).to_int()
    if c == ' '.to_int() || c == '\t'.to_int() || c == '\r'.to_int() {
      start = start + 1
    } else {
      break
    }
  }
  while end > start {
    let c = s.unsafe_get(end - 1).to_int()
    if c == ' '.to_int() || c == '\t'.to_int() || c == '\r'.to_int() {
      end = end - 1
    } else {
      break
    }
  }
  s[start:end].to_string()
}

///|
/// Split a line by runs of ASCII whitespace, dropping empty tokens.
fn gp_split_whitespace(s : String) -> Array[String] {
  let n = s.length()
  let tokens : Array[String] = []
  let mut start = 0
  let mut i = 0
  let mut in_token = false
  while i < n {
    let c = s.unsafe_get(i).to_int()
    let is_ws = c == ' '.to_int() || c == '\t'.to_int() || c == '\r'.to_int()
    if is_ws {
      if in_token {
        tokens.push(s[start:i].to_string())
        in_token = false
      }
    } else if !in_token {
      start = i
      in_token = true
    }
    i = i + 1
  }
  if in_token {
    tokens.push(s[start:n].to_string())
  }
  tokens
}

///|
/// Parse a single locus token (e.g. "0102", "001002", "01") into a genotype.
/// Auto-detects 2-digit vs 3-digit based on token length.
fn gp_parse_locus_token(token : String) -> GenePopGenotype {
  let n = token.length()
  // Missing data shortcuts.
  if n == 1 && token.unsafe_get(0).to_int() == '0'.to_int() {
    return GenePopGenotype::haploid(0)
  }
  if n == 4 && token == "0000" {
    return GenePopGenotype::diploid(0, 0)
  }
  if n == 6 && token == "000000" {
    return GenePopGenotype::diploid(0, 0)
  }
  // 3-digit diploid: "AAABBB"
  if n == 6 {
    let a1_str = token[0:3].to_string()
    let a2_str = token[3:6].to_string()
    let a1 = gp_parse_allele(a1_str)
    let a2 = gp_parse_allele(a2_str)
    return GenePopGenotype::diploid(a1, a2)
  }
  // 2-digit diploid: "AABB"
  if n == 4 {
    let a1_str = token[0:2].to_string()
    let a2_str = token[2:4].to_string()
    let a1 = gp_parse_allele(a1_str)
    let a2 = gp_parse_allele(a2_str)
    return GenePopGenotype::diploid(a1, a2)
  }
  // Haploid: single allele "AA" or "AAA"
  if n == 2 || n == 3 {
    let a1 = gp_parse_allele(token)
    return GenePopGenotype::haploid(a1)
  }
  // Unknown format: treat as missing.
  GenePopGenotype::diploid(0, 0)
}

///|
/// Parse an allele string (e.g. "01", "001", "15", "999") into an Int.
/// Returns 0 for empty/unparseable strings.
fn gp_parse_allele(s : String) -> Int {
  let n = s.length()
  if n == 0 {
    return 0
  }
  let mut result = 0
  let mut i = 0
  let zero = '0'.to_int()
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c < zero || c > zero + 9 {
      // Non-digit encountered; return what we have so far.
      return result
    }
    result = result * 10 + (c - zero)
    i = i + 1
  }
  result
}

///|
/// Determine the per-allele digit size (2 or 3) from a locus token.
/// Returns 0 if it cannot be determined.
fn gp_detect_digit_size(token : String) -> Int {
  let n = token.length()
  if n == 4 {
    return 2
  }
  if n == 6 {
    return 3
  }
  0
}

///|
/// Format an allele as a zero-padded string of the given digit size.
fn gp_format_allele(allele : Int, digit_size : Int) -> String {
  if digit_size == 3 {
    if allele < 10 {
      return "00" + allele.to_string()
    }
    if allele < 100 {
      return "0" + allele.to_string()
    }
    return allele.to_string()
  }
  // Default: 2-digit.
  if allele < 10 {
    return "0" + allele.to_string()
  }
  allele.to_string()
}

// ============================================================================
// Parsing
// ============================================================================

///|
/// Parse a GenePop file into a GenePopRecord.
///
/// The parser:
///   1. Reads the first non-empty line as the title (optional).
///   2. Detects "Pop" markers (case-insensitive) as population separators.
///   3. Parses each individual line as "name, genotype1 genotype2 ...".
///   4. Auto-detects the allele digit size (2 or 3) from the first genotype.
///   5. Locus names are NOT present in the GenePop file format; they are
///      auto-generated as "Locus_1", "Locus_2", ... unless a comment line
///      beginning with "#Loci:" provides them.
pub fn parse_genepop(content : String) -> GenePopRecord {
  let lines = split_by_char(content, '\n'.to_int())

  let mut title = ""
  let mut loci_names : Array[String] = []
  let populations : Array[GenePopPopulation] = []
  let mut cur_individuals : Array[GenePopIndividual] = []
  let mut digit_size = 0
  let mut pop_count = 0
  let mut has_title = false

  let mut i = 0
  while i < lines.length() {
    let raw = lines[i]
    let trimmed = gp_trim(raw)
    if trimmed.length() == 0 {
      i = i + 1
      continue
    }
    // Check for locus name comment line: "#Loci: L1 L2 L3"
    if gp_starts_with(trimmed, "#Loci:") {
      let rest = trimmed[6:trimmed.length()].to_string()
      loci_names = gp_split_whitespace(rest)
      i = i + 1
      continue
    }
    // Skip other comment lines.
    if gp_starts_with(trimmed, "#") {
      i = i + 1
      continue
    }
    // Check for "Pop" marker (case-insensitive).
    if gp_iequals(trimmed, "Pop") {
      // Flush current population if it has individuals.
      if cur_individuals.length() > 0 {
        pop_count = pop_count + 1
        populations.push(
          GenePopPopulation::new("Pop_" + pop_count.to_string(), cur_individuals),
        )
        cur_individuals = []
      } else if populations.length() == 0 {
        // First "Pop" marker with no preceding individuals; start a new pop.
        pop_count = pop_count + 1
      }
      i = i + 1
      continue
    }
    // If no title yet and no "Pop" seen, this line is the title.
    if !has_title && pop_count == 0 {
      title = trimmed
      has_title = true
      i = i + 1
      continue
    }
    // Otherwise, this is an individual line: "name, g1 g2 g3 ..."
    let ind = gp_parse_individual(trimmed)
    if ind.genotypes.length() > 0 && digit_size == 0 {
      // Detect digit size from the first genotype token's length.
      let comma_pos = gp_find_char(trimmed, ',')
      let geno_str = if comma_pos >= 0 {
        gp_trim(trimmed[comma_pos + 1:trimmed.length()].to_string())
      } else {
        trimmed
      }
      let tokens = gp_split_whitespace(geno_str)
      if tokens.length() > 0 {
        digit_size = gp_detect_digit_size(tokens[0])
      }
      if digit_size == 0 {
        digit_size = gp_detect_digit_size_from_ind(ind)
      }
    }
    cur_individuals.push(ind)
    i = i + 1
  }

  // Flush the last population.
  if cur_individuals.length() > 0 {
    pop_count = pop_count + 1
    populations.push(
      GenePopPopulation::new("Pop_" + pop_count.to_string(), cur_individuals),
    )
  }

  // Auto-generate locus names if none were provided.
  if loci_names.length() == 0 && populations.length() > 0 {
    let n_loci = gp_count_loci(populations)
    let mut k = 0
    while k < n_loci {
      loci_names.push("Locus_" + (k + 1).to_string())
      k = k + 1
    }
  }

  GenePopRecord::new(title, loci_names, populations, digit_size)
}

///|
/// Detect allele digit size from the first genotype in an individual.
fn gp_detect_digit_size_from_ind(ind : GenePopIndividual) -> Int {
  if ind.genotypes.length() == 0 {
    return 0
  }
  // We can't recover the original token length from the parsed genotype,
  // so we use the magnitude of the alleles as a heuristic:
  //  * If any allele > 99, it must be 3-digit.
  //  * If all alleles <= 99, assume 2-digit (the common default).
  let mut max_allele = 0
  let mut k = 0
  while k < ind.genotypes.length() {
    let g = ind.genotypes[k]
    if g.allele1 > max_allele {
      max_allele = g.allele1
    }
    match g.allele2 {
      Some(a2) =>
        if a2 > max_allele {
          max_allele = a2
        }
      None => ()
    }
    k = k + 1
  }
  if max_allele > 99 {
    return 3
  }
  2
}

///|
/// Parse an individual line: "name, g1 g2 g3 ..."
fn gp_parse_individual(line : String) -> GenePopIndividual {
  // Split on the first comma to separate name from genotype data.
  let comma_pos = gp_find_char(line, ',')
  let (name_str, geno_str) = if comma_pos >= 0 {
    (gp_trim(line[0:comma_pos].to_string()), gp_trim(line[comma_pos + 1:line.length()].to_string()))
  } else {
    // No comma: treat the whole line as genotype data with empty name.
    ("", line)
  }
  let tokens = gp_split_whitespace(geno_str)
  let genotypes : Array[GenePopGenotype] = []
  let mut t = 0
  while t < tokens.length() {
    genotypes.push(gp_parse_locus_token(tokens[t]))
    t = t + 1
  }
  GenePopIndividual::new(name_str, genotypes)
}

///|
/// Find the first occurrence of a character in a string. Returns -1 if not found.
fn gp_find_char(s : String, ch : Char) -> Int {
  let target = ch.to_int()
  let n = s.length()
  let mut i = 0
  while i < n {
    if s.unsafe_get(i).to_int() == target {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Count the number of loci from the first individual in the first population.
fn gp_count_loci(populations : Array[GenePopPopulation]) -> Int {
  if populations.length() == 0 {
    return 0
  }
  let pop = populations[0]
  if pop.individuals.length() == 0 {
    return 0
  }
  pop.individuals[0].genotypes.length()
}

// ============================================================================
// Serialization
// ============================================================================

///|
/// Serialize a GenePopRecord back to GenePop file format text.
pub fn genepop_to_string(record : GenePopRecord) -> String {
  let buf = StringBuilder::new()
  // Title line.
  if record.title.length() > 0 {
    buf.write_string(record.title)
    buf.write_char('\n')
  }
  let digit = if record.locus_digit_size > 0 {
    record.locus_digit_size
  } else {
    2
  }
  let mut pi = 0
  while pi < record.populations.length() {
    let pop = record.populations[pi]
    buf.write_string("Pop\n")
    let mut ii = 0
    while ii < pop.individuals.length() {
      let ind = pop.individuals[ii]
      buf.write_string(ind.name)
      buf.write_char(',')
      let mut gi = 0
      while gi < ind.genotypes.length() {
        buf.write_char(' ')
        let g = ind.genotypes[gi]
        if g.is_missing {
          if digit == 3 {
            buf.write_string("000000")
          } else {
            buf.write_string("0000")
          }
        } else {
          buf.write_string(gp_format_allele(g.allele1, digit))
          match g.allele2 {
            Some(a2) => buf.write_string(gp_format_allele(a2, digit))
            None => ()
          }
        }
        gi = gi + 1
      }
      buf.write_char('\n')
      ii = ii + 1
    }
    pi = pi + 1
  }
  buf.to_string()
}

// ============================================================================
// Analysis utilities
// ============================================================================

///|
/// Count the total number of individuals across all populations.
pub fn genepop_total_individuals(record : GenePopRecord) -> Int {
  let mut total = 0
  let mut i = 0
  while i < record.populations.length() {
    total = total + record.populations[i].individuals.length()
    i = i + 1
  }
  total
}

///|
/// Count the number of populations in the record.
pub fn genepop_num_populations(record : GenePopRecord) -> Int {
  record.populations.length()
}

///|
/// Count the number of loci.
pub fn genepop_num_loci(record : GenePopRecord) -> Int {
  record.loci_names.length()
}

///|
/// Collect all distinct alleles observed at a given locus index across
/// all populations. Returns a sorted array of allele codes.
pub fn genepop_alleles_at_locus(
  record : GenePopRecord,
  locus_index : Int,
) -> Array[Int] {
  let seen : Map[Int, Unit] = Map([], capacity=0)
  let mut i = 0
  while i < record.populations.length() {
    let pop = record.populations[i]
    let mut j = 0
    while j < pop.individuals.length() {
      let ind = pop.individuals[j]
      if locus_index < ind.genotypes.length() {
        let g = ind.genotypes[locus_index]
        if !g.is_missing {
          seen.set(g.allele1, ())
          match g.allele2 {
            Some(a2) => seen.set(a2, ())
            None => ()
          }
        }
      }
      j = j + 1
    }
    i = i + 1
  }
  let result : Array[Int] = []
  let keys = seen.keys()
  for k in keys {
    result.push(k)
  }
  // Bubble sort (MoonBit's sort has no comparator argument).
  let n = result.length()
  let mut a = 0
  while a < n - 1 {
    let mut b = a + 1
    while b < n {
      if result[a] > result[b] {
        let tmp = result[a]
        result[a] = result[b]
        result[b] = tmp
      }
      b = b + 1
    }
    a = a + 1
  }
  result
}

///|
/// Compute allele frequencies at a locus within a single population.
/// Returns a Map from allele code to frequency (0.0-1.0). Missing
/// genotypes are excluded from the denominator.
pub fn genepop_allele_frequencies(
  population : GenePopPopulation,
  locus_index : Int,
) -> Map[Int, Double] {
  let counts : Map[Int, Int] = Map([], capacity=0)
  let mut total = 0
  let mut i = 0
  while i < population.individuals.length() {
    let ind = population.individuals[i]
    if locus_index < ind.genotypes.length() {
      let g = ind.genotypes[locus_index]
      if !g.is_missing {
        let c1 = counts.get_or_default(g.allele1, 0)
        counts.set(g.allele1, c1 + 1)
        total = total + 1
        match g.allele2 {
          Some(a2) => {
            let c2 = counts.get_or_default(a2, 0)
            counts.set(a2, c2 + 1)
            total = total + 1
          }
          None => ()
        }
      }
    }
    i = i + 1
  }
  let freqs : Map[Int, Double] = Map([], capacity=0)
  if total == 0 {
    return freqs
  }
  let total_d = total.to_double()
  let keys_iter = counts.keys()
  for k in keys_iter {
    let c = counts.get_or_default(k, 0)
    freqs.set(k, c.to_double() / total_d)
  }
  freqs
}

///|
/// Compute the observed heterozygosity at a locus within a population.
/// H_obs = (number of heterozygotes) / (number of non-missing genotypes).
pub fn genepop_observed_heterozygosity(
  population : GenePopPopulation,
  locus_index : Int,
) -> Double {
  let mut het = 0
  let mut total = 0
  let mut i = 0
  while i < population.individuals.length() {
    let ind = population.individuals[i]
    if locus_index < ind.genotypes.length() {
      let g = ind.genotypes[locus_index]
      if !g.is_missing {
        total = total + 1
        match g.allele2 {
          Some(a2) => if a2 != g.allele1 {
            het = het + 1
          }
          None => ()
        }
      }
    }
    i = i + 1
  }
  if total == 0 {
    return 0.0
  }
  het.to_double() / total.to_double()
}

///|
/// Compute the expected heterozygosity (gene diversity) at a locus
/// within a population. H_exp = 1 - sum(p_i^2) where p_i is the
/// frequency of allele i.
pub fn genepop_expected_heterozygosity(
  population : GenePopPopulation,
  locus_index : Int,
) -> Double {
  let freqs = genepop_allele_frequencies(population, locus_index)
  let mut sum_sq = 0.0
  let keys_iter = freqs.keys()
  for k in keys_iter {
    let p = freqs.get_or_default(k, 0.0)
    sum_sq = sum_sq + p * p
  }
  1.0 - sum_sq
}

///|
/// Return a short human-readable summary of the record.
pub fn genepop_summary(record : GenePopRecord) -> String {
  let buf = StringBuilder::new()
  buf.write_string("GenePop Record Summary\n")
  buf.write_string("  Title: " + record.title + "\n")
  buf.write_string("  Populations: " + record.populations.length().to_string() + "\n")
  buf.write_string("  Loci: " + record.loci_names.length().to_string() + "\n")
  buf.write_string("  Individuals: " + genepop_total_individuals(record).to_string() + "\n")
  buf.write_string("  Digit size: " + record.locus_digit_size.to_string() + "\n")
  let mut i = 0
  while i < record.populations.length() {
    let pop = record.populations[i]
    buf.write_string(
      "  " +
      pop.name +
      ": " +
      pop.individuals.length().to_string() +
      " individuals\n",
    )
    i = i + 1
  }
  buf.to_string()
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Generate a small GenePop-format sample with 2 populations, 3 loci,
/// and 2-digit allele encoding, for tests and demos.
pub fn genepop_sample() -> String {
  "Sample GenePop File\n" +
    "Pop\n" +
    "ind1,  0102  0304  0101\n" +
    "ind2,  0202  0303  0102\n" +
    "ind3,  0101  0404  0202\n" +
    "Pop\n" +
    "ind4,  0102  0304  0101\n" +
    "ind5,  0202  0303  0102\n" +
    "ind6,  0101  0000  0202\n"
}
