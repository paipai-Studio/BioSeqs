///|
/// maftools - Cancer genomics MAF file parsing and mutation analysis.
/// Inspired by Bioconductor maftools.
/// Provides MAF file parsing, mutation classification (SNV/Indel/Complex),
/// tumor mutational burden (TMB) calculation, mutation spectrum analysis,
/// mutational signature decomposition, co-occurrence analysis, and oncoplot data generation.

///|
pub struct MAFMutation {
  hugo_symbol : String
  chromosome : String
  start_position : Int
  end_position : Int
  reference_allele : String
  tumor_seq_allele1 : String
  tumor_sample_barcode : String
  variant_classification : String
}

///|
pub struct MAFData {
  mutations : Array[MAFMutation]
  sample_barcodes : Array[String]
  gene_symbols : Array[String]
}

///|
pub struct MutationSpectrum {
  snv_count : Int
  indel_count : Int
  complex_count : Int
  transitions : Array[String]
  transversions : Array[String]
  transition_count : Int
  transversion_count : Int
}

///|
pub struct TMBResult {
  tmb : Double
  coding_region_size : Double
  total_mutations : Int
  per_sample_tmb : @hashmap.HashMap[String, Double]
}

///|
pub struct OncoplotData {
  genes : Array[String]
  samples : Array[String]
  mutation_matrix : Array[Array[String]]
}

///|
pub struct CoOccurrenceResult {
  pair : (String, String)
  co_occurrence : Int
  mutual_exclusivity : Int
  p_value : Double
}

///|
pub struct MutationalSignature {
  name : String
  contribution : Double
}

///|
pub fn MAFMutation::new(
  hugo_symbol : String,
  chromosome : String,
  start_position : Int,
  end_position : Int,
  reference_allele : String,
  tumor_seq_allele1 : String,
  tumor_sample_barcode : String,
  variant_classification : String
) -> MAFMutation {
  MAFMutation::{
    hugo_symbol,
    chromosome,
    start_position,
    end_position,
    reference_allele,
    tumor_seq_allele1,
    tumor_sample_barcode,
    variant_classification
  }
}

///|
pub fn MAFMutation::is_snv(self : MAFMutation) -> Bool {
  let ref_len = self.reference_allele.length()
  let alt_len = self.tumor_seq_allele1.length()
  ref_len == 1 && alt_len == 1
}

///|
pub fn MAFMutation::is_indel(self : MAFMutation) -> Bool {
  let ref_len = self.reference_allele.length()
  let alt_len = self.tumor_seq_allele1.length()
  if self.variant_classification.to_upper().contains("INSERTION") {
    return true
  }
  if self.variant_classification.to_upper().contains("DELETION") {
    return true
  }
  ref_len != alt_len
}

///|
pub fn MAFMutation::is_complex(self : MAFMutation) -> Bool {
  !self.is_snv() && !self.is_indel()
}

///|
pub fn MAFMutation::is_transition(self : MAFMutation) -> Bool {
  if !self.is_snv() {
    return false
  }
  let ref_base = self.reference_allele.to_upper()
  let alt_base = self.tumor_seq_allele1.to_upper()
  (ref_base == "A" && alt_base == "G") ||
  (ref_base == "G" && alt_base == "A") ||
  (ref_base == "C" && alt_base == "T") ||
  (ref_base == "T" && alt_base == "C")
}

///|
pub fn MAFMutation::is_transversion(self : MAFMutation) -> Bool {
  if !self.is_snv() {
    return false
  }
  !self.is_transition()
}

///|
pub fn MAFData::new() -> MAFData {
  MAFData::{
    mutations: Array::new(),
    sample_barcodes: Array::new(),
    gene_symbols: Array::new()
  }
}

///|
pub fn MAFData::add_mutation(self : MAFData, mutation : MAFMutation) -> MAFData {
  let new_mutations = self.mutations
  new_mutations.push(mutation)
  
  let new_samples = self.sample_barcodes
  if !maftools_contains_string(new_samples, mutation.tumor_sample_barcode) {
    new_samples.push(mutation.tumor_sample_barcode)
  }
  
  let new_genes = self.gene_symbols
  if !maftools_contains_string(new_genes, mutation.hugo_symbol) {
    new_genes.push(mutation.hugo_symbol)
  }
  
  MAFData::{
    mutations: new_mutations,
    sample_barcodes: new_samples,
    gene_symbols: new_genes
  }
}

///|
pub fn MAFData::count_mutations(self : MAFData) -> Int {
  self.mutations.length()
}

///|
pub fn MAFData::count_unique_genes(self : MAFData) -> Int {
  self.gene_symbols.length()
}

///|
pub fn MAFData::count_unique_samples(self : MAFData) -> Int {
  self.sample_barcodes.length()
}

///|
pub fn MAFData::filter_genes(self : MAFData, genes : Array[String]) -> MAFData {
  let filtered : Array[MAFMutation] = Array::new()
  let mut i = 0
  while i < self.mutations.length() {
    let mutation = self.mutations[i]
    if maftools_contains_string(genes, mutation.hugo_symbol) {
      filtered.push(mutation)
    }
    i = i + 1
  }
  
  let new_genes : Array[String] = Array::new()
  let mut j = 0
  while j < genes.length() {
    new_genes.push(genes[j])
    j = j + 1
  }
  
  MAFData::{
    mutations: filtered,
    sample_barcodes: self.sample_barcodes,
    gene_symbols: new_genes
  }
}

///|
pub fn MAFData::filter_samples(self : MAFData, samples : Array[String]) -> MAFData {
  let filtered : Array[MAFMutation] = Array::new()
  let mut i = 0
  while i < self.mutations.length() {
    let mutation = self.mutations[i]
    if maftools_contains_string(samples, mutation.tumor_sample_barcode) {
      filtered.push(mutation)
    }
    i = i + 1
  }
  
  let new_samples : Array[String] = Array::new()
  let mut j = 0
  while j < samples.length() {
    new_samples.push(samples[j])
    j = j + 1
  }
  
  MAFData::{
    mutations: filtered,
    sample_barcodes: new_samples,
    gene_symbols: self.gene_symbols
  }
}

///|
pub fn MAFData::filter_variant_type(self : MAFData, variant_type : String) -> MAFData {
  let filtered : Array[MAFMutation] = Array::new()
  let mut i = 0
  while i < self.mutations.length() {
    let mutation = self.mutations[i]
    if mutation.variant_classification.to_upper().contains(variant_type.to_upper()) {
      filtered.push(mutation)
    }
    i = i + 1
  }
  
  MAFData::{
    mutations: filtered,
    sample_barcodes: self.sample_barcodes,
    gene_symbols: self.gene_symbols
  }
}

///|
pub fn MAFData::get_mutation_gene(self : MAFData, index : Int) -> String {
  if index >= 0 && index < self.mutations.length() {
    self.mutations[index].hugo_symbol
  } else {
    ""
  }
}

///|
pub fn MAFData::get_mutation_sample(self : MAFData, index : Int) -> String {
  if index >= 0 && index < self.mutations.length() {
    self.mutations[index].tumor_sample_barcode
  } else {
    ""
  }
}

///|
pub fn get_gene_mutation_counts(maf : MAFData) -> @hashmap.HashMap[String, Int] {
  let counts : @hashmap.HashMap[String, Int] = @hashmap.HashMap([])
  let mut i = 0
  while i < maf.mutations.length() {
    let mutation = maf.mutations[i]
    let current_count = counts.get(mutation.hugo_symbol).unwrap_or(0)
    counts.set(mutation.hugo_symbol, current_count + 1)
    i = i + 1
  }
  counts
}

///|
pub fn get_sample_mutation_counts(maf : MAFData) -> @hashmap.HashMap[String, Int] {
  let counts : @hashmap.HashMap[String, Int] = @hashmap.HashMap([])
  let mut i = 0
  while i < maf.mutations.length() {
    let mutation = maf.mutations[i]
    let current_count = counts.get(mutation.tumor_sample_barcode).unwrap_or(0)
    counts.set(mutation.tumor_sample_barcode, current_count + 1)
    i = i + 1
  }
  counts
}

///|
pub fn calculate_mutation_spectrum(maf : MAFData) -> MutationSpectrum {
  let transitions : Array[String] = ["A>G", "G>A", "C>T", "T>C"]
  let transversions : Array[String] = ["A>C", "A>T", "C>A", "C>G", "G>C", "G>T"]
  
  let mut snv_count = 0
  let mut indel_count = 0
  let mut complex_count = 0
  let mut transition_count = 0
  let mut transversion_count = 0
  
  let mut i = 0
  while i < maf.mutations.length() {
    let mutation = maf.mutations[i]
    if mutation.is_snv() {
      snv_count = snv_count + 1
      if mutation.is_transition() {
        transition_count = transition_count + 1
      } else if mutation.is_transversion() {
        transversion_count = transversion_count + 1
      }
    } else if mutation.is_indel() {
      indel_count = indel_count + 1
    } else {
      complex_count = complex_count + 1
    }
    i = i + 1
  }
  
  MutationSpectrum::{
    snv_count,
    indel_count,
    complex_count,
    transitions,
    transversions,
    transition_count,
    transversion_count
  }
}

///|
pub fn calculate_tmb(maf : MAFData, coding_region_size : Double) -> TMBResult {
  let mutations = maf.mutations
  let mut total_mutations = 0
  let per_sample_tmb : @hashmap.HashMap[String, Double] = @hashmap.HashMap([])
  
  let mut i = 0
  while i < mutations.length() {
    let mutation = mutations[i]
    if mutation.is_snv() {
      total_mutations = total_mutations + 1
      let sample = mutation.tumor_sample_barcode
      let current_count = per_sample_tmb.get(sample).unwrap_or(0.0)
      per_sample_tmb.set(sample, current_count + 1.0)
    }
    i = i + 1
  }
  
  let result_map : @hashmap.HashMap[String, Double] = @hashmap.HashMap([])
  let sample_barcodes = maf.sample_barcodes
  let mut j = 0
  while j < sample_barcodes.length() {
    let sample = sample_barcodes[j]
    let mut_count = per_sample_tmb.get(sample).unwrap_or(0.0)
    let sample_tmb = mut_count / coding_region_size
    result_map.set(sample, sample_tmb)
    j = j + 1
  }
  
  let tmb = total_mutations.to_double() / coding_region_size
  
  TMBResult::{
    tmb,
    coding_region_size,
    total_mutations,
    per_sample_tmb: result_map
  }
}

///|
pub fn generate_oncoplot_data(maf : MAFData, top_n : Int) -> OncoplotData {
  let gene_counts = get_gene_mutation_counts(maf)
  let sorted_genes = maftools_sort_by_count_desc(gene_counts, top_n)
  let samples = maf.sample_barcodes
  
  let mutation_matrix : Array[Array[String]] = Array::new()
  
  let mut i = 0
  while i < sorted_genes.length() {
    let gene = sorted_genes[i]
    let row : Array[String] = Array::new()
    let mut j = 0
    while j < samples.length() {
      let mutation_type = maftools_get_mutation_type_for_sample_gene(maf, samples[j], gene)
      row.push(mutation_type)
      j = j + 1
    }
    mutation_matrix.push(row)
    i = i + 1
  }
  
  OncoplotData::{
    genes: sorted_genes,
    samples,
    mutation_matrix
  }
}

///|
pub fn analyze_co_occurrence(maf : MAFData, gene1 : String, gene2 : String) -> CoOccurrenceResult {
  let samples1 = maftools_get_samples_with_gene(maf, gene1)
  let samples2 = maftools_get_samples_with_gene(maf, gene2)
  
  let mut co_occurrence = 0
  let mut mutual_exclusivity = 0
  
  let all_samples = maf.sample_barcodes
  let mut i = 0
  while i < all_samples.length() {
    let sample = all_samples[i]
    let has_gene1 = maftools_contains_string(samples1, sample)
    let has_gene2 = maftools_contains_string(samples2, sample)
    
    if has_gene1 && has_gene2 {
      co_occurrence = co_occurrence + 1
    } else if has_gene1 != has_gene2 {
      mutual_exclusivity = mutual_exclusivity + 1
    }
    i = i + 1
  }
  
  let total = all_samples.length().to_double()
  let p_value = if total > 0.0 {
    1.0 - (co_occurrence.to_double() / total)
  } else {
    1.0
  }
  
  CoOccurrenceResult::{
    pair: (gene1, gene2),
    co_occurrence,
    mutual_exclusivity,
    p_value
  }
}

///|
pub fn decompose_mutational_signatures(maf : MAFData, reference_signatures : Array[MutationalSignature]) -> Array[MutationalSignature] {
  if reference_signatures.length() == 0 {
    return Array::new()
  }
  
  let spectrum = calculate_mutation_spectrum(maf)
  let total_mutations = spectrum.snv_count.to_double()
  
  if total_mutations == 0.0 {
    return Array::new()
  }
  
  let results : Array[MutationalSignature] = Array::new()
  let mut i = 0
  while i < reference_signatures.length() {
    let sig = reference_signatures[i]
    let contribution = sig.contribution * total_mutations
    results.push(MutationalSignature::{
      name: sig.name,
      contribution
    })
    i = i + 1
  }
  
  results
}

///|
pub fn summarize_maf(maf : MAFData) -> String {
  let mut summary = "\nMAF Summary\n"
  summary = summary + "==============================\n"
  summary = summary + "Total mutations: " + maf.count_mutations().to_string() + "\n"
  summary = summary + "Total genes: " + maf.count_unique_genes().to_string() + "\n"
  summary = summary + "Total samples: " + maf.count_unique_samples().to_string() + "\n"
  
  let spectrum = calculate_mutation_spectrum(maf)
  summary = summary + "SNVs: " + spectrum.snv_count.to_string() + "\n"
  summary = summary + "Indels: " + spectrum.indel_count.to_string() + "\n"
  summary = summary + "Complex: " + spectrum.complex_count.to_string() + "\n"
  
  let tmb = calculate_tmb(maf, 3.0e7)
  summary = summary + "TMB: " + tmb.tmb.to_string() + " mut/Mb\n"
  
  summary
}

///|
pub fn create_example_maf() -> MAFData {
  let mut maf = MAFData::new()
  
  let genes = ["TP53", "BRCA1", "EGFR", "KRAS", "MYC", "PTEN", "APC", "RB1", "VHL", "NF1"]
  let samples = ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6", "Sample7", "Sample8"]
  let variants = ["Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "Splice_Site"]
  
  let mut i = 0
  while i < 80 {
    let gene_idx = i % 10
    let sample_idx = i % 8
    let variant_idx = i % 5
    let pos = 1000000 + i * 10000
    
    let mutation = MAFMutation::new(
      genes[gene_idx],
      "chr" + ((i % 22) + 1).to_string(),
      pos,
      pos,
      if i % 3 == 0 { "A" } else { if i % 3 == 1 { "C" } else { "G" } },
      if i % 2 == 0 { "T" } else { if i % 3 == 0 { "A" } else { "C" } },
      samples[sample_idx],
      variants[variant_idx]
    )
    
    maf = maf.add_mutation(mutation)
    i = i + 1
  }
  
  maf
}

///|
pub fn parse_maf_content(content : String) -> MAFData {
  let lines = split_by_char(content, '\n'.to_int())
  let header = if lines.length() > 0 { split_by_char(lines[0], '\t'.to_int()) } else { Array::new() }
  
  let hugo_idx = maftools_find_column_index(header, "Hugo_Symbol")
  let chrom_idx = maftools_find_column_index(header, "Chromosome")
  let start_idx = maftools_find_column_index(header, "Start_Position")
  let end_idx = maftools_find_column_index(header, "End_Position")
  let ref_idx = maftools_find_column_index(header, "Reference_Allele")
  let alt_idx = maftools_find_column_index(header, "Tumor_Seq_Allele1")
  let sample_idx = maftools_find_column_index(header, "Tumor_Sample_Barcode")
  let variant_idx = maftools_find_column_index(header, "Variant_Classification")
  
  let mut maf = MAFData::new()
  let mut i = 1
  while i < lines.length() {
    let line = lines[i]
    if line.length() > 0 {
      let fields = split_by_char(line, '\t'.to_int())
      if fields.length() >= 8 {
        let mutation = MAFMutation::new(
          if hugo_idx >= 0 { fields[hugo_idx] } else { "" },
          if chrom_idx >= 0 { fields[chrom_idx] } else { "" },
          if start_idx >= 0 { maftools_parse_int(fields[start_idx]) } else { 0 },
          if end_idx >= 0 { maftools_parse_int(fields[end_idx]) } else { 0 },
          if ref_idx >= 0 { fields[ref_idx] } else { "" },
          if alt_idx >= 0 { fields[alt_idx] } else { "" },
          if sample_idx >= 0 { fields[sample_idx] } else { "" },
          if variant_idx >= 0 { fields[variant_idx] } else { "" }
        )
        maf = maf.add_mutation(mutation)
      }
    }
    i = i + 1
  }
  
  maf
}

// Helper functions

fn maftools_contains_string(arr : Array[String], value : String) -> Bool {
  let mut i = 0
  while i < arr.length() {
    if arr[i] == value {
      return true
    }
    i = i + 1
  }
  false
}

fn maftools_find_column_index(header : Array[String], column_name : String) -> Int {
  let mut i = 0
  while i < header.length() {
    if header[i].to_upper() == column_name.to_upper() {
      return i
    }
    i = i + 1
  }
  -1
}

fn maftools_parse_int(s : String) -> Int {
  parse_int(s)
}

fn maftools_sort_by_count_desc(counts : @hashmap.HashMap[String, Int], top_n : Int) -> Array[String] {
  let entries = counts.iter().collect()
  let sorted : Array[(String, Int)] = Array::new()
  
  let mut i = 0
  while i < entries.length() {
    let (key, value) = entries[i]
    sorted.push((key, value))
    i = i + 1
  }
  
  // Simple bubble sort by count descending
  let len = sorted.length()
  let mut j = 0
  while j < len - 1 {
    let mut swapped = false
    let mut k = 0
    while k < len - j - 1 {
      if sorted[k].1 < sorted[k + 1].1 {
        let temp = sorted[k]
        sorted[k] = sorted[k + 1]
        sorted[k + 1] = temp
        swapped = true
      }
      k = k + 1
    }
    if !swapped {
      break
    }
    j = j + 1
  }
  
  let result : Array[String] = Array::new()
  let count = if top_n < sorted.length() { top_n } else { sorted.length() }
  let mut l = 0
  while l < count {
    result.push(sorted[l].0)
    l = l + 1
  }
  
  result
}

fn maftools_get_mutation_type_for_sample_gene(maf : MAFData, sample : String, gene : String) -> String {
  let mut i = 0
  while i < maf.mutations.length() {
    let mutation = maf.mutations[i]
    if mutation.hugo_symbol == gene && mutation.tumor_sample_barcode == sample {
      return mutation.variant_classification
    }
    i = i + 1
  }
  ""
}

fn maftools_get_samples_with_gene(maf : MAFData, gene : String) -> Array[String] {
  let samples : Array[String] = Array::new()
  let mut i = 0
  while i < maf.mutations.length() {
    let mutation = maf.mutations[i]
    if mutation.hugo_symbol == gene {
      if !maftools_contains_string(samples, mutation.tumor_sample_barcode) {
        samples.push(mutation.tumor_sample_barcode)
      }
    }
    i = i + 1
  }
  samples
}
