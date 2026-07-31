///|
/// RaggedExperiment module - Ragged mutation data arrays.
/// Inspired by Bioconductor RaggedExperiment package.
///
/// References:
///   - https://www.bioconductor.org/packages/release/bioc/html/RaggedExperiment.html
///   - Ramos AH et al. (2020) RaggedExperiment: an approach for integrating
///     mutation and copy number alteration data in object-oriented systems.
///     Cancer Research 80(16 Suppl):Abstract nr 3668.

///|
/// Mutation type enumeration.
pub enum MutationType {
  Missense_Mutation
  Nonsense_Mutation
  Frame_Shift_Ins
  Frame_Shift_Del
  In_Frame_Ins
  In_Frame_Del
  Splice_Site
  Translation_Start_Site
  Translation_Stop_Site
  Multi_Hit_Mutation
  Silent_Mutation
  Other
} derive(Eq, Debug)

///|
/// Create a missense mutation type.
pub fn mutation_missense() -> MutationType { MutationType::Missense_Mutation }
///|
/// Create a nonsense mutation type.
pub fn mutation_nonsense() -> MutationType { MutationType::Nonsense_Mutation }
///|
/// Create a frameshift insertion mutation type.
pub fn mutation_fs_ins() -> MutationType { MutationType::Frame_Shift_Ins }
///|
/// Create a frameshift deletion mutation type.
pub fn mutation_fs_del() -> MutationType { MutationType::Frame_Shift_Del }
///|
/// Create a splice site mutation type.
pub fn mutation_splice_site() -> MutationType { MutationType::Splice_Site }
///|
/// Create a silent mutation type.
pub fn mutation_silent() -> MutationType { MutationType::Silent_Mutation }

///|
/// Convert mutation type to string.
pub fn MutationType::to_string(self : MutationType) -> String {
  match self {
    MutationType::Missense_Mutation => "Missense_Mutation"
    MutationType::Nonsense_Mutation => "Nonsense_Mutation"
    MutationType::Frame_Shift_Ins => "Frame_Shift_Ins"
    MutationType::Frame_Shift_Del => "Frame_Shift_Del"
    MutationType::In_Frame_Ins => "In_Frame_Ins"
    MutationType::In_Frame_Del => "In_Frame_Del"
    MutationType::Splice_Site => "Splice_Site"
    MutationType::Translation_Start_Site => "Translation_Start_Site"
    MutationType::Translation_Stop_Site => "Translation_Stop_Site"
    MutationType::Multi_Hit_Mutation => "Multi_Hit_Mutation"
    MutationType::Silent_Mutation => "Silent_Mutation"
    MutationType::Other => "Other"
  }
}

///|
/// Mutation record.
pub struct MutationRecord {
  /// Sample/barcode ID
  sample_id : String
  /// Gene symbol
  gene_symbol : String
  /// Chromosome
  chrom : String
  /// Position
  pos : Int
  /// Reference allele
  ref_allele : String
  /// Alternate allele
  alt_allele : String
  /// Mutation type
  mutation_type : MutationType
  /// Variant allele frequency
  variant_allele_frequency : Double
  /// Read count supporting variant
  t_var_count : Int
  /// Total read count
  t_ref_count : Int
  /// Mean tumor allele frequency
  m_var_count : Int
  /// Mean reference count
  m_ref_count : Int
  /// Scored for damaging effect
  is_damaging : Bool
  /// Sequencer
  sequencer : String
  /// Validation status
  validation : String
}

///|
/// Create a new MutationRecord.
pub fn MutationRecord::new(
  sample_id~ : String,
  gene_symbol~ : String,
  chrom~ : String,
  pos~ : Int,
  ref_allele~ : String,
  alt_allele~ : String,
  mutation_type? : MutationType = MutationType::Missense_Mutation,
  variant_allele_frequency? : Double = 0.0,
) -> MutationRecord {
  MutationRecord::{
    sample_id,
    gene_symbol,
    chrom,
    pos,
    ref_allele,
    alt_allele,
    mutation_type,
    variant_allele_frequency,
    t_var_count: 0,
    t_ref_count: 0,
    m_var_count: 0,
    m_ref_count: 0,
    is_damaging: false,
    sequencer: "",
    validation: "",
  }
}

///|
/// Get the mutation ID (gene + chrom + pos + ref + alt).
pub fn MutationRecord::mutation_id(self : MutationRecord) -> String {
  self.gene_symbol + ":" + self.chrom + ":" + self.pos.to_string() + ":" + self.ref_allele + ">" + self.alt_allele
}

///|
/// Get the mutation key for grouping.
pub fn MutationRecord::mutation_key(self : MutationRecord) -> String {
  self.chrom + ":" + self.pos.to_string() + ":" + self.ref_allele + ">" + self.alt_allele
}

///|
/// Check if this is a loss-of-function mutation.
pub fn MutationRecord::is_loss_of_function(self : MutationRecord) -> Bool {
  match self.mutation_type {
    MutationType::Nonsense_Mutation => true
    MutationType::Frame_Shift_Ins => true
    MutationType::Frame_Shift_Del => true
    MutationType::Splice_Site => true
    MutationType::Translation_Start_Site => true
    _ => false
  }
}

///|
/// RaggedExperiment data structure.
/// Stores ragged arrays where each row (gene) can have different numbers
/// of mutations per sample.
pub struct RaggedExperiment {
  /// Row names (gene symbols)
  rownames : Array[String]
  /// Column names (sample barcodes)
  colnames : Array[String]
  /// Data: ragged array of mutation records per gene per sample
  /// data[row_idx][col_idx] = array of MutationRecord
  data : Array[Array[Array[MutationRecord]]]
  /// Number of mutations per gene per sample (cached)
  mut counts : Array[Array[Int]]
  /// Total mutation counts per sample
  mut sample_tmb : Array[Double]
}

///|
/// Create an empty RaggedExperiment.
pub fn RaggedExperiment::new() -> RaggedExperiment {
  RaggedExperiment::{
    rownames: Array::new(),
    colnames: Array::new(),
    data: Array::new(),
    counts: Array::new(),
    sample_tmb: Array::new(),
  }
}

///|
/// Create a RaggedExperiment from mutation records.
pub fn RaggedExperiment::from_records(records~ : Array[MutationRecord]) -> RaggedExperiment {
  let exp = RaggedExperiment::new()
  let mut i = 0
  while i < records.length() {
    exp.add_record(records[i])
    i = i + 1
  }
  exp.recalculate()
  exp
}

///|
/// Add a mutation record to the experiment.
pub fn RaggedExperiment::add_record(self : RaggedExperiment, record : MutationRecord) -> Unit {
  let mut row_idx = -1
  let mut i = 0
  while i < self.rownames.length() {
    if self.rownames[i] == record.gene_symbol {
      row_idx = i
      break
    }
    i = i + 1
  }

  if row_idx < 0 {
    self.rownames.push(record.gene_symbol)
    let new_row : Array[Array[MutationRecord]] = Array::new()
    let mut c = 0
    while c < self.colnames.length() {
      new_row.push(Array::new())
      c = c + 1
    }
    self.data.push(new_row)
    let count_row : Array[Int] = Array::new()
    let mut c2 = 0
    while c2 < self.colnames.length() {
      count_row.push(0)
      c2 = c2 + 1
    }
    self.counts.push(count_row)
    row_idx = self.rownames.length() - 1
  }

  let mut col_idx = -1
  let mut j = 0
  while j < self.colnames.length() {
    if self.colnames[j] == record.sample_id {
      col_idx = j
      break
    }
    j = j + 1
  }

  if col_idx < 0 {
    self.colnames.push(record.sample_id)
    let mut r = 0
    while r < self.data.length() {
      self.data[r].push(Array::new())
      self.counts[r].push(0)
      r = r + 1
    }
    col_idx = self.colnames.length() - 1
  }

  self.data[row_idx][col_idx].push(record)
  self.counts[row_idx][col_idx] = self.counts[row_idx][col_idx] + 1
}

///|
/// Recalculate sample TMB and counts.
pub fn RaggedExperiment::recalculate(self : RaggedExperiment) -> Unit {
  // Recount from data
  let n_rows = self.rownames.length()
  let n_cols = self.colnames.length()
  let new_counts : Array[Array[Int]] = Array::new()
  let mut r = 0
  while r < n_rows {
    let row_counts : Array[Int] = Array::new()
    let mut c = 0
    while c < n_cols {
      row_counts.push(self.data[r][c].length())
      c = c + 1
    }
    new_counts.push(row_counts)
    r = r + 1
  }
  self.counts = new_counts

  // Calculate TMB per sample
  let tmb : Array[Double] = Array::new()
  let mut c2 = 0
  while c2 < n_cols {
    let mut total = 0
    let mut rr = 0
    while rr < n_rows {
      total = total + self.data[rr][c2].length()
      rr = rr + 1
    }
    tmb.push(total.to_double())
    c2 = c2 + 1
  }
  self.sample_tmb = tmb
}

///|
/// Get number of rows (genes).
pub fn RaggedExperiment::n_rows(self : RaggedExperiment) -> Int {
  self.rownames.length()
}

///|
/// Get number of columns (samples).
pub fn RaggedExperiment::n_cols(self : RaggedExperiment) -> Int {
  self.colnames.length()
}

///|
/// Get row names.
pub fn RaggedExperiment::get_rownames(self : RaggedExperiment) -> Array[String] {
  self.rownames
}

///|
/// Get column names.
pub fn RaggedExperiment::get_colnames(self : RaggedExperiment) -> Array[String] {
  self.colnames
}

///|
/// Get mutation records for a specific gene and sample.
pub fn RaggedExperiment::get_records(
  self : RaggedExperiment,
  gene~ : String,
  sample~ : String,
) -> Array[MutationRecord] {
  let ri = row_index(self, gene)
  let ci = col_index(self, sample)
  if ri < 0 || ci < 0 { return Array::new() }
  self.data[ri][ci]
}

///|
/// Get all records for a gene across all samples.
pub fn RaggedExperiment::get_gene_records(
  self : RaggedExperiment,
  gene~ : String,
) -> Array[MutationRecord] {
  let ri = row_index(self, gene)
  if ri < 0 { return Array::new() }
  let result : Array[MutationRecord] = Array::new()
  let mut c = 0
  while c < self.colnames.length() {
    let recs = self.data[ri][c]
    let mut j = 0
    while j < recs.length() {
      result.push(recs[j])
      j = j + 1
    }
    c = c + 1
  }
  result
}

///|
/// Get all records for a sample across all genes.
pub fn RaggedExperiment::get_sample_records(
  self : RaggedExperiment,
  sample~ : String,
) -> Array[MutationRecord] {
  let ci = col_index(self, sample)
  if ci < 0 { return Array::new() }
  let result : Array[MutationRecord] = Array::new()
  let mut r = 0
  while r < self.rownames.length() {
    let recs = self.data[r][ci]
    let mut j = 0
    while j < recs.length() {
      result.push(recs[j])
      j = j + 1
    }
    r = r + 1
  }
  result
}

///|
/// Get the mutation count matrix (genes x samples).
pub fn RaggedExperiment::get_count_matrix(self : RaggedExperiment) -> Array[Array[Int]] {
  self.counts
}

///|
/// Get TMB (Tumor Mutational Burden) per sample.
/// Returns total mutation count per sample.
pub fn RaggedExperiment::get_tmb(self : RaggedExperiment) -> Array[Double] {
  if self.sample_tmb.length() == 0 {
    self.recalculate()
  }
  self.sample_tmb
}

///|
/// Get TMB as mutations per Mb (requires genome size).
pub fn RaggedExperiment::get_tmb_per_mb(
  self : RaggedExperiment,
  genome_size_mb? : Double = 3000.0,
) -> Array[Double] {
  let tmb = self.get_tmb()
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < tmb.length() {
    result.push(tmb[i] / genome_size_mb)
    i = i + 1
  }
  result
}

///|
/// Summary statistics: number of genes mutated per sample.
pub fn RaggedExperiment::genes_mutated_per_sample(self : RaggedExperiment) -> Array[Int] {
  let n_cols = self.colnames.length()
  let n_rows = self.rownames.length()
  let result : Array[Int] = Array::new()
  let mut c = 0
  while c < n_cols {
    let mut count = 0
    let mut r = 0
    while r < n_rows {
      if self.data[r][c].length() > 0 { count = count + 1 }
      r = r + 1
    }
    result.push(count)
    c = c + 1
  }
  result
}

///|
/// Get summary statistics.
pub fn RaggedExperiment::summary(self : RaggedExperiment) -> String {
  let n_rows = self.n_rows()
  let n_cols = self.n_cols()
  let tmb = self.get_tmb()
  let genes_mutated = self.genes_mutated_per_sample()

  let mut s = "=== RaggedExperiment Summary ===" + "\n"
  s = s + "Genes (rows): " + n_rows.to_string() + "\n"
  s = s + "Samples (cols): " + n_cols.to_string() + "\n"
  s = s + "\nPer-sample TMB and mutated genes:" + "\n"
  let mut i = 0
  while i < n_cols {
    s = s + "  " + self.colnames[i] + ": TMB=" + tmb[i].to_string() + ", mutated_genes=" + genes_mutated[i].to_string() + "\n"
    i = i + 1
  }
  s
}

///|
/// Filter by mutation type.
pub fn RaggedExperiment::filter_by_type(
  self : RaggedExperiment,
  types~ : Array[MutationType],
) -> RaggedExperiment {
  let new_exp = RaggedExperiment::new()
  let mut r = 0
  while r < self.rownames.length() {
    let mut c = 0
    while c < self.colnames.length() {
      let recs = self.data[r][c]
      let mut j = 0
      while j < recs.length() {
        let mut keep = false
        let mut t = 0
        while t < types.length() {
          if recs[j].mutation_type == types[t] { keep = true; break }
          t = t + 1
        }
        if keep {
          new_exp.add_record(recs[j])
        }
        j = j + 1
      }
      c = c + 1
    }
    r = r + 1
  }
  new_exp.recalculate()
  new_exp
}

///|
/// Filter by gene list.
pub fn RaggedExperiment::filter_by_genes(
  self : RaggedExperiment,
  genes~ : Array[String],
) -> RaggedExperiment {
  let new_exp = RaggedExperiment::new()
  let mut r = 0
  while r < self.rownames.length() {
    let mut keep = false
    let mut g = 0
    while g < genes.length() {
      if self.rownames[r] == genes[g] { keep = true; break }
      g = g + 1
    }
    if keep {
      let mut c = 0
      while c < self.colnames.length() {
        let recs = self.data[r][c]
        let mut j = 0
        while j < recs.length() {
          new_exp.add_record(recs[j])
          j = j + 1
        }
        c = c + 1
      }
    }
    r = r + 1
  }
  new_exp.recalculate()
  new_exp
}

///|
/// Subset to specific samples.
pub fn RaggedExperiment::filter_by_samples(
  self : RaggedExperiment,
  samples~ : Array[String],
) -> RaggedExperiment {
  let new_exp = RaggedExperiment::new()
  let mut r = 0
  while r < self.rownames.length() {
    let mut c = 0
    while c < self.colnames.length() {
      let mut keep = false
      let mut s = 0
      while s < samples.length() {
        if self.colnames[c] == samples[s] { keep = true; break }
        s = s + 1
      }
      if keep {
        let recs = self.data[r][c]
        let mut j = 0
        while j < recs.length() {
          new_exp.add_record(recs[j])
          j = j + 1
        }
      }
      c = c + 1
    }
    r = r + 1
  }
  new_exp.recalculate()
  new_exp
}

///|
/// Create sample data for testing.
pub fn ragged_sample_data() -> RaggedExperiment {
  let exp = RaggedExperiment::new()

  // Sample 1: TP53 missense, KRAS mutation
  exp.add_record(MutationRecord::new(
    sample_id="Sample1", gene_symbol="TP53", chrom="chr17", pos=7577121, ref_allele="C", alt_allele="T",
    mutation_type=MutationType::Missense_Mutation, variant_allele_frequency=0.45,
  ))
  exp.add_record(MutationRecord::new(
    sample_id="Sample1", gene_symbol="KRAS", chrom="chr12", pos=25398284, ref_allele="G", alt_allele="T",
    mutation_type=MutationType::Missense_Mutation, variant_allele_frequency=0.32,
  ))
  exp.add_record(MutationRecord::new(
    sample_id="Sample1", gene_symbol="MYC", chrom="chr8", pos=128748315, ref_allele="C", alt_allele="A",
    mutation_type=MutationType::Nonsense_Mutation, variant_allele_frequency=0.28,
  ))

  // Sample 2: TP53 frameshift, EGFR mutation
  exp.add_record(MutationRecord::new(
    sample_id="Sample2", gene_symbol="TP53", chrom="chr17", pos=7578456, ref_allele="AG", alt_allele="A",
    mutation_type=MutationType::Frame_Shift_Del, variant_allele_frequency=0.52,
  ))
  exp.add_record(MutationRecord::new(
    sample_id="Sample2", gene_symbol="EGFR", chrom="chr7", pos=55086714, ref_allele="G", alt_allele="A",
    mutation_type=MutationType::Missense_Mutation, variant_allele_frequency=0.38,
  ))

  // Sample 3: KRAS mutation only
  exp.add_record(MutationRecord::new(
    sample_id="Sample3", gene_symbol="KRAS", chrom="chr12", pos=25398284, ref_allele="G", alt_allele="C",
    mutation_type=MutationType::Missense_Mutation, variant_allele_frequency=0.41,
  ))

  exp.recalculate()
  exp
}

// --- Helper functions ---

fn row_index(exp : RaggedExperiment, gene : String) -> Int {
  let mut i = 0
  while i < exp.rownames.length() {
    if exp.rownames[i] == gene { return i }
    i = i + 1
  }
  -1
}

fn col_index(exp : RaggedExperiment, sample : String) -> Int {
  let mut i = 0
  while i < exp.colnames.length() {
    if exp.colnames[i] == sample { return i }
    i = i + 1
  }
  -1
}
