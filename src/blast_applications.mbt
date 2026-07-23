///|
/// Bio.Blast.Applications module - BLAST command-line tool wrappers.
/// Inspired by Biopython's Bio.Blast.Applications package.
/// Provides wrappers for the NCBI BLAST+ suite: blastp, blastn, blastx,
/// tblastn, tblastx, psiblast, rpsblast, makeblastdb.

// ===== Common parameter definitions =====

///|
/// A single command-line parameter definition.
pub struct BlastParamSpec {
  /// Parameter name (e.g., "-evalue")
  name : String
  /// Human-readable description
  description : String
  /// Whether the parameter takes a value (false for flags)
  takes_value : Bool
  /// Default value (empty if no default)
  default_value : String
  /// Whether parameter is required
  required : Bool
}

///|
/// Create a parameter specification.
pub fn BlastParamSpec::new(name : String, description : String, takes_value : Bool, default_value : String, required : Bool) -> BlastParamSpec {
  BlastParamSpec::{ name, description, takes_value, default_value, required }
}

///|
/// Common BLAST parameters shared by all variants.
fn blastapp_common_params() -> Array[BlastParamSpec] {
  [
    BlastParamSpec::new("-query", "Query sequence file", true, "", false),
    BlastParamSpec::new("-query_loc", "Query location (start-stop)", true, "", false),
    BlastParamSpec::new("-db", "BLAST database name", true, "", false),
    BlastParamSpec::new("-out", "Output file name", true, "", false),
    BlastParamSpec::new("-evalue", "Expectation value threshold", true, "10.0", false),
    BlastParamSpec::new("-word_size", "Word size for initial match", true, "0", false),
    BlastParamSpec::new("-gapopen", "Cost to open a gap", true, "0", false),
    BlastParamSpec::new("-gapextend", "Cost to extend a gap", true, "0", false),
    BlastParamSpec::new("-matrix", "Scoring matrix", true, "", false),
    BlastParamSpec::new("-threshold", "Minimum word score", true, "0", false),
    BlastParamSpec::new("-comp_based_stats", "Composition-based stats", true, "0", false),
    BlastParamSpec::new("-num_descriptions", "Number of descriptions", true, "500", false),
    BlastParamSpec::new("-num_alignments", "Number of alignments", true, "250", false),
    BlastParamSpec::new("-num_threads", "Number of CPU threads", true, "1", false),
    BlastParamSpec::new("-max_target_seqs", "Max target sequences", true, "500", false),
    BlastParamSpec::new("-dust", "DUST filter setting", true, "", false),
    BlastParamSpec::new("-seg", "SEG filter setting", true, "", false),
    BlastParamSpec::new("-soft_masking", "Soft masking", true, "false", false),
    BlastParamSpec::new("-lcase_masking", "Use lowercase masking", false, "", false),
    BlastParamSpec::new("-show_gis", "Show NCBI GIs in output", false, "", false),
    BlastParamSpec::new("-html", "Produce HTML output", false, "", false),
  ]
}

///|
/// Output format options parameters.
fn blastapp_output_params() -> Array[BlastParamSpec] {
  [
    BlastParamSpec::new("-outfmt", "Output format (0-18)", true, "0", false),
    BlastParamSpec::new("-max_hsps", "Max HSPs per subject", true, "0", false),
    BlastParamSpec::new("-max_intron_length", "Max intron length", true, "0", false),
  ]
}

// ===== Generic BLAST commandline =====

///|
/// Generic BLAST command-line wrapper.
pub struct BlastCommandline {
  /// Executable name (e.g., "blastp")
  executable : String
  /// Set parameters as (name, value) pairs
  parameters : Map[String, String]
  /// Flag parameters (true = set)
  flags : Map[String, Bool]
  /// Parameter specifications
  param_specs : Array[BlastParamSpec]
}

///|
/// Create a new BLAST commandline wrapper.
pub fn BlastCommandline::new(executable : String, param_specs : Array[BlastParamSpec]) -> BlastCommandline {
  BlastCommandline::{
    executable,
    parameters: Map([], capacity=32),
    flags: Map([], capacity=16),
    param_specs
  }
}

///|
/// Set a parameter value.
pub fn BlastCommandline::blastapp_set_parameter(self : BlastCommandline, name : String, value : String) -> BlastCommandline {
  let new_params = Map([], capacity=32)
  let keys = self.parameters.keys().collect()
  let mut i = 0
  while i < keys.length() {
    new_params.set(keys[i], self.parameters[keys[i]])
    i = i + 1
  }
  new_params.set(name, value)
  BlastCommandline::{
    executable: self.executable,
    parameters: new_params,
    flags: self.flags,
    param_specs: self.param_specs
  }
}

///|
/// Set a flag (boolean parameter).
pub fn BlastCommandline::blastapp_set_flag(self : BlastCommandline, name : String, on : Bool) -> BlastCommandline {
  let new_flags = Map([], capacity=16)
  let keys = self.flags.keys().collect()
  let mut i = 0
  while i < keys.length() {
    new_flags.set(keys[i], self.flags[keys[i]])
    i = i + 1
  }
  new_flags.set(name, on)
  BlastCommandline::{
    executable: self.executable,
    parameters: self.parameters,
    flags: new_flags,
    param_specs: self.param_specs
  }
}

///|
/// Get a parameter value (returns Option).
pub fn BlastCommandline::blastapp_get_parameter(self : BlastCommandline, name : String) -> String? {
  self.parameters.get(name)
}

///|
/// Get a flag value (returns Option).
pub fn BlastCommandline::blastapp_get_flag(self : BlastCommandline, name : String) -> Bool? {
  self.flags.get(name)
}

///|
/// Build the full command-line string.
pub fn BlastCommandline::blastapp_build_command(self : BlastCommandline) -> String {
  let mut cmd = self.executable
  // Add value parameters
  let keys = self.parameters.keys().collect()
  let mut i = 0
  while i < keys.length() {
    cmd = cmd + " " + keys[i] + " " + self.parameters[keys[i]]
    i = i + 1
  }
  // Add flags (only those that are true)
  let flag_keys = self.flags.keys().collect()
  let mut j = 0
  while j < flag_keys.length() {
    if self.flags[flag_keys[j]] {
      cmd = cmd + " " + flag_keys[j]
    }
    j = j + 1
  }
  cmd
}

///|
/// Validate that all required parameters are set.
pub fn BlastCommandline::blastapp_validate(self : BlastCommandline) -> Array[String] {
  let errors : Array[String] = Array::new()
  let mut i = 0
  while i < self.param_specs.length() {
    let spec = self.param_specs[i]
    if spec.required {
      let val = self.parameters.get(spec.name)
      if val is None {
        errors.push("Missing required parameter: " + spec.name)
      }
    }
    i = i + 1
  }
  errors
}

///|
/// List all available parameter specifications.
pub fn BlastCommandline::blastapp_available_params(self : BlastCommandline) -> Array[BlastParamSpec] {
  self.param_specs.copy()
}

// ===== Specific BLAST variant constructors =====

///|
/// NCBI blastp command line (protein query vs protein database).
pub fn ncbi_blastp_commandline() -> BlastCommandline {
  let specs = blastapp_common_params()
  let out_specs = blastapp_output_params()
  let mut i = 0
  while i < out_specs.length() {
    specs.push(out_specs[i])
    i = i + 1
  }
  // blastp-specific params
  specs.push(BlastParamSpec::new("-remote", "Execute search remotely", false, "", false))
  specs.push(BlastParamSpec::new("-subject", "Subject sequence file", true, "", false))
  specs.push(BlastParamSpec::new("-subject_loc", "Subject location", true, "", false))
  BlastCommandline::new("blastp", specs)
}

///|
/// NCBI blastn command line (nucleotide query vs nucleotide database).
pub fn ncbi_blastn_commandline() -> BlastCommandline {
  let specs = blastapp_common_params()
  let out_specs = blastapp_output_params()
  let mut i = 0
  while i < out_specs.length() {
    specs.push(out_specs[i])
    i = i + 1
  }
  // blastn-specific params
  specs.push(BlastParamSpec::new("-strand", "Query strand (both/minus/plus)", true, "both", false))
  specs.push(BlastParamSpec::new("-task", "Task (blastn/blastn-short/dc-megablast/etc.)", true, "megablast", false))
  specs.push(BlastParamSpec::new("-remote", "Execute search remotely", false, "", false))
  specs.push(BlastParamSpec::new("-subject", "Subject sequence file", true, "", false))
  specs.push(BlastParamSpec::new("-subject_loc", "Subject location", true, "", false))
  BlastCommandline::new("blastn", specs)
}

///|
/// NCBI blastx command line (translated nucleotide query vs protein database).
pub fn ncbi_blastx_commandline() -> BlastCommandline {
  let specs = blastapp_common_params()
  let out_specs = blastapp_output_params()
  let mut i = 0
  while i < out_specs.length() {
    specs.push(out_specs[i])
    i = i + 1
  }
  specs.push(BlastParamSpec::new("-strand", "Query strand (both/minus/plus)", true, "both", false))
  specs.push(BlastParamSpec::new("-query_gencode", "Query genetic code", true, "1", false))
  specs.push(BlastParamSpec::new("-max_intron_length", "Max intron length", true, "0", false))
  specs.push(BlastParamSpec::new("-remote", "Execute search remotely", false, "", false))
  specs.push(BlastParamSpec::new("-subject", "Subject sequence file", true, "", false))
  specs.push(BlastParamSpec::new("-subject_loc", "Subject location", true, "", false))
  BlastCommandline::new("blastx", specs)
}

///|
/// NCBI tblastn command line (protein query vs translated nucleotide database).
pub fn ncbi_tblastn_commandline() -> BlastCommandline {
  let specs = blastapp_common_params()
  let out_specs = blastapp_output_params()
  let mut i = 0
  while i < out_specs.length() {
    specs.push(out_specs[i])
    i = i + 1
  }
  specs.push(BlastParamSpec::new("-db_gencode", "Database genetic code", true, "1", false))
  specs.push(BlastParamSpec::new("-frame_shift_penalty", "Frame shift penalty", true, "0", false))
  specs.push(BlastParamSpec::new("-remote", "Execute search remotely", false, "", false))
  specs.push(BlastParamSpec::new("-subject", "Subject sequence file", true, "", false))
  specs.push(BlastParamSpec::new("-subject_loc", "Subject location", true, "", false))
  BlastCommandline::new("tblastn", specs)
}

///|
/// NCBI tblastx command line (translated nucleotide query vs translated nucleotide database).
pub fn ncbi_tblastx_commandline() -> BlastCommandline {
  let specs = blastapp_common_params()
  let out_specs = blastapp_output_params()
  let mut i = 0
  while i < out_specs.length() {
    specs.push(out_specs[i])
    i = i + 1
  }
  specs.push(BlastParamSpec::new("-strand", "Query strand (both/minus/plus)", true, "both", false))
  specs.push(BlastParamSpec::new("-query_gencode", "Query genetic code", true, "1", false))
  specs.push(BlastParamSpec::new("-db_gencode", "Database genetic code", true, "1", false))
  specs.push(BlastParamSpec::new("-remote", "Execute search remotely", false, "", false))
  specs.push(BlastParamSpec::new("-subject", "Subject sequence file", true, "", false))
  specs.push(BlastParamSpec::new("-subject_loc", "Subject location", true, "", false))
  BlastCommandline::new("tblastx", specs)
}

///|
/// NCBI psiblast command line (position-specific iterated BLAST).
pub fn ncbi_psiblast_commandline() -> BlastCommandline {
  let specs = blastapp_common_params()
  let out_specs = blastapp_output_params()
  let mut i = 0
  while i < out_specs.length() {
    specs.push(out_specs[i])
    i = i + 1
  }
  specs.push(BlastParamSpec::new("-num_iterations", "Number of iterations", true, "1", false))
  specs.push(BlastParamSpec::new("-in_pssm", "Input PSSM file", true, "", false))
  specs.push(BlastParamSpec::new("-out_pssm", "Output PSSM file", true, "", false))
  specs.push(BlastParamSpec::new("-out_ascii_pssm", "Output ASCII PSSM file", true, "", false))
  specs.push(BlastParamSpec::new("-save_pssm_after_last_round", "Save PSSM after last iteration", false, "", false))
  specs.push(BlastParamSpec::new("-save_each_pssm", "Save PSSM after each iteration", false, "", false))
  specs.push(BlastParamSpec::new("-pseudocount", "Pseudocount", true, "0", false))
  specs.push(BlastParamSpec::new("-inclusion_ethresh", "Inclusion e-value threshold", true, "0.002", false))
  specs.push(BlastParamSpec::new("-remote", "Execute search remotely", false, "", false))
  BlastCommandline::new("psiblast", specs)
}

///|
/// NCBI rpsblast command line (reverse position-specific BLAST).
pub fn ncbi_rpsblast_commandline() -> BlastCommandline {
  let specs = blastapp_common_params()
  specs.push(BlastParamSpec::new("-outfmt", "Output format (0-18)", true, "0", false))
  specs.push(BlastParamSpec::new("-max_hsps", "Max HSPs per subject", true, "0", false))
  specs.push(BlastParamSpec::new("-remote", "Execute search remotely", false, "", false))
  BlastCommandline::new("rpsblast", specs)
}

///|
/// NCBI makeblastdb command line (build BLAST database).
pub fn ncbi_makeblastdb_commandline() -> BlastCommandline {
  let specs : Array[BlastParamSpec] = [
    BlastParamSpec::new("-in", "Input FASTA file", true, "", false),
    BlastParamSpec::new("-input_type", "Input file type (fasta/asn1_bin/asn1_txt)", true, "fasta", false),
    BlastParamSpec::new("-dbtype", "Database type (nucl/prot)", true, "nucl", false),
    BlastParamSpec::new("-title", "Database title", true, "", false),
    BlastParamSpec::new("-parse_seqids", "Parse sequence IDs", false, "", false),
    BlastParamSpec::new("-hash_index", "Create hash index", false, "", false),
    BlastParamSpec::new("-mask_data", "Masking algorithm data file", true, "", false),
    BlastParamSpec::new("-mask_algo", "Masking algorithm ID", true, "", false),
    BlastParamSpec::new("-gilist", "GI list file", true, "", false),
    BlastParamSpec::new("-seqidlist", "Sequence ID list file", true, "", false),
    BlastParamSpec::new("-negative_gilist", "Negative GI list file", true, "", false),
    BlastParamSpec::new("-taxid", "Taxonomy ID", true, "", false),
    BlastParamSpec::new("-taxid_map", "Taxonomy ID map file", true, "", false),
    BlastParamSpec::new("-out", "Output database name", true, "", false),
    BlastParamSpec::new("-blastdb_version", "BLAST database version (4 or 5)", true, "5", false),
    BlastParamSpec::new("-max_file_sz", "Maximum file size", true, "1000000000", false),
  ]
  BlastCommandline::new("makeblastdb", specs)
}

// ===== Convenience builders =====

///|
/// Build a typical blastp query against a local database.
pub fn blastapp_quick_blastp(query_file : String, db_name : String, evalue : String) -> BlastCommandline {
  ncbi_blastp_commandline()
    .blastapp_set_parameter("-query", query_file)
    .blastapp_set_parameter("-db", db_name)
    .blastapp_set_parameter("-evalue", evalue)
}

///|
/// Build a typical blastn query against a local database.
pub fn blastapp_quick_blastn(query_file : String, db_name : String, evalue : String) -> BlastCommandline {
  ncbi_blastn_commandline()
    .blastapp_set_parameter("-query", query_file)
    .blastapp_set_parameter("-db", db_name)
    .blastapp_set_parameter("-evalue", evalue)
}

///|
/// Build a typical blastx query against a local protein database.
pub fn blastapp_quick_blastx(query_file : String, db_name : String, evalue : String) -> BlastCommandline {
  ncbi_blastx_commandline()
    .blastapp_set_parameter("-query", query_file)
    .blastapp_set_parameter("-db", db_name)
    .blastapp_set_parameter("-evalue", evalue)
}

///|
/// Build a typical tblastn query against a local nucleotide database.
pub fn blastapp_quick_tblastn(query_file : String, db_name : String, evalue : String) -> BlastCommandline {
  ncbi_tblastn_commandline()
    .blastapp_set_parameter("-query", query_file)
    .blastapp_set_parameter("-db", db_name)
    .blastapp_set_parameter("-evalue", evalue)
}

///|
/// Build a makeblastdb command for a given input file.
pub fn blastapp_quick_makeblastdb(input_file : String, dbtype : String, db_title : String) -> BlastCommandline {
  ncbi_makeblastdb_commandline()
    .blastapp_set_parameter("-in", input_file)
    .blastapp_set_parameter("-dbtype", dbtype)
    .blastapp_set_parameter("-title", db_title)
}

// ===== Example data generators =====

///|
/// Create a sample blastp commandline for demonstration.
pub fn blastapp_create_example_blastp() -> BlastCommandline {
  ncbi_blastp_commandline()
    .blastapp_set_parameter("-query", "query.fasta")
    .blastapp_set_parameter("-db", "nr")
    .blastapp_set_parameter("-evalue", "0.001")
    .blastapp_set_parameter("-out", "results.txt")
    .blastapp_set_parameter("-outfmt", "7")
    .blastapp_set_parameter("-num_threads", "4")
    .blastapp_set_parameter("-max_target_seqs", "10")
    .blastapp_set_flag("-show_gis", true)
}

///|
/// Create a sample blastn commandline for demonstration.
pub fn blastapp_create_example_blastn() -> BlastCommandline {
  ncbi_blastn_commandline()
    .blastapp_set_parameter("-query", "gene.fasta")
    .blastapp_set_parameter("-db", "nt")
    .blastapp_set_parameter("-evalue", "1e-50")
    .blastapp_set_parameter("-task", "megablast")
    .blastapp_set_parameter("-out", "blastn_results.xml")
    .blastapp_set_parameter("-outfmt", "5")
    .blastapp_set_parameter("-num_descriptions", "100")
    .blastapp_set_flag("-html", false)
}

///|
/// Create a sample makeblastdb commandline for demonstration.
pub fn blastapp_create_example_makeblastdb() -> BlastCommandline {
  ncbi_makeblastdb_commandline()
    .blastapp_set_parameter("-in", "genome.fa")
    .blastapp_set_parameter("-dbtype", "nucl")
    .blastapp_set_parameter("-title", "ExampleGenome")
    .blastapp_set_flag("-parse_seqids", true)
    .blastapp_set_flag("-hash_index", true)
}
