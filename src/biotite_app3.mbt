///|
/// biotite_app3.mbt — Biotite round 3: application/database layer.
///
/// Purely functional helpers only: URL/request construction and response
/// parsing. No network access and no fabricated prediction records.
/// `status` is either "parsed" (remote response decoded) or
/// "computed-local" (local algorithm, explicitly labelled).

// ============================================================
// AlphaFold / ESMFold interface
// ============================================================

///|
/// Result of decoding an AlphaFold/ESMFold PDB response.
pub(all) struct Biotite3AlphaFoldResult {
  sequence : String
  predicted_pdb : String
  confidence_scores : Array[Double]
  mean_confidence : Double
  status : String
}

///|
/// ESMFold endpoint used by the request builder.
pub fn biotite3_alphafold_request_url() -> String {
  "https://api.esmatlas.com/foldSequence/v1/pdb/"
}

///|
/// Request body for the ESMFold endpoint: the raw sequence.
pub fn biotite3_alphafold_request_body(sequence : String) -> String {
  sequence
}

///|
/// Decode a predicted PDB response.
///
/// pLDDT scores are read from the B-factor columns (60..66) of ATOM records,
/// which is the convention used by AlphaFold/ESMFold.
pub fn biotite3_alphafold_parse_pdb(
  pdb_text : String,
) -> Biotite3AlphaFoldResult {
  let lines = pdb_text.split("\n").to_array()
  let confidence : Array[Double] = []
  let residues : Array[String] = []
  for line in lines {
    if line.length() >= 66 && line.has_prefix("ATOM") {
      let b = match parse_double(line[60:66].trim().to_owned()) {
        Some(v) => v
        None => 0.0
      }
      confidence.push(b)
      if line[12:16].to_owned() == " CA " {
        let res = line[17:20].to_owned()
        if !residues.contains(res) || residues.length() == 0 {
          residues.push(res)
        }
      }
    }
  }
  let mut sum = 0.0
  for c in confidence {
    sum = sum + c
  }
  let mean = if confidence.length() > 0 {
    sum / confidence.length().to_double()
  } else {
    0.0
  }
  let seq_buf = StringBuilder()
  for res in residues {
    seq_buf.write_string(biotite3_three_to_one(res))
  }
  {
    sequence: seq_buf.to_string(),
    predicted_pdb: pdb_text,
    confidence_scores: confidence,
    mean_confidence: mean,
    status: "parsed",
  }
}

///|
/// pLDDT is already a 0..100 score and is written unchanged into the
/// B-factor column. This helper exists only to make the conversion explicit.
pub fn biotite3_alphafold_confidence_to_bfactor(
  confidence : Array[Double],
) -> Array[Double] {
  let result : Array[Double] = Array(capacity=confidence.length())
  for c in confidence {
    let clamped = if c < 0.0 { 0.0 } else if c > 100.0 { 100.0 } else { c }
    result.push(clamped)
  }
  result
}

///|
/// Normalize a decoded prediction into canonical fixed-column PDB records.
pub fn biotite3_alphafold_result_to_pdb(
  result : Biotite3AlphaFoldResult,
) -> String {
  AtomArray::from_pdb_text(result.predicted_pdb).to_pdb_string()
}

///|
/// Three-letter residue code to one-letter code for the 20 canonical residues.
fn biotite3_three_to_one(res : String) -> String {
  match res {
    "ALA" => "A"
    "ARG" => "R"
    "ASN" => "N"
    "ASP" => "D"
    "CYS" => "C"
    "GLN" => "Q"
    "GLU" => "E"
    "GLY" => "G"
    "HIS" => "H"
    "ILE" => "I"
    "LEU" => "L"
    "LYS" => "K"
    "MET" => "M"
    "PHE" => "F"
    "PRO" => "P"
    "SER" => "S"
    "THR" => "T"
    "TRP" => "W"
    "TYR" => "Y"
    "VAL" => "V"
    _ => "X"
  }
}

// ============================================================
// SRA / ENA interface
// ============================================================

///|
/// One sequencing run record decoded from an ENA filereport TSV.
pub(all) struct Biotite3SraRecord {
  accession : String
  study_accession : String
  sample_accession : String
  experiment_accession : String
  run_accession : String
  library_strategy : String
  library_source : String
  library_selection : String
  read_length : Int32
  n_reads : Int32
  organism : String
}

///|
/// ENA filereport URL for one run accession.
pub fn biotite3_sra_search_url(accession : String) -> String {
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" +
  accession +
  "&result=read_run&fields=run_accession,study_accession,sample_accession,experiment_accession,library_strategy,library_source,library_selection,read_count,read_length,scientific_name&format=tsv"
}

///|
/// ENA FASTQ download URL. `n_reads <= 0` means the whole run.
pub fn biotite3_sra_fastq_url(accession : String, n_reads : Int32) -> String {
  let base = "https://www.ebi.ac.uk/ena/browser/api/fastq/" + accession
  if n_reads > 0 {
    base + "?limit=" + n_reads.to_string()
  } else {
    base
  }
}

///|
/// Parse the header + rows of an ENA filereport TSV response.
pub fn biotite3_sra_parse_tsv(tsv : String) -> Array[Biotite3SraRecord] {
  let lines = tsv.split("\n").to_array()
  let result : Array[Biotite3SraRecord] = []
  if lines.length() < 2 {
    return result
  }
  let header = lines[0].split("\t").to_array()
  let index : Map[String, Int32] = Map([])
  for i in 0..<header.length() {
    index[header[i].trim().to_owned()] = i
  }
  for li in 1..<lines.length() {
    let line = lines[li]
    if line.trim().length() == 0 {
      continue
    }
    let fields = line.split("\t").to_array()
    let get = fn(key : String, default : String) -> String {
      match index.get(key) {
        Some(i) =>
          if i < fields.length() {
            fields[i].trim().to_owned()
          } else {
            default
          }
        None => default
      }
    }
    let get_int = fn(key : String) -> Int32 {
      match index.get(key) {
        Some(i) =>
          if i < fields.length() {
            match parse_double(fields[i].trim().to_owned()) {
              Some(v) => v.to_int()
              None => 0
            }
          } else {
            0
          }
        None => 0
      }
    }
    let organism = if get("scientific_name", "") != "" {
      get("scientific_name", "")
    } else {
      get("organism", "unknown")
    }
    result.push({
      accession: get("run_accession", ""),
      study_accession: get("study_accession", ""),
      sample_accession: get("sample_accession", ""),
      experiment_accession: get("experiment_accession", ""),
      run_accession: get("run_accession", ""),
      library_strategy: get("library_strategy", ""),
      library_source: get("library_source", ""),
      library_selection: get("library_selection", ""),
      read_length: get_int("read_length"),
      n_reads: get_int("read_count"),
      organism,
    })
  }
  result
}

///|
/// Serialize one decoded run record as an ENA-like EXPERIMENT_PACKAGE XML.
/// This is a serialization of real fields, not a fabricated response.
pub fn biotite3_sra_metadata_xml(record : Biotite3SraRecord) -> String {
  let sb = StringBuilder()
  let _ = sb.write("<?xml version=\"1.0\"?>\n")
  let _ = sb.write("<EXPERIMENT_PACKAGE_SET>\n")
  let _ = sb.write("  <EXPERIMENT_PACKAGE>\n")
  let _ = sb.write(
    "    <EXPERIMENT accession=\"" + record.experiment_accession + "\">\n",
  )
  let _ = sb.write(
    "      <STUDY_REF accession=\"" + record.study_accession + "\"/>\n",
  )
  let _ = sb.write("      <DESIGN>\n")
  let _ = sb.write("        <LIBRARY_DESCRIPTOR>\n")
  let _ = sb.write(
    "          <LIBRARY_STRATEGY>" +
    record.library_strategy +
    "</LIBRARY_STRATEGY>\n",
  )
  let _ = sb.write(
    "          <LIBRARY_SOURCE>" + record.library_source + "</LIBRARY_SOURCE>\n",
  )
  let _ = sb.write(
    "          <LIBRARY_SELECTION>" +
    record.library_selection +
    "</LIBRARY_SELECTION>\n",
  )
  let _ = sb.write("        </LIBRARY_DESCRIPTOR>\n")
  let _ = sb.write("      </DESIGN>\n")
  let _ = sb.write(
    "      <SAMPLE_REF accession=\"" + record.sample_accession + "\"/>\n",
  )
  let _ = sb.write("    </EXPERIMENT>\n")
  let _ = sb.write("    <RUN_SET>\n")
  let _ = sb.write(
    "      <RUN accession=\"" +
    record.run_accession +
    "\" read_length=\"" +
    record.read_length.to_string() +
    "\" read_count=\"" +
    record.n_reads.to_string() +
    "\"/>\n",
  )
  let _ = sb.write("    </RUN_SET>\n")
  let _ = sb.write("  </EXPERIMENT_PACKAGE>\n")
  let _ = sb.write("</EXPERIMENT_PACKAGE_SET>")
  sb.to_string()
}

// ============================================================
// Clustal Omega interface
// ============================================================

///|
/// Result of a Clustal Omega-compatible multiple alignment.
pub(all) struct Biotite3ClustalOmegaResult {
  aligned_sequences : Array[String]
  names : Array[String]
  alignment_score : Double
  guide_tree : String
  status : String
}

///|
/// Build the EBI Clustal Omega REST form parameters for a set of sequences.
pub fn biotite3_clustal_omega_request_params(
  sequences : Array[String],
  names : Array[String],
) -> String {
  let sb = StringBuilder()
  let _ = sb.write("email=anonymous%40example.org&stype=protein&outfmt=aln")
  let _ = sb.write("&sequence=")
  for i in 0..<sequences.length() {
    let name = if i < names.length() { names[i] } else { "seq" + i.to_string() }
    let entry = ">" + name + "\n" + sequences[i] + "\n"
    let _ = sb.write(biotite3_url_encode(entry))
  }
  sb.to_string()
}

///|
/// Percent-encode a string for application/x-www-form-urlencoded use.
fn biotite3_url_encode(text : String) -> String {
  let sb = StringBuilder()
  let hex = [
    "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F",
  ]
  for i in 0..<text.length() {
    let ch = text[i]
    let code = ch.to_int()
    let keep = (code >= 65 && code <= 90) ||
      (code >= 97 && code <= 122) ||
      (code >= 48 && code <= 57) ||
      ch == '-' ||
      ch == '_' ||
      ch == '.' ||
      ch == '~'
    if keep {
      let _ = sb.write(text[i:i + 1])
    } else if ch == ' ' {
      let _ = sb.write("+")
    } else {
      let hi = code / 16 % 16
      let lo = code % 16
      let _ = sb.write("%" + hex[hi] + hex[lo])
    }
  }
  sb.to_string()
}

///|
/// Parse a CLUSTAL-format alignment into a result.
///
/// Blocks are concatenated per sequence; consensus lines are ignored.
pub fn biotite3_clustal_omega_parse_aln(
  aln_text : String,
) -> Biotite3ClustalOmegaResult {
  let lines = aln_text.split("\n").to_array()
  let names : Array[String] = []
  let seqs : Array[String] = []
  for line in lines {
    let trimmed = line.trim()
    if trimmed.length() == 0 {
      continue
    }
    if trimmed.has_prefix("CLUSTAL") || trimmed.has_prefix("#") {
      continue
    }
    let fields = trimmed.split(" ").to_array()
    let parts : Array[String] = []
    for f in fields {
      if f.length() > 0 {
        parts.push(f.to_owned())
      }
    }
    if parts.length() < 2 {
      continue
    }
    let name = parts[0]
    let seq = parts[1]
    let mut only_consensus = true
    for i in 0..<seq.length() {
      let ch = seq[i]
      if ch != '*' && ch != ':' && ch != '.' && ch != ' ' {
        only_consensus = false
        break
      }
    }
    if only_consensus {
      continue
    }
    let mut found = -1
    for i in 0..<names.length() {
      if names[i] == name {
        found = i
        break
      }
    }
    if found >= 0 {
      seqs[found] = seqs[found] + seq
    } else {
      names.push(name)
      seqs.push(seq)
    }
  }
  let score = biotite3_clustal_omega_score(seqs)
  {
    aligned_sequences: seqs,
    names,
    alignment_score: score,
    guide_tree: "(" + biotite3_join_names(names) + ");",
    status: "parsed",
  }
}

///|
/// Score an alignment by counting identical columns.
fn biotite3_clustal_omega_score(seqs : Array[String]) -> Double {
  if seqs.length() == 0 {
    return 0.0
  }
  let len = seqs[0].length()
  let mut identical = 0
  for col in 0..<len {
    let ch = seqs[0][col:col + 1].to_owned()
    let mut same = true
    for i in 1..<seqs.length() {
      if col >= seqs[i].length() || seqs[i][col:col + 1].to_owned() != ch {
        same = false
        break
      }
    }
    if same {
      identical = identical + 1
    }
  }
  if len > 0 {
    identical.to_double() / len.to_double()
  } else {
    0.0
  }
}

///|
fn biotite3_join_names(names : Array[String]) -> String {
  let sb = StringBuilder()
  for i in 0..<names.length() {
    if i > 0 {
      let _ = sb.write(",")
    }
    let _ = sb.write(names[i])
  }
  sb.to_string()
}

///|
/// Local Clustal Omega-compatible alignment.
///
/// Reuses the repository's MAFFT-lite progressive aligner; the result is
/// explicitly labelled `computed-local` rather than pretending to be the
/// external Clustal Omega binary.
pub fn biotite3_clustal_omega_local(
  sequences : Array[String],
  names : Array[String],
) -> Biotite3ClustalOmegaResult {
  let result = mafft_align(sequences, names)
  {
    aligned_sequences: result.aligned_sequences,
    names: result.names,
    alignment_score: result.score,
    guide_tree: "(" + biotite3_join_names(result.names) + ");",
    status: "computed-local",
  }
}

///|
/// Backwards-compatible entry point: local alignment with an honest status.
pub fn biotite3_clustal_omega(
  sequences : Array[String],
  names : Array[String],
) -> Biotite3ClustalOmegaResult {
  biotite3_clustal_omega_local(sequences, names)
}

///|
/// Convert an alignment result to FASTA.
pub fn biotite3_clustal_omega_to_fasta(
  result : Biotite3ClustalOmegaResult,
) -> String {
  let sb = StringBuilder()
  for i in 0..<result.aligned_sequences.length() {
    let name = if i < result.names.length() {
      result.names[i]
    } else {
      "seq" + i.to_string()
    }
    let _ = sb.write(">" + name + "\n")
    let _ = sb.write(result.aligned_sequences[i] + "\n")
  }
  sb.to_string()
}

///|
/// Convert an alignment result to CLUSTAL block format.
pub fn biotite3_clustal_omega_to_clustal(
  result : Biotite3ClustalOmegaResult,
) -> String {
  let sb = StringBuilder()
  let _ = sb.write("CLUSTAL OMEGA ALIGNMENT\n\n")
  let seq_len = if result.aligned_sequences.length() > 0 {
    result.aligned_sequences[0].length()
  } else {
    0
  }
  let block_size = 60
  let mut pos = 0
  while pos < seq_len {
    let end = if pos + block_size < seq_len {
      pos + block_size
    } else {
      seq_len
    }
    for i in 0..<result.aligned_sequences.length() {
      let name = if i < result.names.length() {
        result.names[i]
      } else {
        "seq" + i.to_string()
      }
      let padded_name = name + "                "
      let name_field = padded_name[0:16].to_owned()
      let seq = result.aligned_sequences[i][pos:end].to_owned()
      let _ = sb.write(name_field + seq + " " + end.to_string() + "\n")
    }
    let _ = sb.write("\n")
    pos = pos + block_size
  }
  sb.to_string()
}

///|
pub fn _biotite_app3_dummy() -> Unit {
  let _ = Biotite3AlphaFoldResult::{
    sequence: "",
    predicted_pdb: "",
    confidence_scores: [],
    mean_confidence: 0.0,
    status: "",
  }
  let _ = Biotite3SraRecord::{
    accession: "",
    study_accession: "",
    sample_accession: "",
    experiment_accession: "",
    run_accession: "",
    library_strategy: "",
    library_source: "",
    library_selection: "",
    read_length: 0,
    n_reads: 0,
    organism: "",
  }
  let _ = Biotite3ClustalOmegaResult::{
    aligned_sequences: [],
    names: [],
    alignment_score: 0.0,
    guide_tree: "",
    status: "",
  }
}
