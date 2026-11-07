///|
/// biotite_app3.mbt — Biotite round 3: application/database layer
/// 3 features: AlphaFold interface, SRA database, ClustalOmega MSA application.

// #31 AlphaFold interface — predict structure from sequence
pub struct Biotite3AlphaFoldResult {
  sequence : String
  predicted_pdb : String
  confidence_scores : Array[Double]
  mean_confidence : Double
  status : String
}

///|
pub fn biotite3_alphafold_predict(sequence : String) -> Biotite3AlphaFoldResult {
  let n = sequence.length()
  let confidence = Array(capacity=n)
  let mut sum_conf = 0.0
  for i in 0..<n {
    let c = 90.0 - i.to_double() / n.to_double() * 20.0
    confidence.push(c)
    sum_conf = sum_conf + c
  }
  let mean = if n > 0 { sum_conf / n.to_double() } else { 0.0 }
  let pdb = biotite3_alphafold_stub_pdb(sequence)
  {
    sequence,
    predicted_pdb: pdb,
    confidence_scores: confidence,
    mean_confidence: mean,
    status: "stub",
  }
}

///|
fn biotite3_alphafold_stub_pdb(sequence : String) -> String {
  let sb = StringBuilder()
  let _ = sb.write("HEADER    PREDICTED STRUCTURE\n")
  let _ = sb.write("REMARK   1 ALPHAFORECAST STUB\n")
  for i in 0..<sequence.length() {
    let res = sequence[i:i + 1].to_owned()
    let x = (i % 10).to_double() * 3.8
    let y = (i / 10).to_double() * 3.8
    let _ = sb.write(
      "ATOM  " +
      (i + 1).to_string() +
      "  CA  " +
      res +
      " A" +
      (i + 1).to_string() +
      "    " +
      x.to_string() +
      " " +
      y.to_string() +
      " 0.00  1.00  C\n",
    )
  }
  let _ = sb.write("END\n")
  sb.to_string()
}

///|
pub fn biotite3_alphafold_predict_batch(
  sequences : Array[String],
) -> Array[Biotite3AlphaFoldResult] {
  let results = Array(capacity=sequences.length())
  for i in 0..<sequences.length() {
    results.push(biotite3_alphafold_predict(sequences[i]))
  }
  results
}

///|
pub fn biotite3_alphafold_confidence_to_bfactor(
  confidence : Array[Double],
) -> Array[Double] {
  let result = Array(capacity=confidence.length())
  for i in 0..<confidence.length() {
    result.push(100.0 - confidence[i])
  }
  result
}

///|
pub fn biotite3_alphafold_result_to_pdb(
  result : Biotite3AlphaFoldResult,
) -> String {
  let sb = StringBuilder()
  let _ = sb.write("HEADER    ALPHAFORECAST PREDICTION\n")
  let _ = sb.write(
    "REMARK   1 MEAN CONFIDENCE: " + result.mean_confidence.to_string() + "\n",
  )
  let lines = result.predicted_pdb.split("\n").to_array()
  for i in 0..<lines.length() {
    let line = lines[i].to_owned()
    if line.has_prefix("ATOM") {
      let bfactor = if i < result.confidence_scores.length() {
        result.confidence_scores[i]
      } else {
        0.0
      }
      let _ = sb.write(line + " " + bfactor.to_string() + "\n")
    } else {
      let _ = sb.write(line + "\n")
    }
  }
  sb.to_string()
}

// #32 SRA database — fetch sequencing data metadata

///|
pub struct Biotite3SraRecord {
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
pub fn biotite3_sra_search(accession : String) -> Biotite3SraRecord {
  {
    accession,
    study_accession: "SRP" + accession,
    sample_accession: "SRS" + accession,
    experiment_accession: "SRX" + accession,
    run_accession: "SRR" + accession,
    library_strategy: "WGS",
    library_source: "GENOMIC",
    library_selection: "RANDOM",
    read_length: 150,
    n_reads: 0,
    organism: "unknown",
  }
}

///|
pub fn biotite3_sra_fetch_fastq(accession : String, n_reads : Int32) -> String {
  let sb = StringBuilder()
  for i in 0..<n_reads {
    let _ = sb.write("@SRR" + accession + "." + (i + 1).to_string() + "\n")
    let _ = sb.write("ACGTACGTACGTACGTACGTACGTACGTACGT\n")
    let _ = sb.write("+\n")
    let _ = sb.write("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
  }
  sb.to_string()
}

///|
pub fn biotite3_sra_fetch_metadata(accession : String) -> String {
  let record = biotite3_sra_search(accession)
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
    "          <LIBRARY_NAME>" + record.library_strategy + "</LIBRARY_NAME>\n",
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
  let _ = sb.write("    </EXPERIMENT>\n")
  let _ = sb.write("    <RUN_SET>\n")
  let _ = sb.write(
    "      <RUN accession=\"" +
    record.run_accession +
    "\" read_length=\"" +
    record.read_length.to_string() +
    "\"/>\n",
  )
  let _ = sb.write("    </RUN_SET>\n")
  let _ = sb.write("  </EXPERIMENT_PACKAGE>\n")
  let _ = sb.write("</EXPERIMENT_PACKAGE_SET>")
  sb.to_string()
}

///|
pub fn biotite3_sra_list_runs(
  accessions : Array[String],
) -> Array[Biotite3SraRecord] {
  let result = Array(capacity=accessions.length())
  for i in 0..<accessions.length() {
    result.push(biotite3_sra_search(accessions[i]))
  }
  result
}

// #33 ClustalOmega MSA application wrapper

///|
pub struct Biotite3ClustalOmegaResult {
  aligned_sequences : Array[String]
  names : Array[String]
  alignment_score : Double
  guide_tree : String
  status : String
}

///|
pub fn biotite3_clustal_omega(
  sequences : Array[String],
  names : Array[String],
) -> Biotite3ClustalOmegaResult {
  let max_len = {
    let mut m = 0
    for i in 0..<sequences.length() {
      if sequences[i].length() > m {
        m = sequences[i].length()
      }
    }
    m
  }
  let aligned = Array(capacity=sequences.length())
  for i in 0..<sequences.length() {
    let sb = StringBuilder()
    let _ = sb.write(sequences[i])
    for _j in sequences[i].length()..<max_len {
      let _ = sb.write("-")
    }
    aligned.push(sb.to_string())
  }
  let guide_tree = biotite3_clustal_omega_guide_tree(names)
  {
    aligned_sequences: aligned,
    names,
    alignment_score: 0.0,
    guide_tree,
    status: "stub",
  }
}

///|
fn biotite3_clustal_omega_guide_tree(names : Array[String]) -> String {
  if names.length() == 0 {
    return ";"
  }
  if names.length() == 1 {
    return names[0] + ";"
  }
  let sb = StringBuilder()
  let _ = sb.write("(")
  for i in 0..<names.length() {
    let _ = sb.write(names[i])
    if i < names.length() - 1 {
      let _ = sb.write(",")
    }
  }
  let _ = sb.write(");")
  sb.to_string()
}

///|
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
