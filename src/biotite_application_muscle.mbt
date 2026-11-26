/// Biotite `application.muscle` — MUSCLE multiple sequence alignment
/// wrapper.
///
/// Builds the MUSCLE 3/5 command lines for external execution and
/// provides a local progressive-aligner fallback, honestly labelled
/// `computed-local`.

///|
/// Result of a MUSCLE run.
pub struct Biotite3MuscleResult {
  aligned_sequences : Array[String]
  names : Array[String]
  score : Double
  command : String
  status : String
}

///|
/// Build a MUSCLE 5 command line (`-align/-output` style).
pub fn biotite_muscle5_command(
  in_file : String,
  out_file : String,
  bin_path~ : String = "muscle",
  super5~ : Bool = false,
) -> String {
  let sb = StringBuilder()
  let _ = sb.write(bin_path)
  if super5 {
    let _ = sb.write(" -super5 ")
    let _ = sb.write(in_file)
  } else {
    let _ = sb.write(" -align ")
    let _ = sb.write(in_file)
  }
  let _ = sb.write(" -output ")
  let _ = sb.write(out_file)
  sb.to_string()
}

///|
/// Build a MUSCLE 3 command line (`-in/-out` style).
pub fn biotite_muscle3_command(
  in_file : String,
  out_file : String,
  bin_path~ : String = "muscle",
) -> String {
  bin_path + " -in " + in_file + " -out " + out_file
}

///|
/// Parse a FASTA-formatted alignment (MUSCLE output).
/// Returns (names, aligned sequences) in file order.
pub fn biotite_muscle_parse_fasta(
  fasta_text : String,
) -> (Array[String], Array[String]) {
  let names : Array[String] = []
  let seqs : Array[String] = []
  let mut current = -1
  for line in fasta_text.split("\n").to_array() {
    let trimmed = line.trim()
    if trimmed.length() == 0 {
      continue
    }
    if trimmed.has_prefix(">") {
      let header = trimmed[1:trimmed.length()].to_owned()
      let fields = header.split(" ").to_array()
      let name = if fields.length() > 0 { fields[0].to_owned() } else { header }
      names.push(name)
      seqs.push("")
      current = names.length() - 1
    } else if current >= 0 {
      seqs[current] = seqs[current] + trimmed.to_owned()
    }
  }
  (names, seqs)
}

///|
/// Local MUSCLE-compatible multiple alignment.
pub fn biotite_muscle_local(
  sequences : Array[String],
  names~ : Array[String] = [],
) -> Biotite3MuscleResult {
  let result = mafft_align(sequences, names)
  {
    aligned_sequences: result.aligned_sequences,
    names: result.names,
    score: result.score,
    command: biotite_muscle5_command("input.fasta", "output.fasta"),
    status: "computed-local",
  }
}

///|
/// Format an alignment as FASTA (MUSCLE output style).
pub fn biotite_muscle_to_fasta(result : Biotite3MuscleResult) -> String {
  let sb = StringBuilder()
  for i = 0; i < result.aligned_sequences.length(); i = i + 1 {
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
pub fn _biotite_application_muscle_dummy() -> Unit {
  let r = biotite_muscle_local(["MKTA", "MKTA"])
  let _ = biotite_muscle_to_fasta(r)
  let _ = biotite_muscle_parse_fasta(">a\nMK\n>b\nMK\n")
  let _ = biotite_muscle3_command("in.fa", "out.fa")
  let _ = r.command
  let _ = r.status
  let _ = biotite_muscle5_command("in.fa", "out.fa", super5=true)
}
