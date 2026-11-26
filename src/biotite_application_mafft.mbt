/// Biotite `application.mafft` — MAFFT multiple sequence alignment
/// wrapper.
///
/// Builds the MAFFT command line for external execution and provides a
/// local progressive-aligner fallback (the repository's MAFFT-lite
/// implementation), honestly labelled `computed-local`.

///|
/// Result of a MAFFT run.
pub struct Biotite3MafftResult {
  aligned_sequences : Array[String]
  names : Array[String]
  score : Double
  command : String
  status : String
}

///|
/// Build the MAFFT command line.
///
/// Mirrors Biotite's MafftApp defaults: quiet output, any-symbol
/// alphabet, optional global-pair (L-INS-i style) iteration.
pub fn biotite_mafft_command(
  sequence_file : String,
  bin_path~ : String = "mafft",
  quiet~ : Bool = true,
  any_symbol~ : Bool = true,
  global_pair~ : Bool = false,
  max_iterate~ : Int32 = 0,
) -> String {
  let sb = StringBuilder()
  let _ = sb.write(bin_path)
  if quiet {
    let _ = sb.write(" --quiet")
  }
  if any_symbol {
    let _ = sb.write(" --anysymbol")
  }
  if global_pair {
    let _ = sb.write(" --globalpair")
  }
  if global_pair || max_iterate > 0 {
    let _ = sb.write(" --maxiterate ")
    let _ = sb.write(
      if max_iterate > 0 {
        max_iterate.to_string()
      } else {
        "1000"
      },
    )
  }
  let _ = sb.write(" ")
  let _ = sb.write(sequence_file)
  sb.to_string()
}

///|
/// Parse a FASTA-formatted alignment (MAFFT output).
/// Returns (names, aligned sequences) in file order.
pub fn biotite_mafft_parse_fasta(
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
      // Header may carry a trailing description; keep the first token.
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
/// Local MAFFT-compatible multiple alignment.
///
/// Uses the repository's progressive aligner; the status string
/// `computed-local` makes clear the external binary was not called.
pub fn biotite_mafft_local(
  sequences : Array[String],
  names~ : Array[String] = [],
) -> Biotite3MafftResult {
  let result = mafft_align(sequences, names)
  {
    aligned_sequences: result.aligned_sequences,
    names: result.names,
    score: result.score,
    command: biotite_mafft_command("input.fasta"),
    status: "computed-local",
  }
}

///|
/// Format an alignment as FASTA (MAFFT output style).
pub fn biotite_mafft_to_fasta(result : Biotite3MafftResult) -> String {
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
pub fn _biotite_application_mafft_dummy() -> Unit {
  let r = biotite_mafft_local(["ACGU", "ACGC"])
  let _ = biotite_mafft_to_fasta(r)
  let _ = biotite_mafft_parse_fasta(">a\nAC\n>b\nAG\n")
  let _ = r.command
  let _ = r.status
  let _ = biotite_mafft_command("in.fa", global_pair=true, max_iterate=1000)
}
