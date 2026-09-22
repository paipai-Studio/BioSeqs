///|
/// SeqKit `seq` / `fq2fa` / `fx2tab` / `tab2fx` ports — pure sequence
/// transforms and FASTA/FASTQ ↔ tabular conversion, reusing the IUPAC
/// complement table from `complement.mbt` and the FASTA/FASTQ writers.

///|
fn is_gap_char(c : Char) -> Bool {
  let ci = c.to_int()
  ci == 45 || ci == 32 || ci == 46 || ci == 9
}

///|
/// Reverse a sequence string (`seqkit seq -r`).
pub fn seqkit_reverse(seq : String) -> String {
  let buf = StringBuilder::new()
  let mut i = seq.length() - 1
  while i >= 0 {
    buf.write_char(seq.unsafe_get(i).unsafe_to_char())
    i = i - 1
  }
  buf.to_string()
}

///|
/// Complement a sequence (`seqkit seq -p`). IUPAC alphabet, case-preserving;
/// unknown characters pass through unchanged.
pub fn seqkit_complement(seq : String) -> String {
  complement_chars(seq, dna_complement_arr)
}

///|
/// Reverse-complement a sequence (`seqkit seq -r -p`).
pub fn seqkit_reverse_complement(seq : String) -> String {
  reverse_complement_chars(seq, dna_complement_arr)
}

///|
/// Remove gap letters (`seqkit seq -g`): `-`, space, tab and `.`.
pub fn seqkit_remove_gaps(seq : String) -> String {
  let buf = StringBuilder::new()
  for i = 0; i < seq.length(); i = i + 1 {
    let c = seq.unsafe_get(i).unsafe_to_char()
    if not(is_gap_char(c)) {
      buf.write_char(c)
    }
  }
  buf.to_string()
}

///|
/// Convert already-parsed FASTQ records to FASTA text (`seqkit fq2fa`): drops
/// quality. Equivalent to `write_fasta`.
pub fn seqkit_fq2fa(records : Array[SeqRecord]) -> String {
  write_fasta(records)
}

///|
fn phred_to_qual(quals : Array[Int]) -> String {
  let buf = StringBuilder::new()
  for q in quals {
    buf.write_char((q + 33).to_uint16().unsafe_to_char())
  }
  buf.to_string()
}

///|
/// Convert FASTA/FASTQ records to tabular text (`seqkit fx2tab`): one record
/// per line as `id \t seq` (FASTA) or `id \t seq \t qual` (FASTQ).
pub fn seqkit_fx2tab(records : Array[SeqRecord], is_fastq : Bool) -> String {
  let buf = StringBuilder::new()
  let mut first = true
  for rec in records {
    if not(first) {
      buf.write_char('\n')
    }
    buf.write_string(rec.id)
    buf.write_char('\t')
    buf.write_string(rec.seq.to_string())
    if is_fastq {
      buf.write_char('\t')
      match rec.letter_annotations.get("phred_quality") {
        None => ()
        Some(quals) => buf.write_string(phred_to_qual(quals))
      }
    }
    first = false
  }
  buf.to_string()
}

///|
/// Convert tabular text back to FASTA/FASTQ (`seqkit tab2fx`): a 2-column line
/// becomes a FASTA record, a 3-column line becomes a FASTQ record.
pub fn seqkit_tab2fx(content : String) -> String {
  let lines = content.split("\n")
  let buf = StringBuilder::new()
  for line in lines {
    if line.length() == 0 {
      continue
    }
    let cols = split_by_tab(line.to_owned())
    if cols.length() < 2 {
      continue
    }
    let id = cols[0]
    let seq = cols[1]
    if cols.length() >= 3 {
      buf.write_string("@" + id + "\n" + seq + "\n+\n" + cols[2] + "\n")
    } else {
      buf.write_string(">" + id + "\n" + seq + "\n")
    }
  }
  buf.to_string()
}///|
/// Build the `.fai` (FASTA index) text for a FASTA (`seqkit faidx`).
pub fn seqkit_faidx(fasta_content : String) -> String {
  write_fai(build_fai(fasta_content))
}

///|
/// Extract a subsequence by 1-based inclusive coordinates (`seqkit subseq
/// -r start:end`, positive coordinates only). Returns "" for an unknown
/// chromosome or an out-of-range region.
pub fn seqkit_subseq(fasta : Fasta, chrom : String, start : Int, end : Int) -> String {
  if start < 1 || end < start {
    return ""
  }
  match fasta.get_seq(chrom) {
    None => ""
    Some(seq) => {
      let full = seq.to_string()
      let s0 = start - 1
      if s0 >= 0 && end <= full.length() {
        full[s0:end].to_owned()
      } else {
        ""
      }
    }
  }
}