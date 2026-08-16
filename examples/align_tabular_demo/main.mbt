///|
fn parse_document(text : String) -> @src.AlignTabularDocument {
  @src.align_tabular_parse(text) catch {
    AlignTabularError(message) => abort("tabular parsing failed: " + message)
  }
}

///|
fn show_path(values : Array[Int]) -> String {
  let mut result = "["
  for index in 0..<values.length() {
    if index > 0 {
      result = result + ", "
    }
    result = result + values[index].to_string()
  }
  result + "]"
}

///|
fn main {
  println("=== Biopython Bio.Align.tabular Demo ===")

  println("\n1. Reconstruct BLAST outfmt 7 BTOP coordinate paths")
  let blast = "# BLASTP 2.14.1+\n" +
    "# Query: q1 protein query\n" +
    "# Database: proteins\n" +
    "# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, BTOP, query length, subject length, query seq, subject seq\n" +
    "# 2 hits found\n" +
    "q1\ts1\t75.0\t8\t2\t2\t3\t9\t10\t16\t1e-5\t42.0\t2AC1-GT-2\t20\t30\tAAAT-TGG\tAACTG-GG\n" +
    "q1\ts2\t75.0\t8\t2\t2\t3\t9\t4\t10\t0.02\t31.0\t2AC1-GT-2\t20\t25\tAAAT-TGG\tAACTG-GG\n" +
    "# BLAST processed 1 queries\n"
  let blast_document = parse_document(blast)
  let query = blast_document.queries[0]
  println(
    "  program=" +
    query.metadata.program +
    ", query=" +
    query.query_id +
    ", hits=" +
    query.alignments.length().to_string(),
  )
  for alignment in query.alignments {
    println("  " + alignment.summary())
    println(
      "    target path=" +
      show_path(alignment.target_coordinates) +
      ", query path=" +
      show_path(alignment.query_coordinates),
    )
  }
  let best = query.best_by_evalue() catch {
    AlignTabularError(message) => abort("best-hit selection failed: " + message)
  }
  match best {
    Some(alignment) => println("  best E-value target=" + alignment.target_id)
    None => println("  no E-value-bearing hit")
  }
  let coordinate = query.alignments[0].coordinate_alignment() catch {
    AlignTabularError(message) =>
      abort("coordinate conversion failed: " + message)
  }
  println(
    "  coordinate alignment columns=" +
    show_path(coordinate.target_coordinates) +
    " / " +
    show_path(coordinate.query_coordinates),
  )

  println("\n2. Decode FASTA -m 8CC aln_code")
  let fasta = "# fasta36 -q -m 8CC query.aa database.aa\n" +
    "# FASTA 36.3.8h May, 2020\n" +
    "# Query: qf protein query - 218 aa\n" +
    "# Database: database.aa\n" +
    "# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, aln_code\n" +
    "# 1 hits found\n" +
    "qf\tsf\t50.0\t9\t4\t2\t4\t10\t6\t13\t1e-4\t30.0\t3M2I1D3M\n" +
    "# FASTA processed 1 queries\n"
  let fasta_document = parse_document(fasta)
  let fasta_alignment = fasta_document.queries[0].alignments[0]
  println("  command=" + fasta_alignment.metadata.command_line)
  println(
    "  " +
    fasta_alignment.summary() +
    ", target gaps=" +
    fasta_alignment.target_gap_residues().to_string() +
    ", query gaps=" +
    fasta_alignment.query_gap_residues().to_string(),
  )

  println("\n3. Preserve reverse translated TBLASTX coordinates")
  let translated = "# TBLASTX 2.14.1+\n" +
    "# Query: qn translated query\n" +
    "# Database: nucleotides\n" +
    "# Fields: query id, subject id, query length, subject length, q. start, q. end, s. start, s. end, alignment length, BTOP\n" +
    "# 1 hits found\n" +
    "qn\tsn\t60\t90\t30\t19\t45\t34\t4\t4\n"
  let translated_alignment = parse_document(translated).queries[0].alignments[0]
  println(
    "  query units/residue=" +
    translated_alignment.query_units_per_residue.to_string() +
    ", target units/residue=" +
    translated_alignment.target_units_per_residue.to_string(),
  )
  println(
    "  reverse query=" +
    translated_alignment.query_is_reverse().to_string() +
    ", reverse target=" +
    translated_alignment.target_is_reverse().to_string(),
  )
  println(
    "  target path=" +
    show_path(translated_alignment.target_coordinates) +
    ", query path=" +
    show_path(translated_alignment.query_coordinates),
  )

  println("\n=== Demo Complete ===")
}
