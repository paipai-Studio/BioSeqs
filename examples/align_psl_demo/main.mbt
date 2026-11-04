///|
fn parse_psl(text : String) -> @src.PslDocument {
  @src.psl_parse(text) catch {
    PslError(message) => abort("PSL parsing failed: " + message)
  }
}

///|
fn no_header(pslx : Bool) -> @src.PslWriteConfig {
  @src.PslWriteConfig::create(header=false, pslx~) catch {
    PslError(message) => abort("PSL writer configuration failed: " + message)
  }
}

///|
fn main {
  println("=== Biopython Bio.Align.psl Demo ===")

  println("\n1. Write and parse a multi-record PSL document")
  let alignments = @src.psl_example_data() catch {
    PslError(message) => abort("PSL example construction failed: " + message)
  }
  let text = @src.psl_write(alignments) catch {
    PslError(message) => abort("PSL writing failed: " + message)
  }
  let document = parse_psl(text)
  let summary = document.summary()
  println(
    "  records=" +
    summary.alignment_count.to_string() +
    ", nucleotide=" +
    summary.nucleotide_count.to_string() +
    ", translated=" +
    summary.protein_count.to_string(),
  )
  println(
    "  aligned query units=" +
    summary.aligned_query_units.to_string() +
    ", query insert bases=" +
    summary.query_insert_bases.to_string() +
    ", target insert bases=" +
    summary.target_insert_bases.to_string(),
  )

  println("\n2. Map a reverse-strand nucleotide query")
  let reverse = document.query("readReverse")[0]
  println("  strand reverse=" + reverse.is_reverse().to_string())
  match reverse.target_to_query(4) {
    Some(position) =>
      println("  target base 4 maps to query base " + position.to_string())
    None => println("  target base 4 is outside aligned blocks")
  }
  match reverse.query_to_target_interval(7) {
    Some((start, end)) =>
      println(
        "  query base 7 maps to target interval [" +
        start.to_string() +
        ", " +
        end.to_string() +
        ")",
      )
    None => println("  query base 7 is outside aligned blocks")
  }

  println("\n3. Emit and parse PSLX block sequences")
  let protein = alignments[2].recount() catch {
    PslError(message) => abort("translated recount failed: " + message)
  }
  let pslx = protein.format(config=no_header(true)) catch {
    PslError(message) => abort("PSLX formatting failed: " + message)
  }
  let pslx_alignment = parse_psl(pslx).alignments[0]
  println(
    "  PSLX blocks=" +
    pslx_alignment.blocks().length().to_string() +
    ", query fragment=" +
    pslx_alignment.query_block_sequences[0] +
    ", translated target fragment=" +
    pslx_alignment.target_block_sequences[0],
  )

  println("\n4. Recount a reverse-target translated alignment")
  let translated_reverse = @src.PslAlignment::create(
    "codingReverse",
    18,
    "peptide",
    4,
    [18, 12, 9, 3],
    [0, 2, 2, 4],
    sequence_type=@src.PslProtein,
    target_sequence="AAATCCAAAAAAAGCCAT",
    query_sequence="MAFG",
  ) catch {
    PslError(message) =>
      abort("reverse translated alignment failed: " + message)
  }
  let recounted = translated_reverse.recount() catch {
    PslError(message) => abort("reverse translated recount failed: " + message)
  }
  let translated_line = recounted.format(config=no_header(false)) catch {
    PslError(message) => abort("translated PSL formatting failed: " + message)
  }
  println(
    "  strand=+-, matches=" +
    recounted.matches.to_string() +
    ", mismatches=" +
    recounted.mismatches.to_string(),
  )
  println("  " + translated_line)

  println("=== Demo Complete ===")
}
