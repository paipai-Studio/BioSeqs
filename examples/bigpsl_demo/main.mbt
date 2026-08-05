///|
fn main {
  println("=== Biopython Bio.Align.bigpsl Demo ===")
  let (targets, alignments) = @src.bigpsl_example_data() catch {
    BigPslError(message) => abort("failed to create BigPsl data: " + message)
  }
  let config = @src.BigPslWriteConfig::create(
    compress=true,
    block_size=3,
    items_per_slot=1,
  ) catch {
    BigPslError(message) => abort("failed to configure BigPsl: " + message)
  }

  println("\n1. Write and parse a compressed bed12+13 BigPsl stream")
  let bytes = @src.bigpsl_write(targets, alignments, config~) catch {
    BigPslError(message) => abort("failed to write BigPsl: " + message)
  }
  let file = @src.bigpsl_parse(bytes) catch {
    BigPslError(message) => abort("failed to parse BigPsl: " + message)
  }
  let summary = file.summary() catch {
    BigPslError(message) => abort("failed to summarize BigPsl: " + message)
  }
  println(
    "   bytes=\{bytes.length()}, alignments=\{summary.alignment_count}, targets=\{summary.target_count}, nucleotide=\{summary.nucleotide_count}, protein=\{summary.protein_count}, compressed=\{summary.compressed}",
  )

  println("\n2. Inspect the standard bigPsl AutoSQL schema")
  println(@src.bigpsl_schema().to_auto_sql())

  println("3. Query chr1 through the BigBed R-tree")
  let hits = file.search("chr1", start=105, end=106) catch {
    BigPslError(message) => abort("failed to query BigPsl: " + message)
  }
  println("   hits=" + hits.length().to_string())
  for alignment in hits {
    println("   " + alignment.summary())
    println(
      "   target chr1:105 -> query " +
      alignment.target_to_query(105).to_string(),
    )
  }

  println("\n4. Recover reverse-query nucleotide coordinates")
  let reverse = file.find_by_query("rna-reverse") catch {
    BigPslError(message) => abort("failed to find reverse query: " + message)
  }
  println("   " + reverse[0].summary())
  println(
    "   target chr1:200 -> query " + reverse[0].target_to_query(200).to_string(),
  )

  println("\n5. Map amino acids to reverse-target codon intervals")
  let protein = file.find_by_query("protein-reverse") catch {
    BigPslError(message) => abort("failed to find protein query: " + message)
  }
  println("   " + protein[0].summary())
  println(
    "   amino acid 0 -> target codon " +
    protein[0].query_to_target_interval(0).to_string(),
  )

  println("\n6. Export standard PSL text")
  println(
    file.to_psl() catch {
      BigPslError(message) => abort("failed to export PSL: " + message)
    },
  )
  println("=== Demo Complete ===")
}
