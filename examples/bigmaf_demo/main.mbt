///|
fn main {
  println("=== Biopython Bio.Align.bigmaf Demo ===")
  let (reference, targets, blocks) = @src.bigmaf_example_data() catch {
    BigMafError(message) => abort("failed to create BigMaf data: " + message)
  }
  let config = @src.BigMafWriteConfig::create(
    compress=true,
    block_size=3,
    items_per_slot=2,
  ) catch {
    BigMafError(message) => abort("failed to configure BigMaf: " + message)
  }

  println("\n1. Write and parse a compressed bed3+1 BigMaf stream")
  let bytes = @src.bigmaf_write(reference, targets, blocks, config~) catch {
    BigMafError(message) => abort("failed to write BigMaf: " + message)
  }
  let file = @src.bigmaf_parse(bytes) catch {
    BigMafError(message) => abort("failed to parse BigMaf: " + message)
  }
  let summary = file.summary() catch {
    BigMafError(message) => abort("failed to summarize BigMaf: " + message)
  }
  println(
    "   bytes=\{bytes.length()}, alignments=\{summary.alignment_count}, targets=\{summary.target_count}, components=\{summary.component_count}, compressed=\{summary.compressed}",
  )

  println("\n2. Inspect the standard bedMaf AutoSQL schema")
  println(@src.bigmaf_schema().to_auto_sql())

  println("3. Query chr7 through the BigBed R-tree")
  let hits = file.search("hg38.chr7", start=95, end=110) catch {
    BigMafError(message) => abort("failed to query BigMaf: " + message)
  }
  let block = hits[0]
  let reference_component = block.reference_component()
  println(
    "   hits=\{hits.length()}, reference=\{reference_component.source}:\{reference_component.start}-\{reference_component.start + reference_component.size}, columns=\{block.aligned_columns()}",
  )
  println(
    "   score=\{block.score}, pass=\{block.pass_number}, identity=\{block.pairwise_identity()}",
  )

  println("\n4. Map a reference base to a reverse-strand component")
  match block.map_reference_to(104, "mm39.chr5") {
    Some(position) =>
      println("   hg38.chr7:104 -> mm39.chr5:" + position.to_string())
    None => println("   reference base maps to a gap")
  }
  match block.component("mm39.chr5") {
    Some(component) => {
      let (start, end) = component.forward_interval()
      println(
        "   mouse forward interval=\{start}-\{end}, sequence=\{component.forward_sequence()}",
      )
    }
    None => abort("example mouse component is missing")
  }

  println("\n5. Preserve MAF annotations and export standard MAF")
  println(
    "   empty components=\{block.empty_components.length()}, comments=\{block.comments.length()}",
  )
  println(
    file.to_maf() catch {
      BigMafError(message) => abort("failed to export MAF: " + message)
    },
  )

  println("=== Demo Complete ===")
}
