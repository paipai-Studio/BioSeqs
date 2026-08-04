///|
fn format_coordinate_path(path : Array[Array[Int]]) -> String {
  let output = StringBuilder::new()
  output.write_string("[")
  for i in 0..<path.length() {
    if i > 0 {
      output.write_string(", ")
    }
    output.write_string("[")
    for j in 0..<path[i].length() {
      if j > 0 {
        output.write_string(", ")
      }
      output.write_string(path[i][j].to_string())
    }
    output.write_string("]")
  }
  output.write_string("]")
  output.to_string()
}

///|
fn main {
  println("=== Biopython Bio.Align.bigbed Demo ===")
  let targets = [
    @src.BigBedTarget::create("chr1", 1000000) catch {
      _ => abort("failed to create chr1")
    },
    @src.BigBedTarget::create("chr2", 500000) catch {
      _ => abort("failed to create chr2")
    },
  ]
  let records = [
    @src.BigBedRecord::create(
      "chr1",
      100,
      280,
      name="transcript-A",
      score=960,
      strand="+",
      thick_start=120,
      thick_end=260,
      item_rgb="40,120,220",
      block_count=3,
      block_sizes=[50, 40, 30],
      block_starts=[0, 80, 150],
      extra_fields=["GENE1"],
    ) catch {
      _ => abort("failed to create transcript-A")
    },
    @src.BigBedRecord::create(
      "chr1",
      250,
      410,
      name="transcript-B",
      score=870,
      strand="-",
      thick_start=270,
      thick_end=390,
      item_rgb="220,80,80",
      block_count=2,
      block_sizes=[60, 50],
      block_starts=[0, 110],
      extra_fields=["GENE1"],
    ) catch {
      _ => abort("failed to create transcript-B")
    },
    @src.BigBedRecord::create(
      "chr2",
      1000,
      1120,
      name="transcript-C",
      score=700,
      strand="+",
      block_count=2,
      block_sizes=[45, 50],
      block_starts=[0, 70],
      extra_fields=["GENE2"],
    ) catch {
      _ => abort("failed to create transcript-C")
    },
  ]
  let gene_field = @src.BigBedField::create(
    "string",
    "gene",
    comment="Gene identifier",
  ) catch {
    _ => abort("failed to create custom AutoSQL field")
  }
  let config = @src.BigBedWriteConfig::create(
    bed_columns=12,
    items_per_slot=2,
    block_size=2,
    compress=true,
    schema_name="transcripts",
    schema_comment="Indexed transcript annotations",
    custom_fields=[gene_field],
  ) catch {
    _ => abort("failed to create BigBed configuration")
  }

  println("\n1. Write and parse a compressed BigBed v4 stream")
  let bytes = @src.bigbed_write(targets, records, config~) catch {
    _ => abort("failed to write BigBed")
  }
  let file = @src.bigbed_parse(bytes) catch {
    _ => abort("failed to parse BigBed")
  }
  println("   bytes=" + bytes.length().to_string())
  let summary = file.summary() catch {
    _ => abort("failed to summarize BigBed")
  }
  println(
    "   records=\{summary.record_count}, targets=\{summary.target_count}, blocks=\{summary.data_block_count}, covered bases=\{summary.covered_bases}, compressed=\{summary.compressed}",
  )

  println("\n2. Inspect the embedded AutoSQL schema")
  println(file.schema.to_auto_sql())

  println("3. Query an overlapping chromosome interval through the R-tree")
  let hits = file.search("chr1", start=240, end=300) catch {
    _ => abort("failed to search BigBed")
  }
  for hit in hits {
    let gene = match
      hit.annotation(file.schema, "gene", file.header.defined_field_count) {
      Some(value) => value
      None => "missing"
    }
    println(
      "   " +
      hit.name +
      " " +
      hit.chrom +
      ":" +
      hit.chrom_start.to_string() +
      "-" +
      hit.chrom_end.to_string() +
      " gene=" +
      gene,
    )
  }

  println("\n4. Recover exon-aware reverse-strand coordinates")
  let reverse = file.find_by_name("transcript-B") catch {
    _ => abort("failed to query transcript name")
  }
  println(
    "   target/query path=" + format_coordinate_path(reverse[0].coordinates()),
  )

  println("\n5. Export records as BED text")
  println(file.to_bed() catch { _ => abort("failed to export BED") })

  println("\n6. Diagnose a damaged BigBed header")
  let damaged = bytes.copy()
  damaged[0] = 0
  try {
    ignore(@src.bigbed_parse(damaged))
    abort("damaged BigBed unexpectedly parsed")
  } catch {
    BigBedError(message) => println("   rejected: " + message)
  }

  println("\n=== Demo Complete ===")
}
