///|
fn format_optional(value : Double?) -> String {
  match value {
    Some(v) => v.to_string()
    None => "None"
  }
}

///|
fn main {
  println("=== rtracklayer BigWig Demo ===")
  let targets = [
    @src.BigBedTarget::create("chr1", 1000) catch {
      _ => abort("failed to create chr1")
    },
    @src.BigBedTarget::create("chr2", 500) catch {
      _ => abort("failed to create chr2")
    },
  ]
  let intervals = [
    @src.BigWigInterval::create("chr1", 0, 10, 1.0),
    @src.BigWigInterval::create("chr1", 10, 20, 2.5),
    @src.BigWigInterval::create("chr1", 20, 30, -3.0),
    @src.BigWigInterval::create("chr2", 0, 100, 7.5),
  ]

  println("\n1. Write a compressed BigWig v4 stream (zlib blocks)")
  let config = @src.BigWigWriteConfig::create(block_size=64, compress=true) catch {
    _ => abort("failed to create BigWig configuration")
  }
  let bytes = @src.bigwig_write(targets, intervals, config=config) catch {
    _ => abort("failed to write BigWig")
  }
  let file = @src.bigwig_parse(bytes) catch {
    _ => abort("failed to parse BigWig")
  }
  println("   bytes=" + bytes.length().to_string())
  let header = file.header()
  println(
    "   version=\{header.version}, zoom levels=\{file.zoom_levels().length()}, sections=\{file.section_count()}, little endian=\{header.little_endian}",
  )
  let summary = file.total_summary()
  println(
    "   summary: bases=\{summary.valid_count}, min=\{summary.min_val}, max=\{summary.max_val}, sum=\{summary.sum_data}",
  )

  println("\n2. Query overlapping intervals (whole entries, unclipped)")
  let hits = file.intervals("chr1", start=5, end_=Some(25)) catch {
    _ => abort("failed to query intervals")
  }
  for hit in hits {
    println(
      "   " +
      hit.chrom +
      ":" +
      hit.start.to_string() +
      "-" +
      hit.end.to_string() +
      " value=" +
      hit.value.to_string(),
    )
  }

  println("\n3. Per-base values with NaN filling the gaps")
  let values = file.values("chr1", start=0, end_=Some(40)) catch {
    _ => abort("failed to query values")
  }
  let cells : Array[String] = []
  for i in 0..<values.length() {
    if i % 10 == 0 {
      cells.push(if values[i].is_nan() { "NaN" } else { values[i].to_string() })
    }
  }
  println("   value at every 10th base: " + cells.join(", "))

  println("\n4. Binned mean statistics (pyBigWig-compatible bins)")
  let means = file.stats("chr1", start=0, end_=Some(100), n_bins=4) catch {
    _ => abort("failed to compute means")
  }
  for i in 0..<means.length() {
    println("   bin " + i.to_string() + " mean=" + format_optional(means[i]))
  }

  println("\n5. Every statistic over the whole chr1 track")
  let statistics = ["mean", "min", "max", "cov", "sum", "std"]
  for name in statistics {
    let values_one = file.stats("chr1", n_bins=1, stat=name) catch {
      _ => abort("failed to compute " + name)
    }
    println("   " + name + " = " + format_optional(values_one[0]))
  }

  println("\n6. Writer rejects overlapping intervals")
  let bad = [
    @src.BigWigInterval::create("chr1", 0, 20, 1.0),
    @src.BigWigInterval::create("chr1", 10, 30, 2.0),
  ]
  try {
    ignore(@src.bigwig_write(targets, bad))
    abort("overlapping input unexpectedly accepted")
  } catch {
    @src.BigWigError(message) => println("   rejected: " + message)
  }

  println("\n7. Reader rejects a damaged BigWig header")
  let damaged = bytes.copy()
  damaged[0] = 0
  try {
    ignore(@src.bigwig_parse(damaged))
    abort("damaged BigWig unexpectedly parsed")
  } catch {
    @src.BigWigError(message) => println("   rejected: " + message)
  }

  println("\n=== Demo Complete ===")
}
