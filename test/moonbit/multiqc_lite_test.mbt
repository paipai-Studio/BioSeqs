// ===========================================================================
// multiqc_lite_test.mbt — 多样本聚合报告测试
// ===========================================================================

///|
fn multiqc_test_fastq(sample : String) -> String {
  let r1 = "@\{sample}1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"
  let r2 = "@\{sample}2\nGGCCTTAAGGCC\n+\nJJJJJJJJJJJJ\n"
  r1 + r2
}

///|
test "multiqc module add sample" {
  let m = @src.multiqc_module("demo")
  let _ = m.add_sample("s1", [("reads", 100.0), ("gc", 45.5)])
  let _ = m.add_sample("s2", [("reads", 200.0), ("gc", 51.0)])
  assert_eq(m.name, "demo")
  assert_eq(m.samples.length(), 2)
  assert_eq(m.samples.get("s1").unwrap().get("reads").unwrap(), 100.0)
  assert_eq(m.samples.get("s2").unwrap().get("gc").unwrap(), 51.0)
}

///|
test "multiqc collect fqqc" {
  let m = @src.multiqc_collect_fqqc(multiqc_test_fastq("a"), "sampleA")
  assert_eq(m.name, "fqqc")
  let metrics = m.samples.get("sampleA").unwrap()
  assert_eq(metrics.get("reads").unwrap(), 2.0)
  assert_eq(metrics.get("mean_length").unwrap(), 12.0)
  // GC: read1 6/12=50, read2 8/12=66(整数截断) -> 平均 58
  let gc = metrics.get("gc_percent").unwrap()
  assert_true(gc > 57.5 && gc < 58.5)
  assert_true(metrics.get("mean_quality").unwrap() > 0.0)
}

///|
test "multiqc collect kallisto" {
  let tx_names : Array[String] = ["t1", "t2"]
  let tx_seqs : Array[String] = [
    "AAAAACCCCCGGGGGTTTTTATATGCGCGCTATAT", "TTTTTGGGGGCCCCCAAAAAATATGCGCGCTATAT",
  ]
  let quant = @src.kallisto_run(
    tx_names,
    tx_seqs,
    ["AAAAACCCCCGGGGGTTTTTATAT", "TTTTTGGGGGCCCCCAAAAA"],
    k=9,
  )
  let m = @src.multiqc_collect_kallisto(quant, "sampleK")
  assert_eq(m.name, "kallisto-lite")
  let metrics = m.samples.get("sampleK").unwrap()
  assert_eq(metrics.get("reads").unwrap(), 2.0)
  assert_eq(metrics.get("pseudoaligned_percent").unwrap(), 100.0)
  assert_eq(metrics.get("expressed_transcripts").unwrap(), 2.0)
}

///|
test "multiqc collect gatk" {
  let reference = "ATCGGATTCAGCCTAGGATCAAGTTCGGATACCTAGGCATTGCAAGGTCCTAGAAGT"
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < 10; i = i + 1 {
    let sliced = reference[10:30].to_owned()
    let mutated = if i < 5 {
      sliced[0:5].to_owned() + "C" + sliced[6:].to_owned()
    } else {
      sliced
    }
    reads.push({ name: "r" + i.to_string(), seq: mutated, pos: 10, })
  }
  let variants = @src.gatk_call_variants("chr1", reference, reads)
  let m = @src.multiqc_collect_gatk(variants, "sampleG")
  assert_eq(m.name, "gatk-lite")
  let metrics = m.samples.get("sampleG").unwrap()
  assert_eq(metrics.get("snps").unwrap(), 1.0)
  assert_eq(metrics.get("variants").unwrap(), 1.0)
}

///|
test "multiqc json roundtrip" {
  let report = @src.multiqc_report(title="Test Report")
  let _ = report.add_module(
    @src.multiqc_module("m1")
    .add_sample("s1", [("a", 1.0), ("b", 2.5)])
    .add_sample("s2", [("a", 3.0), ("b", 4.0)])
    .with_lower_better(["b"]),
  )
  let _ = report.add_module(
    @src.multiqc_module("m2").add_sample("s1", [("x", 9.0)]),
  )
  let json = @src.multiqc_write_json(report)
  let parsed = @src.multiqc_parse_json(json).unwrap()
  assert_eq(parsed.title, "Test Report")
  assert_eq(parsed.modules.length(), 2)
  assert_eq(parsed.modules[0].name, "m1")
  assert_eq(parsed.modules[0].samples.length(), 2)
  assert_eq(parsed.modules[0].samples.get("s1").unwrap().get("a").unwrap(), 1.0)
  assert_eq(parsed.modules[0].samples.get("s2").unwrap().get("b").unwrap(), 4.0)
  assert_eq(parsed.modules[0].lower_better, ["b"])
  assert_eq(parsed.modules[1].samples.get("s1").unwrap().get("x").unwrap(), 9.0)
}

///|
test "multiqc parse invalid json returns none" {
  let parsed = @src.multiqc_parse_json("not json {")
  assert_true(parsed is None)
}

///|
test "multiqc html output structure" {
  let report = @src.multiqc_report(title="HTML Test")
  let _ = report.add_module(
    @src.multiqc_module("mod1").add_sample("sa", [
      ("reads", 100.0),
      ("gc", 40.0),
    ]),
  )
  let html = @src.multiqc_write_html(report)
  assert_true(html.has_prefix("<!DOCTYPE html>"))
  assert_true(html.contains("<title>HTML Test</title>"))
  assert_true(html.contains("<h2>mod1</h2>"))
  assert_true(html.contains("<th>reads</th>"))
  assert_true(html.contains("<td>100</td>"))
  assert_true(html.contains("<td>40</td>"))
  assert_true(html.contains("multiqc-lite"))
}

///|
test "multiqc html coloring best worst" {
  let report = @src.multiqc_report()
  let _ = report.add_module(
    @src.multiqc_module("m")
    .add_sample("good", [("yield", 99.0)])
    .add_sample("bad", [("yield", 10.0)]),
  )
  let html = @src.multiqc_write_html(report)
  // 两个样本: 最优绿最差红(值不同时)
  assert_true(html.contains("class=\"good\""))
  assert_true(html.contains("class=\"bad\""))
}

///|
test "multiqc html single sample no coloring" {
  let report = @src.multiqc_report()
  let _ = report.add_module(
    @src.multiqc_module("m").add_sample("only", [("x", 5.0)]),
  )
  let html = @src.multiqc_write_html(report)
  assert_false(html.contains("class=\"good\""))
  assert_false(html.contains("class=\"bad\""))
}

///|
test "multiqc markdown output" {
  let report = @src.multiqc_report(title="MD Test")
  let _ = report.add_module(
    @src.multiqc_module("m1")
    .add_sample("s1", [("reads", 10.0)])
    .add_sample("s2", [("reads", 20.0)]),
  )
  let md = @src.multiqc_write_markdown(report)
  assert_true(md.has_prefix("# MD Test"))
  assert_true(md.contains("## m1"))
  assert_true(md.contains("| sample | reads |"))
  assert_true(md.contains("| s1 | 10 |"))
  assert_true(md.contains("| s2 | 20 |"))
}

///|
test "multiqc end to end aggregation" {
  // 两个样本 x 两个模块
  let report = @src.multiqc_report(title="E2E")
  let fqqc_a = @src.multiqc_collect_fqqc(multiqc_test_fastq("a"), "sampleA")
  let fqqc_b = @src.multiqc_collect_fqqc(multiqc_test_fastq("b"), "sampleB")
  let _ = report.add_module(fqqc_a)
  let _ = report.add_module(fqqc_b)
  let kal = @src.kallisto_run(
    ["t1"],
    ["ACGTACGTACGTACGTACGT"],
    ["ACGTACGTACGTACGT"],
    k=9,
  )
  let _ = report.add_module(@src.multiqc_collect_kallisto(kal, "sampleA"))
  assert_eq(report.modules.length(), 3)
  let s = report.summary()
  assert_true(s.contains("multiqc-lite"))
  assert_true(s.contains("3 modules"))
  assert_true(s.contains("2 samples"))
}

///|
test "multiqc summary empty report" {
  let report = @src.multiqc_report()
  let s = report.summary()
  assert_true(s.contains("0 modules"))
  assert_true(s.contains("0 samples"))
}

///|
test "multiqc merge modules same name" {
  let report = @src.multiqc_report()
  let _ = report.add_module(
    @src.multiqc_module("fqqc").add_sample("s1", [("reads", 10.0)]),
  )
  let _ = report.add_module(
    @src.multiqc_module("fqqc").add_sample("s2", [("reads", 20.0)]),
  )
  let _ = report.add_module(
    @src.multiqc_module("gatk").add_sample("s1", [("snps", 3.0)]),
  )
  let merged = report.merge_modules()
  assert_eq(merged.modules.length(), 2)
  assert_eq(merged.modules[0].name, "fqqc")
  assert_eq(merged.modules[0].samples.length(), 2)
  assert_eq(
    merged.modules[0].samples.get("s1").unwrap().get("reads").unwrap(),
    10.0,
  )
  assert_eq(
    merged.modules[0].samples.get("s2").unwrap().get("reads").unwrap(),
    20.0,
  )
  assert_eq(merged.modules[1].name, "gatk")
}
