///|
test "fqqc basic stats" {
  let fastq_ok = "@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nGGGGCCCC\n+\nFFFFFFFF\n"
  let st = @src.fqqc_compute(fastq_ok)
  assert_eq(st.n_reads, 2)
  assert_eq(st.total_bases, 16)
  assert_eq(st.per_base_hist.length(), 8)
}

///|
test "fqqc invalid fastq raises" {
  let fastq = "@r1\nACGTACGT\n+\nIIIIII\n"
  let msg = try {
    ignore(@src.fqqc_compute(fastq))
    "NO_ERROR"
  } catch {
    SeqIOError(m) => m
  }
  assert_true(msg.contains("don't match"))
}

///|
test "fqqc per base histogram counts" {
  let fastq = "@r1\nAA\n+\nII\n@r2\nAA\n+\nHH\n"
  let st = @src.fqqc_compute(fastq)
  // 位置0: I(40)+H(39); 位置1: I(40)+H(39)
  assert_eq(st.per_base_count[0], 2)
  assert_eq(st.per_base_hist[0][40], 1)
  assert_eq(st.per_base_hist[0][39], 1)
  assert_eq(st.per_base_hist[1][40], 1)
}

///|
test "fqqc gc percent" {
  let fastq = "@r1\nGGGGCCCC\n+\nIIIIIIII\n"
  let st = @src.fqqc_compute(fastq)
  assert_eq(st.gc_percents.length(), 1)
  assert_eq(st.gc_percents[0], 100)
  let fastq2 = "@r1\nAAAATTTT\n+\nIIIIIIII\n"
  let st2 = @src.fqqc_compute(fastq2)
  assert_eq(st2.gc_percents[0], 0)
}

///|
test "fqqc adapter detection" {
  let fastq = "@r1\nAGATCGGAAGAGCAAA\n+\nIIIIIIIIIIIIIIII\n@r2\nACGTACGT\n+\nIIIIIIII\n"
  let st = @src.fqqc_compute(fastq)
  assert_eq(st.adapter_hits.get("AGATCGGAAGAGC"), Some(1))
  assert_true(st.adapter_hits.get("GGGGGGGGGGGG") is None)
}

///|
test "fqqc overrepresented" {
  let fastq = "@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nACGTACGT\n+\nIIIIIIII\n@r3\nTTTTAAAA\n+\nIIIIIIII\n"
  let st = @src.fqqc_compute(fastq)
  assert_eq(st.seq_counts.get("ACGTACGT"), Some(2))
  assert_eq(st.seq_counts.get("TTTTAAAA"), Some(1))
}

///|
test "fqqc report json parseable" {
  let fastq_ok = "@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nGGGGCCCCGG\n+\nFFFFFFFFFF\n"
  let report = @src.fqqc_report(fastq_ok)
  let parsed = @src.parse_json(report)
  assert_true(parsed is Some(_))
  match parsed {
    Some(JsonObj(m)) => {
      assert_true(m.contains("basic"))
      assert_true(m.contains("per_base_quality"))
      assert_true(m.contains("per_sequence_quality"))
      assert_true(m.contains("per_sequence_gc"))
      assert_true(m.contains("per_base_gc"))
      assert_true(m.contains("length_distribution"))
      assert_true(m.contains("adapters"))
      assert_true(m.contains("overrepresented"))
    }
    _ => assert_true(false)
  }
}

///|
test "fqqc quantile from histogram" {
  // 4 条单碱基 read, 位置 0 的质量值 [1,3,5,7]: 中位 4
  let fastq = "@r1\nA\n+\n\"\n@r2\nA\n+\n$\n@r3\nA\n+\n&\n@r4\nA\n+\n(\n"
  let st = @src.fqqc_compute(fastq)
  let v = @src.fqqc_report_value(st)
  match v {
    JsonObj(m) =>
      match m.get("per_base_quality") {
        Some(JsonObj(pb)) =>
          match pb.get("median") {
            Some(JsonArr(med)) =>
              match med[0] {
                JsonNum(v) => assert_true((v - 4.0).abs() < 1.0e-9)
                _ => assert_true(false)
              }
            _ => assert_true(false)
          }
        _ => assert_true(false)
      }
    _ => assert_true(false)
  }
}

///|
test "fqqc empty input" {
  let st = @src.fqqc_compute("")
  assert_eq(st.n_reads, 0)
  let report = @src.fqqc_report("")
  assert_true(@src.parse_json(report) is Some(_))
}

///|
test "fqqc length distribution" {
  let fastq = "@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nGGGGCC\n+\nFFFFFF\n"
  let st = @src.fqqc_compute(fastq)
  assert_eq(st.lengths.get(8), Some(1))
  assert_eq(st.lengths.get(6), Some(1))
}
