// ===========================================================================
// read_simulator_test.mbt — 读段模拟器测试
// ===========================================================================

///|
fn read_sim_test_split_fastq(fastq : String) -> Array[Array[String]] {
  let lines = fastq.split("\n").to_array()
  let records : Array[Array[String]] = Array(capacity=0)
  let mut i = 0
  while i + 3 < lines.length() {
    records.push([
      lines[i].to_owned(),
      lines[i + 1].to_owned(),
      lines[i + 2].to_owned(),
      lines[i + 3].to_owned(),
    ])
    i = i + 4
  }
  records
}

///|
fn read_sim_test_reference() -> String {
  let sb = StringBuilder()
  for i = 0; i < 1000; i = i + 1 {
    sb.write_char(['A', 'C', 'G', 'T'][(i * 7 + 3) % 4])
  }
  sb.to_string()
}

///|
test "read simulator config validation" {
  let c = @src.read_simulator_config()
  assert_eq(c.read_length, 150)
  assert_true(c.paired_end)

  let mut raised = false
  try {
    let bad = @src.read_simulator_config(read_length=0)
    ignore(bad)
  } catch {
    _ => raised = true
  }
  assert_true(raised)

  raised = false
  try {
    let bad = @src.read_simulator_config(error_rate=1.5)
    ignore(bad)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "read simulator single end coverage and fastq shape" {
  let config = @src.read_simulator_config(
    read_length=100,
    coverage=25.0,
    paired_end=false,
    error_rate=0.0,
    seed=1UL,
  )
  let result = @src.simulate_reads(read_sim_test_reference(), config~)
  assert_eq(result.n_reads, 250)
  assert_eq(result.n_pairs, 0)
  assert_true(
    result.achieved_coverage > 24.9 && result.achieved_coverage < 25.1,
  )

  let records = read_sim_test_split_fastq(result.fastq)
  assert_eq(records.length(), 250)
  for rec in records {
    assert_true(rec[0].has_prefix("@read"))
    assert_eq(rec[1].length(), rec[3].length())
    assert_eq(rec[1].length(), 100)
    assert_eq(rec[2], "+")
  }
}

///|
test "read simulator paired end produces pairs" {
  let config = @src.read_simulator_config(
    read_length=100,
    coverage=20.0,
    paired_end=true,
    insert_mean=300,
    insert_sd=20,
    error_rate=0.0,
    seed=2UL,
  )
  let result = @src.simulate_reads(read_sim_test_reference(), config~)
  assert_eq(result.n_pairs, 100)
  assert_eq(result.n_reads, 200)
  assert_true(
    result.achieved_coverage > 19.9 && result.achieved_coverage < 20.1,
  )

  let records = read_sim_test_split_fastq(result.fastq)
  assert_eq(records.length(), 200)
  assert_true(records[0][0].has_suffix("/1"))
  assert_true(records[1][0].has_suffix("/2"))
}

///|
test "read simulator injects variants and errors" {
  let reference = read_sim_test_reference()
  let config = @src.read_simulator_config(
    read_length=100,
    coverage=10.0,
    paired_end=false,
    error_rate=0.01,
    snp_rate=0.01,
    indel_rate=0.005,
    seed=3UL,
  )
  let result = @src.simulate_reads(reference, config~)
  assert_true(result.n_variants > 0)
  assert_true(result.n_errors > 0)
}

///|
test "read simulator appends adapter" {
  let adapter = "AGATCGGAAGAG"
  let config = @src.read_simulator_config(
    read_length=100,
    coverage=5.0,
    paired_end=false,
    error_rate=0.0,
    adapter~,
    adapter_rate=1.0,
    seed=4UL,
  )
  let result = @src.simulate_reads(read_sim_test_reference(), config~)
  let records = read_sim_test_split_fastq(result.fastq)
  assert_eq(records.length(), 50)
  for rec in records {
    assert_true(rec[1].has_suffix(adapter))
    assert_eq(rec[1].length(), rec[3].length())
  }
}

///|
test "read simulator is deterministic" {
  let config = @src.read_simulator_config(
    read_length=100,
    coverage=8.0,
    paired_end=false,
    error_rate=0.005,
    seed=99UL,
  )
  let a = @src.simulate_reads(read_sim_test_reference(), config~)
  let b = @src.simulate_reads(read_sim_test_reference(), config~)
  assert_eq(a.fastq, b.fastq)
  assert_eq(a.n_reads, b.n_reads)
  assert_eq(a.n_errors, b.n_errors)
}
