///|
/// Performance measurement tests for FASTA index operations.

fn report_us(name : String, total : Double, iters : Int) -> Unit {
  let per_call = total / iters.to_double() / 1000.0
  println(
    name +
    " x" +
    iters.to_string() +
    ": " +
    (total / 1000000.0).to_string() +
    "ms (" +
    per_call.to_string() +
    " us/call)",
  )
}

fn report_ms(name : String, total : Double, iters : Int) -> Unit {
  let per_call = total / iters.to_double() / 1000000.0
  println(
    name +
    " x" +
    iters.to_string() +
    ": " +
    (total / 1000000.0).to_string() +
    "ms (" +
    per_call.to_string() +
    " ms/call)",
  )
}

fn perf_gen_large_fasta(num_seqs : Int, seq_len : Int) -> String {
  let buf = StringBuilder::new()
  let mut i = 0
  while i < num_seqs {
    buf.write_char('>')
    buf.write_string("seq")
    buf.write_string(i.to_string())
    buf.write_char('\n')
    
    let mut j = 0
    while j < seq_len {
      let c = match j % 4 {
        0 => 'A'
        1 => 'T'
        2 => 'G'
        _ => 'C'
      }
      buf.write_char(c)
      j = j + 1
    }
    buf.write_char('\n')
    
    i = i + 1
  }
  buf.to_string()
}

///|
test "perf_measure_faidx" {
  let fasta_100_1000 = perf_gen_large_fasta(100, 1000)
  let fasta_1000_100 = perf_gen_large_fasta(1000, 100)
  
  let idx_100_1000 = FastaIndex::build(fasta_100_1000)
  
  // Build index x10 for 100 sequences of 1000bp each
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < 10; i = i + 1 {
    let _ = FastaIndex::build(fasta_100_1000)
  }
  report_ms("  faidx_build_100x1000", @bench.monotonic_clock_end(t0), 10)
  
  // Build index x10 for 1000 sequences of 100bp each
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < 10; i = i + 1 {
    let _ = FastaIndex::build(fasta_1000_100)
  }
  report_ms("  faidx_build_1000x100", @bench.monotonic_clock_end(t0), 10)
  
  // Fetch 100bp subsequence x100
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < 100; i = i + 1 {
    let _ = idx_100_1000.fetch("seq50", 1, 100)
  }
  report_us("  faidx_fetch_100", @bench.monotonic_clock_end(t0), 100)
  
  // Fetch full sequence x100
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < 100; i = i + 1 {
    let _ = idx_100_1000.fetch("seq50", 1, 1000)
  }
  report_us("  faidx_fetch_full", @bench.monotonic_clock_end(t0), 100)
  
  // Fetch random 50bp subsequences x100
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < 100; i = i + 1 {
    let _ = idx_100_1000.fetch("seq50", i + 1, i + 50)
  }
  report_us("  faidx_fetch_random", @bench.monotonic_clock_end(t0), 100)
  
  // Fetch from multiple sequences x100
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < 100; i = i + 1 {
    let _ = idx_100_1000.fetch("seq" + i.to_string(), 1, 100)
  }
  report_us("  faidx_fetch_multi", @bench.monotonic_clock_end(t0), 100)
  
  assert_true(true)
}