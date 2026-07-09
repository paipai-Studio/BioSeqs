// Performance benchmark for the MoonBit Seq implementation.
//
// Generates a large synthetic DNA sequence and runs each Seq operation
// repeatedly, printing the elapsed time (us) for each. Run with:
//   moon run cmd/bench

///|
/// Generate a deterministic pseudo-random DNA sequence of length `n`
/// using a simple LCG so results are reproducible across runs.
fn gen_dna(n : Int) -> String {
  let buf = StringBuilder::new()
  let bases = ['A', 'C', 'G', 'T']
  let mut state = 12345
  for i = 0; i < n; i = i + 1 {
    state = (state * 1103515245 + 12345) & 0x7fffffff
    let idx = state % 4
    buf.write_char(bases[idx])
  }
  buf.to_string()
}

///|
/// Measure how long `f` takes when called `iters` times; returns total us.
fn bench(label : String, iters : Int, f : () -> Unit) -> Double {
  let ts = @mb.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    f()
  }
  let elapsed_us = @mb.monotonic_clock_end(ts)
  let per_iter = elapsed_us / iters.to_double()
  println(
    "\{label}: \{elapsed_us} us total, \{per_iter} us/iter (\{iters} iters)",
  )
  elapsed_us
}

///|
fn main {
  let dna = gen_dna(100_000)
  let rna = gen_dna(100_000).replace_all(old="T", new="U")
  let s = @lib.Seq::new(dna)
  let sr = @lib.Seq::new(rna)
  let iters_small = 20
  let iters_find = 200

  println("=== MoonBit Seq benchmark ===")
  println("sequence length: \{dna.char_length()} chars")
  println("iterations (ops): \{iters_small}")
  println("")

  // complement
  ignore(bench("complement", iters_small, fn() { ignore(s.complement()) }))
  ignore(
    bench("complement_rna", iters_small, fn() { ignore(s.complement_rna()) }),
  )
  ignore(
    bench("reverse_complement", iters_small, fn() {
      ignore(s.reverse_complement())
    }),
  )
  ignore(bench("reverse", iters_small, fn() { ignore(s.reverse()) }))

  // transcribe
  ignore(bench("transcribe", iters_small, fn() { ignore(s.transcribe()) }))
  ignore(
    bench("back_transcribe", iters_small, fn() { ignore(sr.back_transcribe()) }),
  )

  // translate
  ignore(
    bench("translate", iters_small, fn() {
      ignore(s.translate()) catch {
        _ => ()
      }
    }),
  )

  // upper / lower
  ignore(bench("upper", iters_small, fn() { ignore(s.upper()) }))
  ignore(bench("lower", iters_small, fn() { ignore(s.lower()) }))

  // count / find
  ignore(bench("count AAA", iters_find, fn() { ignore(s.count("AAA")) }))
  ignore(
    bench("count_overlap AAA", iters_find, fn() {
      ignore(s.count_overlap("AAA"))
    }),
  )
  ignore(bench("find ACGT", iters_find, fn() { ignore(s.find("ACGT")) }))
  ignore(bench("rfind ACGT", iters_find, fn() { ignore(s.rfind("ACGT")) }))

  // startswith / endswith / contains
  ignore(
    bench("startswith ATG", iters_find, fn() { ignore(s.startswith("ATG")) }),
  )
  ignore(bench("endswith TAA", iters_find, fn() { ignore(s.endswith("TAA")) }))
  ignore(bench("contains GGG", iters_find, fn() { ignore(s.contains("GGG")) }))

  // slice
  ignore(
    bench("slice 0..1000", iters_find, fn() {
      ignore(s.slice(start=0, end=1000))
    }),
  )

  // replace
  ignore(
    bench("replace A->N", iters_small, fn() {
      ignore(s.replace(old="A", new="N"))
    }),
  )

  println("")
  println("done.")
}
