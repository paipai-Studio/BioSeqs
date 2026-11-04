// SeqIO test companion: runs the same battery of SeqIO test cases as
// tests/python_seqio_reference.py and emits results between case markers,
// matching BioPython's output exactly so the two can be diffed.
//
// Run:  moon run cmd/seqio_main

///|
fn main {
  seqio_test()
}
