///|
/// pyfaidx Sequence slicing demo.
///
/// Shows the coordinate-tracked slice model: forward/reverse slices, integer
/// indexing, complement/reverse-complement, and orientation.

///|
fn main raise {
  let x = @src.pyfaidx_seq_new("chr1", "ATCGTA", Some(1), Some(6), false)

  println("x          : \{@src.pyfaidx_seq_fancy_name(x)}  \{x.seq}")
  println("x[:3]      : \{@src.pyfaidx_seq_fancy_name(@src.pyfaidx_seq_slice(x, None, Some(3), None))}")
  println("x[3:]      : \{@src.pyfaidx_seq_fancy_name(@src.pyfaidx_seq_slice(x, Some(3), None, None))}")
  println("x[1:-1]    : \{@src.pyfaidx_seq_fancy_name(@src.pyfaidx_seq_slice(x, Some(1), Some(-1), None))}")
  println("x[::-1]    : \{@src.pyfaidx_seq_fancy_name(@src.pyfaidx_seq_slice(x, None, None, Some(-1)))}")
  println("x[::2]     : \{@src.pyfaidx_seq_fancy_name(@src.pyfaidx_seq_slice(x, None, None, Some(2)))}")

  let neg = @src.pyfaidx_seq_neg(x)
  println("-x         : \{@src.pyfaidx_seq_fancy_name(neg)}  \{neg.seq}")
  println("orientation of -x : \{@src.pyfaidx_seq_orientation(neg)}")
}