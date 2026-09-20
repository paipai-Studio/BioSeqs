///|
/// CSI index demo: parse a CSI index and show its binning parameters.

///|
let demo_csi_base64 = "Q1NJAQ4AAAAGAAAAJgAAAAIAAAABAAAAAgAAAAAAAAAjAAAAAAAAAAoAAABjaHIxAGNocjIAAgAAAAIAAABJkgAAeQAAAAAAAAABAAAAeQAAAAAAAACjAAAAAAAAAEqSBAAAAAAAAAAAAAIAAAB5AAAAAAAAAKMAAAAAAAAAAgAAAAAAAAAAAAAAAAAAAAIAAABJkgAAowAAAAAAAAABAAAAowAAAAAAAAAAAJ4AAAAAAEqSBAAAAAAAAAAAAAIAAACjAAAAAAAAAAAAngAAAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

///|
fn main raise {
  let idx = @src.parse_csi(@src.base64_decode(demo_csi_base64))
  println("CSI: min_shift=\{idx.min_shift} n_lvls=\{idx.n_lvls} n_ref=\{idx.references.length()} n_no_coor=\{idx.n_no_coor}")
  for i = 0; i < idx.references.length(); i = i + 1 {
    let r = idx.references[i]
    println("  ref \{i}: \{r.bins.length()} bins")
    for b in r.bins {
      println("    bin=\{b.bin} loffset=\{b.loffset} chunks=\{b.chunks.length()}")
    }
  }
  println("reg2bins(100, 20000, 14, 6) = \{@src.csi_reg2bins(100, 20000, 14, 6)}")
  let q = @src.csi_query(idx, 0, 0, 1000)
  println("query chr1 [0,1000): min_off=\{q.min_off} chunks=\{q.chunks.length()}")
}
