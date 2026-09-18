///|
/// BAM auxiliary tag + flag demo.
///
/// Parses a raw (BGZF-decompressed) BAM and prints each record's CIGAR string,
/// flags, sequence, quality, and auxiliary tags.

///|
fn main raise {
  let b64 = "QkFNAR8AAABASEQJVk46MS42CkBTUQlTTjpjaHIxCUxOOjEwMDAKAQAAAAUAAABjaHIxAOgDAABNAAAAAAAAAGQAAAAGPEkSAQADAAgAAAD//////////wAAAAByZWFkMQCAAAAAEkgSSAAAAAAAAAAATk1pAAAAAE1EWjgAQVNpCAAAAFhTQStXAAAAAAAAAMgAAAAGFEkSAgAABgUAAAD//////////wAAAAByZWFkMgBAAAAAEQAAABJIEAAAAAAATk1pAQAAAEFTZgAAkEBSR1pncnAxAFhBQkMDAAAAAQIDRQAAAAAAAAAsAQAABihJEgEAEAAGAAAA//////////8AAAAAcmVhZDMAYAAAAP///wAAAAAAAE1DYwFTQ0P6U0lz/f9VSVNg6g=="
  let bam = @src.parse_bam_base64(b64)
  for r in bam.records {
    let md = match r.get_tag("MD") {
      Some(@src.BamTagValue::Str(s)) => s
      _ => "-"
    }
    let nm = match r.get_tag("NM") {
      Some(@src.BamTagValue::Int(v)) => v.to_string()
      _ => "-"
    }
    println(
      "\{r.qname}\tflag=\{r.flag}\tcigar=\{r.cigarstring()}\tseq=\{r.seq}\tNM=\{nm}\tMD=\{md}\tdup=\{r.is_duplicate()}\trev=\{r.is_reverse()}",
    )
  }
}