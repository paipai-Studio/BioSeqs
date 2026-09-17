///|
/// BCF (binary VCF) Demo.
///
/// Demonstrates the BCF2 parser: decoding the header name dictionaries and
/// the per-record typed INFO/FORMAT values of a real (htslib-generated) BCF
/// byte stream.

///|
fn main raise {
  let bcf = @src.parse_bcf(
    @src.base64_decode(
      "QkNGAgJLAgAAIyNmaWxlZm9ybWF0PVZDRnY0LjIKIyNGSUxURVI9PElEPVBBU1MsRGVzY3JpcHRpb249IkFsbCBmaWx0ZXJzIHBhc3NlZCIsSURYPTA+CiMjY29udGlnPTxJRD1jaHIxLGxlbmd0aD0xMDAwLElEWD0wPgojI2NvbnRpZz08SUQ9Y2hyMixsZW5ndGg9NTAwLElEWD0xPgojI0ZJTFRFUj08SUQ9cTEwLERlc2NyaXB0aW9uPSJRdWFsaXR5ID4gMTAiLElEWD0xPgojI0lORk89PElEPURQLE51bWJlcj0xLFR5cGU9SW50ZWdlcixEZXNjcmlwdGlvbj0iRGVwdGgiLElEWD0yPgojI0lORk89PElEPUFGLE51bWJlcj1BLFR5cGU9RmxvYXQsRGVzY3JpcHRpb249IkFsbGVsZSBmcmVxIixJRFg9Mz4KIyNJTkZPPTxJRD1GTEFHWCxOdW1iZXI9MCxUeXBlPUZsYWcsRGVzY3JpcHRpb249IkEgZmxhZyIsSURYPTQ+CiMjRk9STUFUPTxJRD1HVCxOdW1iZXI9MSxUeXBlPVN0cmluZyxEZXNjcmlwdGlvbj0iR2Vub3R5cGUiLElEWD01PgojI0ZPUk1BVD08SUQ9QUQsTnVtYmVyPVIsVHlwZT1JbnRlZ2VyLERlc2NyaXB0aW9uPSJBbGxlbGUgZGVwdGhzIixJRFg9Nj4KI0NIUk9NCVBPUwlJRAlSRUYJQUxUCVFVQUwJRklMVEVSCUlORk8JRk9STUFUCVMxCgAwAAAACgAAAAAAAABlAAAAAQAAAAAA9EEDAAIAAQAAAjdyczEXQRdDEQARAhEMEQMVAACAPhEEABEFIQIEEQYhCgIwAAAACwAAAAEAAADHAAAAAQAAAAEAgH8CAAMAAQAAAgcXRxdUF0ERARECEQURAyXNzMw9zcxMPhEFIQQGEQYxAAMC",
    ),
  )

  println("=== header ===")
  println("contigs: \{bcf.header.contig_names[0]}, \{bcf.header.contig_names[1]}")
  println("filters: \{bcf.header.id_name(bcf.header.filter_ids[0])}, \{bcf.header.id_name(bcf.header.filter_ids[1])}")
  println("samples: \{bcf.header.sample_names[0]}")

  println("\n=== records ===")
  for i = 0; i < bcf.records.length(); i = i + 1 {
    let r = bcf.records[i]
    println(
      "\{r.chrom}:\{r.pos} \{r.ref_seq}>\{r.alt[0]} id=\{r.id}",
    )
    for j = 0; j < r.info.length(); j = j + 1 {
      let inf = r.info[j]
      let name = bcf.header.id_name(inf.key_id)
      match inf.values[0] {
        BcfInt(v) => println("  INFO \{name} = \{v} (int)")
        BcfFloat(x) => println("  INFO \{name} = \{x} (float)")
        BcfFlag => println("  INFO \{name} (flag)")
        _ => println("  INFO \{name}")
      }
    }
    // decode GT from the first FORMAT field
    for j = 0; j < r.formats.length(); j = j + 1 {
      let fm = r.formats[j]
      let name = bcf.header.id_name(fm.key_id)
      if name == "GT" {
        match (fm.values[0], fm.values[1]) {
          (BcfInt(a), BcfInt(b)) =>
            println("  FORMAT GT = \{@src.bcf_gt_allele(a)}/\{@src.bcf_gt_allele(b)}")
          _ => println("  FORMAT GT = ?")
        }
      }
    }
  }
}