///|
/// BAM pileup demo: coverage columns over a small region.

///|
let demo_bam_base64 = "QkFNAS0AAABASEQJVk46MS42CVNPOmNvb3JkaW5hdGUKQFNRCVNOOmNocjEJTE46MTAwMAoBAAAABQAAAGNocjEA6AMAADcAAAAAAAAACAAAAAY8SRICAAAABgAAAP//////////AAAAAHJlYWQ1ACQAAABAAAAAEkgS////////MwAAAAAAAAAKAAAABjxJEgEAAAAGAAAA//////////8AAAAAcmVhZDEAYAAAABJIEv///////zAAAAAAAAAACgAAAAY8SRIBAAAABAAAAP//////////AAAAAHJlYWQ2AEAAAAASSP////87AAAAAAAAAAsAAAAGPEkSAwAAAAYAAAD//////////wAAAAByZWFkMwAwAAAAEQAAACAAAAASSBL///////86AAAAAAAAAAwAAAAGPEkSAwAAAAUAAAD//////////wAAAAByZWFkMgAgAAAAEgAAADAAAAASSBD//////zgAAAAAAAAADQAAAAY8SRIDAAAABAAAAP//////////AAAAAHJlYWQ0ACAAAAATAAAAIAAAABJI/////w=="

///|
fn main raise {
  let bam = @src.parse_bam_base64(demo_bam_base64)
  for col in bam.pileup("chr1", 0, 25) {
    let parts = StringBuilder()
    for pr in col.pileups {
      let name = bam.records[pr.record_index].qname
      match pr.query_position {
        Some(q) => {
          parts.write_string("\{name}:q\{q}")
          if pr.indel != 0 {
            parts.write_string("(indel=\{pr.indel})")
          }
        }
        None => parts.write_string("\{name}:del")
      }
      parts.write_string(" ")
    }
    println("pos=\{col.reference_pos} n=\{col.nsegments}  \{parts}")
  }
}