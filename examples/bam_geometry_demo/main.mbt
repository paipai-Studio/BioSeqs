///|
/// BAM record alignment-geometry demo.
///
/// Parses a small BAM (records with M/I/D/N/S/=/X CIGAR operations) and shows
/// each record's reference span and covered reference positions.

///|
let demo_bam_base64 = "QkFNAS0AAABASEQJVk46MS42CVNPOmNvb3JkaW5hdGUKQFNRCVNOOmNocjEJTE46MTAwMAoBAAAABQAAAGNocjEA6AMAAC0AAAAAAAAACgAAAAM8SRIBAAAABAAAAP//////////AAAAAHIxAEAAAAASSP////87AAAAAAAAAAoAAAADPEkSAwAAAAgAAAD//////////wAAAAByMgBAAAAAIQAAACAAAAASSBJI//////////83AAAAAAAAAAoAAAADPEkSAwAAAAUAAAD//////////wAAAAByMwAwAAAAEgAAACAAAAASSED//////zwAAAAAAAAACgAAAAM8SRIEAAAABgAAAP//////////AAAAAHI0ACAAAAARAAAAEAAAACQAAAASSBL///////8+AAAAAAAAAAoAAAADPEkSBAAAAAcAAAD//////////wAAAAByNQAkAAAAMAAAABMAAAAgAAAAEkgSgP////////9AAAAAAAAAAAoAAAADPEkSBQAAAAYAAAD//////////wAAAAByNgAUAAAAIAAAABIAAAAQAAAAJAAAABJIEv///////zsAAAAAAAAACgAAAAM8SRIDAAAACAAAAP//////////AAAAAHI3AFcAAAAYAAAAIAAAABJIEkj//////////w=="

///|
fn main raise {
  let bam = @src.parse_bam_base64(demo_bam_base64)
  for r in bam.records {
    println(
      "span \{r.reference_start()}-\{r.reference_end()}  " +
      "ref_pos=\{r.get_reference_positions(false)}  " +
      "overlap(8,20)=\{r.get_overlap(8, 20)}",
    )
  }
}