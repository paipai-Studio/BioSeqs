///|
/// illuminaio IDAT raw file parsing demo.
///
/// Parses a non-encrypted (Illumina Infinium methylation) IDAT file and prints
/// the per-probe bead type + summary intensities and the metadata strings.

///|
fn main raise {
  let raw = @src.base64_decode(@src.idat_sample_base64())
  let idat = @src.parse_idat(raw)
  println("version      : \{idat.version}")
  println("nFields      : \{idat.n_fields}")
  println("nSNPsRead    : \{idat.n_snps_read}")
  println("barcode      : \{idat.barcode}")
  println("chip type    : \{idat.chip_type}")
  println("probe  illumina_id  mean    sd   n_beads")
  for i = 0; i < idat.n_snps_read; i = i + 1 {
    println(
      "\{i}\t\{idat.illumina_id[i]}\t\{idat.mean[i]}\t\{idat.sd[i]}\t\{idat.n_beads[i]}",
    )
  }
}