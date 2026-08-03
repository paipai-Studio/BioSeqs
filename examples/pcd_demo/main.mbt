///|
/// PCD mass spectrometry format parser demo.
///
/// Demonstrates the PCD (Proteomics Computational Data) module for:
///   - Parsing PCD mass spectrometry data files
///   - Spectrum manipulation and querying
///   - Chromatogram generation (TIC, BPC)
///   - Peak filtering and metadata extraction
///   - Serialization (round-trip write/read)
///
/// Inspired by Biopython's Bio.PCD and related mass spec modules.
/// Run: moon run examples/pcd_demo
fn main {
  println("=== PCD Mass Spectrometry Parser Demo ===")
  println("Mass spectrometry data analysis tools")
  println("Reference: Biopython Bio.PCD")
  println("")

  // ------------------------------------------------------------------
  // 1. Parse PCD sample data
  // ------------------------------------------------------------------
  println("=== 1. Parsing PCD Sample Data ===")
  let content = @bio.pcd_sample()
  println("Input content (first 150 chars):")
  println(content[0:150].to_owned())
  println("")

  let file = @bio.pcd_parse(content)
  println("Parsed PCD file:")
  println("  Version: " + file.version)
  println("  Number of spectra: " + file.num_spectra.to_string())
  println("  Total peaks: " + file.num_peaks.to_string())
  println("")

  // ------------------------------------------------------------------
  // 2. List all spectra
  // ------------------------------------------------------------------
  println("=== 2. Spectrum Details ===")
  let mut i = 1
  while i <= file.num_spectra {
    match file.get_spectrum(i) {
      Some(spec) => {
        println("  Spectrum " + spec.scan_number.to_string() + ":")
        println("    Retention time: " + spec.rt.to_string() + " min")
        println("    Number of peaks: " + spec.num_peaks().to_string())
        println("    Total intensity: " + spec.total_intensity().to_string())
        match spec.precursor_mz {
          Some(mz) => println("    Precursor m/z: " + mz.to_string())
          None => println("    Precursor m/z: N/A")
        }
        match spec.precursor_charge {
          Some(z) => println("    Precursor charge: " + z.to_string())
          None => println("    Precursor charge: N/A")
        }

        // Show top 3 peaks
        println("    Top peaks (m/z, intensity):")
        let mut p = 0
        while p < spec.num_peaks() && p < 3 {
          let mz = spec.mz_values[p]
          let intensity = spec.intensity_values[p]
          println("      " + mz.to_string() + " -> " + intensity.to_string())
          p = p + 1
        }
        println("")
      }
      None => ()
    }
    i = i + 1
  }

  // ------------------------------------------------------------------
  // 3. Chromatogram generation
  // ------------------------------------------------------------------
  println("=== 3. Chromatogram Generation ===")

  // Total Ion Chromatogram (TIC)
  let tic = @bio.pcd_tic(file)
  println("Total Ion Chromatogram (TIC):")
  println("  Data points: " + tic.length().to_string())
  let mut max_tic = tic[0]
  let mut t = 1
  while t < tic.length() {
    if tic[t] > max_tic { max_tic = tic[t] }
    t = t + 1
  }
  println("  Max intensity: " + max_tic.to_string())
  println("")

  // Base Peak Chromatogram (BPC)
  let bpc = @bio.pcd_bpc(file)
  println("Base Peak Chromatogram (BPC):")
  println("  Data points: " + bpc.length().to_string())
  println("")

  // ------------------------------------------------------------------
  // 4. Spectrum querying
  // ------------------------------------------------------------------
  println("=== 4. Spectrum Querying ===")

  // Get spectrum by scan number
  println("Looking up spectrum #1:")
  match file.get_spectrum(1) {
    Some(spec) => {
      println("  Found spectrum with " + spec.num_peaks().to_string() + " peaks")
    }
    None => println("  Not found")
  }

  // Filter by m/z range
  println("")
  println("Filtering spectrum #1 for m/z in range 500-600:")
  match file.get_spectrum(1) {
    Some(spec) => {
      let filtered = @bio.pcd_spectrum_by_mz(spec, 500.0, 600.0)
      println("  Peaks after filtering: " + filtered.num_peaks().to_string())
    }
    None => println("  Spectrum not found")
  }
  println("")

  // Total peaks in file
  println("Total peaks across all spectra: " + @bio.pcd_total_peaks(file).to_string())
  println("Number of spectra: " + @bio.pcd_num_spectra(file).to_string())
  println("")

  // ------------------------------------------------------------------
  // 5. Round-trip serialization
  // ------------------------------------------------------------------
  println("=== 5. Round-trip Serialization ===")
  let written = @bio.pcd_write(file)
  println("Serialized PCD content (first 200 chars):")
  println(written[0:200].to_owned())
  println("")

  // Parse back and verify
  let file2 = @bio.pcd_parse(written)
  println("Re-parsed PCD file:")
  println("  Spectra count match: " + (file.num_spectra == file2.num_spectra).to_string())
  println("  Peak count match: " + (file.num_peaks == file2.num_peaks).to_string())
  println("")

  // ------------------------------------------------------------------
  // 6. Summary report
  // ------------------------------------------------------------------
  println("=== 6. Summary Report ===")
  let summary = @bio.pcd_summary(file)
  println(summary)

  println("=== PCD Demo Complete ===")
}
