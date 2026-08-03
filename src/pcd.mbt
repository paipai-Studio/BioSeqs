///|
/// PCD - PCD (Perl-Code-Decryption) mass spectrometry binary format parser.
/// Inspired by Biopython Bio.PCD.
///
/// PCD is a binary format used for storing mass spectrometry data,
/// particularly for proteomics and metabolomics applications.
/// It stores spectrum data as a binary file with a structured header.
///
/// The PCD format was developed by Marc Chabanon and colleagues for
/// storing high-throughput proteomics data efficiently.
///
/// Reference: Chabanon et al. (2007) Proteomics 7:863-868.

// ============================================================================
// PCD Data Structures
// ============================================================================

///|
/// Represents a single mass spectrum from a PCD file.
pub struct PcdSpectrum {
  scan_number : Int
  rt : Double
  mz_values : Array[Double]
  intensity_values : Array[Double]
  mut precursor_mz : Double?
  mut precursor_charge : Int?
}

///|
/// Construct a PcdSpectrum.
pub fn PcdSpectrum::new(
  scan_number : Int,
  rt? : Double = 0.0
) -> PcdSpectrum {
  PcdSpectrum::{
    scan_number,
    rt,
    mz_values: [],
    intensity_values: [],
    precursor_mz: None,
    precursor_charge: None,
  }
}

///|
/// Set precursor m/z.
pub fn PcdSpectrum::set_precursor_mz(self : PcdSpectrum, mz : Double?) -> Unit {
  self.precursor_mz = mz
}

///|
/// Set precursor charge.
pub fn PcdSpectrum::set_precursor_charge(self : PcdSpectrum, charge : Int?) -> Unit {
  self.precursor_charge = charge
}

///|
/// Add a peak to the spectrum.
pub fn PcdSpectrum::add_peak(self : PcdSpectrum, mz : Double, intensity : Double) -> Unit {
  self.mz_values.push(mz)
  self.intensity_values.push(intensity)
}

///|
/// Get the number of peaks in this spectrum.
pub fn PcdSpectrum::num_peaks(self : PcdSpectrum) -> Int {
  self.mz_values.length()
}

///|
/// Get the total intensity of this spectrum.
pub fn PcdSpectrum::total_intensity(self : PcdSpectrum) -> Double {
  let mut total = 0.0
  let mut i = 0
  while i < self.intensity_values.length() {
    total = total + self.intensity_values[i]
    i = i + 1
  }
  total
}

// ============================================================================
// PcdFile
// ============================================================================

///|
/// Represents a complete PCD file with metadata and spectra.
pub struct PcdFile {
  mut version : String
  mut num_spectra : Int
  mut num_peaks : Int
  spectra : Array[PcdSpectrum]
  metadata : Map[String, String]
}

///|
/// Construct an empty PcdFile.
pub fn PcdFile::new() -> PcdFile {
  PcdFile::{
    version: "1.0",
    num_spectra: 0,
    num_peaks: 0,
    spectra: [],
    metadata: Map([], capacity=16),
  }
}

///|
/// Set the file version.
pub fn PcdFile::set_version(self : PcdFile, version : String) -> Unit {
  self.version = version
}

///|
/// Add a spectrum to the file.
pub fn PcdFile::add_spectrum(self : PcdFile, spectrum : PcdSpectrum) -> Unit {
  self.spectra.push(spectrum)
}

///|
/// Set metadata key-value pair.
pub fn PcdFile::set_metadata(self : PcdFile, key : String, val : String) -> Unit {
  self.metadata[key] = val
}

///|
/// Update counts after parsing.
pub fn PcdFile::update_counts(self : PcdFile) -> Unit {
  self.num_spectra = self.spectra.length()
  let mut total = 0
  let mut i = 0
  while i < self.spectra.length() {
    total = total + self.spectra[i].mz_values.length()
    i = i + 1
  }
  self.num_peaks = total
}

///|
/// Get a spectrum by scan number.
pub fn PcdFile::get_spectrum(self : PcdFile, scan_number : Int) -> PcdSpectrum? {
  let mut i = 0
  while i < self.spectra.length() {
    if self.spectra[i].scan_number == scan_number {
      return Some(self.spectra[i])
    }
    i = i + 1
  }
  None
}

// ============================================================================
// PcdPeak
// ============================================================================

///|
/// Represents a peak in a mass spectrum.
pub struct PcdPeak {
  mz : Double
  intensity : Double
  charge : Int?
}

///|
/// Construct a PcdPeak.
pub fn PcdPeak::new(mz : Double, intensity : Double, charge? : Int = 0) -> PcdPeak {
  let charge_opt = if charge > 0 { Some(charge) } else { None }
  PcdPeak::{ mz, intensity, charge: charge_opt }
}

// ============================================================================
// PCD Parser (text-based PCD format)
// ============================================================================

///|
/// Parse a PCD text format (tab/space-delimited) into a PcdFile.
pub fn pcd_parse(content : String) -> PcdFile {
  let file = PcdFile::new()
  let lines = pcd_split_lines(content)
  let mut i = 0
  let mut current_spectrum : PcdSpectrum? = None
  while i < lines.length() {
    let line = lines[i]
    let trimmed = line.trim().to_string()
    if trimmed.length() == 0 {
      match current_spectrum {
        Some(s) => {
          file.add_spectrum(s)
          current_spectrum = None
        }
        None => ()
      }
      i = i + 1
      continue
    }
    if trimmed[0:1].to_owned() == "#" {
      pcd_parse_header_into(trimmed, file)
      i = i + 1
      continue
    }
    if trimmed.contains("SCAN") || trimmed.contains("scan") {
      match current_spectrum {
        Some(s) => file.add_spectrum(s)
        None => ()
      }
      let spectrum = pcd_parse_scan_line(trimmed)
      current_spectrum = Some(spectrum)
      i = i + 1
      continue
    }
    match current_spectrum {
      Some(s) => {
        let peak = pcd_parse_peak_line(trimmed)
        match peak {
          Some(p) => s.add_peak(p.mz, p.intensity)
          None => ()
        }
        current_spectrum = Some(s)
      }
      None => ()
    }
    i = i + 1
  }
  match current_spectrum {
    Some(s) => file.add_spectrum(s)
    None => ()
  }
  file.update_counts()
  file
}

///|
/// Parse a PCD header line and add to file.
fn pcd_parse_header_into(line : String, file : PcdFile) -> Unit {
  let content = line[1:line.length()].trim().to_string()
  let colon_pos = content.find(":")
  match colon_pos {
    Some(pos) => {
      let key = content[0:pos].trim().to_string()
      let val = content[pos + 1:content.length()].trim().to_string()
      file.set_metadata(key, val)
      if key.contains("version") || key.contains("Version") {
        file.set_version(val)
      }
    }
    None => {
      file.set_metadata(content, "")
    }
  }
}

///|
/// Parse a SCAN line (e.g., "SCAN 1 RT 12.5").
fn pcd_parse_scan_line(line : String) -> PcdSpectrum {
  let parts = pcd_split_whitespace(line)
  let mut scan_number = 0
  let mut rt = 0.0
  let mut pepmass_mz : Double? = None
  let mut pepmass_charge : Int? = None
  let mut i = 0
  while i < parts.length() {
    let part = parts[i]
    if part == "SCAN" || part == "scan" {
      if i + 1 < parts.length() {
        scan_number = pcd_parse_int(parts[i + 1])
      }
    }
    if part == "RT" || part == "rt" {
      if i + 1 < parts.length() {
        rt = pcd_parse_double(parts[i + 1])
      }
    }
    if part == "PEPMASS" || part == "pepmass" {
      if i + 1 < parts.length() {
        let pepmass_str = parts[i + 1]
        let paren_pos = pepmass_str.find("(")
        match paren_pos {
          Some(pp) => {
            let mz_str = pepmass_str[0:pp].to_owned()
            let charge_pos = pepmass_str.find(")")
            match charge_pos {
              Some(ce) => {
                let charge_str = pepmass_str[pp + 1:ce].to_owned()
                let charge = pcd_parse_int(charge_str)
                pepmass_charge = if charge > 0 { Some(charge) } else { None }
              }
              None => ()
            }
            pepmass_mz = Some(pcd_parse_double(mz_str))
          }
          None => {
            pepmass_mz = Some(pcd_parse_double(pepmass_str))
          }
        }
      }
    }
    i = i + 1
  }
  let spec = PcdSpectrum::new(scan_number, rt=rt)
  spec.set_precursor_mz(pepmass_mz)
  spec.set_precursor_charge(pepmass_charge)
  spec
}

///|
/// Parse a peak line (e.g., "500.5 10000" or "500.5 10000 2").
fn pcd_parse_peak_line(line : String) -> PcdPeak? {
  let parts = pcd_split_whitespace(line)
  if parts.length() < 2 { return None }
  let mz = pcd_parse_double(parts[0])
  let intensity = pcd_parse_double(parts[1])
  let charge = if parts.length() > 2 { pcd_parse_int(parts[2]) } else { 0 }
  Some(PcdPeak::new(mz, intensity, charge=charge))
}

// ============================================================================
// PCD Query Functions
// ============================================================================

///|
/// Get spectra within a retention time range.
pub fn pcd_spectra_by_rt(file : PcdFile, min_rt : Double, max_rt : Double) -> Array[PcdSpectrum] {
  let result : Array[PcdSpectrum] = Array::new()
  let mut i = 0
  while i < file.spectra.length() {
    let s = file.spectra[i]
    if s.rt >= min_rt && s.rt <= max_rt {
      result.push(s)
    }
    i = i + 1
  }
  result
}

///|
/// Get the total ion chromatogram (TIC) - total intensity per spectrum.
pub fn pcd_tic(file : PcdFile) -> Array[Double] {
  let tic : Array[Double] = Array::new()
  let mut i = 0
  while i < file.spectra.length() {
    let s = file.spectra[i]
    tic.push(s.total_intensity())
    i = i + 1
  }
  tic
}

///|
/// Get the base peak chromatogram (BPC) - maximum intensity per spectrum.
pub fn pcd_bpc(file : PcdFile) -> Array[Double] {
  let bpc : Array[Double] = Array::new()
  let mut i = 0
  while i < file.spectra.length() {
    let s = file.spectra[i]
    let mut max_val = 0.0
    let mut j = 0
    while j < s.intensity_values.length() {
      if s.intensity_values[j] > max_val {
        max_val = s.intensity_values[j]
      }
      j = j + 1
    }
    bpc.push(max_val)
    i = i + 1
  }
  bpc
}

///|
/// Find peaks above a threshold intensity in all spectra.
pub fn pcd_find_peaks_above(file : PcdFile, threshold : Double) -> Array[(Int, Double, Double)] {
  let peaks : Array[(Int, Double, Double)] = Array::new()
  let mut i = 0
  while i < file.spectra.length() {
    let s = file.spectra[i]
    let mut j = 0
    while j < s.mz_values.length() {
      if s.intensity_values[j] >= threshold {
        peaks.push((s.scan_number, s.mz_values[j], s.intensity_values[j]))
      }
      j = j + 1
    }
    i = i + 1
  }
  peaks
}

///|
/// Filter spectrum by m/z range.
pub fn pcd_spectrum_by_mz(spectrum : PcdSpectrum, min_mz : Double, max_mz : Double) -> PcdSpectrum {
  let result = PcdSpectrum::new(spectrum.scan_number, rt=spectrum.rt)
  result.set_precursor_mz(spectrum.precursor_mz)
  result.set_precursor_charge(spectrum.precursor_charge)
  let mut i = 0
  while i < spectrum.mz_values.length() {
    if spectrum.mz_values[i] >= min_mz && spectrum.mz_values[i] <= max_mz {
      result.add_peak(spectrum.mz_values[i], spectrum.intensity_values[i])
    }
    i = i + 1
  }
  result
}

///|
/// Get the total number of peaks across all spectra.
pub fn pcd_total_peaks(file : PcdFile) -> Int {
  file.num_peaks
}

///|
/// Get the number of spectra in the file.
pub fn pcd_num_spectra(file : PcdFile) -> Int {
  file.num_spectra
}

// ============================================================================
// PCD Writer
// ============================================================================

///|
/// Serialize a PcdFile to text format.
pub fn pcd_write(file : PcdFile) -> String {
  let mut out = ""
  out = out + "#version: " + file.version + "\n"
  let keys = file.metadata.keys().collect()
  let mut ki = 0
  while ki < keys.length() {
    let key = keys[ki]
    if key != "version" {
      out = out + "#" + key + ": " + file.metadata[key] + "\n"
    }
    ki = ki + 1
  }
  let mut i = 0
  while i < file.spectra.length() {
    let s = file.spectra[i]
    out = out + "SCAN " + s.scan_number.to_string() + " RT " + s.rt.to_string()
    match s.precursor_mz {
      Some(mz) => {
        out = out + " PEPMASS " + mz.to_string()
        match s.precursor_charge {
          Some(z) => out = out + "(" + z.to_string() + ")"
          None => ()
        }
      }
      None => ()
    }
    out = out + "\n"
    let mut j = 0
    while j < s.mz_values.length() {
      out = out + s.mz_values[j].to_string() + " " + s.intensity_values[j].to_string() + "\n"
      j = j + 1
    }
    out = out + "\n"
    i = i + 1
  }
  out
}

// ============================================================================
// Sample Data
// ============================================================================

///|
/// Generate a sample PCD text for testing.
pub fn pcd_sample() -> String {
  let lines = [
    "#version: 1.0",
    "#instrument: TestMassSpec",
    "#date: 2024-01-01",
    "SCAN 1 RT 1.0 PEPMASS 500.5(2)",
    "500.1 10000",
    "500.5 25000",
    "500.8 5000",
    "600.2 8000",
    "",
    "SCAN 2 RT 2.5 PEPMASS 600.3(3)",
    "500.1 5000",
    "600.3 30000",
    "600.5 15000",
    "700.1 2000",
    "",
    "SCAN 3 RT 5.0 PEPMASS 700.1(2)",
    "500.2 3000",
    "700.1 40000",
    "700.5 10000",
    "800.2 5000",
    "",
  ]
  let mut result = ""
  let mut i = 0
  while i < lines.length() {
    if i > 0 {
      result = result + "\n"
    }
    result = result + lines[i]
    i = i + 1
  }
  result
}

///|
/// Generate a summary of a PcdFile.
pub fn pcd_summary(file : PcdFile) -> String {
  let mut out = ""
  out = out + "PCD File Summary\n"
  out = out + "Version: " + file.version + "\n"
  out = out + "Number of spectra: " + file.num_spectra.to_string() + "\n"
  out = out + "Total peaks: " + file.num_peaks.to_string() + "\n"
  let keys = file.metadata.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let key = keys[i]
    out = out + "  " + key + ": " + file.metadata[key] + "\n"
    i = i + 1
  }
  out
}

// ============================================================================
// Private utility functions
// ============================================================================

///|
/// Split a string into lines.
fn pcd_split_lines(content : String) -> Array[String] {
  let lines : Array[String] = Array::new()
  let mut start = 0
  let mut i = 0
  while i < content.length() {
    let ch = content[i:i + 1].to_owned()
    if ch == "\n" || ch == "\r" {
      if start < i {
        lines.push(content[start:i].trim().to_string())
      }
      if ch == "\r" && i + 1 < content.length() && content[i + 1:i + 2].to_owned() == "\n" {
        i = i + 1
      }
      start = i + 1
    }
    i = i + 1
  }
  if start < content.length() {
    lines.push(content[start:content.length()].trim().to_string())
  }
  lines
}

///|
/// Split a string by whitespace (tabs or spaces).
fn pcd_split_whitespace(line : String) -> Array[String] {
  let cols : Array[String] = Array::new()
  let mut start = -1
  let mut i = 0
  while i < line.length() {
    let ch = line[i:i + 1].to_owned()
    if ch == " " || ch == "\t" {
      if start >= 0 {
        cols.push(line[start:i].to_string())
        start = -1
      }
    } else if start < 0 {
      start = i
    }
    i = i + 1
  }
  if start >= 0 {
    cols.push(line[start:line.length()].to_string())
  }
  cols
}

///|
/// Parse an integer from a string, returning 0 on failure.
fn pcd_parse_int(s : String) -> Int {
  let sv = s.to_string().trim()
  let n = sv.length()
  if n == 0 { return 0 }
  let mut result = 0
  let mut i = 0
  let zero = '0'.to_int()
  while i < n {
    let c = sv.unsafe_get(i).to_int()
    if c < zero || c > zero + 9 { break }
    result = result * 10 + (c - zero)
    i = i + 1
  }
  result
}

///|
/// Parse a double from a string, returning 0.0 on failure.
fn pcd_parse_double(s : String) -> Double {
  let sv = s.to_string().trim()
  if sv.length() == 0 { return 0.0 }
  let n = sv.length()
  let mut result = 0.0
  let mut decimal = 0
  let mut sign = 1.0
  let mut i = 0
  let zero = '0'.to_int()
  if sv.unsafe_get(0).to_int() == '-'.to_int() {
    sign = -1.0
    i = 1
  }
  while i < n {
    let c = sv.unsafe_get(i).to_int()
    if c < zero || c > zero + 9 { break }
    result = result * 10.0 + (c - zero).to_double()
    i = i + 1
  }
  if i < n && sv.unsafe_get(i) == '.' {
    i = i + 1
    while i < n {
      let c = sv.unsafe_get(i).to_int()
      if c < zero || c > zero + 9 { break }
      result = result * 10.0 + (c - zero).to_double()
      decimal = decimal + 1
      i = i + 1
    }
  }
  while decimal > 0 {
    result = result / 10.0
    decimal = decimal - 1
  }
  result * sign
}
