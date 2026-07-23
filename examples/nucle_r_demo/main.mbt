///|
/// nucleR Demo - Nucleosome Positioning Analysis.
/// Inspired by Bioconductor nucleR package.

fn main {
  println("=== nucleR Demo ===")

  // 1. Create example track
  println("")
  println("1. Create example signal track:")
  let track = @bio.nuc_create_example_track()
  println("   Chromosome: " + track.nuc_get_chrom())
  println("   Length: " + track.nuc_get_length().to_string())
  println("   Signal values (first 10):")
  let mut i = 0
  while i < 10 {
    let pos = track.nuc_get_position(i)
    let val = track.nuc_get_signal(pos)
    println("     pos=" + pos.to_string() + ", signal=" + val.to_string())
    i = i + 1
  }

  // 2. Smooth signal
  println("")
  println("2. Smooth signal (window=100):")
  let smoothed = @bio.nuc_smooth_signal(track, 100)
  println("   Smoothed signal (first 10):")
  let mut j = 0
  while j < 10 {
    let pos = smoothed.nuc_get_position(j)
    let val = smoothed.nuc_get_signal(pos)
    println("     pos=" + pos.to_string() + ", signal=" + val.to_string())
    j = j + 1
  }

  // 3. Call nucleosomes with custom params
  println("")
  println("3. Call nucleosomes (peak_threshold=0.1, min_distance=100):")
  let params = @bio.NucCallParams::new().nuc_set_peak_threshold(0.1)
  let result = @bio.nuc_call_nucleosomes(track, params)
  println("   Nucleosome count: " + result.nuc_n_count().to_string())
  println("   Mean spacing: " + result.nuc_mean_spacing().to_string())
  println("   Mean occupancy: " + result.nuc_mean_occupancy().to_string())
  println("   Fraction well-positioned: " + result.nuc_frac_well_positioned().to_string())
  println("   Nucleosome details:")
  let mut k = 0
  while k < result.nuc_n_count() {
    let nuc = result.nuc_get(k)
    println("     " + k.to_string() +
      ": center=" + nuc.nuc_center().to_string() +
      ", width=" + nuc.nuc_width().to_string() +
      ", occupancy=" + nuc.nuc_occupancy().to_string() +
      ", well_positioned=" + nuc.nuc_is_well_positioned().to_string())
    k = k + 1
  }

  // 4. Compare positioning between two identical tracks
  println("")
  println("4. Compare nucleosome positioning (tolerance=50):")
  let result1 = @bio.nuc_call_nucleosomes(track, params)
  let result2 = @bio.nuc_call_nucleosomes(track, params)
  let comparison = @bio.nuc_compare_positioning(result1, result2, 50)
  println("   Shared: " + comparison.nuc_n_shared().to_string())
  println("   Gained: " + comparison.nuc_n_gained().to_string())
  println("   Lost: " + comparison.nuc_n_lost().to_string())
  println("   Fraction changed: " + comparison.nuc_frac_changed().to_string())
  println("   Direction: " + comparison.nuc_direction().to_string())

  println("")
  println("=== nucleR Demo Completed ===")
}