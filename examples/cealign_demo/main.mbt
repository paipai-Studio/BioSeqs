///|
fn main {
  println("=== Biopython Bio.PDB.cealign Demo ===")
  let (reference, mobile) = @src.cealign_sample_structures()
  let original_mobile_coordinate = mobile.get_atoms()[0].coord

  println("\n1. Configure a CE aligner")
  let aligner = @src.CEAligner::new() catch {
    _ => abort("failed to create CEAligner")
  }
  let configured = aligner.set_reference(reference) catch {
    _ => abort("failed to set CE reference structure")
  }
  println(
    "  Window size / maximum gap: " +
    configured.window_size.to_string() +
    " / " +
    configured.max_gap.to_string(),
  )
  println(
    "  Reference guide atoms: " + configured.reference_length().to_string(),
  )

  println("\n2. Search CE fragment paths and superimpose structures")
  let result = configured.align(mobile, transform=true, final_optimization=true) catch {
    _ => abort("failed to align sample structures")
  }
  println("  " + result.summary())
  println(
    "  Aligned guide atoms / fragments: " +
    result.path.aligned_length().to_string() +
    " / " +
    result.path.fragment_count.to_string(),
  )
  println(
    "  RMSD / CE Z-score: " +
    result.rmsd.to_string() +
    " / " +
    result.path.z_score.to_string(),
  )
  println(
    "  Reference / mobile coverage: " +
    result.reference_coverage().to_string() +
    " / " +
    result.mobile_coverage().to_string(),
  )
  println(
    "  Final local optimization applied: " +
    result.optimization_applied.to_string(),
  )

  println("\n3. Inspect the selected monotonic path")
  let pairs = result.path.aligned_pairs()
  let first = pairs[0]
  let last = pairs[pairs.length() - 1]
  println(
    "  First reference/mobile index: " +
    first.0.to_string() +
    " / " +
    first.1.to_string(),
  )
  println(
    "  Last reference/mobile index: " +
    last.0.to_string() +
    " / " +
    last.1.to_string(),
  )

  println("\n4. Verify immutable all-atom transformation")
  let unchanged_mobile_coordinate = mobile.get_atoms()[0].coord
  let transformed_coordinate = result.structure.get_atoms()[0].coord
  let reference_coordinate = reference.get_atoms()[0].coord
  println(
    "  Input mobile coordinate unchanged: " +
    (unchanged_mobile_coordinate.distance(original_mobile_coordinate) < 1.0e-12).to_string(),
  )
  println(
    "  Transformed first guide distance to reference: " +
    transformed_coordinate.distance(reference_coordinate).to_string(),
  )
  println(
    "  Transformed atom count: " +
    result.structure.get_atoms().length().to_string(),
  )

  println("\n=== Demo Complete ===")
}
