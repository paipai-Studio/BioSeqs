///|
/// quaternion_superimposer demo - Protein structure superposition using Quaternion method.
/// Inspired by Biopython Bio.PDB.Superimposer (Quaternion method).

fn main {
  println("=== quaternion_superimposer Demo: Quaternion Structure Alignment ===")
  println("")

  // Step 1: Two identical structures
  println("--- Step 1: Identical Structures ---")
  let a1 = [
    @src.QsPoint::new(1.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 1.0, 0.0),
    @src.QsPoint::new(0.0, 0.0, 1.0),
  ]
  let b1 = [
    @src.QsPoint::new(1.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 1.0, 0.0),
    @src.QsPoint::new(0.0, 0.0, 1.0),
  ]
  let result1 = @src.qs_superimpose(a1, b1)
  println("RMSD: " + result1.rmsd.to_string())
  println("")

  // Step 2: Rotated structure
  println("--- Step 2: Rotated Structure ---")
  let a2 = [
    @src.QsPoint::new(1.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 1.0, 0.0),
    @src.QsPoint::new(0.0, 0.0, 1.0),
  ]
  let b2 = [
    @src.QsPoint::new(0.0, 1.0, 0.0),
    @src.QsPoint::new(-1.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 0.0, 1.0),
  ]
  let result2 = @src.qs_superimpose(a2, b2)
  println("RMSD after superposition: " + result2.rmsd.to_string())
  println("")

  // Step 3: Translated structure
  println("--- Step 3: Translated Structure ---")
  let a3 = [
    @src.QsPoint::new(1.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 1.0, 0.0),
    @src.QsPoint::new(0.0, 0.0, 1.0),
  ]
  let b3 = [
    @src.QsPoint::new(2.0, 1.0, 1.0),
    @src.QsPoint::new(1.0, 2.0, 1.0),
    @src.QsPoint::new(1.0, 1.0, 2.0),
  ]
  let result3 = @src.qs_superimpose(a3, b3)
  println("RMSD after superposition: " + result3.rmsd.to_string())
  println("")

  // Step 4: Apply transformation
  println("--- Step 4: Apply Transformation ---")
  let points = [
    @src.QsPoint::new(1.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 1.0, 0.0),
  ]
  let transformed = @src.qs_apply_transform(points, result3.rotation, result3.translation)
  println("Transformed point 0: (" + transformed[0].x.to_string() + ", " + transformed[0].y.to_string() + ", " + transformed[0].z.to_string() + ")")
  println("")

  // Step 5: Quaternion operations
  println("--- Step 5: Quaternion Operations ---")
  let q = @src.quaternion_identity()
  println("Identity quaternion: w=" + q.w.to_string() + " x=" + q.x.to_string() + " y=" + q.y.to_string() + " z=" + q.z.to_string())
  let qc = @src.quaternion_conjugate(q)
  println("Conjugate: w=" + qc.w.to_string() + " x=" + qc.x.to_string() + " y=" + qc.y.to_string() + " z=" + qc.z.to_string())
  let mat = @src.quaternion_to_matrix(q)
  println("Rotation matrix from identity quaternion:")
  println("  [" + mat[0][0].to_string() + " " + mat[0][1].to_string() + " " + mat[0][2].to_string() + "]")
  println("  [" + mat[1][0].to_string() + " " + mat[1][1].to_string() + " " + mat[1][2].to_string() + "]")
  println("  [" + mat[2][0].to_string() + " " + mat[2][1].to_string() + " " + mat[2][2].to_string() + "]")
  println("")

  println("=== quaternion_superimposer Demo Complete ===")
}