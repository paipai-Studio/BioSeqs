///|
/// Demo for polypose module: multi-conformer PDB structures (NMR ensembles).

fn main {
  println("=== polypose Demo: Multi-Conformer Structures ===")
  println("")

  // Example 1: Build two conformers
  println("--- Step 1: Build Conformers ---")
  let a1 = @src.PpyAtom::new("CA", 1.0, 0.0, 0.0)
  let a2 = @src.PpyAtom::new("CA", 2.0, 0.0, 0.0)
  let a3 = @src.PpyAtom::new("CA", 3.0, 0.0, 0.0)
  let r1 = @src.PpyResidue::new("ALA", 1, [a1])
  let r2 = @src.PpyResidue::new("GLY", 2, [a2])
  let r3 = @src.PpyResidue::new("LEU", 3, [a3])
  let c1 = @src.PpyConformer::new("model_1", [r1, r2, r3])

  let a4 = @src.PpyAtom::new("CA", 1.1, 0.1, 0.0)
  let a5 = @src.PpyAtom::new("CA", 2.1, 0.0, 0.1)
  let a6 = @src.PpyAtom::new("CA", 2.9, -0.1, 0.0)
  let r4 = @src.PpyResidue::new("ALA", 1, [a4])
  let r5 = @src.PpyResidue::new("GLY", 2, [a5])
  let r6 = @src.PpyResidue::new("LEU", 3, [a6])
  let c2 = @src.PpyConformer::new("model_2", [r4, r5, r6])
  println("Conformer 1: " + c1.conformer_id + ", residues=" + c1.residues.length().to_string())
  println("Conformer 2: " + c2.conformer_id + ", residues=" + c2.residues.length().to_string())
  println("")

  // Example 2: Build a multi-conformer structure
  println("--- Step 2: Build Structure ---")
  let st0 = @src.ppy_new("nmr_bundle")
  let st1 = @src.ppy_add_conformer(st0, c1)
  let st = @src.ppy_add_conformer(st1, c2)
  println("Structure: " + st.name + ", conformers=" + @src.ppy_n_conformers(st).to_string())
  println("")

  // Example 3: Extract Cα coordinates
  println("--- Step 3: Cα Coordinates ---")
  let coords1 = @src.ppy_coords(c1, "CA")
  let coords2 = @src.ppy_coords(c2, "CA")
  println("Conformer 1 Cα coords: " + coords1.length().to_string() + " points")
  for i in 0..<coords1.length() {
    println("  ALA " + (i+1).to_string() + ": (" + coords1[i].x.to_string() + ", " + coords1[i].y.to_string() + ", " + coords1[i].z.to_string() + ")")
  }
  println("")

  // Example 4: Pairwise RMSD matrix
  println("--- Step 4: Pairwise RMSD Matrix ---")
  match @src.ppy_pairwise_rmsds(st) {
    Some(mat) => {
      println("RMSD matrix: " + mat.length().to_string() + " x " + mat[0].length().to_string())
      for i in 0..<mat.length() {
        var line = ""
        for j in 0..<mat[i].length() {
          line = line + " " + mat[i][j].to_string()
        }
        println(line)
      }
    }
    None => println("Coordinate mismatch")
  }
  println("")

  // Example 5: Convergence statistics
  println("--- Step 5: Convergence Statistics ---")
  match @src.ppy_convergence(st) {
    Some(c) => {
      println("Conformers: " + c.n_conformers.to_string())
      println("Cα atoms per conformer: " + c.n_calpha.to_string())
      println("Mean RMSD: " + c.mean_rmsd.to_string())
      println("Min RMSD: " + c.min_rmsd.to_string())
      println("Max RMSD: " + c.max_rmsd.to_string())
      println("Median RMSD: " + c.median_rmsd.to_string())
    }
    None => println("No convergence data")
  }
  println("")

  // Example 6: Superpose and average structure
  println("--- Step 6: Superpose & Average Structure ---")
  match @src.ppy_mean_structure(st) {
    Some(avg) => {
      println("Average structure: " + avg.length().to_string() + " Cα points")
      for i in 0..<avg.length() {
        println("  " + (i+1).to_string() + ": (" + avg[i].x.to_string() + ", " + avg[i].y.to_string() + ", " + avg[i].z.to_string() + ")")
      }
    }
    None => println("Could not compute average structure")
  }
  println("")

  // Example 7: altLoc detection
  println("--- Step 7: altLoc Alternate Conformers ---")
  let alt_a = @src.PpyAtom::new("CA", 1.0, 0.0, 0.0, altloc="A")
  let alt_b = @src.PpyAtom::new("CA", 1.0, 0.1, 0.0, altloc="B")
  let alt_r = @src.PpyResidue::new("MET", 1, [alt_a, alt_b])
  let alt_c = @src.PpyConformer::new("alt_conf", [alt_r])
  let alt_ids = @src.ppy_altloc_ids(alt_r)
  println("altLoc IDs: " + alt_ids.to_string())
  println("Residues with altLoc: " + @src.ppy_altloc_residue_count(alt_c).to_string())
  println("")

  println("=== polypose Demo Complete ===")
}
