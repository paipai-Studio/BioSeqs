///|
/// ResidueDepth demo - Calculate residue depth in protein structures.
/// Inspired by Biopython Bio.PDB.ResidueDepth.
///
/// Run: moon run examples/residue_depth_demo/main.mbt
fn main {
  println("=== ResidueDepth Demo ===")
  println("Inspired by Biopython Bio.PDB.ResidueDepth")
  println("")

  println("1. Create atoms and residues")
  let atom1_1 = @bio.RDAtom::new("CA", @bio.RDPoint3D::new(0.0, 0.0, 0.0), "C")
  let atom1_2 = @bio.RDAtom::new("N", @bio.RDPoint3D::new(-1.4, -0.2, 0.0), "N")
  let atom1_3 = @bio.RDAtom::new("C", @bio.RDPoint3D::new(1.5, 0.0, 0.0), "C")
  let atom1_4 = @bio.RDAtom::new("O", @bio.RDPoint3D::new(2.3, -1.0, 0.0), "O")
  
  let atom2_1 = @bio.RDAtom::new("CA", @bio.RDPoint3D::new(3.8, 0.0, 0.0), "C")
  let atom2_2 = @bio.RDAtom::new("N", @bio.RDPoint3D::new(2.5, 0.2, 0.0), "N")
  let atom2_3 = @bio.RDAtom::new("C", @bio.RDPoint3D::new(5.2, 0.0, 0.0), "C")
  
  let atom3_1 = @bio.RDAtom::new("CA", @bio.RDPoint3D::new(7.8, 0.0, 0.0), "C")
  let atom3_2 = @bio.RDAtom::new("N", @bio.RDPoint3D::new(6.5, -0.2, 0.0), "N")
  let atom3_3 = @bio.RDAtom::new("C", @bio.RDPoint3D::new(9.2, 0.0, 0.0), "C")
  
  let residue1 = @bio.RDResidue::new("ALA", 1, [atom1_1, atom1_2, atom1_3, atom1_4])
  let residue2 = @bio.RDResidue::new("GLY", 2, [atom2_1, atom2_2, atom2_3])
  let residue3 = @bio.RDResidue::new("SER", 3, [atom3_1, atom3_2, atom3_3])
  
  println("  Residue 1: ALA (4 atoms)")
  println("  Residue 2: GLY (3 atoms)")
  println("  Residue 3: SER (3 atoms)")
  println("")

  println("2. Analyze residue depth")
  let residues = [residue1, residue2, residue3]
  let results = @bio.bio_residue_depth_analyze(residues)
  
  println("  Depth Analysis Results:")
  println("  " + "-".repeat(60))
  println("  Residue | CA Depth | COM Depth | CA SASA | COM SASA")
  println("  " + "-".repeat(60))
  
  let mut i = 0
  while i < results.length() {
    let r = results[i]
    println("  " + r.resname + "(" + r.resseq.to_string() + ")" + 
            " | " + r.ca_depth.to_string() + " | " + 
            r.com_depth.to_string() + " | " + 
            r.ca_sasa.to_string() + " | " + 
            r.com_sasa.to_string())
    i = i + 1
  }
  println("")

  println("3. Find surface residues")
  let surface = @bio.bio_residue_depth_find_surface(residues, 3.0)
  println("  Surface residues (depth < 3.0): " + surface.length().to_string())
  
  let mut j = 0
  while j < surface.length() {
    println("    " + surface[j].resname + "(" + surface[j].resseq.to_string() + "): " + surface[j].ca_depth.to_string())
    j = j + 1
  }
  println("")

  println("4. Calculate average depth")
  let (avg_ca, avg_com) = @bio.bio_residue_depth_average(results)
  println("  Average CA depth: " + avg_ca.to_string())
  println("  Average COM depth: " + avg_com.to_string())
  println("")

  println("=== ResidueDepth Demo Complete ===")
}