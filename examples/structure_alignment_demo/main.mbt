///|
/// StructureAlignment demo - Multiple protein structure alignment.
/// Inspired by Biopython Bio.PDB.StructureAlignment.
///
/// Run: moon run examples/structure_alignment_demo/main.mbt
fn main {
  println("=== StructureAlignment Demo ===")
  println("Inspired by Biopython Bio.PDB.StructureAlignment")
  println("")

  println("1. Create test structures")
  let residues1 = [
    @bio.SAResidue::new("ALA", 1, @bio.SAPoint3D::new(0.0, 0.0, 0.0)),
    @bio.SAResidue::new("GLY", 2, @bio.SAPoint3D::new(1.5, 0.0, 0.0)),
    @bio.SAResidue::new("SER", 3, @bio.SAPoint3D::new(3.0, 0.0, 0.0)),
    @bio.SAResidue::new("THR", 4, @bio.SAPoint3D::new(4.5, 0.0, 0.0)),
  ]
  let residues2 = [
    @bio.SAResidue::new("ALA", 1, @bio.SAPoint3D::new(0.0, 0.0, 0.5)),
    @bio.SAResidue::new("GLY", 2, @bio.SAPoint3D::new(1.5, 0.2, 0.5)),
    @bio.SAResidue::new("SER", 3, @bio.SAPoint3D::new(3.0, 0.0, 0.5)),
  ]
  let residues3 = [
    @bio.SAResidue::new("ALA", 1, @bio.SAPoint3D::new(0.0, 1.0, 0.0)),
    @bio.SAResidue::new("GLY", 2, @bio.SAPoint3D::new(1.5, 1.2, 0.0)),
  ]
  
  let struct1 = @bio.SAStructure::new("PDB_1A", residues1)
  let struct2 = @bio.SAStructure::new("PDB_2B", residues2)
  let struct3 = @bio.SAStructure::new("PDB_3C", residues3)
  
  println("  Structure 1 (PDB_1A): " + struct1.residues.length().to_string() + " residues")
  println("  Structure 2 (PDB_2B): " + struct2.residues.length().to_string() + " residues")
  println("  Structure 3 (PDB_3C): " + struct3.residues.length().to_string() + " residues")
  println("")

  println("2. Align two structures")
  let result2 = @bio.bio_structure_align_two(struct1, struct2)
  println("  Aligned structures: " + result2.structure_ids[0] + " vs " + result2.structure_ids[1])
  println("  Aligned length: " + result2.aligned_length.to_string())
  println("  RMSD: " + result2.rmsd.to_string())
  println("  TM-score: " + result2.tm_score.to_string())
  println("")

  println("3. Multiple structure alignment")
  let structures = [struct1, struct2, struct3]
  let result_multi = @bio.bio_structure_align_multiple(structures)
  println("  Number of structures: " + result_multi.structure_ids.length().to_string())
  println("  Aligned length: " + result_multi.aligned_length.to_string())
  println("  Average RMSD: " + result_multi.rmsd.to_string())
  println("")

  println("4. Pairwise RMSD matrix")
  let rmsd_matrix = @bio.bio_structure_pairwise_rmsd(structures)
  println("  RMSD Matrix:")
  println("         " + result_multi.structure_ids.join("     "))
  
  let mut i = 0
  while i < rmsd_matrix.length() {
    let mut row_str = result_multi.structure_ids[i] + "  "
    let mut j = 0
    while j < rmsd_matrix[i].length() {
      row_str = row_str + rmsd_matrix[i][j].to_string() + "  "
      j = j + 1
    }
    println("  " + row_str)
    i = i + 1
  }
  println("")

  println("5. Calculate TM-score directly")
  let tm_score = @bio.bio_structure_calculate_tm_score(struct1, struct2)
  println("  TM-score (PDB_1A vs PDB_2B): " + tm_score.to_string())
  println("")

  println("=== StructureAlignment Demo Complete ===")
}