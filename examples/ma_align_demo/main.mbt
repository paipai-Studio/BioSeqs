///|
/// Multiple Protein Structure Alignment Demo
/// Demonstrates Bio.PDB.MAalign functionality
///
/// Run: moon run examples/ma_align_demo/

fn main {
  println("=== Multiple Structure Alignment Demo ===")

  // Create aligner
  let aligner = @bio.ma_aligner_new()
  println("Aligner parameters:")
  println("  Gap opening: \{aligner.get_gap_opening().to_string()}")
  println("  Gap extension: \{aligner.get_gap_extension().to_string()}")
  println("  Max iterations: \{aligner.get_max_iterations().to_string()}")
  println("  RMSD threshold: \{aligner.get_rmsd_threshold().to_string()} A")

  // Create mock protein structures (CA atoms only)
  // Structure 1: Simple helix-like arrangement
  let coords1 : Array[(Double, Double, Double)] = [
    (0.0, 0.0, 0.0),
    (1.5, 0.5, 0.2),
    (2.8, 1.2, 0.5),
    (3.9, 2.1, 1.0),
    (4.8, 3.2, 1.8),
    (5.5, 4.5, 2.8),
    (6.0, 6.0, 4.0),
    (6.2, 7.7, 5.4),
    (6.1, 9.5, 6.8),
    (5.8, 11.3, 8.2),
  ]
  let seq1 = "ALAGLYALA"
  let residues1 = ["ALA", "GLY", "ALA", "VAL", "LEU", "ILE", "PHE", "MET", "TRP", "TYR"]

  // Structure 2: Same shape but translated and rotated
  let coords2 : Array[(Double, Double, Double)] = [
    (10.0, 10.0, 5.0),
    (11.5, 10.5, 5.2),
    (12.8, 11.2, 5.5),
    (13.9, 12.1, 6.0),
    (14.8, 13.2, 6.8),
    (15.5, 14.5, 7.8),
    (16.0, 16.0, 9.0),
    (16.2, 17.7, 10.4),
    (16.1, 19.5, 11.8),
    (15.8, 21.3, 13.2),
  ]
  let seq2 = "ALAGLYALA"
  let residues2 = ["ALA", "GLY", "ALA", "VAL", "LEU", "ILE", "PHE", "MET", "TRP", "TYR"]

  // Structure 3: Different conformation
  let coords3 : Array[(Double, Double, Double)] = [
    (0.0, 0.0, 0.0),
    (1.5, 0.0, 0.0),
    (3.0, 0.0, 0.0),
    (4.5, 0.0, 0.0),
    (6.0, 0.0, 0.0),
    (7.5, 0.0, 0.0),
    (9.0, 0.0, 0.0),
    (10.5, 0.0, 0.0),
    (12.0, 0.0, 0.0),
    (13.5, 0.0, 0.0),
  ]
  let seq3 = "ALAGLYALA"
  let residues3 = ["ALA", "GLY", "ALA", "VAL", "LEU", "ILE", "PHE", "MET", "TRP", "TYR"]

  // Create structures
  println("")
  println("Creating structures...")
  let s1 = @bio.add_structure(aligner, "helix1", coords1, seq1, residues1)
  let s2 = @bio.add_structure(aligner, "helix2", coords2, seq2, residues2)
  let s3 = @bio.add_structure(aligner, "extended", coords3, seq3, residues3)

  println("Structure 1: \{s1.get_id()}, \{s1.length().to_string()} CA atoms")
  println("Structure 2: \{s2.get_id()}, \{s2.length().to_string()} CA atoms")
  println("Structure 3: \{s3.get_id()}, \{s3.length().to_string()} CA atoms")

  // Pairwise alignment demo
  println("")
  println("=== Pairwise Alignment ===")
  let (rot1, rmsd1) = @bio.ma_align_pairwise(coords1, coords2)
  println("Alignment helix1 vs helix2:")
  println("  RMSD: \{rmsd1.to_string()} A")
  println("  Rotation matrix (Euler angles):")
  println("    x-axis: \{rot1.0.to_string()}, y-axis: \{rot1.1.to_string()}, z-axis: \{rot1.2.to_string()}")

  let (rot2, rmsd2) = @bio.ma_align_pairwise(coords1, coords3)
  println("Alignment helix1 vs extended:")
  println("  RMSD: \{rmsd2.to_string()} A")
  println("  Rotation matrix (Euler angles):")
  println("    x-axis: \{rot2.0.to_string()}, y-axis: \{rot2.1.to_string()}, z-axis: \{rot2.2.to_string()}")
  println("  Higher RMSD = more different structures")

  // Multiple structure alignment
  println("")
  println("=== Multiple Structure Alignment ===")
  let structures = [s1, s2, s3]
  let result = @bio.align_structures(aligner, structures)

  println("Alignment result:")
  println("  Structures aligned: \{result.get_structure_ids().length().to_string()}")
  println("  Total RMSD: \{result.get_total_rmsd().to_string()} A")

  // Display aligned sequences
  println("")
  println("Aligned sequences:")
  let aligned_seqs = result.get_aligned_sequences()
  let struct_ids = result.get_structure_ids()
  for i = 0; i < struct_ids.length(); i = i + 1 {
    println("  \{struct_ids[i]}: \{aligned_seqs[i]}")
  }

  // Conservation profile
  println("")
  println("Conservation profile:")
  let conservation = @bio.get_conservation(result)
  for i = 0; i < conservation.length(); i = i + 1 {
    let score = conservation[i]
    let bar_count = (score * 10.0).to_int()
    let mut bars = ""
    for j = 0; j < bar_count; j = j + 1 {
      bars = bars + "|"
    }
    println("  Col \{(i + 1).to_string()}: \{score.to_string()} \{bars}")
  }

  // RMSD matrix
  println("")
  println("Pairwise RMSD matrix (A):")
  let rmsd_mat = @bio.get_rmsd_matrix(result)
  let mut header = "      "
  for id in struct_ids {
    header = header + id + "  "
  }
  println(header)
  for i = 0; i < rmsd_mat.length(); i = i + 1 {
    let mut row = struct_ids[i] + "  "
    for j = 0; j < rmsd_mat[i].length(); j = j + 1 {
      row = row + rmsd_mat[i][j].to_string() + "  "
    }
    println(row)
  }

  // Aligned coordinates
  println("")
  println("Aligned CA coordinates (structure 1):")
  let aligned_coords = result.get_aligned_coordinates()
  let coords_set1 = aligned_coords[0]
  for i = 0; i < coords_set1.length(); i = i + 1 {
    let (x, y, z) = coords_set1[i]
    println("  Atom \{(i + 1).to_string()}: (\{x.to_string()}, \{y.to_string()}, \{z.to_string()})")
  }

  println("=== Demo complete ===")
}