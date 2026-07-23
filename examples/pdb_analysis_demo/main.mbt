///|
/// PDB Structure Analysis Demo - Advanced protein structure analysis.
///
/// Demonstrates dihedral angles, distance matrices, contact maps,
/// hydrogen bond detection, secondary structure assignment,
/// radius of gyration, and Ramachandran analysis.
///
/// Run: moon run examples/pdb_analysis_demo/main.mbt
fn format_angle_opt(angle : Double?) -> String {
  match angle {
    Some(a) => a.to_string()
    None => "N/A"
  }
}

///|
fn main {
  println("=== PDB Structure Analysis Demo ===\n")

  // === 1. Create demo structure ===
  println("1. Create Demo Protein Structure:")
  println("================================")
  let structure = @src.create_demo_structure()
  println("  Structure ID: " + structure.get_id())
  println("  Models: " + structure.get_num_models().to_string())
  println("  Chains: " + structure.get_num_chains().to_string())
  println("  Residues: " + structure.get_num_residues().to_string())
  println("  Atoms: " + structure.get_num_atoms().to_string())

  // Get the chain for analysis
  let chain_opt = structure.find_model(0).unwrap().find_chain('A')
  let chain = chain_opt.unwrap()
  println("  Chain ID: " + chain.get_id().to_string())
  println("  Chain residues: " + chain.get_num_residues().to_string())

  // === 2. Calculate dihedral angles ===
  println("\n2. Calculate Dihedral Angles (phi, psi, omega):")
  println("================================================")
  let dihedrals = @src.calc_chain_dihedrals(chain)
  let mut i = 0
  while i < dihedrals.length() {
    let (resseq, phi, psi, omega) = dihedrals[i]
    println(
      "  Residue " +
      resseq.to_string() +
      ": phi=" +
      format_angle_opt(phi) +
      ", psi=" +
      format_angle_opt(psi) +
      ", omega=" +
      format_angle_opt(omega),
    )
    i = i + 1
  }

  // === 3. CA distance matrix ===
  println("\n3. CA Distance Matrix:")
  println("======================")
  let dist_matrix = @src.calc_ca_distance_matrix(chain)
  println(
    "  Matrix size: " +
    dist_matrix.length().to_string() +
    "x" +
    dist_matrix.length().to_string(),
  )
  let mut row_i = 0
  while row_i < dist_matrix.length() {
    let row = dist_matrix[row_i]
    let mut line = "  Res " + (row_i + 1).to_string() + ": "
    let mut col_i = 0
    while col_i < row.length() {
      let d = row[col_i]
      let s = d.to_string()
      let len = if s.length() < 6 { s.length() } else { 6 }
      line = line + s[0:len].to_owned() + "\t"
      col_i = col_i + 1
    }
    println(line)
    row_i = row_i + 1
  }

  // === 4. Contact map ===
  println("\n4. Contact Map (threshold=8.0 Å):")
  println("==================================")
  let contact_map = @src.calc_contact_map(chain)
  let mut ci = 0
  while ci < contact_map.length() {
    let row = contact_map[ci]
    let mut line = "  Res " + (ci + 1).to_string() + ": "
    let mut cj = 0
    while cj < row.length() {
      line = line + (if row[cj] { "1" } else { "0" }) + " "
      cj = cj + 1
    }
    println(line)
    ci = ci + 1
  }
  // Count contacts
  let contact_map2 = @src.calc_contact_map(chain, threshold=20.0)
  let mut contact_count = 0
  let mut cci = 0
  while cci < contact_map2.length() {
    let mut ccj = 0
    while ccj < contact_map2[cci].length() {
      if contact_map2[cci][ccj] {
        contact_count = contact_count + 1
      }
      ccj = ccj + 1
    }
    cci = cci + 1
  }
  println("  Total contacts (threshold=20.0 Å): " + contact_count.to_string())

  // === 5. Hydrogen bond detection ===
  println("\n5. Hydrogen Bond Detection:")
  println("============================")
  let hbonds = @src.detect_hydrogen_bonds(chain, max_dist=5.0, min_angle=0.0)
  println("  Number of hydrogen bonds detected: " + hbonds.length().to_string())
  let mut hi = 0
  while hi < hbonds.length() {
    let hb = hbonds[hi]
    println(
      "  [" +
      hi.to_string() +
      "] " +
      hb.donor_resname +
      hb.donor_resseq.to_string() +
      " -> " +
      hb.acceptor_resname +
      hb.acceptor_resseq.to_string() +
      " (dist=" +
      hb.distance.to_string() +
      ", angle=" +
      hb.angle.to_string() +
      ")",
    )
    hi = hi + 1
  }

  // === 6. Secondary structure assignment ===
  println("\n6. Secondary Structure Assignment:")
  println("===================================")
  let ss_assignments = @src.assign_secondary_structure(chain)
  let mut si = 0
  while si < ss_assignments.length() {
    let (resseq, ss) = ss_assignments[si]
    println("  Residue " + resseq.to_string() + ": " + ss.to_code().to_string())
    si = si + 1
  }

  // === 7. Secondary structure statistics ===
  println("\n7. Secondary Structure Statistics:")
  println("==================================")
  let ss_counts = @src.count_secondary_structure(ss_assignments)
  let ss_keys = ["H", "G", "I", "E", "B", "T", "S", "-"]
  let mut ski = 0
  while ski < ss_keys.length() {
    let key = ss_keys[ski]
    match ss_counts.get(key) {
      Some(count) => println("  " + key + ": " + count.to_string())
      None => ()
    }
    ski = ski + 1
  }

  // === 8. Radius of gyration ===
  println("\n8. Radius of Gyration:")
  println("=======================")
  let rg = @src.calc_radius_of_gyration(chain)
  match rg {
    Some(value) => println("  Radius of gyration: " + value.to_string() + " Å")
    None => println("  Radius of gyration: N/A")
  }

  // === 9. Ramachandran plot data ===
  println("\n9. Ramachandran Plot Data:")
  println("===========================")
  let rama_data = @src.ramachandran_plot(chain)
  println("  Residues with phi/psi data: " + rama_data.length().to_string())
  let mut ri = 0
  while ri < rama_data.length() {
    let (resseq, resname, phi, psi) = rama_data[ri]
    println(
      "  Residue " +
      resseq.to_string() +
      " (" +
      resname +
      "): phi=" +
      phi.to_string() +
      ", psi=" +
      psi.to_string(),
    )
    ri = ri + 1
  }

  // === 10. Ramachandran region classification ===
  println("\n10. Ramachandran Region Classification:")
  println("=======================================")
  let mut rc_i = 0
  while rc_i < rama_data.length() {
    let (resseq, resname, phi, psi) = rama_data[rc_i]
    let region = @src.classify_ramachandran(phi, psi)
    println(
      "  Residue " +
      resseq.to_string() +
      " (" +
      resname +
      "): " +
      region.to_string(),
    )
    rc_i = rc_i + 1
  }

  // Demo classification for standard regions
  println("\n  Standard region classifications:")
  let test_cases = [(-60.0, -45.0), (-120.0, 130.0), (60.0, 30.0), (0.0, 0.0)]
  let mut tc_i = 0
  while tc_i < test_cases.length() {
    let (phi, psi) = test_cases[tc_i]
    let region = @src.classify_ramachandran(phi, psi)
    println(
      "  phi=" +
      phi.to_string() +
      ", psi=" +
      psi.to_string() +
      " -> " +
      region.to_string(),
    )
    tc_i = tc_i + 1
  }

  // === 11. Ramachandran quality assessment ===
  println("\n11. Ramachandran Quality Assessment:")
  println("====================================")
  let (favored_pct, allowed_pct, generous_pct) = @src.ramachandran_quality(chain)
  println(
    "  Favored: " + (favored_pct * 100.0).to_string() + "%"
  )
  println(
    "  Allowed: " + (allowed_pct * 100.0).to_string() + "%"
  )
  println(
    "  Generous: " + (generous_pct * 100.0).to_string() + "%"
  )

  // === 12. Solvent Accessible Surface Area (SASA) ===
  println("\n12. Solvent Accessible Surface Area (SASA):")
  println("==========================================")
  let total_sasa = @src.calc_total_sasa(chain)
  println("  Total SASA: " + total_sasa.to_string() + " Å²")
  
  let atom_sasa = @src.calc_sasa(chain, probe_radius=1.4, n_points=50)
  println("  Number of atoms analyzed: " + atom_sasa.length().to_string())
  // Print top 5 atoms by SASA
  println("  Atom SASA values (first 5):")
  let mut ai = 0
  while ai < atom_sasa.length() && ai < 5 {
    let entry = atom_sasa[ai]
    println(
      "    " + entry.0 + ": " + entry.1.to_string() + " Å²"
    )
    ai = ai + 1
  }

  // === 13. Per-residue SASA ===
  println("\n13. Per-Residue SASA:")
  println("=====================")
  let res_sasa = @src.calc_residue_sasa(chain)
  let mut ri_sasa = 0
  while ri_sasa < res_sasa.length() {
    let (resseq, resname, sasa_val) = res_sasa[ri_sasa]
    println(
      "  Residue " +
      resseq.to_string() +
      " (" +
      resname +
      "): " +
      sasa_val.to_string() +
      " Å²"
    )
    ri_sasa = ri_sasa + 1
  }

  // === 14. Hydrophobicity analysis ===
  println("\n14. Hydrophobicity Analysis (Kyte-Doolittle):")
  println("===============================================")
  let hydrophobic_residues = ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO"]
  let hydrophilic_residues = ["ARG", "ASN", "ASP", "GLN", "GLU", "LYS", "SER", "THR", "CYS", "TYR", "HIS", "GLY"]
  
  println("  Hydrophobic residues:")
  let mut hi_i = 0
  while hi_i < hydrophobic_residues.length() {
    let res = hydrophobic_residues[hi_i]
    let hydro = @src.get_hydrophobicity(res)
    println("    " + res + ": " + hydro.to_string())
    hi_i = hi_i + 1
  }
  
  println("  Hydrophilic residues:")
  let mut hyi = 0
  while hyi < hydrophilic_residues.length() {
    let res = hydrophilic_residues[hyi]
    let hydro = @src.get_hydrophobicity(res)
    println("    " + res + ": " + hydro.to_string())
    hyi = hyi + 1
  }

  // === 15. Sequence property matrix ===
  println("\n15. Sequence Property Matrix (Hydrophobicity distance):")
  println("=======================================================")
  let prop_matrix = @src.calc_sequence_property_matrix(chain)
  println("  Matrix size: " + prop_matrix.length().to_string() + "x" + prop_matrix.length().to_string())
  let mut pmi = 0
  while pmi < prop_matrix.length() {
    let row = prop_matrix[pmi]
    let mut line = "  Res " + (pmi + 1).to_string() + ": "
    let mut pmj = 0
    while pmj < row.length() {
      let val = row[pmj]
      let s = val.to_string()
      let len = if s.length() < 6 { s.length() } else { 6 }
      line = line + s[0:len].to_owned() + "  "
      pmj = pmj + 1
    }
    println(line)
    pmi = pmi + 1
  }

  // === 16. Van der Waals radii ===
  println("\n16. Van der Waals Radii:")
  println("========================")
  let atom_types = ["C", "CA", "CB", "N", "NZ", "O", "OG", "S", "SG", "H"]
  let mut vwi = 0
  while vwi < atom_types.length() {
    let atom = atom_types[vwi]
    let radius = @src.vdw_radius(atom)
    println("  " + atom + ": " + radius.to_string() + " Å")
    vwi = vwi + 1
  }

  println("\n=== PDB Structure Analysis Demo Complete ===")
}
