///|
/// MMTF (Macromolecular Transmission Format) demo.
///
/// Demonstrates loading an MMTF structure, inspecting the hierarchy
/// (models/chains/groups/atoms), run-length and delta encoding,
/// text serialization round-trip, and conversion to a PDB Structure.
///
/// Run: moon run examples/mmtf_demo --target wasm-gc
fn main {
  // Load the sample MMTF structure
  let structure = @bio.mmtf_sample()

  println("=== MMTF Structure Demo ===")
  println("")

  // === Structure info ===
  println("1. Structure Info")
  let header = structure.get_header()
  let xtal = structure.get_xtal_info()
  println("   ID: \{structure.structure_id}")
  println("   Title: \{header.title}")
  println("   Resolution: \{header.resolution} Å")
  println("   R-Work: \{header.r_work}")
  println("   R-Free: \{header.r_free}")
  println("   Deposition Date: \{header.deposition_date}")
  println("   Release Date: \{header.release_date}")
  println("   Experimental Methods: \{header.experimental_methods.length()}")
  for m in header.experimental_methods {
    println("     - \{m}")
  }
  println("   Space Group: \{xtal.space_group}")
  println(
    "   Unit Cell: a=\{xtal.unit_cell[0]}, b=\{xtal.unit_cell[1]}, c=\{xtal.unit_cell[2]}",
  )
  println("")

  // === Models, chains, groups, atoms ===
  println("2. Hierarchy")
  println("   Models: \{structure.get_num_models()}")
  println("   Chains (all models): \{structure.count_chains()}")
  println("   Groups (all models): \{structure.count_groups()}")
  println("   Atoms (all models): \{structure.count_atoms()}")
  println("   Bonds: \{structure.bonds.length()}")
  println("   Entities: \{structure.entities.length()}")
  println("")

  // === Detail per model/chain/group ===
  println("3. Detail")
  let models = structure.get_models()
  for mi = 0; mi < models.length(); mi = mi + 1 {
    let model = models[mi]
    println("   Model \{model.model_id}:")
    let chains = model.get_chains()
    for ci = 0; ci < chains.length(); ci = ci + 1 {
      let chain = chains[ci]
      let seq = @bio.mmtf_chain_sequence(chain)
      println(
        "     Chain \{chain.chain_id} (name=\{chain.chain_name}): \{chain.get_num_groups()} groups, sequence='\{seq}'",
      )
      let groups = chain.get_groups()
      for gi = 0; gi < groups.length(); gi = gi + 1 {
        let group = groups[gi]
        println(
          "       Group \{group.group_name} \{group.group_number} (\{group.single_letter_code}) ss=\{group.sec_struct.to_string()} atoms=\{group.get_num_atoms()}",
        )
        let atoms = group.get_atoms()
        for ai = 0; ai < atoms.length(); ai = ai + 1 {
          let atom = atoms[ai]
          println(
            "         Atom \{atom.name} serial=\{atom.serial} (\{atom.x}, \{atom.y}, \{atom.z}) element=\{atom.element}",
          )
        }
      }
    }
  }
  println("")

  // === Run-length encoding/decoding ===
  println("4. Run-Length Encoding")
  let chain_id_list = [0, 0, 0, 1, 1, 2, 2, 2, 2]
  let encoded = @bio.mmtf_encode_run_length(chain_id_list)
  println("   Original: \{chain_id_list.length()} values")
  println("   Encoded: \{encoded.length()} values (pairs of value,count)")
  let decoded = @bio.mmtf_decode_run_length(encoded)
  println("   Decoded: \{decoded.length()} values")
  // Verify round-trip
  let mut rle_ok = true
  if decoded.length() != chain_id_list.length() {
    rle_ok = false
  } else {
    for i = 0; i < chain_id_list.length(); i = i + 1 {
      if decoded[i] != chain_id_list[i] {
        rle_ok = false
      }
    }
  }
  println("   Round-trip OK: \{rle_ok}")
  println("")

  // === Delta encoding/decoding ===
  println("5. Delta Encoding")
  let coords = [100, 101, 103, 106, 110, 115]
  let delta_encoded = @bio.mmtf_encode_delta(coords)
  println("   Original: \{coords}")
  println("   Delta-encoded: \{delta_encoded}")
  let delta_decoded = @bio.mmtf_decode_delta(delta_encoded)
  println("   Decoded: \{delta_decoded}")
  println("")

  // === Text serialization round-trip ===
  println("6. Text Serialization")
  let text = @bio.mmtf_write_text(structure)
  println("   Written \{text.length()} chars of text")
  let parsed = @bio.mmtf_read_text(text)
  println("   Parsed back: ID=\{parsed.structure_id}")
  println("   Parsed models: \{parsed.get_num_models()}")
  println("   Parsed chains: \{parsed.count_chains()}")
  println("   Parsed groups: \{parsed.count_groups()}")
  println("   Parsed atoms: \{parsed.count_atoms()}")
  println("   Parsed bonds: \{parsed.bonds.length()}")
  println("")

  // === Convert to PDB Structure ===
  println("7. Convert to PDB Structure")
  let pdb = @bio.mmtf_to_pdb_structure(structure)
  println("   PDB ID: \{pdb.get_id()}")
  println("   PDB Models: \{pdb.get_num_models()}")
  println("   PDB Chains: \{pdb.get_num_chains()}")
  println("   PDB Residues: \{pdb.get_num_residues()}")
  println("   PDB Atoms: \{pdb.get_num_atoms()}")
  // Show first few residues
  let residues = pdb.get_residues()
  let limit = if residues.length() < 3 { residues.length() } else { 3 }
  for i = 0; i < limit; i = i + 1 {
    let r = residues[i]
    println(
      "   Residue[\{i}]: \{r.resname} \{r.chainid}\{r.resseq} atoms=\{r.get_num_atoms()} het=\{r.is_het()}",
    )
  }
  println("")

  // === Summary ===
  println("8. Summary")
  let summary = @bio.mmtf_summary(structure)
  println(summary)

  println("=== Done ===")
}
