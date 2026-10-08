test "gro_roundtrip" {
  // Build a small GRO file with proper fixed-width columns (44 chars per atom line)
  // 5(resseq) + 5(resname) + 5(atomname) + 5(atomseq) + 8(x) + 8(y) + 8(z)
  let line1 = "    1ALA  CA       1   0.500   0.500   0.500"
  let line2 = "    1ALA  N        2   0.400   0.500   0.500"
  let gro_str = "Test GRO\n2\n" + line1 + "\n" + line2 + "\n1.00000 1.00000 1.00000\n"
  let parsed = @src.parse_gro(gro_str)
  assert_true(parsed is Some(_))
  match parsed {
    Some(g) => {
      assert_eq(g.atoms.length(), 2)
      assert_eq(g.atoms[0].resname, "ALA")
      assert_eq(g.atoms[0].atom_name, "CA")
      assert_true(g.atoms[0].x > 0.49 && g.atoms[0].x < 0.51)
      assert_true(g.atoms[1].atom_name == "N")
      // Box
      assert_true(g.box.0 > 0.99 && g.box.0 < 1.01)
      // Write back
      let written = @src.write_gro(g)
      let reparsed = @src.parse_gro(written)
      assert_true(reparsed is Some(_))
      match reparsed {
        Some(g2) => {
          assert_eq(g2.atoms.length(), 2)
          assert_eq(g2.atoms[0].resname, "ALA")
        }
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
}

test "gro_to_pdb_conversion" {
  // 0.5 nm = 5 Å
  let atoms : Array[@src.GroAtom] = [
    @src.GroAtom::{
      resseq: 1,
      resname: "GLY",
      atom_name: "CA",
      atom_seq: 1,
      x: 0.5,
      y: 0.0,
      z: 0.0,
      vx: None,
      vy: None,
      vz: None,
    },
  ]
  let gro = @src.GroFile::{ title: "test", atoms, box: (1.0, 1.0, 1.0) }
  let pdb_atoms = @src.gro_to_pdb_atoms(gro)
  assert_eq(pdb_atoms.length(), 1)
  let coord = pdb_atoms[0].get_coord()
  assert_true(coord.x > 4.9 && coord.x < 5.1)
}
