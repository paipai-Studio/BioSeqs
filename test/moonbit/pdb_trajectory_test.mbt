///|
/// Tests for PDB trajectory I/O — multi-model PDB as trajectory.

///|
fn make_multi_model_pdb() -> String {
  let buf = StringBuilder()
  // Model 1
  buf.write_string("MODEL     1\n")
  buf.write_string("ATOM      1  N   ALA A   1       1.000   0.000   0.000  1.00 20.00           N\n")
  buf.write_string("ATOM      2  CA  ALA A   1       2.000   0.000   0.000  1.00 20.00           C\n")
  buf.write_string("ATOM      3  N   GLY A   2       3.000   0.000   0.000  1.00 15.00           N\n")
  buf.write_string("ATOM      4  CA  GLY A   2       4.000   0.000   0.000  1.00 15.00           C\n")
  buf.write_string("ENDMDL\n")
  // Model 2 (shifted by 1 in x)
  buf.write_string("MODEL     2\n")
  buf.write_string("ATOM      1  N   ALA A   1       2.000   0.000   0.000  1.00 20.00           N\n")
  buf.write_string("ATOM      2  CA  ALA A   1       3.000   0.000   0.000  1.00 20.00           C\n")
  buf.write_string("ATOM      3  N   GLY A   2       4.000   0.000   0.000  1.00 15.00           N\n")
  buf.write_string("ATOM      4  CA  GLY A   2       5.000   0.000   0.000  1.00 15.00           C\n")
  buf.write_string("ENDMDL\n")
  // Model 3 (shifted by 2 in x)
  buf.write_string("MODEL     3\n")
  buf.write_string("ATOM      1  N   ALA A   1       3.000   0.000   0.000  1.00 20.00           N\n")
  buf.write_string("ATOM      2  CA  ALA A   1       4.000   0.000   0.000  1.00 20.00           C\n")
  buf.write_string("ATOM      3  N   GLY A   2       5.000   0.000   0.000  1.00 15.00           N\n")
  buf.write_string("ATOM      4  CA  GLY A   2       6.000   0.000   0.000  1.00 15.00           C\n")
  buf.write_string("ENDMDL\n")
  buf.write_string("END\n")
  buf.to_string()
}

test "parse_pdb_trajectory basic" {
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      assert_eq(t.num_frames(), 3)
      assert_eq(t.num_atoms(), 4)
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory get_frame" {
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      let frame0 = t.get_frame(0)
      assert_eq(frame0.length(), 4)
      assert_true((frame0[0].x - 1.0).abs() < 1.0e-6)
      let frame1 = t.get_frame(1)
      assert_true((frame1[0].x - 2.0).abs() < 1.0e-6)
      let frame2 = t.get_frame(2)
      assert_true((frame2[0].x - 3.0).abs() < 1.0e-6)
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory get_atom_coords" {
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      // Atom 0 across frames: x = 1, 2, 3
      let coords = t.get_atom_coords(0)
      assert_eq(coords.length(), 3)
      assert_true((coords[0].x - 1.0).abs() < 1.0e-6)
      assert_true((coords[1].x - 2.0).abs() < 1.0e-6)
      assert_true((coords[2].x - 3.0).abs() < 1.0e-6)
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory get_coords" {
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      let coords = t.get_coords()
      assert_eq(coords.length(), 3)
      assert_eq(coords[0].length(), 4)
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory get_coords_flat" {
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      let flat = t.get_coords_flat()
      assert_eq(flat.length(), 3)
      assert_eq(flat[0].length(), 4)
      assert_eq(flat[0][0].length(), 3)
      assert_true((flat[0][0][0] - 1.0).abs() < 1.0e-6)
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory slice" {
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      let sub = t.slice(1, 3, 1)
      assert_eq(sub.num_frames(), 2)
      let frame0 = sub.get_frame(0)
      assert_true((frame0[0].x - 2.0).abs() < 1.0e-6)
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory slice with step" {
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      let sub = t.slice(0, 3, 2)
      assert_eq(sub.num_frames(), 2) // frames 0 and 2
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory mean_coords" {
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      let mean = t.mean_coords()
      assert_eq(mean.length(), 4)
      // Atom 0: mean of (1,2,3) = 2.0
      assert_true((mean[0].x - 2.0).abs() < 1.0e-6)
      // Atom 1: mean of (2,3,4) = 3.0
      assert_true((mean[1].x - 3.0).abs() < 1.0e-6)
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory rmsd_matrix" {
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      let mat = t.rmsd_matrix()
      assert_eq(mat.length(), 3)
      // Diagonal should be 0
      assert_true(mat[0][0].abs() < 1.0e-10)
      assert_true(mat[1][1].abs() < 1.0e-10)
      assert_true(mat[2][2].abs() < 1.0e-10)
      // Symmetric
      assert_true((mat[0][1] - mat[1][0]).abs() < 1.0e-10)
      // Frame 0 vs 1: all atoms shifted by 1 in x → RMSD = 1.0
      assert_true((mat[0][1] - 1.0).abs() < 1.0e-6)
      // Frame 0 vs 2: shifted by 2 → RMSD = 2.0
      assert_true((mat[0][2] - 2.0).abs() < 1.0e-6)
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory rmsf" {
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      let r = t.rmsf()
      assert_eq(r.length(), 4)
      // All atoms have same displacement pattern → same RMSF
      // Atom 0: values 1,2,3, mean=2, deviations -1,0,1, RMSF = sqrt(2/3)
      let expected = (2.0 / 3.0).sqrt()
      assert_true((r[0] - expected).abs() < 1.0e-6)
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory deltag" {
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      let d = t.deltag()
      assert_eq(d.length(), 3)
      assert_eq(d[0].length(), 4)
      // Frame 0, atom 0: 1.0 - 2.0 = -1.0
      assert_true((d[0][0].x - (-1.0)).abs() < 1.0e-6)
      // Frame 1, atom 0: 2.0 - 2.0 = 0.0
      assert_true(d[1][0].x.abs() < 1.0e-6)
      // Frame 2, atom 0: 3.0 - 2.0 = 1.0
      assert_true((d[2][0].x - 1.0).abs() < 1.0e-6)
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory superpose_frames identical" {
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      // Superpose to first frame — frame 0 should be unchanged
      let aligned = t.superpose_frames()
      let frame0 = aligned.get_frame(0)
      // First frame is the reference, should be ~identical
      let orig = t.get_frame(0)
      for i = 0; i < 4; i = i + 1 {
        assert_true((frame0[i].x - orig[i].x).abs() < 1.0e-4)
      }
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory superpose_frames translation" {
  // Pure translation should be fully removed by superposition
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      let aligned = t.superpose_frames()
      // After superposition, all frames should have ~same centroid
      let c0 = aligned.get_frame(0)
      let c1 = aligned.get_frame(1)
      let mut centroid0_x = 0.0
      let mut centroid1_x = 0.0
      for i = 0; i < 4; i = i + 1 {
        centroid0_x = centroid0_x + c0[i].x
        centroid1_x = centroid1_x + c1[i].x
      }
      centroid0_x = centroid0_x / 4.0
      centroid1_x = centroid1_x / 4.0
      assert_true((centroid0_x - centroid1_x).abs() < 1.0e-4)
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory write_pdb" {
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      let written = t.write_pdb()
      // Should contain MODEL and ENDMDL
      assert_true(written.contains("MODEL"))
      assert_true(written.contains("ENDMDL"))
      assert_true(written.contains("END"))
      // Re-parse should give same trajectory
      let reparsed = @src.parse_pdb_trajectory(written)
      match reparsed {
        Some(t2) => {
          assert_eq(t2.num_frames(), 3)
          assert_eq(t2.num_atoms(), 4)
          // Coordinates should match
          for f = 0; f < 3; f = f + 1 {
            let f1 = t.get_frame(f)
            let f2 = t2.get_frame(f)
            for a = 0; a < 4; a = a + 1 {
              assert_true((f1[a].x - f2[a].x).abs() < 0.01)
            }
          }
        }
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory to_structure" {
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      let s = t.to_structure()
      assert_eq(s.get_num_models(), 3)
      assert_eq(s.get_models()[0].get_atoms().length(), 4)
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory from_structure" {
  let pdb = make_multi_model_pdb()
  let structure = @src.parse_pdb(pdb)
  let traj = @src.PDBTrajectory::from_structure(structure)
  match traj {
    Some(t) => {
      assert_eq(t.num_frames(), 3)
      assert_eq(t.num_atoms(), 4)
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory from_structure empty" {
  let s = @src.Structure::new(id="empty")
  let traj = @src.PDBTrajectory::from_structure(s)
  match traj {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

test "PDBTrajectory single model" {
  let pdb = "MODEL     1\nATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nENDMDL\nEND\n"
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      assert_eq(t.num_frames(), 1)
      assert_eq(t.num_atoms(), 1)
      let frame = t.get_frame(0)
      assert_true((frame[0].x - 1.0).abs() < 1.0e-6)
      assert_true((frame[0].y - 2.0).abs() < 1.0e-6)
      assert_true((frame[0].z - 3.0).abs() < 1.0e-6)
      // RMSF of single frame = 0
      let r = t.rmsf()
      assert_true(r[0].abs() < 1.0e-10)
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory round_trip" {
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      let s = t.to_structure()
      let t2 = @src.PDBTrajectory::from_structure(s)
      match t2 {
        Some(t2v) => {
          assert_eq(t2v.num_frames(), t.num_frames())
          assert_eq(t2v.num_atoms(), t.num_atoms())
          for f = 0; f < 3; f = f + 1 {
            let f1 = t.get_frame(f)
            let f2 = t2v.get_frame(f)
            for a = 0; a < 4; a = a + 1 {
              assert_true((f1[a].x - f2[a].x).abs() < 1.0e-6)
              assert_true((f1[a].y - f2[a].y).abs() < 1.0e-6)
              assert_true((f1[a].z - f2[a].z).abs() < 1.0e-6)
            }
          }
        }
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
}

test "PDBTrajectory rmsd_matrix symmetric" {
  let pdb = make_multi_model_pdb()
  let traj = @src.parse_pdb_trajectory(pdb)
  match traj {
    Some(t) => {
      let mat = t.rmsd_matrix()
      for i = 0; i < 3; i = i + 1 {
        for j = 0; j < 3; j = j + 1 {
          assert_true((mat[i][j] - mat[j][i]).abs() < 1.0e-10)
        }
      }
    }
    None => assert_true(false)
  }
}
