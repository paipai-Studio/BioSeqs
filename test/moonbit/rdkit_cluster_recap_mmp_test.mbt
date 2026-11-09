///|
/// Tests for rdkit_cluster_butina.mbt, rdkit_recap.mbt and rdkit_mmp.mbt.

///|
test "butina_two_pair_cluster" {
  // 4 points: {0,1} tight, {2,3} tight, cross-distance far
  let d : Array[Array[Double]] = [
    [0.0, 0.1, 0.9, 0.9],
    [0.1, 0.0, 0.9, 0.9],
    [0.9, 0.9, 0.0, 0.1],
    [0.9, 0.9, 0.1, 0.0],
  ]
  let clusters = @src.rdkit_butina_cluster(d, 0.3)
  assert_eq(clusters.length(), 2)
  let all : Array[Int32] = []
  for c in clusters {
    for m in c {
      all.push(m)
    }
  }
  all.sort()
  assert_eq(all, [0, 1, 2, 3])
  // members of one cluster are within cutoff of each other
  for c in clusters {
    for a in c {
      for b in c {
        assert_true(d[a][b] <= 0.3)
      }
    }
  }
}

///|
test "butina_all_singletons_when_cutoff_small" {
  let d : Array[Array[Double]] = [
    [0.0, 0.5, 0.5],
    [0.5, 0.0, 0.5],
    [0.5, 0.5, 0.0],
  ]
  let clusters = @src.rdkit_butina_cluster(d, 0.1)
  assert_eq(clusters.length(), 3)
}

///|
test "butina_cluster_smiles_similar_group" {
  // three ethanols + one decane: identical fingerprints must land in
  // the same cluster
  let clusters = @src.rdkit_butina_cluster_smiles(
    ["CCO", "CCO", "CCO", "CCCCCCCCCC"],
    0.4,
  )
  assert_true(clusters.length() >= 1)
  // find the cluster holding the ethanols
  let mut ethanol_cluster = -1
  for ci = 0; ci < clusters.length(); ci = ci + 1 {
    for m in clusters[ci] {
      if m == 0 {
        ethanol_cluster = ci
      }
    }
  }
  for m in clusters[ethanol_cluster] {
    assert_true(m == 0 || m == 1 || m == 2)
  }
}

///|
test "recap_decompose_ester_ethyl_propionate" {
  let frags = @src.rdkit_recap_decompose(@src.parse_smiles("CCC(=O)OCC"))
  assert_true(frags.length() > 0)
  let mut has_ester_ethanol = false
  let mut has_ester_acyl = false
  for f in frags {
    if f.rule == "ester" {
      // canonical form may order the fragment either way
      if f.smiles == "OCC" || f.smiles == "CCO" {
        has_ester_ethanol = true
      }
      if f.smiles == "O=CCC" || f.smiles == "CCC=O" {
        has_ester_acyl = true
      }
    }
  }
  assert_true(has_ester_ethanol)
  assert_true(has_ester_acyl)
}

///|
test "recap_decompose_amine_diethylamine" {
  let frags = @src.rdkit_recap_decompose(@src.parse_smiles("CCNCC"))
  let mut found_ethyl = false
  for f in frags {
    if f.rule == "amine" && f.smiles == "CC" {
      found_ethyl = true
    }
  }
  assert_true(found_ethyl)
}

///|
test "recap_no_fragment_without_rules" {
  // methane: no RECAP rule applies
  let frags = @src.rdkit_recap_decompose(@src.parse_smiles("C"))
  assert_eq(frags.length(), 0)
}

///|
test "mmp_single_cut_propane" {
  let cuts = @src.rdkit_mmp_cuts(@src.parse_smiles("CCC"), max_cuts=1)
  assert_eq(cuts.length(), 2)
  // cutting an end bond: fragment C, context CC
  let mut found_end = false
  for c in cuts {
    if c.fragments.length() == 1 && c.fragments[0] == "C" {
      assert_eq(c.context, "CC")
      found_end = true
    }
  }
  assert_true(found_end)
}

///|
test "mmp_single_cut_butane_middle" {
  let cuts = @src.rdkit_mmp_cuts(@src.parse_smiles("CCCC"), max_cuts=1)
  let mut found_middle = false
  for c in cuts {
    if c.fragments[0] == "CC" {
      assert_eq(c.context, "CC")
      found_middle = true
    }
  }
  assert_true(found_middle)
}

///|
test "mmp_double_cut_butane" {
  let cuts = @src.rdkit_mmp_cuts(@src.parse_smiles("CCCC"), max_cuts=2)
  let mut found_two_frag = false
  for c in cuts {
    if c.fragments.length() == 2 {
      found_two_frag = true
      for f in c.fragments {
        assert_eq(f, "C")
      }
    }
  }
  assert_true(found_two_frag)
}

///|
test "mmp_ignores_ring_bonds" {
  // benzene: no non-ring single bonds -> no cuts
  let cuts = @src.rdkit_mmp_cuts(@src.parse_smiles("C1=CC=CC=C1"), max_cuts=2)
  assert_eq(cuts.length(), 0)
}
