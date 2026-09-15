///|
/// Tests for abstract_property_map module (Bio.PDB.AbstractPropertyMap).

test "abstract_property_map basic add and get" {
  let pm = @src.AbstractPropertyMap::new(["SS", "ASA"])
  let props : Map[String, String] = {
    let m : Map[String, String] = Map([], capacity=2)
    m["SS"] = "H"
    m["ASA"] = "75.5"
    m
  }
  pm.add("A", 1, props)
  assert_eq(pm.n_residues, 1)
  let ss = pm.get_property("A", 1, "SS")
  assert_true(ss is Some(_))
  assert_eq(ss.unwrap(), "H")
}

test "abstract_property_map has_residue" {
  let pm = @src.AbstractPropertyMap::new(["SS"])
  let props : Map[String, String] = {
    let m : Map[String, String] = Map([], capacity=1)
    m["SS"] = "H"
    m
  }
  pm.add("A", 1, props)
  assert_true(pm.has_residue("A", 1))
  assert_false(pm.has_residue("A", 2))
  assert_false(pm.has_residue("B", 1))
}

test "abstract_property_map icode" {
  let pm = @src.AbstractPropertyMap::new(["SS"])
  let props : Map[String, String] = {
    let m : Map[String, String] = Map([], capacity=1)
    m["SS"] = "H"
    m
  }
  pm.add("A", 1, props, icode="A")
  assert_true(pm.has_residue("A", 1, icode="A"))
  assert_false(pm.has_residue("A", 1))
}

test "abstract_property_map keys" {
  let pm = @src.abstract_property_map_example()
  let keys = pm.keys()
  assert_eq(keys.length(), 3)
}

test "abstract_property_map all_records" {
  let pm = @src.abstract_property_map_example()
  let records = pm.all_records()
  assert_eq(records.length(), 3)
}

test "abstract_property_map filter_by_property" {
  let pm = @src.abstract_property_map_example()
  // Filter residues where SS == "H"
  let filtered = pm.filter_by_property("SS", fn(s) { s == "H" })
  assert_eq(filtered.n_residues, 1)
}

test "abstract_property_map stats_property" {
  let pm = @src.abstract_property_map_example()
  let (count, min_v, max_v, mean, sum) = pm.stats_property("ASA")
  assert_eq(count, 3)
  assert_true((min_v - 45.3).abs() < 1.0e-9)
  assert_true((max_v - 120.0).abs() < 1.0e-9)
  assert_true((sum - 240.8).abs() < 1.0e-9)
  assert_true((mean - 80.26666666666667).abs() < 1.0e-6)
}

test "abstract_property_map to_tsv" {
  let pm = @src.abstract_property_map_example()
  let tsv = pm.to_tsv()
  assert_true(tsv.contains("chain_id"))
  assert_true(tsv.contains("resseq"))
  assert_true(tsv.contains("SS"))
  assert_true(tsv.contains("ASA"))
  assert_true(tsv.contains("A\t1"))
}

test "abstract_property_map merge" {
  let pm1 = @src.AbstractPropertyMap::new(["SS"])
  let props1 : Map[String, String] = {
    let m : Map[String, String] = Map([], capacity=1)
    m["SS"] = "H"
    m
  }
  pm1.add("A", 1, props1)
  let pm2 = @src.AbstractPropertyMap::new(["ASA"])
  let props2 : Map[String, String] = {
    let m : Map[String, String] = Map([], capacity=1)
    m["ASA"] = "80.0"
    m
  }
  pm2.add("A", 1, props2)
  let merged = @src.AbstractPropertyMap::merge(pm1, pm2)
  assert_true(merged.has_residue("A", 1))
  let ss = merged.get_property("A", 1, "SS")
  let asa = merged.get_property("A", 1, "ASA")
  assert_eq(ss.unwrap(), "H")
  assert_eq(asa.unwrap(), "80.0")
}

test "abstract_property_map get_residue" {
  let pm = @src.abstract_property_map_example()
  let props = pm.get_residue("A", 1)
  assert_true(props is Some(_))
  let p = props.unwrap()
  assert_eq(p["SS"], "H")
  assert_eq(p["ASA"], "75.5")
}

test "abstract_property_map example" {
  let pm = @src.abstract_property_map_example()
  assert_eq(pm.n_residues, 3)
  assert_eq(pm.property_names.length(), 3)
}

test "residue_key_to_string" {
  let key = @src.ResidueKey::new("A", 42, icode="B")
  let s = @src.residue_key_to_string(key)
  assert_eq(s, "A:42:B")
}

test "abstract_property_map missing property" {
  let pm = @src.abstract_property_map_example()
  let missing = pm.get_property("A", 1, "NonExistent")
  assert_true(missing is None)
}