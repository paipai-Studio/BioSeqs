///|
/// Tests for Bio.SCOP structural classification of proteins parser.

// ---------------------------------------------------------------------------
// Node type helpers
// ---------------------------------------------------------------------------

test "scop_node_type_name" {
  assert_eq(@src.scop_node_type_name("ro"), "root")
  assert_eq(@src.scop_node_type_name("cl"), "class")
  assert_eq(@src.scop_node_type_name("cf"), "fold")
  assert_eq(@src.scop_node_type_name("sf"), "superfamily")
  assert_eq(@src.scop_node_type_name("fa"), "family")
  assert_eq(@src.scop_node_type_name("dm"), "protein")
  assert_eq(@src.scop_node_type_name("sp"), "species")
  assert_eq(@src.scop_node_type_name("px"), "domain")
  assert_eq(@src.scop_node_type_name("xx"), "unknown")
}

test "scop_node_type_code" {
  assert_eq(@src.scop_node_type_code("root"), "ro")
  assert_eq(@src.scop_node_type_code("class"), "cl")
  assert_eq(@src.scop_node_type_code("fold"), "cf")
  assert_eq(@src.scop_node_type_code("superfamily"), "sf")
  assert_eq(@src.scop_node_type_code("family"), "fa")
  assert_eq(@src.scop_node_type_code("protein"), "dm")
  assert_eq(@src.scop_node_type_code("species"), "sp")
  assert_eq(@src.scop_node_type_code("domain"), "px")
}

test "scop_node_type_order" {
  assert_true(@src.scop_node_type_order("ro") < @src.scop_node_type_order("cl"))
  assert_true(@src.scop_node_type_order("cl") < @src.scop_node_type_order("cf"))
  assert_true(@src.scop_node_type_order("cf") < @src.scop_node_type_order("sf"))
  assert_true(@src.scop_node_type_order("sf") < @src.scop_node_type_order("fa"))
  assert_true(@src.scop_node_type_order("fa") < @src.scop_node_type_order("dm"))
  assert_true(@src.scop_node_type_order("dm") < @src.scop_node_type_order("sp"))
  assert_true(@src.scop_node_type_order("sp") < @src.scop_node_type_order("px"))
}

// ---------------------------------------------------------------------------
// Residues parsing
// ---------------------------------------------------------------------------

test "scop_parse_residues_dash" {
  let r = @src.scop_parse_residues("-")
  assert_eq(r.pdbid(), "")
  assert_eq(r.fragments().length(), 0)
}

test "scop_parse_residues_empty" {
  let r = @src.scop_parse_residues("")
  assert_eq(r.pdbid(), "")
  assert_eq(r.fragments().length(), 0)
}

test "scop_parse_residues_paren_dash" {
  let r = @src.scop_parse_residues("(-)")
  assert_eq(r.pdbid(), "")
  assert_eq(r.fragments().length(), 1)
  assert_eq(r.fragments()[0].chain(), "")
  assert_eq(r.fragments()[0].start(), "")
  assert_eq(r.fragments()[0].end_(), "")
}

test "scop_parse_residues_with_pdbid" {
  let r = @src.scop_parse_residues("1bba A:10-20,B:")
  assert_eq(r.pdbid(), "1bba")
  assert_eq(r.fragments().length(), 2)
  // First fragment: A:10-20
  assert_eq(r.fragments()[0].chain(), "A")
  assert_eq(r.fragments()[0].start(), "10")
  assert_eq(r.fragments()[0].end_(), "20")
  // Second fragment: B: (whole chain)
  assert_eq(r.fragments()[1].chain(), "B")
  assert_eq(r.fragments()[1].start(), "")
  assert_eq(r.fragments()[1].end_(), "")
}

test "scop_parse_residues_single_chain" {
  let r = @src.scop_parse_residues("A:1-141")
  assert_eq(r.pdbid(), "")
  assert_eq(r.fragments().length(), 1)
  assert_eq(r.fragments()[0].chain(), "A")
  assert_eq(r.fragments()[0].start(), "1")
  assert_eq(r.fragments()[0].end_(), "141")
}

test "scop_residues_to_string_roundtrip" {
  let r = @src.scop_parse_residues("1hba A:1-141")
  let s = @src.scop_residues_to_string(r)
  let r2 = @src.scop_parse_residues(s)
  assert_eq(r2.pdbid(), "1hba")
  assert_eq(r2.fragments().length(), 1)
  assert_eq(r2.fragments()[0].chain(), "A")
  assert_eq(r2.fragments()[0].start(), "1")
  assert_eq(r2.fragments()[0].end_(), "141")
}

test "scop_residues_to_string_dash" {
  let r = @src.ScopResidues::new("", [])
  let s = @src.scop_residues_to_string(r)
  assert_eq(s, "-")
}

// ---------------------------------------------------------------------------
// Record creation and accessors
// ---------------------------------------------------------------------------

test "scop_cla_record_creation" {
  let res = @src.scop_parse_residues("1hba A:1-141")
  let hier : Map[String, Int] = Map([], capacity=8)
  hier.set("cl", 100)
  hier.set("cf", 200)
  let r = @src.ClaRecord::new("d1hba_", "1hba", res, "a.1.1.1.1.1.1", 1000, hier)
  assert_eq(r.sid(), "d1hba_")
  assert_eq(r.pdbid(), "1hba")
  assert_eq(r.sccs(), "a.1.1.1.1.1.1")
  assert_eq(r.sunid(), 1000)
  assert_eq(r.hierarchy().get("cl").unwrap(), 100)
}

test "scop_des_record_creation" {
  let r = @src.DesRecord::new(
    1000, "px", "a.1.1.1.1.1.1", "d1hba_", "1hba Hemoglobin alpha chain",
  )
  assert_eq(r.sunid(), 1000)
  assert_eq(r.nodetype(), "px")
  assert_eq(r.sccs(), "a.1.1.1.1.1.1")
  assert_eq(r.name(), "d1hba_")
  assert_eq(r.description(), "1hba Hemoglobin alpha chain")
}

test "scop_hie_record_creation" {
  let r = @src.HieRecord::new(600, 500, [1000, 1001])
  assert_eq(r.sunid(), 600)
  assert_eq(r.parent(), 500)
  assert_eq(r.children().length(), 2)
  assert_eq(r.children()[0], 1000)
  assert_eq(r.children()[1], 1001)
}

// ---------------------------------------------------------------------------
// File parsers
// ---------------------------------------------------------------------------

test "scop_parse_cla_line" {
  let line = "d1hba_\t1hba\tA:1-141\ta.1.1.1.1.1.1\t1000\tcl=100,cf=200,sf=300,fa=400,dm=500,sp=600,px=1000"
  let r = @src.scop_parse_cla_line(line)
  assert_true(r.is_some())
  let rec = r.unwrap()
  assert_eq(rec.sid(), "d1hba_")
  assert_eq(rec.pdbid(), "1hba")
  assert_eq(rec.sccs(), "a.1.1.1.1.1.1")
  assert_eq(rec.sunid(), 1000)
  assert_eq(rec.hierarchy().get("cl").unwrap(), 100)
  assert_eq(rec.hierarchy().get("px").unwrap(), 1000)
}

test "scop_parse_cla_line_comment" {
  let r = @src.scop_parse_cla_line("# comment line")
  assert_true(r.is_none())
}

test "scop_parse_cla_line_empty" {
  let r = @src.scop_parse_cla_line("")
  assert_true(r.is_none())
}

test "scop_parse_des_line" {
  let line = "1000\tpx\ta.1.1.1.1.1.1\td1hba_\t1hba Hemoglobin alpha chain"
  let r = @src.scop_parse_des_line(line)
  assert_true(r.is_some())
  let rec = r.unwrap()
  assert_eq(rec.sunid(), 1000)
  assert_eq(rec.nodetype(), "px")
  assert_eq(rec.name(), "d1hba_")
}

test "scop_parse_des_line_class" {
  let line = "100\tcl\ta\t-\tAll alpha proteins"
  let r = @src.scop_parse_des_line(line)
  assert_true(r.is_some())
  let rec = r.unwrap()
  assert_eq(rec.nodetype(), "cl")
  assert_eq(rec.sccs(), "a")
  assert_eq(rec.name(), "-")
  assert_eq(rec.description(), "All alpha proteins")
}

test "scop_parse_hie_line" {
  let line = "600\t500\t1000,1001"
  let r = @src.scop_parse_hie_line(line)
  assert_true(r.is_some())
  let rec = r.unwrap()
  assert_eq(rec.sunid(), 600)
  assert_eq(rec.parent(), 500)
  assert_eq(rec.children().length(), 2)
}

test "scop_parse_hie_line_root" {
  let line = "0\t-\t100"
  let r = @src.scop_parse_hie_line(line)
  assert_true(r.is_some())
  let rec = r.unwrap()
  assert_eq(rec.sunid(), 0)
  assert_eq(rec.parent(), -1)
  assert_eq(rec.children().length(), 1)
}

test "scop_parse_hie_line_leaf" {
  let line = "1000\t600\t-"
  let r = @src.scop_parse_hie_line(line)
  assert_true(r.is_some())
  let rec = r.unwrap()
  assert_eq(rec.sunid(), 1000)
  assert_eq(rec.parent(), 600)
  assert_eq(rec.children().length(), 0)
}

test "scop_parse_cla_multiple" {
  let content = @src.scop_sample_cla()
  let records = @src.scop_parse_cla(content)
  assert_eq(records.length(), 2)
  assert_eq(records[0].sid(), "d1hba_")
  assert_eq(records[1].sid(), "d1hbb_")
}

test "scop_parse_des_multiple" {
  let content = @src.scop_sample_des()
  let records = @src.scop_parse_des(content)
  // root + 6 hierarchy nodes (cl/cf/sf/fa/dm/sp) + 2 domains = 9 records
  assert_eq(records.length(), 9)
}

test "scop_parse_hie_multiple" {
  let content = @src.scop_sample_hie()
  let records = @src.scop_parse_hie(content)
  assert_eq(records.length(), 9)
}

// ---------------------------------------------------------------------------
// Scop hierarchy construction and queries
// ---------------------------------------------------------------------------

test "scop_build_hierarchy" {
  let cla = @src.scop_parse_cla(@src.scop_sample_cla())
  let des = @src.scop_parse_des(@src.scop_sample_des())
  let hie = @src.scop_parse_hie(@src.scop_sample_hie())
  let scop = @src.Scop::new(cla, des, hie)
  // Root node
  assert_eq(scop.root().sunid(), 0)
  assert_eq(scop.root().nodetype(), "ro")
  // Check domain count
  assert_eq(scop.get_domains().length(), 2)
}

test "scop_get_node_by_sunid" {
  let cla = @src.scop_parse_cla(@src.scop_sample_cla())
  let des = @src.scop_parse_des(@src.scop_sample_des())
  let hie = @src.scop_parse_hie(@src.scop_sample_hie())
  let scop = @src.Scop::new(cla, des, hie)
  let node = scop.get_node_by_sunid(100)
  assert_true(node.is_some())
  assert_eq(node.unwrap().nodetype(), "cl")
  assert_eq(node.unwrap().description(), "All alpha proteins")
  // Non-existent
  let none = scop.get_node_by_sunid(99999)
  assert_true(none.is_none())
}

test "scop_get_domain_by_sid" {
  let cla = @src.scop_parse_cla(@src.scop_sample_cla())
  let des = @src.scop_parse_des(@src.scop_sample_des())
  let hie = @src.scop_parse_hie(@src.scop_sample_hie())
  let scop = @src.Scop::new(cla, des, hie)
  let d = scop.get_domain_by_sid("d1hba_")
  assert_true(d.is_some())
  assert_eq(d.unwrap().sunid(), 1000)
  assert_eq(d.unwrap().sid(), "d1hba_")
  assert_eq(d.unwrap().pdbid(), "1hba")
  assert_eq(d.unwrap().sccs(), "a.1.1.1.1.1.1")
  // Non-existent
  let none = scop.get_domain_by_sid("d9zzzz_")
  assert_true(none.is_none())
}

test "scop_get_parent" {
  let cla = @src.scop_parse_cla(@src.scop_sample_cla())
  let des = @src.scop_parse_des(@src.scop_sample_des())
  let hie = @src.scop_parse_hie(@src.scop_sample_hie())
  let scop = @src.Scop::new(cla, des, hie)
  // Parent of domain 1000 is species 600
  let d = scop.get_domain_by_sid("d1hba_").unwrap()
  let p = scop.get_parent(d)
  assert_true(p.is_some())
  assert_eq(p.unwrap().sunid(), 600)
  assert_eq(p.unwrap().nodetype(), "sp")
  // Root has no parent
  let root_parent = scop.get_parent(scop.root())
  assert_true(root_parent.is_none())
}

test "scop_get_children" {
  let cla = @src.scop_parse_cla(@src.scop_sample_cla())
  let des = @src.scop_parse_des(@src.scop_sample_des())
  let hie = @src.scop_parse_hie(@src.scop_sample_hie())
  let scop = @src.Scop::new(cla, des, hie)
  // Children of species 600 are domains 1000 and 1001
  let sp = scop.get_node_by_sunid(600).unwrap()
  let children = scop.get_children(sp)
  assert_eq(children.length(), 2)
  assert_eq(children[0].sunid(), 1000)
  assert_eq(children[1].sunid(), 1001)
  // Leaf node has no children
  let leaf = scop.get_domain_by_sid("d1hba_").unwrap()
  assert_eq(scop.get_children(leaf).length(), 0)
}

test "scop_get_ascendent_by_code" {
  let cla = @src.scop_parse_cla(@src.scop_sample_cla())
  let des = @src.scop_parse_des(@src.scop_sample_des())
  let hie = @src.scop_parse_hie(@src.scop_sample_hie())
  let scop = @src.Scop::new(cla, des, hie)
  let d = scop.get_domain_by_sid("d1hba_").unwrap()
  // Ascendent of type "cl" (class)
  let cl = scop.get_ascendent(d, "cl")
  assert_true(cl.is_some())
  assert_eq(cl.unwrap().sunid(), 100)
  assert_eq(cl.unwrap().description(), "All alpha proteins")
  // Ascendent of type "fa" (family)
  let fa = scop.get_ascendent(d, "fa")
  assert_true(fa.is_some())
  assert_eq(fa.unwrap().sunid(), 400)
}

test "scop_get_ascendent_by_long_name" {
  let cla = @src.scop_parse_cla(@src.scop_sample_cla())
  let des = @src.scop_parse_des(@src.scop_sample_des())
  let hie = @src.scop_parse_hie(@src.scop_sample_hie())
  let scop = @src.Scop::new(cla, des, hie)
  let d = scop.get_domain_by_sid("d1hba_").unwrap()
  // Using long name "fold" instead of "cf"
  let fold = scop.get_ascendent(d, "fold")
  assert_true(fold.is_some())
  assert_eq(fold.unwrap().sunid(), 200)
  assert_eq(fold.unwrap().description(), "Globin-like")
}

test "scop_get_descendents" {
  let cla = @src.scop_parse_cla(@src.scop_sample_cla())
  let des = @src.scop_parse_des(@src.scop_sample_des())
  let hie = @src.scop_parse_hie(@src.scop_sample_hie())
  let scop = @src.Scop::new(cla, des, hie)
  // Descendents of root (type "px") should be all 2 domains
  let domains = scop.get_descendents(scop.root(), "px")
  assert_eq(domains.length(), 2)
  // Descendents of class 100 (type "px") should also be 2 domains
  let cl = scop.get_node_by_sunid(100).unwrap()
  let domains2 = scop.get_descendents(cl, "px")
  assert_eq(domains2.length(), 2)
  // Descendents of root (type "fa") should be 1 family
  let families = scop.get_descendents(scop.root(), "fa")
  assert_eq(families.length(), 1)
  assert_eq(families[0].sunid(), 400)
}

test "scop_get_descendents_long_name" {
  let cla = @src.scop_parse_cla(@src.scop_sample_cla())
  let des = @src.scop_parse_des(@src.scop_sample_des())
  let hie = @src.scop_parse_hie(@src.scop_sample_hie())
  let scop = @src.Scop::new(cla, des, hie)
  let domains = scop.get_descendents(scop.root(), "domain")
  assert_eq(domains.length(), 2)
}

test "scop_node_is_domain" {
  let cla = @src.scop_parse_cla(@src.scop_sample_cla())
  let des = @src.scop_parse_des(@src.scop_sample_des())
  let hie = @src.scop_parse_hie(@src.scop_sample_hie())
  let scop = @src.Scop::new(cla, des, hie)
  let d = scop.get_domain_by_sid("d1hba_").unwrap()
  assert_true(d.is_domain())
  let cl = scop.get_node_by_sunid(100).unwrap()
  assert_true(!cl.is_domain())
}

// ---------------------------------------------------------------------------
// SCCS comparison
// ---------------------------------------------------------------------------

test "scop_cmp_sccs_equal" {
  assert_eq(@src.scop_cmp_sccs("a.1.1.1", "a.1.1.1"), 0)
}

test "scop_cmp_sccs_letter_diff" {
  assert_true(@src.scop_cmp_sccs("a.1.1.1", "b.1.1.1") < 0)
  assert_true(@src.scop_cmp_sccs("b.1.1.1", "a.1.1.1") > 0)
}

test "scop_cmp_sccs_numeric_diff" {
  assert_true(@src.scop_cmp_sccs("a.1.1.1", "a.1.1.2") < 0)
  assert_true(@src.scop_cmp_sccs("a.1.1.2", "a.1.1.1") > 0)
  assert_true(@src.scop_cmp_sccs("a.1.1.1", "a.1.2.1") < 0)
}

test "scop_cmp_sccs_length_diff" {
  // Shorter prefix sorts first when all compared components are equal
  assert_true(@src.scop_cmp_sccs("a.1", "a.1.1") < 0)
  assert_true(@src.scop_cmp_sccs("a.1.1", "a.1") > 0)
}

test "scop_cmp_sccs_numeric_not_lexical" {
  // Numerically, 2 < 11, so a.1.2 < a.1.11 (NOT lexical where "2" > "11")
  assert_true(@src.scop_cmp_sccs("a.1.2", "a.1.11") < 0)
  assert_true(@src.scop_cmp_sccs("a.1.11", "a.1.2") > 0)
}

// ---------------------------------------------------------------------------
// Serialization
// ---------------------------------------------------------------------------

test "scop_write_hie" {
  let cla = @src.scop_parse_cla(@src.scop_sample_cla())
  let des = @src.scop_parse_des(@src.scop_sample_des())
  let hie = @src.scop_parse_hie(@src.scop_sample_hie())
  let scop = @src.Scop::new(cla, des, hie)
  let output = scop.write_hie()
  // Should contain root line
  assert_true(output.contains("0\t-\t100"))
  // Should contain a domain line
  assert_true(output.contains("1000\t600\t-"))
}

test "scop_write_des" {
  let cla = @src.scop_parse_cla(@src.scop_sample_cla())
  let des = @src.scop_parse_des(@src.scop_sample_des())
  let hie = @src.scop_parse_hie(@src.scop_sample_hie())
  let scop = @src.Scop::new(cla, des, hie)
  let output = scop.write_des()
  assert_true(output.contains("All alpha proteins"))
  assert_true(output.contains("d1hba_"))
}

test "scop_write_cla" {
  let cla = @src.scop_parse_cla(@src.scop_sample_cla())
  let des = @src.scop_parse_des(@src.scop_sample_des())
  let hie = @src.scop_parse_hie(@src.scop_sample_hie())
  let scop = @src.Scop::new(cla, des, hie)
  let output = scop.write_cla()
  assert_true(output.contains("d1hba_"))
  assert_true(output.contains("1hba"))
  assert_true(output.contains("cl=100"))
}

// ---------------------------------------------------------------------------
// Empty Scop
// ---------------------------------------------------------------------------

test "scop_empty" {
  let scop = @src.scop_empty()
  assert_eq(scop.root().sunid(), 0)
  assert_eq(scop.root().nodetype(), "ro")
  assert_eq(scop.get_domains().length(), 0)
  assert_true(scop.get_node_by_sunid(100).is_none())
  assert_true(scop.get_domain_by_sid("d1hba_").is_none())
}
