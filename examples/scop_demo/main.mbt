///|
fn main {
  println("=== SCOP (Structural Classification of Proteins) Demo ===")
  println("Parsing SCOP database files and building hierarchy.")
  println("")

  // 1. Parse CLA, DES, HIE files
  println("1. Parse SCOP Files:")
  println("====================")
  let cla_content = @bio.scop_sample_cla()
  let des_content = @bio.scop_sample_des()
  let hie_content = @bio.scop_sample_hie()
  let cla_records = @bio.scop_parse_cla(cla_content)
  let des_records = @bio.scop_parse_des(des_content)
  let hie_records = @bio.scop_parse_hie(hie_content)
  println("CLA records: " + cla_records.length().to_string())
  println("DES records: " + des_records.length().to_string())
  println("HIE records: " + hie_records.length().to_string())

  // 2. Build hierarchy
  println("\n2. Build SCOP Hierarchy:")
  println("========================")
  let scop = @bio.Scop::new(cla_records, des_records, hie_records)
  println("Root sunid: " + scop.root().sunid().to_string())
  println("Root type: " + scop.root().nodetype() + " (" + @bio.scop_node_type_name(scop.root().nodetype()) + ")")
  println("Total domains: " + scop.get_domains().length().to_string())

  // 3. Look up domain by sid
  println("\n3. Look Up Domain by SID:")
  println("=========================")
  let d = scop.get_domain_by_sid("d1hba_")
  if d.is_some() {
    let domain = d.unwrap()
    println("SID: " + domain.sid())
    println("Sunid: " + domain.sunid().to_string())
    println("SCCS: " + domain.sccs())
    println("PDB ID: " + domain.pdbid())
    println("Description: " + domain.description())
    let res = domain.residues()
    if res.is_some() {
      let r = res.unwrap()
      println("PDB ID (residues): " + r.pdbid())
      println("Fragments: " + r.fragments().length().to_string())
      for frag in r.fragments() {
        println("  Chain " + frag.chain() + ": " + frag.start() + "-" + frag.end_())
      }
    }
  }

  // 4. Traverse hierarchy (get ascendents)
  println("\n4. Get Ascendents of Domain:")
  println("============================")
  if d.is_some() {
    let domain = d.unwrap()
    let types = ["class", "fold", "superfamily", "family", "protein", "species"]
    for t in types {
      let asc = scop.get_ascendent(domain, t)
      if asc.is_some() {
        let n = asc.unwrap()
        println(t + ": sunid=" + n.sunid().to_string() + " desc=\"" + n.description() + "\"")
      }
    }
  }

  // 5. Get descendents of root
  println("\n5. Get Descendents of Root:")
  println("===========================")
  let domains = scop.get_descendents(scop.root(), "domain")
  println("Domain descendents of root: " + domains.length().to_string())
  for dom in domains {
    println("  " + dom.sid() + " (sunid=" + dom.sunid().to_string() + ")")
  }

  // 6. SCCS comparison
  println("\n6. SCCS Comparison:")
  println("===================")
  let sccs1 = "a.1.1.1"
  let sccs2 = "a.1.1.2"
  let cmp = @bio.scop_cmp_sccs(sccs1, sccs2)
  println(sccs1 + " vs " + sccs2 + ": " + cmp.to_string())
  let cmp2 = @bio.scop_cmp_sccs("b.1.1.1", "a.1.1.1")
  println("b.1.1.1 vs a.1.1.1: " + cmp2.to_string())

  // 7. Serialize back to files
  println("\n7. Serialize SCOP Database:")
  println("===========================")
  let hie_out = scop.write_hie()
  let des_out = scop.write_des()
  let cla_out = scop.write_cla()
  println("HIE output lines: " + count_lines(hie_out).to_string())
  println("DES output lines: " + count_lines(des_out).to_string())
  println("CLA output lines: " + count_lines(cla_out).to_string())

  println("\n=== Demo Complete ===")
}

fn count_lines(s : String) -> Int {
  let mut count = 1
  for i = 0; i < s.length(); i = i + 1 {
    if s.unsafe_get(i).to_int() == '\n'.to_int() {
      count = count + 1
    }
  }
  count
}
