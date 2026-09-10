///|
/// Bio.Restriction demo: a REBASE-style restriction map of a small linear
/// construct and a circular mini-plasmid.
///
/// Showcases: the 1088-enzyme database, batch search / RestrictionAnalysis
/// filters, single and double digest fragment tables (linear and circular),
/// overhang classification, compatible ends and isoschizomer relationships.
///
/// Run with: moon run examples/restriction_demo
fn fmt_cuts(cuts : Array[Int]) -> String {
  if cuts.length() == 0 {
    "-"
  } else {
    let mut s = ""
    for i in 0..<cuts.length() {
      if i > 0 {
        s += ", "
      }
      s += cuts[i].to_string()
    }
    s
  }
}

///|
fn fmt_names(enzymes : Array[@bio.BioRestrictionEnzyme]) -> String {
  if enzymes.length() == 0 {
    "-"
  } else {
    let mut s = ""
    for i in 0..<enzymes.length() {
      if i > 0 {
        s += ", "
      }
      s += enzymes[i].name
    }
    s
  }
}

///|
fn fmt_strs(xs : Array[String]) -> String {
  if xs.length() == 0 {
    "-"
  } else {
    let mut s = ""
    for i in 0..<xs.length() {
      if i > 0 {
        s += ", "
      }
      s += xs[i]
    }
    s
  }
}

///|
/// Fragment lengths of a digest from sorted cut positions.
fn digest_lengths(cuts : Array[Int], length : Int, linear : Bool) -> Array[Int] {
  let c = cuts.copy()
  c.sort()
  let r : Array[Int] = []
  if c.length() == 0 {
    r.push(length)
    return r
  }
  if linear {
    let mut prev = 0
    for pos in c {
      r.push(pos - prev)
      prev = pos
    }
    r.push(length - prev)
  } else {
    for i in 1..<c.length() {
      r.push(c[i] - c[i - 1])
    }
    // bridge the last site through the origin back to the first site
    r.push(length - c[c.length() - 1] + c[0])
  }
  r
}

///|
fn main raise {
  // ------------------------------------------------------------------
  // 1. Database overview
  // ------------------------------------------------------------------
  println("== Bio.Restriction: REBASE enzyme database ==")
  println(
    "enzymes in database: " + @bio.bio_restriction_enzyme_count().to_string(),
  )
  let ecori = @bio.require_bio_restriction_enzyme("EcoRI")
  println(
    "EcoRI: site=\{ecori.site}  elucidate: \{ecori.elucidate()}  " +
    "palindromic=\{ecori.is_palindromic()}  " +
    "theoretical frequency: one cut per \{ecori.frequency().to_string()} bp",
  )
  println(
    "end type: \{ecori.overhang_type()}  suppliers (" +
    ecori.supplier_list().length().to_string() +
    ")",
  )
  let bbsi = @bio.require_bio_restriction_enzyme("BbsI")
  println("BbsI (cuts outside its recognition site): \{bbsi.elucidate()}")
  println("")

  // ------------------------------------------------------------------
  // 2. Restriction map of a 96 bp linear construct
  // ------------------------------------------------------------------
  let construct = "ACGTACGTAC" +
    "GAATTC" +
    "TTTT" +
    "CCCGGG" +
    "TTTT" +
    "GGATCC" +
    "TTTT" +
    "GGTACC" +
    "TTTT" +
    "TCTAGA" +
    "AAAATTTTGGAAAATTTTGGAAAATTTTGGAAAATTTTGG"
  let batch = @bio.RestrictionBatch::new(names=[
    "EcoRI", "SmaI", "BamHI", "KpnI", "XbaI", "BstYI", "AvaI", "HindIII",
  ])
  println("== restriction map of a 96 bp linear construct ==")
  println(construct)
  println("length: " + construct.length().to_string() + " bp")
  println("")
  let sites = batch.search(construct)
  for name in batch.elements() {
    let e = @bio.require_bio_restriction_enzyme(name)
    let cuts = sites.get(name).unwrap()
    let mut label = name
    while label.length() < 8 {
      label += " "
    }
    println(
      "\{label} \{e.site}  cuts at [\{fmt_cuts(cuts)}]  " +
      "(\{e.overhang_type()})",
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Single digest: EcoRI fragments, linear and circular
  // ------------------------------------------------------------------
  println("== single digest with EcoRI ==")
  let lin_frag = ecori.catalyse(construct)
  println("linear fragments: " + fmt_strs(lin_frag))
  let ring = "GGGGAATTCAAAAAAGAATTCGGG"
  let circ_frag = ecori.catalyse(ring, linear=false)
  println("circular 24 bp mini-plasmid: " + ring)
  let mut lengths = ""
  for f in circ_frag {
    if lengths != "" {
      lengths += ", "
    }
    lengths += f.length().to_string()
  }
  println("circular fragments (origin-bridging digest): " + lengths + " bp")
  println("")

  // ------------------------------------------------------------------
  // 4. Double digest: EcoRI + BamHI fragment table
  // ------------------------------------------------------------------
  println("== double digest (EcoRI + BamHI), linear ==")
  let eco_cuts = ecori.search(construct)
  let bam_cuts = @bio.require_bio_restriction_enzyme("BamHI").search(construct)
  let both : Array[Int] = []
  for x in eco_cuts {
    both.push(x)
  }
  for x in bam_cuts {
    both.push(x)
  }
  let double_lengths = digest_lengths(both, construct.length(), true)
  println(
    "cut positions: " +
    fmt_cuts(both) +
    "  fragment lengths: " +
    fmt_cuts(double_lengths) +
    " bp",
  )
  println("")

  // ------------------------------------------------------------------
  // 5. End chemistry and compatible (sticky) ends
  // ------------------------------------------------------------------
  println("== ends and isoschizomers ==")
  println(
    "SmaI  CCC^GGG  is blunt: " +
    @bio.require_bio_restriction_enzyme("SmaI").is_blunt().to_string(),
  )
  println(
    "KpnI  GGTAC^C is 3' overhang: " +
    @bio.require_bio_restriction_enzyme("KpnI").is_3overhang().to_string(),
  )
  let mfei = @bio.require_bio_restriction_enzyme("MfeI")
  let xhoi = @bio.require_bio_restriction_enzyme("XhoI")
  println(
    "enzymes producing ends compatible with EcoRI: " +
    ecori.compatible_ends().length().to_string(),
  )
  println("  EcoRI % MfeI: " + ecori.compatible_overhang(mfei).to_string())
  println("  EcoRI % XhoI: " + ecori.compatible_overhang(xhoi).to_string())
  let saci = @bio.require_bio_restriction_enzyme("SacI")
  let ssti = @bio.require_bio_restriction_enzyme("SstI")
  println(
    "SacI/SstI are isoschizomers with the same cut: " +
    saci.is_equischizomer(ssti).to_string(),
  )
  let smai = @bio.require_bio_restriction_enzyme("SmaI")
  let xmai = @bio.require_bio_restriction_enzyme("XmaI")
  println(
    "SmaI/XmaI share the site but cut differently (neoschizomers): " +
    smai.is_neoschizomer(xmai).to_string(),
  )
  println("equischizomers of SacI: " + fmt_names(saci.equischizomers()))
  println("neoschizomers of SmaI: " + fmt_names(smai.neoschizomers()))
  println("")

  // ------------------------------------------------------------------
  // 6. RestrictionAnalysis: digest summaries and region filters
  // ------------------------------------------------------------------
  println("== RestrictionAnalysis summaries ==")
  let analysis = @bio.RestrictionAnalysis::new(batch~, seq=construct)
  println(
    "enzymes with exactly one site: " +
    analysis.with_n_sites(1).length().to_string(),
  )
  println(
    "enzymes with no site: " + analysis.without_site().length().to_string(),
  )
  println(
    "blunt cutters in the batch: " + analysis.blunt().length().to_string(),
  )
  println("5' overhang cutters: " + analysis.overhang5().length().to_string())
  println("3' overhang cutters: " + analysis.overhang3().length().to_string())
  let window = analysis.between(20, 40)
  let mut window_names = ""
  for name in batch.elements() {
    match window.get(name) {
      Some(cuts) =>
        if cuts.length() > 0 {
          if window_names != "" {
            window_names += ", "
          }
          window_names += name + " [" + fmt_cuts(cuts) + "]"
        }
      None => ()
    }
  }
  println("enzymes cutting in region 20..40: " + window_names)
  println(
    "enzymes that do not cut anywhere: " +
      analysis.do_not_cut(1, construct.length() + 1).length().to_string(),
  )
  println("")

  // ------------------------------------------------------------------
  // 7. Supplier catalogue
  // ------------------------------------------------------------------
  println("== supplier catalogue ==")
  println(
    "enzymes sold by supplier code 'N' (e.g. NEB): " +
    @bio.bio_restriction_enzymes_by_supplier("N").length().to_string(),
  )
}
