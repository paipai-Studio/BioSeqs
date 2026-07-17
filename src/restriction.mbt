///|
/// Restriction - Restriction enzyme analysis inspired by Biopython Bio.Restriction.

///|
pub struct RestrictionEnzyme {
  name : String
  recognition_site : String
  cutting_position : Int
  overhang_5 : String
  overhang_3 : String
}

///|
pub fn restriction_enzyme(
  name : String,
  recognition_site : String,
  cutting_position : Int,
  overhang_5 : String,
  overhang_3 : String,
) -> RestrictionEnzyme {
  RestrictionEnzyme::{
    name,
    recognition_site,
    cutting_position,
    overhang_5,
    overhang_3,
  }
}

///|
pub fn builtin_enzymes() -> Map[String, RestrictionEnzyme] {
  let enzymes : Map[String, RestrictionEnzyme] = Map([], capacity=50)

  enzymes["EcoRI"] = RestrictionEnzyme::{
    name: "EcoRI",
    recognition_site: "GAATTC",
    cutting_position: 1,
    overhang_5: "AATT",
    overhang_3: "",
  }
  enzymes["HindIII"] = RestrictionEnzyme::{
    name: "HindIII",
    recognition_site: "AAGCTT",
    cutting_position: 1,
    overhang_5: "AGCT",
    overhang_3: "",
  }
  enzymes["BamHI"] = RestrictionEnzyme::{
    name: "BamHI",
    recognition_site: "GGATCC",
    cutting_position: 1,
    overhang_5: "GATC",
    overhang_3: "",
  }
  enzymes["NotI"] = RestrictionEnzyme::{
    name: "NotI",
    recognition_site: "GCGGCCGC",
    cutting_position: 2,
    overhang_5: "GGCC",
    overhang_3: "",
  }
  enzymes["XhoI"] = RestrictionEnzyme::{
    name: "XhoI",
    recognition_site: "CTCGAG",
    cutting_position: 1,
    overhang_5: "TCGA",
    overhang_3: "",
  }
  enzymes["SalI"] = RestrictionEnzyme::{
    name: "SalI",
    recognition_site: "GTCGAC",
    cutting_position: 1,
    overhang_5: "TCGA",
    overhang_3: "",
  }
  enzymes["PstI"] = RestrictionEnzyme::{
    name: "PstI",
    recognition_site: "CTGCAG",
    cutting_position: 4,
    overhang_5: "",
    overhang_3: "TGCA",
  }
  enzymes["KpnI"] = RestrictionEnzyme::{
    name: "KpnI",
    recognition_site: "GGTACC",
    cutting_position: 5,
    overhang_5: "",
    overhang_3: "GTAC",
  }
  enzymes["SmaI"] = RestrictionEnzyme::{
    name: "SmaI",
    recognition_site: "CCCGGG",
    cutting_position: 3,
    overhang_5: "",
    overhang_3: "",
  }
  enzymes["EcoRV"] = RestrictionEnzyme::{
    name: "EcoRV",
    recognition_site: "GATATC",
    cutting_position: 3,
    overhang_5: "",
    overhang_3: "",
  }
  enzymes["BglII"] = RestrictionEnzyme::{
    name: "BglII",
    recognition_site: "AGATCT",
    cutting_position: 1,
    overhang_5: "GATC",
    overhang_3: "",
  }
  enzymes["XbaI"] = RestrictionEnzyme::{
    name: "XbaI",
    recognition_site: "TCTAGA",
    cutting_position: 1,
    overhang_5: "CTAG",
    overhang_3: "",
  }
  enzymes["SacI"] = RestrictionEnzyme::{
    name: "SacI",
    recognition_site: "GAGCTC",
    cutting_position: 5,
    overhang_5: "",
    overhang_3: "AGCT",
  }
  enzymes["SpeI"] = RestrictionEnzyme::{
    name: "SpeI",
    recognition_site: "ACTAGT",
    cutting_position: 1,
    overhang_5: "CTAG",
    overhang_3: "",
  }
  enzymes["NdeI"] = RestrictionEnzyme::{
    name: "NdeI",
    recognition_site: "CATATG",
    cutting_position: 1,
    overhang_5: "ATAT",
    overhang_3: "",
  }
  enzymes["NheI"] = RestrictionEnzyme::{
    name: "NheI",
    recognition_site: "GCTAGC",
    cutting_position: 1,
    overhang_5: "CTAG",
    overhang_3: "",
  }
  enzymes["ApaI"] = RestrictionEnzyme::{
    name: "ApaI",
    recognition_site: "GGGCCC",
    cutting_position: 5,
    overhang_5: "",
    overhang_3: "GGCC",
  }
  enzymes["BstXI"] = RestrictionEnzyme::{
    name: "BstXI",
    recognition_site: "CCANNNNNNNTGG",
    cutting_position: 6,
    overhang_5: "",
    overhang_3: "",
  }

  enzymes
}

///|
pub fn find_cut_sites(seq : String, enzyme : RestrictionEnzyme) -> Array[Int] {
  let sites : Array[Int] = Array::new()
  let pattern = enzyme.recognition_site
  let seq_len = seq.length()
  let pattern_len = pattern.length()

  for i = 0; i <= seq_len - pattern_len; i = i + 1 {
    let mut match_found = true
    for j = 0; j < pattern_len; j = j + 1 {
      let seq_char = seq.unsafe_get(i + j).to_string()
      let pattern_char = pattern.unsafe_get(j).to_string()

      if pattern_char != "N" && seq_char != pattern_char {
        match_found = false
        break
      }
    }

    if match_found {
      sites.push(i + enzyme.cutting_position)
    }
  }

  sites
}

///|
pub fn digest(seq : String, enzyme : RestrictionEnzyme) -> Array[String] {
  let sites = find_cut_sites(seq, enzyme)
  let fragments : Array[String] = Array::new()

  let mut prev_pos = 0
  for i = 0; i < sites.length(); i = i + 1 {
    fragments.push(substring(seq, prev_pos, sites[i] - prev_pos))
    prev_pos = sites[i]
  }

  if prev_pos < seq.length() {
    fragments.push(substring(seq, prev_pos, seq.length() - prev_pos))
  }

  fragments
}

///|
pub fn digest_with_multiple(
  seq : String,
  enzymes : Array[RestrictionEnzyme],
) -> Array[String] {
  let all_sites : Array[Int] = Array::new()

  for i = 0; i < enzymes.length(); i = i + 1 {
    let sites = find_cut_sites(seq, enzymes[i])
    for j = 0; j < sites.length(); j = j + 1 {
      all_sites.push(sites[j])
    }
  }

  for i = 0; i < all_sites.length() - 1; i = i + 1 {
    for j = i + 1; j < all_sites.length(); j = j + 1 {
      if all_sites[i] > all_sites[j] {
        let temp = all_sites[i]
        all_sites[i] = all_sites[j]
        all_sites[j] = temp
      }
    }
  }

  let fragments : Array[String] = Array::new()
  let mut prev_pos = 0

  for i = 0; i < all_sites.length(); i = i + 1 {
    if all_sites[i] > prev_pos {
      fragments.push(substring(seq, prev_pos, all_sites[i] - prev_pos))
      prev_pos = all_sites[i]
    }
  }

  if prev_pos < seq.length() {
    fragments.push(substring(seq, prev_pos, seq.length() - prev_pos))
  }

  fragments
}

///|
pub fn fragment_sizes(fragments : Array[String]) -> Array[Int] {
  let sizes : Array[Int] = Array::make(fragments.length(), 0)
  for i = 0; i < fragments.length(); i = i + 1 {
    sizes[i] = fragments[i].length()
  }
  sizes
}

///|
pub fn is_compatible_overhang(
  enzyme1 : RestrictionEnzyme,
  enzyme2 : RestrictionEnzyme,
) -> Bool {
  enzyme1.overhang_5 == enzyme2.overhang_5 ||
  enzyme1.overhang_3 == enzyme2.overhang_3
}

///|
pub fn get_enzyme_by_name(name : String) -> RestrictionEnzyme? {
  let enzymes = builtin_enzymes()
  enzymes.get(name)
}

///|
pub fn list_enzymes() -> Array[String] {
  builtin_enzymes().keys().collect()
}

///|
pub fn enzymes_by_overhang(overhang : String) -> Array[RestrictionEnzyme] {
  let enzymes = builtin_enzymes()
  let keys = enzymes.keys().collect()
  let result : Array[RestrictionEnzyme] = Array::new()

  for i = 0; i < keys.length(); i = i + 1 {
    let enzyme = enzymes[keys[i]]
    if enzyme.overhang_5 == overhang || enzyme.overhang_3 == overhang {
      result.push(enzyme)
    }
  }

  result
}

///|
pub fn enzymes_by_site_length(
  min_len : Int,
  max_len : Int,
) -> Array[RestrictionEnzyme] {
  let enzymes = builtin_enzymes()
  let keys = enzymes.keys().collect()
  let result : Array[RestrictionEnzyme] = Array::new()

  for i = 0; i < keys.length(); i = i + 1 {
    let enzyme = enzymes[keys[i]]
    let len = enzyme.recognition_site.replace(old="N", new="").length()
    if len >= min_len && len <= max_len {
      result.push(enzyme)
    }
  }

  result
}
