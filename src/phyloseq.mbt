///|
/// phyloseq - Microbiome analysis.

///|
pub struct OTU {
  id : String
  taxonomy : Array[String]
  abundance : Array[Double]
}

///|
pub fn OTU::new(id : String, taxonomy : Array[String], abundance : Array[Double]) -> OTU {
  OTU::{ id, taxonomy, abundance }
}

///|
pub struct Sample {
  id : String
  metadata : Map[String, String]
}

///|
pub fn Sample::new(id : String, metadata : Map[String, String]) -> Sample {
  Sample::{ id, metadata }
}

///|
pub struct TaxRank {
  kingdom : String
  phylum : String
  class : String
  order : String
  family : String
  genus : String
  species : String
}

///|
pub fn TaxRank::new(kingdom : String, phylum : String, class : String, order : String, family : String, genus : String, species : String) -> TaxRank {
  TaxRank::{ kingdom, phylum, class, order, family, genus, species }
}

///|
pub struct Phyloseq {
  otus : Array[OTU]
  samples : Array[Sample]
  taxonomy : Array[TaxRank]
  metadata : Map[String, String]
}

///|
pub fn Phyloseq::new() -> Phyloseq {
  Phyloseq::{
    otus: Array::new(),
    samples: Array::new(),
    taxonomy: Array::new(),
    metadata: Map([], capacity=0)
  }
}

///|
pub fn ps_add_otu(ps : Phyloseq, otu : OTU) -> Phyloseq {
  ps.otus.push(otu)
  ps
}

///|
pub fn ps_add_sample(ps : Phyloseq, sample : Sample) -> Phyloseq {
  ps.samples.push(sample)
  ps
}

///|
pub fn ps_add_taxonomy(ps : Phyloseq, tax : TaxRank) -> Phyloseq {
  ps.taxonomy.push(tax)
  ps
}

///|
pub fn ps_set_metadata(ps : Phyloseq, key : String, value : String) -> Phyloseq {
  ps.metadata.set(key, value)
  ps
}

///|
pub fn ps_num_otus(ps : Phyloseq) -> Int {
  ps.otus.length()
}

///|
pub fn ps_num_samples(ps : Phyloseq) -> Int {
  ps.samples.length()
}

///|
pub fn ps_total_abundance(ps : Phyloseq) -> Double {
  let mut total = 0.0
  let mut i = 0
  while i < ps.otus.length() {
    let mut j = 0
    while j < ps.otus[i].abundance.length() {
      total = total + ps.otus[i].abundance[j]
      j = j + 1
    }
    i = i + 1
  }
  total
}

///|
pub fn ps_filter_by_abundance(ps : Phyloseq, min_abundance : Double) -> Phyloseq {
  let filtered = Phyloseq::new()
  
  let mut i = 0
  while i < ps.samples.length() {
    filtered.samples.push(ps.samples[i])
    i = i + 1
  }
  
  let mut j = 0
  while j < ps.taxonomy.length() {
    filtered.taxonomy.push(ps.taxonomy[j])
    j = j + 1
  }
  
  let mut k = 0
  while k < ps.otus.length() {
    let mut total = 0.0
    let mut m = 0
    while m < ps.otus[k].abundance.length() {
      total = total + ps.otus[k].abundance[m]
      m = m + 1
    }
    if total >= min_abundance {
      filtered.otus.push(ps.otus[k])
    }
    k = k + 1
  }
  
  filtered
}

///|
pub fn ps_filter_by_taxonomy(ps : Phyloseq, rank : String, value : String) -> Phyloseq {
  let filtered = Phyloseq::new()
  
  let mut i = 0
  while i < ps.samples.length() {
    filtered.samples.push(ps.samples[i])
    i = i + 1
  }
  
  let mut j = 0
  while j < ps.taxonomy.length() {
    let matched = match rank {
      "kingdom" => ps.taxonomy[j].kingdom == value
      "phylum" => ps.taxonomy[j].phylum == value
      "class" => ps.taxonomy[j].class == value
      "order" => ps.taxonomy[j].order == value
      "family" => ps.taxonomy[j].family == value
      "genus" => ps.taxonomy[j].genus == value
      "species" => ps.taxonomy[j].species == value
      _ => false
    }
    if matched {
      filtered.taxonomy.push(ps.taxonomy[j])
      if j < ps.otus.length() {
        filtered.otus.push(ps.otus[j])
      }
    }
    j = j + 1
  }
  
  filtered
}

///|
pub fn ps_get_taxa_summary(ps : Phyloseq, rank : String) -> Map[String, Double] {
  let summary = Map([], capacity=0)
  
  let mut i = 0
  while i < ps.taxonomy.length() {
    let tax_value = match rank {
      "kingdom" => ps.taxonomy[i].kingdom
      "phylum" => ps.taxonomy[i].phylum
      "class" => ps.taxonomy[i].class
      "order" => ps.taxonomy[i].order
      "family" => ps.taxonomy[i].family
      "genus" => ps.taxonomy[i].genus
      "species" => ps.taxonomy[i].species
      _ => "unknown"
    }
    
    let mut abundance = 0.0
    if i < ps.otus.length() {
      let mut j = 0
      while j < ps.otus[i].abundance.length() {
        abundance = abundance + ps.otus[i].abundance[j]
        j = j + 1
      }
    }
    
    if summary.contains(tax_value) {
      summary.set(tax_value, summary.get(tax_value).unwrap() + abundance)
    } else {
      summary.set(tax_value, abundance)
    }
    
    i = i + 1
  }
  
  summary
}

///|
pub fn ps_relative_abundance(ps : Phyloseq) -> Phyloseq {
  let total = ps_total_abundance(ps)
  let ra = Phyloseq::new()
  
  let mut i = 0
  while i < ps.samples.length() {
    ra.samples.push(ps.samples[i])
    i = i + 1
  }
  
  let mut j = 0
  while j < ps.taxonomy.length() {
    ra.taxonomy.push(ps.taxonomy[j])
    j = j + 1
  }
  
  let mut k = 0
  while k < ps.otus.length() {
    let ra_abundance : Array[Double] = Array::new()
    let mut m = 0
    while m < ps.otus[k].abundance.length() {
      ra_abundance.push(ps.otus[k].abundance[m] / total)
      m = m + 1
    }
    let ra_otu = OTU::new(ps.otus[k].id, ps.otus[k].taxonomy, ra_abundance)
    ra.otus.push(ra_otu)
    k = k + 1
  }
  
  ra
}

///|
pub fn ps_rarefy(ps : Phyloseq, depth : Int) -> Phyloseq {
  let rarefied = Phyloseq::new()
  
  let mut i = 0
  while i < ps.samples.length() {
    rarefied.samples.push(ps.samples[i])
    i = i + 1
  }
  
  let mut j = 0
  while j < ps.taxonomy.length() {
    rarefied.taxonomy.push(ps.taxonomy[j])
    j = j + 1
  }
  
  let mut k = 0
  while k < ps.otus.length() {
    let rarefied_abundance : Array[Double] = Array::new()
    let mut m = 0
    while m < ps.otus[k].abundance.length() {
      let proportion = ps.otus[k].abundance[m] / ps_total_abundance(ps)
      rarefied_abundance.push(proportion * (depth.to_double()))
      m = m + 1
    }
    let rarefied_otu = OTU::new(ps.otus[k].id, ps.otus[k].taxonomy, rarefied_abundance)
    rarefied.otus.push(rarefied_otu)
    k = k + 1
  }
  
  rarefied
}

///|
pub fn create_example_phyloseq() -> Phyloseq {
  let ps = Phyloseq::new()
  
  let otu1_tax : Array[String] = Array::new()
  otu1_tax.push("Bacteria")
  otu1_tax.push("Proteobacteria")
  otu1_tax.push("Gammaproteobacteria")
  let otu1_abund : Array[Double] = Array::new()
  otu1_abund.push(100.0)
  otu1_abund.push(150.0)
  otu1_abund.push(200.0)
  let otu1 = OTU::new("OTU001", otu1_tax, otu1_abund)
  ps.otus.push(otu1)
  
  let otu2_tax : Array[String] = Array::new()
  otu2_tax.push("Bacteria")
  otu2_tax.push("Firmicutes")
  otu2_tax.push("Clostridia")
  let otu2_abund : Array[Double] = Array::new()
  otu2_abund.push(50.0)
  otu2_abund.push(75.0)
  otu2_abund.push(100.0)
  let otu2 = OTU::new("OTU002", otu2_tax, otu2_abund)
  ps.otus.push(otu2)
  
  let otu3_tax : Array[String] = Array::new()
  otu3_tax.push("Bacteria")
  otu3_tax.push("Bacteroidetes")
  otu3_tax.push("Bacteroidia")
  let otu3_abund : Array[Double] = Array::new()
  otu3_abund.push(75.0)
  otu3_abund.push(100.0)
  otu3_abund.push(125.0)
  let otu3 = OTU::new("OTU003", otu3_tax, otu3_abund)
  ps.otus.push(otu3)
  
  let sample1_meta = Map([("group", "control")], capacity=1)
  let sample1 = Sample::new("Sample001", sample1_meta)
  ps.samples.push(sample1)
  
  let sample2_meta = Map([("group", "treatment")], capacity=1)
  let sample2 = Sample::new("Sample002", sample2_meta)
  ps.samples.push(sample2)
  
  let sample3_meta = Map([("group", "treatment")], capacity=1)
  let sample3 = Sample::new("Sample003", sample3_meta)
  ps.samples.push(sample3)
  
  let tax1 = TaxRank::new("Bacteria", "Proteobacteria", "Gammaproteobacteria", "Enterobacteriales", "Enterobacteriaceae", "Escherichia", "coli")
  ps.taxonomy.push(tax1)
  
  let tax2 = TaxRank::new("Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Clostridiaceae", "Clostridium", "perfringens")
  ps.taxonomy.push(tax2)
  
  let tax3 = TaxRank::new("Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Bacteroidaceae", "Bacteroides", "fragilis")
  ps.taxonomy.push(tax3)
  
  ps.metadata.set("project", "microbiome_study")
  ps.metadata.set("platform", "Illumina")
  
  ps
}