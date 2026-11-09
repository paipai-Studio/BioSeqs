/// Scanpy second-round: Extended gene queries
/// Inspired by scanpy.queries: gene_coordinates, gene_symbols, gene_biotype, miwi

///|
/// Gene coordinate record.
pub struct GeneCoordinate {
  chromosome : String
  start : Int32
  end : Int32
  strand : String // "+" or "-"
}

///|
/// Look up gene coordinates from a reference table.
/// In real scanpy this queries Ensembl/UCSC; here we use a local lookup.
pub fn scanpy_gene_coordinates(
  gene_names : Array[String],
  ref_table : Map[String, GeneCoordinate],
) -> Map[String, GeneCoordinate] {
  let result = Map([])
  for g in gene_names {
    let coord = ref_table.get(g)
    match coord {
      Some(c) => result[g] = c
      None => ()
    }
  }
  result
}

///|
/// Build a gene coordinate reference table from GTF-like text.
/// Each line: gene_name \t chromosome \t start \t end \t strand
pub fn scanpy_gene_coordinates_from_table(
  lines : Array[String],
) -> Map[String, GeneCoordinate] {
  let result = Map([])
  for line in lines {
    let parts = line.split("\t").to_array()
    if parts.length() >= 5 {
      let name = parts[0].to_owned()
      let chrom = parts[1].to_owned()
      let start = parse_int(parts[2].to_owned())
      let end = parse_int(parts[3].to_owned())
      let strand = parts[4].to_owned()
      result[name] = GeneCoordinate::{ chromosome: chrom, start, end, strand, }
    }
  }
  result
}

///|
/// Gene ID <-> Gene Symbol mapping.
pub struct GeneSymbolMapping {
  id_to_symbol : Map[String, String]
  symbol_to_id : Map[String, String]
}

///|
/// Build bidirectional gene ID/symbol mapping from a table.
/// Each entry: (gene_id, gene_symbol)
pub fn scanpy_gene_symbols_build(
  pairs : Array[(String, String)],
) -> GeneSymbolMapping {
  let id_to_symbol = Map([])
  let symbol_to_id = Map([])
  for p in pairs {
    id_to_symbol[p.0] = p.1
    symbol_to_id[p.1] = p.0
  }
  GeneSymbolMapping::{ id_to_symbol, symbol_to_id, }
}

///|
/// Convert gene IDs to symbols using a mapping.
pub fn scanpy_gene_ids_to_symbols(
  gene_ids : Array[String],
  mapping : GeneSymbolMapping,
) -> Array[String] {
  let result = Array(capacity=gene_ids.length())
  for id in gene_ids {
    let sym = mapping.id_to_symbol.get(id)
    match sym {
      Some(s) => result.push(s)
      None => result.push(id) // keep original if not found
    }
  }
  result
}

///|
/// Convert gene symbols to IDs using a mapping.
pub fn scanpy_gene_symbols_to_ids(
  gene_symbols : Array[String],
  mapping : GeneSymbolMapping,
) -> Array[String] {
  let result = Array(capacity=gene_symbols.length())
  for sym in gene_symbols {
    let id = mapping.symbol_to_id.get(sym)
    match id {
      Some(s) => result.push(s)
      None => result.push(sym)
    }
  }
  result
}

///|
/// Gene biotype classification.
pub enum GeneBiotype {
  ProteinCoding
  LncRNA
  MiRNA
  SnRNA
  Pseudogene
  RRNA
  TRNA
  Other(String)
}

///|
/// Classify genes by biotype based on a reference lookup.
pub fn scanpy_gene_biotype(
  gene_names : Array[String],
  biotype_table : Map[String, String],
) -> Array[(String, GeneBiotype)] {
  let result = Array(capacity=gene_names.length())
  for g in gene_names {
    let bt = biotype_table.get(g).unwrap_or("unknown")
    let biotype = match bt {
      "protein_coding" => ProteinCoding
      "lncRNA" => LncRNA
      "miRNA" => MiRNA
      "snRNA" => SnRNA
      "pseudogene" => Pseudogene
      "rRNA" => RRNA
      "tRNA" => TRNA
      _ => Other(bt)
    }
    result.push((g, biotype))
  }
  result
}

///|
/// Filter genes by biotype.
pub fn scanpy_filter_by_biotype(
  gene_names : Array[String],
  biotypes : Array[(String, GeneBiotype)],
  target~ : GeneBiotype = ProteinCoding,
) -> Array[String] {
  let result = Array(capacity=gene_names.length())
  for entry in biotypes {
    let match_target = match target {
      ProteinCoding =>
        match entry.1 {
          ProteinCoding => true
          _ => false
        }
      LncRNA =>
        match entry.1 {
          LncRNA => true
          _ => false
        }
      MiRNA =>
        match entry.1 {
          MiRNA => true
          _ => false
        }
      SnRNA =>
        match entry.1 {
          SnRNA => true
          _ => false
        }
      Pseudogene =>
        match entry.1 {
          Pseudogene => true
          _ => false
        }
      RRNA =>
        match entry.1 {
          RRNA => true
          _ => false
        }
      TRNA =>
        match entry.1 {
          TRNA => true
          _ => false
        }
      Other(_) => false
    }
    if match_target {
      result.push(entry.0)
    }
  }
  ignore(gene_names)
  result
}

///|
/// Query microRNA (miwi) genes from gene names.
/// Identifies miRNA genes by prefix pattern.
pub fn scanpy_miwi_genes(gene_names : Array[String]) -> Array[String] {
  let result = Array(capacity=0)
  for g in gene_names {
    let upper = g.to_upper()
    if upper.has_prefix("MIR") || upper.has_prefix("MIR-") {
      result.push(g)
    }
  }
  result
}

// Dummy functions to suppress unused constructor warnings

///|
pub fn _scanpy_queries2_ignore_unused() -> Unit {
  ignore(ProteinCoding)
  ignore(LncRNA)
  ignore(MiRNA)
  ignore(SnRNA)
  ignore(Pseudogene)
  ignore(RRNA)
  ignore(TRNA)
  ignore(GeneBiotype::Other("x"))
}
