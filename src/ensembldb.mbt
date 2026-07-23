///|
/// ensembldb - Ensembl annotation database interface.

///|
pub struct EnsDbExon {
  exon_id : String
  gene_id : String
  tx_id : String
  exon_rank : Int
  seq_start : Int
  seq_end : Int
  strand : Int
}

///|
pub fn EnsDbExon::new(exon_id : String, gene_id : String, tx_id : String, exon_rank : Int, seq_start : Int, seq_end : Int, strand : Int) -> EnsDbExon {
  EnsDbExon::{ exon_id, gene_id, tx_id, exon_rank, seq_start, seq_end, strand }
}

///|
pub struct EnsDbCDS {
  cds_id : String
  gene_id : String
  tx_id : String
  cds_rank : Int
  seq_start : Int
  seq_end : Int
  strand : Int
  codon_start : Int
}

///|
pub fn EnsDbCDS::new(cds_id : String, gene_id : String, tx_id : String, cds_rank : Int, seq_start : Int, seq_end : Int, strand : Int, codon_start : Int) -> EnsDbCDS {
  EnsDbCDS::{ cds_id, gene_id, tx_id, cds_rank, seq_start, seq_end, strand, codon_start }
}

///|
pub struct EnsDbTranscript {
  tx_id : String
  gene_id : String
  tx_name : String
  biotype : String
  seq_start : Int
  seq_end : Int
  strand : Int
  chromosome : String
}

///|
pub fn EnsDbTranscript::new(tx_id : String, gene_id : String, tx_name : String, biotype : String, seq_start : Int, seq_end : Int, strand : Int, chromosome : String) -> EnsDbTranscript {
  EnsDbTranscript::{ tx_id, gene_id, tx_name, biotype, seq_start, seq_end, strand, chromosome }
}

///|
pub struct EnsDbGene {
  gene_id : String
  gene_name : String
  biotype : String
  seq_start : Int
  seq_end : Int
  strand : Int
  chromosome : String
}

///|
pub fn EnsDbGene::new(gene_id : String, gene_name : String, biotype : String, seq_start : Int, seq_end : Int, strand : Int, chromosome : String) -> EnsDbGene {
  EnsDbGene::{ gene_id, gene_name, biotype, seq_start, seq_end, strand, chromosome }
}

///|
pub struct EnsDb {
  version : String
  species : String
  genes : Array[EnsDbGene]
  transcripts : Array[EnsDbTranscript]
  exons : Array[EnsDbExon]
  cds : Array[EnsDbCDS]
  metadata : Map[String, String]
}

///|
pub fn EnsDb::new(version : String, species : String) -> EnsDb {
  EnsDb::{
    version,
    species,
    genes: Array::new(),
    transcripts: Array::new(),
    exons: Array::new(),
    cds: Array::new(),
    metadata: Map([], capacity=0)
  }
}

///|
pub fn edb_add_gene(db : EnsDb, gene : EnsDbGene) -> EnsDb {
  db.genes.push(gene)
  db
}

///|
pub fn edb_add_transcript(db : EnsDb, tx : EnsDbTranscript) -> EnsDb {
  db.transcripts.push(tx)
  db
}

///|
pub fn edb_add_exon(db : EnsDb, exon : EnsDbExon) -> EnsDb {
  db.exons.push(exon)
  db
}

///|
pub fn edb_add_cds(db : EnsDb, cds : EnsDbCDS) -> EnsDb {
  db.cds.push(cds)
  db
}

///|
pub fn edb_set_metadata(db : EnsDb, key : String, value : String) -> EnsDb {
  db.metadata.set(key, value)
  db
}

///|
pub fn edb_num_genes(db : EnsDb) -> Int {
  db.genes.length()
}

///|
pub fn edb_num_transcripts(db : EnsDb) -> Int {
  db.transcripts.length()
}

///|
pub fn edb_num_exons(db : EnsDb) -> Int {
  db.exons.length()
}

///|
pub fn edb_num_cds(db : EnsDb) -> Int {
  db.cds.length()
}

///|
pub fn edb_get_gene_by_id(db : EnsDb, gene_id : String) -> EnsDbGene? {
  let mut i = 0
  while i < db.genes.length() {
    if db.genes[i].gene_id == gene_id {
      return Some(db.genes[i])
    }
    i = i + 1
  }
  None
}

///|
pub fn edb_get_gene_by_name(db : EnsDb, gene_name : String) -> Array[EnsDbGene] {
  let results : Array[EnsDbGene] = Array::new()
  
  let mut i = 0
  while i < db.genes.length() {
    if db.genes[i].gene_name == gene_name {
      results.push(db.genes[i])
    }
    i = i + 1
  }
  
  results
}

///|
pub fn edb_get_transcripts_by_gene(db : EnsDb, gene_id : String) -> Array[EnsDbTranscript] {
  let results : Array[EnsDbTranscript] = Array::new()
  
  let mut i = 0
  while i < db.transcripts.length() {
    if db.transcripts[i].gene_id == gene_id {
      results.push(db.transcripts[i])
    }
    i = i + 1
  }
  
  results
}

///|
pub fn edb_get_exons_by_transcript(db : EnsDb, tx_id : String) -> Array[EnsDbExon] {
  let results : Array[EnsDbExon] = Array::new()
  
  let mut i = 0
  while i < db.exons.length() {
    if db.exons[i].tx_id == tx_id {
      results.push(db.exons[i])
    }
    i = i + 1
  }
  
  results
}

///|
pub fn edb_get_cds_by_transcript(db : EnsDb, tx_id : String) -> Array[EnsDbCDS] {
  let results : Array[EnsDbCDS] = Array::new()
  
  let mut i = 0
  while i < db.cds.length() {
    if db.cds[i].tx_id == tx_id {
      results.push(db.cds[i])
    }
    i = i + 1
  }
  
  results
}

///|
pub fn edb_filter_by_chromosome(db : EnsDb, chromosome : String) -> EnsDb {
  let filtered = EnsDb::new(db.version, db.species)
  
  let mut i = 0
  while i < db.genes.length() {
    if db.genes[i].chromosome == chromosome {
      filtered.genes.push(db.genes[i])
    }
    i = i + 1
  }
  
  let mut j = 0
  while j < db.transcripts.length() {
    if db.transcripts[j].chromosome == chromosome {
      filtered.transcripts.push(db.transcripts[j])
    }
    j = j + 1
  }
  
  filtered
}

///|
pub fn edb_filter_by_biotype(db : EnsDb, biotype : String) -> EnsDb {
  let filtered = EnsDb::new(db.version, db.species)
  
  let mut i = 0
  while i < db.genes.length() {
    if db.genes[i].biotype == biotype {
      filtered.genes.push(db.genes[i])
    }
    i = i + 1
  }
  
  let mut j = 0
  while j < db.transcripts.length() {
    if db.transcripts[j].biotype == biotype {
      filtered.transcripts.push(db.transcripts[j])
    }
    j = j + 1
  }
  
  filtered
}

///|
pub fn edb_get_gene_length(db : EnsDb, gene_id : String) -> Int {
  let gene_opt = edb_get_gene_by_id(db, gene_id)
  if gene_opt is Some(_) {
    gene_opt.unwrap().seq_end - gene_opt.unwrap().seq_start + 1
  } else {
    0
  }
}

///|
pub fn edb_get_transcript_length(db : EnsDb, tx_id : String) -> Int {
  let mut i = 0
  while i < db.transcripts.length() {
    if db.transcripts[i].tx_id == tx_id {
      return db.transcripts[i].seq_end - db.transcripts[i].seq_start + 1
    }
    i = i + 1
  }
  0
}

///|
pub fn create_example_ensdb() -> EnsDb {
  let db = EnsDb::new("109", "Homo sapiens")
  let _ = db.metadata.set("source", "Ensembl")
  let _ = db.metadata.set("build", "GRCh38")
  
  let gene1 = EnsDbGene::new("ENSG000001", "ACTB", "protein_coding", 5566659, 5572915, 1, "17")
  let gene2 = EnsDbGene::new("ENSG000002", "GAPDH", "protein_coding", 653207, 674323, 1, "12")
  let gene3 = EnsDbGene::new("ENSG000003", "TP53", "protein_coding", 7661779, 7687550, -1, "17")
  db.genes.push(gene1)
  db.genes.push(gene2)
  db.genes.push(gene3)
  
  let tx1 = EnsDbTranscript::new("ENST000001", "ENSG000001", "ACTB-201", "protein_coding", 5566659, 5572915, 1, "17")
  let tx2 = EnsDbTranscript::new("ENST000002", "ENSG000001", "ACTB-202", "protein_coding", 5566659, 5572915, 1, "17")
  let tx3 = EnsDbTranscript::new("ENST000003", "ENSG000002", "GAPDH-201", "protein_coding", 653207, 674323, 1, "12")
  let tx4 = EnsDbTranscript::new("ENST000004", "ENSG000003", "TP53-201", "protein_coding", 7661779, 7687550, -1, "17")
  db.transcripts.push(tx1)
  db.transcripts.push(tx2)
  db.transcripts.push(tx3)
  db.transcripts.push(tx4)
  
  let exon1 = EnsDbExon::new("ENSE000001", "ENSG000001", "ENST000001", 1, 5566659, 5566743, 1)
  let exon2 = EnsDbExon::new("ENSE000002", "ENSG000001", "ENST000001", 2, 5567108, 5567231, 1)
  let exon3 = EnsDbExon::new("ENSE000003", "ENSG000001", "ENST000001", 3, 5567332, 5567420, 1)
  let exon4 = EnsDbExon::new("ENSE000004", "ENSG000002", "ENST000003", 1, 653207, 653322, 1)
  let exon5 = EnsDbExon::new("ENSE000005", "ENSG000002", "ENST000003", 2, 654417, 654536, 1)
  db.exons.push(exon1)
  db.exons.push(exon2)
  db.exons.push(exon3)
  db.exons.push(exon4)
  db.exons.push(exon5)
  
  let cds1 = EnsDbCDS::new("ENSC000001", "ENSG000001", "ENST000001", 1, 5566680, 5566743, 1, 1)
  let cds2 = EnsDbCDS::new("ENSC000002", "ENSG000001", "ENST000001", 2, 5567108, 5567231, 1, 1)
  let cds3 = EnsDbCDS::new("ENSC000003", "ENSG000001", "ENST000001", 3, 5567332, 5567420, 1, 1)
  let cds4 = EnsDbCDS::new("ENSC000004", "ENSG000002", "ENST000003", 1, 653247, 653322, 1, 1)
  let cds5 = EnsDbCDS::new("ENSC000005", "ENSG000002", "ENST000003", 2, 654417, 654536, 1, 1)
  db.cds.push(cds1)
  db.cds.push(cds2)
  db.cds.push(cds3)
  db.cds.push(cds4)
  db.cds.push(cds5)
  
  db
}