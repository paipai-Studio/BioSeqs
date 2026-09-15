///|
/// msigdbr - Molecular Signatures Database (MSigDB) gene set collections.
/// Inspired by Bioconductor {msigdbr} (Dolgalev I., 2022).
///
/// Provides a curated catalog of representative MSigDB gene sets (Hallmark H,
/// Curated canonical C2, Ontology C5), GMT round-tripping and gene identifier
/// mapping. Interoperates with GSEABase structures (`GmtGeneSet` /
/// `GmtGeneSetCollection`) so results feed directly into fgsea / GSEA /
/// clusterProfiler style enrichment.

///|
/// MSigDB major collection categories.
pub(all) enum MsigCategory {
  MsigHallmark
  MsigPositional
  MsigCurated
  MsigRegulatory
  MsigComputational
  MsigOntology
  MsigOncogenic
  MsigImmunologic
  MsigCellType
  MsigOther
}

///|
/// Official MSigDB collection code (H, C1..C8).
pub fn MsigCategory::code(self : MsigCategory) -> String {
  match self {
    MsigHallmark => "H"
    MsigPositional => "C1"
    MsigCurated => "C2"
    MsigRegulatory => "C3"
    MsigComputational => "C4"
    MsigOntology => "C5"
    MsigOncogenic => "C6"
    MsigImmunologic => "C7"
    MsigCellType => "C8"
    MsigOther => "other"
  }
}

///|
/// Human readable label for a collection category.
pub fn MsigCategory::label(self : MsigCategory) -> String {
  match self {
    MsigHallmark => "Hallmark gene sets"
    MsigPositional => "Positional gene sets"
    MsigCurated => "Curated gene sets"
    MsigRegulatory => "Regulatory target gene sets"
    MsigComputational => "Computational gene sets"
    MsigOntology => "Ontology gene sets"
    MsigOncogenic => "Oncogenic signatures"
    MsigImmunologic => "Immunologic signatures"
    MsigCellType => "Cell type signature gene sets"
    MsigOther => "Unclassified collection"
  }
}

// ---- internal helpers -----------------------------------------------------

///|
fn msig_trim(s : String) -> String {
  let mut start = 0
  let n = s.length()
  while start < n {
    let c = s.unsafe_get(start)
    if c == 32 || c == 9 || c == 10 || c == 13 {
      start = start + 1
    } else {
      break
    }
  }
  let mut end = n
  while end > start {
    let c = s.unsafe_get(end - 1)
    if c == 32 || c == 9 || c == 10 || c == 13 {
      end = end - 1
    } else {
      break
    }
  }
  if start == 0 && end == n {
    s
  } else {
    s[start:end].to_owned()
  }
}

///|
/// Parse a category from an MSigDB code ("H", "C2", ...) or a name.
pub fn msig_category_from_string(s : String) -> MsigCategory {
  let t = msig_trim(s).to_lower()
  if t == "h" || t == "hallmark" || t.contains("hallmark") {
    MsigHallmark
  } else if t == "c1" || t.contains("positional") {
    MsigPositional
  } else if t == "c2" || t.contains("curated") {
    MsigCurated
  } else if t == "c3" || t.contains("regulatory") || t.contains("motif") {
    MsigRegulatory
  } else if t == "c4" || t.contains("computational") {
    MsigComputational
  } else if t == "c5" || t.contains("ontology") || t.contains("go_") {
    MsigOntology
  } else if t == "c6" || t.contains("oncogenic") {
    MsigOncogenic
  } else if t == "c7" || t.contains("immunologic") {
    MsigImmunologic
  } else if t == "c8" || t.contains("cell type") || t.contains("celltype") {
    MsigCellType
  } else {
    MsigOther
  }
}

///|
/// A named MSigDB collection: one category, an optional subcategory (e.g.
/// "CP:KEGG", "GO:BP") and the gene sets it contains.
pub struct MsigCollection {
  pub category : MsigCategory
  pub subcategory : String
  pub sets : Array[GmtGeneSet]
}

///|
/// Construct a collection from parts.
pub fn MsigCollection::new(
  category : MsigCategory,
  subcategory : String,
  sets : Array[GmtGeneSet],
) -> MsigCollection {
  MsigCollection::{ category, subcategory, sets }
}

///|
/// Number of gene sets in the collection.
pub fn MsigCollection::size(self : MsigCollection) -> Int {
  self.sets.length()
}

///|
/// Convert to a GSEABase-style collection for downstream enrichment.
pub fn MsigCollection::to_gmt_collection(
  self : MsigCollection,
) -> GmtGeneSetCollection {
  GmtGeneSetCollection::new(self.sets)
}

///|
fn msig_set(
  name : String,
  description : String,
  genes : Array[String],
  collection_type : GeneSetCollectionType,
  subcategory : String,
) -> GmtGeneSet {
  let id = "\{name}@\{subcategory}"
  GmtGeneSet::with_annotation(
    name,
    description,
    genes,
    collection_type,
    "Homo sapiens",
    id,
  )
}

// ---- built-in catalog -----------------------------------------------------

///|
/// Curated Hallmark (H) gene sets with representative real members.
pub fn msig_hallmark() -> MsigCollection {
  let sets : Array[GmtGeneSet] = []
  sets.push(
    msig_set(
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "Genes regulated by NF-kB in response to TNF",
      [
        "NFKB1", "NFKB2", "RELA", "RELB", "TNF", "IL6", "CXCL8", "CCL2",
        "ICAM1", "BIRC3", "NFKBIA", "TNFAIP3", "JUN", "FOS", "ATF3", "CEBPB",
      ],
      GeneSetCollectionType::Hallmark,
      "H",
    ),
  )
  sets.push(
    msig_set(
      "HALLMARK_P53_PATHWAY",
      "Genes involved in p53 pathways and networks",
      [
        "TP53", "CDKN1A", "MDM2", "BAX", "GADD45A", "BBC3", "RRM2B", "SESN1",
        "TIGAR", "ZMAT3", "DDB2", "FAS", "TNFRSF10B",
      ],
      GeneSetCollectionType::Hallmark,
      "H",
    ),
  )
  sets.push(
    msig_set(
      "HALLMARK_APOPTOSIS",
      "Genes mediating programmed cell death",
      [
        "BAX", "BAK1", "CASP3", "CASP7", "CASP8", "CASP9", "BID", "BCL2L11",
        "PMAIP1", "BBC3", "FAS", "TNFRSF10B", "DIABLO",
      ],
      GeneSetCollectionType::Hallmark,
      "H",
    ),
  )
  sets.push(
    msig_set(
      "HALLMARK_E2F_TARGETS",
      "Genes regulated by E2F transcription factors",
      [
        "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "PCNA", "CDC6", "CDT1",
        "ORC1", "CCNA2", "CCNE1", "E2F1", "E2F2", "MCM8",
      ],
      GeneSetCollectionType::Hallmark,
      "H",
    ),
  )
  sets.push(
    msig_set(
      "HALLMARK_G2M_CHECKPOINT",
      "Genes involved in the G2/M checkpoint",
      [
        "CCNB1", "CCNB2", "CDK1", "CDC20", "CDC25C", "AURKA", "AURKB", "BUB1",
        "MAD2L1", "PLK1", "TOP2A", "CCNA2", "MKI67",
      ],
      GeneSetCollectionType::Hallmark,
      "H",
    ),
  )
  sets.push(
    msig_set(
      "HALLMARK_INFLAMMATORY_RESPONSE",
      "Genes defining the inflammatory response",
      [
        "IL6", "CXCL8", "CCL2", "IL1B", "TNF", "CXCL1", "CXCL2", "CCL20",
        "IL1A", "PTGS2", "CXCL6", "CCL5", "IL6R",
      ],
      GeneSetCollectionType::Hallmark,
      "H",
    ),
  )
  MsigCollection::new(MsigHallmark, "H", sets)
}

///|
/// Curated canonical (C2) gene sets - a compact KEGG / Reactome selection.
pub fn msig_curated_canonical() -> MsigCollection {
  let sets : Array[GmtGeneSet] = []
  sets.push(
    msig_set(
      "KEGG_CELL_CYCLE",
      "KEGG cell cycle pathway",
      [
        "CDK1", "CDK2", "CDK4", "CDK6", "CDKN1A", "CDKN1B", "CCNA2", "CCNB1",
        "CCND1", "CCNE1", "CDC25A", "CDC25B", "E2F1", "RB1", "MCM2", "PCNA",
      ],
      GeneSetCollectionType::Canonical,
      "CP:KEGG",
    ),
  )
  sets.push(
    msig_set(
      "KEGG_APOPTOSIS",
      "KEGG apoptosis pathway",
      [
        "TP53", "BAX", "BAK1", "BID", "CASP3", "CASP6", "CASP7", "CASP8",
        "CASP9", "CASP10", "CYCS", "APAF1", "BCL2", "BCL2L1", "FAS", "TNFRSF10B",
      ],
      GeneSetCollectionType::Canonical,
      "CP:KEGG",
    ),
  )
  sets.push(
    msig_set(
      "REACTOME_DNA_REPLICATION",
      "Reactome DNA replication",
      [
        "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "PCNA", "RPA1", "RPA2",
        "RPA3", "POLA1", "POLD1", "POLE", "RFC1", "RFC2", "LIG1",
      ],
      GeneSetCollectionType::Canonical,
      "CP:REACTOME",
    ),
  )
  sets.push(
    msig_set(
      "REACTOME_CELL_CYCLE",
      "Reactome cell cycle",
      [
        "CDK1", "CDK2", "CCNA2", "CCNB1", "CCNE1", "CDC20", "CDC25A", "CDC25C",
        "PLK1", "AURKA", "BUB1", "MAD2L1", "ESPL1", "ANAPC1", "SKP1",
      ],
      GeneSetCollectionType::Canonical,
      "CP:REACTOME",
    ),
  )
  MsigCollection::new(MsigCurated, "CP", sets)
}

///|
/// Ontology (C5) gene sets - a compact GO biological process selection.
pub fn msig_ontology_go_bp() -> MsigCollection {
  let sets : Array[GmtGeneSet] = []
  sets.push(
    msig_set(
      "GO_APOPTOTIC_PROCESS",
      "GO:0006915 apoptotic process",
      [
        "TP53", "BAX", "BAK1", "CASP3", "CASP8", "CASP9", "BID", "BBC3",
        "PMAIP1", "FAS", "CYCS", "APAF1",
      ],
      GeneSetCollectionType::GO_BiologicalProcess,
      "GO:BP",
    ),
  )
  sets.push(
    msig_set(
      "GO_DNA_REPAIR",
      "GO:0006281 DNA repair",
      [
        "TP53", "BRCA1", "BRCA2", "RAD51", "ATM", "ATR", "CHEK1", "CHEK2",
        "XRCC1", "PARP1", "MLH1", "MSH2", "MSH6",
      ],
      GeneSetCollectionType::GO_BiologicalProcess,
      "GO:BP",
    ),
  )
  sets.push(
    msig_set(
      "GO_IMMUNE_RESPONSE",
      "GO:0006955 immune response",
      [
        "IL6", "TNF", "IFNG", "IL1B", "CXCL8", "CCL2", "CD4", "CD8A", "CD19",
        "TLR4", "MYD88", "NFKB1", "STAT1",
      ],
      GeneSetCollectionType::GO_BiologicalProcess,
      "GO:BP",
    ),
  )
  MsigCollection::new(MsigOntology, "GO:BP", sets)
}

///|
/// The full built-in MSigDB catalog (all bundled collections).
pub fn msig_catalog() -> Array[MsigCollection] {
  [msig_hallmark(), msig_curated_canonical(), msig_ontology_go_bp()]
}

///|
/// Total number of gene sets across a catalog.
pub fn msig_catalog_size(catalog : Array[MsigCollection]) -> Int {
  let mut total = 0
  for col in catalog {
    total = total + col.size()
  }
  total
}

///|
/// Look up a single gene set by name across a catalog.
pub fn msig_lookup(
  catalog : Array[MsigCollection],
  name : String,
) -> GmtGeneSet? {
  for col in catalog {
    for s in col.sets {
      if s.name == name {
        return Some(s)
      }
    }
  }
  None
}

///|
/// Keep only gene sets whose size lies in [min_size, max_size].
pub fn msig_filter_by_size(
  col : MsigCollection,
  min_size : Int,
  max_size : Int,
) -> MsigCollection {
  let kept : Array[GmtGeneSet] = []
  for s in col.sets {
    let sz = s.size()
    if sz >= min_size && sz <= max_size {
      kept.push(s)
    }
  }
  MsigCollection::new(col.category, col.subcategory, kept)
}

// ---- GMT round-tripping ---------------------------------------------------

///|
/// Serialize a collection to GMT text (one set per line).
pub fn msig_to_gmt(col : MsigCollection) -> String {
  write_gmt(col.to_gmt_collection())
}

///|
/// Parse GMT text into a collection, tagging it with a category.
pub fn msig_from_gmt(
  content : String,
  category : MsigCategory,
) -> MsigCollection {
  let col = parse_gmt(content)
  MsigCollection::new(category, "imported", col.sets)
}

// ---- gene identifier mapping ---------------------------------------------

///|
/// Build an identifier mapping table from (from_id, to_id) pairs.
pub fn msig_id_table(pairs : Array[(String, String)]) -> Map[String, String] {
  let table : Map[String, String] = Map([], capacity=pairs.length() * 2 + 1)
  for p in pairs {
    let k = msig_trim(p.0)
    let v = msig_trim(p.1)
    table.set(k, v)
  }
  table
}

///|
/// Remap the gene identifiers of one gene set. Unmapped identifiers are kept
/// unchanged so that no membership information is silently lost.
pub fn msig_map_ids(gs : GmtGeneSet, table : Map[String, String]) -> GmtGeneSet {
  let mapped : Array[String] = []
  for g in gs.gene_ids {
    match table.get(g) {
      Some(v) => mapped.push(v)
      None => mapped.push(g)
    }
  }
  GmtGeneSet::with_annotation(
    gs.name,
    gs.description,
    mapped,
    gs.collection_type,
    gs.organism,
    gs.id,
  )
}

///|
/// Remap every gene set in a collection.
pub fn msig_map_collection(
  col : MsigCollection,
  table : Map[String, String],
) -> MsigCollection {
  let mapped : Array[GmtGeneSet] = []
  for s in col.sets {
    mapped.push(msig_map_ids(s, table))
  }
  MsigCollection::new(col.category, col.subcategory, mapped)
}

///|
/// Convenience: map gene symbols to Entrez ids using a supplied table.
pub fn msig_symbol_to_entrez(
  col : MsigCollection,
  symbol_to_entrez : Map[String, String],
) -> MsigCollection {
  msig_map_collection(col, symbol_to_entrez)
}

///|
/// Convenience: map Entrez ids back to gene symbols.
pub fn msig_entrez_to_symbol(
  col : MsigCollection,
  entrez_to_symbol : Map[String, String],
) -> MsigCollection {
  msig_map_collection(col, entrez_to_symbol)
}

// ---- collection analytics ------------------------------------------------

///|
/// How many sets each gene belongs to (useful for identifying hub genes).
pub fn msig_gene_frequency(col : MsigCollection) -> Map[String, Int] {
  let freq : Map[String, Int] = Map([], capacity=256)
  for s in col.sets {
    for g in s.gene_ids {
      let cur = freq.get(g).unwrap_or(0)
      freq.set(g, cur + 1)
    }
  }
  freq
}

///|
/// Genes shared by two gene sets.
pub fn msig_shared_genes(a : GmtGeneSet, b : GmtGeneSet) -> Array[String] {
  gene_set_intersect_arr(a, b)
}

///|
/// Jaccard similarity between two gene sets.
pub fn msig_jaccard(a : GmtGeneSet, b : GmtGeneSet) -> Double {
  gene_set_jaccard(a, b)
}

///|
/// Human readable summary of a collection.
pub fn msig_summary(col : MsigCollection) -> String {
  let buf = StringBuilder()
  buf.write_string(
    "MsigCollection[\{col.category.code()}] \{col.subcategory} - \{col.size()} sets\n",
  )
  let mut total_genes = 0
  for s in col.sets {
    total_genes = total_genes + s.size()
  }
  buf.write_string("  total gene memberships: \{total_genes}\n")
  for s in col.sets {
    buf.write_string("  - \{s.name} (\{s.size()} genes)\n")
  }
  buf.to_string()
}
