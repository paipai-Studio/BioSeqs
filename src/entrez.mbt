///|
/// Bio.Entrez - NCBI Entrez database access (simulated).
/// Inspired by Biopython's Bio.Entrez module.

///|
/// Entrez database record type.
pub enum EntrezDB {
  Pubmed // PubMed citations
  Nucleotide // Nucleotide sequences
  Protein // Protein sequences
  Gene // Gene information
  SNP // Single nucleotide polymorphisms
  Taxonomy // Taxonomy data
} derive(Eq, Show)

///|
/// Construct Pubmed EntrezDB.
pub fn entrez_db_pubmed() -> EntrezDB {
  Pubmed
}

///|
/// Construct Nucleotide EntrezDB.
pub fn entrez_db_nucleotide() -> EntrezDB {
  Nucleotide
}

///|
/// Construct Protein EntrezDB.
pub fn entrez_db_protein() -> EntrezDB {
  Protein
}

///|
/// Construct Gene EntrezDB.
pub fn entrez_db_gene() -> EntrezDB {
  Gene
}

///|
/// Construct SNP EntrezDB.
pub fn entrez_db_snp() -> EntrezDB {
  SNP
}

///|
/// Construct Taxonomy EntrezDB.
pub fn entrez_db_taxonomy() -> EntrezDB {
  Taxonomy
}

///|
/// Convert EntrezDB to string name.
pub fn EntrezDB::to_string(self : EntrezDB) -> String {
  match self {
    Pubmed => "pubmed"
    Nucleotide => "nuccore"
    Protein => "protein"
    Gene => "gene"
    SNP => "snp"
    Taxonomy => "taxonomy"
  }
}

///|
/// Construct EntrezDB from string name.
pub fn EntrezDB::from_string(name : String) -> EntrezDB {
  let lower = name
  if lower == "pubmed" || lower == "PUBMED" || lower == "Pubmed" {
    Pubmed
  } else if lower == "nucleotide" ||
    lower == "NUCLEOTIDE" ||
    lower == "Nucleotide" {
    Nucleotide
  } else if lower == "protein" || lower == "PROTEIN" || lower == "Protein" {
    Protein
  } else if lower == "gene" || lower == "GENE" || lower == "Gene" {
    Gene
  } else if lower == "snp" || lower == "SNP" {
    SNP
  } else if lower == "taxonomy" || lower == "TAXONOMY" || lower == "Taxonomy" {
    Taxonomy
  } else {
    Pubmed
  }
}

///|
/// Search result item.
pub struct SearchResult {
  id : String
  title : String
  score : Double
}

///|
/// Construct SearchResult.
pub fn SearchResult::new(
  id~ : String,
  title~ : String,
  score~ : Double,
) -> SearchResult {
  SearchResult::{ id, title, score }
}

///|
/// PubMed article record.
pub struct PubmedArticle {
  pmid : String
  title : String
  abstract : String
  authors : Array[String]
  journal : String
  publication_date : String
  keywords : Array[String]
}

///|
/// Construct PubmedArticle.
pub fn PubmedArticle::new(
  pmid~ : String,
  title~ : String,
  abstract~ : String,
  authors~ : Array[String],
  journal~ : String,
  publication_date~ : String,
  keywords~ : Array[String],
) -> PubmedArticle {
  PubmedArticle::{
    pmid,
    title,
    abstract,
    authors,
    journal,
    publication_date,
    keywords,
  }
}

///|
/// Gene record.
pub struct GeneRecord {
  gene_id : String
  symbol : String
  description : String
  organism : String
  chromosome : String
  start : Int
  end : Int
  strand : String
}

///|
/// Construct GeneRecord.
pub fn GeneRecord::new(
  gene_id~ : String,
  symbol~ : String,
  description~ : String,
  organism~ : String,
  chromosome~ : String,
  start~ : Int,
  end~ : Int,
  strand~ : String,
) -> GeneRecord {
  GeneRecord::{
    gene_id,
    symbol,
    description,
    organism,
    chromosome,
    start,
    end,
    strand,
  }
}

///|
/// Taxonomy record.
pub struct TaxonomyRecord {
  tax_id : String
  scientific_name : String
  common_name : String
  rank : String
  lineage : Array[String]
  genetic_code : Int
}

///|
/// Construct TaxonomyRecord.
pub fn TaxonomyRecord::new(
  tax_id~ : String,
  scientific_name~ : String,
  common_name~ : String,
  rank~ : String,
  lineage~ : Array[String],
  genetic_code~ : Int,
) -> TaxonomyRecord {
  TaxonomyRecord::{
    tax_id,
    scientific_name,
    common_name,
    rank,
    lineage,
    genetic_code,
  }
}

///|
/// Simulated ESearch: Search the Entrez databases.
pub fn esearch(
  db : EntrezDB,
  term : String,
  retmax? : Int = 10,
) -> Array[SearchResult] {
  let results : Array[SearchResult] = Array::new()
  match db {
    Pubmed =>
      if term.contains("CRISPR") || term.contains("crispr") {
        results.push(
          SearchResult::new(
            id="37566892",
            title="CRISPR-Cas9 gene editing in human cells",
            score=98.5,
          ),
        )
        results.push(
          SearchResult::new(
            id="36421234",
            title="Applications of CRISPR technology in medicine",
            score=92.3,
          ),
        )
        results.push(
          SearchResult::new(
            id="35123456",
            title="CRISPR-based diagnostics for viral infections",
            score=87.1,
          ),
        )
      } else if term.contains("COVID") || term.contains("covid") {
        results.push(
          SearchResult::new(
            id="38123456",
            title="SARS-CoV-2 genome sequencing and analysis",
            score=95.7,
          ),
        )
        results.push(
          SearchResult::new(
            id="37234567",
            title="COVID-19 vaccine development strategies",
            score=91.2,
          ),
        )
      } else if term.contains("cancer") || term.contains("Cancer") {
        results.push(
          SearchResult::new(
            id="36890123",
            title="Novel biomarkers for cancer diagnosis",
            score=89.4,
          ),
        )
        results.push(
          SearchResult::new(
            id="35678901",
            title="Immunotherapy approaches in cancer treatment",
            score=86.8,
          ),
        )
      } else {
        results.push(
          SearchResult::new(
            id="12345678",
            title="General biology research article",
            score=50.0,
          ),
        )
      }
    Gene =>
      if term.contains("TP53") || term.contains("p53") {
        results.push(
          SearchResult::new(
            id="7157",
            title="TP53 - tumor protein p53",
            score=99.0,
          ),
        )
        results.push(
          SearchResult::new(
            id="7161",
            title="TP53BP1 - TP53 binding protein 1",
            score=85.0,
          ),
        )
      } else if term.contains("BRCA1") {
        results.push(
          SearchResult::new(
            id="672",
            title="BRCA1 - breast cancer 1",
            score=98.5,
          ),
        )
      } else {
        results.push(
          SearchResult::new(id="1", title="Gene search result", score=50.0),
        )
      }
    Nucleotide =>
      if term.contains("MT-ND1") || term.contains("mitochondrial") {
        results.push(
          SearchResult::new(
            id="NC_012920.1",
            title="Homo sapiens mitochondrion",
            score=99.5,
          ),
        )
        results.push(
          SearchResult::new(
            id="NC_001607.1",
            title="Mus musculus mitochondrion",
            score=95.0,
          ),
        )
      } else {
        results.push(
          SearchResult::new(
            id="NM_000546.6",
            title="TP53 mRNA sequence",
            score=88.0,
          ),
        )
      }
    Protein =>
      if term.contains("P53") || term.contains("p53") {
        results.push(
          SearchResult::new(id="P04637", title="Tumor protein p53", score=99.0),
        )
      } else {
        results.push(
          SearchResult::new(id="P03372", title="BRCA1 protein", score=90.0),
        )
      }
    Taxonomy =>
      if term.contains("Homo sapiens") || term.contains("human") {
        results.push(
          SearchResult::new(id="9606", title="Homo sapiens", score=99.9),
        )
      } else if term.contains("Mus musculus") || term.contains("mouse") {
        results.push(
          SearchResult::new(id="10090", title="Mus musculus", score=99.5),
        )
      } else {
        results.push(SearchResult::new(id="2", title="Bacteria", score=70.0))
      }
    SNP =>
      results.push(
        SearchResult::new(
          id="rs123456",
          title="Common SNP in human genome",
          score=80.0,
        ),
      )
  }
  if results.length() > retmax {
    results[0:retmax].to_owned()
  } else {
    results
  }
}

///|
/// Simulated EFetch: Retrieve full records from Entrez.
pub fn efetch(db : EntrezDB, id : String) -> String {
  match db {
    Pubmed =>
      match id {
        "37566892" =>
          "<PubmedArticle><PMID>37566892</PMID><ArticleTitle>CRISPR-Cas9 gene editing in human cells</ArticleTitle><Abstract>CRISPR-Cas9 is a revolutionary gene editing technology that allows precise modification of genomic DNA. This study demonstrates efficient gene knockout and knock-in in various human cell lines.</Abstract><AuthorList><Author>Zhang F</Author><Author>Li W</Author></AuthorList><Journal>Nature Biotechnology</Journal><PublicationDate>2023-05-15</PublicationDate></PubmedArticle>"
        "36421234" =>
          "<PubmedArticle><PMID>36421234</PMID><ArticleTitle>Applications of CRISPR technology in medicine</ArticleTitle><Abstract>Review of current applications of CRISPR-Cas systems in therapeutic development.</Abstract><Journal>Cell</Journal><PublicationDate>2023-01-20</PublicationDate></PubmedArticle>"
        _ =>
          "<PubmedArticle><PMID>" +
          id +
          "</PMID><ArticleTitle>Article not found</ArticleTitle></PubmedArticle>"
      }
    Gene =>
      match id {
        "7157" =>
          "<Gene><GeneID>7157</GeneID><Symbol>TP53</Symbol><Description>Tumor protein p53</Description><Organism>Homo sapiens</Organism><Chromosome>17</Chromosome><Start>7661779</Start><End>7687550</End><Strand>-1</Strand></Gene>"
        "672" =>
          "<Gene><GeneID>672</GeneID><Symbol>BRCA1</Symbol><Description>Breast cancer type 1 susceptibility protein</Description><Organism>Homo sapiens</Organism><Chromosome>17</Chromosome><Start>43044295</Start><End>43125370</End><Strand>-1</Strand></Gene>"
        _ =>
          "<Gene><GeneID>" +
          id +
          "</GeneID><Symbol>Gene not found</Symbol></Gene>"
      }
    Taxonomy =>
      match id {
        "9606" =>
          "<Taxonomy><TaxId>9606</TaxId><ScientificName>Homo sapiens</ScientificName><CommonName>human</CommonName><Rank>species</Rank><Lineage>cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; Dipnotetrapodomorpha; Tetrapoda; Amniota; Mammalia; Theria; Eutheria; Boreoeutheria; Euarchontoglires; Primates; Haplorrhini; Simiiformes; Catarrhini; Hominoidea; Hominidae; Homininae; Homo</Lineage><GeneticCode>1</GeneticCode></Taxonomy>"
        "10090" =>
          "<Taxonomy><TaxId>10090</TaxId><ScientificName>Mus musculus</ScientificName><CommonName>house mouse</CommonName><Rank>species</Rank><Lineage>cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; Dipnotetrapodomorpha; Tetrapoda; Amniota; Mammalia; Theria; Eutheria; Boreoeutheria; Euarchontoglires; Glires; Rodentia; Sciurognathi; Muroidea; Muridae; Murinae; Mus</Lineage><GeneticCode>1</GeneticCode></Taxonomy>"
        _ =>
          "<Taxonomy><TaxId>" +
          id +
          "</TaxId><ScientificName>Not found</ScientificName></Taxonomy>"
      }
    Nucleotide =>
      "<GBSeq><GBSeq_locus>NM_" +
      id +
      "</GBSeq_locus><GBSeq_sequence>ATGCGTACGTACGTACGT</GBSeq_sequence><GBSeq_definition>Nucleotide sequence</GBSeq_definition></GBSeq>"
    Protein =>
      "<GBSeq><GBSeq_locus>NP_" +
      id +
      "</GBSeq_locus><GBSeq_sequence>MRSLLLRFLFLLLPPLP</GBSeq_sequence><GBSeq_definition>Protein sequence</GBSeq_definition></GBSeq>"
    SNP => "<SNP><RsId>" + id + "</RsId><SNPType>single</SNPType></SNP>"
  }
}

///|
/// Parse a PubMed article from XML string.
pub fn parse_pubmed(xml : String) -> PubmedArticle? {
  let pmid = extract_tag(xml, "PMID")
  let title = extract_tag(xml, "ArticleTitle")
  let abstract_text = extract_tag(xml, "Abstract")
  let authors = extract_tags(xml, "Author")
  let journal = extract_tag(xml, "Journal")
  let pub_date = extract_tag(xml, "PublicationDate")
  let keywords = []
  match (pmid, title) {
    (Some(p), Some(t)) =>
      Some(
        PubmedArticle::new(
          pmid=p,
          title=t,
          abstract=abstract_text.unwrap_or(""),
          authors~,
          journal=journal.unwrap_or(""),
          publication_date=pub_date.unwrap_or(""),
          keywords~,
        ),
      )
    _ => None
  }
}

///|
/// Parse a gene record from XML string.
pub fn parse_gene(xml : String) -> GeneRecord? {
  let gene_id = extract_tag(xml, "GeneID")
  let symbol = extract_tag(xml, "Symbol")
  let description = extract_tag(xml, "Description")
  let organism = extract_tag(xml, "Organism")
  let chromosome = extract_tag(xml, "Chromosome")
  let start_str = extract_tag(xml, "Start")
  let end_str = extract_tag(xml, "End")
  let strand = extract_tag(xml, "Strand")
  match (gene_id, symbol) {
    (Some(g), Some(s)) =>
      Some(
        GeneRecord::new(
          gene_id=g,
          symbol=s,
          description=description.unwrap_or(""),
          organism=organism.unwrap_or(""),
          chromosome=chromosome.unwrap_or(""),
          start=parse_int(start_str.unwrap_or("0")),
          end=parse_int(end_str.unwrap_or("0")),
          strand=strand.unwrap_or(""),
        ),
      )
    _ => None
  }
}

///|
/// Parse a taxonomy record from XML string.
pub fn parse_taxonomy(xml : String) -> TaxonomyRecord? {
  let tax_id = extract_tag(xml, "TaxId")
  let scientific_name = extract_tag(xml, "ScientificName")
  let common_name = extract_tag(xml, "CommonName")
  let rank = extract_tag(xml, "Rank")
  let lineage_str = extract_tag(xml, "Lineage")
  let genetic_code_str = extract_tag(xml, "GeneticCode")
  let lineage = match lineage_str {
    Some(s) => split_by_char(s, ';'.to_int())
    None => []
  }
  match (tax_id, scientific_name) {
    (Some(t), Some(s)) =>
      Some(
        TaxonomyRecord::new(
          tax_id=t,
          scientific_name=s,
          common_name=common_name.unwrap_or(""),
          rank=rank.unwrap_or(""),
          lineage~,
          genetic_code=parse_int(genetic_code_str.unwrap_or("1")),
        ),
      )
    _ => None
  }
}

///|
/// Extract the content of a single XML tag.
fn extract_tag(xml : String, tag : String) -> String? {
  let start_tag = "<" + tag + ">"
  let end_tag = "</" + tag + ">"
  let start_pos = xml.find(start_tag)
  if start_pos is None {
    return None
  }
  let start = start_pos.unwrap() + start_tag.length()
  let remaining = xml[start:].to_owned()
  let end_pos = remaining.find(end_tag)
  if end_pos is None {
    return None
  }
  Some(remaining[0:end_pos.unwrap()].to_owned())
}

///|
/// Extract the content of multiple XML tags.
fn extract_tags(xml : String, tag : String) -> Array[String] {
  let result : Array[String] = Array::new()
  let start_tag = "<" + tag + ">"
  let end_tag = "</" + tag + ">"
  let mut remaining = xml
  while remaining.length() > 0 {
    let start_pos = remaining.find(start_tag)
    if start_pos is None {
      break
    }
    let start = start_pos.unwrap() + start_tag.length()
    let after_start = remaining[start:].to_owned()
    let end_pos = after_start.find(end_tag)
    if end_pos is None {
      break
    }
    let end = end_pos.unwrap()
    result.push(after_start[0:end].to_owned())
    remaining = after_start[end + end_tag.length():].to_owned()
  }
  result
}

///|
/// Simulated EInfo: Get information about available databases.
pub fn einfo() -> Map[String, String] {
  let result : Map[String, String] = Map([], capacity=10)
  result["pubmed"] = "PubMed bibliographic records"
  result["nuccore"] = "Nucleotide sequence records"
  result["protein"] = "Protein sequence records"
  result["gene"] = "Gene records"
  result["snp"] = "Single nucleotide polymorphisms"
  result["taxonomy"] = "Taxonomy records"
  result
}

///|
/// Simulated EGQuery: Global query for multiple databases.
pub fn egquery(term : String) -> Array[(String, Int)] {
  let results : Array[(String, Int)] = Array::new()
  let mut count = 100
  if term.contains("CRISPR") || term.contains("crispr") {
    count = 5000
  } else if term.contains("cancer") || term.contains("Cancer") {
    count = 15000
  }
  results.push(("pubmed", count))
  results.push(("nuccore", count * 2))
  results.push(("protein", count * 3))
  results.push(("gene", count / 10))
  results.push(("snp", count / 5))
  results.push(("taxonomy", count / 100))
  results
}

///|
/// Get gene information from NCBI (simulated).
pub fn get_gene_info(
  gene_symbol : String,
  organism? : String = "Homo sapiens",
) -> GeneRecord? {
  let lower = gene_symbol
  if lower == "tp53" || lower == "p53" || lower == "TP53" || lower == "P53" {
    Some(
      GeneRecord::new(
        gene_id="7157",
        symbol="TP53",
        description="Tumor protein p53",
        organism="Homo sapiens",
        chromosome="17",
        start=7661779,
        end=7687550,
        strand="-1",
      ),
    )
  } else if lower == "brca1" || lower == "BRCA1" {
    Some(
      GeneRecord::new(
        gene_id="672",
        symbol="BRCA1",
        description="Breast cancer type 1 susceptibility protein",
        organism="Homo sapiens",
        chromosome="17",
        start=43044295,
        end=43125370,
        strand="-1",
      ),
    )
  } else if lower == "egfr" || lower == "EGFR" {
    Some(
      GeneRecord::new(
        gene_id="1956",
        symbol="EGFR",
        description="Epidermal growth factor receptor",
        organism="Homo sapiens",
        chromosome="7",
        start=55019019,
        end=55211670,
        strand="-1",
      ),
    )
  } else {
    None
  }
}

///|
/// Get taxonomy information for an organism (simulated).
pub fn get_taxonomy(scientific_name : String) -> TaxonomyRecord? {
  let lower = scientific_name
  if lower.contains("homo sapiens") ||
    lower.contains("Homo sapiens") ||
    lower == "human" ||
    lower == "Human" {
    Some(
      TaxonomyRecord::new(
        tax_id="9606",
        scientific_name="Homo sapiens",
        common_name="human",
        rank="species",
        lineage=[
          "cellular organisms", "Eukaryota", "Opisthokonta", "Metazoa", "Chordata",
          "Craniata", "Vertebrata", "Mammalia", "Primates", "Homo",
        ],
        genetic_code=1,
      ),
    )
  } else if lower.contains("mus musculus") ||
    lower.contains("Mus musculus") ||
    lower == "mouse" ||
    lower == "Mouse" {
    Some(
      TaxonomyRecord::new(
        tax_id="10090",
        scientific_name="Mus musculus",
        common_name="house mouse",
        rank="species",
        lineage=[
          "cellular organisms", "Eukaryota", "Opisthokonta", "Metazoa", "Chordata",
          "Craniata", "Vertebrata", "Mammalia", "Rodentia", "Mus",
        ],
        genetic_code=1,
      ),
    )
  } else {
    None
  }
}
