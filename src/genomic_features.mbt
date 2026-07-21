///|
pub struct GFGene {
  gene_id : String
  gene_name : String
  gene_biotype : String
  chromosome : String
  start : Int
  end : Int
  strand : String
  description : String
  transcripts : Array[GFTranscript]
}

///|
pub fn GFGene::new(
  gene_id : String,
  gene_name : String,
  gene_biotype : String,
  chromosome : String,
  start : Int,
  end : Int,
  strand : String,
  description : String,
) -> GFGene {
  GFGene::{
    gene_id,
    gene_name,
    gene_biotype,
    chromosome,
    start,
    end,
    strand,
    description,
    transcripts: Array::new(),
  }
}

///|
pub fn GFGene::to_string(self : GFGene) -> String {
  let mut result = self.gene_id + " (" + self.gene_name + ")\n"
  result = result + "  Biotype: " + self.gene_biotype + "\n"
  result = result +
    "  Location: " +
    self.chromosome +
    ":" +
    self.start.to_string() +
    "-" +
    self.end.to_string() +
    " " +
    self.strand +
    "\n"
  result = result +
    "  Transcripts: " +
    self.transcripts.length().to_string() +
    "\n"
  result = result + "  Description: " + self.description
  result
}

///|
pub struct GFTranscript {
  transcript_id : String
  transcript_name : String
  gene_id : String
  chromosome : String
  start : Int
  end : Int
  strand : String
  exons : Array[GFExon]
  cds : Array[GFExon]
  utr5 : Array[GFExon]
  utr3 : Array[GFExon]
}

///|
pub fn GFTranscript::new(
  transcript_id : String,
  transcript_name : String,
  gene_id : String,
  chromosome : String,
  start : Int,
  end : Int,
  strand : String,
) -> GFTranscript {
  GFTranscript::{
    transcript_id,
    transcript_name,
    gene_id,
    chromosome,
    start,
    end,
    strand,
    exons: Array::new(),
    cds: Array::new(),
    utr5: Array::new(),
    utr3: Array::new(),
  }
}

///|
pub struct GFExon {
  exon_id : String
  exon_number : Int
  chromosome : String
  start : Int
  end : Int
  strand : String
}

///|
pub fn GFExon::new(
  exon_id : String,
  exon_number : Int,
  chromosome : String,
  start : Int,
  end : Int,
  strand : String,
) -> GFExon {
  GFExon::{ exon_id, exon_number, chromosome, start, end, strand }
}

///|
pub struct GenomeAnnotation {
  genes : Array[GFGene]
  transcripts : Array[GFTranscript]
  genome : String
  version : String
  source : String
}

///|
pub fn GenomeAnnotation::new(
  genome : String,
  version : String,
  source : String,
) -> GenomeAnnotation {
  GenomeAnnotation::{
    genes: Array::new(),
    transcripts: Array::new(),
    genome,
    version,
    source,
  }
}

///|
pub fn genome_annotation_add_gene(
  annotation : GenomeAnnotation,
  gene : GFGene,
) -> Unit {
  annotation.genes.push(gene)
}

///|
pub fn genome_annotation_add_transcript(
  annotation : GenomeAnnotation,
  transcript : GFTranscript,
) -> Unit {
  annotation.transcripts.push(transcript)
}

///|
pub fn genome_annotation_get_gene_by_id(
  annotation : GenomeAnnotation,
  gene_id : String,
) -> GFGene? {
  let mut i = 0
  while i < annotation.genes.length() {
    let gene = annotation.genes[i]
    if gene.gene_id == gene_id {
      return Some(gene)
    }
    i = i + 1
  }
  None
}

///|
pub fn genome_annotation_get_gene_by_name(
  annotation : GenomeAnnotation,
  gene_name : String,
) -> Array[GFGene] {
  let results : Array[GFGene] = Array::new()
  let mut i = 0
  while i < annotation.genes.length() {
    let gene = annotation.genes[i]
    if gene.gene_name == gene_name {
      results.push(gene)
    }
    i = i + 1
  }
  results
}

///|
pub fn genome_annotation_get_transcripts_by_gene(
  annotation : GenomeAnnotation,
  gene_id : String,
) -> Array[GFTranscript] {
  let results : Array[GFTranscript] = Array::new()
  let mut i = 0
  while i < annotation.transcripts.length() {
    let transcript = annotation.transcripts[i]
    if transcript.gene_id == gene_id {
      results.push(transcript)
    }
    i = i + 1
  }
  results
}

///|
pub fn genome_annotation_genes_on_chromosome(
  annotation : GenomeAnnotation,
  chromosome : String,
) -> Array[GFGene] {
  let results : Array[GFGene] = Array::new()
  let mut i = 0
  while i < annotation.genes.length() {
    let gene = annotation.genes[i]
    if gene.chromosome == chromosome {
      results.push(gene)
    }
    i = i + 1
  }
  results
}

///|
pub fn genome_annotation_genes_in_region(
  annotation : GenomeAnnotation,
  chromosome : String,
  start : Int,
  end : Int,
) -> Array[GFGene] {
  let results : Array[GFGene] = Array::new()
  let mut i = 0
  while i < annotation.genes.length() {
    let gene = annotation.genes[i]
    if gene.chromosome == chromosome && gene.start <= end && gene.end >= start {
      results.push(gene)
    }
    i = i + 1
  }
  results
}

///|
pub fn genome_annotation_count_genes(annotation : GenomeAnnotation) -> Int {
  annotation.genes.length()
}

///|
pub fn genome_annotation_count_transcripts(
  annotation : GenomeAnnotation,
) -> Int {
  annotation.transcripts.length()
}

///|
pub fn genome_annotation_summarize(annotation : GenomeAnnotation) -> String {
  let mut result = "GenomeAnnotation Summary\n"
  result = result + "========================\n"
  result = result + "Genome: " + annotation.genome + "\n"
  result = result + "Version: " + annotation.version + "\n"
  result = result + "Source: " + annotation.source + "\n\n"
  result = result + "Statistics:\n"
  result = result +
    "  - Genes: " +
    genome_annotation_count_genes(annotation).to_string() +
    "\n"
  result = result +
    "  - Transcripts: " +
    genome_annotation_count_transcripts(annotation).to_string() +
    "\n"

  let biotype_counts : Map[String, Int] = Map([], capacity=10)
  let biotype_names : Array[String] = Array::new()
  let mut i = 0
  while i < annotation.genes.length() {
    let gene = annotation.genes[i]
    let key = gene.gene_biotype
    let current = biotype_counts.get(key)
    if current.is_some() {
      biotype_counts[key] = current.unwrap() + 1
    } else {
      biotype_counts[key] = 1
      biotype_names.push(key)
    }
    i = i + 1
  }

  result = result + "\nGene biotypes:\n"
  let mut j = 0
  while j < biotype_names.length() {
    let key = biotype_names[j]
    let value = biotype_counts.get(key).unwrap()
    result = result + "  - " + key + ": " + value.to_string() + "\n"
    j = j + 1
  }
  result
}

///|
pub fn parse_gtf_line(
  line : String,
) -> (String, String, String, Int, Int, String, Map[String, String])? {
  if line.length() == 0 || line.starts_with("#") {
    return None
  }

  let parts = split_by_char(line, 9)
  if parts.length() < 9 {
    return None
  }

  let chromosome = parts[0]
  let source = parts[1]
  let feature_type = parts[2]

  let start = match parse_double(parts[3]) {
    Some(v) => v.to_int()
    None => 0
  }

  let end = match parse_double(parts[4]) {
    Some(v) => v.to_int()
    None => 0
  }

  let strand = parts[6]
  let attributes_str = parts[8]

  let attributes : Map[String, String] = Map([], capacity=5)
  let attr_pairs = split_by_char(attributes_str, ';'.to_int())
  let mut j = 0
  while j < attr_pairs.length() {
    let pair = attr_pairs[j]
    let eq_idx = pair.find("=")
    if eq_idx.is_some() {
      let pos = eq_idx.unwrap()
      let key = trim(pair[0:pos].to_string())
      let mut value = trim(pair[pos + 1:].to_string())
      if value.length() >= 2 &&
        value.unsafe_get(0).to_int() == '\"'.to_int() &&
        value.unsafe_get(value.length() - 1).to_int() == '\"'.to_int() {
        value = value[1:value.length() - 1].to_string()
      }
      attributes[key] = value
    }
    j = j + 1
  }

  Some((chromosome, source, feature_type, start, end, strand, attributes))
}

///|
pub fn gf_parse_gtf(content : String) -> GenomeAnnotation {
  let annotation = GenomeAnnotation::new("Unknown", "Unknown", "GTF")
  let lines = split_by_char(content, '\n'.to_int())

  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    let parsed = parse_gtf_line(line)
    if parsed.is_some() {
      let (chromosome, _, feature_type, start, end, strand, attributes) = parsed.unwrap()

      if feature_type == "gene" {
        let gene_id = attributes.get("gene_id")
        let gene_name = attributes.get("gene_name")
        let gene_biotype = attributes.get("gene_biotype")
        let description = attributes.get("gene_description")

        let gene = GFGene::new(
          if gene_id.is_some() {
            gene_id.unwrap()
          } else {
            "unknown"
          },
          if gene_name.is_some() {
            gene_name.unwrap()
          } else {
            "unknown"
          },
          if gene_biotype.is_some() {
            gene_biotype.unwrap()
          } else {
            "protein_coding"
          },
          chromosome,
          start,
          end,
          strand,
          if description.is_some() {
            description.unwrap()
          } else {
            ""
          },
        )
        genome_annotation_add_gene(annotation, gene)
      } else if feature_type == "transcript" {
        let transcript_id = attributes.get("transcript_id")
        let transcript_name = attributes.get("transcript_name")
        let gene_id = attributes.get("gene_id")

        let transcript = GFTranscript::new(
          if transcript_id.is_some() {
            transcript_id.unwrap()
          } else {
            "unknown"
          },
          if transcript_name.is_some() {
            transcript_name.unwrap()
          } else {
            "unknown"
          },
          if gene_id.is_some() {
            gene_id.unwrap()
          } else {
            "unknown"
          },
          chromosome,
          start,
          end,
          strand,
        )
        genome_annotation_add_transcript(annotation, transcript)
      }
    }
    i = i + 1
  }

  annotation
}

///|
pub fn create_example_annotation() -> GenomeAnnotation {
  let annotation = GenomeAnnotation::new("Homo sapiens", "GRCh38", "Ensembl")

  let tp53 = GFGene::new(
    "ENSG00000130203", "TP53", "protein_coding", "chr17", 7661779, 7687550, "-",
    "Tumor protein p53",
  )

  let tp53_tx1 = GFTranscript::new(
    "ENST00000357654", "TP53-201", "ENSG00000130203", "chr17", 7661779, 7687550,
    "-",
  )

  let tp53_exon1 = GFExon::new(
    "ENSE00001234567", 1, "chr17", 7661779, 7661933, "-",
  )
  let tp53_exon2 = GFExon::new(
    "ENSE00001234568", 2, "chr17", 7662125, 7662220, "-",
  )
  tp53_tx1.exons.push(tp53_exon1)
  tp53_tx1.exons.push(tp53_exon2)
  tp53.transcripts.push(tp53_tx1)

  genome_annotation_add_gene(annotation, tp53)
  genome_annotation_add_transcript(annotation, tp53_tx1)

  let brca1 = GFGene::new(
    "ENSG00000012048", "BRCA1", "protein_coding", "chr17", 43044295, 43125370, "-",
    "Breast cancer type 1 susceptibility protein",
  )

  let brca1_tx1 = GFTranscript::new(
    "ENST00000357655", "BRCA1-201", "ENSG00000012048", "chr17", 43044295, 43125370,
    "-",
  )

  genome_annotation_add_gene(annotation, brca1)
  genome_annotation_add_transcript(annotation, brca1_tx1)

  let brca2 = GFGene::new(
    "ENSG00000139618", "BRCA2", "protein_coding", "chr13", 32315086, 32400266, "+",
    "Breast cancer type 2 susceptibility protein",
  )

  let brca2_tx1 = GFTranscript::new(
    "ENST00000357656", "BRCA2-201", "ENSG00000139618", "chr13", 32315086, 32400266,
    "+",
  )

  genome_annotation_add_gene(annotation, brca2)
  genome_annotation_add_transcript(annotation, brca2_tx1)

  let atm = GFGene::new(
    "ENSG00000140469", "ATM", "protein_coding", "chr11", 108187468, 108251873, "+",
    "Ataxia telangiectasia mutated",
  )

  let atm_tx1 = GFTranscript::new(
    "ENST00000357657", "ATM-201", "ENSG00000140469", "chr11", 108187468, 108251873,
    "+",
  )

  genome_annotation_add_gene(annotation, atm)
  genome_annotation_add_transcript(annotation, atm_tx1)

  let mir123 = GFGene::new(
    "ENSG00000284332", "MIR123", "miRNA", "chr1", 1000000, 1000100, "+", "MicroRNA 123",
  )

  genome_annotation_add_gene(annotation, mir123)

  annotation
}
