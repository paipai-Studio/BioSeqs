///|
pub struct Annotation {
  gene_id : String
  symbol : String
  description : String
  chromosome : String
  start : Int
  end : Int
  strand : String
  biotype : String
}

///|
pub fn Annotation::new(
  gene_id : String,
  symbol : String,
  description : String,
  chromosome : String,
  start : Int,
  end : Int,
  strand : String,
  biotype : String,
) -> Annotation {
  Annotation::{
    gene_id,
    symbol,
    description,
    chromosome,
    start,
    end,
    strand,
    biotype,
  }
}

///|
pub struct FilterCondition {
  field : String
  operator : String
  value : String
}

///|
pub fn FilterCondition::new(
  field : String,
  operator : String,
  value : String,
) -> FilterCondition {
  FilterCondition::{ field, operator, value }
}

///|
pub fn af_evaluate_condition(ann : Annotation, cond : FilterCondition) -> Bool {
  let field_value = match cond.field {
    "gene_id" => ann.gene_id
    "symbol" => ann.symbol
    "description" => ann.description
    "chromosome" => ann.chromosome
    "strand" => ann.strand
    "biotype" => ann.biotype
    "start" => ann.start.to_string()
    "end" => ann.end.to_string()
    _ => ""
  }

  match cond.operator {
    "==" => field_value == cond.value
    "!=" => field_value != cond.value
    "contains" => field_value.contains(cond.value)
    "starts_with" => field_value.has_prefix(cond.value)
    "ends_with" => field_value.has_suffix(cond.value)
    _ => false
  }
}

///|
pub fn af_filter(
  annotations : Array[Annotation],
  conditions : Array[FilterCondition],
) -> Array[Annotation] {
  let result : Array[Annotation] = Array::new()
  let mut i = 0
  while i < annotations.length() {
    let mut match_all = true
    let mut j = 0
    while j < conditions.length() && match_all {
      if !af_evaluate_condition(annotations[i], conditions[j]) {
        match_all = false
      }
      j = j + 1
    }
    if match_all {
      result.push(annotations[i])
    }
    i = i + 1
  }
  result
}

///|
pub fn af_filter_chromosome(
  annotations : Array[Annotation],
  chromosome : String,
) -> Array[Annotation] {
  af_filter(annotations, [FilterCondition::new("chromosome", "==", chromosome)])
}

///|
pub fn af_filter_biotype(
  annotations : Array[Annotation],
  biotype : String,
) -> Array[Annotation] {
  af_filter(annotations, [FilterCondition::new("biotype", "==", biotype)])
}

///|
pub fn af_filter_strand(
  annotations : Array[Annotation],
  strand : String,
) -> Array[Annotation] {
  af_filter(annotations, [FilterCondition::new("strand", "==", strand)])
}

///|
pub fn af_filter_region(
  annotations : Array[Annotation],
  chromosome : String,
  start : Int,
  end : Int,
) -> Array[Annotation] {
  let result : Array[Annotation] = Array::new()
  let mut i = 0
  while i < annotations.length() {
    let ann = annotations[i]
    if ann.chromosome == chromosome && ann.end > start && ann.start < end {
      result.push(ann)
    }
    i = i + 1
  }
  result
}

///|
pub fn af_filter_description_contains(
  annotations : Array[Annotation],
  pattern : String,
) -> Array[Annotation] {
  af_filter(annotations, [
    FilterCondition::new("description", "contains", pattern),
  ])
}

///|
pub fn af_filter_symbol_pattern(
  annotations : Array[Annotation],
  pattern : String,
) -> Array[Annotation] {
  af_filter(annotations, [FilterCondition::new("symbol", "contains", pattern)])
}

///|
pub fn af_filter_biotype_in(
  annotations : Array[Annotation],
  biotypes : Array[String],
) -> Array[Annotation] {
  let result : Array[Annotation] = Array::new()
  let mut i = 0
  while i < annotations.length() {
    let mut match_biotype = false
    let mut j = 0
    while j < biotypes.length() && !match_biotype {
      if annotations[i].biotype == biotypes[j] {
        match_biotype = true
      }
      j = j + 1
    }
    if match_biotype {
      result.push(annotations[i])
    }
    i = i + 1
  }
  result
}

///|
pub fn af_filter_by_multiple_conditions(
  annotations : Array[Annotation],
  conditions : Array[FilterCondition],
  logic : String,
) -> Array[Annotation] {
  let result : Array[Annotation] = Array::new()
  let mut i = 0
  while i < annotations.length() {
    let mut match_result = match logic {
      "AND" => true
      "OR" => false
      _ => true
    }

    let mut j = 0
    while j < conditions.length() {
      let matched = af_evaluate_condition(annotations[i], conditions[j])
      match logic {
        "AND" => match_result = match_result && matched
        "OR" => match_result = match_result || matched
        _ => match_result = match_result && matched
      }
      j = j + 1
    }

    if match_result {
      result.push(annotations[i])
    }
    i = i + 1
  }
  result
}

///|
pub fn af_count_by_biotype(annotations : Array[Annotation]) -> Map[String, Int] {
  let counts : Map[String, Int] = Map([], capacity=10)
  let mut i = 0
  while i < annotations.length() {
    let biotype = annotations[i].biotype
    let current = counts.get(biotype)
    counts[biotype] = if current.is_some() { current.unwrap() + 1 } else { 1 }
    i = i + 1
  }
  counts
}

///|
pub fn af_count_by_chromosome(
  annotations : Array[Annotation],
) -> Map[String, Int] {
  let counts : Map[String, Int] = Map([], capacity=25)
  let mut i = 0
  while i < annotations.length() {
    let chr = annotations[i].chromosome
    let current = counts.get(chr)
    counts[chr] = if current.is_some() { current.unwrap() + 1 } else { 1 }
    i = i + 1
  }
  counts
}

///|
pub fn af_filter_summary(annotations : Array[Annotation]) -> String {
  let mut result = "Annotation Filter Summary\n"
  result = result + "========================\n"
  result = result +
    "Total annotations: " +
    annotations.length().to_string() +
    "\n"

  let biotype_counts = af_count_by_biotype(annotations)
  result = result + "\nBiotype distribution:\n"
  for key in biotype_counts.keys() {
    let value = biotype_counts.get(key)
    if value.is_some() {
      result = result + "  " + key + ": " + value.unwrap().to_string() + "\n"
    }
  }

  let chr_counts = af_count_by_chromosome(annotations)
  result = result + "\nChromosome distribution:\n"
  for key in chr_counts.keys() {
    let value = chr_counts.get(key)
    if value.is_some() {
      result = result + "  " + key + ": " + value.unwrap().to_string() + "\n"
    }
  }

  result
}

///|
pub fn af_create_example_annotations() -> Array[Annotation] {
  [
    Annotation::new(
      "ENSG000001", "TP53", "Tumor protein p53", "chr17", 7661779, 7687550, "+",
      "protein_coding",
    ),
    Annotation::new(
      "ENSG000002", "BRCA1", "Breast cancer type 1", "chr17", 43044294, 43125370,
      "-", "protein_coding",
    ),
    Annotation::new(
      "ENSG000003", "BRCA2", "Breast cancer type 2", "chr13", 32315086, 32400266,
      "+", "protein_coding",
    ),
    Annotation::new(
      "ENSG000004", "ATM", "Ataxia telangiectasia", "chr11", 108176558, 108217876,
      "-", "protein_coding",
    ),
    Annotation::new(
      "ENSG000005", "CDK4", "Cyclin-dependent kinase 4", "chr12", 58147027, 58153446,
      "+", "protein_coding",
    ),
    Annotation::new(
      "ENSG000006", "MYC", "Myc proto-oncogene", "chr8", 128748310, 128750027, "+",
      "protein_coding",
    ),
    Annotation::new(
      "ENSG000007", "EGFR", "Epidermal growth factor receptor", "chr7", 55019034,
      55201588, "+", "protein_coding",
    ),
    Annotation::new(
      "ENSG000008", "ACTB", "Actin beta", "chr7", 55348172, 55356365, "+", "protein_coding",
    ),
    Annotation::new(
      "ENSG000009", "GAPDH", "Glyceraldehyde-3-phosphate", "chr12", 65363620, 65374338,
      "+", "protein_coding",
    ),
    Annotation::new(
      "ENSG000010", "MIR123", "MicroRNA 123", "chr1", 1000000, 1000022, "+", "miRNA",
    ),
    Annotation::new(
      "ENSG000011", "MIR456", "MicroRNA 456", "chr2", 2000000, 2000025, "-", "miRNA",
    ),
    Annotation::new(
      "ENSG000012", "LINC001", "Long intergenic non-protein coding", "chr3", 3000000,
      3001000, "+", "lincRNA",
    ),
  ]
}
