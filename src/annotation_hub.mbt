///|
pub enum ResourceType {
  GeneModel
  Variation
  Epigenome
  Sequence
  Alignment
  Expression
  Other(String)
}

pub fn ResourceType::to_string(self : ResourceType) -> String {
  match self {
    GeneModel => "GeneModel"
    Variation => "Variation"
    Epigenome => "Epigenome"
    Sequence => "Sequence"
    Alignment => "Alignment"
    Expression => "Expression"
    Other(s) => s
  }
}

pub enum DataProvider {
  Ensembl
  UCSC
  NCBI
  EBI
  ENCODE
  GEO
  Other(String)
}

pub fn DataProvider::to_string(self : DataProvider) -> String {
  match self {
    Ensembl => "Ensembl"
    UCSC => "UCSC"
    NCBI => "NCBI"
    EBI => "EBI"
    ENCODE => "ENCODE"
    GEO => "GEO"
    Other(s) => s
  }
}

pub struct HubResource {
  resource_id : String
  resource_type : ResourceType
  provider : DataProvider
  genome : String
  version : String
  description : String
  source_url : String
  file_size : Int
  cache_path : String
  metadata : Map[String, String]
}

pub fn HubResource::new(
  resource_id : String,
  resource_type : ResourceType,
  provider : DataProvider,
  genome : String,
  version : String,
  description : String,
  source_url : String,
  file_size : Int,
) -> HubResource {
  HubResource::{
    resource_id,
    resource_type,
    provider,
    genome,
    version,
    description,
    source_url,
    file_size,
    cache_path: "",
    metadata: Map([], capacity=10),
  }
}

pub fn HubResource::to_string(self : HubResource) -> String {
  let mut result = self.resource_id + " [" + self.resource_type.to_string() + "]\n";
  result = result + "  Provider: " + self.provider.to_string() + "\n";
  result = result + "  Genome: " + self.genome + "\n";
  result = result + "  Version: " + self.version + "\n";
  result = result + "  Size: " + self.file_size.to_string() + " bytes\n";
  result = result + "  Description: " + self.description;
  result
}

pub struct AnnotationHub {
  resources : Array[HubResource]
  cache_dir : String
  default_genome : String
}

pub fn AnnotationHub::new(cache_dir : String, default_genome : String) -> AnnotationHub {
  AnnotationHub::{
    resources: Array::new(),
    cache_dir,
    default_genome,
  }
}

pub fn annotation_hub_add_resource(hub : AnnotationHub, resource : HubResource) -> Unit {
  hub.resources.push(resource)
}

pub fn annotation_hub_search_by_type(hub : AnnotationHub, resource_type_str : String) -> Array[HubResource] {
  let results : Array[HubResource] = Array::new()
  let mut i = 0
  while i < hub.resources.length() {
    let r = hub.resources[i]
    if r.resource_type.to_string() == resource_type_str {
      results.push(r)
    }
    i = i + 1
  }
  results
}

pub fn annotation_hub_search_by_provider(hub : AnnotationHub, provider_str : String) -> Array[HubResource] {
  let results : Array[HubResource] = Array::new()
  let mut i = 0
  while i < hub.resources.length() {
    let r = hub.resources[i]
    if r.provider.to_string() == provider_str {
      results.push(r)
    }
    i = i + 1
  }
  results
}

pub fn annotation_hub_search_by_genome(hub : AnnotationHub, genome : String) -> Array[HubResource] {
  let results : Array[HubResource] = Array::new()
  let mut i = 0
  while i < hub.resources.length() {
    let r = hub.resources[i]
    if r.genome == genome {
      results.push(r)
    }
    i = i + 1
  }
  results
}

pub fn annotation_hub_search(hub : AnnotationHub, query : String) -> Array[HubResource] {
  let results : Array[HubResource] = Array::new()
  let query_lower = to_lowercase(query)
  let mut i = 0
  while i < hub.resources.length() {
    let r = hub.resources[i]
    if contains(to_lowercase(r.resource_id), query_lower) ||
       contains(to_lowercase(r.description), query_lower) ||
       contains(to_lowercase(r.genome), query_lower) {
      results.push(r)
    }
    i = i + 1
  }
  results
}

pub fn annotation_hub_get_resource(hub : AnnotationHub, resource_id : String) -> HubResource? {
  let mut i = 0
  while i < hub.resources.length() {
    let r = hub.resources[i]
    if r.resource_id == resource_id {
      return Some(r)
    }
    i = i + 1
  }
  None
}

pub fn annotation_hub_count(hub : AnnotationHub) -> Int {
  hub.resources.length()
}

pub fn annotation_hub_summarize(hub : AnnotationHub) -> String {
  let mut result = "AnnotationHub Summary\n"
  result = result + "=====================\n"
  result = result + "Total resources: " + annotation_hub_count(hub).to_string() + "\n"
  result = result + "Cache directory: " + hub.cache_dir + "\n"
  result = result + "Default genome: " + hub.default_genome + "\n\n"
  
  let type_counts : Map[String, Int] = Map([], capacity=7)
  let type_names : Array[String] = Array::new()
  let mut i = 0
  while i < hub.resources.length() {
    let r = hub.resources[i]
    let key = r.resource_type.to_string()
    let current = type_counts.get(key)
    if current.is_some() {
      type_counts[key] = current.unwrap() + 1
    } else {
      type_counts[key] = 1
      type_names.push(key)
    }
    i = i + 1
  }
  
  result = result + "Resource types:\n"
  let mut j = 0
  while j < type_names.length() {
    let key = type_names[j]
    let value = type_counts.get(key).unwrap()
    result = result + "  - " + key + ": " + value.to_string() + "\n"
    j = j + 1
  }
  result
}

pub fn create_example_hub() -> AnnotationHub {
  let hub = AnnotationHub::new("/tmp/annotation_hub", "Homo sapiens")
  
  annotation_hub_add_resource(hub, HubResource::new(
    "AH12345",
    ResourceType::GeneModel,
    DataProvider::Ensembl,
    "Homo sapiens",
    "GRCh38",
    "Ensembl gene models for human GRCh38",
    "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/",
    15000000,
  ))
  
  annotation_hub_add_resource(hub, HubResource::new(
    "AH12346",
    ResourceType::Variation,
    DataProvider::NCBI,
    "Homo sapiens",
    "GRCh38",
    "dbSNP variation data for human",
    "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b155/",
    500000000,
  ))
  
  annotation_hub_add_resource(hub, HubResource::new(
    "AH12347",
    ResourceType::Epigenome,
    DataProvider::ENCODE,
    "Homo sapiens",
    "GRCh38",
    "ENCODE ChIP-seq peaks for H3K4me3",
    "https://www.encodeproject.org/search/?type=Experiment&assay_title=ChIP-seq/",
    20000000,
  ))
  
  annotation_hub_add_resource(hub, HubResource::new(
    "AH12348",
    ResourceType::Sequence,
    DataProvider::UCSC,
    "Mus musculus",
    "GRCm39",
    "Mouse reference genome sequence",
    "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/",
    300000000,
  ))
  
  annotation_hub_add_resource(hub, HubResource::new(
    "AH12349",
    ResourceType::Expression,
    DataProvider::GEO,
    "Homo sapiens",
    "GPL96",
    "GEO expression data for human",
    "https://www.ncbi.nlm.nih.gov/geo/",
    100000000,
  ))
  
  annotation_hub_add_resource(hub, HubResource::new(
    "AH12350",
    ResourceType::GeneModel,
    DataProvider::Ensembl,
    "Mus musculus",
    "GRCm39",
    "Ensembl gene models for mouse GRCm39",
    "https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/",
    12000000,
  ))
  
  annotation_hub_add_resource(hub, HubResource::new(
    "AH12351",
    ResourceType::Epigenome,
    DataProvider::EBI,
    "Homo sapiens",
    "GRCh38",
    "Roadmap Epigenomics histone marks",
    "https://ftp.ebi.ac.uk/pub/databases/blueprint/data/",
    80000000,
  ))
  
  annotation_hub_add_resource(hub, HubResource::new(
    "AH12352",
    ResourceType::Alignment,
    DataProvider::ENCODE,
    "Homo sapiens",
    "GRCh38",
    "ENCODE RNA-seq alignments",
    "https://www.encodeproject.org/search/?type=Experiment&assay_title=RNA-seq/",
    500000000,
  ))
  
  hub
}
