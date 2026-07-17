///|
/// GenomeInfoDb - Genome information management inspired by Bioconductor GenomeInfoDb.
/// Manages chromosome names, lengths, centromere positions, and genome build versions.

///|
/// Chromosome information.
pub struct ChromInfo {
  name : String
  length : Int
  is_circular : Bool
  genome : String
}

///|
/// Construct ChromInfo.
pub fn ChromInfo::new(
  name~ : String,
  length~ : Int,
  is_circular? : Bool = false,
  genome? : String = "",
) -> ChromInfo {
  ChromInfo::{ name, length, is_circular, genome }
}

///|
/// Centromere position information.
pub struct Centromere {
  chrom : String
  start : Int
  end : Int
}

///|
/// Construct Centromere.
pub fn Centromere::new(chrom~ : String, start~ : Int, end~ : Int) -> Centromere {
  Centromere::{ chrom, start, end }
}

///|
/// Genome information database.
pub struct GenomeInfoDb {
  genome_build : String
  organism : String
  chrom_infos : Array[ChromInfo]
  centromeres : Array[Centromere]
}

///|
/// Construct GenomeInfoDb.
pub fn GenomeInfoDb::new(
  genome_build~ : String,
  organism? : String = "",
) -> GenomeInfoDb {
  GenomeInfoDb::{ genome_build, organism, chrom_infos: [], centromeres: [] }
}

///|
/// Add chromosome information.
pub fn GenomeInfoDb::add_chrom(
  self : GenomeInfoDb,
  name : String,
  length : Int,
  is_circular? : Bool = false,
) -> GenomeInfoDb {
  let new_info = ChromInfo::new(
    name~,
    length~,
    is_circular~,
    genome=self.genome_build,
  )
  GenomeInfoDb::{
    genome_build: self.genome_build,
    organism: self.organism,
    chrom_infos: self.chrom_infos + [new_info],
    centromeres: self.centromeres,
  }
}

///|
/// Add centromere position.
pub fn GenomeInfoDb::add_centromere(
  self : GenomeInfoDb,
  chrom : String,
  start : Int,
  end : Int,
) -> GenomeInfoDb {
  let c = Centromere::new(chrom~, start~, end~)
  GenomeInfoDb::{
    genome_build: self.genome_build,
    organism: self.organism,
    chrom_infos: self.chrom_infos,
    centromeres: self.centromeres + [c],
  }
}

///|
/// Get all chromosome names.
pub fn GenomeInfoDb::seqnames(self : GenomeInfoDb) -> Array[String] {
  self.chrom_infos.map(fn(ci) { ci.name })
}

///|
/// Get total genome length.
pub fn GenomeInfoDb::genome_length(self : GenomeInfoDb) -> Int {
  let mut total = 0
  for ci in self.chrom_infos {
    total = total + ci.length
  }
  total
}

///|
/// Get chromosome length by name.
pub fn GenomeInfoDb::seqLength(self : GenomeInfoDb, chrom : String) -> Int? {
  let mut i = 0
  while i < self.chrom_infos.length() {
    if self.chrom_infos[i].name == chrom {
      return Some(self.chrom_infos[i].length)
    }
    i = i + 1
  }
  None
}

///|
/// Check if chromosome is circular.
pub fn GenomeInfoDb::is_circular(self : GenomeInfoDb, chrom : String) -> Bool? {
  let mut i = 0
  while i < self.chrom_infos.length() {
    if self.chrom_infos[i].name == chrom {
      return Some(self.chrom_infos[i].is_circular)
    }
    i = i + 1
  }
  None
}

///|
/// Get centromere position for a chromosome.
pub fn GenomeInfoDb::centromere(
  self : GenomeInfoDb,
  chrom : String,
) -> Centromere? {
  let mut i = 0
  while i < self.centromeres.length() {
    if self.centromeres[i].chrom == chrom {
      return Some(self.centromeres[i])
    }
    i = i + 1
  }
  None
}

///|
/// Get standard chromosomes (excluding random, unmapped, alt, etc.).
pub fn GenomeInfoDb::standard_chromosomes(self : GenomeInfoDb) -> Array[String] {
  self.chrom_infos
  .filter(fn(ci) {
    let name = ci.name
    if name.contains("_") ||
      name.contains("random") ||
      name.contains("Un") ||
      name.contains("alt") {
      false
    } else {
      true
    }
  })
  .map(fn(ci) { ci.name })
}

///|
/// Get non-standard chromosomes.
pub fn GenomeInfoDb::non_standard_chromosomes(
  self : GenomeInfoDb,
) -> Array[String] {
  self.chrom_infos
  .filter(fn(ci) {
    let name = ci.name
    if name.contains("_") ||
      name.contains("random") ||
      name.contains("Un") ||
      name.contains("alt") {
      true
    } else {
      false
    }
  })
  .map(fn(ci) { ci.name })
}

///|
/// Get number of chromosomes.
pub fn GenomeInfoDb::num_chromosomes(self : GenomeInfoDb) -> Int {
  self.chrom_infos.length()
}

///|
/// Create a pre-built human genome (hg38).
pub fn genome_hg38() -> GenomeInfoDb {
  GenomeInfoDb::{
    genome_build: "hg38",
    organism: "Homo sapiens",
    chrom_infos: [
      ChromInfo::new(name="chr1", length=248956422, genome="hg38"),
      ChromInfo::new(name="chr2", length=242193529, genome="hg38"),
      ChromInfo::new(name="chr3", length=198295559, genome="hg38"),
      ChromInfo::new(name="chr4", length=190214555, genome="hg38"),
      ChromInfo::new(name="chr5", length=181538259, genome="hg38"),
      ChromInfo::new(name="chr6", length=170805979, genome="hg38"),
      ChromInfo::new(name="chr7", length=159345973, genome="hg38"),
      ChromInfo::new(name="chr8", length=145138636, genome="hg38"),
      ChromInfo::new(name="chr9", length=138394717, genome="hg38"),
      ChromInfo::new(name="chr10", length=133797422, genome="hg38"),
      ChromInfo::new(name="chr11", length=135086622, genome="hg38"),
      ChromInfo::new(name="chr12", length=133275309, genome="hg38"),
      ChromInfo::new(name="chr13", length=114364328, genome="hg38"),
      ChromInfo::new(name="chr14", length=107043718, genome="hg38"),
      ChromInfo::new(name="chr15", length=101991189, genome="hg38"),
      ChromInfo::new(name="chr16", length=90338345, genome="hg38"),
      ChromInfo::new(name="chr17", length=83257441, genome="hg38"),
      ChromInfo::new(name="chr18", length=80373285, genome="hg38"),
      ChromInfo::new(name="chr19", length=58617616, genome="hg38"),
      ChromInfo::new(name="chr20", length=64444167, genome="hg38"),
      ChromInfo::new(name="chr21", length=46709983, genome="hg38"),
      ChromInfo::new(name="chr22", length=50818468, genome="hg38"),
      ChromInfo::new(name="chrX", length=156040895, genome="hg38"),
      ChromInfo::new(name="chrY", length=57227415, genome="hg38"),
      ChromInfo::new(name="chrM", length=16569, is_circular=true, genome="hg38"),
    ],
    centromeres: [
      Centromere::new(chrom="chr1", start=121535434, end=124595842),
      Centromere::new(chrom="chr2", start=91191172, end=94733720),
      Centromere::new(chrom="chr3", start=87589075, end=91540110),
      Centromere::new(chrom="chr4", start=48716757, end=52716757),
      Centromere::new(chrom="chr5", start=46142139, end=50142139),
      Centromere::new(chrom="chr6", start=57809333, end=61809333),
      Centromere::new(chrom="chr7", start=54495348, end=58495348),
      Centromere::new(chrom="chr8", start=44037835, end=48037835),
      Centromere::new(chrom="chr9", start=42000000, end=46000000),
      Centromere::new(chrom="chr10", start=38215317, end=42215317),
      Centromere::new(chrom="chr11", start=51000000, end=55000000),
      Centromere::new(chrom="chr12", start=33272291, end=37272291),
      Centromere::new(chrom="chr13", start=13553000, end=17553000),
      Centromere::new(chrom="chr14", start=13553000, end=17553000),
      Centromere::new(chrom="chr15", start=13553000, end=17553000),
      Centromere::new(chrom="chr16", start=34898202, end=38898202),
      Centromere::new(chrom="chr17", start=10672816, end=14672816),
      Centromere::new(chrom="chr18", start=11700000, end=15700000),
      Centromere::new(chrom="chr19", start=14500000, end=18500000),
      Centromere::new(chrom="chr20", start=25300000, end=29300000),
      Centromere::new(chrom="chr21", start=6900000, end=10900000),
      Centromere::new(chrom="chr22", start=10500000, end=14500000),
      Centromere::new(chrom="chrX", start=56000000, end=60000000),
      Centromere::new(chrom="chrY", start=10200000, end=14200000),
    ],
  }
}

///|
/// Create a pre-built human genome (hg19).
pub fn genome_hg19() -> GenomeInfoDb {
  GenomeInfoDb::{
    genome_build: "hg19",
    organism: "Homo sapiens",
    chrom_infos: [
      ChromInfo::new(name="chr1", length=249250621, genome="hg19"),
      ChromInfo::new(name="chr2", length=243199373, genome="hg19"),
      ChromInfo::new(name="chr3", length=198022430, genome="hg19"),
      ChromInfo::new(name="chr4", length=191154276, genome="hg19"),
      ChromInfo::new(name="chr5", length=180915260, genome="hg19"),
      ChromInfo::new(name="chr6", length=171115067, genome="hg19"),
      ChromInfo::new(name="chr7", length=159138663, genome="hg19"),
      ChromInfo::new(name="chr8", length=146364022, genome="hg19"),
      ChromInfo::new(name="chr9", length=141213431, genome="hg19"),
      ChromInfo::new(name="chr10", length=135534747, genome="hg19"),
      ChromInfo::new(name="chr11", length=135006516, genome="hg19"),
      ChromInfo::new(name="chr12", length=133851895, genome="hg19"),
      ChromInfo::new(name="chr13", length=115169878, genome="hg19"),
      ChromInfo::new(name="chr14", length=107349540, genome="hg19"),
      ChromInfo::new(name="chr15", length=102531392, genome="hg19"),
      ChromInfo::new(name="chr16", length=90354753, genome="hg19"),
      ChromInfo::new(name="chr17", length=81195210, genome="hg19"),
      ChromInfo::new(name="chr18", length=78077248, genome="hg19"),
      ChromInfo::new(name="chr19", length=59128983, genome="hg19"),
      ChromInfo::new(name="chr20", length=63025520, genome="hg19"),
      ChromInfo::new(name="chr21", length=48129895, genome="hg19"),
      ChromInfo::new(name="chr22", length=51304566, genome="hg19"),
      ChromInfo::new(name="chrX", length=155270560, genome="hg19"),
      ChromInfo::new(name="chrY", length=59373566, genome="hg19"),
      ChromInfo::new(name="chrM", length=16569, is_circular=true, genome="hg19"),
    ],
    centromeres: [],
  }
}

///|
/// Create a pre-built mouse genome (mm10).
pub fn genome_mm10() -> GenomeInfoDb {
  GenomeInfoDb::{
    genome_build: "mm10",
    organism: "Mus musculus",
    chrom_infos: [
      ChromInfo::new(name="chr1", length=195471971, genome="mm10"),
      ChromInfo::new(name="chr2", length=182113224, genome="mm10"),
      ChromInfo::new(name="chr3", length=160039680, genome="mm10"),
      ChromInfo::new(name="chr4", length=156508116, genome="mm10"),
      ChromInfo::new(name="chr5", length=151834684, genome="mm10"),
      ChromInfo::new(name="chr6", length=149736546, genome="mm10"),
      ChromInfo::new(name="chr7", length=145441459, genome="mm10"),
      ChromInfo::new(name="chr8", length=129401213, genome="mm10"),
      ChromInfo::new(name="chr9", length=124595110, genome="mm10"),
      ChromInfo::new(name="chr10", length=130694993, genome="mm10"),
      ChromInfo::new(name="chr11", length=122082543, genome="mm10"),
      ChromInfo::new(name="chr12", length=120129022, genome="mm10"),
      ChromInfo::new(name="chr13", length=120421639, genome="mm10"),
      ChromInfo::new(name="chr14", length=124902244, genome="mm10"),
      ChromInfo::new(name="chr15", length=104043685, genome="mm10"),
      ChromInfo::new(name="chr16", length=98207768, genome="mm10"),
      ChromInfo::new(name="chr17", length=94987271, genome="mm10"),
      ChromInfo::new(name="chr18", length=90702639, genome="mm10"),
      ChromInfo::new(name="chr19", length=61431566, genome="mm10"),
      ChromInfo::new(name="chrX", length=171031299, genome="mm10"),
      ChromInfo::new(name="chrY", length=91744698, genome="mm10"),
      ChromInfo::new(name="chrM", length=16299, is_circular=true, genome="mm10"),
    ],
    centromeres: [],
  }
}

///|
/// Compare two genome builds.
pub fn compare_genome_builds(build1 : String, build2 : String) -> Bool {
  build1 == build2
}

///|
/// Get supported genome builds.
pub fn supported_genome_builds() -> Array[String] {
  ["hg38", "hg19", "mm10", "mm9", "rn6", "danRer11", "dm6", "ce11"]
}

///|
/// Get organism from genome build.
pub fn organism_from_build(build : String) -> String {
  if build == "hg38" || build == "hg19" || build == "hg18" {
    "Homo sapiens"
  } else if build == "mm10" || build == "mm9" {
    "Mus musculus"
  } else if build == "rn6" || build == "rn5" {
    "Rattus norvegicus"
  } else if build == "danRer11" || build == "danRer10" {
    "Danio rerio"
  } else if build == "dm6" || build == "dm3" {
    "Drosophila melanogaster"
  } else if build == "ce11" || build == "ce10" {
    "Caenorhabditis elegans"
  } else {
    "Unknown"
  }
}

///|
/// Get arm (p or q) for a position on a chromosome.
pub fn GenomeInfoDb::get_arm(
  self : GenomeInfoDb,
  chrom : String,
  position : Int,
) -> String {
  match self.centromere(chrom) {
    Some(c) =>
      if position < c.start {
        "p"
      } else if position > c.end {
        "q"
      } else {
        "centromere"
      }
    None => "unknown"
  }
}

///|
/// Filter chromosomes by pattern.
pub fn GenomeInfoDb::filter_chromosomes(
  self : GenomeInfoDb,
  pattern : String,
) -> Array[String] {
  self.chrom_infos
  .filter(fn(ci) { ci.name.contains(pattern) })
  .map(fn(ci) { ci.name })
}

///|
/// Get genome summary as string.
pub fn GenomeInfoDb::summary(self : GenomeInfoDb) -> String {
  let n = self.chrom_infos.length()
  let total = self.genome_length()
  let build = self.genome_build
  let org = self.organism
  build +
  " (" +
  org +
  "): " +
  n.to_string() +
  " chromosomes, " +
  total.to_string() +
  " bp total"
}
