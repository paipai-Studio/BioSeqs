///|
/// VariantAnnotation module, inspired by Bioconductor VariantAnnotation.

///|
pub enum VariantType {
  SNP
  Insertion
  Deletion
  Indel
}

///|
pub fn VariantType::show(self : VariantType) -> String {
  match self {
    SNP => "SNP"
    Insertion => "Insertion"
    Deletion => "Deletion"
    Indel => "Indel"
  }
}

///|
pub struct VariantAnnotation {
  chrom : String
  pos : Int
  ref_allele : String
  alt_allele : String
  variant_type : VariantType
  gene_id : String
  gene_name : String
  transcript_id : String
  exon_id : String
  cds_position : Int
  aa_change : String
  codon_change : String
  consequence : String
}

///|
pub fn VariantAnnotation::new(
  chrom : String,
  pos : Int,
  ref_allele : String,
  alt_allele : String,
  variant_type : VariantType,
  gene_id : String,
  gene_name : String,
  transcript_id : String,
  exon_id : String,
  cds_position : Int,
  aa_change : String,
  codon_change : String,
  consequence : String,
) -> VariantAnnotation {
  VariantAnnotation::{
    chrom,
    pos,
    ref_allele,
    alt_allele,
    variant_type,
    gene_id,
    gene_name,
    transcript_id,
    exon_id,
    cds_position,
    aa_change,
    codon_change,
    consequence,
  }
}

///|
pub fn detect_variant_type(
  ref_allele : String,
  alt_allele : String,
) -> VariantType {
  let ref_len = ref_allele.length()
  let alt_len = alt_allele.length()

  if ref_len == 1 && alt_len == 1 {
    VariantType::SNP
  } else if ref_len == 1 && alt_len > 1 {
    VariantType::Insertion
  } else if ref_len > 1 && alt_len == 1 {
    VariantType::Deletion
  } else {
    VariantType::Indel
  }
}

///|
pub fn locate_variants(
  vcf_records : Array[VcfRecord],
  txdb : TxDb,
) -> Array[VariantAnnotation] {
  let result : Array[VariantAnnotation] = Array::new()

  let mut i = 0
  while i < vcf_records.length() {
    let rec = vcf_records[i]

    let mut j = 0
    while j < txdb.genes.length() {
      let gene = txdb.genes[j]
      if rec.chrom == gene.seqname &&
        rec.pos >= gene.start &&
        rec.pos <= gene.end {
        let mut k = 0
        while k < txdb.transcripts.length() {
          let tx = txdb.transcripts[k]
          if tx.gene_id == gene.gene_id {
            let mut l = 0
            while l < txdb.exons.length() {
              let exon = txdb.exons[l]
              if exon.transcript_id == tx.transcript_id &&
                rec.pos >= exon.start &&
                rec.pos <= exon.end {
                let alt = if rec.alt.length() > 0 { rec.alt[0] } else { "" }
                let va = VariantAnnotation::new(
                  rec.chrom,
                  rec.pos,
                  rec.ref,
                  alt,
                  detect_variant_type(rec.ref, alt),
                  gene.gene_id,
                  gene.gene_id,
                  tx.transcript_id,
                  exon.exon_id,
                  0,
                  "",
                  "",
                  "exon",
                )
                result.push(va)
              }
              l = l + 1
            }
          }
          k = k + 1
        }
      }
      j = j + 1
    }

    let mut m = 0
    while m < txdb.cds.length() {
      let cds = txdb.cds[m]
      if rec.chrom == cds.seqname && rec.pos >= cds.start && rec.pos <= cds.end {
        let alt = if rec.alt.length() > 0 { rec.alt[0] } else { "" }
        let va = VariantAnnotation::new(
          rec.chrom,
          rec.pos,
          rec.ref,
          alt,
          detect_variant_type(rec.ref, alt),
          cds.gene_id,
          cds.gene_id,
          cds.transcript_id,
          "",
          rec.pos - cds.start + 1,
          "",
          "",
          "coding",
        )
        result.push(va)
      }
      m = m + 1
    }

    i = i + 1
  }

  result
}

///|
pub struct CodingPrediction {
  chrom : String
  pos : Int
  ref_allele : String
  alt_allele : String
  gene_id : String
  transcript_id : String
  cds_pos : Int
  codon_pos : Int
  original_codon : String
  mutated_codon : String
  original_aa : String
  mutated_aa : String
  aa_change : String
  consequence : String
}

///|
pub fn CodingPrediction::new(
  chrom : String,
  pos : Int,
  ref_allele : String,
  alt_allele : String,
  gene_id : String,
  transcript_id : String,
  cds_pos : Int,
  codon_pos : Int,
  original_codon : String,
  mutated_codon : String,
  original_aa : String,
  mutated_aa : String,
  aa_change : String,
  consequence : String,
) -> CodingPrediction {
  CodingPrediction::{
    chrom,
    pos,
    ref_allele,
    alt_allele,
    gene_id,
    transcript_id,
    cds_pos,
    codon_pos,
    original_codon,
    mutated_codon,
    original_aa,
    mutated_aa,
    aa_change,
    consequence,
  }
}

///|
pub fn translate_codon(codon : String) -> String {
  let codon_upper = codon.to_upper()
  if codon_upper == "TTT" || codon_upper == "TTC" {
    return "F"
  }
  if codon_upper == "TTA" ||
    codon_upper == "TTG" ||
    codon_upper == "CTT" ||
    codon_upper == "CTC" ||
    codon_upper == "CTA" ||
    codon_upper == "CTG" {
    return "L"
  }
  if codon_upper == "TCT" ||
    codon_upper == "TCC" ||
    codon_upper == "TCA" ||
    codon_upper == "TCG" ||
    codon_upper == "AGT" ||
    codon_upper == "AGC" {
    return "S"
  }
  if codon_upper == "TAT" || codon_upper == "TAC" {
    return "Y"
  }
  if codon_upper == "TAA" || codon_upper == "TAG" || codon_upper == "TGA" {
    return "*"
  }
  if codon_upper == "TGT" || codon_upper == "TGC" {
    return "C"
  }
  if codon_upper == "TGG" {
    return "W"
  }
  if codon_upper == "CCT" ||
    codon_upper == "CCC" ||
    codon_upper == "CCA" ||
    codon_upper == "CCG" {
    return "P"
  }
  if codon_upper == "CAT" || codon_upper == "CAC" {
    return "H"
  }
  if codon_upper == "CAA" || codon_upper == "CAG" {
    return "Q"
  }
  if codon_upper == "CGT" ||
    codon_upper == "CGC" ||
    codon_upper == "CGA" ||
    codon_upper == "CGG" ||
    codon_upper == "AGA" ||
    codon_upper == "AGG" {
    return "R"
  }
  if codon_upper == "ATT" || codon_upper == "ATC" || codon_upper == "ATA" {
    return "I"
  }
  if codon_upper == "ATG" {
    return "M"
  }
  if codon_upper == "ACT" ||
    codon_upper == "ACC" ||
    codon_upper == "ACA" ||
    codon_upper == "ACG" {
    return "T"
  }
  if codon_upper == "AAT" || codon_upper == "AAC" {
    return "N"
  }
  if codon_upper == "AAA" || codon_upper == "AAG" {
    return "K"
  }
  if codon_upper == "GTT" ||
    codon_upper == "GTC" ||
    codon_upper == "GTA" ||
    codon_upper == "GTG" {
    return "V"
  }
  if codon_upper == "GCT" ||
    codon_upper == "GCC" ||
    codon_upper == "GCA" ||
    codon_upper == "GCG" {
    return "A"
  }
  if codon_upper == "GAT" || codon_upper == "GAC" {
    return "D"
  }
  if codon_upper == "GAA" || codon_upper == "GAG" {
    return "E"
  }
  if codon_upper == "GGT" ||
    codon_upper == "GGC" ||
    codon_upper == "GGA" ||
    codon_upper == "GGG" {
    return "G"
  }
  if codon.length() != 3 {
    return "?"
  }
  "?"
}

///|
pub fn predict_coding(
  vcf_records : Array[VcfRecord],
  txdb : TxDb,
  reference : @hashmap.HashMap[String, String],
) -> Array[CodingPrediction] {
  let result : Array[CodingPrediction] = Array::new()

  let mut i = 0
  while i < vcf_records.length() {
    let rec = vcf_records[i]

    let ref_seq = match reference.get(rec.chrom) {
      Some(s) => s
      None => ""
    }

    if ref_seq == "" {
      i = i + 1
      continue
    }

    let mut k = 0
    while k < txdb.cds.length() {
      let cds = txdb.cds[k]
      if rec.chrom == cds.seqname && rec.pos >= cds.start && rec.pos <= cds.end {
        let cds_pos = rec.pos - cds.start + 1
        let codon_pos = (cds_pos - 1) % 3
        let codon_start = cds_pos - codon_pos

        let mut original_codon = ""
        if codon_start >= 1 && codon_start + 2 <= cds.end - cds.start + 1 {
          let ref_start = cds.start + codon_start - 1
          if ref_start >= 0 && ref_start + 2 < ref_seq.length() {
            original_codon = ref_seq[ref_start:ref_start + 3].to_string()
          }
        }

        let mut mutated_codon = original_codon
        let alt = if rec.alt.length() > 0 { rec.alt[0] } else { "" }
        if original_codon.length() == 3 &&
          codon_pos >= 0 &&
          codon_pos < 3 &&
          alt.length() == 1 {
          let c0 = original_codon[0:1].to_string()
          let c1 = original_codon[1:2].to_string()
          let c2 = original_codon[2:3].to_string()
          if codon_pos == 0 {
            mutated_codon = alt + c1 + c2
          } else if codon_pos == 1 {
            mutated_codon = c0 + alt + c2
          } else {
            mutated_codon = c0 + c1 + alt
          }
        }

        let original_aa = translate_codon(original_codon)
        let mutated_aa = translate_codon(mutated_codon)

        let aa_change = if original_aa == mutated_aa {
          "synonymous"
        } else if mutated_aa == "*" {
          "nonsense"
        } else {
          "missense"
        }

        let consequence = if mutated_aa == "*" {
          "stop_gained"
        } else if original_aa != mutated_aa {
          "missense_variant"
        } else {
          "synonymous_variant"
        }

        let cp = CodingPrediction::new(
          rec.chrom,
          rec.pos,
          rec.ref,
          alt,
          cds.gene_id,
          cds.transcript_id,
          cds_pos,
          codon_pos,
          original_codon,
          mutated_codon,
          original_aa,
          mutated_aa,
          aa_change,
          consequence,
        )
        result.push(cp)
      }
      k = k + 1
    }

    i = i + 1
  }

  result
}

///|
pub struct VariantSummary {
  mut total_variants : Int
  mut snp_count : Int
  mut insertion_count : Int
  mut deletion_count : Int
  mut indel_count : Int
  mut coding_variants : Int
  mut non_coding_variants : Int
  mut synonymous_count : Int
  mut missense_count : Int
  mut nonsense_count : Int
  gene_ids : Array[String]
  gene_counts : Array[Int]
}

///|
pub fn VariantSummary::new() -> VariantSummary {
  VariantSummary::{
    total_variants: 0,
    snp_count: 0,
    insertion_count: 0,
    deletion_count: 0,
    indel_count: 0,
    coding_variants: 0,
    non_coding_variants: 0,
    synonymous_count: 0,
    missense_count: 0,
    nonsense_count: 0,
    gene_ids: [],
    gene_counts: [],
  }
}

///|
pub fn variant_summary(
  vcf_records : Array[VcfRecord],
  annotations : Array[VariantAnnotation],
) -> VariantSummary {
  let summary = VariantSummary::new()

  summary.total_variants = vcf_records.length()

  let mut i = 0
  while i < vcf_records.length() {
    let rec = vcf_records[i]
    let alt = if rec.alt.length() > 0 { rec.alt[0] } else { "" }
    let vt = detect_variant_type(rec.ref, alt)

    if vt is VariantType::SNP {
      summary.snp_count = summary.snp_count + 1
    }
    if vt is VariantType::Insertion {
      summary.insertion_count = summary.insertion_count + 1
    }
    if vt is VariantType::Deletion {
      summary.deletion_count = summary.deletion_count + 1
    }
    if vt is VariantType::Indel {
      summary.indel_count = summary.indel_count + 1
    }

    i = i + 1
  }

  let mut j = 0
  while j < annotations.length() {
    let va = annotations[j]
    if va.consequence == "coding" {
      summary.coding_variants = summary.coding_variants + 1
    } else {
      summary.non_coding_variants = summary.non_coding_variants + 1
    }

    let mut found = false
    let mut k = 0
    while k < summary.gene_ids.length() {
      if summary.gene_ids[k] == va.gene_id {
        summary.gene_counts[k] = summary.gene_counts[k] + 1
        found = true
        k = summary.gene_ids.length()
      }
      k = k + 1
    }
    if !found {
      summary.gene_ids.push(va.gene_id)
      summary.gene_counts.push(1)
    }

    j = j + 1
  }

  summary
}

///|
pub fn variant_summary_with_coding(
  vcf_records : Array[VcfRecord],
  predictions : Array[CodingPrediction],
) -> VariantSummary {
  let summary = variant_summary(vcf_records, [])

  let mut i = 0
  while i < predictions.length() {
    let cp = predictions[i]

    if cp.consequence == "synonymous_variant" {
      summary.synonymous_count = summary.synonymous_count + 1
    }
    if cp.consequence == "missense_variant" {
      summary.missense_count = summary.missense_count + 1
    }
    if cp.consequence == "stop_gained" {
      summary.nonsense_count = summary.nonsense_count + 1
    }

    let mut found = false
    let mut j = 0
    while j < summary.gene_ids.length() {
      if summary.gene_ids[j] == cp.gene_id {
        summary.gene_counts[j] = summary.gene_counts[j] + 1
        found = true
        j = summary.gene_ids.length()
      }
      j = j + 1
    }
    if !found {
      summary.gene_ids.push(cp.gene_id)
      summary.gene_counts.push(1)
    }

    i = i + 1
  }

  summary
}
