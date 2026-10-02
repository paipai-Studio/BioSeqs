// ============================================================================
// SnpEff-lite: Variant annotation and effect prediction
// ============================================================================
// A lite reimplementation of SnpEff for variant effect prediction.
// Annotates VCF variants with Sequence Ontology terms, impact levels,
// HGVS notation, and standardized ANN fields.
//
// Pipeline:
//   1. Parse VCF + genome annotation (GFF/GTF)
//   2. For each variant, locate overlapping genes/transcripts
//   3. Predict molecular effect (missense, nonsense, frameshift, etc.)
//   4. Generate HGVS c./p. notation
//   5. Output annotated VCF with ANN info field
//
// Reuses: GenomeAnnotation, VcfRecord/VcfFile, translate_codon

// ============================================================================
// 1. Sequence Ontology terms + impact levels
// ============================================================================

///|
/// Sequence Ontology (SO) variant effect terms.
pub(all) enum SnpEffSO {
  // HIGH impact
  StopGained              // stop_gained
  StopLost                // stop_lost
  StartLost               // start_lost
  FrameshiftVariant       // frameshift_variant
  SpliceDonorVariant      // splice_donor_variant
  SpliceAcceptorVariant   // splice_acceptor_variant
  MissenseVariant         // missense_variant (MODERATE)
  // MODERATE impact
  SpliceRegionVariant     // splice_region_variant
  // LOW impact
  SynonymousVariant       // synonymous_variant
  // MODIFIER impact
  UpstreamVariant         // upstream_gene_variant
  DownstreamVariant       // downstream_gene_variant
  FivePrimeUTRVariant     // 5_prime_UTR_variant
  ThreePrimeUTRVariant    // 3_prime_UTR_variant
  IntronVariant           // intron_variant
  IntergenicVariant       // intergenic_variant
  NonCodingTranscriptVariant // non_coding_transcript_variant
}

///|
/// Impact level.
pub(all) enum SnpEffImpact {
  High
  Moderate
  Low
  Modifier
}

///|
/// Convert SO term to string.
pub fn snpeff_so_string(so : SnpEffSO) -> String {
  match so {
    StopGained => "stop_gained"
    StopLost => "stop_lost"
    StartLost => "start_lost"
    FrameshiftVariant => "frameshift_variant"
    SpliceDonorVariant => "splice_donor_variant"
    SpliceAcceptorVariant => "splice_acceptor_variant"
    MissenseVariant => "missense_variant"
    SpliceRegionVariant => "splice_region_variant"
    SynonymousVariant => "synonymous_variant"
    UpstreamVariant => "upstream_gene_variant"
    DownstreamVariant => "downstream_gene_variant"
    FivePrimeUTRVariant => "5_prime_UTR_variant"
    ThreePrimeUTRVariant => "3_prime_UTR_variant"
    IntronVariant => "intron_variant"
    IntergenicVariant => "intergenic_variant"
    NonCodingTranscriptVariant => "non_coding_transcript_variant"
  }
}

///|
/// Convert impact level to string.
pub fn snpeff_impact_string(impact : SnpEffImpact) -> String {
  match impact {
    High => "HIGH"
    Moderate => "MODERATE"
    Low => "LOW"
    Modifier => "MODIFIER"
  }
}

///|
/// Map SO term to impact level.
pub fn snpeff_so_impact(so : SnpEffSO) -> SnpEffImpact {
  match so {
    StopGained => High
    StopLost => High
    StartLost => High
    FrameshiftVariant => High
    SpliceDonorVariant => High
    SpliceAcceptorVariant => High
    MissenseVariant => Moderate
    SpliceRegionVariant => Moderate
    SynonymousVariant => Low
    UpstreamVariant => Modifier
    DownstreamVariant => Modifier
    FivePrimeUTRVariant => Modifier
    ThreePrimeUTRVariant => Modifier
    IntronVariant => Modifier
    IntergenicVariant => Modifier
    NonCodingTranscriptVariant => Modifier
  }
}

// ============================================================================
// 2. Effect data structures
// ============================================================================

///|
/// A single variant effect prediction.
pub(all) struct SnpEffEffect {
  pub allele : String
  pub so_term : SnpEffSO
  pub impact : SnpEffImpact
  pub gene_name : String
  pub gene_id : String
  pub feature_type : String     // "transcript" or "gene" or "intergenic"
  pub feature_id : String       // transcript_id or gene_id
  pub transcript_biotype : String
  pub rank : Int                // exon rank (0 if not applicable)
  pub hgvs_c : String           // HGVS c. notation
  pub hgvs_p : String           // HGVS p. notation
  pub cdna_pos : Int            // cDNA position
  pub cdna_len : Int            // cDNA length
  pub cds_pos : Int             // CDS position
  pub cds_len : Int             // CDS length
  pub aa_pos : Int              // amino acid position
  pub aa_len : Int              // protein length
  pub distance : Int            // distance to feature (for upstream/downstream)
  pub original_codon : String
  pub mutated_codon : String
  pub original_aa : String
  pub mutated_aa : String
}

///|
/// Annotation result for a single variant.
pub struct SnpEffVariantResult {
  pub chrom : String
  pub pos : Int
  pub ref_allele : String
  pub alt_allele : String
  pub effects : Array[SnpEffEffect]
}

///|
/// Full annotation result.
pub struct SnpEffResult {
  pub seqid : String
  pub num_variants : Int
  pub num_annotated : Int
  pub variant_results : Array[SnpEffVariantResult]
  pub summary : Map[String, Int]  // SO term -> count
}

// ============================================================================
// 3. Helper functions
// ============================================================================

///|
/// Amino acid 3-letter codes for HGVS p. notation.
fn aa_three_letter(aa : String) -> String {
  match aa {
    "A" => "Ala"
    "R" => "Arg"
    "N" => "Asn"
    "D" => "Asp"
    "C" => "Cys"
    "E" => "Glu"
    "Q" => "Gln"
    "G" => "Gly"
    "H" => "His"
    "I" => "Ile"
    "L" => "Leu"
    "K" => "Lys"
    "M" => "Met"
    "F" => "Phe"
    "P" => "Pro"
    "S" => "Ser"
    "T" => "Thr"
    "W" => "Trp"
    "Y" => "Tyr"
    "V" => "Val"
    "*" => "Ter"
    _ => aa
  }
}

///|
/// Extract a subsequence from a reference string.
fn ref_subseq(ref_seq : String, start : Int, end_ : Int) -> String {
  // 1-based inclusive coordinates
  let buf : StringBuilder = StringBuilder()
  for i = start - 1; i < end_; i = i + 1 {
    if i >= 0 && i < ref_seq.length() {
      buf.write_char(ref_seq[i].unsafe_to_char())
    }
  }
  buf.to_string()
}

///|
/// Check if a position is within a range (1-based inclusive).
fn in_range(pos : Int, start : Int, end_ : Int) -> Bool {
  pos >= start && pos <= end_
}

///|
/// Check if a position is in splice donor region (first 2 nt of intron).
fn is_splice_donor(pos : Int, exon_start : Int, exon_end : Int, strand : String) -> Bool {
  if strand == "+" {
    // Donor = first 2 nt after exon end (intron start)
    pos == exon_end + 1 || pos == exon_end + 2
  } else {
    // On minus strand, donor = last 2 nt before exon start
    pos == exon_start - 1 || pos == exon_start - 2
  }
}

///|
/// Check if a position is in splice acceptor region (last 2 nt of intron).
fn is_splice_acceptor(pos : Int, exon_start : Int, exon_end : Int, strand : String) -> Bool {
  if strand == "+" {
    // Acceptor = last 2 nt before exon start (intron end)
    pos == exon_start - 1 || pos == exon_start - 2
  } else {
    // On minus strand, acceptor = first 2 nt after exon end
    pos == exon_end + 1 || pos == exon_end + 2
  }
}

///|
/// Check if a position is in splice region (±1 to ±3 around splice sites).
fn is_splice_region(pos : Int, exon_start : Int, exon_end : Int) -> Bool {
  let dist_to_start = (pos - exon_start).abs()
  let dist_to_end = (pos - exon_end).abs()
  dist_to_start <= 3 || dist_to_end <= 3
}

// ============================================================================
// 4. Effect prediction
// ============================================================================

///|
/// Predict the effect of a single variant on a transcript.
fn predict_effect_on_transcript(
  chrom : String,
  pos : Int,          // 1-based VCF position
  ref_allele : String,
  alt_allele : String,
  gene : GFGene,
  tx : GFTranscript,
  ref_seq : String,   // reference sequence for this chrom
) -> SnpEffEffect {
  let is_snp = ref_allele.length() == 1 && alt_allele.length() == 1
  let is_frameshift = (ref_allele.length() - alt_allele.length()) % 3 != 0 ||
    (alt_allele.length() - ref_allele.length()) % 3 != 0

  // Default intergenic effect
  let empty_effect : SnpEffEffect = {
    allele: alt_allele,
    so_term: IntergenicVariant,
    impact: Modifier,
    gene_name: gene.gene_name,
    gene_id: gene.gene_id,
    feature_type: "transcript",
    feature_id: tx.transcript_id,
    transcript_biotype: "",
    rank: 0,
    hgvs_c: "",
    hgvs_p: "",
    cdna_pos: 0,
    cdna_len: 0,
    cds_pos: 0,
    cds_len: 0,
    aa_pos: 0,
    aa_len: 0,
    distance: 0,
    original_codon: "",
    mutated_codon: "",
    original_aa: "",
    mutated_aa: "",
  }

  // Check splice sites first (before exon/CDS)
  for exon in tx.exons {
    // Splice donor
    if is_splice_donor(pos, exon.start, exon.end, tx.strand) {
      return { ..empty_effect, so_term: SpliceDonorVariant, impact: High }
    }
    // Splice acceptor
    if is_splice_acceptor(pos, exon.start, exon.end, tx.strand) {
      return { ..empty_effect, so_term: SpliceAcceptorVariant, impact: High }
    }
  }

  // Check splice region (broader, ±3)
  for exon in tx.exons {
    if is_splice_region(pos, exon.start, exon.end) {
      // Only if not already donor/acceptor
      let in_exon = in_range(pos, exon.start, exon.end)
      if not(in_exon) {
        return { ..empty_effect, so_term: SpliceRegionVariant, impact: Moderate }
      }
    }
  }

  // Check CDS (coding sequence)
  for cds in tx.cds {
    if in_range(pos, cds.start, cds.end) {
      // Compute CDS position relative to transcript
      let cds_pos = pos - cds.start + 1
      let codon_offset = (cds_pos - 1) % 3
      let codon_start_cds = cds_pos - codon_offset
      let codon_start_genome = cds.start + codon_start_cds - 1

      // Extract original codon from reference
      let original_codon = ref_subseq(
        ref_seq,
        codon_start_genome,
        codon_start_genome + 2,
      )

      // Build mutated codon
      let mutated_codon = if is_snp && original_codon.length() == 3 {
        let mut mc = ""
        for k = 0; k < 3; k = k + 1 {
          if k == codon_offset {
            mc += alt_allele
          } else {
            mc += original_codon[k].unsafe_to_char().to_string()
          }
        }
        mc
      } else {
        original_codon
      }

      let original_aa = translate_codon(original_codon)
      let mutated_aa = translate_codon(mutated_codon)
      let aa_pos = (codon_start_cds + 2) / 3

      // Determine effect
      let (so_term, impact) = if is_frameshift {
        (FrameshiftVariant, High)
      } else if is_snp && original_codon.length() == 3 {
        if mutated_aa == "*" && original_aa != "*" {
          (StopGained, High)
        } else if original_aa == "*" && mutated_aa != "*" {
          (StopLost, High)
        } else if original_aa == "M" && aa_pos == 1 {
          (StartLost, High)
        } else if original_aa == mutated_aa {
          (SynonymousVariant, Low)
        } else {
          (MissenseVariant, Moderate)
        }
      } else {
        (MissenseVariant, Moderate)
      }

      // HGVS c. notation
      let hgvs_c = if is_snp && original_codon.length() == 3 {
        "c." + cds_pos.to_string() + ref_allele + ">" + alt_allele
      } else if ref_allele.length() > alt_allele.length() {
        // Deletion
        let del_len = ref_allele.length() - alt_allele.length()
        "c." + cds_pos.to_string() + "_" + (cds_pos + del_len - 1).to_string() + "del"
      } else if alt_allele.length() > ref_allele.length() {
        // Insertion
        "c." + (cds_pos - 1).to_string() + "_" + cds_pos.to_string() + "ins" + alt_allele[1:].to_string()
      } else {
        "c." + cds_pos.to_string() + "delins" + alt_allele
      }

      // HGVS p. notation
      let hgvs_p = if is_snp && original_codon.length() == 3 && original_aa != "?" && mutated_aa != "?" {
        "p." + aa_three_letter(original_aa) + aa_pos.to_string() + aa_three_letter(mutated_aa)
      } else if is_frameshift {
        "p." + aa_three_letter(original_aa) + aa_pos.to_string() + "fs"
      } else {
        ""
      }

      // Find exon rank
      let mut exon_rank = 0
      for ei = 0; ei < tx.exons.length(); ei = ei + 1 {
        if in_range(pos, tx.exons[ei].start, tx.exons[ei].end) {
          exon_rank = tx.exons[ei].exon_number
          break
        }
      }

      return {
        allele: alt_allele,
        so_term,
        impact,
        gene_name: gene.gene_name,
        gene_id: gene.gene_id,
        feature_type: "transcript",
        feature_id: tx.transcript_id,
        transcript_biotype: gene.gene_biotype,
        rank: exon_rank,
        hgvs_c,
        hgvs_p,
        cdna_pos: cds_pos,
        cdna_len: tx.end - tx.start + 1,
        cds_pos,
        cds_len: 0, // computed from total CDS length
        aa_pos,
        aa_len: 0,
        distance: 0,
        original_codon,
        mutated_codon,
        original_aa,
        mutated_aa,
      }
    }
  }

  // Check 5' UTR
  for utr in tx.utr5 {
    if in_range(pos, utr.start, utr.end) {
      return { ..empty_effect, so_term: FivePrimeUTRVariant, impact: Modifier }
    }
  }

  // Check 3' UTR
  for utr in tx.utr3 {
    if in_range(pos, utr.start, utr.end) {
      return { ..empty_effect, so_term: ThreePrimeUTRVariant, impact: Modifier }
    }
  }

  // Check introns (in exon range but not in exon itself)
  let in_transcript = in_range(pos, tx.start, tx.end)
  let mut in_exon = false
  for exon in tx.exons {
    if in_range(pos, exon.start, exon.end) {
      in_exon = true
    }
  }
  if in_transcript && not(in_exon) && tx.exons.length() > 1 {
    return { ..empty_effect, so_term: IntronVariant, impact: Modifier }
  }

  // Non-coding transcript
  if in_transcript && tx.cds.length() == 0 {
    return { ..empty_effect, so_term: NonCodingTranscriptVariant, impact: Modifier }
  }

  empty_effect
}

///|
/// Predict effect for a variant with no overlapping gene (intergenic).
fn predict_intergenic(
  chrom : String,
  pos : Int,
  ref_allele : String,
  alt_allele : String,
  annotation : GenomeAnnotation,
) -> Array[SnpEffEffect] {
  let effects : Array[SnpEffEffect] = Array(capacity=0)
  // Find nearest genes upstream and downstream
  let mut nearest_up : GFGene? = None
  let mut nearest_up_dist = 999999
  let mut nearest_down : GFGene? = None
  let mut nearest_down_dist = 999999
  for gene in annotation.genes {
    if gene.chromosome != chrom {
      continue
    }
    if gene.end < pos {
      // Downstream (variant is 3' of gene)
      let dist = pos - gene.end
      if dist < nearest_down_dist {
        nearest_down_dist = dist
        nearest_down = Some(gene)
      }
    } else if gene.start > pos {
      // Upstream (variant is 5' of gene)
      let dist = gene.start - pos
      if dist < nearest_up_dist {
        nearest_up_dist = dist
        nearest_up = Some(gene)
      }
    }
  }
  match nearest_up {
    Some(g) =>
      effects.push({
        allele: alt_allele,
        so_term: UpstreamVariant,
        impact: Modifier,
        gene_name: g.gene_name,
        gene_id: g.gene_id,
        feature_type: "gene",
        feature_id: g.gene_id,
        transcript_biotype: "",
        rank: 0,
        hgvs_c: "",
        hgvs_p: "",
        cdna_pos: 0,
        cdna_len: 0,
        cds_pos: 0,
        cds_len: 0,
        aa_pos: 0,
        aa_len: 0,
        distance: nearest_up_dist,
        original_codon: "",
        mutated_codon: "",
        original_aa: "",
        mutated_aa: "",
      })
    None => ()
  }
  match nearest_down {
    Some(g) =>
      effects.push({
        allele: alt_allele,
        so_term: DownstreamVariant,
        impact: Modifier,
        gene_name: g.gene_name,
        gene_id: g.gene_id,
        feature_type: "gene",
        feature_id: g.gene_id,
        transcript_biotype: "",
        rank: 0,
        hgvs_c: "",
        hgvs_p: "",
        cdna_pos: 0,
        cdna_len: 0,
        cds_pos: 0,
        cds_len: 0,
        aa_pos: 0,
        aa_len: 0,
        distance: nearest_down_dist,
        original_codon: "",
        mutated_codon: "",
        original_aa: "",
        mutated_aa: "",
      })
    None => ()
  }
  if effects.length() == 0 {
    effects.push({
      allele: alt_allele,
      so_term: IntergenicVariant,
      impact: Modifier,
      gene_name: "",
      gene_id: "",
      feature_type: "intergenic",
      feature_id: "",
      transcript_biotype: "",
      rank: 0,
      hgvs_c: "",
      hgvs_p: "",
      cdna_pos: 0,
      cdna_len: 0,
      cds_pos: 0,
      cds_len: 0,
      aa_pos: 0,
      aa_len: 0,
      distance: 0,
      original_codon: "",
      mutated_codon: "",
      original_aa: "",
      mutated_aa: "",
    })
  }
  effects
}

// ============================================================================
// 5. Main annotation pipeline
// ============================================================================

///|
/// Annotate a VCF file with variant effects.
pub fn snpeff_annotate(
  vcf : VcfFile,
  annotation : GenomeAnnotation,
  reference : Map[String, String],
) -> SnpEffResult {
  let variant_results : Array[SnpEffVariantResult] = Array(capacity=0)
  let summary : Map[String, Int] = Map([], capacity=0)
  let mut num_annotated = 0
  for rec in vcf.records {
    let ref_seq = match reference.get(rec.chrom) {
      Some(s) => s
      None => ""
    }
    // Process each alt allele
    for alt in rec.alt {
      let effects = snpeff_predict_variant(
        rec.chrom,
        rec.pos,
        rec.ref_seq,
        alt,
        annotation,
        ref_seq,
      )
      if effects.length() > 0 {
        num_annotated += 1
      }
      for eff in effects {
        let so_str = snpeff_so_string(eff.so_term)
        match summary.get(so_str) {
          Some(c) => summary[so_str] = c + 1
          None => summary[so_str] = 1
        }
      }
      variant_results.push({
        chrom: rec.chrom,
        pos: rec.pos,
        ref_allele: rec.ref_seq,
        alt_allele: alt,
        effects,
      })
    }
  }
  {
    seqid: annotation.genome,
    num_variants: variant_results.length(),
    num_annotated,
    variant_results,
    summary,
  }
}

///|
/// Predict effects for a single variant.
pub fn snpeff_predict_variant(
  chrom : String,
  pos : Int,
  ref_allele : String,
  alt_allele : String,
  annotation : GenomeAnnotation,
  ref_seq : String,
) -> Array[SnpEffEffect] {
  let effects : Array[SnpEffEffect] = Array(capacity=0)
  let mut found_gene = false
  // Find overlapping genes
  for gene in annotation.genes {
    if gene.chromosome != chrom {
      continue
    }
    if in_range(pos, gene.start, gene.end) {
      found_gene = true
      // Check each transcript
      for tx in gene.transcripts {
        if in_range(pos, tx.start, tx.end) {
          let eff = predict_effect_on_transcript(
            chrom,
            pos,
            ref_allele,
            alt_allele,
            gene,
            tx,
            ref_seq,
          )
          effects.push(eff)
        }
      }
      // If in gene but not in any transcript, still annotate
      if effects.length() == 0 {
        for tx in gene.transcripts {
          let eff = predict_effect_on_transcript(
            chrom,
            pos,
            ref_allele,
            alt_allele,
            gene,
            tx,
            ref_seq,
          )
          effects.push(eff)
        }
      }
    }
  }
  // If no overlapping gene, find intergenic/upstream/downstream
  if not(found_gene) || effects.length() == 0 {
    let intergenic = predict_intergenic(
      chrom,
      pos,
      ref_allele,
      alt_allele,
      annotation,
    )
    for eff in intergenic {
      effects.push(eff)
    }
  }
  effects
}

// ============================================================================
// 6. ANN field output
// ============================================================================

///|
/// Generate ANN field string for a single effect.
pub fn snpeff_to_ann(effect : SnpEffEffect) -> String {
  // ANN format:
  // Allele|Annotation|Annotation_Impact|Gene_Name|Gene_ID|
  // Feature_Type|Feature_ID|Transcript_BioType|Rank|HGVS.c|HGVS.p|
  // cDNA.pos|cDNA.length|CDS.pos|CDS.length|AA.pos|AA.length|Distance
  let parts : Array[String] = Array(capacity=17)
  parts.push(effect.allele)
  parts.push(snpeff_so_string(effect.so_term))
  parts.push(snpeff_impact_string(effect.impact))
  parts.push(effect.gene_name)
  parts.push(effect.gene_id)
  parts.push(effect.feature_type)
  parts.push(effect.feature_id)
  parts.push(effect.transcript_biotype)
  parts.push(if effect.rank > 0 { effect.rank.to_string() } else { "" })
  parts.push(effect.hgvs_c)
  parts.push(effect.hgvs_p)
  parts.push(if effect.cdna_pos > 0 { effect.cdna_pos.to_string() + "/" + effect.cdna_len.to_string() } else { "" })
  parts.push(if effect.cds_pos > 0 { effect.cds_pos.to_string() + "/" + effect.cds_len.to_string() } else { "" })
  parts.push(if effect.aa_pos > 0 { effect.aa_pos.to_string() + "/" + effect.aa_len.to_string() } else { "" })
  parts.push(if effect.distance > 0 { effect.distance.to_string() } else { "" })
  parts.push("")  // ERRORS / WARNINGS / INFO field
  // Join with |
  let buf : StringBuilder = StringBuilder()
  for i = 0; i < parts.length(); i = i + 1 {
    if i > 0 {
      buf.write_char('|')
    }
    buf.write_string(parts[i])
  }
  buf.to_string()
}

///|
/// Generate ANN field for all effects of a variant.
pub fn snpeff_to_ann_field(effects : Array[SnpEffEffect]) -> String {
  if effects.length() == 0 {
    return ""
  }
  let buf : StringBuilder = StringBuilder()
  for i = 0; i < effects.length(); i = i + 1 {
    if i > 0 {
      buf.write_char(',')
    }
    buf.write_string(snpeff_to_ann(effects[i]))
  }
  buf.to_string()
}

// ============================================================================
// 7. VCF output with ANN field
// ============================================================================

///|
/// Annotate a VCF file and return annotated VCF string.
pub fn snpeff_annotate_vcf(
  vcf : VcfFile,
  annotation : GenomeAnnotation,
  reference : Map[String, String],
) -> String {
  let result = snpeff_annotate(vcf, annotation, reference)
  let buf : StringBuilder = StringBuilder()
  // Write meta headers
  buf.write_string("##fileformat=VCFv4.2\n")
  buf.write_string("##INFO=<ID=ANN,Number=1,Type=String,Description=\"Functional annotations: 'Allele|Annotation|Annotation_Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|Transcript_BioType|Rank|HGVS.c|HGVS.p|cDNA.pos|cDNA.length|CDS.pos|CDS.length|AA.pos|AA.length|Distance'\">\n")
  // Copy existing meta headers
  for header in vcf.meta_headers {
    if header.key != "INFO" || not(header.value.contains("ANN")) {
      buf.write_string("##")
      buf.write_string(header.key)
      buf.write_string("=")
      buf.write_string(header.value)
      buf.write_string("\n")
    }
  }
  // Write column header
  buf.write_string("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
  match vcf.fixed_header {
    Some(fh) =>
      for i = 8; i < fh.columns.length(); i = i + 1 {
        buf.write_string("\t")
        buf.write_string(fh.columns[i])
      }
    None => ()
  }
  buf.write_string("\n")
  // Write records with ANN
  let mut ri = 0
  for rec in vcf.records {
    buf.write_string(rec.chrom)
    buf.write_string("\t")
    buf.write_string(rec.pos.to_string())
    buf.write_string("\t")
    // ID
    if rec.id.length() > 0 {
      for i = 0; i < rec.id.length(); i = i + 1 {
        if i > 0 {
          buf.write_char(';')
        }
        buf.write_string(rec.id[i])
      }
    } else {
      buf.write_char('.')
    }
    buf.write_string("\t")
    buf.write_string(rec.ref_seq)
    buf.write_string("\t")
    // ALT
    for i = 0; i < rec.alt.length(); i = i + 1 {
      if i > 0 {
        buf.write_char(',')
      }
      buf.write_string(rec.alt[i])
    }
    buf.write_string("\t")
    // QUAL
    match rec.qual {
      Some(q) => buf.write_string(q.to_string())
      None => buf.write_char('.')
    }
    buf.write_string("\t")
    // FILTER
    if rec.filter.length() > 0 {
      for i = 0; i < rec.filter.length(); i = i + 1 {
        if i > 0 {
          buf.write_char(';')
        }
        buf.write_string(rec.filter[i])
      }
    } else {
      buf.write_char('.')
    }
    buf.write_string("\t")
    // INFO with ANN
    let mut has_info = false
    for (k, v) in rec.info {
      if k != "ANN" {
        if has_info {
          buf.write_char(';')
        }
        buf.write_string(k)
        buf.write_string("=")
        buf.write_string(v)
        has_info = true
      }
    }
    // Add ANN field
    if ri < result.variant_results.length() {
      let vr = result.variant_results[ri]
      if vr.effects.length() > 0 {
        if has_info {
          buf.write_char(';')
        }
        buf.write_string("ANN=")
        buf.write_string(snpeff_to_ann_field(vr.effects))
        has_info = true
      }
    }
    if not(has_info) {
      buf.write_char('.')
    }
    // FORMAT + samples
    match rec.format {
      Some(fmt) => {
        buf.write_string("\t")
        for i = 0; i < fmt.length(); i = i + 1 {
          if i > 0 {
            buf.write_char(':')
          }
          buf.write_string(fmt[i])
        }
        for (sample_name, vals) in rec.samples {
          buf.write_string("\t")
          for i = 0; i < vals.length(); i = i + 1 {
            if i > 0 {
              buf.write_char(':')
            }
            buf.write_string(vals[i])
          }
          ignore(sample_name)
        }
      }
      None => ()
    }
    buf.write_string("\n")
    ri += 1
  }
  buf.to_string()
}

// ============================================================================
// 8. Summary and utilities
// ============================================================================

///|
/// Generate a summary of annotation results.
pub fn snpeff_summary_text(result : SnpEffResult) -> String {
  let buf : StringBuilder = StringBuilder()
  buf.write_string("SnpEff annotation summary\n")
  buf.write_string("=========================\n")
  buf.write_string("Genome: ")
  buf.write_string(result.seqid)
  buf.write_string("\n")
  buf.write_string("Total variants: ")
  buf.write_string(result.num_variants.to_string())
  buf.write_string("\n")
  buf.write_string("Annotated variants: ")
  buf.write_string(result.num_annotated.to_string())
  buf.write_string("\n\n")
  buf.write_string("Effect counts:\n")
  // Sort by count descending
  let entries : Array[(String, Int)] = Array(capacity=0)
  for (k, v) in result.summary {
    entries.push((k, v))
  }
  for i = 0; i < entries.length(); i = i + 1 {
    for j = i + 1; j < entries.length(); j = j + 1 {
      if entries[j].1 > entries[i].1 {
        let tmp = entries[i]
        entries[i] = entries[j]
        entries[j] = tmp
      }
    }
  }
  for (term, count) in entries {
    buf.write_string("  ")
    buf.write_string(term)
    buf.write_string(": ")
    buf.write_string(count.to_string())
    buf.write_string("\n")
  }
  buf.to_string()
}

///|
/// Split multi-allelic VCF records into single-allelic.
pub fn snpeff_split_multiallelic(vcf : VcfFile) -> VcfFile {
  let result : VcfFile = {
    meta_headers: vcf.meta_headers,
    fixed_header: vcf.fixed_header,
    records: Array(capacity=0),
  }
  for rec in vcf.records {
    if rec.alt.length() <= 1 {
      result.records.push(rec)
    } else {
      for alt in rec.alt {
        let new_rec : VcfRecord = {
          chrom: rec.chrom,
          pos: rec.pos,
          id: rec.id,
          ref_seq: rec.ref_seq,
          alt: [alt],
          qual: rec.qual,
          filter: rec.filter,
          info: rec.info,
          format: rec.format,
          samples: rec.samples,
        }
        result.records.push(new_rec)
      }
    }
  }
  result
}

///|
/// Count effects by impact level.
pub fn snpeff_count_by_impact(
  result : SnpEffResult,
) -> Map[String, Int] {
  let counts : Map[String, Int] = Map([], capacity=0)
  for vr in result.variant_results {
    for eff in vr.effects {
      let impact_str = snpeff_impact_string(eff.impact)
      match counts.get(impact_str) {
        Some(c) => counts[impact_str] = c + 1
        None => counts[impact_str] = 1
      }
    }
  }
  counts
}
