///|
/// pysam_bcf_call — `bcftools call` 的内存版实现 (PL + GT + GQ)
///
/// 在 gatk_lite.mbt 的 pileup site 核心之上，额外提供：
/// - diploid genotype call with GQ / PL (Phred-scaled likelihoods)
/// - -c (consensus) and -m (multiallelic) mode 枚举
/// - convenience GT/PL 字符串生成
///
/// 这是 bcftools call 在 BioSeqs 里的基因型层包装。

///|
/// Genotype call mode, matching bcftools call -c / -m flags.
pub enum BcfCallMode {
  /// `-c` : consensus calling — one REF/ALT pair per site
  Consensus
  /// `-m` : multiallelic calling — up to 3 alleles per site
  Multiallelic
}

///|
/// Configuration for one call run.
pub struct BcfCallConfig {
  min_depth : Int32
  min_af : Double
  min_gq : Double
  mode : BcfCallMode
  sample_name : String
  write_pl : Bool
}

///|
/// Default config — matches `bcftools call -m -v`.
pub fn BcfCallConfig::default() -> BcfCallConfig {
  BcfCallConfig::{
    min_depth: 4,
    min_af: 0.2,
    min_gq: 20.0,
    mode: BcfCallMode::Multiallelic,
    sample_name: "SAMPLE",
    write_pl: true,
  }
}

///|
/// Compute genotype quality GQ and PL scores for a REF/ALT site.
///
/// Returns (gt_index, gq, [pl_rr, pl_ra, pl_aa]).
/// gt_index: 0 = 0/0, 1 = 0/1, 2 = 1/1.
pub fn call_genotype_pl(
  ref_count : Int32,
  alt_count : Int32,
) -> (Int32, Double, Array[String]) {
  let total = ref_count + alt_count
  if total == 0 {
    return (0, 0.0, ["99", "99", "99"])
  }
  let alt_f = alt_count.to_double() / total.to_double()
  // Diploid prior: 0.25 / 0.5 / 0.25
  let prior_rr : Double = 0.25
  let prior_ra : Double = 0.5
  let prior_aa : Double = 0.25
  // Error rate from base quality assumption (Phred 30 = 0.001)
  let err : Double = 0.001
  // Likelihoods:  binomial-like
  let mut lik_rr : Double = 1.0
  let mut lik_ra : Double = 1.0
  let mut lik_aa : Double = 1.0
  let mut i = 0
  while i < ref_count {
    lik_rr = lik_rr * (1.0 - err)
    lik_ra = lik_ra * 0.5
    lik_aa = lik_aa * err
    i = i + 1
  }
  i = 0
  while i < alt_count {
    lik_rr = lik_rr * err
    lik_ra = lik_ra * 0.5
    lik_aa = lik_aa * (1.0 - err)
    i = i + 1
  }
  lik_rr = lik_rr * prior_rr
  lik_ra = lik_ra * prior_ra
  lik_aa = lik_aa * prior_aa
  // Normalize
  let sum = lik_rr + lik_ra + lik_aa
  let p_rr = lik_rr / sum
  let p_ra = lik_ra / sum
  let p_aa = lik_aa / sum
  // Call best genotype
  let mut gt : Int32 = 0
  let mut gt_prob : Double = p_rr
  if p_ra > gt_prob {
    gt = 1
    gt_prob = p_ra
  }
  if p_aa > gt_prob {
    gt = 2
    gt_prob = p_aa
  }
  // GQ = -10 * log10(1 - posterior)
  let one_minus : Double = 1.0 - gt_prob
  let tiny : Double = 0.000000000000000000000000000001
  let log_arg : Double = if one_minus < tiny { tiny } else { one_minus }
  let gq : Double = -4.34294481903 * @math.ln(log_arg)
  // PL scores: -10 * log10(lik), normalized to best
  let log10_factor : Double = -4.34294481903
  let lr : Double = if lik_rr < tiny { tiny } else { lik_rr }
  let la : Double = if lik_ra < tiny { tiny } else { lik_ra }
  let laa : Double = if lik_aa < tiny { tiny } else { lik_aa }
  let pl_rr_d : Double = log10_factor * @math.ln(lr)
  let pl_ra_d : Double = log10_factor * @math.ln(la)
  let pl_aa_d : Double = log10_factor * @math.ln(laa)
  let best_pl : Double = if pl_rr_d < pl_ra_d {
    if pl_rr_d < pl_aa_d {
      pl_rr_d
    } else {
      pl_aa_d
    }
  } else if pl_ra_d < pl_aa_d {
    pl_ra_d
  } else {
    pl_aa_d
  }
  let pl_rr_s = (pl_rr_d - best_pl).to_int().to_string()
  let pl_ra_s = (pl_ra_d - best_pl).to_int().to_string()
  let pl_aa_s = (pl_aa_d - best_pl).to_int().to_string()
  (gt, gq, [pl_rr_s, pl_ra_s, pl_aa_s])
}

///|
/// Genotype index → diploid GT string ("0/0", "0/1", "1/1").
pub fn gt_str(gt_index : Int32) -> String {
  match gt_index {
    0 => "0/0"
    1 => "0/1"
    2 => "1/1"
    _ => "./."
  }
}

///|
/// Build INFO/FORMAT meta header lines suitable for bcftools-call output.
pub fn bcf_call_meta_headers() -> Array[String] {
  [
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=Total Depth>", "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele Frequency>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>", "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Normalized Genotype Likelihoods>",
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype Quality>", "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read Depth>",
  ]
}
