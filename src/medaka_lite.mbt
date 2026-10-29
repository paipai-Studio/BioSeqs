///|
/// medaka_lite.mbt — 长读长变异检出(高错误率共识模型, medaka 风格).
///
/// 对齐 medaka 的共识变异检出语义(简化本地实现):
/// - 长读错误模型: ONT ~10% 一致性(epsilon=0.10), indel 偏置
/// - 逐位点多数投票 + 长读自适应基因型似然(epsilon 放宽)
/// - active region 扫描 + 插入/缺失事件多数投票
///
/// 与 minimap_align_preset(map-ont) 组合完成"长读比对 -> 变异检出"
/// 本地闭环; 输出复用 GatkVariant 与 gatk_write_vcf 进入 VCF 栈.
///
/// 已知限制: 逐位点共识模型不做读段级 indel 感知比对, 单条读段的
/// 缺失会在其下游产生连锁错位假阳性; 适用单变异位点/低密度场景,
/// 高 indet 密度建议分 active region 处理.
pub(all) struct MedakaRead {
  name : String
  seq : String
  pos : Int32
}

///|
/// ONT 长读错误率(共识一致性 ~90%).
let medaka_ont_epsilon : Double = 0.10

///|
/// HiFi 长读错误率(一致性 ~99.9%).
let medaka_hifi_epsilon : Double = 0.005

///|
/// 长读错误模型.
pub(all) enum MedakaModel {
  Ont
  Hifi
} derive(Eq, Debug)

///|
pub extend MedakaModel with Eq::{equal, not_equal}

///|
pub extend MedakaModel with @moonbitlang/core/debug.Debug::{to_repr}

///|
fn medaka_epsilon(model : MedakaModel) -> Double {
  match model {
    Ont => medaka_ont_epsilon
    Hifi => medaka_hifi_epsilon
  }
}

///|
fn medaka_to_gatk_reads(reads : Array[MedakaRead]) -> Array[GatkRead] {
  let out : Array[GatkRead] = Array(capacity=0)
  for r in reads {
    out.push({ name: r.name, seq: r.seq, pos: r.pos, })
  }
  out
}

///|
/// 长读变异检出: 高错误率共识模型 + 多数投票 indel.
/// 与 gatk-lite 的差异: epsilon 放宽(长读错误模型), min_depth 降低
/// (长读覆盖深但条数少), indel 证据容忍度提高.
pub fn medaka_call_variants(
  chrom : String,
  reference : String,
  reads : Array[MedakaRead],
  model~ : MedakaModel = Ont,
  min_depth~ : Int32 = 3,
  min_af~ : Double = 0.3,
  min_gq~ : Double = 10.0,
) -> Array[GatkVariant] {
  let gatk_reads = medaka_to_gatk_reads(reads)
  let eps = medaka_epsilon(model)
  let variants : Array[GatkVariant] = Array(capacity=0)
  let n = reference.length()
  for p = 0; p < n; p = p + 1 {
    let site = gatk_pileup_site(reference, gatk_reads, p)
    if site.depth < min_depth {
      continue
    }
    let mut alt : Char = 'N'
    let mut alt_count = 0
    for pair in site.base_counts {
      if pair.0 != site.ref_base {
        alt = pair.0
        alt_count = pair.1
        break
      }
    }
    let af = alt_count.to_double() / site.depth.to_double()
    if alt == 'N' || af < min_af {
      continue
    }
    let gts = gatk_genotype_likelihoods(site, alt, epsilon=eps)
    let best = gts[0]
    if best.gt == "0/0" {
      continue
    }
    if best.log10_error < min_gq {
      continue
    }
    let ref_count = site.depth - alt_count
    variants.push({
      chrom,
      pos: p + 1,
      ref: site.ref_base.to_string(),
      alt: alt.to_string(),
      qual: best.log10_error,
      filter: "PASS",
      gt: best.gt,
      dp: site.depth,
      ad_ref: ref_count,
      ad_alt: alt_count,
      gq: best.log10_error.to_int(),
      af,
    })
  }
  // active region indel(长读 indel 偏置: 错配密度更宽容)
  let regions = gatk_active_regions(
    reference,
    gatk_reads,
    window=25,
    min_mismatches=2,
  )
  for region in regions {
    let (rs, re) = region
    let mut best_ins : (Int32, String, Int32, Int32) = (-1, "", 0, 0)
    let mut best_del : (Int32, Int32, Int32, Int32) = (-1, 0, 0, 0)
    for p = rs; p < re && p < n; p = p + 1 {
      let site = gatk_pileup_site(reference, gatk_reads, p)
      if site.depth < min_depth {
        continue
      }
      if site.insertions.length() > 0 {
        let (ins_seq, cnt) = site.insertions[0]
        if cnt > best_ins.2 {
          best_ins = (p, ins_seq, cnt, site.depth)
        }
      }
      if site.deletions.length() > 0 {
        let (dl, cnt) = site.deletions[0]
        if cnt > best_del.2 {
          best_del = (p, dl, cnt, site.depth)
        }
      }
    }
    if best_del.2 >= best_ins.2 && best_del.2 > 0 {
      let (p, dl, cnt, depth) = best_del
      let af = cnt.to_double() / depth.to_double()
      if af >= min_af && p >= 0 && p + 1 + dl <= n {
        let deleted = reference[p:p + 1 + dl].to_owned()
        let kept = reference[p:p + 1].to_owned()
        variants.push({
          chrom,
          pos: p + 1,
          ref: deleted,
          alt: kept,
          qual: 20.0,
          filter: "PASS",
          gt: "0/1",
          dp: depth,
          ad_ref: depth - cnt,
          ad_alt: cnt,
          gq: 20,
          af,
        })
      }
    } else if best_ins.2 > 0 {
      let (p, ins_seq, cnt, depth) = best_ins
      let af = cnt.to_double() / depth.to_double()
      if af >= min_af && p >= 0 && p < n {
        let anchor = reference[p:p + 1].to_owned()
        variants.push({
          chrom,
          pos: p + 1,
          ref: anchor,
          alt: anchor + ins_seq,
          qual: 20.0,
          filter: "PASS",
          gt: "0/1",
          dp: depth,
          ad_ref: depth - cnt,
          ad_alt: cnt,
          gq: 20,
          af,
        })
      }
    }
  }
  variants.sort_by(fn(a, b) { a.pos - b.pos })
  variants
}

///|
/// 长读检出摘要.
pub fn medaka_variants_summary(variants : Array[GatkVariant]) -> String {
  let mut n_snp = 0
  let mut n_indel = 0
  for v in variants {
    if v.ref.length() == 1 && v.alt.length() == 1 {
      n_snp = n_snp + 1
    } else {
      n_indel = n_indel + 1
    }
  }
  "medaka-lite: \{variants.length()} variants (\{n_snp} SNP, \{n_indel} indels)"
}
