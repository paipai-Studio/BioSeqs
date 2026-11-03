///|
/// qtlseq_demo — QTL-seq / MutMap 风格 BSA 混池分离定位演示:
/// 合成双混池 VCF(注入峰值区间) -> 逐位点 SNP-index -> 滑窗
/// ΔSNP-index -> 候选区间 -> TSV/BED/JSON.

///|
fn qsd_vcf() -> String {
  let sb = StringBuilder()
  sb.write_string("##fileformat=VCFv4.2\n")
  sb.write_string(
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMUT\tWT\n",
  )
  let bases = ["A", "C", "G", "T"]
  let mut pos : Int32 = 500
  let mut i : Int32 = 0
  while pos <= 20500 {
    let bi : Int32 = i % 4
    let alt : String = if bases[bi] == "A" {
      "G"
    } else if bases[bi] == "C" {
      "T"
    } else if bases[bi] == "G" {
      "A"
    } else {
      "C"
    }
    // 峰值区间 [10000, 12500): mut 0.9 / wt 0.1; 其余中性 0.5/0.5
    let mut_ad : Int32 = if pos >= 10000 && pos < 12500 { 2 } else { 10 }
    let mut_alt : Int32 = if pos >= 10000 && pos < 12500 { 18 } else { 10 }
    let wt_ad : Int32 = if pos >= 10000 && pos < 12500 { 18 } else { 10 }
    let wt_alt : Int32 = if pos >= 10000 && pos < 12500 { 2 } else { 10 }
    sb.write_string(
      "chr1\t\{pos}\t.\t\{bases[bi]}\t\{alt}\t.\tPASS\t.\tGT:AD\t0/1:\{mut_ad},\{mut_alt}\t0/1:\{wt_ad},\{wt_alt}\n",
    )
    pos = pos + 500
    i = i + 1
  }
  sb.to_string()
}

///|
fn main {
  println("1. QTL-seq / MutMap 风格 BSA 混池分离定位演示")
  let vcf = @src.parse_vcf(qsd_vcf())
  println(
    "2. 双混池 VCF: \{vcf.records.length()} 双等位 SNP(MUT + WT 池, 注入峰值区间 [10000, 12500))",
  )
  // 3. SNP-index
  let params : @src.QtlseqParams = {
    min_depth: 5,
    window: 2000,
    step: 500,
    delta_threshold: 0.3,
  }
  let sites = @src.qtlseq_snp_index(vcf.records, "MUT", wt_bulk="WT", params~)
  println("3. 逐位点 SNP-index: \{sites.length()} 位点(深度过滤后)")
  // 4. 滑窗
  let windows = @src.qtlseq_windows(sites, params~)
  println(
    "4. 滑窗聚合(\{params.window}bp 窗 / \{params.step}bp 步): \{windows.length()} 窗",
  )
  // 展示峰值附近窗口
  let mut shown = 0
  for w in windows {
    if w.start >= 9000 && w.start <= 13000 && shown < 6 {
      println(
        "  \{w.start}-\{w.end_}: n=\{w.n_snps}, mut=\{w.mean_mut_index}, wt=\{w.mean_wt_index}, delta=\{w.mean_delta}",
      )
      shown = shown + 1
    }
  }
  // 5. 候选区间
  let intervals = @src.qtlseq_intervals(windows, params~)
  println("5. 候选区间(|delta| >= \{params.delta_threshold}):")
  for iv in intervals {
    println(
      "  \{iv.chrom}:\{iv.start}-\{iv.end_} 峰值=\{iv.peak_delta} 方向=\{iv.direction} SNPs=\{iv.n_snps}",
    )
  }
  // 6. 导出
  let tsv = @src.qtlseq_windows_tsv(windows)
  let bed = @src.qtlseq_intervals_bed(intervals)
  let json = @src.qtlseq_json(sites, windows, intervals)
  println(
    "6. TSV \{tsv.length()} 字节, BED \{bed.length()} 字节, JSON \{json.length()} 字节",
  )
  println(@src.qtlseq_summary_line(sites, windows, intervals))
  println(
    "演示完成: 双池 VCF -> SNP-index -> ΔSNP-index 滑窗 -> 候选区间 BED 闭环",
  )
}
