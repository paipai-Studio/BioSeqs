///|
/// fusion.mbt — 融合基因 caller(split-read 断点簇聚类, STAR-Fusion 风格).
///
/// 对齐 STAR-Fusion / Arriba 的核心检出语义(本地实现):
/// - split-read 提取: SAM SA 标签(supplementary alignment)解析,
///   主定位 + 每条补充分割比对构成断点对(junction read)
/// - 断点簇聚类: (chromA, chromB) 有序对内断点位置邻近(容差 bp)聚类
/// - 候选判定: 支持读数 >= min_support, 跨基因/跨染色体或大跨距
/// - 过滤: 同链 read-through(同向且距离 < readthrough_max)
/// - 输出: STAR-Fusion 风格 TSV + BEDPE 记录(与 bedpe.mbt 闭环)
///
/// 断点近似: pos + 比对长度(CIGAR M/D 段和), 忽略 soft-clip 精修.
/// 数据不出本机.

///|
/// 一条 split-read 断点对(来自主定位 + 补充定位).
pub(all) struct FusionJunction {
  qname : String
  chrom_a : String
  break_a : Int32
  strand_a : Int32
  chrom_b : String
  break_b : Int32
  strand_b : Int32
}

///|
/// 一个融合候选(断点簇).
pub(all) struct FusionCandidate {
  chrom_a : String
  break_a : Int32
  chrom_b : String
  break_b : Int32
  strand_a : Int32
  strand_b : Int32
  support : Int32
  /// 参与读段名
  read_names : Array[String]
  is_interchromosomal : Bool
  fusion_name : String
  span : Int32
  filter : String
}

///|
/// 从 SAM CIGAR 计算比对跨度(M/D 段和).
fn fusion_cigar_span(cigar : Array[CigarElement]) -> Int32 {
  let mut span = 0
  for el in cigar {
    match el.op {
      CigarOp::Match | CigarOp::Equal | CigarOp::Diff | CigarOp::Delete =>
        span = span + el.length
      _ => ()
    }
  }
  span
}

///|
/// 解析 SA 标签: "rname,pos,strand,CIGAR,mapQ,NM;..." -> 每条补充分割比对.
fn fusion_parse_sa(sa_value : String) -> Array[(String, Int32, Int32, String)] {
  let entries : Array[(String, Int32, Int32, String)] = Array(capacity=0)
  let parts = sa_value.split(";").to_array()
  for part in parts {
    let p = part.to_owned()
    if p.length() == 0 {
      continue
    }
    let fields = p.split(",").to_array()
    if fields.length() < 4 {
      continue
    }
    let rname = fields[0].to_owned()
    let pos = parse_int(fields[1].to_owned())
    let strand = if fields[2].to_owned() == "+" { 0 } else { 1 }
    let cigar = fields[3].to_owned()
    entries.push((rname, pos, strand, cigar))
  }
  entries
}

///|
/// 从 SAM 文本提取 split-read 断点对.
/// 主定位 + 每条 SA 补充定位构成一个 junction.
pub fn fusion_junctions_from_sam(content : String) -> Array[FusionJunction] {
  let sam = parse_sam(content)
  let out : Array[FusionJunction] = Array(capacity=0)
  for rec in sam.records {
    let sa = rec.tags.get("SA")
    match sa {
      None => ()
      Some(sa_value) => {
        let supp = fusion_parse_sa(sa_value)
        for entry in supp {
          // 主定位断点: pos + span(5' 方向近似)
          let main_span = fusion_cigar_span(rec.cigar)
          let main_break = rec.pos + main_span
          let supp_break = entry.1 + fusion_cigar_span_str(entry.3)
          // 有序对: 按字典序排 chrom, 保证 (A,B) 与 (B,A) 同簇
          let j = if rec.rname <= entry.0 {
            {
              qname: rec.qname,
              chrom_a: rec.rname,
              break_a: main_break,
              strand_a: if rec.flag.land(16) != 0 {
                1
              } else {
                0
              },
              chrom_b: entry.0,
              break_b: supp_break,
              strand_b: entry.2,
            }
          } else {
            {
              qname: rec.qname,
              chrom_a: entry.0,
              break_a: supp_break,
              strand_a: entry.2,
              chrom_b: rec.rname,
              break_b: main_break,
              strand_b: if rec.flag.land(16) != 0 {
                1
              } else {
                0
              },
            }
          }
          out.push(j)
        }
      }
    }
  }
  out
}

///|
/// CIGAR 字符串解析的跨度(轻量, 避免 SAM CIGAR 结构依赖).
fn fusion_cigar_span_str(cigar : String) -> Int32 {
  let mut span = 0
  let mut num = 0
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar.unsafe_get(i).unsafe_to_char()
    if c >= '0' && c <= '9' {
      num = num * 10 + (c.to_int() - '0'.to_int())
    } else {
      match c {
        'M' | 'D' | '=' | 'X' => span = span + num
        _ => ()
      }
      num = 0
    }
  }
  span
}

///|
/// 融合检出: 断点簇聚类 + 支持数过滤.
/// cluster_bp: 断点聚类容差(默认 50bp); min_support: 最小支持读数;
/// readthrough_max: 同链同染色体且跨距小于该值视为 read-through 过滤.
pub fn fusion_call(
  junctions : Array[FusionJunction],
  cluster_bp~ : Int32 = 50,
  min_support~ : Int32 = 2,
  readthrough_max~ : Int32 = 200000,
) -> Array[FusionCandidate] {
  // 按 (chromA, chromB, strandA, strandB) 分组(入口先归一有序对)
  let groups : Map[String, Array[FusionJunction]] = Map([], capacity=0)
  for j0 in junctions {
    // 若 chrom_a > chrom_b 则交换, 保证 (A,B) 与 (B,A) 同簇
    let j : FusionJunction = if j0.chrom_a > j0.chrom_b {
      {
        qname: j0.qname,
        chrom_a: j0.chrom_b,
        break_a: j0.break_b,
        strand_a: j0.strand_b,
        chrom_b: j0.chrom_a,
        break_b: j0.break_a,
        strand_b: j0.strand_a,
      }
    } else {
      j0
    }
    let key = j.chrom_a +
      "|" +
      j.chrom_b +
      "|" +
      j.strand_a.to_string() +
      "|" +
      j.strand_b.to_string()
    match groups.get(key) {
      Some(arr) => arr.push(j)
      None => {
        let arr : Array[FusionJunction] = Array(capacity=0)
        arr.push(j)
        groups.set(key, arr)
      }
    }
  }
  let candidates : Array[FusionCandidate] = Array(capacity=0)
  for key in groups.keys() {
    let js = groups.get(key).unwrap()
    if js.length() == 0 {
      continue
    }
    // 断点对聚类: 先按 break_a 排序, 滑动合并 break_a 与 break_b 均在容差内
    let sorted : Array[FusionJunction] = Array(capacity=0)
    for j in js {
      sorted.push(j)
    }
    sorted.sort_by(fn(a, b) { a.break_a - b.break_a })
    let clusters : Array[Array[FusionJunction]] = Array(capacity=0)
    let mut current : Array[FusionJunction] = Array(capacity=0)
    current.push(sorted[0])
    for i = 1; i < sorted.length(); i = i + 1 {
      let prev = current[current.length() - 1]
      let j = sorted[i]
      let da = j.break_a - prev.break_a
      let db = if j.break_b > prev.break_b {
        j.break_b - prev.break_b
      } else {
        prev.break_b - j.break_b
      }
      let da_abs = if da >= 0 { da } else { -da }
      if da_abs <= cluster_bp && db <= cluster_bp {
        current.push(j)
      } else {
        clusters.push(current)
        current = Array(capacity=0)
        current.push(j)
      }
    }
    clusters.push(current)
    // 每簇产出候选
    for cluster in clusters {
      if cluster.length() < min_support {
        continue
      }
      // 簇代表断点: 中位数
      let ba : Array[Int32] = Array(capacity=0)
      let bb : Array[Int32] = Array(capacity=0)
      let names : Array[String] = Array(capacity=0)
      for j in cluster {
        ba.push(j.break_a)
        bb.push(j.break_b)
        names.push(j.qname)
      }
      ba.sort()
      bb.sort()
      let med_a = ba[ba.length() / 2]
      let med_b = bb[bb.length() / 2]
      let chrom_a = cluster[0].chrom_a
      let chrom_b = cluster[0].chrom_b
      let strand_a = cluster[0].strand_a
      let strand_b = cluster[0].strand_b
      let inter = chrom_a != chrom_b
      let span = if inter {
        -1
      } else {
        let s = med_b - med_a
        if s < 0 {
          -s
        } else {
          s
        }
      }
      // read-through 过滤: 同链同染色体近距
      let mut filter = "PASS"
      if !inter && strand_a == strand_b && span < readthrough_max {
        filter = "readthrough"
      }
      candidates.push({
        chrom_a,
        break_a: med_a,
        chrom_b,
        break_b: med_b,
        strand_a,
        strand_b,
        support: cluster.length(),
        read_names: names,
        is_interchromosomal: inter,
        fusion_name: chrom_a + "--" + chrom_b,
        span,
        filter,
      })
    }
  }
  // 排序: 支持数降序
  candidates.sort_by(fn(a, b) {
    if a.support != b.support {
      b.support - a.support
    } else {
      compare(a.fusion_name, b.fusion_name)
    }
  })
  candidates
}

///|
/// 候选摘要行.
pub fn FusionCandidate::summary(self : FusionCandidate) -> String {
  "\{self.fusion_name} support=\{self.support} \{self.loc()} filter=\{self.filter}"
}

///|
/// 位置描述.
pub fn FusionCandidate::loc(self : FusionCandidate) -> String {
  if self.is_interchromosomal {
    "\{self.chrom_a}:\{self.break_a}--\{self.chrom_b}:\{self.break_b} (interchromosomal)"
  } else {
    "\{self.chrom_a}:\{self.break_a}--\{self.chrom_b}:\{self.break_b} (span \{self.span}bp)"
  }
}

///|
/// STAR-Fusion 风格 TSV 输出.
pub fn fusion_write_tsv(candidates : Array[FusionCandidate]) -> String {
  let sb = StringBuilder()
  sb.write_string(
    "#FusionName\tJunctionReadCount\tChromA\tBreakA\tStrandA\tChromB\tBreakB\tStrandB\tInterchromosomal\tSpan\tFilter\n",
  )
  for c in candidates {
    sb.write_string(c.fusion_name)
    sb.write_char('\t')
    sb.write_string(c.support.to_string())
    sb.write_char('\t')
    sb.write_string(c.chrom_a)
    sb.write_char('\t')
    sb.write_string(c.break_a.to_string())
    sb.write_char('\t')
    sb.write_string(if c.strand_a == 0 { "+" } else { "-" })
    sb.write_char('\t')
    sb.write_string(c.chrom_b)
    sb.write_char('\t')
    sb.write_string(c.break_b.to_string())
    sb.write_char('\t')
    sb.write_string(if c.strand_b == 0 { "+" } else { "-" })
    sb.write_char('\t')
    sb.write_string(if c.is_interchromosomal { "Y" } else { "N" })
    sb.write_char('\t')
    sb.write_string(c.span.to_string())
    sb.write_char('\t')
    sb.write_string(c.filter)
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// BEDPE 输出(与 bedpe.mbt 的解析闭环).
pub fn fusion_write_bedpe(candidates : Array[FusionCandidate]) -> String {
  let sb = StringBuilder()
  sb.write_string(
    "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\n",
  )
  for c in candidates {
    sb.write_string(c.chrom_a)
    sb.write_char('\t')
    sb.write_string((c.break_a - 1).to_string())
    sb.write_char('\t')
    sb.write_string(c.break_a.to_string())
    sb.write_char('\t')
    sb.write_string(c.chrom_b)
    sb.write_char('\t')
    sb.write_string((c.break_b - 1).to_string())
    sb.write_char('\t')
    sb.write_string(c.break_b.to_string())
    sb.write_char('\t')
    sb.write_string(c.fusion_name)
    sb.write_char('\t')
    sb.write_string(c.support.to_string())
    sb.write_char('\t')
    sb.write_string(if c.strand_a == 0 { "+" } else { "-" })
    sb.write_char('\t')
    sb.write_string(if c.strand_b == 0 { "+" } else { "-" })
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// 检出结果摘要.
pub fn fusion_candidates_summary(candidates : Array[FusionCandidate]) -> String {
  let mut n_pass = 0
  let mut n_inter = 0
  for c in candidates {
    if c.filter == "PASS" {
      n_pass = n_pass + 1
    }
    if c.is_interchromosomal {
      n_inter = n_inter + 1
    }
  }
  "fusion: \{candidates.length()} candidates (\{n_pass} PASS, \{n_inter} interchromosomal)"
}
