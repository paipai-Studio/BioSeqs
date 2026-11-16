// ============================================================================
// pysam_vcf_serialize — VcfFile → VCF 文本 序列化（支持 round-trip）
//
// 核心函数：vcf_to_text(vcf) → String
// 反向函数：vcf_from_text = parse_vcf（已存在于 vcf.mbt）
//
// round-trip 验证：vcf_from_text(vcf_to_text(vcf)) 应该等于原 vcf
// ============================================================================

///|
fn tags_to_string(tags : Map[String, String]) -> String {
  // Map 没有固定顺序，我们按键名排序后输出以保证确定性
  let keys = tags.keys()
  // 简单插入排序
  let ks : Array[String] = Array(capacity=0)
  for k in keys {
    ks.push(k)
  }
  let n = ks.length()
  let mut p = 0
  while p < n {
    let mut q = p + 1
    while q < n {
      if ks[p] > ks[q] {
        let tmp = ks[p]
        ks[p] = ks[q]
        ks[q] = tmp
      }
      q = q + 1
    }
    p = p + 1
  }
  let mut out = ""
  let mut i = 0
  while i < ks.length() {
    if i > 0 {
      out = out + ","
    }
    let v = tags.get(ks[i])
    match v {
      Some(val) => out = out + ks[i] + "=" + val
      None => out = out + ks[i]
    }
    i = i + 1
  }
  out
}

///|
fn _join_str(arr : Array[String], sep : String) -> String {
  let mut out = ""
  let mut i = 0
  while i < arr.length() {
    if i > 0 {
      out = out + sep
    }
    out = out + arr[i]
    i = i + 1
  }
  out
}

///|
fn alt_to_string(alt : Array[String]) -> String {
  if alt.length() == 0 {
    return "."
  }
  _join_str(alt, ",")
}

///|
fn filter_to_string(f : Array[String]) -> String {
  if f.length() == 0 {
    return "PASS"
  }
  _join_str(f, ";")
}

///|
fn id_to_string(id : Array[String]) -> String {
  if id.length() == 0 {
    return "."
  }
  _join_str(id, ";")
}

///|
fn qual_to_string(q : Double?) -> String {
  match q {
    Some(v) => {
      // 如果是整数值，输出不带小数点
      let rounded : Double = v.to_int().to_double()
      if v == rounded {
        v.to_int().to_string()
      } else {
        v.to_string()
      }
    }
    None => "."
  }
}

///|
fn info_to_string(info : Map[String, String]) -> String {
  if info.size() == 0 {
    return "."
  }
  // Map entries key=value; 按键排序
  let keys = info.keys()
  let ks : Array[String] = Array(capacity=0)
  for k in keys {
    ks.push(k)
  }
  let n = ks.length()
  let mut p = 0
  while p < n {
    let mut q = p + 1
    while q < n {
      if ks[p] > ks[q] {
        let tmp = ks[p]
        ks[p] = ks[q]
        ks[q] = tmp
      }
      q = q + 1
    }
    p = p + 1
  }
  let mut out = ""
  let mut i = 0
  while i < ks.length() {
    if i > 0 {
      out = out + ";"
    }
    match info.get(ks[i]) {
      Some(val) =>
        // 值为空时只输出 key
        if val.length() == 0 {
          out = out + ks[i]
        } else {
          out = out + ks[i] + "=" + val
        }
      None => out = out + ks[i]
    }
    i = i + 1
  }
  out
}

///|
/// 序列化一条 VcfRecord。`columns` 来自 fixed_header 决定哪些列存在。
pub fn vcf_record_to_line(rec : VcfRecord, columns : Array[String]) -> String {
  let samples_present = columns.length() > 8 // FORMAT + 至少一个 SAMPLE
  // 1-8 列固定
  let mut out = rec.chrom +
    "\t" +
    rec.pos.to_string() +
    "\t" +
    id_to_string(rec.id) +
    "\t" +
    rec.ref_seq +
    "\t" +
    alt_to_string(rec.alt) +
    "\t" +
    qual_to_string(rec.qual) +
    "\t" +
    filter_to_string(rec.filter) +
    "\t" +
    info_to_string(rec.info)
  // FORMAT + SAMPLES
  if samples_present {
    out = out + "\t"
    match rec.format {
      Some(fmt) => out = out + _join_str(fmt, ":")
      None => out = out + "."
    }
    // SAMPLE columns
    let fmt_keys = columns
    let mut fmt_idx_val : Int32 = -1
    let mut _i = 0
    while _i < fmt_keys.length() {
      if fmt_keys[_i] == "FORMAT" {
        fmt_idx_val = _i
        break
      }
      _i = _i + 1
    }
    let fmt_idx : Int32? = if fmt_idx_val >= 0 {
      Some(fmt_idx_val)
    } else {
      None
    }
    let sample_names : Array[String] = match fmt_idx {
      Some(idx) => {
        let out2 : Array[String] = Array(capacity=fmt_keys.length() - idx - 1)
        let mut si = idx + 1
        while si < fmt_keys.length() {
          out2.push(fmt_keys[si])
          si = si + 1
        }
        out2
      }
      None => []
    }
    let mut i = 0
    while i < sample_names.length() {
      out = out + "\t"
      match rec.samples.get(sample_names[i]) {
        Some(arr) => out = out + _join_str(arr, ":")
        None => out = out + "."
      }
      i = i + 1
    }
  }
  out
}

///|
/// 把 VcfFile 完整序列化为 VCF 文本字符串。
pub fn vcf_to_text(vcf : VcfFile) -> String {
  let mut out = ""
  // meta headers
  let mut mi = 0
  while mi < vcf.meta_headers.length() {
    let mh = vcf.meta_headers[mi]
    out = out + "##" + mh.key + "=" + mh.value + "\n"
    mi = mi + 1
  }
  // fixed header
  let columns = match vcf.fixed_header {
    Some(fh) => fh.columns
    None => ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
  }
  out = out + "#" + _join_str(columns, "\t") + "\n"
  // records
  let mut ri = 0
  while ri < vcf.records.length() {
    out = out + vcf_record_to_line(vcf.records[ri], columns) + "\n"
    ri = ri + 1
  }
  out
}

///|
/// Convenience: write the serialized VCF to a String. Alias of `vcf_to_text`.
pub fn write_vcf(vcf : VcfFile) -> String {
  vcf_to_text(vcf)
}
