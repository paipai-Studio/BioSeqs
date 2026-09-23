///|
/// scikit-allel `allel.io` bridge — convert a parsed VCF (`VcfFile`) into the
/// genotype / haplotype / allele-count arrays consumed by the `skallel`
/// statistics (LD, PCA, iHS/nSL, Fst, …).

///|
/// Parse a GT field string ("0/1", "1|0", "./.", "2", …) into allele indices
/// (missing calls become −1).
pub fn skallel_parse_gt(s : String) -> Array[Int] {
  let out : Array[Int] = Array(capacity=0)
  let mut cur = StringBuilder::new()
  let mut i = 0
  while i <= s.length() {
    let is_sep = if i == s.length() {
      true
    } else {
      let c = s.unsafe_get(i).unsafe_to_char()
      c == '/' || c == '|'
    }
    if is_sep {
      let token = cur.to_string()
      if token.length() > 0 {
        if token == "." {
          out.push(-1)
        } else {
          out.push(parse_int(token))
        }
      }
      cur = StringBuilder::new()
    } else {
      cur.write_char(s.unsafe_get(i).unsafe_to_char())
    }
    i = i + 1
  }
  out
}

///|
fn skallel_gt_index(format : Array[String]) -> Int {
  for i = 0; i < format.length(); i = i + 1 {
    if format[i] == "GT" {
      return i
    }
  }
  0
}

///|
/// Extract the GT FORMAT field as a genotype array of shape
/// `(n_variants × n_samples × ploidy)`, with missing alleles set to −1.
pub fn skallel_vcf_genotypes(
  vcf : VcfFile,
  sample_names : Array[String],
) -> Array[Array[Array[Int]]] {
  let out : Array[Array[Array[Int]]] = Array(capacity=0)
  for rec in vcf.records {
    let gt_idx = match rec.format {
      None => 0
      Some(fmt) => skallel_gt_index(fmt)
    }
    let row : Array[Array[Int]] = Array(capacity=0)
    for name in sample_names {
      match rec.samples.get(name) {
        None => row.push([-1, -1])
        Some(data) => {
          let gt_str = if gt_idx < data.length() {
            data[gt_idx]
          } else {
            "."
          }
          row.push(skallel_parse_gt(gt_str))
        }
      }
    }
    out.push(row)
  }
  out
}

///|
/// Transform a genotype array into the number of non-reference alleles per
/// call (0/1/2 dosage matrix), matching `GenotypeArray.to_n_alt()`. Missing
/// calls become `fill` (default 0, the same as homozygous reference).
pub fn skallel_gt_to_n_alt(
  g : Array[Array[Array[Int]]],
  fill? : Int = 0,
) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = Array(capacity=0)
  for row in g {
    let out_row : Array[Int] = Array(capacity=0)
    for call in row {
      let mut n_alt = 0
      let mut has_missing = false
      for a in call {
        if a < 0 {
          has_missing = true
        } else if a > 0 {
          n_alt = n_alt + 1
        }
      }
      out_row.push(if fill != 0 && has_missing { fill } else { n_alt })
    }
    out.push(out_row)
  }
  out
}

///|
/// Reshape a genotype array to haplotypes (`GenotypeArray.to_haplotypes()`): a
/// `(n_variants × n_samples·ploidy)` array flattening each sample's alleles.
pub fn skallel_gt_to_haplotypes(g : Array[Array[Array[Int]]]) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = Array(capacity=0)
  for row in g {
    let out_row : Array[Int] = Array(capacity=0)
    for call in row {
      for a in call {
        out_row.push(a)
      }
    }
    out.push(out_row)
  }
  out
}

///|
/// Count alleles (0 … `max_allele`) across all samples of a genotype array,
/// returning a `(n_variants × (max_allele+1))` matrix.
pub fn skallel_gt_count_alleles(
  g : Array[Array[Array[Int]]],
  max_allele : Int,
) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = Array(capacity=0)
  for row in g {
    let counts = Array::make(max_allele + 1, 0)
    for call in row {
      for a in call {
        if a >= 0 && a <= max_allele {
          counts[a] = counts[a] + 1
        }
      }
    }
    out.push(counts)
  }
  out
}

///|
/// Count alleles across a subset of samples (given by 0-based sample indices),
/// returning a `(n_variants × (max_allele+1))` matrix.
pub fn skallel_gt_count_alleles_subpop(
  g : Array[Array[Array[Int]]],
  sample_indices : Array[Int],
  max_allele : Int,
) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = Array(capacity=0)
  for row in g {
    let counts = Array::make(max_allele + 1, 0)
    for si in sample_indices {
      if si >= 0 && si < row.length() {
        for a in row[si] {
          if a >= 0 && a <= max_allele {
            counts[a] = counts[a] + 1
          }
        }
      }
    }
    out.push(counts)
  }
  out
}