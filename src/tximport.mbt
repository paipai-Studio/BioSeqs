///|
pub struct TxImportData {
  abundance : Array[Array[Double]]
  counts : Array[Array[Double]]
  length : Array[Array[Double]]
  tx_names : Array[String]
  sample_names : Array[String]
  gene_names : Array[String]
}

///|
pub struct Tx2Gene {
  tx_id : String
  gene_id : String
  gene_name : String
}

///|
pub fn tx_import_data_new(
  abundance : Array[Array[Double]],
  counts : Array[Array[Double]],
  length : Array[Array[Double]],
  tx_names : Array[String],
  sample_names : Array[String],
  gene_names : Array[String],
) -> TxImportData {
  TxImportData::{
    abundance,
    counts,
    length,
    tx_names,
    sample_names,
    gene_names,
  }
}

///|
pub fn tx2gene_new(
  tx_id : String,
  gene_id : String,
  gene_name : String,
) -> Tx2Gene {
  Tx2Gene::{ tx_id, gene_id, gene_name }
}

///|
pub fn parse_quant_sf(
  content : String,
) -> (Array[String], Array[Double], Array[Double], Array[Double]) {
  let lines = split_by_char(content, '\n'.to_int())

  let tx_names : Array[String] = Array::new()
  let abundance : Array[Double] = Array::new()
  let counts : Array[Double] = Array::new()
  let length : Array[Double] = Array::new()

  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    if line.length() == 0 || line.starts_with("#") {
      i = i + 1
      continue
    }

    if i == 0 {
      i = i + 1
      continue
    }

    let parts = split_by_char(line, 9)
    if parts.length() >= 5 {
      tx_names.push(parts[0])
      let len_val = match parse_double(parts[1]) {
        Some(v) => v
        None => 0.0
      }
      length.push(len_val)
      let ab_val = match parse_double(parts[3]) {
        Some(v) => v
        None => 0.0
      }
      abundance.push(ab_val)
      let ct_val = match parse_double(parts[4]) {
        Some(v) => v
        None => 0.0
      }
      counts.push(ct_val)
    }
    i = i + 1
  }

  (tx_names, abundance, counts, length)
}

///|
pub fn tx_import(
  files : Array[String],
  tx2gene : Array[Tx2Gene],
  ignore_tx_version : Bool,
) -> TxImportData {
  let mut all_tx_names : Array[String] = Array::new()
  let all_abundance : Array[Array[Double]] = Array::new()
  let all_counts : Array[Array[Double]] = Array::new()
  let all_length : Array[Array[Double]] = Array::new()

  let sample_names : Array[String] = Array::new()

  let mut i = 0
  while i < files.length() {
    let file_content = files[i]
    let (tx_names, abundance, counts, length) = parse_quant_sf(file_content)

    if i == 0 {
      all_tx_names = tx_names.clone()
    }

    all_abundance.push(abundance)
    all_counts.push(counts)
    all_length.push(length)

    let name = extract_sample_name(file_content)
    if name.length() > 0 {
      sample_names.push(name)
    } else {
      sample_names.push("sample_" + i.to_string())
    }

    i = i + 1
  }

  let gene_names = map_tx_to_gene(all_tx_names, tx2gene, ignore_tx_version)

  tx_import_data_new(
    all_abundance, all_counts, all_length, all_tx_names, sample_names, gene_names,
  )
}

///|
fn extract_sample_name(content : String) -> String {
  let lines = split_by_char(content, '\n'.to_int())
  if lines.length() > 0 {
    let first_line = lines[0]
    if first_line.starts_with("#") {
      let idx = first_line.find(":")
      if idx.is_some() {
        return first_line[idx.unwrap() + 1:].to_string()
      }
    }
  }
  ""
}

///|
fn map_tx_to_gene(
  tx_names : Array[String],
  tx2gene : Array[Tx2Gene],
  ignore_tx_version : Bool,
) -> Array[String] {
  let gene_map : Map[String, String] = Map([], capacity=tx2gene.length())

  let mut i = 0
  while i < tx2gene.length() {
    let tx = tx2gene[i]
    let tx_id = if ignore_tx_version {
      remove_tx_version(tx.tx_id)
    } else {
      tx.tx_id
    }
    gene_map[tx_id] = tx.gene_name
    i = i + 1
  }

  let gene_names : Array[String] = Array::new()
  i = 0
  while i < tx_names.length() {
    let tx = if ignore_tx_version {
      remove_tx_version(tx_names[i])
    } else {
      tx_names[i]
    }
    if gene_map.contains(tx) {
      gene_names.push(gene_map[tx])
    } else {
      gene_names.push(tx_names[i])
    }
    i = i + 1
  }

  gene_names
}

///|
fn remove_tx_version(tx_id : String) -> String {
  let idx = tx_id.find(".")
  if idx.is_some() {
    tx_id[0:idx.unwrap()].to_string()
  } else {
    tx_id
  }
}

///|
pub fn summarize_to_gene(
  txi : TxImportData,
  tx2gene : Array[Tx2Gene],
  ignore_tx_version : Bool,
) -> TxImportData {
  let gene_to_tx : Map[String, Array[Int]] = Map([], capacity=100)

  let mut i = 0
  while i < txi.tx_names.length() {
    let tx = if ignore_tx_version {
      remove_tx_version(txi.tx_names[i])
    } else {
      txi.tx_names[i]
    }

    let mut gene_name = txi.gene_names[i]
    let mut j = 0
    while j < tx2gene.length() {
      let txg = tx2gene[j]
      let txg_id = if ignore_tx_version {
        remove_tx_version(txg.tx_id)
      } else {
        txg.tx_id
      }
      if txg_id == tx && txg.gene_name.length() > 0 {
        gene_name = txg.gene_name
        break
      }
      j = j + 1
    }

    if !gene_to_tx.contains(gene_name) {
      gene_to_tx[gene_name] = []
    }
    gene_to_tx[gene_name].push(i)
    i = i + 1
  }

  let gene_list = gene_to_tx.keys().collect()
  let ngene = gene_list.length()
  let nsample = txi.sample_names.length()

  let gene_abundance : Array[Array[Double]] = Array::new()
  let gene_counts : Array[Array[Double]] = Array::new()
  let gene_length : Array[Array[Double]] = Array::new()

  let mut k = 0
  while k < nsample {
    let row_ab : Array[Double] = Array::new()
    let row_ct : Array[Double] = Array::new()
    let row_len : Array[Double] = Array::new()
    let mut g = 0
    while g < ngene {
      row_ab.push(0.0)
      row_ct.push(0.0)
      row_len.push(0.0)
      g = g + 1
    }
    gene_abundance.push(row_ab)
    gene_counts.push(row_ct)
    gene_length.push(row_len)
    k = k + 1
  }

  let mut g = 0
  while g < ngene {
    let gene_name = gene_list[g]
    let tx_indices = gene_to_tx[gene_name]

    let mut k = 0
    while k < nsample {
      let mut total_abundance = 0.0
      let mut total_counts = 0.0
      let mut total_length = 0.0

      let mut j = 0
      while j < tx_indices.length() {
        let tx_idx = tx_indices[j]
        if tx_idx < txi.abundance[k].length() {
          total_abundance = total_abundance + txi.abundance[k][tx_idx]
        }
        if tx_idx < txi.counts[k].length() {
          total_counts = total_counts + txi.counts[k][tx_idx]
        }
        if tx_idx < txi.length[k].length() {
          total_length = total_length + txi.length[k][tx_idx]
        }
        j = j + 1
      }

      if tx_indices.length() > 0 {
        gene_abundance[k][g] = total_abundance
        gene_counts[k][g] = total_counts
        gene_length[k][g] = total_length / tx_indices.length().to_double()
      }
      k = k + 1
    }
    g = g + 1
  }

  tx_import_data_new(
    gene_abundance,
    gene_counts,
    gene_length,
    gene_list,
    txi.sample_names,
    gene_list,
  )
}

///|
pub fn tx_import_data_to_expression_set(
  txi : TxImportData,
  type_ : String,
) -> ExpressionSet {
  let assay_data = if type_ == "abundance" {
    transpose(txi.abundance)
  } else if type_ == "counts" {
    transpose(txi.counts)
  } else {
    transpose(txi.length)
  }

  let pheno_data : Map[String, Array[String]] = Map([], capacity=1)
  pheno_data["sample_name"] = txi.sample_names

  let feature_data : Map[String, Array[String]] = Map([], capacity=2)
  feature_data["feature_name"] = txi.gene_names
  feature_data["transcript_name"] = txi.tx_names

  ExpressionSet::new(assay_data, pheno_data, feature_data, "tximport", Map([]))
}

///|
pub fn tx_import_data_summary(txi : TxImportData) -> String {
  let ntx = txi.tx_names.length()
  let nsample = txi.sample_names.length()

  let mut result = "TxImportData Summary\n"
  result = result + "--------------------\n"
  result = result + "Transcripts: " + ntx.to_string() + "\n"
  result = result + "Samples: " + nsample.to_string() + "\n"

  if txi.sample_names.length() > 0 {
    result = result + "\nSample names:\n"
    let mut i = 0
    while i < txi.sample_names.length() && i < 5 {
      result = result + "  - " + txi.sample_names[i] + "\n"
      i = i + 1
    }
    if txi.sample_names.length() > 5 {
      result = result +
        "  ... and " +
        (txi.sample_names.length() - 5).to_string() +
        " more\n"
    }
  }

  if txi.gene_names.length() > 0 {
    let unique_genes = count_unique(txi.gene_names)
    result = result + "\nUnique genes: " + unique_genes.to_string() + "\n"
  }

  result
}

///|
fn count_unique(arr : Array[String]) -> Int {
  let seen : Map[String, Bool] = Map([], capacity=arr.length())
  let mut count = 0
  let mut i = 0
  while i < arr.length() {
    if !seen.contains(arr[i]) {
      seen[arr[i]] = true
      count = count + 1
    }
    i = i + 1
  }
  count
}

///|
pub fn tx_import_filter_low_expression(
  txi : TxImportData,
  min_counts : Double,
) -> TxImportData {
  let nsample = txi.counts.length()
  let ntx = if nsample > 0 { txi.counts[0].length() } else { 0 }

  let keep_tx : Array[Int] = Array::new()
  let mut i = 0
  while i < ntx {
    let mut total = 0.0
    let mut j = 0
    while j < nsample {
      if i < txi.counts[j].length() {
        total = total + txi.counts[j][i]
      }
      j = j + 1
    }
    if total >= min_counts {
      keep_tx.push(i)
    }
    i = i + 1
  }

  let new_abundance : Array[Array[Double]] = Array::new()
  let new_counts : Array[Array[Double]] = Array::new()
  let new_length : Array[Array[Double]] = Array::new()
  let new_tx_names : Array[String] = Array::new()
  let new_gene_names : Array[String] = Array::new()

  let mut k = 0
  while k < nsample {
    let row_ab : Array[Double] = Array::new()
    let row_ct : Array[Double] = Array::new()
    let row_len : Array[Double] = Array::new()

    let mut j = 0
    while j < keep_tx.length() {
      let tx_idx = keep_tx[j]
      if tx_idx < txi.abundance[k].length() {
        row_ab.push(txi.abundance[k][tx_idx])
      } else {
        row_ab.push(0.0)
      }
      if tx_idx < txi.counts[k].length() {
        row_ct.push(txi.counts[k][tx_idx])
      } else {
        row_ct.push(0.0)
      }
      if tx_idx < txi.length[k].length() {
        row_len.push(txi.length[k][tx_idx])
      } else {
        row_len.push(0.0)
      }
      j = j + 1
    }
    new_abundance.push(row_ab)
    new_counts.push(row_ct)
    new_length.push(row_len)
    k = k + 1
  }

  let mut j = 0
  while j < keep_tx.length() {
    let tx_idx = keep_tx[j]
    if tx_idx < txi.tx_names.length() {
      new_tx_names.push(txi.tx_names[tx_idx])
    }
    if tx_idx < txi.gene_names.length() {
      new_gene_names.push(txi.gene_names[tx_idx])
    }
    j = j + 1
  }

  tx_import_data_new(
    new_abundance,
    new_counts,
    new_length,
    new_tx_names,
    txi.sample_names,
    new_gene_names,
  )
}
