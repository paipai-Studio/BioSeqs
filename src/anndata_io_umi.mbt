/// anndata_io_umi.mbt — Read UMI-tools output format
/// Item 5: read_umi_tools

///|
/// UMI-tools output: tab-delimited with columns like:
/// gene_id, cell_1, cell_2, ... (or barcode, gene, count long format)
pub fn anndata_read_umi_tools(
  content : String,
  format~ : String = "wide",
) -> AnnData? {
  let lines = content.trim().split("\n").to_array()
  if lines.length() < 2 {
    return None
  }
  if format == "long" {
    // Long format: barcode\tgene\tcount
    let barcodes : Array[String] = Array(capacity=0)
    let genes : Array[String] = Array(capacity=0)
    let bc_set : Set[String] = Set([])
    let gene_set : Set[String] = Set([])
    let bc_idx : Map[String, Int32] = Map([], capacity=16)
    let gene_idx : Map[String, Int32] = Map([], capacity=16)
    let values : Array[(Int32, Int32, Double)] = Array(capacity=0)
    for r = 1; r < lines.length(); r = r + 1 {
      let fields = lines[r].split("\t").to_array()
      if fields.length() < 3 {
        continue
      }
      let bc = fields[0].to_owned()
      let gene = fields[1].to_owned()
      let count = parse_double(fields[2].to_owned()).unwrap_or(0.0)
      if not(bc_set.contains(bc)) {
        bc_set.add(bc)
        bc_idx[bc] = barcodes.length()
        barcodes.push(bc)
      }
      if not(gene_set.contains(gene)) {
        gene_set.add(gene)
        gene_idx[gene] = genes.length()
        genes.push(gene)
      }
      values.push((bc_idx[bc], gene_idx[gene], count))
    }
    // Build matrix
    let n_obs = barcodes.length()
    let n_vars = genes.length()
    let x : Array[Array[Double]] = Array(capacity=n_obs)
    for _i = 0; _i < n_obs; _i = _i + 1 {
      let row : Array[Double] = Array(capacity=n_vars)
      for _j = 0; _j < n_vars; _j = _j + 1 {
        row.push(0.0)
      }
      x.push(row)
    }
    for v in values {
      x[v.0][v.1] = v.2
    }
    return Some(AnnData::new(x, barcodes, genes))
  }
  // Wide format: gene_id, cell_1, cell_2, ...
  let header = lines[0].split("\t").to_array()
  if header.length() < 2 {
    return None
  }
  let obs_names : Array[String] = Array(capacity=header.length() - 1)
  for c = 1; c < header.length(); c = c + 1 {
    obs_names.push(header[c].to_owned())
  }
  let var_names : Array[String] = Array(capacity=lines.length() - 1)
  let x : Array[Array[Double]] = Array(capacity=lines.length() - 1)
  // Transpose: input is var × obs, AnnData is obs × var
  let n_obs = obs_names.length()
  let n_vars = lines.length() - 1
  for _i = 0; _i < n_obs; _i = _i + 1 {
    let row : Array[Double] = Array(capacity=n_vars)
    for _j = 0; _j < n_vars; _j = _j + 1 {
      row.push(0.0)
    }
    x.push(row)
  }
  for r = 1; r < lines.length(); r = r + 1 {
    let fields = lines[r].split("\t").to_array()
    if fields.length() < 1 {
      continue
    }
    var_names.push(fields[0].to_owned())
    for c = 1; c < fields.length() && c - 1 < n_obs; c = c + 1 {
      x[c - 1][r - 1] = parse_double(fields[c].to_owned()).unwrap_or(0.0)
    }
  }
  Some(AnnData::new(x, obs_names, var_names))
}
