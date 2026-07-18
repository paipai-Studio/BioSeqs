///|
pub struct GEOData {
  series_id : String
  sample_ids : Array[String]
  platform_id : String
  data : Array[Array[Double]]
  gene_names : Array[String]
  sample_info : Map[String, Array[String]]
}

///|
pub struct GEOSeries {
  gse_id : String
  title : String
  description : String
  platform_ids : Array[String]
  sample_ids : Array[String]
}

///|
pub struct GEOSample {
  gsm_id : String
  title : String
  description : String
  platform_id : String
  characteristics : Map[String, String]
}

///|
pub struct GEOPlatform {
  gpl_id : String
  title : String
  organism : String
  probe_count : Int
}

///|
pub fn geo_data_new(
  series_id : String,
  sample_ids : Array[String],
  platform_id : String,
  data : Array[Array[Double]],
  gene_names : Array[String],
  sample_info : Map[String, Array[String]],
) -> GEOData {
  GEOData::{ series_id, sample_ids, platform_id, data, gene_names, sample_info }
}

///|
pub fn geo_series_new(
  gse_id : String,
  title : String,
  description : String,
  platform_ids : Array[String],
  sample_ids : Array[String],
) -> GEOSeries {
  GEOSeries::{ gse_id, title, description, platform_ids, sample_ids }
}

///|
pub fn geo_sample_new(
  gsm_id : String,
  title : String,
  description : String,
  platform_id : String,
  characteristics : Map[String, String],
) -> GEOSample {
  GEOSample::{ gsm_id, title, description, platform_id, characteristics }
}

///|
pub fn geo_platform_new(
  gpl_id : String,
  title : String,
  organism : String,
  probe_count : Int,
) -> GEOPlatform {
  GEOPlatform::{ gpl_id, title, organism, probe_count }
}

///|
fn parse_geo_header_value(line : String) -> String {
  let idx = line.find("=")
  if idx.is_none() {
    return ""
  }
  let pos = idx.unwrap() + 1
  let result = line[pos:].to_string()
  let mut start = 0
  while start < result.length() &&
        result.unsafe_get(start).to_int() == ' '.to_int() {
    start = start + 1
  }
  result[start:].to_string()
}

///|
fn parse_geo_header_key_value(line : String) -> (String, String)? {
  let idx = line.find("=")
  if idx.is_none() {
    return None
  }
  let pos = idx.unwrap()
  let key = line[1:pos].to_string()
  let value = line[pos + 1:].to_string()
  let mut start = 0
  while start < value.length() &&
        value.unsafe_get(start).to_int() == ' '.to_int() {
    start = start + 1
  }
  Some((key, value[start:].to_string()))
}

///|
pub fn parse_geo_series_matrix(content : String) -> GEOData {
  let lines = split_by_char(content, '\n'.to_int())

  let mut series_id = ""
  let mut platform_id = ""
  let sample_ids : Array[String] = Array::new()
  let gene_names : Array[String] = Array::new()
  let data : Array[Array[Double]] = Array::new()
  let sample_info : Map[String, Array[String]] = Map([], capacity=10)

  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    if line.starts_with("!Series_geo_accession") {
      series_id = parse_geo_header_value(line)
    } else if line.starts_with("!Platform_geo_accession") {
      platform_id = parse_geo_header_value(line)
    } else if line.starts_with("!Sample_geo_accession") {
      sample_ids.push(parse_geo_header_value(line))
    } else if line.starts_with("!") &&
      !line.starts_with("!series_matrix_table_begin") {
      let key_value = parse_geo_header_key_value(line)
      if key_value.is_some() {
        let (key, value) = key_value.unwrap()
        if !sample_info.contains(key) {
          sample_info[key] = []
        }
        sample_info[key].push(value)
      }
    } else if line.starts_with("ID_REF") {
      sample_ids.clear()
      let header_parts = split_by_char(line, 9)
      if header_parts.length() > 1 {
        let mut j = 1
        while j < header_parts.length() {
          let id = header_parts[j]
          if id.length() > 0 {
            sample_ids.push(id)
          }
          j = j + 1
        }
      }
    } else if line.length() > 0 && !line.starts_with("!") {
      let parts = split_by_char(line, 9)
      if parts.length() > 1 {
        let gene_name = parts[0]
        gene_names.push(gene_name)

        let row : Array[Double] = []
        let mut j = 1
        while j < parts.length() {
          let val = parts[j]
          if val.length() > 0 {
            let parsed = match parse_double(val) {
              Some(v) => v
              None => 0.0
            }
            row.push(parsed)
          } else {
            row.push(0.0)
          }
          j = j + 1
        }
        data.push(row)
      }
    }
    i = i + 1
  }

  geo_data_new(
    series_id, sample_ids, platform_id, data, gene_names, sample_info,
  )
}

///|
pub fn parse_soft_format(content : String) -> GEOSeries {
  let lines = split_by_char(content, '\n'.to_int())

  let mut gse_id = ""
  let mut title = ""
  let mut description = ""
  let platform_ids : Array[String] = Array::new()
  let sample_ids : Array[String] = Array::new()

  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    if line.starts_with("^SAMPLE=") {
      let idx = line.find("=")
      if idx.is_some() {
        sample_ids.push(line[idx.unwrap() + 1:].to_string())
      }
    } else if line.starts_with("!Sample_platform_id") {
      let val = parse_geo_header_value(line)
      if val.length() > 0 && !platform_ids.contains(val) {
        platform_ids.push(val)
      }
    } else if line.starts_with("!Series_geo_accession") {
      gse_id = parse_geo_header_value(line)
    } else if line.starts_with("!Series_title") {
      title = parse_geo_header_value(line)
    } else if line.starts_with("!Series_summary") {
      description = parse_geo_header_value(line)
    }
    i = i + 1
  }

  geo_series_new(gse_id, title, description, platform_ids, sample_ids)
}

///|
pub fn geo_data_to_expression_set(geo : GEOData) -> ExpressionSet {
  let pheno_data : Map[String, Array[String]] = Map(
    [],
    capacity=geo.sample_info.length(),
  )
  let pheno_keys = geo.sample_info.keys().collect()
  let mut i = 0
  while i < pheno_keys.length() {
    let key = pheno_keys[i]
    pheno_data[key] = geo.sample_info[key]
    i = i + 1
  }
  pheno_data["sample_name"] = geo.sample_ids

  let feature_data : Map[String, Array[String]] = Map([], capacity=1)
  feature_data["feature_name"] = geo.gene_names

  ExpressionSet::new(
    geo.data,
    pheno_data,
    feature_data,
    geo.platform_id,
    Map([]),
  )
}

///|
pub fn geo_data_filter_genes(
  geo : GEOData,
  min_expression : Double,
  max_expression : Double,
) -> GEOData {
  let filtered_data : Array[Array[Double]] = Array::new()
  let filtered_genes : Array[String] = Array::new()

  let mut i = 0
  while i < geo.data.length() {
    let row = geo.data[i]
    let mut keep = true
    let mut j = 0
    while j < row.length() {
      if row[j] < min_expression || row[j] > max_expression {
        keep = false
      }
      j = j + 1
    }
    if keep {
      filtered_data.push(row)
      if i < geo.gene_names.length() {
        filtered_genes.push(geo.gene_names[i])
      }
    }
    i = i + 1
  }

  geo_data_new(
    geo.series_id,
    geo.sample_ids,
    geo.platform_id,
    filtered_data,
    filtered_genes,
    geo.sample_info,
  )
}

///|
pub fn geo_data_log2_transform(geo : GEOData) -> GEOData {
  let transformed_data : Array[Array[Double]] = Array::new()

  let mut i = 0
  while i < geo.data.length() {
    let row = geo.data[i]
    let new_row : Array[Double] = Array::new()
    let mut j = 0
    while j < row.length() {
      if row[j] > 0.0 {
        new_row.push(@math.log2(row[j]))
      } else {
        new_row.push(0.0)
      }
      j = j + 1
    }
    transformed_data.push(new_row)
    i = i + 1
  }

  geo_data_new(
    geo.series_id,
    geo.sample_ids,
    geo.platform_id,
    transformed_data,
    geo.gene_names,
    geo.sample_info,
  )
}

///|
pub fn geo_data_summary(geo : GEOData) -> String {
  let ngenes = geo.gene_names.length()
  let nsamples = geo.sample_ids.length()
  let nrows = geo.data.length()

  let mut result = "GEOData Summary\n"
  result = result + "---------------\n"
  result = result + "Series ID: " + geo.series_id + "\n"
  result = result + "Platform: " + geo.platform_id + "\n"
  result = result + "Genes: " + ngenes.to_string() + "\n"
  result = result + "Samples: " + nsamples.to_string() + "\n"
  result = result + "Data rows: " + nrows.to_string() + "\n"

  if geo.sample_ids.length() > 0 {
    result = result + "\nSample IDs:\n"
    let mut i = 0
    while i < geo.sample_ids.length() && i < 5 {
      result = result + "  - " + geo.sample_ids[i] + "\n"
      i = i + 1
    }
    if geo.sample_ids.length() > 5 {
      result = result +
        "  ... and " +
        (geo.sample_ids.length() - 5).to_string() +
        " more\n"
    }
  }

  result
}

///|
pub fn geo_series_summary(series : GEOSeries) -> String {
  let mut result = "GEOSeries Summary\n"
  result = result + "-----------------\n"
  result = result + "GSE ID: " + series.gse_id + "\n"
  result = result + "Title: " + series.title + "\n"
  result = result + "Description: " + series.description + "\n"
  result = result +
    "Platforms: " +
    series.platform_ids.length().to_string() +
    "\n"
  result = result + "Samples: " + series.sample_ids.length().to_string() + "\n"

  if series.platform_ids.length() > 0 {
    result = result + "\nPlatform IDs:\n"
    let mut i = 0
    while i < series.platform_ids.length() && i < 3 {
      result = result + "  - " + series.platform_ids[i] + "\n"
      i = i + 1
    }
  }

  if series.sample_ids.length() > 0 {
    result = result + "\nSample IDs:\n"
    let mut i = 0
    while i < series.sample_ids.length() && i < 5 {
      result = result + "  - " + series.sample_ids[i] + "\n"
      i = i + 1
    }
    if series.sample_ids.length() > 5 {
      result = result +
        "  ... and " +
        (series.sample_ids.length() - 5).to_string() +
        " more\n"
    }
  }

  result
}
