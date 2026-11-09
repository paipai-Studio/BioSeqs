/// Scanpy-inspired: read_excel — read Excel-like tabular data
/// Since MoonBit doesn't have native Excel parsing, this reads TSV/CSV
/// exported from Excel, or a simple JSON representation.

///|
/// Read an Excel-exported TSV file as AnnData.
/// Excel data is typically exported as tab-separated values.
///
/// Parameters:
/// - `content`: TSV file content (tab-separated)
/// - `sheet_name`: logical sheet name (for metadata, default "Sheet1")
/// - `first_col_is_index`: if true, first column is obs names (default true)
pub fn scanpy_read_excel_tsv(
  content : String,
  sheet_name~ : String = "Sheet1",
  _first_col_is_index~ : Bool = true,
) -> AnnData {
  let tab : String = "\t"
  let adata = scanpy_read_csv(content, delimiter=tab, first_col_is_index=true)
  ignore(adata.set_uns("excel_sheet", JsonStr(sheet_name)))
  adata
}

///|
/// Read multi-sheet Excel data from JSON representation.
/// Each sheet becomes a layer or obs/var annotation.
///
/// Parameters:
/// - `json_str`: JSON with sheets: {"sheet1": [[...]], "sheet2": [[...]]}
/// - `main_sheet`: name of the main expression sheet (default "Sheet1")
pub fn scanpy_read_excel_multi(
  json_str : String,
  main_sheet~ : String = "Sheet1",
) -> AnnData {
  let json = parse_json(json_str)
  match json {
    Some(j) => scanpy_parse_excel_json(j, main_sheet)
    None => AnnData::new([], [], [])
  }
}

///|
fn scanpy_parse_excel_json(j : JsonValue, main_sheet : String) -> AnnData {
  let x : Array[Array[Double]] = []
  let obs_names : Array[String] = []
  let var_names : Array[String] = []
  let layers : Map[String, Array[Array[Double]]] = Map([], capacity=0)
  match j {
    JsonObj(m) => {
      match m.get(main_sheet) {
        Some(JsonArr(rows)) =>
          for row in rows {
            match row {
              JsonArr(cols) => {
                let r : Array[Double] = []
                for v in cols {
                  match v {
                    JsonNum(n) => r.push(n)
                    _ => r.push(0.0)
                  }
                }
                x.push(r)
              }
              _ => ()
            }
          }
        _ => ()
      }
      for sheet_name in m.keys() {
        if sheet_name != main_sheet {
          match m.get(sheet_name) {
            Some(JsonArr(rows)) => {
              let mat : Array[Array[Double]] = []
              for row in rows {
                match row {
                  JsonArr(cols) => {
                    let r : Array[Double] = []
                    for v in cols {
                      match v {
                        JsonNum(n) => r.push(n)
                        _ => r.push(0.0)
                      }
                    }
                    mat.push(r)
                  }
                  _ => ()
                }
              }
              if mat.length() > 0 {
                layers[sheet_name] = mat
              }
            }
            _ => ()
          }
        }
      }
    }
    _ => ()
  }
  if x.length() > 0 {
    for i in 0..<x.length() {
      obs_names.push("cell_" + i.to_string())
    }
    for j in 0..<x[0].length() {
      var_names.push("gene_" + j.to_string())
    }
  }
  let adata = AnnData::new(x, obs_names, var_names)
  for layer_name in layers.keys() {
    match layers.get(layer_name) {
      Some(mat) => ignore(adata.set_layer(layer_name, mat))
      None => ()
    }
  }
  adata
}
