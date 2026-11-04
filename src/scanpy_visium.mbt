/// Scanpy-inspired: read_10x_visium — read 10x Visium spatial transcriptomics data

///|
/// Spatial transcriptomics data container.
pub struct VisiumData {
  adata : AnnData // expression matrix
  spatial_coords : Array[Array[Double]] // spot coordinates (n_spots × 2)
  image_coords : Array[Array[Double]] // image pixel coordinates (n_spots × 2)
  spot_diameter : Double // spot diameter in pixels
  scale_factor : Double // scale factor for image
}

///|
/// Read 10x Visium spatial data from JSON representation.
/// Expects: "matrix" (counts), "barcodes", "features",
/// "coords" (array_x, array_y), "image_coords" (pxl_row, pxl_col),
/// "diameter", "scalef".
///
/// Parameters:
/// - `json_str`: JSON string with Visium data
pub fn scanpy_read_10x_visium(json_str : String) -> VisiumData {
  let json = parse_json(json_str)
  match json {
    Some(j) => scanpy_parse_visium_json(j)
    None =>
      VisiumData::{
        adata: AnnData::new([], [], []),
        spatial_coords: [],
        image_coords: [],
        spot_diameter: 0.0,
        scale_factor: 1.0,
      }
  }
}

///|
fn scanpy_parse_visium_json(j : JsonValue) -> VisiumData {
  // Parse expression matrix
  let x = Array(capacity=0)
  let obs_names = Array(capacity=0)
  let var_names = Array(capacity=0)
  let spatial_coords = Array(capacity=0)
  let image_coords = Array(capacity=0)
  let mut diameter = 0.0
  let mut scalef = 1.0
  match j {
    JsonObj(m) => {
      // Parse matrix
      match m.get("matrix") {
        Some(JsonArr(rows)) =>
          for row in rows {
            match row {
              JsonArr(cols) => {
                let r = Array(capacity=0)
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
      // Parse barcodes
      match m.get("barcodes") {
        Some(JsonArr(arr)) =>
          for item in arr {
            match item {
              JsonStr(s) => obs_names.push(s)
              _ => obs_names.push("")
            }
          }
        _ => ()
      }
      // Parse features
      match m.get("features") {
        Some(JsonArr(arr)) =>
          for item in arr {
            match item {
              JsonStr(s) => var_names.push(s)
              _ => var_names.push("")
            }
          }
        _ => ()
      }
      // Parse spatial coordinates
      match m.get("coords") {
        Some(JsonArr(arr)) =>
          for item in arr {
            match item {
              JsonArr(c) => {
                let coord = Array(capacity=2)
                for v in c {
                  match v {
                    JsonNum(n) => coord.push(n)
                    _ => coord.push(0.0)
                  }
                }
                spatial_coords.push(coord)
              }
              _ => ()
            }
          }
        _ => ()
      }
      // Parse image coordinates
      match m.get("image_coords") {
        Some(JsonArr(arr)) =>
          for item in arr {
            match item {
              JsonArr(c) => {
                let coord = Array(capacity=2)
                for v in c {
                  match v {
                    JsonNum(n) => coord.push(n)
                    _ => coord.push(0.0)
                  }
                }
                image_coords.push(coord)
              }
              _ => ()
            }
          }
        _ => ()
      }
      // Parse diameter
      match m.get("diameter") {
        Some(JsonNum(n)) => diameter = n
        _ => ()
      }
      // Parse scale factor
      match m.get("scalef") {
        Some(JsonNum(n)) => scalef = n
        _ => ()
      }
    }
    _ => ()
  }
  // Fill names if empty
  if obs_names.length() == 0 {
    for i in 0..<x.length() {
      obs_names.push("spot_" + i.to_string())
    }
  }
  if var_names.length() == 0 && x.length() > 0 {
    for j in 0..<x[0].length() {
      var_names.push("gene_" + j.to_string())
    }
  }
  VisiumData::{
    adata: AnnData::new(x, obs_names, var_names),
    spatial_coords,
    image_coords,
    spot_diameter: diameter,
    scale_factor: scalef,
  }
}

///|
/// Filter Visium spots by tissue area.
/// Parameters:
/// - `visium`: Visium data
/// - `in_tissue`: boolean mask (true = in tissue)
pub fn scanpy_visium_filter_tissue(
  visium : VisiumData,
  in_tissue : Array[Bool],
) -> VisiumData {
  let n = visium.adata.n_obs()
  let keep = Array(capacity=0)
  for i in 0..<n {
    if i < in_tissue.length() && in_tissue[i] {
      keep.push(i)
    }
  }
  // Subset expression
  let x = Array(capacity=keep.length())
  let obs_names = Array(capacity=keep.length())
  let coords = Array(capacity=keep.length())
  let img_coords = Array(capacity=keep.length())
  for k in keep {
    x.push(visium.adata.x[k])
    obs_names.push(visium.adata.obs_names[k])
    if k < visium.spatial_coords.length() {
      coords.push(visium.spatial_coords[k])
    }
    if k < visium.image_coords.length() {
      img_coords.push(visium.image_coords[k])
    }
  }
  VisiumData::{
    adata: AnnData::new(x, obs_names, visium.adata.get_var_names()),
    spatial_coords: coords,
    image_coords: img_coords,
    spot_diameter: visium.spot_diameter,
    scale_factor: visium.scale_factor,
  }
}
