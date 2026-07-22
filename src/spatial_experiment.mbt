///|
/// SpatialExperiment - Spatial transcriptomics data structure.

///|
pub struct SpatialCoord {
  x : Double
  y : Double
  z : Double
}

///|
pub fn SpatialCoord::new(x : Double, y : Double, z : Double) -> SpatialCoord {
  SpatialCoord::{ x, y, z }
}

///|
pub fn SpatialCoord::new_2d(x : Double, y : Double) -> SpatialCoord {
  SpatialCoord::{ x, y, z: 0.0 }
}

///|
pub struct SpatialImage {
  id : String
  data : Array[Array[Double]]
  scale_factor : Double
}

///|
pub fn SpatialImage::new(id : String, data : Array[Array[Double]], scale_factor : Double) -> SpatialImage {
  SpatialImage::{ id, data, scale_factor }
}

///|
pub struct SpatialExperiment {
  assay : Map[String, Array[Array[Double]]]
  row_data : Array[Map[String, String]]
  col_data : Array[Map[String, String]]
  spatial_coords : Array[SpatialCoord]
  images : Array[SpatialImage]
  metadata : Map[String, String]
}

///|
pub fn SpatialExperiment::new() -> SpatialExperiment {
  SpatialExperiment::{
    assay: Map([], capacity=0),
    row_data: Array::new(),
    col_data: Array::new(),
    spatial_coords: Array::new(),
    images: Array::new(),
    metadata: Map([], capacity=0)
  }
}

///|
pub fn se_add_assay(se : SpatialExperiment, name : String, data : Array[Array[Double]]) -> SpatialExperiment {
  se.assay.set(name, data)
  se
}

///|
pub fn se_add_row(se : SpatialExperiment, row : Map[String, String]) -> SpatialExperiment {
  se.row_data.push(row)
  se
}

///|
pub fn se_add_col(se : SpatialExperiment, col : Map[String, String]) -> SpatialExperiment {
  se.col_data.push(col)
  se
}

///|
pub fn se_add_spatial_coord(se : SpatialExperiment, coord : SpatialCoord) -> SpatialExperiment {
  se.spatial_coords.push(coord)
  se
}

///|
pub fn se_add_image(se : SpatialExperiment, image : SpatialImage) -> SpatialExperiment {
  se.images.push(image)
  se
}

///|
pub fn se_set_metadata(se : SpatialExperiment, key : String, value : String) -> SpatialExperiment {
  se.metadata.set(key, value)
  se
}

///|
pub fn se_num_rows(se : SpatialExperiment) -> Int {
  se.row_data.length()
}

///|
pub fn se_num_cols(se : SpatialExperiment) -> Int {
  se.col_data.length()
}

///|
pub fn se_num_assays(se : SpatialExperiment) -> Int {
  se.assay.keys().length()
}

///|
pub fn se_get_spatial_range(se : SpatialExperiment) -> (Double, Double, Double, Double) {
  if se.spatial_coords.length() == 0 {
    return (0.0, 0.0, 0.0, 0.0)
  }
  
  let mut min_x = se.spatial_coords[0].x
  let mut max_x = se.spatial_coords[0].x
  let mut min_y = se.spatial_coords[0].y
  let mut max_y = se.spatial_coords[0].y
  
  let mut i = 1
  while i < se.spatial_coords.length() {
    let coord = se.spatial_coords[i]
    if coord.x < min_x {
      min_x = coord.x
    }
    if coord.x > max_x {
      max_x = coord.x
    }
    if coord.y < min_y {
      min_y = coord.y
    }
    if coord.y > max_y {
      max_y = coord.y
    }
    i = i + 1
  }
  
  (min_x, max_x, min_y, max_y)
}

///|
pub fn se_filter_spots_by_range(se : SpatialExperiment, min_x : Double, max_x : Double, min_y : Double, max_y : Double) -> SpatialExperiment {
  let filtered = SpatialExperiment::new()
  
  let mut i = 0
  while i < se.row_data.length() {
    filtered.row_data.push(se.row_data[i])
    i = i + 1
  }
  
  let mut j = 0
  while j < se.images.length() {
    filtered.images.push(se.images[j])
    j = j + 1
  }
  
  let mut k = 0
  while k < se.spatial_coords.length() {
    let coord = se.spatial_coords[k]
    if coord.x >= min_x && coord.x <= max_x && coord.y >= min_y && coord.y <= max_y {
      filtered.spatial_coords.push(coord)
      if k < se.col_data.length() {
        filtered.col_data.push(se.col_data[k])
      }
    }
    k = k + 1
  }
  
  filtered
}

///|
pub fn create_example_spatial_experiment() -> SpatialExperiment {
  let se = SpatialExperiment::new()
  let _ = se_set_metadata(se, "platform", "Visium")
  let _ = se_set_metadata(se, "sample_id", "sample001")
  
  se.row_data.push(Map([("gene_id", "ENSG000001"), ("gene_name", "ACTB")], capacity=2))
  se.row_data.push(Map([("gene_id", "ENSG000002"), ("gene_name", "GAPDH")], capacity=2))
  se.row_data.push(Map([("gene_id", "ENSG000003"), ("gene_name", "TP53")], capacity=2))
  
  se.col_data.push(Map([("spot_id", "AAACAAGTATCTCCCA"), ("tissue", "tumor")], capacity=2))
  se.col_data.push(Map([("spot_id", "AAACAGCTTTCAGAAG"), ("tissue", "tumor")], capacity=2))
  se.col_data.push(Map([("spot_id", "AAACAGGGTCTATATT"), ("tissue", "normal")], capacity=2))
  se.col_data.push(Map([("spot_id", "AAACATTTCCCGGATT"), ("tissue", "normal")], capacity=2))
  se.col_data.push(Map([("spot_id", "AAACCGGGTGTATGGC"), ("tissue", "tumor")], capacity=2))
  se.col_data.push(Map([("spot_id", "AAACCGTTCGTCCACA"), ("tissue", "normal")], capacity=2))
  
  se.spatial_coords.push(SpatialCoord::new_2d(100.0, 200.0))
  se.spatial_coords.push(SpatialCoord::new_2d(150.0, 250.0))
  se.spatial_coords.push(SpatialCoord::new_2d(200.0, 150.0))
  se.spatial_coords.push(SpatialCoord::new_2d(250.0, 200.0))
  se.spatial_coords.push(SpatialCoord::new_2d(300.0, 250.0))
  se.spatial_coords.push(SpatialCoord::new_2d(350.0, 150.0))
  
  let image_data : Array[Array[Double]] = Array::new()
  let image = SpatialImage::new("HE", image_data, 0.1)
  se.images.push(image)
  
  se
}