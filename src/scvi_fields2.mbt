///|
/// Data fields: NumericalJointObsField, NumericalJointVarField,
/// CategoricalJointObsField, CategoricalJointVarField,
/// ObsmField, VarmField, ProteinObsmField,
/// LabelsWithUnlabeledObsField, MuDataWrapper.
pub struct JointNumericalField {
  name : String
  attr_key : String
  columns : Array[String]
  is_obs : Bool
  data : Array[Array[Double]]
}

///|
pub fn joint_numerical_field_init(
  name : String,
  attr_key : String,
  columns : Array[String],
  is_obs : Bool,
) -> JointNumericalField {
  { name, attr_key, columns, is_obs, data: [], }
}

///|
pub fn joint_numerical_field_concat(
  field : JointNumericalField,
  obs_data : Map[String, Array[Double]],
) -> Array[Array[Double]] {
  let n = {
    let mut max_n = 0
    for col in field.columns {
      match obs_data.get(col) {
        Some(arr) => if arr.length() > max_n { max_n = arr.length() }
        None => ()
      }
    }
    max_n
  }
  let result = Array(capacity=n)
  for i in 0..<n {
    let row = Array(capacity=field.columns.length())
    for col in field.columns {
      match obs_data.get(col) {
        Some(arr) => row.push(if i < arr.length() { arr[i] } else { 0.0 })
        None => row.push(0.0)
      }
    }
    result.push(row)
  }
  result
}

///|
pub struct JointCategoricalField {
  name : String
  attr_key : String
  columns : Array[String]
  is_obs : Bool
  mappings : Map[String, Map[String, Int32]]
}

///|
pub fn joint_categorical_field_init(
  name : String,
  attr_key : String,
  columns : Array[String],
  is_obs : Bool,
) -> JointCategoricalField {
  { name, attr_key, columns, is_obs, mappings: Map([]), }
}

///|
pub fn joint_categorical_field_encode(
  field : JointCategoricalField,
  data : Map[String, Array[String]],
) -> Array[Array[Int32]] {
  let f = field
  let n = {
    let mut max_n = 0
    for col in f.columns {
      match data.get(col) {
        Some(arr) => if arr.length() > max_n { max_n = arr.length() }
        None => ()
      }
    }
    max_n
  }
  for col in f.columns {
    let col_map = Map([])
    let mut idx = 0
    match data.get(col) {
      Some(arr) =>
        for i in 0..<arr.length() {
          match col_map.get(arr[i]) {
            Some(_) => ()
            None => {
              col_map[arr[i]] = idx
              idx = idx + 1
            }
          }
        }
      None => ()
    }
    f.mappings[col] = col_map
  }
  let result = Array(capacity=n)
  for i in 0..<n {
    let row = Array(capacity=f.columns.length())
    for col in f.columns {
      match data.get(col) {
        Some(arr) =>
          if i < arr.length() {
            match f.mappings.get(col) {
              Some(col_map) =>
                match col_map.get(arr[i]) {
                  Some(v) => row.push(v)
                  None => row.push(-1)
                }
              None => row.push(-1)
            }
          } else {
            row.push(-1)
          }
        None => row.push(-1)
      }
    }
    result.push(row)
  }
  result
}

///|
pub struct ObsmFieldData {
  name : String
  key : String
  data : Array[Array[Double]]
  n_vars : Int32
}

///|
pub fn obsm_field_init(
  name : String,
  key : String,
  n_vars : Int32,
) -> ObsmFieldData {
  { name, key, data: [], n_vars, }
}

///|
pub fn obsm_field_set_data(
  field : ObsmFieldData,
  data : Array[Array[Double]],
) -> ObsmFieldData {
  { ..field, data, }
}

///|
pub struct VarmFieldData {
  name : String
  key : String
  data : Array[Array[Double]]
  n_obs : Int32
}

///|
pub fn varm_field_init(
  name : String,
  key : String,
  n_obs : Int32,
) -> VarmFieldData {
  { name, key, data: [], n_obs, }
}

///|
pub fn varm_field_set_data(
  field : VarmFieldData,
  data : Array[Array[Double]],
) -> VarmFieldData {
  { ..field, data, }
}

///|
pub struct ProteinObsmFieldData {
  name : String
  key : String
  protein_names : Array[String]
  data : Array[Array[Double]]
  background : Array[Array[Double]]
}

///|
pub fn protein_obsm_field_init(
  name : String,
  key : String,
  protein_names : Array[String],
) -> ProteinObsmFieldData {
  { name, key, protein_names, data: [], background: [], }
}

///|
pub fn protein_obsm_field_set_data(
  field : ProteinObsmFieldData,
  data : Array[Array[Double]],
) -> ProteinObsmFieldData {
  { ..field, data, }
}

///|
pub fn protein_obsm_field_set_background(
  field : ProteinObsmFieldData,
  background : Array[Array[Double]],
) -> ProteinObsmFieldData {
  { ..field, background, }
}

///|
pub fn protein_obsm_field_correct_background(
  field : ProteinObsmFieldData,
) -> Array[Array[Double]] {
  let result = Array(capacity=field.data.length())
  for i in 0..<field.data.length() {
    let row = Array(capacity=field.data[i].length())
    for j in 0..<field.data[i].length() {
      let bg = if i < field.background.length() &&
        j < field.background[i].length() {
        field.background[i][j]
      } else {
        0.0
      }
      let corrected = field.data[i][j] - bg
      row.push(if corrected < 0.0 { 0.0 } else { corrected })
    }
    result.push(row)
  }
  result
}

///|
pub struct LabelsWithUnlabeledField {
  name : String
  attr_key : String
  unlabeled_category : String
  labels : Array[String]
  is_labeled : Array[Bool]
}

///|
pub fn labels_with_unlabeled_init(
  name : String,
  attr_key : String,
  unlabeled_category : String,
) -> LabelsWithUnlabeledField {
  { name, attr_key, unlabeled_category, labels: [], is_labeled: [], }
}

///|
pub fn labels_with_unlabeled_set_labels(
  field : LabelsWithUnlabeledField,
  raw_labels : Array[String],
) -> LabelsWithUnlabeledField {
  let labels = Array(capacity=raw_labels.length())
  let is_labeled = Array(capacity=raw_labels.length())
  for i in 0..<raw_labels.length() {
    labels.push(raw_labels[i])
    is_labeled.push(raw_labels[i] != field.unlabeled_category)
  }
  { ..field, labels, is_labeled, }
}

///|
pub fn labels_with_unlabeled_get_labeled_mask(
  field : LabelsWithUnlabeledField,
) -> Array[Bool] {
  field.is_labeled
}

///|
pub fn labels_with_unlabeled_get_unique_labels(
  field : LabelsWithUnlabeledField,
) -> Array[String] {
  let seen = Set([])
  let result = Array(capacity=0)
  for i in 0..<field.labels.length() {
    if field.is_labeled[i] {
      if not(seen.contains(field.labels[i])) {
        seen.add(field.labels[i])
        result.push(field.labels[i])
      }
    }
  }
  result
}

///|
pub struct MuDataWrapper {
  modality_keys : Array[String]
  modality_data : Map[String, (Int32, Int32)]
  shared_obs : Array[String]
  shared_var : Array[String]
}

///|
pub fn mudata_wrapper_init(modality_keys : Array[String]) -> MuDataWrapper {
  { modality_keys, modality_data: Map([]), shared_obs: [], shared_var: [], }
}

///|
pub fn mudata_wrapper_register_modality(
  wrapper : MuDataWrapper,
  key : String,
  n_obs : Int32,
  n_var : Int32,
) -> MuDataWrapper {
  let w = wrapper
  w.modality_data[key] = (n_obs, n_var)
  w
}

///|
pub fn mudata_wrapper_get_shape(
  wrapper : MuDataWrapper,
  key : String,
) -> (Int32, Int32) {
  match wrapper.modality_data.get(key) {
    Some(shape) => shape
    None => (0, 0)
  }
}

///|
pub fn mudata_wrapper_get_total_obs(wrapper : MuDataWrapper) -> Int32 {
  let mut max_obs = 0
  for key in wrapper.modality_data.keys() {
    match wrapper.modality_data.get(key) {
      Some((n_obs, _)) => if n_obs > max_obs { max_obs = n_obs }
      None => ()
    }
  }
  max_obs
}

///|
pub fn _scvi_fields2_dummy() -> Unit {
  let _ = JointNumericalField::{
    name: "",
    attr_key: "",
    columns: [],
    is_obs: false,
    data: [],
  }
  let _ = JointCategoricalField::{
    name: "",
    attr_key: "",
    columns: [],
    is_obs: false,
    mappings: Map([]),
  }
  let _ = ObsmFieldData::{ name: "", key: "", data: [], n_vars: 0, }
  let _ = VarmFieldData::{ name: "", key: "", data: [], n_obs: 0, }
  let _ = ProteinObsmFieldData::{
    name: "",
    key: "",
    protein_names: [],
    data: [],
    background: [],
  }
  let _ = LabelsWithUnlabeledField::{
    name: "",
    attr_key: "",
    unlabeled_category: "",
    labels: [],
    is_labeled: [],
  }
  let _ = MuDataWrapper::{
    modality_keys: [],
    modality_data: Map([]),
    shared_obs: [],
    shared_var: [],
  }
}
