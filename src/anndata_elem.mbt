/// anndata_elem.mbt — read_elem / write_elem + encoding type system
/// Items 15-17: Low-level element I/O and encoding registry

///|
/// Encoding type identifiers
pub enum AnnDataEncoding {
  EncArray // dense array
  CSR // compressed sparse row
  CSC // compressed sparse column
  Scalar // single value
  Dict // dictionary
  StringArray // array of strings
  NullableStringArray // nullable strings
  DataFrame // table-like
}

///|
/// Encoding spec: type + version
pub struct AnnDataIOSpec {
  encoding_type : AnnDataEncoding
  encoding_version : String
}

///|
/// Get encoding type name as string
pub fn anndata_encoding_name(enc : AnnDataEncoding) -> String {
  match enc {
    EncArray => "array"
    CSR => "csr_matrix"
    CSC => "csc_matrix"
    Scalar => "scalar"
    Dict => "dict"
    StringArray => "string-array"
    NullableStringArray => "nullable-string-array"
    DataFrame => "dataframe"
  }
}

///|
/// Detect encoding type from a value.
pub fn anndata_detect_encoding(value : JsonValue) -> AnnDataIOSpec {
  match value {
    JsonNum(_) =>
      AnnDataIOSpec::{ encoding_type: Scalar, encoding_version: "0.1.0", }
    JsonStr(_) =>
      AnnDataIOSpec::{ encoding_type: Scalar, encoding_version: "0.1.0", }
    JsonArr(_) =>
      AnnDataIOSpec::{ encoding_type: EncArray, encoding_version: "0.2.0", }
    JsonObj(_) =>
      AnnDataIOSpec::{ encoding_type: Dict, encoding_version: "0.1.0", }
    JsonNull =>
      AnnDataIOSpec::{ encoding_type: Scalar, encoding_version: "0.1.0", }
    JsonBool(_) =>
      AnnDataIOSpec::{ encoding_type: Scalar, encoding_version: "0.1.0", }
  }
}

///|
/// read_elem: read a single element from a JSON-like storage group.
pub fn anndata_read_elem(value : JsonValue) -> JsonValue {
  // In a full implementation, this would decode based on encoding metadata.
  // Here we just return the value as-is.
  value
}

///|
/// write_elem: write a single element to a JSON-like storage with encoding.
pub fn anndata_write_elem(
  key : String,
  value : JsonValue,
) -> (String, JsonValue, AnnDataIOSpec) {
  let spec = anndata_detect_encoding(value)
  (key, value, spec)
}

///|
/// read_elem from a 2D array (dense matrix).
pub fn anndata_read_elem_array(value : JsonValue) -> Array[Array[Double]]? {
  match value {
    JsonArr(arr) => {
      let result : Array[Array[Double]] = Array(capacity=arr.length())
      for row in arr {
        match row {
          JsonArr(inner) => {
            let r : Array[Double] = Array(capacity=inner.length())
            for v in inner {
              match v {
                JsonNum(d) => r.push(d)
                _ => r.push(0.0)
              }
            }
            result.push(r)
          }
          _ => result.push(Array(capacity=0))
        }
      }
      Some(result)
    }
    _ => None
  }
}

///|
/// write_elem for a 2D array.
pub fn anndata_write_elem_array(
  key : String,
  matrix : Array[Array[Double]],
) -> (String, JsonValue) {
  let rows : Array[JsonValue] = Array(capacity=matrix.length())
  for row in matrix {
    let cols : Array[JsonValue] = Array(capacity=row.length())
    for v in row {
      cols.push(JsonNum(v))
    }
    rows.push(JsonArr(cols))
  }
  (key, JsonArr(rows))
}

///|
/// read_elem for a scalar.
pub fn anndata_read_elem_scalar(value : JsonValue) -> Double? {
  match value {
    JsonNum(d) => Some(d)
    _ => None
  }
}

///|
/// write_elem for a scalar.
pub fn anndata_write_elem_scalar(
  key : String,
  value : Double,
) -> (String, JsonValue) {
  (key, JsonNum(value))
}

///|
/// read_elem for a string array.
pub fn anndata_read_elem_string_array(value : JsonValue) -> Array[String]? {
  match value {
    JsonArr(arr) => {
      let result : Array[String] = Array(capacity=arr.length())
      for v in arr {
        match v {
          JsonStr(s) => result.push(s)
          _ => result.push("")
        }
      }
      Some(result)
    }
    _ => None
  }
}

///|
/// write_elem for a string array.
pub fn anndata_write_elem_string_array(
  key : String,
  values : Array[String],
) -> (String, JsonValue) {
  let arr : Array[JsonValue] = Array(capacity=values.length())
  for v in values {
    arr.push(JsonStr(v))
  }
  (key, JsonArr(arr))
}

///|
/// Encoding registry: map type name to encoding.
pub fn anndata_encoding_registry(type_name : String) -> AnnDataEncoding? {
  if type_name == "array" {
    Some(EncArray)
  } else if type_name == "csr_matrix" {
    Some(CSR)
  } else if type_name == "csc_matrix" {
    Some(CSC)
  } else if type_name == "scalar" {
    Some(Scalar)
  } else if type_name == "dict" {
    Some(Dict)
  } else if type_name == "string-array" {
    Some(StringArray)
  } else if type_name == "nullable-string-array" {
    Some(NullableStringArray)
  } else if type_name == "dataframe" {
    Some(DataFrame)
  } else {
    None
  }
}

///|
pub fn anndata_elem_dummy() -> Unit {
  ignore(AnnDataIOSpec::{ encoding_type: EncArray, encoding_version: "0.1.0", })
  ignore(EncArray)
  ignore(CSR)
  ignore(CSC)
  ignore(Scalar)
  ignore(Dict)
  ignore(StringArray)
  ignore(NullableStringArray)
  ignore(DataFrame)
}
