///|
/// rhdf5 - HDF5 file format support for BioSeqs.
/// Inspired by Bioconductor rhdf5 package.
/// Provides HDF5 file I/O operations for genomic data.

///|
pub struct HDF5Attribute {
  name : String
  data_type : String
  data : String
}

///|
pub fn HDF5Attribute::new(name : String, data_type : String, data : String) -> HDF5Attribute {
  HDF5Attribute::{ name, data_type, data }
}

///|
pub struct HDF5Dataset {
  name : String
  data_type : String
  dimensions : Array[Int]
  data : String
  attributes : @hashmap.HashMap[String, HDF5Attribute]
}

///|
pub fn HDF5Dataset::new(name : String, data_type : String, dimensions : Array[Int]) -> HDF5Dataset {
  HDF5Dataset::{
    name,
    data_type,
    dimensions,
    data: "",
    attributes: @hashmap.HashMap([])
  }
}

///|
pub fn HDF5Dataset::set_data(self : HDF5Dataset, data : String) -> HDF5Dataset {
  HDF5Dataset::{
    name: self.name,
    data_type: self.data_type,
    dimensions: self.dimensions,
    data: data,
    attributes: self.attributes
  }
}

///|
pub fn HDF5Dataset::add_attribute(self : HDF5Dataset, attr : HDF5Attribute) -> HDF5Dataset {
  let new_attrs = self.attributes
  new_attrs.set(attr.name, attr)
  HDF5Dataset::{
    name: self.name,
    data_type: self.data_type,
    dimensions: self.dimensions,
    data: self.data,
    attributes: new_attrs
  }
}

///|
pub fn HDF5Dataset::get_attribute(self : HDF5Dataset, name : String) -> HDF5Attribute? {
  if self.attributes.contains(name) {
    Some(self.attributes.get(name).unwrap())
  } else {
    None
  }
}

///|
pub struct HDF5Group {
  name : String
  groups : @hashmap.HashMap[String, HDF5Group]
  datasets : @hashmap.HashMap[String, HDF5Dataset]
  attributes : @hashmap.HashMap[String, HDF5Attribute]
}

///|
pub fn HDF5Group::new(name : String) -> HDF5Group {
  HDF5Group::{
    name,
    groups: @hashmap.HashMap([]),
    datasets: @hashmap.HashMap([]),
    attributes: @hashmap.HashMap([])
  }
}

///|
pub fn HDF5Group::find_group(self : HDF5Group, name : String) -> HDF5Group? {
  if self.groups.contains(name) {
    Some(self.groups.get(name).unwrap())
  } else {
    None
  }
}

///|
pub fn HDF5Group::find_dataset(self : HDF5Group, name : String) -> HDF5Dataset? {
  if self.datasets.contains(name) {
    Some(self.datasets.get(name).unwrap())
  } else {
    None
  }
}

///|
pub struct HDF5File {
  filename : String
  mode : String
  root : HDF5Group
}

///|
pub fn HDF5File::new(filename : String, mode : String) -> HDF5File {
  HDF5File::{
    filename,
    mode,
    root: HDF5Group::new("/")
  }
}

///|
fn split_path(path : String) -> Array[String] {
  let parts : Array[String] = Array::new()
  let mut current = ""
  let mut i = 0
  while i < path.length() {
    let c = path[i]
    if c == '/' {
      if current != "" {
        parts.push(current)
        current = ""
      }
    } else {
      current = current + c.to_string()
    }
    i = i + 1
  }
  if current != "" {
    parts.push(current)
  }
  parts
}

///|
fn split_path_and_name(path : String) -> (Array[String], String) {
  let parts = split_path(path)
  if parts.length() == 0 {
    (Array::new(), "")
  } else if parts.length() == 1 {
    (Array::new(), parts[0])
  } else {
    let name = parts[parts.length() - 1]
    let group_parts : Array[String] = Array::new()
    let mut i = 0
    while i < parts.length() - 1 {
      group_parts.push(parts[i])
      i = i + 1
    }
    (group_parts, name)
  }
}

///|
pub fn HDF5File::create_group(self : HDF5File, path : String) -> HDF5File {
  let parts = split_path(path)
  let _ = ensure_path(self.root, parts)
  self
}

///|
fn ensure_path(group : HDF5Group, parts : Array[String]) -> Unit {
  let mut current = group
  let mut i = 0
  while i < parts.length() {
    let part = parts[i]
    let g = current.groups
    if g.contains(part) {
      current = g.get(part).unwrap()
    } else {
      let new_group = HDF5Group::new(part)
      g.set(part, new_group)
      current = new_group
    }
    i = i + 1
  }
}

///|
pub fn HDF5File::write_dataset(self : HDF5File, path : String, dataset : HDF5Dataset) -> HDF5File {
  let (group_parts, ds_name) = split_path_and_name(path)
  
  let _ = ensure_path(self.root, group_parts)
  
  let mut current = self.root
  let mut i = 0
  while i < group_parts.length() {
    current = current.groups.get(group_parts[i]).unwrap()
    i = i + 1
  }
  
  let new_ds = dataset.set_name(ds_name)
  current.datasets.set(ds_name, new_ds)
  
  self
}

///|
pub fn HDF5File::read_dataset(self : HDF5File, path : String) -> HDF5Dataset? {
  let (group_parts, ds_name) = split_path_and_name(path)
  
  let mut current = self.root
  let mut i = 0
  while i < group_parts.length() {
    if current.groups.contains(group_parts[i]) {
      current = current.groups.get(group_parts[i]).unwrap()
    } else {
      return None
    }
    i = i + 1
  }
  
  if current.datasets.contains(ds_name) {
    Some(current.datasets.get(ds_name).unwrap())
  } else {
    None
  }
}

///|
pub fn HDF5Dataset::set_name(self : HDF5Dataset, name : String) -> HDF5Dataset {
  HDF5Dataset::{
    name: name,
    data_type: self.data_type,
    dimensions: self.dimensions,
    data: self.data,
    attributes: self.attributes
  }
}

///|
pub fn h5create_file(filename : String) -> HDF5File {
  HDF5File::new(filename, "w")
}

///|
pub fn h5create_group(file : HDF5File, path : String) -> HDF5File {
  file.create_group(path)
}

///|
pub fn h5create_dataset(file : HDF5File, path : String, data_type : String, dimensions : Array[Int]) -> HDF5File {
  let ds = HDF5Dataset::new("", data_type, dimensions)
  file.write_dataset(path, ds)
}

///|
pub fn h5write_dataset(file : HDF5File, path : String, data : String, data_type : String, dimensions : Array[Int]) -> HDF5File {
  let ds = HDF5Dataset::new("", data_type, dimensions)
  let ds_with_data = ds.set_data(data)
  file.write_dataset(path, ds_with_data)
}

///|
pub fn h5read_dataset(file : HDF5File, path : String) -> HDF5Dataset? {
  file.read_dataset(path)
}

///|
pub fn h5ls(file : HDF5File) -> Array[String] {
  let result : Array[String] = Array::new()
  ignore(list_group(file.root, "", result))
  result
}

///|
fn list_group(group : HDF5Group, prefix : String, result : Array[String]) -> Unit {
  let ds_entries = group.datasets.iter().collect()
  
  let mut i = 0
  while i < ds_entries.length() {
    let (name, ds) = ds_entries[i]
    let dims_str = dimensions_to_string(ds.dimensions)
    result.push(prefix + "/" + name + " [" + dims_str + "]")
    i = i + 1
  }
  
  let grp_entries = group.groups.iter().collect()
  
  let mut j = 0
  while j < grp_entries.length() {
    let (name, grp) = grp_entries[j]
    let new_prefix = if prefix == "" { name } else { prefix + "/" + name }
    result.push(new_prefix + "/")
    ignore(list_group(grp, new_prefix, result))
    j = j + 1
  }
}

///|
fn dimensions_to_string(dims : Array[Int]) -> String {
  let mut result = ""
  let mut i = 0
  while i < dims.length() {
    if i > 0 {
      result = result + " x "
    }
    result = result + dims[i].to_string()
    i = i + 1
  }
  if result == "" {
    result = "scalar"
  }
  result
}

///|
pub fn h5create_example_file() -> HDF5File {
  let file = h5create_file("example.h5")
  let file2 = h5create_group(file, "data")
  let file3 = h5create_dataset(file2, "/data/matrix", "H5T_NATIVE_DOUBLE", [100, 1000])
  let file4 = h5write_dataset(file3, "/data/matrix", "1.0,2.0,3.0", "H5T_NATIVE_DOUBLE", [3])
  file4
}

///|
pub fn h5create_attribute(name : String, data_type : String, data : String) -> HDF5Attribute {
  HDF5Attribute::new(name, data_type, data)
}
