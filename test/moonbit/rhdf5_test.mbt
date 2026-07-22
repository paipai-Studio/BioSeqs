///|
/// Tests for rhdf5 module.

test "HDF5Attribute creation" {
  let attr = @src.HDF5Attribute::new("species", "H5T_NATIVE_STRING", "Homo sapiens")
  assert_eq(attr.name, "species")
}

test "HDF5Dataset creation" {
  let ds = @src.HDF5Dataset::new("expression", "H5T_NATIVE_DOUBLE", [100, 10])
  assert_eq(ds.name, "expression")
  assert_eq(ds.dimensions.length(), 2)
}

test "HDF5Group creation" {
  let group = @src.HDF5Group::new("/genome")
  assert_eq(group.name, "/genome")
}

test "HDF5File creation" {
  let file = @src.h5create_file("test.h5")
  assert_eq(file.filename, "test.h5")
}

test "h5create_dataset simple" {
  let file = @src.h5create_file("test.h5")
  let file2 = @src.h5create_dataset(file, "/data/matrix", "H5T_NATIVE_DOUBLE", [3, 3])
  let result = @src.h5read_dataset(file2, "/data/matrix")
  assert_true(result is Some(_))
}

test "h5write_dataset simple" {
  let file = @src.h5create_file("test.h5")
  let file2 = @src.h5write_dataset(file, "/seq/dna", "ACGTACGT", "H5T_NATIVE_STRING", [8])
  let result = @src.h5read_dataset(file2, "/seq/dna")
  assert_true(result is Some(_))
}

test "h5read_dataset" {
  let file = @src.h5create_file("test.h5")
  let file2 = @src.h5write_dataset(file, "/data/values", "1.0,2.0,3.0", "H5T_NATIVE_DOUBLE", [3])
  let result = @src.h5read_dataset(file2, "/data/values")
  assert_true(result is Some(_))
}

test "h5ls" {
  let file = @src.h5create_example_file()
  let listing = @src.h5ls(file)
  assert_true(listing.length() > 0)
}

test "h5create_example_file" {
  let file = @src.h5create_example_file()
  assert_eq(file.filename, "example.h5")
}

test "HDF5Dataset add_attribute" {
  let ds = @src.HDF5Dataset::new("test", "H5T_NATIVE_INT", [10])
  let attr = @src.HDF5Attribute::new("unit", "H5T_NATIVE_STRING", "counts")
  let ds_with_attr = ds.add_attribute(attr)
  assert_true(ds_with_attr.get_attribute("unit") is Some(_))
}

test "HDF5Group find_group" {
  let group = @src.HDF5Group::new("test")
  let result = group.find_group("data")
  assert_true(result is None)
}

test "HDF5Group find_dataset" {
  let group = @src.HDF5Group::new("test")
  let result = group.find_dataset("data")
  assert_true(result is None)
}
