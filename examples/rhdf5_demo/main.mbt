///|
/// rhdf5 Demo - HDF5 file operations demo.

pub fn main() {
  println("=== rhdf5 Demo ===")
  
  let file = @src.h5create_file("genomic_data.h5")
  
  let file1 = @src.h5write_dataset(file, "/genome/chr1", "ACGTACGTACGT", @src.H5T_NATIVE_STRING, [12])
  let file2 = @src.h5write_dataset(file1, "/genome/chr2", "GCTAGCTAGCTA", @src.H5T_NATIVE_STRING, [12])
  let file3 = @src.h5write_dataset(file2, "/expression/data", "1.2,3.4,5.6,7.8,9.0", @src.H5T_NATIVE_DOUBLE, [5])
  
  println("HDF5 file created: genomic_data.h5")
  println("Contents:")
  
  let listing = @src.h5ls(file3)
  let mut i = 0
  while i < listing.length() {
    println("  " + listing[i])
    i = i + 1
  }
  
  let read_result = @src.h5read_dataset(file3, "/expression/data")
  if read_result is Some(_) {
    println("\nRead expression data: " + read_result.unwrap())
  }
  
  let example_file = @src.h5create_example_file()
  println("\nExample file created with genome, expression, and variants data")
}
