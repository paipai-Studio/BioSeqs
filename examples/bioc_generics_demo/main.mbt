///|
/// BiocGenerics Demo - Generic functions demo.

pub fn main() {
  println("=== BiocGenerics Demo ===")
  
  println("\n1. NA handling:")
  let na_i = @src.na_int()
  println("NA int: " + na_i.to_string())
  println("Is NA: " + @src.is_na_int(na_i).to_string())
  
  println("\n2. Sorting:")
  let sorted = @src.sort_int([5, 2, 8, 1, 9])
  println("Sorted [5,2,8,1,9]: " + int_array_to_string(sorted))
  
  println("\n3. Unique values:")
  let unique = @src.unique_int([1, 2, 1, 3, 2])
  println("Unique [1,2,1,3,2]: " + int_array_to_string(unique))
  
  println("\n4. Set operations:")
  let inter = @src.intersect_int([1, 2, 3], [2, 3, 4])
  println("Intersection: " + int_array_to_string(inter))
  
  let uni = @src.union_int([1, 2], [3, 4])
  println("Union: " + int_array_to_string(uni))
  
  let diff = @src.setdiff_int([1, 2, 3], [2])
  println("Setdiff: " + int_array_to_string(diff))
  
  println("\n5. Ordering:")
  let indices = @src.order_int([3, 1, 4, 2])
  println("Order indices: " + int_array_to_string(indices))
  
  println("\n6. Table:")
  let table = @src.table_int([1, 2, 1, 3, 1, 2])
  println("Table count for '1': " + table.get("1").unwrap().to_string())
  
  println("\n7. Sequence generation:")
  let seq = @src.seq_int(1, 10, by=2)
  println("seq(1, 10, by=2): " + int_array_to_string(seq))
  
  println("\n8. Head/Tail:")
  let arr = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
  let hd = @src.head_int(arr, 3)
  let tl = @src.tail_int(arr, 3)
  println("Head: " + int_array_to_string(hd))
  println("Tail: " + int_array_to_string(tl))
  
  println("\n9. Reverse:")
  let rev_arr = @src.rev_int([1, 2, 3])
  println("Reverse [1,2,3]: " + int_array_to_string(rev_arr))
  
  println("\n10. Repeat:")
  let rep_arr = @src.rep_int(5, 4)
  println("Repeat 5, 4 times: " + int_array_to_string(rep_arr))
  
  println("\n11. Match:")
  let matched = @src.match_int([2, 4], [1, 2, 3])
  println("Match [2,4] in [1,2,3]: " + int_array_to_string(matched))
  
  println("\n12. NA omit:")
  let cleaned = @src.na_omit_int([1, -2147483648, 3, 4])
  println("NA omit [1, NA, 3, 4]: " + int_array_to_string(cleaned))
}

fn int_array_to_string(arr : Array[Int]) -> String {
  let mut result = "["
  let mut i = 0
  while i < arr.length() {
    if i > 0 {
      result = result + ", "
    }
    result = result + arr[i].to_string()
    i = i + 1
  }
  result = result + "]"
  result
}
