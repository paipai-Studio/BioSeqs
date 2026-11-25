///|
/// biotite_structure_annotation.mbt — residue/SSE annotation helpers.

///|
/// Per-residue secondary structure string (PSEA 3-state: H/E/C).
pub fn biotite_annotate_sse(arr : AtomArray) -> String {
  arr.annotate_sse().psea_ss_string()
}

///|
/// Per-atom flag: true when the atom belongs to a canonical amino acid.
pub fn biotite_annotate_amino_acids(arr : AtomArray) -> Array[Bool] {
  let result : Array[Bool] = Array(capacity=arr.len())
  for i in 0..<arr.len() {
    result.push(biotite_is_amino_acid_name(arr.get_res_name(i)))
  }
  result
}

///|
/// Per-atom residue descriptor: (chain ID, residue ID, residue name).
pub fn biotite_annotate_residues(
  arr : AtomArray,
) -> Array[(String, Int32, String)] {
  let result : Array[(String, Int32, String)] = Array(capacity=arr.len())
  for i in 0..<arr.len() {
    result.push(
      (arr.get_chain_id(i).to_string(), arr.res_ids[i], arr.get_res_name(i)),
    )
  }
  result
}

///|
/// SSE assignment per residue: (residue number, residue name, SSE code).
pub fn biotite_annotate_sse_by_residue(
  arr : AtomArray,
) -> Array[(Int32, String, String)] {
  let sse = arr.annotate_sse()
  let result : Array[(Int32, String, String)] = Array(
    capacity=sse.psea_n_residues(),
  )
  for i in 0..<sse.psea_n_residues() {
    result.push(
      (sse.psea_get_res_num(i), sse.psea_get_res_name(i), sse.psea_get_ss(i)),
    )
  }
  result
}

///|
fn biotite_is_amino_acid_name(res_name : String) -> Bool {
  match res_name {
    "ALA"
    | "ARG"
    | "ASN"
    | "ASP"
    | "CYS"
    | "GLN"
    | "GLU"
    | "GLY"
    | "HIS"
    | "ILE"
    | "LEU"
    | "LYS"
    | "MET"
    | "PHE"
    | "PRO"
    | "SER"
    | "THR"
    | "TRP"
    | "TYR"
    | "VAL" => true
    _ => false
  }
}
