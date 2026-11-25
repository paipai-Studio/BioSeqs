///|
/// biotite_sequence_bundles.mbt — a simple named bundle container for
/// sequences and structures (Biotite Bundles analogue).

///|
pub(all) struct BiotiteBundles {
  sequence_keys : Array[String]
  sequences : Array[String]
  structure_keys : Array[String]
  structures : Array[AtomArray]
}

///|
pub fn biotite_bundles_new() -> BiotiteBundles {
  { sequence_keys: [], sequences: [], structure_keys: [], structures: [], }
}

///|
/// Add or replace a named sequence.
pub fn biotite_bundles_add_sequence(
  bundles : BiotiteBundles,
  name : String,
  sequence : String,
) -> BiotiteBundles {
  let keys : Array[String] = Array(capacity=bundles.sequence_keys.length() + 1)
  let seqs : Array[String] = Array(capacity=bundles.sequences.length() + 1)
  let mut replaced = false
  for i in 0..<bundles.sequence_keys.length() {
    keys.push(bundles.sequence_keys[i])
    if bundles.sequence_keys[i] == name {
      seqs.push(sequence)
      replaced = true
    } else {
      seqs.push(bundles.sequences[i])
    }
  }
  if !replaced {
    keys.push(name)
    seqs.push(sequence)
  }
  {
    sequence_keys: keys,
    sequences: seqs,
    structure_keys: bundles.structure_keys,
    structures: bundles.structures,
  }
}

///|
/// Add or replace a named structure.
pub fn biotite_bundles_add_structure(
  bundles : BiotiteBundles,
  name : String,
  structure : AtomArray,
) -> BiotiteBundles {
  let keys : Array[String] = Array(capacity=bundles.structure_keys.length() + 1)
  let values : Array[AtomArray] = Array(
    capacity=bundles.structures.length() + 1,
  )
  let mut replaced = false
  for i in 0..<bundles.structure_keys.length() {
    keys.push(bundles.structure_keys[i])
    if bundles.structure_keys[i] == name {
      values.push(structure)
      replaced = true
    } else {
      values.push(bundles.structures[i])
    }
  }
  if !replaced {
    keys.push(name)
    values.push(structure)
  }
  {
    sequence_keys: bundles.sequence_keys,
    sequences: bundles.sequences,
    structure_keys: keys,
    structures: values,
  }
}

///|
pub fn biotite_bundles_get_sequence(
  bundles : BiotiteBundles,
  name : String,
) -> String? {
  for i in 0..<bundles.sequence_keys.length() {
    if bundles.sequence_keys[i] == name {
      return Some(bundles.sequences[i])
    }
  }
  None
}

///|
pub fn biotite_bundles_get_structure(
  bundles : BiotiteBundles,
  name : String,
) -> AtomArray? {
  for i in 0..<bundles.structure_keys.length() {
    if bundles.structure_keys[i] == name {
      return Some(bundles.structures[i])
    }
  }
  None
}

///|
pub fn biotite_bundles_length(bundles : BiotiteBundles) -> Int32 {
  bundles.sequence_keys.length()
}

///|
pub fn biotite_bundles_structure_length(bundles : BiotiteBundles) -> Int32 {
  bundles.structure_keys.length()
}
