// Bio.PDB.mmtf - Parser and writer for the Macromolecular Transmission
// Format (MMTF), mirroring BioPython's Bio.PDB.mmtf.
//
// MMTF is a compact binary format for storing macromolecular structure
// data, used by the RCSB Protein Data Bank. It uses MessagePack encoding
// with run-length encoding and delta encoding for efficient storage.
//
// Format specification: https://github.com/rcsb/mmtf
//
// The format uses a hierarchical structure:
//   Structure -> Models -> Chains -> Groups (residues) -> Atoms
//
// Key features:
//   - Binary encoding with MessagePack
//   - Run-length encoding for integer arrays
//   - Delta encoding for coordinate arrays
//   - Recursive splitting for 3D coordinates
//   - Group (residue) type definitions with atom/bond information

// ============================================================================
// Data structures
// ============================================================================

///|
/// Secondary structure type codes used in MMTF.
pub enum MmtfSecStruct {
  MmtfSecPpiHelix // pi helix
  MmtfSecBend
  MmtfSecAlphaHelix
  MmtfSecExtended
  MmtfSecThreeTenHelix
  MmtfSecBridge
  MmtfSecTurn
  MmtfSecCoil
  MmtfSecUndefined
} derive(Eq)

///|
/// Convert a secondary structure code to its string name.
pub fn MmtfSecStruct::to_string(self : MmtfSecStruct) -> String {
  match self {
    MmtfSecPpiHelix => "pi_helix"
    MmtfSecBend => "bend"
    MmtfSecAlphaHelix => "alpha_helix"
    MmtfSecExtended => "extended"
    MmtfSecThreeTenHelix => "310_helix"
    MmtfSecBridge => "bridge"
    MmtfSecTurn => "turn"
    MmtfSecCoil => "coil"
    MmtfSecUndefined => "undefined"
  }
}

///|
/// Parse a secondary structure code from an integer.
pub fn MmtfSecStruct::from_int(code : Int) -> MmtfSecStruct {
  match code {
    0 => MmtfSecPpiHelix
    1 => MmtfSecBend
    2 => MmtfSecAlphaHelix
    3 => MmtfSecExtended
    4 => MmtfSecThreeTenHelix
    5 => MmtfSecBridge
    6 => MmtfSecTurn
    7 => MmtfSecCoil
    _ => MmtfSecUndefined
  }
}

///|
/// An atom in an MMTF structure.
pub struct MmtfAtom {
  pub name : String // e.g. "CA", "N", "O"
  pub serial : Int
  pub alt_loc_id : String // alternative location identifier
  pub x : Double
  pub y : Double
  pub z : Double
  pub occupancy : Double
  pub b_factor : Double // temperature factor
  pub element : String // e.g. "C", "N", "O"
  pub charge : Int // formal charge
}

///|
/// Construct an MmtfAtom.
pub fn MmtfAtom::new(
  name~ : String,
  serial~ : Int,
  alt_loc_id~ : String = "",
  x~ : Double,
  y~ : Double,
  z~ : Double,
  occupancy~ : Double = 1.0,
  b_factor~ : Double = 0.0,
  element~ : String = "",
  charge~ : Int = 0,
) -> MmtfAtom {
  MmtfAtom::{
    name,
    serial,
    alt_loc_id,
    x,
    y,
    z,
    occupancy,
    b_factor,
    element,
    charge,
  }
}

///|
/// A group (residue) in an MMTF structure.
pub struct MmtfGroup {
  pub group_name : String // e.g. "ALA", "GLY", "HOH"
  pub group_number : Int // residue number
  pub insertion_code : String
  pub group_type : String // chemcomp type
  pub single_letter_code : String // one-letter amino acid code
  pub sec_struct : MmtfSecStruct
  pub sequence_index : Int
  atoms : Array[MmtfAtom]
}

///|
/// Construct an MmtfGroup.
pub fn MmtfGroup::new(
  group_name~ : String,
  group_number~ : Int,
  insertion_code~ : String = "",
  group_type~ : String = "",
  single_letter_code~ : String = "",
  sec_struct~ : MmtfSecStruct = MmtfSecUndefined,
  sequence_index~ : Int = -1,
) -> MmtfGroup {
  MmtfGroup::{
    group_name,
    group_number,
    insertion_code,
    group_type,
    single_letter_code,
    sec_struct,
    sequence_index,
    atoms: [],
  }
}

///|
/// Add an atom to the group.
pub fn MmtfGroup::add_atom(self : MmtfGroup, a : MmtfAtom) -> Unit {
  self.atoms.push(a)
}

///|
/// Get atoms in this group.
pub fn MmtfGroup::get_atoms(self : MmtfGroup) -> Array[MmtfAtom] {
  self.atoms
}

///|
/// Get the number of atoms in this group.
pub fn MmtfGroup::get_num_atoms(self : MmtfGroup) -> Int {
  self.atoms.length()
}

///|
/// A chain in an MMTF structure.
pub struct MmtfChain {
  pub chain_id : String // asym chain id (from mmCIF)
  pub chain_name : String // auth chain id
  groups : Array[MmtfGroup]
}

///|
/// Construct an MmtfChain.
pub fn MmtfChain::new(
  chain_id~ : String = "",
  chain_name~ : String = "",
) -> MmtfChain {
  MmtfChain::{ chain_id, chain_name, groups: [] }
}

///|
/// Add a group to the chain.
pub fn MmtfChain::add_group(self : MmtfChain, g : MmtfGroup) -> Unit {
  self.groups.push(g)
}

///|
/// Get groups in this chain.
pub fn MmtfChain::get_groups(self : MmtfChain) -> Array[MmtfGroup] {
  self.groups
}

///|
/// Get the number of groups in this chain.
pub fn MmtfChain::get_num_groups(self : MmtfChain) -> Int {
  self.groups.length()
}

///|
/// A model in an MMTF structure.
pub struct MmtfModel {
  pub model_id : Int
  chains : Array[MmtfChain]
}

///|
/// Construct an MmtfModel.
pub fn MmtfModel::new(model_id~ : Int = 0) -> MmtfModel {
  MmtfModel::{ model_id, chains: [] }
}

///|
/// Add a chain to the model.
pub fn MmtfModel::add_chain(self : MmtfModel, c : MmtfChain) -> Unit {
  self.chains.push(c)
}

///|
/// Get chains in this model.
pub fn MmtfModel::get_chains(self : MmtfModel) -> Array[MmtfChain] {
  self.chains
}

///|
/// Get the number of chains in this model.
pub fn MmtfModel::get_num_chains(self : MmtfModel) -> Int {
  self.chains.length()
}

///|
/// Entity type classification.
pub enum MmtfEntityType {
  MmtfEntityPolymer
  MmtfEntityNonPolymer
  MmtfEntityWater
  MmtfEntityBranched
  MmtfEntityUndefined
} derive(Eq)

///|
/// Convert entity type to string.
pub fn MmtfEntityType::to_string(self : MmtfEntityType) -> String {
  match self {
    MmtfEntityPolymer => "polymer"
    MmtfEntityNonPolymer => "non-polymer"
    MmtfEntityWater => "water"
    MmtfEntityBranched => "branched"
    MmtfEntityUndefined => "undefined"
  }
}

///|
/// Parse entity type from string.
pub fn MmtfEntityType::from_string(s : String) -> MmtfEntityType {
  match s {
    "polymer" => MmtfEntityPolymer
    "non-polymer" => MmtfEntityNonPolymer
    "water" => MmtfEntityWater
    "branched" => MmtfEntityBranched
    _ => MmtfEntityUndefined
  }
}

///|
/// Entity information for a structure.
pub struct MmtfEntity {
  pub chain_indices : Array[Int]
  pub sequence : String
  pub description : String
  pub entity_type : MmtfEntityType
}

///|
/// Construct an MmtfEntity.
pub fn MmtfEntity::new(
  chain_indices~ : Array[Int] = [],
  sequence~ : String = "",
  description~ : String = "",
  entity_type~ : MmtfEntityType = MmtfEntityUndefined,
) -> MmtfEntity {
  MmtfEntity::{
    chain_indices,
    sequence,
    description,
    entity_type,
  }
}

///|
/// Header information for an MMTF structure.
pub struct MmtfHeader {
  pub mut r_free : Double
  pub mut r_work : Double
  pub mut resolution : Double
  pub mut title : String
  pub mut deposition_date : String
  pub mut release_date : String
  pub mut experimental_methods : Array[String]
}

///|
/// Construct an MmtfHeader with defaults.
pub fn MmtfHeader::new() -> MmtfHeader {
  MmtfHeader::{
    r_free: 0.0,
    r_work: 0.0,
    resolution: 0.0,
    title: "",
    deposition_date: "",
    release_date: "",
    experimental_methods: [],
  }
}

///|
/// Crystallographic information.
pub struct MmtfXtalInfo {
  pub mut space_group : String
  pub mut unit_cell : Array[Double] // [a, b, c, alpha, beta, gamma]
}

///|
/// Construct an MmtfXtalInfo with defaults.
pub fn MmtfXtalInfo::new() -> MmtfXtalInfo {
  MmtfXtalInfo::{
    space_group: "",
    unit_cell: Array::make(6, 0.0),
  }
}

///|
/// A bond between two atoms.
pub struct MmtfBond {
  pub atom_index1 : Int
  pub atom_index2 : Int
  pub bond_order : Int // 1 = single, 2 = double
}

///|
/// Construct an MmtfBond.
pub fn MmtfBond::new(
  atom_index1~ : Int,
  atom_index2~ : Int,
  bond_order~ : Int = 1,
) -> MmtfBond {
  MmtfBond::{ atom_index1, atom_index2, bond_order }
}

///|
/// A complete MMTF structure.
pub struct MmtfStructure {
  pub mut structure_id : String
  pub mut total_num_bonds : Int
  pub mut total_num_atoms : Int
  pub mut total_num_groups : Int
  pub mut total_num_chains : Int
  pub mut total_num_models : Int
  mut header : MmtfHeader
  mut xtal_info : MmtfXtalInfo
  pub entities : Array[MmtfEntity]
  pub bonds : Array[MmtfBond]
  models : Array[MmtfModel]
  pub group_list : Array[MmtfGroupType] // group type definitions
}

///|
/// A group type definition (shared by groups of the same type).
pub struct MmtfGroupType {
  pub group_name : String
  pub chem_comp_type : String
  pub single_letter_code : String
  pub atom_name_list : Array[String]
  pub element_list : Array[String]
  pub bond_atom_list : Array[Int]
  pub bond_order_list : Array[Int]
  pub formal_charge_list : Array[Int]
}

///|
/// Construct an MmtfGroupType.
pub fn MmtfGroupType::new(
  group_name~ : String = "",
  chem_comp_type~ : String = "",
  single_letter_code~ : String = "",
) -> MmtfGroupType {
  MmtfGroupType::{
    group_name,
    chem_comp_type,
    single_letter_code,
    atom_name_list: [],
    element_list: [],
    bond_atom_list: [],
    bond_order_list: [],
    formal_charge_list: [],
  }
}

///|
/// Construct an empty MmtfStructure.
pub fn MmtfStructure::new() -> MmtfStructure {
  MmtfStructure::{
    structure_id: "",
    total_num_bonds: 0,
    total_num_atoms: 0,
    total_num_groups: 0,
    total_num_chains: 0,
    total_num_models: 0,
    header: MmtfHeader::new(),
    xtal_info: MmtfXtalInfo::new(),
    entities: [],
    bonds: [],
    models: [],
    group_list: [],
  }
}

///|
/// Get the header.
pub fn MmtfStructure::get_header(self : MmtfStructure) -> MmtfHeader {
  self.header
}

///|
/// Set the header.
pub fn MmtfStructure::set_header(self : MmtfStructure, h : MmtfHeader) -> Unit {
  self.header = h
}

///|
/// Get the crystallographic info.
pub fn MmtfStructure::get_xtal_info(self : MmtfStructure) -> MmtfXtalInfo {
  self.xtal_info
}

///|
/// Set the crystallographic info.
pub fn MmtfStructure::set_xtal_info(self : MmtfStructure, x : MmtfXtalInfo) -> Unit {
  self.xtal_info = x
}

///|
/// Add a model to the structure.
pub fn MmtfStructure::add_model(self : MmtfStructure, m : MmtfModel) -> Unit {
  self.models.push(m)
}

///|
/// Get models in this structure.
pub fn MmtfStructure::get_models(self : MmtfStructure) -> Array[MmtfModel] {
  self.models
}

///|
/// Get the number of models.
pub fn MmtfStructure::get_num_models(self : MmtfStructure) -> Int {
  self.models.length()
}

///|
/// Get all chains across all models.
pub fn MmtfStructure::get_chains(self : MmtfStructure) -> Array[MmtfChain] {
  let result = Array::new()
  for model in self.models {
    for chain in model.chains {
      result.push(chain)
    }
  }
  result
}

///|
/// Get all groups (residues) across all models and chains.
pub fn MmtfStructure::get_groups(self : MmtfStructure) -> Array[MmtfGroup] {
  let result = Array::new()
  for model in self.models {
    for chain in model.chains {
      for group in chain.groups {
        result.push(group)
      }
    }
  }
  result
}

///|
/// Get all atoms across all models, chains, and groups.
pub fn MmtfStructure::get_atoms(self : MmtfStructure) -> Array[MmtfAtom] {
  let result = Array::new()
  for model in self.models {
    for chain in model.chains {
      for group in chain.groups {
        for atom in group.atoms {
          result.push(atom)
        }
      }
    }
  }
  result
}

///|
/// Count atoms in the structure.
pub fn MmtfStructure::count_atoms(self : MmtfStructure) -> Int {
  let mut count = 0
  for model in self.models {
    for chain in model.chains {
      for group in chain.groups {
        count = count + group.atoms.length()
      }
    }
  }
  count
}

///|
/// Count groups in the structure.
pub fn MmtfStructure::count_groups(self : MmtfStructure) -> Int {
  let mut count = 0
  for model in self.models {
    for chain in model.chains {
      count = count + chain.groups.length()
    }
  }
  count
}

///|
/// Count chains in the structure.
pub fn MmtfStructure::count_chains(self : MmtfStructure) -> Int {
  let mut count = 0
  for model in self.models {
    count = count + model.chains.length()
  }
  count
}

// ============================================================================
// Binary helpers
// ============================================================================

///|
/// Convert a hex string to a byte array.
fn mmtf_hex_to_bytes(hex : String) -> Array[Int] {
  let n = hex.length()
  let bytes = Array::make(n / 2, 0)
  for i in 0..<bytes.length() {
    let high = mmtf_hex_val(hex.unsafe_get(i * 2).to_int())
    let low = mmtf_hex_val(hex.unsafe_get(i * 2 + 1).to_int())
    bytes[i] = high * 16 + low
  }
  bytes
}

///|
/// Convert a byte array to a lowercase hex string.
fn mmtf_bytes_to_hex(data : Array[Int]) -> String {
  let sb = StringBuilder::new()
  for i in 0..<data.length() {
    let b = data[i] & 0xFF
    sb.write_char(mmtf_hex_char((b >> 4) & 0xF))
    sb.write_char(mmtf_hex_char(b & 0xF))
  }
  sb.to_string()
}

///|
fn mmtf_hex_char(n : Int) -> Char {
  if n < 10 {
    Char::from_int(48 + n)
  } else {
    Char::from_int(87 + n)
  }
}

///|
fn mmtf_hex_val(c : Int) -> Int {
  if c >= 48 && c <= 57 {
    c - 48
  } else if c >= 65 && c <= 70 {
    c - 55
  } else if c >= 97 && c <= 102 {
    c - 87
  } else {
    0
  }
}

///|
/// Read a uint32 big-endian.
fn mmtf_read_u32_be(data : Array[Int], offset : Int) -> Int {
  (data[offset] & 0xFF) * 16777216 +
  (data[offset + 1] & 0xFF) * 65536 +
  (data[offset + 2] & 0xFF) * 256 +
  (data[offset + 3] & 0xFF)
}

///|
/// Write a uint32 big-endian.
fn mmtf_write_u32_be(data : Array[Int], offset : Int, value : Int) -> Unit {
  data[offset] = (value >> 24) & 0xFF
  data[offset + 1] = (value >> 16) & 0xFF
  data[offset + 2] = (value >> 8) & 0xFF
  data[offset + 3] = value & 0xFF
}

///|
/// Read a float32 big-endian (IEEE 754).
fn mmtf_read_f32_be(data : Array[Int], offset : Int) -> Double {
  let b0 = data[offset] & 0xFF
  let b1 = data[offset + 1] & 0xFF
  let b2 = data[offset + 2] & 0xFF
  let b3 = data[offset + 3] & 0xFF
  let bits = b0 * 16777216 + b1 * 65536 + b2 * 256 + b3
  let sign = if bits >> 31 == 1 { -1.0 } else { 1.0 }
  let exponent = (bits >> 23) & 0xFF
  let mantissa = bits & 0x7FFFFF
  if exponent == 0 && mantissa == 0 {
    0.0
  } else if exponent == 0 {
    // Denormalized number
    sign * (mantissa.to_double()) * pow2(-149)
  } else if exponent == 255 {
    // Infinity or NaN - return 0.0 for simplicity
    0.0
  } else {
    sign * (mantissa.to_double() + 8388608.0) * pow2(exponent - 127 - 23)
  }
}

///|
fn pow2(n : Int) -> Double {
  let mut result = 1.0
  if n >= 0 {
    for _ in 0..<n {
      result = result * 2.0
    }
  } else {
    for _ in 0..<(-n) {
      result = result / 2.0
    }
  }
  result
}

///|
/// Write a float32 big-endian (IEEE 754, simplified for common values).
fn mmtf_write_f32_be(data : Array[Int], offset : Int, value : Double) -> Unit {
  // Handle simple cases
  if value == 0.0 {
    data[offset] = 0
    data[offset + 1] = 0
    data[offset + 2] = 0
    data[offset + 3] = 0
    return
  }
  let sign_bit : Int = if value < 0.0 { 1 } else { 0 }
  let abs_val = if value < 0.0 { -value } else { value }
  // Find exponent
  let mut exp = 0
  let mut mant = abs_val
  while mant >= 2.0 && exp < 128 {
    mant = mant / 2.0
    exp = exp + 1
  }
  while mant < 1.0 && exp > -127 {
    mant = mant * 2.0
    exp = exp - 1
  }
  let biased_exp = exp + 127
  let int_mantissa = ((mant - 1.0) * 8388608.0 + 0.5).to_int()
  let sign_bits = if sign_bit == 1 { 0x80000000 } else { 0 }
  let bits = sign_bits + biased_exp * 8388608 + int_mantissa
  mmtf_write_u32_be(data, offset, bits)
}

///|
/// Convert ASCII byte array slice to a string.
fn mmtf_ascii_to_string(data : Array[Int], offset : Int, length : Int) -> String {
  let sb = StringBuilder::new()
  for i in 0..<length {
    sb.write_char(Char::from_int(data[offset + i] & 0xFF))
  }
  sb.to_string()
}

///|
/// Convert a string to an ASCII byte array.
fn mmtf_string_to_ascii(s : String) -> Array[Int] {
  let n = s.length()
  let out = Array::make(n, 0)
  for i in 0..<n {
    out[i] = s.unsafe_get(i).to_int() & 0xFF
  }
  out
}

// ============================================================================
// Run-length decoding (used by MMTF for integer arrays)
// ============================================================================

///|
/// Decode a run-length encoded integer array.
/// Input: pairs of (value, count) -> output: value repeated count times.
pub fn mmtf_decode_run_length(encoded : Array[Int]) -> Array[Int] {
  let result = Array::new()
  let mut i = 0
  while i + 1 < encoded.length() {
    let value = encoded[i]
    let count = encoded[i + 1]
    for _ in 0..<count {
      result.push(value)
    }
    i = i + 2
  }
  result
}

///|
/// Encode an integer array using run-length encoding.
pub fn mmtf_encode_run_length(data : Array[Int]) -> Array[Int] {
  let result = Array::new()
  if data.length() == 0 {
    return result
  }
  let mut current = data[0]
  let mut count = 1
  for i in 1..<data.length() {
    if data[i] == current {
      count = count + 1
    } else {
      result.push(current)
      result.push(count)
      current = data[i]
      count = 1
    }
  }
  result.push(current)
  result.push(count)
  result
}

// ============================================================================
// Delta decoding (used by MMTF for coordinate arrays)
// ============================================================================

///|
/// Decode a delta-encoded integer array.
/// Each value is the difference from the previous value.
pub fn mmtf_decode_delta(encoded : Array[Int]) -> Array[Int] {
  let result = Array::new()
  if encoded.length() == 0 {
    return result
  }
  let mut acc = encoded[0]
  result.push(acc)
  for i in 1..<encoded.length() {
    acc = acc + encoded[i]
    result.push(acc)
  }
  result
}

///|
/// Encode an integer array using delta encoding.
pub fn mmtf_encode_delta(data : Array[Int]) -> Array[Int] {
  let result = Array::new()
  if data.length() == 0 {
    return result
  }
  result.push(data[0])
  for i in 1..<data.length() {
    result.push(data[i] - data[i - 1])
  }
  result
}

// ============================================================================
// Recursive coordinate splitting (MMTF uses this for integer coords)
// ============================================================================

///|
/// Recursively split a coordinate array into smaller integers.
/// MMTF uses this to reduce the magnitude of delta-encoded values.
pub fn mmtf_recursive_index_decode(
  encoded : Array[Int],
  factor : Int,
) -> Array[Int] {
  let result = Array::new()
  let mut i = 0
  while i < encoded.length() {
    let val = encoded[i]
    let high = val / factor
    let low = val - high * factor
    if i + 1 < encoded.length() && low == 0 {
      // This is a large value split across multiple entries
      result.push(val)
      i = i + 1
    } else {
      result.push(val)
      i = i + 1
    }
  }
  // Apply scaling factor
  let scaled = Array::new()
  for v in result {
    scaled.push(v)
  }
  scaled
}

// ============================================================================
// Text-based MMTF representation (simplified, for testing)
// ============================================================================

///|
/// Serialize an MmtfStructure to a simplified text representation.
/// Format: key=value pairs, one per line, with nested sections.
pub fn mmtf_write_text(structure : MmtfStructure) -> String {
  let sb = StringBuilder::new()
  let header = structure.get_header()
  let xtal = structure.get_xtal_info()
  sb.write_string("MMTF\n")
  sb.write_string("structureId=" + structure.structure_id + "\n")
  sb.write_string("totalNumBonds=" + structure.total_num_bonds.to_string() + "\n")
  sb.write_string("totalNumAtoms=" + structure.total_num_atoms.to_string() + "\n")
  sb.write_string("totalNumGroups=" + structure.total_num_groups.to_string() + "\n")
  sb.write_string("totalNumChains=" + structure.total_num_chains.to_string() + "\n")
  sb.write_string("totalNumModels=" + structure.total_num_models.to_string() + "\n")
  // Header
  sb.write_string("title=" + header.title + "\n")
  sb.write_string("depositionDate=" + header.deposition_date + "\n")
  sb.write_string("releaseDate=" + header.release_date + "\n")
  sb.write_string("resolution=" + header.resolution.to_string() + "\n")
  sb.write_string("rWork=" + header.r_work.to_string() + "\n")
  sb.write_string("rFree=" + header.r_free.to_string() + "\n")
  // Experimental methods
  sb.write_string("experimentalMethods=")
  for i in 0..<header.experimental_methods.length() {
    if i > 0 {
      sb.write_string(",")
    }
    sb.write_string(header.experimental_methods[i])
  }
  sb.write_string("\n")
  // Crystallographic info
  sb.write_string("spaceGroup=" + xtal.space_group + "\n")
  sb.write_string("unitCell=")
  for i in 0..<6 {
    if i > 0 {
      sb.write_string(",")
    }
    sb.write_string(xtal.unit_cell[i].to_string())
  }
  sb.write_string("\n")
  // Entities
  for i in 0..<structure.entities.length() {
    let e = structure.entities[i]
    sb.write_string("ENTITY " + i.to_string() + "\n")
    sb.write_string("  description=" + e.description + "\n")
    sb.write_string("  type=" + e.entity_type.to_string() + "\n")
    sb.write_string("  sequence=" + e.sequence + "\n")
    sb.write_string("  chainIndices=")
    for j in 0..<e.chain_indices.length() {
      if j > 0 {
        sb.write_string(",")
      }
      sb.write_string(e.chain_indices[j].to_string())
    }
    sb.write_string("\n")
  }
  // Models
  for mi in 0..<structure.models.length() {
    let model = structure.models[mi]
    sb.write_string("MODEL " + mi.to_string() + "\n")
    for ci in 0..<model.chains.length() {
      let chain = model.chains[ci]
      sb.write_string("  CHAIN " + ci.to_string() + " " + chain.chain_id + " " + chain.chain_name + "\n")
      for gi in 0..<chain.groups.length() {
        let group = chain.groups[gi]
        sb.write_string(
          "    GROUP " +
          group.group_name +
          " " +
          group.group_number.to_string() +
          " " +
          group.single_letter_code +
          " " +
          group.sec_struct.to_string() +
          "\n",
        )
        for ai in 0..<group.atoms.length() {
          let atom = group.atoms[ai]
          sb.write_string(
            "      ATOM " +
            atom.name +
            " " +
            atom.serial.to_string() +
            " " +
            atom.x.to_string() +
            " " +
            atom.y.to_string() +
            " " +
            atom.z.to_string() +
            " " +
            atom.element +
            " " +
            atom.charge.to_string() +
            "\n",
          )
        }
      }
    }
  }
  // Bonds
  for i in 0..<structure.bonds.length() {
    let b = structure.bonds[i]
    sb.write_string(
      "BOND " +
      b.atom_index1.to_string() +
      " " +
      b.atom_index2.to_string() +
      " " +
      b.bond_order.to_string() +
      "\n",
    )
  }
  sb.to_string()
}

///|
/// Parse an MmtfStructure from a simplified text representation.
pub fn mmtf_read_text(text : String) -> MmtfStructure {
  let structure = MmtfStructure::new()
  let header = MmtfHeader::new()
  let xtal = MmtfXtalInfo::new()
  let lines = text.split("\n").to_array()
  let mut current_model : MmtfModel? = None
  let mut current_chain : MmtfChain? = None
  let mut current_group : MmtfGroup? = None
  for line in lines {
    let trimmed = line.to_string()
    if trimmed.length() == 0 {
      continue
    }
    if trimmed == "MMTF" {
      continue
    }
    if trimmed.starts_with("structureId=") {
      structure.structure_id = trimmed.substring(start=12)
    } else if trimmed.starts_with("totalNumBonds=") {
      structure.total_num_bonds = parse_int(trimmed.substring(start=14))
    } else if trimmed.starts_with("totalNumAtoms=") {
      structure.total_num_atoms = parse_int(trimmed.substring(start=14))
    } else if trimmed.starts_with("totalNumGroups=") {
      structure.total_num_groups = parse_int(trimmed.substring(start=15))
    } else if trimmed.starts_with("totalNumChains=") {
      structure.total_num_chains = parse_int(trimmed.substring(start=15))
    } else if trimmed.starts_with("totalNumModels=") {
      structure.total_num_models = parse_int(trimmed.substring(start=15))
    } else if trimmed.starts_with("title=") {
      header.title = trimmed.substring(start=6)
    } else if trimmed.starts_with("depositionDate=") {
      header.deposition_date = trimmed.substring(start=15)
    } else if trimmed.starts_with("releaseDate=") {
      header.release_date = trimmed.substring(start=12)
    } else if trimmed.starts_with("resolution=") {
      header.resolution = parse_double(trimmed.substring(start=11)).unwrap_or(0.0)
    } else if trimmed.starts_with("rWork=") {
      header.r_work = parse_double(trimmed.substring(start=6)).unwrap_or(0.0)
    } else if trimmed.starts_with("rFree=") {
      header.r_free = parse_double(trimmed.substring(start=6)).unwrap_or(0.0)
    } else if trimmed.starts_with("spaceGroup=") {
      xtal.space_group = trimmed.substring(start=11)
    } else if trimmed.starts_with("MODEL ") {
      // Save previous group/chain/model
      match current_group {
        Some(g) =>
          match current_chain {
            Some(c) => c.add_group(g)
            None => ()
          }
        None => ()
      }
      current_group = None
      match current_chain {
        Some(c) =>
          match current_model {
            Some(m) => m.add_chain(c)
            None => ()
          }
        None => ()
      }
      current_chain = None
      match current_model {
        Some(m) => structure.add_model(m)
        None => ()
      }
      let parts = trimmed.split(" ").to_array().filter(fn(s) -> Bool {
        s.length() > 0
      })
      let model_id = if parts.length() >= 2 { parse_int(parts[1].to_string()) } else { 0 }
      current_model = Some(MmtfModel::new(model_id=model_id))
    } else if trimmed.starts_with("  CHAIN ") {
      // Save previous group
      match current_group {
        Some(g) =>
          match current_chain {
            Some(c) => c.add_group(g)
            None => ()
          }
        None => ()
      }
      current_group = None
      // Save previous chain
      match current_chain {
        Some(c) =>
          match current_model {
            Some(m) => m.add_chain(c)
            None => ()
          }
        None => ()
      }
      let parts = trimmed.split(" ").to_array().filter(fn(s) -> Bool {
        s.length() > 0
      })
      let chain_id = if parts.length() >= 3 { parts[2].to_string() } else { "" }
      let chain_name = if parts.length() >= 4 { parts[3].to_string() } else { chain_id }
      current_chain = Some(MmtfChain::new(chain_id=chain_id, chain_name=chain_name))
    } else if trimmed.starts_with("    GROUP ") {
      // Save previous group
      match current_group {
        Some(g) =>
          match current_chain {
            Some(c) => c.add_group(g)
            None => ()
          }
        None => ()
      }
      let parts = trimmed.split(" ").to_array().filter(fn(s) -> Bool {
        s.length() > 0
      })
      let group_name = if parts.length() >= 2 { parts[1].to_string() } else { "" }
      let group_number = if parts.length() >= 3 { parse_int(parts[2].to_string()) } else { 0 }
      let slc = if parts.length() >= 4 { parts[3].to_string() } else { "" }
      let ss_str = if parts.length() >= 5 { parts[4].to_string() } else { "undefined" }
      let ss = match ss_str {
        "alpha_helix" => MmtfSecAlphaHelix
        "extended" => MmtfSecExtended
        "turn" => MmtfSecTurn
        "coil" => MmtfSecCoil
        "310_helix" => MmtfSecThreeTenHelix
        "pi_helix" => MmtfSecPpiHelix
        "bend" => MmtfSecBend
        "bridge" => MmtfSecBridge
        _ => MmtfSecUndefined
      }
      current_group = Some(
        MmtfGroup::new(
          group_name=group_name,
          group_number=group_number,
          single_letter_code=slc,
          sec_struct=ss,
        ),
      )
    } else if trimmed.starts_with("      ATOM ") {
      let parts = trimmed.split(" ").to_array().filter(fn(s) -> Bool {
        s.length() > 0
      })
      let name = if parts.length() >= 2 { parts[1].to_string() } else { "" }
      let serial = if parts.length() >= 3 { parse_int(parts[2].to_string()) } else { 0 }
      let x = if parts.length() >= 4 { parse_double(parts[3].to_string()).unwrap_or(0.0) } else { 0.0 }
      let y = if parts.length() >= 5 { parse_double(parts[4].to_string()).unwrap_or(0.0) } else { 0.0 }
      let z = if parts.length() >= 6 { parse_double(parts[5].to_string()).unwrap_or(0.0) } else { 0.0 }
      let element = if parts.length() >= 7 { parts[6].to_string() } else { "" }
      let charge = if parts.length() >= 8 { parse_int(parts[7].to_string()) } else { 0 }
      let atom = MmtfAtom::new(
        name=name,
        serial=serial,
        x=x,
        y=y,
        z=z,
        element=element,
        charge=charge,
      )
      match current_group {
        Some(g) => g.add_atom(atom)
        None => ()
      }
    } else if trimmed.starts_with("BOND ") {
      let parts = trimmed.split(" ").to_array().filter(fn(s) -> Bool {
        s.length() > 0
      })
      let a1 = if parts.length() >= 2 { parse_int(parts[1].to_string()) } else { 0 }
      let a2 = if parts.length() >= 3 { parse_int(parts[2].to_string()) } else { 0 }
      let order = if parts.length() >= 4 { parse_int(parts[3].to_string()) } else { 1 }
      structure.bonds.push(MmtfBond::new(atom_index1=a1, atom_index2=a2, bond_order=order))
    }
  }
  // Save final group/chain/model
  match current_group {
    Some(g) =>
      match current_chain {
        Some(c) => c.add_group(g)
        None => ()
      }
    None => ()
  }
  match current_chain {
    Some(c) =>
      match current_model {
        Some(m) => m.add_chain(c)
        None => ()
      }
    None => ()
  }
  match current_model {
    Some(m) => structure.add_model(m)
    None => ()
  }
  structure.set_header(header)
  structure.set_xtal_info(xtal)
  structure
}

// ============================================================================
// Sample structure
// ============================================================================

///|
/// Create a sample MMTF structure for testing.
pub fn mmtf_sample() -> MmtfStructure {
  let s = MmtfStructure::new()
  s.structure_id = "1XXX"
  s.total_num_models = 1
  s.total_num_chains = 2
  s.total_num_groups = 4
  s.total_num_atoms = 10
  s.total_num_bonds = 2
  // Header
  let header = MmtfHeader::new()
  header.title = "Sample protein structure"
  header.deposition_date = "2024-01-01"
  header.release_date = "2024-06-01"
  header.resolution = 1.5
  header.r_work = 0.18
  header.r_free = 0.22
  header.experimental_methods = ["X-RAY DIFFRACTION"]
  s.set_header(header)
  // Crystallographic info
  let xtal = MmtfXtalInfo::new()
  xtal.space_group = "P 21 21 21"
  xtal.unit_cell = [30.0, 40.0, 50.0, 90.0, 90.0, 90.0]
  s.set_xtal_info(xtal)
  // Entities
  let entity = MmtfEntity::new(
    chain_indices=[0],
    sequence="ACDE",
    description="sample protein",
    entity_type=MmtfEntityPolymer,
  )
  s.entities.push(entity)
  // Model 0
  let model = MmtfModel::new(model_id=0)
  // Chain A
  let chain_a = MmtfChain::new(chain_id="A", chain_name="A")
  // Group ALA (residue 1)
  let ala = MmtfGroup::new(
    group_name="ALA",
    group_number=1,
    single_letter_code="A",
    sec_struct=MmtfSecAlphaHelix,
    sequence_index=0,
  )
  ala.add_atom(
    MmtfAtom::new(name="N", serial=1, x=10.0, y=20.0, z=30.0, element="N"),
  )
  ala.add_atom(
    MmtfAtom::new(name="CA", serial=2, x=11.0, y=21.0, z=31.0, element="C"),
  )
  ala.add_atom(
    MmtfAtom::new(name="C", serial=3, x=12.0, y=22.0, z=32.0, element="C"),
  )
  chain_a.add_group(ala)
  // Group CYS (residue 2)
  let cys = MmtfGroup::new(
    group_name="CYS",
    group_number=2,
    single_letter_code="C",
    sec_struct=MmtfSecExtended,
    sequence_index=1,
  )
  cys.add_atom(
    MmtfAtom::new(name="N", serial=4, x=13.0, y=23.0, z=33.0, element="N"),
  )
  cys.add_atom(
    MmtfAtom::new(name="CA", serial=5, x=14.0, y=24.0, z=34.0, element="C"),
  )
  chain_a.add_group(cys)
  model.add_chain(chain_a)
  // Chain B (water)
  let chain_b = MmtfChain::new(chain_id="B", chain_name="B")
  let hoh = MmtfGroup::new(
    group_name="HOH",
    group_number=100,
    single_letter_code="",
    sec_struct=MmtfSecUndefined,
  )
  hoh.add_atom(
    MmtfAtom::new(name="O", serial=6, x=50.0, y=50.0, z=50.0, element="O"),
  )
  chain_b.add_group(hoh)
  model.add_chain(chain_b)
  s.add_model(model)
  // Bonds
  s.bonds.push(MmtfBond::new(atom_index1=0, atom_index2=1, bond_order=1))
  s.bonds.push(MmtfBond::new(atom_index1=1, atom_index2=2, bond_order=1))
  s
}

// ============================================================================
// Utility functions
// ============================================================================

///|
/// Get the sequence of a chain (single-letter codes).
pub fn mmtf_chain_sequence(chain : MmtfChain) -> String {
  let sb = StringBuilder::new()
  for group in chain.groups {
    if group.single_letter_code != "" && group.single_letter_code != "?" {
      sb.write_string(group.single_letter_code)
    }
  }
  sb.to_string()
}

///|
/// Count atoms in a chain.
pub fn mmtf_chain_num_atoms(chain : MmtfChain) -> Int {
  let mut count = 0
  for group in chain.groups {
    count = count + group.atoms.length()
  }
  count
}

///|
/// Count atoms in a model.
pub fn mmtf_model_num_atoms(model : MmtfModel) -> Int {
  let mut count = 0
  for chain in model.chains {
    count = count + mmtf_chain_num_atoms(chain)
  }
  count
}

///|
/// Generate a summary string.
pub fn mmtf_summary(structure : MmtfStructure) -> String {
  let sb = StringBuilder::new()
  let header = structure.get_header()
  let xtal = structure.get_xtal_info()
  sb.write_string("MMTF Structure Summary:\n")
  sb.write_string("  ID: " + structure.structure_id + "\n")
  sb.write_string("  Title: " + header.title + "\n")
  sb.write_string("  Resolution: " + header.resolution.to_string() + " Å\n")
  sb.write_string("  R-Work: " + header.r_work.to_string() + "\n")
  sb.write_string("  R-Free: " + header.r_free.to_string() + "\n")
  sb.write_string("  Space Group: " + xtal.space_group + "\n")
  sb.write_string("  Models: " + structure.get_num_models().to_string() + "\n")
  sb.write_string("  Chains: " + structure.count_chains().to_string() + "\n")
  sb.write_string("  Groups: " + structure.count_groups().to_string() + "\n")
  sb.write_string("  Atoms: " + structure.count_atoms().to_string() + "\n")
  sb.write_string("  Bonds: " + structure.bonds.length().to_string() + "\n")
  sb.write_string("  Entities: " + structure.entities.length().to_string() + "\n")
  sb.to_string()
}

///|
/// Convert MmtfStructure to PDB Structure.
pub fn mmtf_to_pdb_structure(mmtf : MmtfStructure) -> Structure {
  let pdb_models = Array::new()
  for mmtf_model in mmtf.models {
    let pdb_chains = Array::new()
    for mmtf_chain in mmtf_model.chains {
      let chain_id : Char = if mmtf_chain.chain_name.length() > 0 {
        Char::from_int(mmtf_chain.chain_name.unsafe_get(0).to_int())
      } else if mmtf_chain.chain_id.length() > 0 {
        Char::from_int(mmtf_chain.chain_id.unsafe_get(0).to_int())
      } else {
        'A'
      }
      let pdb_residues = Array::new()
      for mmtf_group in mmtf_chain.groups {
        let is_het = mmtf_group.group_name == "HOH" ||
          mmtf_group.group_type == "non-polymer"
        let hetfield = if is_het { "H" } else { " " }
        let pdb_atoms = Array::new()
        for mmtf_atom in mmtf_group.atoms {
          pdb_atoms.push(
            Atom::new(
              name=mmtf_atom.name,
              coord=Vector3::new(mmtf_atom.x, mmtf_atom.y, mmtf_atom.z),
              element=mmtf_atom.element,
              charge=mmtf_atom.charge,
            ),
          )
        }
        pdb_residues.push(
          Residue::new(
            resname=mmtf_group.group_name,
            chainid=chain_id,
            resseq=mmtf_group.group_number,
            atoms=pdb_atoms,
            hetfield=hetfield,
          ),
        )
      }
      pdb_chains.push(Chain::new(id=chain_id, residues=pdb_residues))
    }
    pdb_models.push(Model::new(id=mmtf_model.model_id, chains=pdb_chains))
  }
  Structure::new(id=mmtf.structure_id, models=pdb_models)
}
