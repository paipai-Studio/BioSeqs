// Bio.PDB.binary_cif - BinaryCIF decoding and PDB structure construction.
//
// BinaryCIF stores CIF data blocks in MessagePack and applies a reversible
// encoding pipeline to each column. The encoding list is decoded in reverse.

///|
pub suberror BinaryCifError {
  BinaryCifError(String)
}

///|
pub enum BinaryCifMask {
  BinaryCifPresent
  BinaryCifNotPresent
  BinaryCifUnknown
} derive(Eq, Debug)
pub extend BinaryCifMask with Eq::{equal, not_equal}
pub extend BinaryCifMask with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub enum BinaryCifColumnKind {
  BinaryCifInteger
  BinaryCifFloat
  BinaryCifText
} derive(Eq, Debug)
pub extend BinaryCifColumnKind with Eq::{equal, not_equal}
pub extend BinaryCifColumnKind with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub enum BinaryCifColumnData {
  BinaryCifIntegers(Array[Int])
  BinaryCifFloats(Array[Double])
  BinaryCifStrings(Array[String])
}

///|
pub struct BinaryCifColumn {
  name : String
  data : BinaryCifColumnData
  mask : Array[Int]
}

///|
pub struct BinaryCifCategory {
  name : String
  row_count : Int
  columns : Array[BinaryCifColumn]
}

///|
pub struct BinaryCifDataBlock {
  header : String
  categories : Array[BinaryCifCategory]
}

///|
pub struct BinaryCifFile {
  version : String
  encoder : String
  data_blocks : Array[BinaryCifDataBlock]
}

///|
pub fn BinaryCifColumn::kind(self : BinaryCifColumn) -> BinaryCifColumnKind {
  match self.data {
    BinaryCifIntegers(_) => BinaryCifInteger
    BinaryCifFloats(_) => BinaryCifFloat
    BinaryCifStrings(_) => BinaryCifText
  }
}

///|
pub fn BinaryCifColumn::length(self : BinaryCifColumn) -> Int {
  match self.data {
    BinaryCifIntegers(values) => values.length()
    BinaryCifFloats(values) => values.length()
    BinaryCifStrings(values) => values.length()
  }
}

///|
pub fn BinaryCifColumn::mask_at(
  self : BinaryCifColumn,
  index : Int,
) -> BinaryCifMask raise BinaryCifError {
  if index < 0 || index >= self.length() {
    raise BinaryCifError("BinaryCIF column index is out of bounds")
  }
  let value = if self.mask.length() == 0 { 0 } else { self.mask[index] }
  match value {
    0 => BinaryCifPresent
    1 => BinaryCifNotPresent
    2 => BinaryCifUnknown
    _ => raise BinaryCifError("BinaryCIF mask value must be 0, 1, or 2")
  }
}

///|
pub fn BinaryCifColumn::int_at(
  self : BinaryCifColumn,
  index : Int,
) -> Int? raise BinaryCifError {
  match self.mask_at(index) {
    BinaryCifPresent =>
      match self.data {
        BinaryCifIntegers(values) => Some(values[index])
        _ => raise BinaryCifError("BinaryCIF column is not an integer column")
      }
    _ => None
  }
}

///|
pub fn BinaryCifColumn::double_at(
  self : BinaryCifColumn,
  index : Int,
) -> Double? raise BinaryCifError {
  match self.mask_at(index) {
    BinaryCifPresent =>
      match self.data {
        BinaryCifIntegers(values) => Some(values[index].to_double())
        BinaryCifFloats(values) => Some(values[index])
        _ => raise BinaryCifError("BinaryCIF column is not a numeric column")
      }
    _ => None
  }
}

///|
pub fn BinaryCifColumn::string_at(
  self : BinaryCifColumn,
  index : Int,
) -> String? raise BinaryCifError {
  match self.mask_at(index) {
    BinaryCifPresent =>
      match self.data {
        BinaryCifStrings(values) => Some(values[index])
        _ => raise BinaryCifError("BinaryCIF column is not a string column")
      }
    _ => None
  }
}

///|
pub fn BinaryCifColumn::cif_token_at(
  self : BinaryCifColumn,
  index : Int,
) -> String raise BinaryCifError {
  match self.mask_at(index) {
    BinaryCifNotPresent => "."
    BinaryCifUnknown => "?"
    BinaryCifPresent =>
      match self.data {
        BinaryCifIntegers(values) => values[index].to_string()
        BinaryCifFloats(values) => values[index].to_string()
        BinaryCifStrings(values) => values[index]
      }
  }
}

///|
pub fn BinaryCifCategory::get_column(
  self : BinaryCifCategory,
  name : String,
) -> BinaryCifColumn? {
  for column in self.columns {
    if column.name == name || self.name + "." + column.name == name {
      return Some(column)
    }
  }
  None
}

///|
pub fn BinaryCifDataBlock::get_category(
  self : BinaryCifDataBlock,
  name : String,
) -> BinaryCifCategory? {
  for category in self.categories {
    if category.name == name ||
      (
        name.length() > 0 &&
        name.unsafe_get(0).to_int() != '_'.to_int() &&
        category.name == "_" + name
      ) {
      return Some(category)
    }
  }
  None
}

///|
pub fn BinaryCifFile::get_block(
  self : BinaryCifFile,
  index : Int,
) -> BinaryCifDataBlock? {
  if index < 0 || index >= self.data_blocks.length() {
    None
  } else {
    Some(self.data_blocks[index])
  }
}

///|
pub fn BinaryCifFile::summary(self : BinaryCifFile) -> String {
  let mut category_count = 0
  let mut row_count = 0
  for block in self.data_blocks {
    category_count = category_count + block.categories.length()
    for category in block.categories {
      row_count = row_count + category.row_count
    }
  }
  "BinaryCIF(version=" +
  self.version +
  ", blocks=" +
  self.data_blocks.length().to_string() +
  ", categories=" +
  category_count.to_string() +
  ", rows=" +
  row_count.to_string() +
  ")"
}

// MessagePack model and decoder.

///|
priv enum BcifMessage {
  BcifNil
  BcifBoolean(Bool)
  BcifIntegerValue(Int)
  BcifFloatValue(Double)
  BcifTextValue(String)
  BcifBinaryValue(Array[Int])
  BcifArrayValue(Array[BcifMessage])
  BcifMapValue(Map[String, BcifMessage])
}

///|
priv struct BcifReader {
  data : Array[Int]
  mut position : Int
}

///|
fn bcif_read_byte(reader : BcifReader) -> Int raise BinaryCifError {
  if reader.position >= reader.data.length() {
    raise BinaryCifError("Truncated MessagePack input")
  }
  let value = reader.data[reader.position]
  reader.position = reader.position + 1
  value
}

///|
fn bcif_read_unsigned(
  reader : BcifReader,
  byte_count : Int,
) -> Int raise BinaryCifError {
  if byte_count <= 0 || byte_count > 4 {
    raise BinaryCifError("Unsupported MessagePack unsigned integer width")
  }
  if reader.position + byte_count > reader.data.length() {
    raise BinaryCifError("Truncated MessagePack integer")
  }
  let first = reader.data[reader.position]
  if byte_count == 4 && first >= 128 {
    raise BinaryCifError("MessagePack unsigned integer exceeds MoonBit Int")
  }
  let mut value = 0
  for _ in 0..<byte_count {
    value = value * 256 + bcif_read_byte(reader)
  }
  value
}

///|
fn bcif_read_signed(
  reader : BcifReader,
  byte_count : Int,
) -> Int raise BinaryCifError {
  if byte_count <= 0 || byte_count > 4 {
    raise BinaryCifError("Unsupported MessagePack signed integer width")
  }
  if reader.position + byte_count > reader.data.length() {
    raise BinaryCifError("Truncated MessagePack integer")
  }
  let first = bcif_read_byte(reader)
  let mut value = if first >= 128 { first - 256 } else { first }
  for _ in 1..<byte_count {
    value = value * 256 + bcif_read_byte(reader)
  }
  value
}

///|
fn bcif_read_u64(reader : BcifReader) -> Int raise BinaryCifError {
  for _ in 0..<4 {
    if bcif_read_byte(reader) != 0 {
      raise BinaryCifError("MessagePack uint64 exceeds MoonBit Int")
    }
  }
  bcif_read_unsigned(reader, 4)
}

///|
fn bcif_read_i64(reader : BcifReader) -> Int raise BinaryCifError {
  if reader.position + 8 > reader.data.length() {
    raise BinaryCifError("Truncated MessagePack int64")
  }
  let sign_byte = reader.data[reader.position]
  let expected = if sign_byte >= 128 { 255 } else { 0 }
  for _ in 0..<4 {
    if bcif_read_byte(reader) != expected {
      raise BinaryCifError("MessagePack int64 exceeds MoonBit Int")
    }
  }
  let next = reader.data[reader.position]
  if (expected == 0 && next >= 128) || (expected == 255 && next < 128) {
    raise BinaryCifError("MessagePack int64 exceeds MoonBit Int")
  }
  bcif_read_signed(reader, 4)
}

///|
fn bcif_pow2(exponent : Int) -> Double {
  let mut power = if exponent < 0 { -exponent } else { exponent }
  let mut base = 2.0
  let mut result = 1.0
  while power > 0 {
    if power % 2 == 1 {
      result = result * base
    }
    base = base * base
    power = power / 2
  }
  if exponent < 0 {
    1.0 / result
  } else {
    result
  }
}

///|
fn bcif_float32_from_bytes(
  b0 : Int,
  b1 : Int,
  b2 : Int,
  b3 : Int,
) -> Double raise BinaryCifError {
  let sign = if b0 >= 128 { -1.0 } else { 1.0 }
  let exponent = (b0 & 0x7F) * 2 + (b1 >> 7)
  let mantissa = (b1 & 0x7F) * 65536 + b2 * 256 + b3
  if exponent == 255 {
    raise BinaryCifError("BinaryCIF does not accept non-finite floats")
  }
  if exponent == 0 {
    sign * mantissa.to_double() * bcif_pow2(-149)
  } else {
    sign * (1.0 + mantissa.to_double() / 8388608.0) * bcif_pow2(exponent - 127)
  }
}

///|
fn bcif_float64_from_bytes(
  bytes : Array[Int],
  offset : Int,
  little_endian : Bool,
) -> Double raise BinaryCifError {
  let ordered = Array::make(8, 0)
  for i in 0..<8 {
    ordered[i] = if little_endian {
      bytes[offset + 7 - i]
    } else {
      bytes[offset + i]
    }
  }
  let sign = if ordered[0] >= 128 { -1.0 } else { 1.0 }
  let exponent = (ordered[0] & 0x7F) * 16 + (ordered[1] >> 4)
  if exponent == 2047 {
    raise BinaryCifError("BinaryCIF does not accept non-finite doubles")
  }
  let mut fraction = (ordered[1] & 0x0F).to_double() / 16.0
  let mut scale = 1.0 / 4096.0
  for i in 2..<8 {
    fraction = fraction + ordered[i].to_double() * scale
    scale = scale / 256.0
  }
  if exponent == 0 {
    sign * fraction * bcif_pow2(-1022)
  } else {
    sign * (1.0 + fraction) * bcif_pow2(exponent - 1023)
  }
}

///|
fn bcif_read_float32(reader : BcifReader) -> Double raise BinaryCifError {
  let b0 = bcif_read_byte(reader)
  let b1 = bcif_read_byte(reader)
  let b2 = bcif_read_byte(reader)
  let b3 = bcif_read_byte(reader)
  bcif_float32_from_bytes(b0, b1, b2, b3)
}

///|
fn bcif_read_float64(reader : BcifReader) -> Double raise BinaryCifError {
  if reader.position + 8 > reader.data.length() {
    raise BinaryCifError("Truncated MessagePack float64")
  }
  let value = bcif_float64_from_bytes(reader.data, reader.position, false)
  reader.position = reader.position + 8
  value
}

///|
fn bcif_read_utf8(
  reader : BcifReader,
  length : Int,
) -> String raise BinaryCifError {
  if length < 0 || reader.position + length > reader.data.length() {
    raise BinaryCifError("Truncated MessagePack string")
  }
  let end = reader.position + length
  let output = StringBuilder()
  while reader.position < end {
    let first = bcif_read_byte(reader)
    if first < 0x80 {
      output.write_char(first.unsafe_to_char())
    } else if first >= 0xC2 && first <= 0xDF {
      if reader.position >= end {
        raise BinaryCifError("Truncated UTF-8 sequence")
      }
      let second = bcif_read_byte(reader)
      if second < 0x80 || second > 0xBF {
        raise BinaryCifError("Invalid UTF-8 continuation byte")
      }
      output.write_char(
        ((first & 0x1F) * 64 + (second & 0x3F)).unsafe_to_char(),
      )
    } else if first >= 0xE0 && first <= 0xEF {
      if reader.position + 2 > end {
        raise BinaryCifError("Truncated UTF-8 sequence")
      }
      let second = bcif_read_byte(reader)
      let third = bcif_read_byte(reader)
      if second < 0x80 ||
        second > 0xBF ||
        third < 0x80 ||
        third > 0xBF ||
        (first == 0xE0 && second < 0xA0) ||
        (first == 0xED && second >= 0xA0) {
        raise BinaryCifError("Invalid UTF-8 sequence")
      }
      output.write_char(
        ((first & 0x0F) * 4096 + (second & 0x3F) * 64 + (third & 0x3F)).unsafe_to_char(),
      )
    } else if first >= 0xF0 && first <= 0xF4 {
      if reader.position + 3 > end {
        raise BinaryCifError("Truncated UTF-8 sequence")
      }
      let second = bcif_read_byte(reader)
      let third = bcif_read_byte(reader)
      let fourth = bcif_read_byte(reader)
      if second < 0x80 ||
        second > 0xBF ||
        third < 0x80 ||
        third > 0xBF ||
        fourth < 0x80 ||
        fourth > 0xBF ||
        (first == 0xF0 && second < 0x90) ||
        (first == 0xF4 && second >= 0x90) {
        raise BinaryCifError("Invalid UTF-8 sequence")
      }
      output.write_char(
        ((first & 0x07) * 262144 +
        (second & 0x3F) * 4096 +
        (third & 0x3F) * 64 +
        (fourth & 0x3F)).unsafe_to_char(),
      )
    } else {
      raise BinaryCifError("Invalid UTF-8 leading byte")
    }
  }
  output.to_string()
}

///|
fn bcif_read_binary(
  reader : BcifReader,
  length : Int,
) -> Array[Int] raise BinaryCifError {
  if length < 0 || reader.position + length > reader.data.length() {
    raise BinaryCifError("Truncated MessagePack binary value")
  }
  let result = Array::make(length, 0)
  for i in 0..<length {
    result[i] = bcif_read_byte(reader)
  }
  result
}

///|
fn bcif_parse_message(
  reader : BcifReader,
  depth : Int,
) -> BcifMessage raise BinaryCifError {
  if depth > 128 {
    raise BinaryCifError("MessagePack nesting is too deep")
  }
  let marker = bcif_read_byte(reader)
  if marker <= 0x7F {
    return BcifIntegerValue(marker)
  }
  if marker >= 0xE0 {
    return BcifIntegerValue(marker - 256)
  }
  if marker >= 0xA0 && marker <= 0xBF {
    return BcifTextValue(bcif_read_utf8(reader, marker & 0x1F))
  }
  if marker >= 0x90 && marker <= 0x9F {
    let length = marker & 0x0F
    let values = Array(capacity=length)
    for _ in 0..<length {
      values.push(bcif_parse_message(reader, depth + 1))
    }
    return BcifArrayValue(values)
  }
  if marker >= 0x80 && marker <= 0x8F {
    return bcif_parse_map(reader, marker & 0x0F, depth + 1)
  }
  match marker {
    0xC0 => BcifNil
    0xC2 => BcifBoolean(false)
    0xC3 => BcifBoolean(true)
    0xC4 =>
      BcifBinaryValue(bcif_read_binary(reader, bcif_read_unsigned(reader, 1)))
    0xC5 =>
      BcifBinaryValue(bcif_read_binary(reader, bcif_read_unsigned(reader, 2)))
    0xC6 =>
      BcifBinaryValue(bcif_read_binary(reader, bcif_read_unsigned(reader, 4)))
    0xCA => BcifFloatValue(bcif_read_float32(reader))
    0xCB => BcifFloatValue(bcif_read_float64(reader))
    0xCC => BcifIntegerValue(bcif_read_unsigned(reader, 1))
    0xCD => BcifIntegerValue(bcif_read_unsigned(reader, 2))
    0xCE => BcifIntegerValue(bcif_read_unsigned(reader, 4))
    0xCF => BcifIntegerValue(bcif_read_u64(reader))
    0xD0 => BcifIntegerValue(bcif_read_signed(reader, 1))
    0xD1 => BcifIntegerValue(bcif_read_signed(reader, 2))
    0xD2 => BcifIntegerValue(bcif_read_signed(reader, 4))
    0xD3 => BcifIntegerValue(bcif_read_i64(reader))
    0xD9 => BcifTextValue(bcif_read_utf8(reader, bcif_read_unsigned(reader, 1)))
    0xDA => BcifTextValue(bcif_read_utf8(reader, bcif_read_unsigned(reader, 2)))
    0xDB => BcifTextValue(bcif_read_utf8(reader, bcif_read_unsigned(reader, 4)))
    0xDC => {
      let length = bcif_read_unsigned(reader, 2)
      let values = Array(capacity=length)
      for _ in 0..<length {
        values.push(bcif_parse_message(reader, depth + 1))
      }
      BcifArrayValue(values)
    }
    0xDD => {
      let length = bcif_read_unsigned(reader, 4)
      let values = Array(capacity=length)
      for _ in 0..<length {
        values.push(bcif_parse_message(reader, depth + 1))
      }
      BcifArrayValue(values)
    }
    0xDE => bcif_parse_map(reader, bcif_read_unsigned(reader, 2), depth + 1)
    0xDF => bcif_parse_map(reader, bcif_read_unsigned(reader, 4), depth + 1)
    _ =>
      raise BinaryCifError(
        "Unsupported MessagePack marker " + marker.to_string(),
      )
  }
}

///|
fn bcif_parse_map(
  reader : BcifReader,
  length : Int,
  depth : Int,
) -> BcifMessage raise BinaryCifError {
  let values : Map[String, BcifMessage] = Map([], capacity=length)
  for _ in 0..<length {
    let key = match bcif_parse_message(reader, depth) {
      BcifTextValue(value) => value
      _ => raise BinaryCifError("MessagePack map keys must be strings")
    }
    values[key] = bcif_parse_message(reader, depth)
  }
  BcifMapValue(values)
}

///|
fn bcif_as_map(
  value : BcifMessage,
  context : String,
) -> Map[String, BcifMessage] raise BinaryCifError {
  match value {
    BcifMapValue(result) => result
    _ => raise BinaryCifError(context + " must be a MessagePack map")
  }
}

///|
fn bcif_as_array(
  value : BcifMessage,
  context : String,
) -> Array[BcifMessage] raise BinaryCifError {
  match value {
    BcifArrayValue(result) => result
    _ => raise BinaryCifError(context + " must be a MessagePack array")
  }
}

///|
fn bcif_as_string(
  value : BcifMessage,
  context : String,
) -> String raise BinaryCifError {
  match value {
    BcifTextValue(result) => result
    _ => raise BinaryCifError(context + " must be a string")
  }
}

///|
fn bcif_as_int(
  value : BcifMessage,
  context : String,
) -> Int raise BinaryCifError {
  match value {
    BcifIntegerValue(result) => result
    _ => raise BinaryCifError(context + " must be an integer")
  }
}

///|
fn bcif_as_double(
  value : BcifMessage,
  context : String,
) -> Double raise BinaryCifError {
  match value {
    BcifIntegerValue(result) => result.to_double()
    BcifFloatValue(result) => result
    _ => raise BinaryCifError(context + " must be numeric")
  }
}

///|
fn bcif_as_bool(
  value : BcifMessage,
  context : String,
) -> Bool raise BinaryCifError {
  match value {
    BcifBoolean(result) => result
    _ => raise BinaryCifError(context + " must be a boolean")
  }
}

///|
fn bcif_as_binary(
  value : BcifMessage,
  context : String,
) -> Array[Int] raise BinaryCifError {
  match value {
    BcifBinaryValue(result) => result
    _ => raise BinaryCifError(context + " must be binary data")
  }
}

///|
fn bcif_required(
  values : Map[String, BcifMessage],
  key : String,
  context : String,
) -> BcifMessage raise BinaryCifError {
  match values.get(key) {
    Some(value) => value
    None => raise BinaryCifError(context + " is missing '" + key + "'")
  }
}

///|
fn bcif_encoding_int(
  values : Map[String, BcifMessage],
  key : String,
  context : String,
) -> Int raise BinaryCifError {
  bcif_as_int(bcif_required(values, key, context), context + "." + key)
}

///|
fn bcif_encoding_double(
  values : Map[String, BcifMessage],
  key : String,
  context : String,
) -> Double raise BinaryCifError {
  bcif_as_double(bcif_required(values, key, context), context + "." + key)
}

// BinaryCIF encoding decoders.

///|
fn bcif_read_i16_le(data : Array[Int], offset : Int) -> Int {
  let high = data[offset + 1]
  (if high >= 128 { high - 256 } else { high }) * 256 + data[offset]
}

///|
fn bcif_read_u16_le(data : Array[Int], offset : Int) -> Int {
  data[offset] + data[offset + 1] * 256
}

///|
fn bcif_read_i32_le(data : Array[Int], offset : Int) -> Int {
  let high = data[offset + 3]
  (if high >= 128 { high - 256 } else { high }) * 16777216 +
  data[offset + 2] * 65536 +
  data[offset + 1] * 256 +
  data[offset]
}

///|
fn bcif_read_u32_le(
  data : Array[Int],
  offset : Int,
) -> Int raise BinaryCifError {
  if data[offset + 3] >= 128 {
    raise BinaryCifError("BinaryCIF UInt32 exceeds MoonBit Int")
  }
  data[offset + 3] * 16777216 +
  data[offset + 2] * 65536 +
  data[offset + 1] * 256 +
  data[offset]
}

///|
fn bcif_validate_bytes(data : Array[Int]) -> Unit raise BinaryCifError {
  for value in data {
    if value < 0 || value > 255 {
      raise BinaryCifError("BinaryCIF byte array contains a non-byte value")
    }
  }
}

///|
pub fn binary_cif_decode_int_bytes(
  data : Array[Int],
  type_code : Int,
) -> Array[Int] raise BinaryCifError {
  bcif_validate_bytes(data)
  let width = match type_code {
    1 | 4 => 1
    2 | 5 => 2
    3 | 6 => 4
    _ => raise BinaryCifError("BinaryCIF type is not an integer byte array")
  }
  if data.length() % width != 0 {
    raise BinaryCifError("BinaryCIF byte array has an invalid length")
  }
  let result = Array(capacity=data.length() / width)
  let mut offset = 0
  while offset < data.length() {
    let value = match type_code {
      1 => {
        let raw = data[offset]
        if raw >= 128 {
          raw - 256
        } else {
          raw
        }
      }
      2 => bcif_read_i16_le(data, offset)
      3 => bcif_read_i32_le(data, offset)
      4 => data[offset]
      5 => bcif_read_u16_le(data, offset)
      6 => bcif_read_u32_le(data, offset)
      _ => 0
    }
    result.push(value)
    offset = offset + width
  }
  result
}

///|
pub fn binary_cif_decode_float_bytes(
  data : Array[Int],
  type_code : Int,
) -> Array[Double] raise BinaryCifError {
  bcif_validate_bytes(data)
  let width = match type_code {
    32 => 4
    33 => 8
    _ =>
      raise BinaryCifError("BinaryCIF type is not a floating-point byte array")
  }
  if data.length() % width != 0 {
    raise BinaryCifError(
      "BinaryCIF floating-point byte array has an invalid length",
    )
  }
  let result = Array(capacity=data.length() / width)
  let mut offset = 0
  while offset < data.length() {
    if type_code == 32 {
      result.push(
        bcif_float32_from_bytes(
          data[offset + 3],
          data[offset + 2],
          data[offset + 1],
          data[offset],
        ),
      )
    } else {
      result.push(bcif_float64_from_bytes(data, offset, true))
    }
    offset = offset + width
  }
  result
}

///|
pub fn binary_cif_decode_integer_packing(
  data : Array[Int],
  byte_count : Int,
  is_unsigned : Bool,
  source_size : Int,
) -> Array[Int] raise BinaryCifError {
  if byte_count != 1 && byte_count != 2 {
    raise BinaryCifError("BinaryCIF integer packing uses 1- or 2-byte values")
  }
  if source_size < 0 {
    raise BinaryCifError("BinaryCIF integer packing source size is negative")
  }
  let upper = if is_unsigned {
    if byte_count == 1 {
      255
    } else {
      65535
    }
  } else if byte_count == 1 {
    127
  } else {
    32767
  }
  let lower = if is_unsigned {
    0
  } else if byte_count == 1 {
    -128
  } else {
    -32768
  }
  let result = Array(capacity=source_size)
  let mut accumulator = 0
  for value in data {
    if value < lower || value > upper {
      raise BinaryCifError(
        "BinaryCIF integer packing value is outside its byte range",
      )
    }
    accumulator = accumulator + value
    let continuation = value == upper || (!is_unsigned && value == lower)
    if !continuation {
      result.push(accumulator)
      accumulator = 0
    }
  }
  if accumulator != 0 || result.length() != source_size {
    raise BinaryCifError("BinaryCIF integer packing source size does not match")
  }
  result
}

///|
pub fn binary_cif_decode_run_length(
  data : Array[Int],
  source_size : Int,
) -> Array[Int] raise BinaryCifError {
  if source_size < 0 || data.length() % 2 != 0 {
    raise BinaryCifError("Invalid BinaryCIF run-length encoding")
  }
  let result = Array(capacity=source_size)
  let mut index = 0
  while index < data.length() {
    let value = data[index]
    let count = data[index + 1]
    if count < 0 || result.length() + count > source_size {
      raise BinaryCifError("Invalid BinaryCIF run length")
    }
    for _ in 0..<count {
      result.push(value)
    }
    index = index + 2
  }
  if result.length() != source_size {
    raise BinaryCifError("BinaryCIF run-length source size does not match")
  }
  result
}

///|
pub fn binary_cif_decode_delta(data : Array[Int], origin : Int) -> Array[Int] {
  let result = Array(capacity=data.length())
  let mut current = origin
  for value in data {
    current = current + value
    result.push(current)
  }
  result
}

///|
pub fn binary_cif_decode_fixed_point(
  data : Array[Int],
  factor : Double,
) -> Array[Double] raise BinaryCifError {
  if factor == 0.0 || factor.abs() > 1.0e300 {
    raise BinaryCifError(
      "BinaryCIF fixed-point factor must be finite and non-zero",
    )
  }
  let result = Array(capacity=data.length())
  for value in data {
    result.push(value.to_double() / factor)
  }
  result
}

///|
pub fn binary_cif_decode_interval_quantization(
  data : Array[Int],
  minimum : Double,
  maximum : Double,
  number_of_steps : Int,
) -> Array[Double] raise BinaryCifError {
  if number_of_steps < 2 ||
    minimum.abs() > 1.0e300 ||
    maximum.abs() > 1.0e300 ||
    maximum < minimum {
    raise BinaryCifError("Invalid BinaryCIF interval quantization parameters")
  }
  let step = (maximum - minimum) / (number_of_steps - 1).to_double()
  let result = Array(capacity=data.length())
  for value in data {
    if value < 0 || value >= number_of_steps {
      raise BinaryCifError("BinaryCIF quantized value is outside its interval")
    }
    result.push(minimum + value.to_double() * step)
  }
  result
}

///|
pub fn binary_cif_decode_string_array(
  indices : Array[Int],
  offsets : Array[Int],
  string_data : String,
) -> Array[String] raise BinaryCifError {
  if offsets.length() == 0 || offsets[0] != 0 {
    raise BinaryCifError("BinaryCIF string offsets must start at zero")
  }
  for i in 1..<offsets.length() {
    if offsets[i] < offsets[i - 1] || offsets[i] > string_data.length() {
      raise BinaryCifError("BinaryCIF string offsets are invalid")
    }
  }
  if offsets[offsets.length() - 1] != string_data.length() {
    raise BinaryCifError(
      "BinaryCIF string offsets do not cover the string data",
    )
  }
  let dictionary = Array(capacity=offsets.length() - 1)
  for i in 0..<(offsets.length() - 1) {
    dictionary.push(string_data[offsets[i]:offsets[i + 1]].to_owned())
  }
  let result = Array(capacity=indices.length())
  for index in indices {
    if index < 0 || index >= dictionary.length() {
      raise BinaryCifError("BinaryCIF string lookup index is out of bounds")
    }
    result.push(dictionary[index])
  }
  result
}

///|
fn bcif_decoded_length(data : BinaryCifColumnData) -> Int {
  match data {
    BinaryCifIntegers(values) => values.length()
    BinaryCifFloats(values) => values.length()
    BinaryCifStrings(values) => values.length()
  }
}

///|
fn bcif_decode_raw(
  bytes : Array[Int],
  encodings : Array[BcifMessage],
) -> BinaryCifColumnData raise BinaryCifError {
  let mut raw_bytes : Array[Int]? = Some(bytes)
  let mut decoded : BinaryCifColumnData? = None
  let mut index = encodings.length() - 1
  while index >= 0 {
    let encoding = bcif_as_map(encodings[index], "BinaryCIF encoding")
    let kind = bcif_as_string(
      bcif_required(encoding, "kind", "BinaryCIF encoding"),
      "BinaryCIF encoding.kind",
    )
    match kind {
      "ByteArray" => {
        let source = match raw_bytes {
          Some(value) => value
          None =>
            raise BinaryCifError(
              "ByteArray must be the final BinaryCIF encoding",
            )
        }
        let type_code = bcif_encoding_int(encoding, "type", "ByteArray")
        decoded = Some(
          if type_code == 32 || type_code == 33 {
            BinaryCifFloats(binary_cif_decode_float_bytes(source, type_code))
          } else {
            BinaryCifIntegers(binary_cif_decode_int_bytes(source, type_code))
          },
        )
        raw_bytes = None
      }
      "IntegerPacking" => {
        let values = match decoded {
          Some(BinaryCifIntegers(value)) => value
          _ => raise BinaryCifError("IntegerPacking requires integer input")
        }
        let byte_count = bcif_encoding_int(
          encoding, "byteCount", "IntegerPacking",
        )
        let source_size = bcif_encoding_int(
          encoding, "srcSize", "IntegerPacking",
        )
        let is_unsigned = bcif_as_bool(
          bcif_required(encoding, "isUnsigned", "IntegerPacking"),
          "IntegerPacking.isUnsigned",
        )
        decoded = Some(
          BinaryCifIntegers(
            binary_cif_decode_integer_packing(
              values, byte_count, is_unsigned, source_size,
            ),
          ),
        )
      }
      "RunLength" => {
        let values = match decoded {
          Some(BinaryCifIntegers(value)) => value
          _ => raise BinaryCifError("RunLength requires integer input")
        }
        decoded = Some(
          BinaryCifIntegers(
            binary_cif_decode_run_length(
              values,
              bcif_encoding_int(encoding, "srcSize", "RunLength"),
            ),
          ),
        )
      }
      "Delta" => {
        let values = match decoded {
          Some(BinaryCifIntegers(value)) => value
          _ => raise BinaryCifError("Delta requires integer input")
        }
        decoded = Some(
          BinaryCifIntegers(
            binary_cif_decode_delta(
              values,
              bcif_encoding_int(encoding, "origin", "Delta"),
            ),
          ),
        )
      }
      "FixedPoint" => {
        let values = match decoded {
          Some(BinaryCifIntegers(value)) => value
          _ => raise BinaryCifError("FixedPoint requires integer input")
        }
        decoded = Some(
          BinaryCifFloats(
            binary_cif_decode_fixed_point(
              values,
              bcif_encoding_double(encoding, "factor", "FixedPoint"),
            ),
          ),
        )
      }
      "IntervalQuantization" => {
        let values = match decoded {
          Some(BinaryCifIntegers(value)) => value
          _ =>
            raise BinaryCifError("IntervalQuantization requires integer input")
        }
        let steps = match encoding.get("numSteps") {
          Some(value) => bcif_as_int(value, "IntervalQuantization.numSteps")
          None =>
            bcif_as_int(
              bcif_required(encoding, "num_steps", "IntervalQuantization"),
              "IntervalQuantization.num_steps",
            )
        }
        decoded = Some(
          BinaryCifFloats(
            binary_cif_decode_interval_quantization(
              values,
              bcif_encoding_double(encoding, "min", "IntervalQuantization"),
              bcif_encoding_double(encoding, "max", "IntervalQuantization"),
              steps,
            ),
          ),
        )
      }
      "StringArray" => {
        let source = match raw_bytes {
          Some(value) => value
          None => raise BinaryCifError("StringArray requires raw byte input")
        }
        let data_encodings = bcif_as_array(
          bcif_required(encoding, "dataEncoding", "StringArray"),
          "StringArray.dataEncoding",
        )
        let offset_encodings = bcif_as_array(
          bcif_required(encoding, "offsetEncoding", "StringArray"),
          "StringArray.offsetEncoding",
        )
        let offsets_bytes = bcif_as_binary(
          bcif_required(encoding, "offsets", "StringArray"),
          "StringArray.offsets",
        )
        let indices = match bcif_decode_raw(source, data_encodings) {
          BinaryCifIntegers(value) => value
          _ =>
            raise BinaryCifError(
              "StringArray lookup data must decode to integers",
            )
        }
        let offsets = match bcif_decode_raw(offsets_bytes, offset_encodings) {
          BinaryCifIntegers(value) => value
          _ =>
            raise BinaryCifError("StringArray offsets must decode to integers")
        }
        decoded = Some(
          BinaryCifStrings(
            binary_cif_decode_string_array(
              indices,
              offsets,
              bcif_as_string(
                bcif_required(encoding, "stringData", "StringArray"),
                "StringArray.stringData",
              ),
            ),
          ),
        )
        raw_bytes = None
      }
      _ => raise BinaryCifError("Unsupported BinaryCIF encoding '" + kind + "'")
    }
    index = index - 1
  }
  match decoded {
    Some(value) => value
    None => raise BinaryCifError("BinaryCIF data has no encoding")
  }
}

///|
fn bcif_decode_data(
  value : BcifMessage,
  context : String,
) -> BinaryCifColumnData raise BinaryCifError {
  let data = bcif_as_map(value, context)
  let bytes = bcif_as_binary(
    bcif_required(data, "data", context),
    context + ".data",
  )
  let encodings = bcif_as_array(
    bcif_required(data, "encoding", context),
    context + ".encoding",
  )
  bcif_decode_raw(bytes, encodings)
}

///|
fn bcif_parse_column(
  value : BcifMessage,
  row_count : Int,
) -> BinaryCifColumn raise BinaryCifError {
  let column = bcif_as_map(value, "BinaryCIF column")
  let name = bcif_as_string(
    bcif_required(column, "name", "BinaryCIF column"),
    "BinaryCIF column.name",
  )
  let data = bcif_decode_data(
    bcif_required(column, "data", "BinaryCIF column"),
    "BinaryCIF column '" + name + "' data",
  )
  if bcif_decoded_length(data) != row_count {
    raise BinaryCifError(
      "BinaryCIF column '" + name + "' length does not match rowCount",
    )
  }
  let mask = match column.get("mask") {
    None | Some(BcifNil) => Array::make(row_count, 0)
    Some(mask_value) =>
      match bcif_decode_data(mask_value, "BinaryCIF column mask") {
        BinaryCifIntegers(values) => {
          if values.length() != row_count {
            raise BinaryCifError(
              "BinaryCIF mask length does not match rowCount",
            )
          }
          for item in values {
            if item < 0 || item > 2 {
              raise BinaryCifError("BinaryCIF mask value must be 0, 1, or 2")
            }
          }
          values
        }
        _ => raise BinaryCifError("BinaryCIF mask must decode to integers")
      }
  }
  BinaryCifColumn::{ name, data, mask, }
}

///|
fn bcif_parse_category(
  value : BcifMessage,
) -> BinaryCifCategory raise BinaryCifError {
  let category = bcif_as_map(value, "BinaryCIF category")
  let name = bcif_as_string(
    bcif_required(category, "name", "BinaryCIF category"),
    "BinaryCIF category.name",
  )
  let row_count = bcif_as_int(
    bcif_required(category, "rowCount", "BinaryCIF category"),
    "BinaryCIF category.rowCount",
  )
  if row_count < 0 {
    raise BinaryCifError("BinaryCIF category rowCount is negative")
  }
  let raw_columns = bcif_as_array(
    bcif_required(category, "columns", "BinaryCIF category"),
    "BinaryCIF category.columns",
  )
  let columns = Array(capacity=raw_columns.length())
  let seen : Map[String, Bool] = Map([], capacity=raw_columns.length())
  for raw_column in raw_columns {
    let column = bcif_parse_column(raw_column, row_count)
    if seen.contains(column.name) {
      raise BinaryCifError("Duplicate BinaryCIF column '" + column.name + "'")
    }
    seen[column.name] = true
    columns.push(column)
  }
  BinaryCifCategory::{ name, row_count, columns, }
}

///|
fn bcif_parse_block(
  value : BcifMessage,
) -> BinaryCifDataBlock raise BinaryCifError {
  let block = bcif_as_map(value, "BinaryCIF data block")
  let header = bcif_as_string(
    bcif_required(block, "header", "BinaryCIF data block"),
    "BinaryCIF data block.header",
  )
  let raw_categories = bcif_as_array(
    bcif_required(block, "categories", "BinaryCIF data block"),
    "BinaryCIF data block.categories",
  )
  let categories = Array(capacity=raw_categories.length())
  let seen : Map[String, Bool] = Map([], capacity=raw_categories.length())
  for raw_category in raw_categories {
    let category = bcif_parse_category(raw_category)
    if seen.contains(category.name) {
      raise BinaryCifError(
        "Duplicate BinaryCIF category '" + category.name + "'",
      )
    }
    seen[category.name] = true
    categories.push(category)
  }
  BinaryCifDataBlock::{ header, categories, }
}

///|
pub fn binary_cif_parse(
  input : Array[Int],
) -> BinaryCifFile raise BinaryCifError {
  if input.length() == 0 {
    raise BinaryCifError("Empty BinaryCIF input")
  }
  for byte in input {
    if byte < 0 || byte > 255 {
      raise BinaryCifError("BinaryCIF input contains a non-byte value")
    }
  }
  if input.length() >= 2 && input[0] == 0x1F && input[1] == 0x8B {
    raise BinaryCifError(
      "Gzip-compressed BinaryCIF input must be decompressed first",
    )
  }
  let reader = BcifReader::{ data: input, position: 0, }
  let root = bcif_as_map(bcif_parse_message(reader, 0), "BinaryCIF document")
  if reader.position != input.length() {
    raise BinaryCifError("Additional bytes follow the BinaryCIF document")
  }
  let version = bcif_as_string(
    bcif_required(root, "version", "BinaryCIF document"),
    "BinaryCIF version",
  )
  let encoder = bcif_as_string(
    bcif_required(root, "encoder", "BinaryCIF document"),
    "BinaryCIF encoder",
  )
  let raw_blocks = bcif_as_array(
    bcif_required(root, "dataBlocks", "BinaryCIF document"),
    "BinaryCIF dataBlocks",
  )
  if raw_blocks.length() == 0 {
    raise BinaryCifError("BinaryCIF document has no data blocks")
  }
  let data_blocks = Array(capacity=raw_blocks.length())
  for raw_block in raw_blocks {
    data_blocks.push(bcif_parse_block(raw_block))
  }
  BinaryCifFile::{ version, encoder, data_blocks, }
}

// PDB Structure conversion.

///|
priv struct BcifResidueBuilder {
  resname : String
  chain_id : Char
  sequence_id : Int
  insertion_code : Char
  hetero_field : String
  atoms : Array[Atom]
}

///|
priv struct BcifChainBuilder {
  full_id : String
  chain_id : Char
  residues : Array[BcifResidueBuilder]
}

///|
priv struct BcifModelBuilder {
  model_number : Int
  chains : Array[BcifChainBuilder]
}

///|
fn bcif_require_column(
  category : BinaryCifCategory,
  name : String,
) -> BinaryCifColumn raise BinaryCifError {
  match category.get_column(name) {
    Some(column) => column
    None =>
      raise BinaryCifError(
        "BinaryCIF category '" +
        category.name +
        "' is missing column '" +
        name +
        "'",
      )
  }
}

///|
fn bcif_required_string(
  column : BinaryCifColumn,
  row : Int,
) -> String raise BinaryCifError {
  match column.string_at(row) {
    Some(value) => value
    None =>
      raise BinaryCifError(
        "Required BinaryCIF value '" + column.name + "' is missing",
      )
  }
}

///|
fn bcif_required_int(
  column : BinaryCifColumn,
  row : Int,
) -> Int raise BinaryCifError {
  match column.int_at(row) {
    Some(value) => value
    None =>
      raise BinaryCifError(
        "Required BinaryCIF value '" + column.name + "' is missing",
      )
  }
}

///|
fn bcif_required_double(
  column : BinaryCifColumn,
  row : Int,
) -> Double raise BinaryCifError {
  match column.double_at(row) {
    Some(value) => value
    None =>
      raise BinaryCifError(
        "Required BinaryCIF value '" + column.name + "' is missing",
      )
  }
}

///|
fn bcif_optional_string(
  category : BinaryCifCategory,
  name : String,
  row : Int,
  fallback : String,
) -> String raise BinaryCifError {
  match category.get_column(name) {
    None => fallback
    Some(column) =>
      match column.string_at(row) {
        Some(value) => value
        None => fallback
      }
  }
}

///|
fn bcif_optional_int(
  category : BinaryCifCategory,
  name : String,
  row : Int,
  fallback : Int,
) -> Int raise BinaryCifError {
  match category.get_column(name) {
    None => fallback
    Some(column) =>
      match column.int_at(row) {
        Some(value) => value
        None => fallback
      }
  }
}

///|
fn bcif_optional_double(
  category : BinaryCifCategory,
  name : String,
  row : Int,
  fallback : Double,
) -> Double raise BinaryCifError {
  match category.get_column(name) {
    None => fallback
    Some(column) =>
      match column.double_at(row) {
        Some(value) => value
        None => fallback
      }
  }
}

///|
fn bcif_first_char_or_space(value : String) -> Char {
  if value.length() == 0 || value == "." || value == "?" {
    ' '
  } else {
    value.unsafe_get(0).unsafe_to_char()
  }
}

///|
fn bcif_hetero_field(group : String, component : String) -> String {
  if group == "HETATM" {
    if component == "HOH" || component == "WAT" {
      "W"
    } else {
      "H"
    }
  } else {
    " "
  }
}

///|
fn bcif_entry_id(file : BinaryCifFile) -> String raise BinaryCifError {
  let block = file.data_blocks[0]
  match block.get_category("_entry") {
    Some(category) =>
      match category.get_column("id") {
        Some(column) =>
          match column.string_at(0) {
            Some(value) => value
            None => block.header
          }
        None => block.header
      }
    None => block.header
  }
}

///|
pub fn binary_cif_to_structure(
  file : BinaryCifFile,
  structure_id? : String = "",
) -> Structure raise BinaryCifError {
  if file.data_blocks.length() == 0 {
    raise BinaryCifError("BinaryCIF file has no data blocks")
  }
  let block = file.data_blocks[0]
  let atom_site = match block.get_category("_atom_site") {
    Some(category) => category
    None => raise BinaryCifError("BinaryCIF file has no _atom_site category")
  }
  let names = bcif_require_column(atom_site, "label_atom_id")
  let components = bcif_require_column(atom_site, "label_comp_id")
  let chains = bcif_require_column(atom_site, "label_asym_id")
  let sequence_ids = bcif_require_column(atom_site, "auth_seq_id")
  let xs = bcif_require_column(atom_site, "Cartn_x")
  let ys = bcif_require_column(atom_site, "Cartn_y")
  let zs = bcif_require_column(atom_site, "Cartn_z")
  let builders : Array[BcifModelBuilder] = []
  for row in 0..<atom_site.row_count {
    let model_number = bcif_optional_int(
      atom_site, "pdbx_PDB_model_num", row, 1,
    )
    if builders.length() == 0 ||
      builders[builders.length() - 1].model_number != model_number {
      builders.push(BcifModelBuilder::{ model_number, chains: [], })
    }
    let model_builder = builders[builders.length() - 1]
    let full_chain_id = bcif_required_string(chains, row)
    let chain_id = bcif_first_char_or_space(full_chain_id)
    if model_builder.chains.length() == 0 ||
      model_builder.chains[model_builder.chains.length() - 1].full_id !=
      full_chain_id {
      for existing in model_builder.chains {
        if existing.chain_id == chain_id && existing.full_id != full_chain_id {
          raise BinaryCifError(
            "BinaryCIF chain IDs '" +
            existing.full_id +
            "' and '" +
            full_chain_id +
            "' collide in the character-based PDB model",
          )
        }
      }
      model_builder.chains.push(BcifChainBuilder::{
        full_id: full_chain_id,
        chain_id,
        residues: [],
      })
    }
    let chain_builder = model_builder.chains[model_builder.chains.length() - 1]
    let component = bcif_required_string(components, row)
    let sequence_id = bcif_required_int(sequence_ids, row)
    let insertion_code = bcif_first_char_or_space(
      bcif_optional_string(atom_site, "pdbx_PDB_ins_code", row, ""),
    )
    let hetero_field = bcif_hetero_field(
      bcif_optional_string(atom_site, "group_PDB", row, "ATOM"),
      component,
    )
    let new_residue = chain_builder.residues.length() == 0 ||
      ({
        let previous = chain_builder.residues[chain_builder.residues.length() -
          1]
        previous.sequence_id != sequence_id ||
        previous.insertion_code != insertion_code ||
        previous.hetero_field != hetero_field ||
        previous.resname != component
      })
    if new_residue {
      chain_builder.residues.push(BcifResidueBuilder::{
        resname: component,
        chain_id,
        sequence_id,
        insertion_code,
        hetero_field,
        atoms: [],
      })
    }
    let residue_builder = chain_builder.residues[chain_builder.residues.length() -
      1]
    let atom_name = bcif_required_string(names, row)
    let altloc = bcif_first_char_or_space(
      bcif_optional_string(atom_site, "label_alt_id", row, ""),
    )
    residue_builder.atoms.push(Atom::{
      name: atom_name,
      altloc,
      resname: component,
      chainid: chain_id,
      resseq: sequence_id,
      icode: insertion_code,
      coord: Vector3::new(
        bcif_required_double(xs, row),
        bcif_required_double(ys, row),
        bcif_required_double(zs, row),
      ),
      occupancy: bcif_optional_double(atom_site, "occupancy", row, 1.0),
      bfactor: bcif_optional_double(atom_site, "B_iso_or_equiv", row, 0.0),
      element: bcif_optional_string(atom_site, "type_symbol", row, ""),
      charge: bcif_optional_int(atom_site, "pdbx_formal_charge", row, 0),
    })
  }
  let models = Array(capacity=builders.length())
  for model_index in 0..<builders.length() {
    let chain_values = Array(capacity=builders[model_index].chains.length(),
    )
    for chain_builder in builders[model_index].chains {
      let residues = Array(capacity=chain_builder.residues.length())
      for residue_builder in chain_builder.residues {
        residues.push(Residue::{
          resname: residue_builder.resname,
          chainid: residue_builder.chain_id,
          resseq: residue_builder.sequence_id,
          icode: residue_builder.insertion_code,
          atoms: residue_builder.atoms,
          hetfield: residue_builder.hetero_field,
        })
      }
      chain_values.push(Chain::{ id: chain_builder.chain_id, residues, })
    }
    models.push(Model::{ id: model_index, chains: chain_values, })
  }
  Structure::new(
    id=if structure_id.length() > 0 {
      structure_id
    } else {
      bcif_entry_id(file)
    },
    models~,
  )
}

///|
pub fn binary_cif_parse_structure(
  input : Array[Int],
  structure_id? : String = "",
) -> Structure raise BinaryCifError {
  binary_cif_to_structure(binary_cif_parse(input), structure_id~)
}

// A deterministic BinaryCIF fixture used by tests and examples.

///|
fn bcif_message_map(entries : Array[(String, BcifMessage)]) -> BcifMessage {
  let values : Map[String, BcifMessage] = Map([], capacity=entries.length())
  for entry in entries {
    let (key, value) = entry
    values[key] = value
  }
  BcifMapValue(values)
}

///|
fn bcif_message_array(values : Array[BcifMessage]) -> BcifMessage {
  BcifArrayValue(values)
}

///|
fn bcif_write_u16_be(output : Array[Int], value : Int) -> Unit {
  output.push((value >> 8) & 0xFF)
  output.push(value & 0xFF)
}

///|
fn bcif_write_u32_be(output : Array[Int], value : Int) -> Unit {
  output.push((value >> 24) & 0xFF)
  output.push((value >> 16) & 0xFF)
  output.push((value >> 8) & 0xFF)
  output.push(value & 0xFF)
}

///|
fn bcif_write_i32_le(output : Array[Int], value : Int) -> Unit {
  output.push(value & 0xFF)
  output.push((value >> 8) & 0xFF)
  output.push((value >> 16) & 0xFF)
  output.push((value >> 24) & 0xFF)
}

///|
fn bcif_pack_length(
  output : Array[Int],
  small_base : Int,
  marker16 : Int,
  marker32 : Int,
  length : Int,
) -> Unit {
  if length < 16 {
    output.push(small_base + length)
  } else if length <= 65535 {
    output.push(marker16)
    bcif_write_u16_be(output, length)
  } else {
    output.push(marker32)
    bcif_write_u32_be(output, length)
  }
}

///|
fn bcif_pack_string(output : Array[Int], value : String) -> Unit {
  let length = value.length()
  if length < 32 {
    output.push(0xA0 + length)
  } else if length <= 255 {
    output.push(0xD9)
    output.push(length)
  } else if length <= 65535 {
    output.push(0xDA)
    bcif_write_u16_be(output, length)
  } else {
    output.push(0xDB)
    bcif_write_u32_be(output, length)
  }
  for i in 0..<length {
    output.push(value.unsafe_get(i).to_int() & 0xFF)
  }
}

///|
fn bcif_pack_float64(output : Array[Int], value : Double) -> Unit {
  output.push(0xCB)
  if value == 0.0 {
    for _ in 0..<8 {
      output.push(0)
    }
    return
  }
  let sign = if value < 0.0 { 0x80 } else { 0 }
  let mut normalized = value.abs()
  let mut exponent = 0
  while normalized >= 2.0 {
    normalized = normalized / 2.0
    exponent = exponent + 1
  }
  while normalized < 1.0 {
    normalized = normalized * 2.0
    exponent = exponent - 1
  }
  let biased = exponent + 1023
  let mut fraction = normalized - 1.0
  fraction = fraction * 16.0
  let first_nibble = fraction.to_int()
  fraction = fraction - first_nibble.to_double()
  output.push(sign + (biased >> 4))
  output.push(((biased & 0x0F) << 4) + first_nibble)
  for _ in 0..<6 {
    fraction = fraction * 256.0
    let byte = fraction.to_int()
    output.push(byte)
    fraction = fraction - byte.to_double()
  }
}

///|
fn bcif_pack_message(output : Array[Int], value : BcifMessage) -> Unit {
  match value {
    BcifNil => output.push(0xC0)
    BcifBoolean(false) => output.push(0xC2)
    BcifBoolean(true) => output.push(0xC3)
    BcifIntegerValue(number) =>
      if number >= 0 && number <= 127 {
        output.push(number)
      } else if number >= -32 && number < 0 {
        output.push(number & 0xFF)
      } else {
        output.push(0xD2)
        bcif_write_u32_be(output, number)
      }
    BcifFloatValue(number) => bcif_pack_float64(output, number)
    BcifTextValue(text) => bcif_pack_string(output, text)
    BcifBinaryValue(bytes) => {
      if bytes.length() <= 255 {
        output.push(0xC4)
        output.push(bytes.length())
      } else if bytes.length() <= 65535 {
        output.push(0xC5)
        bcif_write_u16_be(output, bytes.length())
      } else {
        output.push(0xC6)
        bcif_write_u32_be(output, bytes.length())
      }
      for byte in bytes {
        output.push(byte)
      }
    }
    BcifArrayValue(values) => {
      bcif_pack_length(output, 0x90, 0xDC, 0xDD, values.length())
      for item in values {
        bcif_pack_message(output, item)
      }
    }
    BcifMapValue(values) => {
      let keys = values.keys().collect()
      bcif_pack_length(output, 0x80, 0xDE, 0xDF, keys.length())
      for key in keys {
        bcif_pack_string(output, key)
        bcif_pack_message(output, values[key])
      }
    }
  }
}

///|
fn bcif_byte_array_encoding(type_code : Int) -> BcifMessage {
  bcif_message_map([
    ("kind", BcifTextValue("ByteArray")),
    ("type", BcifIntegerValue(type_code)),
  ])
}

///|
fn bcif_sample_int_column(name : String, values : Array[Int]) -> BcifMessage {
  let bytes = Array(capacity=values.length() * 4)
  for value in values {
    bcif_write_i32_le(bytes, value)
  }
  bcif_message_map([
    ("name", BcifTextValue(name)),
    (
      "data",
      bcif_message_map([
        ("data", BcifBinaryValue(bytes)),
        ("encoding", bcif_message_array([bcif_byte_array_encoding(3)])),
      ]),
    ),
  ])
}

///|
fn bcif_sample_fixed_column(
  name : String,
  scaled_values : Array[Int],
  factor : Int,
) -> BcifMessage {
  let deltas = Array(capacity=scaled_values.length())
  let mut previous = 0
  for value in scaled_values {
    deltas.push(value - previous)
    previous = value
  }
  let bytes = Array(capacity=deltas.length() * 4)
  for value in deltas {
    bcif_write_i32_le(bytes, value)
  }
  bcif_message_map([
    ("name", BcifTextValue(name)),
    (
      "data",
      bcif_message_map([
        ("data", BcifBinaryValue(bytes)),
        (
          "encoding",
          bcif_message_array([
            bcif_message_map([
              ("kind", BcifTextValue("FixedPoint")),
              ("factor", BcifIntegerValue(factor)),
              ("srcType", BcifIntegerValue(33)),
            ]),
            bcif_message_map([
              ("kind", BcifTextValue("Delta")),
              ("origin", BcifIntegerValue(0)),
              ("srcType", BcifIntegerValue(3)),
            ]),
            bcif_byte_array_encoding(3),
          ]),
        ),
      ]),
    ),
  ])
}

///|
fn bcif_sample_string_column(
  name : String,
  values : Array[String],
  mask? : Array[Int] = [],
) -> BcifMessage {
  let dictionary : Array[String] = []
  let indices : Array[Int] = []
  let lookup : Map[String, Int] = Map([], capacity=values.length())
  for value in values {
    let index = match lookup.get(value) {
      Some(existing) => existing
      None => {
        let created = dictionary.length()
        dictionary.push(value)
        lookup[value] = created
        created
      }
    }
    indices.push(index)
  }
  let mut string_data = ""
  let offsets = [0]
  for value in dictionary {
    string_data = string_data + value
    offsets.push(string_data.length())
  }
  let offset_bytes = Array(capacity=offsets.length() * 4)
  for value in offsets {
    bcif_write_i32_le(offset_bytes, value)
  }
  let index_bytes = Array(capacity=indices.length())
  for value in indices {
    index_bytes.push(value)
  }
  let entries : Array[(String, BcifMessage)] = [
    ("name", BcifTextValue(name)),
    (
      "data",
      bcif_message_map([
        ("data", BcifBinaryValue(index_bytes)),
        (
          "encoding",
          bcif_message_array([
            bcif_message_map([
              ("kind", BcifTextValue("StringArray")),
              (
                "dataEncoding",
                bcif_message_array([bcif_byte_array_encoding(4)]),
              ),
              ("stringData", BcifTextValue(string_data)),
              (
                "offsetEncoding",
                bcif_message_array([bcif_byte_array_encoding(3)]),
              ),
              ("offsets", BcifBinaryValue(offset_bytes)),
            ]),
          ]),
        ),
      ]),
    ),
  ]
  if mask.length() > 0 {
    entries.push(
      (
        "mask",
        bcif_message_map([
          ("data", BcifBinaryValue(mask)),
          ("encoding", bcif_message_array([bcif_byte_array_encoding(4)])),
        ]),
      ),
    )
  }
  bcif_message_map(entries)
}

///|
pub fn binary_cif_sample_bytes() -> Array[Int] {
  let atom_count = 8
  let atom_columns = [
    bcif_sample_string_column("group_PDB", [
      "ATOM", "ATOM", "ATOM", "ATOM", "ATOM", "HETATM", "ATOM", "ATOM",
    ]),
    bcif_sample_int_column("id", [1, 2, 3, 4, 5, 6, 7, 8]),
    bcif_sample_string_column("type_symbol", [
      "N", "C", "C", "N", "C", "O", "N", "C",
    ]),
    bcif_sample_string_column("label_atom_id", [
      "N", "CA", "C", "N", "CA", "O", "N", "CA",
    ]),
    bcif_sample_string_column(
      "label_alt_id",
      ["", "", "", "A", "A", "", "", ""],
      mask=[1, 1, 1, 0, 0, 1, 2, 1],
    ),
    bcif_sample_string_column("label_comp_id", [
      "GLY", "GLY", "GLY", "ALA", "ALA", "HOH", "SER", "SER",
    ]),
    bcif_sample_string_column("label_asym_id", [
      "A", "A", "A", "A", "A", "B", "A", "A",
    ]),
    bcif_sample_int_column("auth_seq_id", [1, 1, 1, 2, 2, 10, 1, 1]),
    bcif_sample_string_column(
      "pdbx_PDB_ins_code",
      ["", "", "", "A", "A", "", "", ""],
      mask=[1, 1, 1, 0, 0, 1, 1, 1],
    ),
    bcif_sample_fixed_column(
      "Cartn_x",
      [100, 220, 340, 460, 580, 700, 820, 940],
      100,
    ),
    bcif_sample_fixed_column(
      "Cartn_y",
      [200, 300, 400, 500, 600, 700, 800, 900],
      100,
    ),
    bcif_sample_fixed_column(
      "Cartn_z",
      [-100, 0, 100, 200, 300, 400, 500, 600],
      100,
    ),
    bcif_sample_fixed_column(
      "occupancy",
      [100, 100, 100, 50, 50, 100, 100, 100],
      100,
    ),
    bcif_sample_fixed_column(
      "B_iso_or_equiv",
      [1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700],
      100,
    ),
    bcif_sample_int_column("pdbx_PDB_model_num", [1, 1, 1, 1, 1, 1, 2, 2]),
  ]
  let document = bcif_message_map([
    ("version", BcifTextValue("0.3.0")),
    ("encoder", BcifTextValue("BioSeqs MoonBit")),
    (
      "dataBlocks",
      bcif_message_array([
        bcif_message_map([
          ("header", BcifTextValue("BCIF")),
          (
            "categories",
            bcif_message_array([
              bcif_message_map([
                ("name", BcifTextValue("_entry")),
                ("rowCount", BcifIntegerValue(1)),
                (
                  "columns",
                  bcif_message_array([bcif_sample_string_column("id", ["BCIF"])]),
                ),
              ]),
              bcif_message_map([
                ("name", BcifTextValue("_atom_site")),
                ("rowCount", BcifIntegerValue(atom_count)),
                ("columns", bcif_message_array(atom_columns)),
              ]),
            ]),
          ),
        ]),
      ]),
    ),
  ])
  let output : Array[Int] = []
  bcif_pack_message(output, document)
  output
}
