///|
/// GROMACS TRR (full-precision trajectory) reader/writer.
///
/// TRR frames are XDR (big-endian) encoded. Every frame stores: magic 1993,
/// a 13-byte version string ("GMX_trn_file\0"), ten int32 size fields
/// (ir/eg/box/vir/pres/top/sym/x/v/f), int32 natoms and step, float64 time
/// and lambda, then optional float32 blocks: a 9-value box (3x3 row major,
/// nm), 9-value virial, 9-value pressure, natoms*3 coordinates, velocities
/// and forces.
///
/// NOTE: this writer never emits virial/pressure blocks; the reader parses
/// and skips them when present.

///|
/// Error type for TRR operations.
pub suberror TrrError {
  TrrError(String)
}

///|
/// A single frame (snapshot) from a TRR trajectory.
/// `box` holds 9 doubles (3x3, row major, nm) and is empty when absent.
/// `coords`/`velocities`/`forces` are flat 3*natoms arrays (empty if absent).
pub(all) struct TrrFrame {
  step : Int32
  time : Double
  lambda : Double
  box : Array[Double]
  coords : Array[Double]
  velocities : Array[Double]
  forces : Array[Double]
}

// ---------------------------------------------------------------------------
// IEEE 754 bit conversion (math decomposition, see dcd_f32_bits)
// ---------------------------------------------------------------------------

///|
fn trr_pow2(e : Int32) -> Double {
  if e >= 0 {
    let mut result = 1.0
    for i = 0; i < e; i = i + 1 {
      result = result * 2.0
    }
    result
  } else {
    let mut result = 1.0
    for i = 0; i < -e; i = i + 1 {
      result = result / 2.0
    }
    result
  }
}

///|
fn trr_f32_bits(v : Double) -> Int32 {
  if v != v {
    return 0x7FC00000
  }
  if v == 0.0 {
    return if 1.0 / v < 0.0 { 0x80000000 } else { 0 }
  }
  let sign = if v < 0.0 { 1 } else { 0 }
  let a = if v < 0.0 { 0.0 - v } else { v }
  if a > 3.4028234663852886e38 {
    return (sign << 31) | 0x7F800000
  }
  let mut e = 0
  let mut m = a
  if m >= 1.0 {
    while m >= 2.0 && e < 128 {
      m = m / 2.0
      e = e + 1
    }
  } else {
    while m < 1.0 && e > -127 {
      m = m * 2.0
      e = e - 1
    }
  }
  let mut exp = e + 127
  if exp <= 0 {
    return sign << 31
  }
  let mut mant = ((m - 1.0) * 8388608.0 + 0.5).to_int()
  if mant >= 8388608 {
    mant = 0
    exp = exp + 1
  }
  if exp >= 255 {
    return (sign << 31) | 0x7F800000
  }
  (sign << 31) | (exp << 23) | (mant & 0x7FFFFF)
}

///|
fn trr_f32_from_bits(bits : Int32) -> Double {
  let sign = if ((bits >> 31) & 1) == 1 { -1.0 } else { 1.0 }
  let exp = (bits >> 23) & 0xFF
  let mant = bits & 0x7FFFFF
  if exp == 0 {
    if mant == 0 {
      return sign * 0.0
    }
    return sign * mant.to_double() * trr_pow2(-149)
  }
  if exp == 255 {
    if mant == 0 {
      return sign * 1.0e308 * 10.0
    }
    return 0.0 / 0.0
  }
  sign * (1.0 + mant.to_double() / 8388608.0) * trr_pow2(exp - 127)
}

///|
fn trr_f64_bits(v : Double) -> UInt64 {
  if v != v {
    return 0x7FF8000000000000UL
  }
  if v == 0.0 {
    return if 1.0 / v < 0.0 { 0x8000000000000000UL } else { 0UL }
  }
  let sign = if v < 0.0 { 1UL } else { 0UL }
  let a = if v < 0.0 { 0.0 - v } else { v }
  if a > 1.7976931348623157e308 {
    return (sign << 63) | 0x7FF0000000000000UL
  }
  let mut e = 0
  let mut m = a
  if m >= 1.0 {
    while m >= 2.0 && e < 1024 {
      m = m / 2.0
      e = e + 1
    }
  } else {
    while m < 1.0 && e > -1023 {
      m = m * 2.0
      e = e - 1
    }
  }
  let mut exp = e + 1023
  if exp <= 0 {
    return sign << 63
  }
  let mut mant = ((m - 1.0) * 4503599627370496.0).to_uint64()
  if mant >= 4503599627370496UL {
    mant = 0UL
    exp = exp + 1
  }
  if exp >= 2047 {
    return (sign << 63) | 0x7FF0000000000000UL
  }
  (sign << 63) | (exp.to_uint64() << 52) | (mant & 0xFFFFFFFFFFFFFUL)
}

///|
fn trr_f64_from_bits(bits : UInt64) -> Double {
  let sign = if ((bits >> 63) & 1UL) == 1UL { -1.0 } else { 1.0 }
  let exp = ((bits >> 52) & 0x7FFUL).to_int()
  let mant = bits & 0xFFFFFFFFFFFFFUL
  if exp == 0 {
    if mant == 0UL {
      return sign * 0.0
    }
    return sign * mant.to_double() * trr_pow2(-1074)
  }
  if exp == 2047 {
    if mant == 0UL {
      return sign * 1.0e308 * 10.0
    }
    return 0.0 / 0.0
  }
  sign * (1.0 + mant.to_double() / 4503599627370496.0) * trr_pow2(exp - 1023)
}

// ---------------------------------------------------------------------------
// XDR read helpers (big-endian, byte buffer = Array[Int32], values 0-255)
// ---------------------------------------------------------------------------

///|
fn trr_read_i32(
  data : Array[Int32],
  pos : Int32,
) -> (Int32, Int32) raise TrrError {
  if pos + 4 > data.length() {
    raise TrrError("unexpected end of data reading int at pos \{pos}")
  }
  let b0 = data[pos]
  let b1 = data[pos + 1]
  let b2 = data[pos + 2]
  let b3 = data[pos + 3]
  ((b0 << 24) | (b1 << 16) | (b2 << 8) | b3, pos + 4)
}

///|
fn trr_read_f32(
  data : Array[Int32],
  pos : Int32,
) -> (Double, Int32) raise TrrError {
  if pos + 4 > data.length() {
    raise TrrError("unexpected end of data reading float32 at pos \{pos}")
  }
  let b0 = data[pos]
  let b1 = data[pos + 1]
  let b2 = data[pos + 2]
  let b3 = data[pos + 3]
  let bits = (b0 << 24) | (b1 << 16) | (b2 << 8) | b3
  (trr_f32_from_bits(bits), pos + 4)
}

///|
fn trr_read_f64(
  data : Array[Int32],
  pos : Int32,
) -> (Double, Int32) raise TrrError {
  if pos + 8 > data.length() {
    raise TrrError("unexpected end of data reading float64 at pos \{pos}")
  }
  let b0 = data[pos].to_uint64()
  let b1 = data[pos + 1].to_uint64()
  let b2 = data[pos + 2].to_uint64()
  let b3 = data[pos + 3].to_uint64()
  let b4 = data[pos + 4].to_uint64()
  let b5 = data[pos + 5].to_uint64()
  let b6 = data[pos + 6].to_uint64()
  let b7 = data[pos + 7].to_uint64()
  let bits = (b0 << 56) |
    (b1 << 48) |
    (b2 << 40) |
    (b3 << 32) |
    (b4 << 24) |
    (b5 << 16) |
    (b6 << 8) |
    b7
  (trr_f64_from_bits(bits), pos + 8)
}

// ---------------------------------------------------------------------------
// XDR write helpers
// ---------------------------------------------------------------------------

///|
fn trr_write_i32(buf : Array[Int32], v : Int32) -> Unit {
  buf.push((v >> 24) & 0xFF)
  buf.push((v >> 16) & 0xFF)
  buf.push((v >> 8) & 0xFF)
  buf.push(v & 0xFF)
}

///|
fn trr_write_f32(buf : Array[Int32], v : Double) -> Unit {
  trr_write_i32(buf, trr_f32_bits(v))
}

///|
fn trr_write_f64(buf : Array[Int32], v : Double) -> Unit {
  let bits = trr_f64_bits(v)
  buf.push(((bits >> 56) & 0xFFUL).to_int())
  buf.push(((bits >> 48) & 0xFFUL).to_int())
  buf.push(((bits >> 40) & 0xFFUL).to_int())
  buf.push(((bits >> 32) & 0xFFUL).to_int())
  buf.push(((bits >> 24) & 0xFFUL).to_int())
  buf.push(((bits >> 16) & 0xFFUL).to_int())
  buf.push(((bits >> 8) & 0xFFUL).to_int())
  buf.push((bits & 0xFFUL).to_int())
}

// ---------------------------------------------------------------------------
// Write / read
// ---------------------------------------------------------------------------

///|
/// Serialize one TRR frame to XDR bytes (big-endian).
pub fn trr_write_frame(frame : TrrFrame) -> Array[Int32] {
  let buf : Array[Int32] = Array(capacity=0)
  let has_box = frame.box.length() >= 9
  let box_size = if has_box { 72 } else { 0 }
  let x_size = frame.coords.length() * 4
  let v_size = frame.velocities.length() * 4
  let f_size = frame.forces.length() * 4
  let natoms = if frame.coords.length() > 0 {
    frame.coords.length() / 3
  } else if frame.velocities.length() > 0 {
    frame.velocities.length() / 3
  } else if frame.forces.length() > 0 {
    frame.forces.length() / 3
  } else {
    0
  }
  trr_write_i32(buf, 1993)
  trr_write_i32(buf, 13)
  let version = "GMX_trn_file"
  for i = 0; i < 12; i = i + 1 {
    buf.push(version.unsafe_get(i).to_int() & 0xFF)
  }
  buf.push(0)
  trr_write_i32(buf, 0) // ir_size
  trr_write_i32(buf, 0) // eg_size
  trr_write_i32(buf, box_size)
  trr_write_i32(buf, 0) // vir_size
  trr_write_i32(buf, 0) // pres_size
  trr_write_i32(buf, 0) // top_size
  trr_write_i32(buf, 0) // sym_size
  trr_write_i32(buf, x_size)
  trr_write_i32(buf, v_size)
  trr_write_i32(buf, f_size)
  trr_write_i32(buf, natoms)
  trr_write_i32(buf, frame.step)
  trr_write_f64(buf, frame.time)
  trr_write_f64(buf, frame.lambda)
  if has_box {
    for i = 0; i < 9; i = i + 1 {
      trr_write_f32(buf, frame.box[i])
    }
  }
  for i = 0; i < frame.coords.length(); i = i + 1 {
    trr_write_f32(buf, frame.coords[i])
  }
  for i = 0; i < frame.velocities.length(); i = i + 1 {
    trr_write_f32(buf, frame.velocities[i])
  }
  for i = 0; i < frame.forces.length(); i = i + 1 {
    trr_write_f32(buf, frame.forces[i])
  }
  buf
}

///|
/// Parse a TRR byte stream (possibly containing multiple frames) and return
/// all frames. Raises TrrError on malformed data (bad magic, truncated
/// frame, invalid version string size).
pub fn trr_read(data : Array[Int32]) -> Array[TrrFrame] raise TrrError {
  let frames : Array[TrrFrame] = Array(capacity=0)
  let len = data.length()
  let mut pos = 0
  while pos < len {
    let (magic, p1) = trr_read_i32(data, pos)
    if magic != 1993 {
      raise TrrError("invalid TRR magic \{magic}, expected 1993")
    }
    pos = p1
    let (vsz, p2) = trr_read_i32(data, pos)
    if vsz != 13 {
      raise TrrError("invalid TRR version string size \{vsz}, expected 13")
    }
    pos = p2
    pos = pos + 13
    let (ir_size, p3) = trr_read_i32(data, pos)
    ignore(ir_size)
    pos = p3
    let (eg_size, p4) = trr_read_i32(data, pos)
    ignore(eg_size)
    pos = p4
    let (box_size, p5) = trr_read_i32(data, pos)
    pos = p5
    let (vir_size, p6) = trr_read_i32(data, pos)
    pos = p6
    let (pres_size, p7) = trr_read_i32(data, pos)
    pos = p7
    let (top_size, p8) = trr_read_i32(data, pos)
    ignore(top_size)
    pos = p8
    let (sym_size, p9) = trr_read_i32(data, pos)
    ignore(sym_size)
    pos = p9
    let (x_size, p10) = trr_read_i32(data, pos)
    pos = p10
    let (v_size, p11) = trr_read_i32(data, pos)
    pos = p11
    let (f_size, p12) = trr_read_i32(data, pos)
    pos = p12
    let (natoms, p13) = trr_read_i32(data, pos)
    ignore(natoms)
    pos = p13
    let (step, p14) = trr_read_i32(data, pos)
    pos = p14
    let (time, p15) = trr_read_f64(data, pos)
    pos = p15
    let (lambda, p16) = trr_read_f64(data, pos)
    pos = p16
    let box : Array[Double] = Array(capacity=0)
    if box_size > 0 {
      for i = 0; i < 9; i = i + 1 {
        let (v, pb) = trr_read_f32(data, pos)
        pos = pb
        box.push(v)
      }
    }
    if vir_size > 0 {
      pos = pos + vir_size
    }
    if pres_size > 0 {
      pos = pos + pres_size
    }
    let coords : Array[Double] = Array(capacity=0)
    if x_size > 0 {
      let nvals = x_size / 4
      for i = 0; i < nvals; i = i + 1 {
        let (v, px) = trr_read_f32(data, pos)
        pos = px
        coords.push(v)
      }
    }
    let velocities : Array[Double] = Array(capacity=0)
    if v_size > 0 {
      let nvals = v_size / 4
      for i = 0; i < nvals; i = i + 1 {
        let (v, pv) = trr_read_f32(data, pos)
        pos = pv
        velocities.push(v)
      }
    }
    let forces : Array[Double] = Array(capacity=0)
    if f_size > 0 {
      let nvals = f_size / 4
      for i = 0; i < nvals; i = i + 1 {
        let (v, pf) = trr_read_f32(data, pos)
        pos = pf
        forces.push(v)
      }
    }
    frames.push(TrrFrame::{
      step,
      time,
      lambda,
      box,
      coords,
      velocities,
      forces,
    })
  }
  frames
}

///|
/// Convert TRR frames into a Trajectory. Frames whose `coords` array is
/// empty are skipped; every kept frame must contain 3*natoms coordinates.
pub fn trr_to_trajectory(frames : Array[TrrFrame]) -> Trajectory {
  let out : Array[Array[Array[Double]]] = Array(capacity=0)
  for f in frames {
    if f.coords.length() >= 3 {
      let natoms = f.coords.length() / 3
      let frame : Array[Array[Double]] = Array(capacity=natoms)
      for a = 0; a < natoms; a = a + 1 {
        frame.push([f.coords[a * 3], f.coords[a * 3 + 1], f.coords[a * 3 + 2]])
      }
      out.push(frame)
    }
  }
  Trajectory::new(out)
}

///|
/// Convert a Trajectory into TRR frames. Frame `i` gets step
/// `step + i` and time `time_step * i`. No box, velocities or forces
/// are attached.
pub fn trajectory_to_trr(
  traj : Trajectory,
  step~ : Int32 = 0,
  time_step~ : Double = 0.0,
) -> Array[TrrFrame] {
  let out : Array[TrrFrame] = Array(capacity=0)
  let mut idx = 0
  for frame in traj.frames {
    let coords : Array[Double] = Array(capacity=frame.length() * 3)
    for atom in frame {
      for d = 0; d < 3; d = d + 1 {
        if d < atom.length() {
          coords.push(atom[d])
        } else {
          coords.push(0.0)
        }
      }
    }
    out.push(TrrFrame::{
      step: step + idx,
      time: time_step * idx.to_double(),
      lambda: 0.0,
      box: [],
      coords,
      velocities: [],
      forces: [],
    })
    idx = idx + 1
  }
  out
}
