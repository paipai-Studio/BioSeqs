///|
/// DCD (CHARMM/NAMD/X-PLOR binary trajectory) + a generic Trajectory type.
///
/// Writes and reads the X-PLOR DCD variant (the common subset every reader —
/// VMD, MDAnalysis, gromacs — accepts): a `CORD` header, two 80-byte title
/// records, a NATOM record, then per-frame float32 little-endian `xyz`.

///|
pub struct Trajectory {
  frames : Array[Array[Array[Double]]]
}

///|
pub fn Trajectory::new(frames : Array[Array[Array[Double]]]) -> Trajectory {
  Trajectory::{ frames }
}

///|
pub fn Trajectory::num_frames(self : Trajectory) -> Int {
  self.frames.length()
}

///|
pub fn Trajectory::num_atoms(self : Trajectory) -> Int {
  if self.frames.length() == 0 { 0 } else { self.frames[0].length() }
}

///|
pub fn Trajectory::get_frame(self : Trajectory, i : Int) -> Array[Array[Double]] {
  self.frames[i]
}

// ---------------------------------------------------------------------------
// Little-endian binary helpers
// ---------------------------------------------------------------------------

///|
fn dcd_read_int32_le(data : String, pos : Int) -> (Int, Int) {
  let b0 = data.unsafe_get(pos).to_int() & 0xFF
  let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
  let b2 = data.unsafe_get(pos + 2).to_int() & 0xFF
  let b3 = data.unsafe_get(pos + 3).to_int() & 0xFF
  (b0 | (b1 << 8) | (b2 << 16) | (b3 << 24), pos + 4)
}

///|
fn dcd_read_f32_le(data : String, pos : Int) -> (Double, Int) {
  let b0 = data.unsafe_get(pos).to_int() & 0xFF
  let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
  let b2 = data.unsafe_get(pos + 2).to_int() & 0xFF
  let b3 = data.unsafe_get(pos + 3).to_int() & 0xFF
  let sign = if b3 >= 128 { -1.0 } else { 1.0 }
  let exponent = (b3 & 0x7F) * 2 + (b2 >> 7)
  let mantissa = (b2 & 0x7F) * 65536 + b1 * 256 + b0
  let v = if exponent == 255 {
    0.0 / 0.0
  } else if exponent == 0 {
    sign * mantissa.to_double() * @math.pow(2.0, -149.0)
  } else {
    sign * (1.0 + mantissa.to_double() / 8388608.0) * @math.pow(2.0, (exponent - 127).to_double())
  }
  (v, pos + 4)
}

///|
fn dcd_write_i32(sb : StringBuilder, v : Int) -> Unit {
  sb.write_char(((v & 0xFF).to_uint16()).unsafe_to_char())
  sb.write_char((((v >> 8) & 0xFF).to_uint16()).unsafe_to_char())
  sb.write_char((((v >> 16) & 0xFF).to_uint16()).unsafe_to_char())
  sb.write_char((((v >> 24) & 0xFF).to_uint16()).unsafe_to_char())
}

///|
fn dcd_f32_bits(v : Double) -> Int {
  if v != v {
    return 0x7FC00000
  }
  if v == 0.0 {
    return if (1.0 / v) < 0.0 { 0x80000000 } else { 0 }
  }
  if v > 3.4028234663852886e38 {
    return 0x7F800000
  }
  let sign = if v < 0.0 { 1 } else { 0 }
  let a = if v < 0.0 { 0.0 - v } else { v }
  let mut e = 0
  let mut m = a
  if m >= 1.0 {
    while m >= 2.0 {
      m = m / 2.0
      e = e + 1
    }
  } else {
    while m < 1.0 {
      m = m * 2.0
      e = e - 1
    }
  }
  let mut exp = e + 127
  if exp <= 0 {
    return 0
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
fn dcd_write_f32(sb : StringBuilder, v : Double) -> Unit {
  dcd_write_i32(sb, dcd_f32_bits(v))
}

///|
fn dcd_write_title(sb : StringBuilder, n : Int) -> Unit {
  let title = "BioSeqs DCD"
  let mut i = 0
  while i < n {
    if i < title.length() {
      sb.write_char(title.unsafe_get(i).unsafe_to_char())
    } else {
      sb.write_char('\u{0000}')
    }
    i = i + 1
  }
}

// ---------------------------------------------------------------------------
// Write / read
// ---------------------------------------------------------------------------

///|
/// Serialize a trajectory to X-PLOR DCD bytes.
pub fn write_dcd(traj : Trajectory) -> String {
  let sb = StringBuilder()
  let n_frames = traj.frames.length()
  let n_atoms = traj.num_atoms()

  // header record: [84]["CORD" + 20 int32][84]
  dcd_write_i32(sb, 84)
  sb.write_string("CORD")
  dcd_write_i32(sb, n_frames) // NSET
  dcd_write_i32(sb, 0)        // ISTART
  dcd_write_i32(sb, 1)        // NSAVC
  dcd_write_i32(sb, 0)        // NAMNF
  let mut k = 0
  while k < 5 {
    dcd_write_i32(sb, 0)
    k = k + 1
  }
  dcd_write_i32(sb, 0)        // NDEGF
  let mut k2 = 0
  while k2 < 4 {
    dcd_write_i32(sb, 0)
    k2 = k2 + 1
  }
  dcd_write_i32(sb, n_atoms)  // NATOM
  dcd_write_i32(sb, n_atoms)  // NFREAT
  let mut k3 = 0
  while k3 < 3 {
    dcd_write_i32(sb, 0)
    k3 = k3 + 1
  }
  dcd_write_i32(sb, 0)        // delt
  dcd_write_i32(sb, 84)

  // two title records (80 bytes each)
  let mut t = 0
  while t < 2 {
    dcd_write_i32(sb, 80)
    dcd_write_title(sb, 80)
    dcd_write_i32(sb, 80)
    t = t + 1
  }

  // NATOM record
  dcd_write_i32(sb, 4)
  dcd_write_i32(sb, n_atoms)
  dcd_write_i32(sb, 4)

  // frames
  for frame in traj.frames {
    let frame_len = 4 * n_atoms * 3
    dcd_write_i32(sb, frame_len)
    for atom in frame {
      for coord in atom {
        dcd_write_f32(sb, coord)
      }
    }
    dcd_write_i32(sb, frame_len)
  }

  sb.to_string()
}

///|
/// Parse X-PLOR DCD bytes into a Trajectory.
pub fn parse_dcd(content : String) -> Trajectory {
  let mut pos = 0
  let n = content.length()

  // header record
  let (hdr_len, p1) = dcd_read_int32_le(content, pos)
  pos = p1
  // magic "CORD"
  pos = pos + 4
  // 20 int32 header fields（读回 NSET 和 NATOM）
  let (n_set, p2) = dcd_read_int32_le(content, pos)
  pos = p2
  // 跳到 NATOM：字段 [1]..[14] 共 14 个 int32（位置 16..72）
  pos = pos + 13 * 4
  let (n_atoms, p3) = dcd_read_int32_le(content, pos)
  pos = p3
  // 剩下的 header 字段 + 后缀（delt 等）
  pos = pos + 5 * 4 + 4
  // 跳过 title records ×2 + NATOM record
  let (t1, _) = dcd_read_int32_le(content, pos)
  pos = pos + 4 + t1 + 4
  let (t2, _) = dcd_read_int32_le(content, pos)
  pos = pos + 4 + t2 + 4
  let (_nrec, _) = dcd_read_int32_le(content, pos)
  pos = pos + 4 + 4 + 4

  // frames
  let frames : Array[Array[Array[Double]]] = Array(capacity=0)
  let mut f = 0
  while f < n_set {
    let (flen, p4) = dcd_read_int32_le(content, pos)
    pos = p4
    let frame : Array[Array[Double]] = Array(capacity=0)
    let mut a = 0
    while a < n_atoms {
      let coord : Array[Double] = Array(capacity=0)
      let mut d = 0
      while d < 3 {
        let (v, p5) = dcd_read_f32_le(content, pos)
        pos = p5
        coord.push(v)
        d = d + 1
      }
      frame.push(coord)
      a = a + 1
    }
    frames.push(frame)
    pos = pos + 4 // frame size suffix
    f = f + 1
  }
  ignore(hdr_len)
  Trajectory::new(frames)
}
