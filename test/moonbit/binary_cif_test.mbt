///|
/// Tests for the Biopython Bio.PDB.binary_cif-compatible decoder.

///|
test "binary cif decodes signed and unsigned 8-bit byte arrays" {
  let bytes = [0, 127, 128, 255]
  assert_eq(@src.binary_cif_decode_int_bytes(bytes, 1), [0, 127, -128, -1])
  assert_eq(@src.binary_cif_decode_int_bytes(bytes, 4), [0, 127, 128, 255])
}

///|
test "binary cif decodes signed and unsigned 16-bit byte arrays" {
  let bytes = [0x34, 0x12, 0x00, 0x80, 0xFF, 0xFF]
  assert_eq(
    @src.binary_cif_decode_int_bytes(bytes, 2),
    [4660, -32768, -1],
  )
  assert_eq(
    @src.binary_cif_decode_int_bytes(bytes, 5),
    [4660, 32768, 65535],
  )
}

///|
test "binary cif decodes signed and unsigned 32-bit byte arrays" {
  assert_eq(
    @src.binary_cif_decode_int_bytes(
      [0x78, 0x56, 0x34, 0x12, 0xFF, 0xFF, 0xFF, 0xFF],
      3,
    ),
    [305419896, -1],
  )
  assert_eq(
    @src.binary_cif_decode_int_bytes([0xFF, 0xFF, 0xFF, 0x7F], 6),
    [2147483647],
  )
}

///|
test "binary cif decodes little-endian float32 byte arrays" {
  let values = @src.binary_cif_decode_float_bytes([
    0x00, 0x00, 0x80, 0x3F, 0x00, 0x00, 0x20, 0xC0,
  ], 32)
  assert_eq(values.length(), 2)
  assert_true((values[0] - 1.0).abs() < 1.0e-12)
  assert_true((values[1] + 2.5).abs() < 1.0e-12)
}

///|
test "binary cif decodes little-endian float64 byte arrays" {
  let values = @src.binary_cif_decode_float_bytes([
    0x00,
    0x00,
    0x00,
    0x00,
    0x00,
    0x00,
    0xF0,
    0x3F,
    0x00,
    0x00,
    0x00,
    0x00,
    0x00,
    0x00,
    0x04,
    0xC0,
  ], 33)
  assert_eq(values.length(), 2)
  assert_true((values[0] - 1.0).abs() < 1.0e-12)
  assert_true((values[1] + 2.5).abs() < 1.0e-12)
}

///|
test "binary cif byte array decoding rejects invalid inputs" {
  let bad_length = try {
    ignore(@src.binary_cif_decode_int_bytes([1], 2))
    false
  } catch {
    BinaryCifError(_) => true
  }
  let bad_type = try {
    ignore(@src.binary_cif_decode_float_bytes([0, 0, 0, 0], 31))
    false
  } catch {
    BinaryCifError(_) => true
  }
  let non_byte = try {
    ignore(@src.binary_cif_decode_int_bytes([256], 4))
    false
  } catch {
    BinaryCifError(_) => true
  }
  assert_true(bad_length)
  assert_true(bad_type)
  assert_true(non_byte)
}

///|
test "binary cif decodes unsigned integer packing" {
  assert_eq(
    @src.binary_cif_decode_integer_packing([255, 5, 1], 1, true, 2),
    [260, 1],
  )
}

///|
test "binary cif decodes signed integer packing" {
  assert_eq(
    @src.binary_cif_decode_integer_packing(
      [127, 5, -128, -2, 10],
      1,
      false,
      3,
    ),
    [132, -130, 10],
  )
}

///|
test "binary cif decodes two-byte integer packing" {
  assert_eq(
    @src.binary_cif_decode_integer_packing([65535, 100], 2, true, 1),
    [65635],
  )
  assert_eq(
    @src.binary_cif_decode_integer_packing([-32768, -5], 2, false, 1),
    [-32773],
  )
}

///|
test "binary cif integer packing validates metadata and values" {
  let bad_width = try {
    ignore(@src.binary_cif_decode_integer_packing([1], 4, true, 1))
    false
  } catch {
    BinaryCifError(_) => true
  }
  let bad_range = try {
    ignore(@src.binary_cif_decode_integer_packing([256], 1, true, 1))
    false
  } catch {
    BinaryCifError(_) => true
  }
  let bad_size = try {
    ignore(@src.binary_cif_decode_integer_packing([1], 1, true, 2))
    false
  } catch {
    BinaryCifError(_) => true
  }
  assert_true(bad_width)
  assert_true(bad_range)
  assert_true(bad_size)
}

///|
test "binary cif decodes run-length encoding" {
  assert_eq(
    @src.binary_cif_decode_run_length([4, 3, -1, 2], 5),
    [4, 4, 4, -1, -1],
  )
}

///|
test "binary cif run-length decoding rejects malformed data" {
  let odd = try {
    ignore(@src.binary_cif_decode_run_length([4, 3, 1], 3))
    false
  } catch {
    BinaryCifError(_) => true
  }
  let overflow = try {
    ignore(@src.binary_cif_decode_run_length([4, 4], 3))
    false
  } catch {
    BinaryCifError(_) => true
  }
  let mismatch = try {
    ignore(@src.binary_cif_decode_run_length([4, 2], 3))
    false
  } catch {
    BinaryCifError(_) => true
  }
  assert_true(odd)
  assert_true(overflow)
  assert_true(mismatch)
}

///|
test "binary cif decodes delta encoding" {
  assert_eq(@src.binary_cif_decode_delta([1, 2, -1, 5], 10), [11, 13, 12, 17])
}

///|
test "binary cif decodes fixed-point encoding" {
  let values = @src.binary_cif_decode_fixed_point([100, -25, 0], 10.0)
  assert_eq(values.length(), 3)
  assert_true((values[0] - 10.0).abs() < 1.0e-12)
  assert_true((values[1] + 2.5).abs() < 1.0e-12)
  assert_true(values[2].abs() < 1.0e-12)
}

///|
test "binary cif fixed-point decoding rejects zero factor" {
  let raised = try {
    ignore(@src.binary_cif_decode_fixed_point([1], 0.0))
    false
  } catch {
    BinaryCifError(_) => true
  }
  assert_true(raised)
}

///|
test "binary cif decodes interval quantization" {
  let values = @src.binary_cif_decode_interval_quantization(
    [0, 2, 4],
    0.0,
    1.0,
    5,
  )
  assert_eq(values.length(), 3)
  assert_true(values[0].abs() < 1.0e-12)
  assert_true((values[1] - 0.5).abs() < 1.0e-12)
  assert_true((values[2] - 1.0).abs() < 1.0e-12)
}

///|
test "binary cif interval quantization validates parameters and indices" {
  let bad_steps = try {
    ignore(
      @src.binary_cif_decode_interval_quantization([0], 0.0, 1.0, 1),
    )
    false
  } catch {
    BinaryCifError(_) => true
  }
  let bad_bounds = try {
    ignore(
      @src.binary_cif_decode_interval_quantization([0], 2.0, 1.0, 3),
    )
    false
  } catch {
    BinaryCifError(_) => true
  }
  let bad_index = try {
    ignore(
      @src.binary_cif_decode_interval_quantization([3], 0.0, 1.0, 3),
    )
    false
  } catch {
    BinaryCifError(_) => true
  }
  assert_true(bad_steps)
  assert_true(bad_bounds)
  assert_true(bad_index)
}

///|
test "binary cif decodes dictionary string arrays" {
  assert_eq(
    @src.binary_cif_decode_string_array([1, 0, 1], [0, 3, 6], "ALAGLY"),
    ["GLY", "ALA", "GLY"],
  )
}

///|
test "binary cif string arrays validate offsets and indices" {
  let bad_start = try {
    ignore(@src.binary_cif_decode_string_array([0], [1, 3], "ALA"))
    false
  } catch {
    BinaryCifError(_) => true
  }
  let bad_end = try {
    ignore(@src.binary_cif_decode_string_array([0], [0, 2], "ALA"))
    false
  } catch {
    BinaryCifError(_) => true
  }
  let bad_index = try {
    ignore(@src.binary_cif_decode_string_array([1], [0, 3], "ALA"))
    false
  } catch {
    BinaryCifError(_) => true
  }
  assert_true(bad_start)
  assert_true(bad_end)
  assert_true(bad_index)
}

///|
test "binary cif parses sample document metadata" {
  let file = @src.binary_cif_parse(@src.binary_cif_sample_bytes())
  assert_eq(file.version, "0.3.0")
  assert_eq(file.encoder, "BioSeqs MoonBit")
  assert_eq(file.data_blocks.length(), 1)
  assert_eq(file.data_blocks[0].header, "BCIF")
}

///|
test "binary cif queries blocks and categories" {
  let file = @src.binary_cif_parse(@src.binary_cif_sample_bytes())
  assert_true(file.get_block(-1) is None)
  assert_true(file.get_block(1) is None)
  let block = file.get_block(0).unwrap()
  assert_true(block.get_category("_atom_site") is Some(_))
  assert_true(block.get_category("atom_site") is Some(_))
  assert_true(block.get_category("missing") is None)
  let atom_site = block.get_category("_atom_site").unwrap()
  assert_eq(atom_site.name, "_atom_site")
  assert_eq(atom_site.row_count, 8)
}

///|
test "binary cif queries columns by short and qualified names" {
  let file = @src.binary_cif_parse(@src.binary_cif_sample_bytes())
  let atom_site = file.data_blocks[0].get_category("atom_site").unwrap()
  assert_true(atom_site.get_column("label_atom_id") is Some(_))
  assert_true(atom_site.get_column("_atom_site.label_atom_id") is Some(_))
  assert_true(atom_site.get_column("missing") is None)
}

///|
test "binary cif reports decoded column kinds and lengths" {
  let atom_site = @src.binary_cif_parse(@src.binary_cif_sample_bytes())
    .data_blocks[0]
    .get_category("atom_site")
    .unwrap()
  assert_true(
    atom_site.get_column("id").unwrap().kind() is
    @src.BinaryCifColumnKind::BinaryCifInteger,
  )
  assert_true(
    atom_site.get_column("Cartn_x").unwrap().kind() is
    @src.BinaryCifColumnKind::BinaryCifFloat,
  )
  assert_true(
    atom_site.get_column("label_atom_id").unwrap().kind() is
    @src.BinaryCifColumnKind::BinaryCifText,
  )
  assert_eq(atom_site.get_column("Cartn_x").unwrap().length(), 8)
}

///|
test "binary cif reads integer columns" {
  let atom_site = @src.binary_cif_parse(@src.binary_cif_sample_bytes())
    .data_blocks[0]
    .get_category("atom_site")
    .unwrap()
  let ids = atom_site.get_column("id").unwrap()
  assert_eq(ids.int_at(0), Some(1))
  assert_eq(ids.int_at(7), Some(8))
  assert_eq(ids.double_at(3), Some(4.0))
}

///|
test "binary cif reverses fixed-point and delta pipelines" {
  let atom_site = @src.binary_cif_parse(@src.binary_cif_sample_bytes())
    .data_blocks[0]
    .get_category("atom_site")
    .unwrap()
  let xs = atom_site.get_column("Cartn_x").unwrap()
  assert_true((xs.double_at(0).unwrap() - 1.0).abs() < 1.0e-12)
  assert_true((xs.double_at(1).unwrap() - 2.2).abs() < 1.0e-12)
  assert_true((xs.double_at(7).unwrap() - 9.4).abs() < 1.0e-12)
}

///|
test "binary cif reverses string array pipelines" {
  let atom_site = @src.binary_cif_parse(@src.binary_cif_sample_bytes())
    .data_blocks[0]
    .get_category("atom_site")
    .unwrap()
  let components = atom_site.get_column("label_comp_id").unwrap()
  assert_eq(components.string_at(0), Some("GLY"))
  assert_eq(components.string_at(3), Some("ALA"))
  assert_eq(components.string_at(5), Some("HOH"))
  assert_eq(components.string_at(7), Some("SER"))
}

///|
test "binary cif exposes present not-present and unknown masks" {
  let atom_site = @src.binary_cif_parse(@src.binary_cif_sample_bytes())
    .data_blocks[0]
    .get_category("atom_site")
    .unwrap()
  let altloc = atom_site.get_column("label_alt_id").unwrap()
  assert_true(
    altloc.mask_at(0) is @src.BinaryCifMask::BinaryCifNotPresent,
  )
  assert_true(altloc.mask_at(3) is @src.BinaryCifMask::BinaryCifPresent)
  assert_true(altloc.mask_at(6) is @src.BinaryCifMask::BinaryCifUnknown)
  assert_true(altloc.string_at(0) is None)
  assert_eq(altloc.string_at(3), Some("A"))
  assert_true(altloc.string_at(6) is None)
  assert_eq(altloc.cif_token_at(0), ".")
  assert_eq(altloc.cif_token_at(3), "A")
  assert_eq(altloc.cif_token_at(6), "?")
}

///|
test "binary cif column access validates type and bounds" {
  let atom_site = @src.binary_cif_parse(@src.binary_cif_sample_bytes())
    .data_blocks[0]
    .get_category("atom_site")
    .unwrap()
  let names = atom_site.get_column("label_atom_id").unwrap()
  let wrong_type = try {
    ignore(names.int_at(0))
    false
  } catch {
    BinaryCifError(_) => true
  }
  let out_of_bounds = try {
    ignore(names.string_at(8))
    false
  } catch {
    BinaryCifError(_) => true
  }
  assert_true(wrong_type)
  assert_true(out_of_bounds)
}

///|
test "binary cif summarizes document contents" {
  let file = @src.binary_cif_parse(@src.binary_cif_sample_bytes())
  assert_eq(
    file.summary(),
    "BinaryCIF(version=0.3.0, blocks=1, categories=2, rows=9)",
  )
}

///|
test "binary cif builds pdb structure hierarchy" {
  let structure = @src.binary_cif_parse_structure(@src.binary_cif_sample_bytes())
  assert_eq(structure.get_id(), "BCIF")
  assert_eq(structure.get_num_models(), 2)
  assert_eq(structure.get_num_chains(), 3)
  assert_eq(structure.get_num_residues(), 4)
  assert_eq(structure.get_num_atoms(), 8)
}

///|
test "binary cif preserves model and chain grouping" {
  let structure = @src.binary_cif_parse_structure(@src.binary_cif_sample_bytes())
  let models = structure.get_models()
  assert_eq(models[0].get_id(), 0)
  assert_eq(models[1].get_id(), 1)
  assert_eq(models[0].get_num_chains(), 2)
  assert_eq(models[1].get_num_chains(), 1)
  assert_eq(models[0].get_chains()[0].get_id(), 'A')
  assert_eq(models[0].get_chains()[1].get_id(), 'B')
}

///|
test "binary cif preserves residues insertion codes and water" {
  let structure = @src.binary_cif_parse_structure(@src.binary_cif_sample_bytes())
  let residues = structure.get_residues()
  assert_eq(residues[0].resname, "GLY")
  assert_eq(residues[0].resseq, 1)
  assert_false(residues[0].is_het())
  assert_eq(residues[1].resname, "ALA")
  assert_eq(residues[1].icode, 'A')
  assert_eq(residues[2].resname, "HOH")
  assert_true(residues[2].is_het())
  assert_eq(residues[2].chainid, 'B')
  assert_eq(residues[3].resname, "SER")
}

///|
test "binary cif preserves atom annotations and coordinates" {
  let structure = @src.binary_cif_parse_structure(@src.binary_cif_sample_bytes())
  let atoms = structure.get_atoms()
  assert_eq(atoms[0].name, "N")
  assert_eq(atoms[0].element, "N")
  let first = atoms[0].get_coord()
  assert_true((first.x - 1.0).abs() < 1.0e-12)
  assert_true((first.y - 2.0).abs() < 1.0e-12)
  assert_true((first.z + 1.0).abs() < 1.0e-12)
  assert_eq(atoms[3].altloc, 'A')
  assert_true((atoms[3].occupancy - 0.5).abs() < 1.0e-12)
  assert_true((atoms[3].bfactor - 13.0).abs() < 1.0e-12)
}

///|
test "binary cif structure conversion supports id override" {
  let file = @src.binary_cif_parse(@src.binary_cif_sample_bytes())
  let structure = @src.binary_cif_to_structure(file, structure_id="custom")
  assert_eq(structure.get_id(), "custom")
  assert_eq(structure.get_num_atoms(), 8)
}

///|
test "binary cif rejects empty gzip and non-byte inputs" {
  let empty = try {
    ignore(@src.binary_cif_parse([]))
    false
  } catch {
    BinaryCifError(_) => true
  }
  let gzip = try {
    ignore(@src.binary_cif_parse([0x1F, 0x8B, 0x08]))
    false
  } catch {
    BinaryCifError(_) => true
  }
  let non_byte = try {
    ignore(@src.binary_cif_parse([300]))
    false
  } catch {
    BinaryCifError(_) => true
  }
  assert_true(empty)
  assert_true(gzip)
  assert_true(non_byte)
}

///|
test "binary cif rejects truncated and trailing MessagePack data" {
  let bytes = @src.binary_cif_sample_bytes()
  let truncated : Array[Int] = []
  for i in 0..<(bytes.length() - 1) {
    truncated.push(bytes[i])
  }
  let truncated_raised = try {
    ignore(@src.binary_cif_parse(truncated))
    false
  } catch {
    BinaryCifError(_) => true
  }
  bytes.push(0)
  let trailing_raised = try {
    ignore(@src.binary_cif_parse(bytes))
    false
  } catch {
    BinaryCifError(_) => true
  }
  assert_true(truncated_raised)
  assert_true(trailing_raised)
}

///|
test "binary cif rejects unsupported MessagePack markers and roots" {
  let marker = try {
    ignore(@src.binary_cif_parse([0xC1]))
    false
  } catch {
    BinaryCifError(_) => true
  }
  let root = try {
    ignore(@src.binary_cif_parse([0x90]))
    false
  } catch {
    BinaryCifError(_) => true
  }
  assert_true(marker)
  assert_true(root)
}
