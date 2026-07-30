///|
/// Test file for seq_location module.

test "exact_position" {
  let pos = @src.exact_position(42)
  assert_eq!(pos.get(), 42)
  assert_eq!(pos.to_string(), "42")
  assert_eq!(pos.is_exact(), true)
}

test "before_position" {
  let pos = @src.before_position(50)
  assert_eq!(pos.get(), 50)
  assert_eq!(pos.to_string(), "<50")
  assert_eq!(pos.is_exact(), false)
}

test "after_position" {
  let pos = @src.after_position(100)
  assert_eq!(pos.get(), 100)
  assert_eq!(pos.to_string(), ">100")
  assert_eq!(pos.is_exact(), false)
}

test "one_of_position" {
  let opts : Array[@src.Pos] = Array::new()
  opts.push(@src.exact_position(5))
  opts.push(@src.exact_position(7))
  opts.push(@src.exact_position(9))
  let pos = @src.one_of_position(7, opts)
  assert_eq!(pos.get(), 7)
  assert_eq!(pos.is_exact(), false)
  let s = pos.to_string()
  assert_eq!(s.starts_with("{"), true)
  assert_eq!(s.contains("5"), true)
  assert_eq!(s.contains("7"), true)
  assert_eq!(s.contains("9"), true)
}

test "within_position" {
  let pos = @src.within_position(50, @src.exact_position(40), @src.exact_position(60))
  assert_eq!(pos.get(), 50)
  assert_eq!(pos.is_exact(), false)
  let s = pos.to_string()
  assert_eq!(s.starts_with("("), true)
  assert_eq!(s.contains("40"), true)
  assert_eq!(s.contains("60"), true)
}

test "simple_location_basic" {
  let loc = @src.simple_location(@src.exact_position(10), @src.exact_position(50), strand="+")
  assert_eq!(loc.start(), 10)
  assert_eq!(loc.end(), 50)
  assert_eq!(loc.strand(), "+")
  assert_eq!(loc.len(), 40)
  assert_eq!(loc.is_compound(), false)
  assert_eq!(loc.to_string(), "10..50")
}

test "simple_location_reverse_strand" {
  let loc = @src.simple_location(@src.exact_position(100), @src.exact_position(200), strand="-")
  assert_eq!(loc.strand(), "-")
  assert_eq!(loc.len(), 100)
}

test "simple_location_overlaps" {
  let loc1 = @src.simple_location(@src.exact_position(10), @src.exact_position(50))
  let loc2 = @src.simple_location(@src.exact_position(30), @src.exact_position(70))
  let loc3 = @src.simple_location(@src.exact_position(100), @src.exact_position(200))
  assert_eq!(loc1.overlaps(loc2), true)
  assert_eq!(loc1.overlaps(loc3), false)
}

test "simple_location_contains" {
  let loc = @src.simple_location(@src.exact_position(10), @src.exact_position(50))
  assert_eq!(loc.contains(10), true)
  assert_eq!(loc.contains(25), true)
  assert_eq!(loc.contains(49), true)
  assert_eq!(loc.contains(5), false)
  assert_eq!(loc.contains(50), false)
  assert_eq!(loc.contains(100), false)
}

test "compound_location_basic" {
  let exon1 = @src.simple_location(@src.exact_position(0), @src.exact_position(100))
  let exon2 = @src.simple_location(@src.exact_position(200), @src.exact_position(350))
  let locs : Array[@src.Loc] = Array::new()
  locs.push(exon1)
  locs.push(exon2)
  let cl = @src.compound_location(locs, strand="+")
  assert_eq!(cl.is_compound(), true)
  assert_eq!(cl.parts().length(), 2)
  assert_eq!(cl.len(), 250)
}

test "compound_from_simple" {
  let starts : Array[Int] = [0, 200, 500]
  let ends : Array[Int] = [100, 350, 600]
  let cl = @src.compound_from_simple(starts, ends, strand="+")
  assert_eq!(cl.is_compound(), true)
  assert_eq!(cl.parts().length(), 3)
  assert_eq!(cl.len(), 350)  // 100 + 150 + 100
}

test "compound_location_overlaps" {
  let exon1 = @src.simple_location(@src.exact_position(0), @src.exact_position(100))
  let exon2 = @src.simple_location(@src.exact_position(200), @src.exact_position(350))
  let locs : Array[@src.Loc] = Array::new()
  locs.push(exon1)
  locs.push(exon2)
  let cl = @src.compound_location(locs)

  let overlapping = @src.simple_location(@src.exact_position(50), @src.exact_position(150))
  let not_overlapping = @src.simple_location(@src.exact_position(400), @src.exact_position(500))
  assert_eq!(cl.overlaps(overlapping), true)
  assert_eq!(cl.overlaps(not_overlapping), false)
}

test "compound_location_contains" {
  let exon1 = @src.simple_location(@src.exact_position(0), @src.exact_position(100))
  let exon2 = @src.simple_location(@src.exact_position(200), @src.exact_position(350))
  let locs : Array[@src.Loc] = Array::new()
  locs.push(exon1)
  locs.push(exon2)
  let cl = @src.compound_location(locs)
  assert_eq!(cl.contains(50), true)
  assert_eq!(cl.contains(250), true)
  assert_eq!(cl.contains(150), false)
  assert_eq!(cl.contains(400), false)
}

test "location_genbank_format" {
  let loc = @src.simple_location(@src.exact_position(9), @src.exact_position(50))
  let gb = @src.location_to_genbank(loc)
  assert_eq!(gb.contains("10"), true)
  assert_eq!(gb.contains("50"), true)
}

test "from_genbank_coords" {
  let loc = @src.from_genbank_coords(10, 50, strand="+")
  assert_eq!(loc.start(), 9)
  assert_eq!(loc.end(), 50)
  assert_eq!(loc.strand(), "+")
  assert_eq!(loc.len(), 41)
}

test "parse_genbank_location_simple" {
  let loc = @src.parse_genbank_location("10..50")
  assert_eq!(loc.start(), 9)  // 0-based
  assert_eq!(loc.end(), 50)
  assert_eq!(loc.len(), 41)
}

test "seq_location_sample" {
  let samples = @src.seq_location_sample()
  assert_eq!(samples.length(), 5)
  // First is simple forward
  assert_eq!(samples[0].is_compound(), false)
  // Third is compound
  assert_eq!(samples[2].is_compound(), true)
}
