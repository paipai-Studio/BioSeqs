///|
test "cooler basic construction" {
  let c = @src.Cooler::new(1000)
  ignore(c.add_chrom("chr1", 2500))
  ignore(c.add_chrom("chr2", 1200))
  assert_eq(c.bin_size(), 1000)
  assert_eq(c.bin_count(), 3 + 2)
  assert_eq(c.chrom_names(), ["chr1", "chr2"])
  assert_eq(c.chrom_length("chr1"), Some(2500))
  assert_true(c.chrom_length("chr3") is None)
}

///|
test "cooler bin index" {
  let c = @src.Cooler::new(1000)
  ignore(c.add_chrom("chr1", 2500))
  ignore(c.add_chrom("chr2", 1200))
  assert_eq(c.bin_index("chr1", 0), Some(0))
  assert_eq(c.bin_index("chr1", 999), Some(0))
  assert_eq(c.bin_index("chr1", 1000), Some(1))
  assert_eq(c.bin_index("chr1", 2499), Some(2))
  assert_eq(c.bin_index("chr2", 0), Some(3))
  assert_eq(c.bin_index("chr2", 1199), Some(4))
  assert_true(c.bin_index("chr2", 1200) is None)
  assert_true(c.bin_index("chr3", 0) is None)
  assert_true(c.bin_index("chr1", -1) is None)
}

///|
test "cooler bin location roundtrip" {
  let c = @src.Cooler::new(1000)
  ignore(c.add_chrom("chr1", 2500))
  ignore(c.add_chrom("chr2", 1200))
  assert_eq(c.bin_location(0), Some(("chr1", 0)))
  assert_eq(c.bin_location(2), Some(("chr1", 2000)))
  assert_eq(c.bin_location(3), Some(("chr2", 0)))
  assert_eq(c.bin_location(4), Some(("chr2", 1000)))
  assert_true(c.bin_location(5) is None)
}

///|
test "cooler add contact and accumulate" {
  let c = @src.Cooler::new(1000)
  ignore(c.add_chrom("chr1", 10000))
  ignore(c.add_contact("chr1", 1500, 2600, 3))
  ignore(c.add_contact("chr1", 2500, 1200, 2))
  assert_eq(c.pixel_count(), 1)
  assert_eq(c.total_counts(), 5)
  assert_eq(c.get(1, 2), Some(5))
  assert_eq(c.get(2, 1), Some(5))
  assert_true(c.get(0, 1) is None)
}

///|
test "cooler add pixel out of range raises" {
  let c = @src.Cooler::new(1000)
  ignore(c.add_chrom("chr1", 2500))
  let msg = try {
    ignore(c.add_pixel(0, 5, 1))
    "NO_ERROR"
  } catch {
    CoolerError(m) => m
  }
  assert_true(msg.contains("out of range"))
}

///|
test "cooler add contact invalid coordinate raises" {
  let c = @src.Cooler::new(1000)
  ignore(c.add_chrom("chr1", 2500))
  let msg = try {
    ignore(c.add_contact("chr1", 0, 99999, 1))
    "NO_ERROR"
  } catch {
    CoolerError(m) => m
  }
  assert_true(msg.contains("invalid coordinate"))
}

///|
test "cooler pixels sorted" {
  let c = @src.Cooler::new(1000)
  ignore(c.add_chrom("chr1", 5000))
  ignore(c.add_pixel(3, 1, 1))
  ignore(c.add_pixel(0, 2, 1))
  ignore(c.add_pixel(1, 4, 1))
  let px = c.pixels()
  assert_eq(px.length(), 3)
  assert_eq(px[0].bin1(), 0)
  assert_eq(px[0].bin2(), 2)
  assert_eq(px[1].bin1(), 1)
  assert_eq(px[1].bin2(), 3)
  assert_eq(px[2].bin1(), 1)
  assert_eq(px[2].bin2(), 4)
}

///|
test "cooler balance equalizes coverage" {
  // 均匀数据: 平衡后所有 bin 权重应相等
  let c = @src.Cooler::new(1000)
  ignore(c.add_chrom("chr1", 5000))
  for i in 0..<5 {
    for j in i..<5 {
      ignore(c.add_pixel(i, j, 10))
    }
  }
  let w = c.balance()
  assert_eq(w.length(), 5)
  for v in w {
    assert_true((v - 1.0).abs() < 0.1)
  }
}

///|
test "cooler balance zero coverage bin" {
  let c = @src.Cooler::new(1000)
  ignore(c.add_chrom("chr1", 5000))
  ignore(c.add_pixel(0, 1, 5))
  let w = c.balance(iterations=5)
  assert_eq(w[2], 0.0)
  assert_true(w[0] > 0.0)
}

///|
test "cooler coarsen merges counts" {
  let c = @src.Cooler::new(1000)
  ignore(c.add_chrom("chr1", 10000))
  // bin 0..3 的接触粗化后应合并到新 bin 0/1
  ignore(c.add_pixel(0, 0, 4))
  ignore(c.add_pixel(1, 1, 6))
  ignore(c.add_pixel(2, 2, 10))
  ignore(c.add_pixel(3, 3, 14))
  let coarse = c.coarsen(2)
  assert_eq(coarse.bin_size(), 2000)
  assert_eq(coarse.get(0, 0), Some(10))
  assert_eq(coarse.get(1, 1), Some(24))
  assert_eq(coarse.total_counts(), 34)
  assert_eq(coarse.chrom_names(), ["chr1"])
}

///|
test "cooler coarsen factor one copies" {
  let c = @src.Cooler::new(1000)
  ignore(c.add_chrom("chr1", 3000))
  ignore(c.add_pixel(0, 1, 7))
  let copy = c.coarsen(1)
  assert_eq(copy.bin_size(), 1000)
  assert_eq(copy.get(0, 1), Some(7))
}

///|
test "cooler to hicdc contacts" {
  let c = @src.Cooler::new(1000)
  ignore(c.add_chrom("chr1", 3000))
  ignore(c.add_chrom("chr2", 2000))
  ignore(c.add_pixel(0, 1, 5))
  ignore(c.add_pixel(3, 4, 2))
  let contacts = c.to_hicdc_contacts()
  assert_eq(contacts.length(), 2)
  assert_eq(contacts[0].chr, "chr1")
  assert_eq(contacts[0].count, 5)
  assert_eq(contacts[1].chr, "chr2")
}

///|
test "cooler bins as bed" {
  let c = @src.Cooler::new(1000)
  ignore(c.add_chrom("chr1", 2500))
  ignore(c.add_pixel(0, 0, 3))
  let bed = c.bins_as_bed("chr1")
  assert_eq(bed, "chr1\t0\t1000\nchr1\t1000\t2000\nchr1\t2000\t2500\n")
  let bed_nz = c.bins_as_bed("chr1", nonzero_only=true)
  assert_eq(bed_nz, "chr1\t0\t1000\n")
  assert_eq(c.bins_as_bed("chrX"), "")
}

///|
test "cooler json roundtrip" {
  let c = @src.Cooler::new(1000)
  ignore(c.add_chrom("chr1", 3000))
  ignore(c.add_chrom("chr2", 1500))
  ignore(c.add_pixel(0, 1, 5))
  ignore(c.add_pixel(3, 3, 2))
  let json = c.to_json()
  let back = @src.cooler_from_json(json).unwrap()
  assert_eq(back.bin_size(), 1000)
  assert_eq(back.bin_count(), c.bin_count())
  assert_eq(back.get(0, 1), Some(5))
  assert_eq(back.get(3, 3), Some(2))
  assert_eq(back.total_counts(), 7)
  assert_eq(back.chrom_names(), ["chr1", "chr2"])
}

///|
test "cooler json invalid returns none" {
  assert_true(@src.cooler_from_json("not json") is None)
  assert_true(@src.cooler_from_json("{\"bin_size\":0}") is None)
}

///|
test "cooler mcool multiresolution pipeline" {
  // 合成 Hi-C 数据 -> 粗化两级 -> 总计数守恒
  let c = @src.Cooler::new(1000)
  ignore(c.add_chrom("chr1", 20000))
  let rng = @src.Rng::from_int(7)
  for _i in 0..<500 {
    let a = rng.next_int(20)
    let b = rng.next_int(20)
    ignore(c.add_pixel(a, b, 1))
  }
  let total = c.total_counts()
  let mid = c.coarsen(2)
  let coarse = c.coarsen(4)
  assert_eq(mid.bin_size(), 2000)
  assert_eq(coarse.bin_size(), 4000)
  assert_eq(mid.total_counts(), total)
  assert_eq(coarse.total_counts(), total)
}
