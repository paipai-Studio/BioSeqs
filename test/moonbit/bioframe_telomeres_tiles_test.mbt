///|
/// Tests for bioframe_telomeres_tiles.mbt — telomere intervals + tiled regions.

// ---- bf_telomeres ---------------------------------------------------------

test "bf_telomeres: empty chromsizes → empty" {
  let out = @src.bf_telomeres([])
  assert_eq(out.length(), 0)
}

test "bf_telomeres: single chrom pad=10000 → p and q arms" {
  let cs = [("chr1", 248956422)]
  let out = @src.bf_telomeres(cs, pad=10000)
  assert_eq(out.length(), 2)
  assert_eq(out[0].chrom, "chr1")
  assert_eq(out[0].start, 0)
  assert_eq(out[0].end, 10000)
  assert_eq(out[1].chrom, "chr1")
  assert_eq(out[1].start, 248956422 - 10000)
  assert_eq(out[1].end, 248956422)
}

test "bf_telomeres: tiny chromosome → skipped" {
  let cs = [("chr1", 5000)]  // 5000 <= 2*10000 → skip
  let out = @src.bf_telomeres(cs, pad=10000)
  assert_eq(out.length(), 0)
}

test "bf_telomeres: custom pad" {
  let cs = [("chr1", 1000)]
  let out = @src.bf_telomeres(cs, pad=200)
  assert_eq(out.length(), 2)
  assert_eq(out[0].end, 200)
  assert_eq(out[1].start, 800)
}

// ---- bf_telomeres_labelled -----------------------------------------------

test "bf_telomeres_labelled: side char 'p' and 'q'" {
  let cs = [("chr1", 1000)]
  let out = @src.bf_telomeres_labelled(cs, pad=100)
  assert_eq(out.length(), 2)
  assert_eq(out[0].side, 'p')
  assert_eq(out[1].side, 'q')
}

// ---- bf_tiled_regions ----------------------------------------------------

// NOTE on source logic:
//   while cur + ts <= length  → emit full tile (is_last=false)
//   if cur < length           → emit partial tile (is_last=true)
// So the LAST FULL tile is still is_last=false; only a final clipping
// gets is_last=true.

test "bf_tiled_regions: length divisible by tile → all full, last=false" {
  let cs = [("chr1", 100)]
  let out = @src.bf_tiled_regions(cs, tile_size=50, stride=50)
  // [0,50) full, [50,100) full, no partial (cur=100 not < length)
  assert_eq(out.length(), 2)
  assert_eq(out[0].start, 0); assert_eq(out[0].end, 50)
  assert_eq(out[1].start, 50); assert_eq(out[1].end, 100)
  // Source: full tiles all is_last=false
  assert_true(!out[0].is_last)
  assert_true(!out[1].is_last)
}

test "bf_tiled_regions: remainder creates partial is_last=true" {
  let cs = [("chr1", 100)]
  let out = @src.bf_tiled_regions(cs, tile_size=30, stride=30)
  // cur=0  : 0+30<=100 → [0,30) full, cur=30
  // cur=30 : 30+30<=100 → [30,60) full, cur=60
  // cur=60 : 60+30<=100 → [60,90) full, cur=90
  // cur=90 : 90+30=120>100 → while exits; cur<100 → [90,100) partial
  assert_eq(out.length(), 4)
  assert_eq(out[0].end, 30)
  assert_eq(out[3].start, 90); assert_eq(out[3].end, 100)
  assert_true(out[3].is_last)
  assert_true(!out[0].is_last)
  assert_true(!out[1].is_last)
  assert_true(!out[2].is_last)
}

test "bf_tiled_regions: overlapping tiles (stride < tile_size)" {
  let cs = [("chr1", 100)]
  let out = @src.bf_tiled_regions(cs, tile_size=40, stride=20)
  // cur=0 : 0+40<=100 → [0,40) full, cur=20
  // cur=20: 20+40<=100 → [20,60) full, cur=40
  // cur=40: 40+40<=100 → [40,80) full, cur=60
  // cur=60: 60+40<=100 → [60,100) full, cur=80
  // cur=80: 80+40=120>100 → while exits; cur=80<100 → [80,100) partial
  assert_eq(out.length(), 5)
  assert_eq(out[0].start, 0);  assert_eq(out[0].end, 40)
  assert_eq(out[1].start, 20); assert_eq(out[1].end, 60)
  assert_eq(out[2].start, 40); assert_eq(out[2].end, 80)
  assert_eq(out[3].start, 60); assert_eq(out[3].end, 100)
  assert_eq(out[4].start, 80); assert_eq(out[4].end, 100)
  assert_true(out[4].is_last)
}

test "bf_tiled_regions: zero-length chrom → no tiles" {
  let cs = [("chr1", 0)]
  let out = @src.bf_tiled_regions(cs, tile_size=100, stride=100)
  assert_eq(out.length(), 0)
}

// ---- bf_tiled_as_bed -----------------------------------------------------

// NOTE: source hardcodes tile_size=1000000 / stride=1000000 internally
// (ignores the parameters passed in — a known limitation).
test "bf_tiled_as_bed: uses default 1_000_000 tile_size" {
  let cs = [("chr1", 100)]
  // tile_size=1_000_000 by default → only partial [0,100)
  let out = @src.bf_tiled_as_bed(cs)
  assert_eq(out.length(), 1)
  assert_eq(out[0].chrom, "chr1")
  assert_eq(out[0].start, 0)
  assert_eq(out[0].end, 100)
}
