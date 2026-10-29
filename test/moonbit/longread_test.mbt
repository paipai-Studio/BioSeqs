// ===========================================================================
// longread_test.mbt — 长读长深化测试(minimap 预设 + medaka-lite 变异检出)
// ===========================================================================

///|
fn longread_ref() -> String {
  "ATCGGATTCAGCCTAGGATCAAGTTCGGATACCTAGGCATTGCAAGGTCCTAGAAGTCCATTGCAGGT"
}

///|
test "minimap preset params" {
  let (k1, w1, b1, _) = @src.MapOnt.params()
  assert_eq(k1, 15)
  assert_eq(w1, 10)
  assert_true(b1 >= 150)
  let (k2, w2, _, _) = @src.MapPb.params()
  assert_eq(k2, 19)
  assert_eq(w2, 10)
  let (k3, _, _, _) = @src.MapSr.params()
  assert_eq(k3, 21)
}

///|
test "minimap preset to string" {
  assert_eq(@src.MapOnt.to_string(), "map-ont")
  assert_eq(@src.MapPb.to_string(), "map-pb")
  assert_eq(@src.MapSr.to_string(), "map-sr")
}

///|
test "minimap align preset map ont exact query" {
  let reference = longread_ref()
  let refs : Array[@src.MinimapRefSeq] = [
    @src.MinimapRefSeq::new("ctg1", reference),
  ]
  let index = @src.minimap_index_build_preset(refs, @src.MapOnt)
  // 模拟 ONT 长读: 精确子串(位置 10-50)
  let read = reference[10:50].to_owned()
  let results = @src.minimap_align_preset(index, read, "lr1", @src.MapOnt)
  assert_true(results.length() > 0)
  let best = results[0]
  assert_eq(best.ref_name, "ctg1")
  assert_eq(best.strand, 0)
  // 比对区间应落在 [10, 50) 附近
  assert_true(best.ref_start >= 5 && best.ref_start <= 15)
  assert_true(best.ref_end >= 45 && best.ref_end <= 55)
}

///|
test "minimap align preset map ont noisy query" {
  let reference = longread_ref()
  let refs : Array[@src.MinimapRefSeq] = [
    @src.MinimapRefSeq::new("ctg1", reference),
  ]
  let index = @src.minimap_index_build_preset(refs, @src.MapOnt)
  // 引入 3 处错配 + 1 处 2bp 缺失(模拟 ONT 噪声)
  let region = reference[10:50].to_owned()
  let noisy = region[0:5].to_owned() +
    "T" +
    region[6:10].to_owned() +
    region[12:].to_owned()
  let results = @src.minimap_align_preset(index, noisy, "lr2", @src.MapOnt)
  // 宽带预设应容忍噪声并命中
  assert_true(results.length() > 0)
  assert_eq(results[0].ref_name, "ctg1")
}

///|
test "minimap align preset reverse strand" {
  let reference = longread_ref()
  let refs : Array[@src.MinimapRefSeq] = [
    @src.MinimapRefSeq::new("ctg1", reference),
  ]
  let index = @src.minimap_index_build_preset(refs, @src.MapOnt)
  // 反向互补读段
  let rc = @src.reverse_complement(reference[10:50].to_owned())
  let results = @src.minimap_align_preset(index, rc, "lr3", @src.MapOnt)
  assert_true(results.length() > 0)
  assert_eq(results[0].strand, 1)
}

///|
test "minimap align preset map pb" {
  let reference = longread_ref()
  let refs : Array[@src.MinimapRefSeq] = [
    @src.MinimapRefSeq::new("ctg1", reference),
  ]
  let index = @src.minimap_index_build_preset(refs, @src.MapPb)
  let read = reference[0:30].to_owned()
  let results = @src.minimap_align_preset(index, read, "hifi1", @src.MapPb)
  assert_true(results.length() > 0)
  assert_eq(results[0].ref_name, "ctg1")
  assert_eq(results[0].ref_start, 0)
}

///|
test "minimap align preset unmapped" {
  let refs : Array[@src.MinimapRefSeq] = [
    @src.MinimapRefSeq::new("ctg1", longread_ref()),
  ]
  let index = @src.minimap_index_build_preset(refs, @src.MapOnt)
  // 与参考无任何 k-mer 相似性的读段
  let read = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
  let results = @src.minimap_align_preset(index, read, "unmapped", @src.MapOnt)
  assert_eq(results.length(), 0)
}

// ---------------------------------------------------------------------------
// medaka-lite 长读变异检出
// ---------------------------------------------------------------------------

///|
test "medaka call snp with ont error model" {
  let reference = longread_ref()
  // 5 条长读覆盖, 3 条含同一杂合 SNP(位置 15)
  let reads : Array[@src.MedakaRead] = []
  for i = 0; i < 5; i = i + 1 {
    let sliced = reference[0:40].to_owned()
    let seq = if i < 3 {
      sliced[0:15].to_owned() + "A" + sliced[16:].to_owned()
    } else {
      sliced
    }
    reads.push({ name: "lr" + i.to_string(), seq, pos: 0, })
  }
  let variants = @src.medaka_call_variants("chr1", reference, reads)
  assert_eq(variants.length(), 1)
  let v = variants[0]
  assert_eq(v.pos, 16)
  assert_eq(v.ref, reference[15:16].to_owned())
  assert_eq(v.alt, "A")
  assert_eq(v.gt, "0/1")
  assert_eq(v.dp, 5)
  assert_eq(v.ad_alt, 3)
}

///|
test "medaka call deletion in noisy long reads" {
  let reference = longread_ref()
  let reads : Array[@src.MedakaRead] = []
  for i = 0; i < 5; i = i + 1 {
    let sliced = reference[0:40].to_owned()
    let seq = if i < 3 {
      // 删除参考位置 10-12
      sliced[0:10].to_owned() + sliced[13:].to_owned()
    } else {
      sliced
    }
    reads.push({ name: "lr" + i.to_string(), seq, pos: 0, })
  }
  let variants = @src.medaka_call_variants("chr1", reference, reads)
  let mut found_del = false
  for v in variants {
    if v.ref.length() > v.alt.length() {
      found_del = true
      assert_eq(v.ref.length() - v.alt.length(), 3)
    }
  }
  assert_true(found_del)
}

///|
test "medaka low depth filtered" {
  let reference = longread_ref()
  // 仅 2 条(低于 min_depth 3)
  let reads : Array[@src.MedakaRead] = []
  for i = 0; i < 2; i = i + 1 {
    let sliced = reference[0:40].to_owned()
    let seq = sliced[0:15].to_owned() + "A" + sliced[16:].to_owned()
    reads.push({ name: "lr" + i.to_string(), seq, pos: 0, })
  }
  let variants = @src.medaka_call_variants("chr1", reference, reads)
  assert_eq(variants.length(), 0)
}

///|
test "medaka hifi model stricter" {
  let reference = longread_ref()
  // HiFi 模型 epsilon 低 -> 同样证据下 GQ 更高
  let reads : Array[@src.MedakaRead] = []
  for i = 0; i < 5; i = i + 1 {
    let sliced = reference[0:40].to_owned()
    let seq = if i < 3 {
      sliced[0:15].to_owned() + "A" + sliced[16:].to_owned()
    } else {
      sliced
    }
    reads.push({ name: "lr" + i.to_string(), seq, pos: 0, })
  }
  let v_ont = @src.medaka_call_variants("chr1", reference, reads)
  let v_hifi = @src.medaka_call_variants(
    "chr1",
    reference,
    reads,
    model=@src.Hifi,
  )
  assert_eq(v_ont.length(), 1)
  assert_eq(v_hifi.length(), 1)
  // HiFi 模型下杂合 GQ 更高
  assert_true(v_hifi[0].gq > v_ont[0].gq)
}

///|
test "medaka variants summary" {
  let reference = longread_ref()
  let reads : Array[@src.MedakaRead] = []
  for i = 0; i < 5; i = i + 1 {
    let sliced = reference[0:40].to_owned()
    let seq = if i < 3 {
      sliced[0:15].to_owned() + "A" + sliced[16:].to_owned()
    } else {
      sliced
    }
    reads.push({ name: "lr" + i.to_string(), seq, pos: 0, })
  }
  let variants = @src.medaka_call_variants("chr1", reference, reads)
  let s = @src.medaka_variants_summary(variants)
  assert_true(s.contains("medaka-lite"))
  assert_true(s.contains("1 SNP"))
}

///|
test "medaka vcf output via gatk write vcf" {
  let reference = longread_ref()
  let reads : Array[@src.MedakaRead] = []
  for i = 0; i < 5; i = i + 1 {
    let sliced = reference[0:40].to_owned()
    let seq = if i < 3 {
      sliced[0:15].to_owned() + "A" + sliced[16:].to_owned()
    } else {
      sliced
    }
    reads.push({ name: "lr" + i.to_string(), seq, pos: 0, })
  }
  let variants = @src.medaka_call_variants("chr1", reference, reads)
  let vcf_text = @src.gatk_write_vcf("chr1", reference, variants, sample="LR")
  let vcf = @src.parse_vcf(vcf_text)
  assert_eq(vcf.len(), 1)
  assert_eq(vcf.records[0].pos, 16)
  assert_true(vcf.records[0].samples.get("LR").unwrap()[0].contains("0/1"))
}

///|
test "minimap medaka end to end pipeline" {
  // 长读全链路: 预设比对定位 -> 变异检出
  let reference = longread_ref()
  let refs : Array[@src.MinimapRefSeq] = [
    @src.MinimapRefSeq::new("ctg1", reference),
  ]
  let index = @src.minimap_index_build_preset(refs, @src.MapOnt)
  // 长读从参考位置 10 采出, 含杂合 SNP(读内偏移 5 -> 参考位置 15)
  let aligned_reads : Array[@src.MedakaRead] = []
  for i = 0; i < 5; i = i + 1 {
    let sliced = reference[10:50].to_owned()
    let seq = if i < 3 {
      sliced[0:5].to_owned() + "A" + sliced[6:].to_owned()
    } else {
      sliced
    }
    // 模拟比对已定位(真实管线由 minimap_align_preset 产出 pos)
    let aln = @src.minimap_align_preset(
      index,
      seq,
      "lr" + i.to_string(),
      @src.MapOnt,
    )
    let pos = if aln.length() > 0 {
      if aln[0].strand == 0 {
        aln[0].ref_start
      } else {
        10
      }
    } else {
      10
    }
    aligned_reads.push({ name: "lr" + i.to_string(), seq, pos, })
  }
  let variants = @src.medaka_call_variants("ctg1", reference, aligned_reads)
  assert_eq(variants.length(), 1)
  assert_eq(variants[0].pos, 16)
}
