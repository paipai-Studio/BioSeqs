///|
/// de novo motif 发现演示: 植入 motif -> EM/Gibbs 恢复 -> 与 motif 栈闭环.
fn main raise {
  println("=== de novo motif 发现 (EM / Gibbs) ===")

  // 构造含植入 motif 的序列集
  let planted = "ACGTTGCA"
  let rng = @src.Rng::from_int(12345)
  let bg : Array[Char] = ['A', 'C', 'G', 'T']
  let seqs : Array[String] = Array(capacity=0)
  for i = 0; i < 10; i = i + 1 {
    let sb = StringBuilder()
    let left_len = 8 + rng.next_int(5)
    for j = 0; j < left_len; j = j + 1 {
      sb.write_char(bg[rng.next_int(4)])
    }
    sb.write_string(planted)
    let right_len = 8 + rng.next_int(5)
    for j = 0; j < right_len; j = j + 1 {
      sb.write_char(bg[rng.next_int(4)])
    }
    ignore(i)
    seqs.push(sb.to_string())
  }
  println("植入 motif: \{planted}, 序列 \{seqs.length()} 条")

  // 1. EM 发现
  println("\n1. EM 发现(MEME 风格多起点)")
  let em = @src.motif_em(seqs, width=8, iterations=60)
  println("  " + em.summary())
  println("  位点指认:")
  for i = 0; i < em.site_positions.length(); i = i + 1 {
    let p = em.site_positions[i]
    println("    seq\{i}: pos=\{p.0} strand=\{p.1} site=\{em.sites[i]}")
  }

  // 2. Gibbs 发现
  println("\n2. Gibbs 采样发现(AlignACE 风格, 确定性 RNG)")
  let gibbs = @src.motif_gibbs(seqs, width=8, restarts=3, sweeps=40, seed=42)
  println("  " + gibbs.summary())

  // 3. 双链一致性
  println("\n3. 双链一致性检查")
  let cons_em = em.consensus()
  let cons_gibbs = gibbs.consensus()
  let rc_gibbs = @src.reverse_complement(cons_gibbs)
  println(
    "  EM consensus=\{cons_em}, Gibbs consensus=\{cons_gibbs} (rc=\{rc_gibbs})",
  )
  let same = cons_em == cons_gibbs || cons_em == rc_gibbs
  println("  双算法一致性: \{same}")

  // 4. 与 motif 栈闭环
  println("\n4. 与 motif 栈闭环")
  let m = em.to_motif()
  println(
    "  to_motif: name=\{m.name} pwm_rows=\{m.pwm.length()} sites=\{m.sites.length()}",
  )
  println(
    "  consensus(库内 get_consensus 风格): \{@src.get_consensus(m.pwm[0], "ACGT")}",
  )

  // 5. 随机背景对照
  println("\n5. 随机背景对照(应低 IC)")
  let rng2 = @src.Rng::from_int(99)
  let null_seqs : Array[String] = Array(capacity=0)
  for i = 0; i < 10; i = i + 1 {
    let sb = StringBuilder()
    for j = 0; j < 30; j = j + 1 {
      sb.write_char(bg[rng2.next_int(4)])
    }
    ignore(i)
    null_seqs.push(sb.to_string())
  }
  let null_em = @src.motif_em(null_seqs, width=8, iterations=30)
  println(
    "  null IC=\{null_em.information_content} (signal IC=\{em.information_content})",
  )
  println(
    "  信噪分离: \{em.information_content > null_em.information_content}",
  )

  println(
    "\n完成: 多起点 EM + Gibbs 双算法 de novo 发现本地完成.",
  )
}
