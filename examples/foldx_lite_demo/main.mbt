///|
/// foldx_lite_demo — FoldX 风格蛋白稳定性 ΔΔG 与饱和突变扫描演示:
/// 合成 α-螺旋链 -> 单点 ΔΔG -> 饱和突变扫描(最稳/最不稳/熵) ->
/// 突变文库排名 -> TSV/JSON.

///|
fn fxd_chain() -> @src.Chain {
  let seq = [
    "ALA", "GLY", "SER", "LEU", "ILE", "VAL", "PHE", "MET", "LEU", "SER", "ILE",
    "VAL",
  ]
  let residues : Array[@src.Residue] = Array(capacity=0)
  for r = 0; r < seq.length(); r = r + 1 {
    let angle = r.to_double() * 1.05
    let x = 4.0 * @math.cos(angle)
    let y = 4.0 * @math.sin(angle)
    let z = r.to_double() * 2.1
    let ca = @src.Atom::new(
      name="CA",
      coord=@src.Vector3::new(x, y, z),
      resname=seq[r],
      chainid='A',
      resseq=r + 1,
      element="C",
    )
    let res = @src.Residue::new(resname=seq[r], chainid='A', resseq=r + 1, atoms=[
      ca,
    ])
    residues.push(res)
  }
  @src.Chain::new(id='A', residues~)
}

///|
fn main {
  let chain = fxd_chain()
  println("1. 合成链: 12 残基 α-螺旋(Cα 模型, DOPE 统计势)")
  // 2. 单点 ΔΔG
  let e_wt = @src.foldx_energy(chain).unwrap()
  let ddg_gly = @src.foldx_ddg(chain, 7, "GLY").unwrap()
  let ddg_trp = @src.foldx_ddg(chain, 7, "TRP").unwrap()
  println("2. 单点 ΔΔG(pos7 PHE):")
  println("  wt 能量 = \{e_wt}")
  println("  PHE->GLY ΔΔG = \{ddg_gly} (小侧链)")
  println("  PHE->TRP ΔΔG = \{ddg_trp} (大侧链)")
  // 3. 饱和突变扫描(核心位点 4/7/10)
  let positions : Array[Int32] = [4, 7, 10]
  let sites = @src.foldx_saturation_scan(chain, positions=Some(positions))
  println("3. 饱和突变扫描(位点 4/7/10):")
  for s in sites {
    println(
      "  pos \{s.pos} \{s.wt}: best=\{s.best.mt}(\{s.best.ddg}), worst=\{s.worst.mt}(\{s.worst.ddg}), entropy=\{s.entropy}",
    )
  }
  // 4. 汇总
  let summary = @src.foldx_scan_summary(sites).unwrap()
  println("4. 位点数=\{summary.n_sites}, 平均熵=\{summary.mean_entropy}")
  println(@src.foldx_summary_line(summary))
  // 5. 文库排名
  let lib = @src.foldx_mutant_library(sites)
  println("5. 突变文库(\{lib.length()} 条, 最稳前 3):")
  let mut shown = 0
  for m in lib {
    if shown < 3 {
      println(
        "  \{m.label} | ΔΔG=\{m.total_ddg} | fitness=\{m.mutations[0].fitness}",
      )
      shown = shown + 1
    }
  }
  // 6. TSV/JSON
  let tsv = @src.foldx_scan_tsv(sites)
  let json = @src.foldx_scan_json(sites)
  println("6. TSV \{tsv.length()} 字节, JSON \{json.length()} 字节")
  println(
    "演示完成: 链能量 -> 单点 ΔΔG -> 饱和扫描 -> 文库排名闭环",
  )
}
