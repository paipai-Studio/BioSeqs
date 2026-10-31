///|
/// 分子钟与祖先序列重建演示(TreeTime 风格):
/// 采样日期 + Newick 树 -> root-to-tip 回归 -> 时间树 -> 祖先状态重建.
fn main raise {
  println("=== 分子钟标定与祖先序列重建 (TreeTime-lite) ===")

  // 一棵 4 叶超度量树的 Newick(单位: 每位点替换数).
  let newick = "(((A:0.02,B:0.02):0.03,C:0.05):0.05,D:0.15);"
  let tree = @src.skbio_tree_read(newick)
  println("\nNewick: " + newick)

  // 采样日期(年).
  let dates : Map[String, Double] = Map([
    ("A", 2020.0),
    ("B", 2021.0),
    ("C", 2022.0),
    ("D", 2023.0),
  ])

  // 1. root-to-tip 回归.
  println("\n1. root-to-tip 回归")
  let fit = @src.phylo_clock_root_to_tip(tree, dates)
  println("  估计替换速率 = \{fit.rate} 替换/位点/年")
  println("  根节点日期   = \{fit.root_date}")
  println("  R²           = \{fit.r_squared}")
  println("  逐叶残差:")
  for tip in fit.tips {
    println(
      "  \{tip.name}: 日期=\{tip.date}, root-to-tip=\{tip.root_to_tip}, 残差=\{tip.residual}",
    )
  }

  // 2. 时间树.
  println("\n2. 时间树节点日期")
  let dated = @src.phylo_clock_dated_tree(tree, fit)
  print_dated_tree(dated, 0)

  // 3. 祖先状态重建(Jukes-Cantor + Felsenstein 剪枝).
  println("\n3. 祖先状态重建")
  let sequences : Map[String, String] = Map([
    ("A", "ACGTACGTAC"),
    ("B", "ACGTACGTAC"),
    ("C", "ACGTACGTAG"),
    ("D", "ACGTACGTAG"),
  ])
  let ancestral = @src.phylo_clock_ancestral(tree, sequences)
  println("  整树对数似然 = \{ancestral.log_likelihood}")
  for taxon in ancestral.taxa {
    let lowest = lowest_probability(taxon.sites)
    println(
      "  \{taxon.name}: \{taxon.sequence} (最低位点置信度 \{lowest})",
    )
  }

  // 4. JSON 输出(衔接 Web/报告架构).
  println("\n4. JSON 输出(节选)")
  let json = @src.phylo_clock_ancestral_json(ancestral)
  let lines = json.split("\n").to_array()
  for i = 0; i < lines.length() && i < 6; i = i + 1 {
    println("  " + lines[i].to_owned())
  }
}

///|
fn print_dated_tree(node : @src.DatedNode, depth : Int32) -> Unit {
  let indent = "  ".repeat(depth)
  let name = match node.name {
    Some(n) => n
    None => "(内部节点)"
  }
  println("  \{indent}\{name}: \{node.time}")
  for child in node.children {
    print_dated_tree(child, depth + 1)
  }
}

///|
fn lowest_probability(sites : Array[@src.AncestralSite]) -> Double {
  let mut lowest = 1.0
  for site in sites {
    if site.probability < lowest {
      lowest = site.probability
    }
  }
  lowest
}
