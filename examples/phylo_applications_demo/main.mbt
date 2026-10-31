///|
/// Bio.Phylo.Applications Demo.
///
/// Demonstrates the command-line wrappers for the FastTree, PhyML and
/// RAxML phylogenetics tools. The wrappers are pure parameter collectors:
/// they serialize to the exact command line that would be executed.

///|
fn main {
  println("=== 1. FastTree: approximately maximum-likelihood trees ===")
  let fasttree = @src.FastTreeCommandline::new()
  fasttree.set_nt()
  fasttree.set_quote()
  fasttree.set_boot(1000)
  fasttree.set_gtr()
  fasttree.set_gamma()
  fasttree.set_spr(rounds=4)
  fasttree.set_input("alignment.fna")
  fasttree.set_out("tree.nwk")
  println(fasttree.build())

  let fasttree_prot = @src.FastTreeCommandline::new(executable="FastTree")
  fasttree_prot.set_wag()
  fasttree_prot.set_pseudo(weight=2.5)
  fasttree_prot.set_mlnni(rounds=2)
  fasttree_prot.set_input("alignment.faa")
  println(fasttree_prot.build())

  println("\n=== 2. PhyML: maximum likelihood with bootstrap ===")
  let phyml = @src.PhymlCommandline::new()
  phyml.set_input("alignment.phy")
  phyml.set_datatype("nt")
  phyml.set_model("GTR")
  phyml.set_bootstrap(100)
  phyml.set_nclasses(4)
  phyml.set_alpha(1.5)
  phyml.set_search("SPR")
  phyml.set_optimize("tlr")
  phyml.set_rand_start()
  phyml.set_n_rand_starts(5)
  phyml.set_r_seed(12345)
  phyml.set_quiet()
  println(phyml.build())

  println("\n=== 3. RAxML: rapid bootstrapping ===")
  let raxml = @src.RaxmlCommandline::new(executable="raxmlHPC-AVX")
  raxml.set_alignment("alignment.phy")
  raxml.set_name("T1")
  raxml.set_model("PROTGAMMAWAG")
  raxml.set_algorithm("a")
  raxml.set_rapid_bootstrap_seed(12345)
  raxml.set_num_replicates(100)
  raxml.set_bootstrap_branch_lengths()
  raxml.set_threads(4)
  raxml.set_outgroup("OutgroupTaxon")
  raxml.set_working_dir("/tmp/raxml-run")
  println(raxml.build())

  println("\n=== 4. RAxML: partitioned model with constraint tree ===")
  let raxml2 = @src.RaxmlCommandline::new()
  raxml2.set_alignment("alignment.phy")
  raxml2.set_name("partitioned")
  raxml2.set_model("GTRCAT")
  raxml2.set_partition_filename("partitions.model")
  raxml2.set_starting_tree("starting.tre")
  raxml2.set_binary_constraint("constraint.tre")
  raxml2.set_epsilon(0.001)
  raxml2.set_num_categories(25)
  println(raxml2.build())

  println("\nDone.")
}
