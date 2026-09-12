///|
/// Tests for the Bio.Phylo.Applications command-line wrappers (FastTree,
/// PhyML and RAxML). The wrappers are pure parameter collectors built on
/// Bio.Application's AbstractCommandline, so the tests verify the exact
/// serialized command lines without executing any external tool.

///|
test "fasttree_defaults" {
  let cmd = @src.FastTreeCommandline::new()
  assert_eq(cmd.build(), "fasttree")
  let cmd2 = @src.FastTreeCommandline::new(executable="FastTree")
  assert_eq(cmd2.build(), "FastTree")
}

///|
test "fasttree_nucleotide_basic" {
  // positional input alignment is appended last
  let cmd = @src.FastTreeCommandline::new()
  cmd.set_nt()
  cmd.set_input("alignment.fsa")
  cmd.set_out("tree.newick")
  assert_eq(cmd.build(), "fasttree -nt -out tree.newick alignment.fsa")
}

///|
test "fasttree_support_and_search_options" {
  let cmd = @src.FastTreeCommandline::new()
  cmd.set_nt()
  cmd.set_quote()
  cmd.set_boot(1000)
  cmd.set_seed(12345)
  cmd.set_spr(rounds=4)
  cmd.set_mlacc(level=2)
  cmd.set_input("align.fsa")
  assert_eq(
    cmd.build(),
    "fasttree -nt -quote -boot 1000 -seed 12345 -spr 4 -mlacc 2 align.fsa",
  )
}

///|
test "fasttree_flag_without_value" {
  // optional-value flags can be set without a value
  let cmd = @src.FastTreeCommandline::new()
  cmd.set_pseudo()
  cmd.set_nni()
  cmd.set_cat()
  cmd.set_input("aln.fsa")
  assert_eq(cmd.build(), "fasttree -pseudo -nni -cat aln.fsa")
  // ...or with one
  let cmd2 = @src.FastTreeCommandline::new()
  cmd2.set_pseudo(weight=2.5)
  cmd2.set_nni(rounds=1)
  cmd2.set_cat(categories=20)
  cmd2.set_input("aln.fsa")
  assert_eq(cmd2.build(), "fasttree -pseudo 2.5 -nni 1 -cat 20 aln.fsa")
}

///|
test "fasttree_model_and_ml_options" {
  let cmd = @src.FastTreeCommandline::new()
  cmd.set_gtr()
  cmd.set_gamma()
  cmd.set_gtrrates("1 1 1 1 1 1")
  cmd.set_gtrfreq("0.25 0.25 0.25 0.25")
  cmd.set_intree("starting.tre")
  cmd.set_log("fasttree.log")
  cmd.set_quiet()
  cmd.set_input("aln.fna")
  assert_eq(
    cmd.build(),
    "fasttree -gtr -gamma -gtrrates 1 1 1 1 1 1 -gtrfreq 0.25 0.25 0.25 0.25 -intree starting.tre -log fasttree.log -quiet aln.fna",
  )
}

///|
test "fasttree_off_switches" {
  let cmd = @src.FastTreeCommandline::new()
  cmd.set_noml()
  cmd.set_mllen()
  cmd.set_nome()
  cmd.set_nosupport()
  cmd.set_nomatrix()
  cmd.set_nopr()
  cmd.set_input("aln.fsa")
  assert_eq(
    cmd.build(),
    "fasttree -noml -mllen -nome -nosupport -nomatrix -nopr aln.fsa",
  )
}

///|
test "fasttree_misc_flags" {
  let cmd = @src.FastTreeCommandline::new()
  cmd.set_bionj()
  cmd.set_slownni()
  cmd.set_wag()
  cmd.set_slow()
  cmd.set_fastest()
  cmd.set_second()
  cmd.set_top()
  cmd.set_sprlength(100)
  cmd.set_n(2)
  cmd.set_input("aln.fsa")
  assert_eq(
    cmd.build(),
    "fasttree -bionj -slownni -wag -slow -fastest -2nd -top -sprlength 100 -n 2 aln.fsa",
  )
  // the negative switches are separate flags
  let cmd2 = @src.FastTreeCommandline::new()
  cmd2.set_no2nd()
  cmd2.set_notop()
  cmd2.set_nocat()
  cmd2.set_nj()
  cmd2.set_input("aln.fsa")
  assert_eq(cmd2.build(), "fasttree -no2nd -notop -nocat -nj aln.fsa")
}

///|
test "phyml_basic_protein" {
  let cmd = @src.PhymlCommandline::new()
  cmd.set_input("align.phy")
  cmd.set_datatype("aa")
  cmd.set_model("WAG")
  cmd.set_quiet()
  assert_eq(cmd.build(), "phyml -i align.phy -d aa -m WAG --quiet")
}

///|
test "phyml_bootstrap_and_gamma" {
  let cmd = @src.PhymlCommandline::new(executable="phyml")
  cmd.set_input("seqs.phy")
  cmd.set_bootstrap(100)
  cmd.set_nclasses(4)
  cmd.set_alpha(1.5)
  cmd.set_prop_invar(0.25)
  cmd.set_ts_tv_ratio(4.0)
  cmd.set_search("SPR")
  cmd.set_optimize("tlr")
  assert_eq(
    cmd.build(),
    "phyml -i seqs.phy -b 100 -c 4 -a 1.5 -v 0.25 -t 4 -s SPR -o tlr",
  )
}

///|
test "phyml_random_starts_and_output" {
  let cmd = @src.PhymlCommandline::new()
  cmd.set_input("aln.phy")
  cmd.set_rand_start()
  cmd.set_n_rand_starts(5)
  cmd.set_r_seed(7)
  cmd.set_print_site_lnl()
  cmd.set_print_trace()
  cmd.set_run_id("run1")
  cmd.set_frequencies("e")
  cmd.set_multiple(1)
  cmd.set_sequential()
  cmd.set_pars(0.5)
  cmd.set_input_tree("tree.phy")
  assert_eq(
    cmd.build(),
    "phyml -i aln.phy --rand_start --n_rand_starts 5 --r_seed 7 --print_site_lnl --print_trace --run_id run1 -f e -n 1 -q -p 0.5 -u tree.phy",
  )
}

///|
test "raxml_basic_mla" {
  let cmd = @src.RaxmlCommandline::new()
  cmd.set_alignment("align.phy")
  cmd.set_name("run1")
  cmd.set_model("GTRGAMMA")
  cmd.set_algorithm("a")
  cmd.set_threads(4)
  assert_eq(cmd.build(), "raxmlHPC -s align.phy -n run1 -m GTRGAMMA -f a -T 4")
}

///|
test "raxml_rapid_bootstrap" {
  let cmd = @src.RaxmlCommandline::new(executable="raxmlHPC-AVX")
  cmd.set_alignment("aln.phy")
  cmd.set_name("T1")
  cmd.set_model("PROTGAMMAWAG")
  cmd.set_algorithm("a")
  cmd.set_rapid_bootstrap_seed(12345)
  cmd.set_num_replicates(100)
  cmd.set_bootstrap_branch_lengths()
  cmd.set_outgroup("OutSpec")
  cmd.set_working_dir("/tmp/tree")
  assert_eq(
    cmd.build(),
    "raxmlHPC-AVX -s aln.phy -n T1 -m PROTGAMMAWAG -f a -x 12345 -N 100 -k -o OutSpec -w /tmp/tree",
  )
}

///|
test "raxml_tree_and_constraint_files" {
  let cmd = @src.RaxmlCommandline::new()
  cmd.set_alignment("aln.phy")
  cmd.set_name("r")
  cmd.set_starting_tree("start.tre")
  cmd.set_binary_constraint("const.tre")
  cmd.set_grouping_constraint("multi.tre")
  cmd.set_bipartition_filename("bipart.tre")
  cmd.set_partition_filename("parts.model")
  cmd.set_exclude_filename("exclude.dat")
  cmd.set_weight_filename("weights.dat")
  cmd.set_protein_model("JTT")
  assert_eq(
    cmd.build(),
    "raxmlHPC -s aln.phy -n r -t start.tre -r const.tre -g multi.tre -z bipart.tre -q parts.model -E exclude.dat -a weights.dat -P JTT",
  )
}

///|
test "raxml_numeric_and_flag_options" {
  let cmd = @src.RaxmlCommandline::new()
  cmd.set_alignment("aln.phy")
  cmd.set_name("r")
  cmd.set_bootstrap_seed(9)
  cmd.set_parsimony_seed(8)
  cmd.set_num_bootstrap_searches(10)
  cmd.set_num_categories(25)
  cmd.set_epsilon(0.001)
  cmd.set_cluster_threshold(0.9)
  cmd.set_cluster_threshold_fast(0.95)
  cmd.set_rearrangements(25)
  cmd.set_random_starting_tree()
  cmd.set_checkpoints()
  cmd.set_partition_branch_lengths()
  cmd.set_parsimony()
  cmd.set_version()
  assert_eq(
    cmd.build(),
    "raxmlHPC -s aln.phy -n r -b 9 -p 8 -u 10 -c 25 -e 0.001 -l 0.9 -L 0.95 -i 25 -d -j -M -y -v",
  )
}
