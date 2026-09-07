///|
/// Bio.Phylo.Applications - wrappers for phylogenetics command-line tools.
///
/// Port of Biopython's Bio.Phylo.Applications package (FastTree, PhyML and
/// RAxML command-line wrappers), built on the Bio.Application
/// AbstractCommandline base. Each wrapper collects the tool's main
/// parameters in place and serializes them to the canonical command line,
/// e.g.:
///
///   let cmd = FastTreeCommandline::new()
///   cmd.set_nt()
///   cmd.set_input("alignment.fsa")
///   cmd.set_out("tree.newick")
///   cmd.build() == "fasttree -nt -out tree.newick alignment.fsa"
///
/// The input alignment is a positional argument and is appended last.

///|
/// Command-line wrapper for the FastTree tree inference program.
///
/// Reference: Price, M.N., Dehal, P.S., and Arkin, A.P. (2010) FastTree 2 --
/// Approximately Maximum-Likelihood Trees for Large Alignments.
/// PLoS ONE, 5(3):e9490.
pub struct FastTreeCommandline {
  commandline : AbstractCommandline
  priv mut input : String
}

///|
/// Create a FastTree wrapper; `executable` defaults to "fasttree".
pub fn FastTreeCommandline::new(executable? : String) -> FastTreeCommandline {
  FastTreeCommandline::{
    commandline: AbstractCommandline::new(executable.unwrap_or("fasttree")),
    input: "",
  }
}

///|
/// Append a flag with an optional value to the underlying command line.
fn FastTreeCommandline::flag(
  self : FastTreeCommandline,
  name : String,
  value~ : String?,
) -> Unit {
  self.commandline.arguments.push(name)
  if value is Some(v) {
    self.commandline.arguments.push(v)
  }
}

///|
/// Set the input alignment (positional argument, fasta or phylip).
pub fn FastTreeCommandline::set_input(
  self : FastTreeCommandline,
  infile : String,
) -> Unit {
  self.input = infile
}

///|
/// Set the output Newick tree file (-out).
pub fn FastTreeCommandline::set_out(
  self : FastTreeCommandline,
  outfile : String,
) -> Unit {
  self.flag("-out", value=Some(outfile))
}

///|
/// Treat the alignment as nucleotides (-nt); default is protein.
pub fn FastTreeCommandline::set_nt(self : FastTreeCommandline) -> Unit {
  self.flag("-nt", value=None)
}

///|
/// Read N multiple alignments in (-n, phylip interleaved).
pub fn FastTreeCommandline::set_n(self : FastTreeCommandline, count : Int) -> Unit {
  self.flag("-n", value=Some(count.to_string()))
}

///|
/// Quote sequence names in output (-quote).
pub fn FastTreeCommandline::set_quote(self : FastTreeCommandline) -> Unit {
  self.flag("-quote", value=None)
}

///|
/// Use pseudocounts for distance estimation (-pseudo [weight]).
pub fn FastTreeCommandline::set_pseudo(
  self : FastTreeCommandline,
  weight? : Double,
) -> Unit {
  self.flag("-pseudo", value=weight.map(fn(w) { w.to_string() }))
}

///|
/// Set the number of resamples for support values (-boot).
pub fn FastTreeCommandline::set_boot(
  self : FastTreeCommandline,
  resamples : Int,
) -> Unit {
  self.flag("-boot", value=Some(resamples.to_string()))
}

///|
/// Turn off support values (-nosupport).
pub fn FastTreeCommandline::set_nosupport(self : FastTreeCommandline) -> Unit {
  self.flag("-nosupport", value=None)
}

///|
/// Read an initial tree from a file (-intree).
pub fn FastTreeCommandline::set_intree(
  self : FastTreeCommandline,
  file : String,
) -> Unit {
  self.flag("-intree", value=Some(file))
}

///|
/// Suppress printing output except for the final tree (-quiet).
pub fn FastTreeCommandline::set_quiet(self : FastTreeCommandline) -> Unit {
  self.flag("-quiet", value=None)
}

///|
/// Turn off progress reports (-nopr).
pub fn FastTreeCommandline::set_nopr(self : FastTreeCommandline) -> Unit {
  self.flag("-nopr", value=None)
}

///|
/// Rounds of minimum-evolution nearest-neighbor interchanges (-nni [n]).
pub fn FastTreeCommandline::set_nni(
  self : FastTreeCommandline,
  rounds? : Int,
) -> Unit {
  self.flag("-nni", value=rounds.map(fn(n) { n.to_string() }))
}

///|
/// Rounds of subtree-prune-regraft moves (-spr [n]).
pub fn FastTreeCommandline::set_spr(
  self : FastTreeCommandline,
  rounds? : Int,
) -> Unit {
  self.flag("-spr", value=rounds.map(fn(n) { n.to_string() }))
}

///|
/// Turn off maximum-likelihood based refinement (-noml).
pub fn FastTreeCommandline::set_noml(self : FastTreeCommandline) -> Unit {
  self.flag("-noml", value=None)
}

///|
/// Optimize branch lengths but not topology (-mllen).
pub fn FastTreeCommandline::set_mllen(self : FastTreeCommandline) -> Unit {
  self.flag("-mllen", value=None)
}

///|
/// Turn off minimum-evolution support values (-nome).
pub fn FastTreeCommandline::set_nome(self : FastTreeCommandline) -> Unit {
  self.flag("-nome", value=None)
}

///|
/// Use maximum-likelihood NNIs (-mlnni [n]).
pub fn FastTreeCommandline::set_mlnni(
  self : FastTreeCommandline,
  rounds? : Int,
) -> Unit {
  self.flag("-mlnni", value=rounds.map(fn(n) { n.to_string() }))
}

///|
/// Set the level of ML-based topology improvement (-mlacc [n]).
pub fn FastTreeCommandline::set_mlacc(
  self : FastTreeCommandline,
  level? : Int,
) -> Unit {
  self.flag("-mlacc", value=level.map(fn(n) { n.to_string() }))
}

///|
/// Use the slower but accurate NNI implementation (-slownni).
pub fn FastTreeCommandline::set_slownni(self : FastTreeCommandline) -> Unit {
  self.flag("-slownni", value=None)
}

///|
/// Use the WAG protein substitution model (-wag).
pub fn FastTreeCommandline::set_wag(self : FastTreeCommandline) -> Unit {
  self.flag("-wag", value=None)
}

///|
/// Use the GTR nucleotide substitution model (-gtr).
pub fn FastTreeCommandline::set_gtr(self : FastTreeCommandline) -> Unit {
  self.flag("-gtr", value=None)
}

///|
/// Use a discrete-gamma model with n rate categories (-cat n).
pub fn FastTreeCommandline::set_cat(
  self : FastTreeCommandline,
  categories? : Int,
) -> Unit {
  self.flag("-cat", value=categories.map(fn(n) { n.to_string() }))
}

///|
/// Turn off rate categories (-nocat).
pub fn FastTreeCommandline::set_nocat(self : FastTreeCommandline) -> Unit {
  self.flag("-nocat", value=None)
}

///|
/// Use AIC-optimized gamma settings (-gamma).
pub fn FastTreeCommandline::set_gamma(self : FastTreeCommandline) -> Unit {
  self.flag("-gamma", value=None)
}

///|
/// Use the slower NNI/SPR search from version 1 (-slow).
pub fn FastTreeCommandline::set_slow(self : FastTreeCommandline) -> Unit {
  self.flag("-slow", value=None)
}

///|
/// Use the fastest search setting (-fastest).
pub fn FastTreeCommandline::set_fastest(self : FastTreeCommandline) -> Unit {
  self.flag("-fastest", value=None)
}

///|
/// Do a second round of ML NNIs (-2nd).
pub fn FastTreeCommandline::set_second(self : FastTreeCommandline) -> Unit {
  self.flag("-2nd", value=None)
}

///|
/// Do not do a second round of ML NNIs (-no2nd).
pub fn FastTreeCommandline::set_no2nd(self : FastTreeCommandline) -> Unit {
  self.flag("-no2nd", value=None)
}

///|
/// Set the seed for tied ranks and bootstrap (-seed).
pub fn FastTreeCommandline::set_seed(
  self : FastTreeCommandline,
  seed : Int,
) -> Unit {
  self.flag("-seed", value=Some(seed.to_string()))
}

///|
/// Turn off the heuristic used to avoid (-top).
pub fn FastTreeCommandline::set_top(self : FastTreeCommandline) -> Unit {
  self.flag("-top", value=None)
}

///|
/// Turn off "topological" improvements entirely (-notop).
pub fn FastTreeCommandline::set_notop(self : FastTreeCommandline) -> Unit {
  self.flag("-notop", value=None)
}

///|
/// Log output to a file (-log).
pub fn FastTreeCommandline::set_log(
  self : FastTreeCommandline,
  file : String,
) -> Unit {
  self.flag("-log", value=Some(file))
}

///|
/// Do not write a distance matrix (-nomatrix).
pub fn FastTreeCommandline::set_nomatrix(self : FastTreeCommandline) -> Unit {
  self.flag("-nomatrix", value=None)
}

///|
/// Use neighbor-joining (default) instead of BioNJ (-nj).
pub fn FastTreeCommandline::set_nj(self : FastTreeCommandline) -> Unit {
  self.flag("-nj", value=None)
}

///|
/// Use BioNJ instead of neighbor-joining (-bionj).
pub fn FastTreeCommandline::set_bionj(self : FastTreeCommandline) -> Unit {
  self.flag("-bionj", value=None)
}

///|
/// Specify GTR rates (-gtrrates, e.g. "1 1 1 1 1 1").
pub fn FastTreeCommandline::set_gtrrates(
  self : FastTreeCommandline,
  rates : String,
) -> Unit {
  self.flag("-gtrrates", value=Some(rates))
}

///|
/// Specify GTR nucleotide frequencies (-gtrfreq, e.g. "0.25 0.25 0.25 0.25").
pub fn FastTreeCommandline::set_gtrfreq(
  self : FastTreeCommandline,
  freqs : String,
) -> Unit {
  self.flag("-gtrfreq", value=Some(freqs))
}

///|
/// Set the maximum SPR move length (-sprlength).
pub fn FastTreeCommandline::set_sprlength(
  self : FastTreeCommandline,
  length : Int,
) -> Unit {
  self.flag("-sprlength", value=Some(length.to_string()))
}

///|
/// Serialize the command line; the input alignment is appended last as a
/// positional argument.
pub fn FastTreeCommandline::build(self : FastTreeCommandline) -> String {
  let base = self.commandline.build_command()
  if self.input.length() > 0 {
    base + " " + self.input
  } else {
    base
  }
}

///|
/// Command-line wrapper for the PhyML maximum-likelihood tree program.
///
/// Reference: Guindon, S. & Gascuel, O. (2003) "A simple, fast, and accurate
/// algorithm to estimate large phylogenies by maximum likelihood".
/// Systematic Biology 52(5):696-704.
pub struct PhymlCommandline {
  commandline : AbstractCommandline
}

///|
/// Create a PhyML wrapper; `executable` defaults to "phyml".
pub fn PhymlCommandline::new(executable? : String) -> PhymlCommandline {
  PhymlCommandline::{
    commandline: AbstractCommandline::new(executable.unwrap_or("phyml")),
  }
}

///|
/// Append a flag with an optional value to the underlying command line.
fn PhymlCommandline::flag(
  self : PhymlCommandline,
  name : String,
  value~ : String?,
) -> Unit {
  self.commandline.arguments.push(name)
  if value is Some(v) {
    self.commandline.arguments.push(v)
  }
}

///|
/// Set the PHYLIP format input file (-i, required).
pub fn PhymlCommandline::set_input(
  self : PhymlCommandline,
  infile : String,
) -> Unit {
  self.flag("-i", value=Some(infile))
}

///|
/// Set the datatype: "nt" (default) or "aa" (-d).
pub fn PhymlCommandline::set_datatype(
  self : PhymlCommandline,
  datatype : String,
) -> Unit {
  self.flag("-d", value=Some(datatype))
}

///|
/// Change interleaved format to sequential format (-q).
pub fn PhymlCommandline::set_sequential(self : PhymlCommandline) -> Unit {
  self.flag("-q", value=None)
}

///|
/// Set the number of data sets to analyse (-n).
pub fn PhymlCommandline::set_multiple(self : PhymlCommandline, count : Int) -> Unit {
  self.flag("-n", value=Some(count.to_string()))
}

///|
/// Set the parsimony starting tree value (-p).
pub fn PhymlCommandline::set_pars(self : PhymlCommandline, pars : Double) -> Unit {
  self.flag("-p", value=Some(pars.to_string()))
}

///|
/// Set the number of bootstrap replicates (-b).
pub fn PhymlCommandline::set_bootstrap(
  self : PhymlCommandline,
  replicates : Int,
) -> Unit {
  self.flag("-b", value=Some(replicates.to_string()))
}

///|
/// Set the substitution model (-m).
pub fn PhymlCommandline::set_model(
  self : PhymlCommandline,
  model : String,
) -> Unit {
  self.flag("-m", value=Some(model))
}

///|
/// Set the base frequency method, e.g. "e", "m" or custom values (-f).
pub fn PhymlCommandline::set_frequencies(
  self : PhymlCommandline,
  frequencies : String,
) -> Unit {
  self.flag("-f", value=Some(frequencies))
}

///|
/// Set the transition/transversion ratio (-t).
pub fn PhymlCommandline::set_ts_tv_ratio(
  self : PhymlCommandline,
  ratio : Double,
) -> Unit {
  self.flag("-t", value=Some(ratio.to_string()))
}

///|
/// Set the proportion of invariable sites (-v).
pub fn PhymlCommandline::set_prop_invar(
  self : PhymlCommandline,
  prop : Double,
) -> Unit {
  self.flag("-v", value=Some(prop.to_string()))
}

///|
/// Set the number of relative substitution rate categories (-c).
pub fn PhymlCommandline::set_nclasses(
  self : PhymlCommandline,
  classes : Int,
) -> Unit {
  self.flag("-c", value=Some(classes.to_string()))
}

///|
/// Set the gamma distribution parameter (-a).
pub fn PhymlCommandline::set_alpha(
  self : PhymlCommandline,
  alpha : Double,
) -> Unit {
  self.flag("-a", value=Some(alpha.to_string()))
}

///|
/// Set the tree-search operation, e.g. "NNI" or "SPR" (-s).
pub fn PhymlCommandline::set_search(
  self : PhymlCommandline,
  search : String,
) -> Unit {
  self.flag("-s", value=Some(search))
}

///|
/// Set an input tree file (-u).
pub fn PhymlCommandline::set_input_tree(
  self : PhymlCommandline,
  treefile : String,
) -> Unit {
  self.flag("-u", value=Some(treefile))
}

///|
/// Set the optimisation targets, e.g. "tlr", "tl", "lr" (-o).
pub fn PhymlCommandline::set_optimize(
  self : PhymlCommandline,
  optimize : String,
) -> Unit {
  self.flag("-o", value=Some(optimize))
}

///|
/// Enable random starting trees (--rand_start).
pub fn PhymlCommandline::set_rand_start(self : PhymlCommandline) -> Unit {
  self.flag("--rand_start", value=None)
}

///|
/// Set the number of random starting trees (--n_rand_starts).
pub fn PhymlCommandline::set_n_rand_starts(
  self : PhymlCommandline,
  count : Int,
) -> Unit {
  self.flag("--n_rand_starts", value=Some(count.to_string()))
}

///|
/// Set the random seed (--r_seed).
pub fn PhymlCommandline::set_r_seed(self : PhymlCommandline, seed : Int) -> Unit {
  self.flag("--r_seed", value=Some(seed.to_string()))
}

///|
/// Print site log-likelihoods (--print_site_lnl).
pub fn PhymlCommandline::set_print_site_lnl(self : PhymlCommandline) -> Unit {
  self.flag("--print_site_lnl", value=None)
}

///|
/// Print the trace of the search (--print_trace).
pub fn PhymlCommandline::set_print_trace(self : PhymlCommandline) -> Unit {
  self.flag("--print_trace", value=None)
}

///|
/// Set a run identifier for output files (--run_id).
pub fn PhymlCommandline::set_run_id(
  self : PhymlCommandline,
  run_id : String,
) -> Unit {
  self.flag("--run_id", value=Some(run_id))
}

///|
/// Run quietly (--quiet).
pub fn PhymlCommandline::set_quiet(self : PhymlCommandline) -> Unit {
  self.flag("--quiet", value=None)
}

///|
/// Serialize the command line.
pub fn PhymlCommandline::build(self : PhymlCommandline) -> String {
  self.commandline.build_command()
}

///|
/// Command-line wrapper for the RAxML tree inference program (raxmlHPC).
///
/// Reference: Stamatakis, A. (2014) "RAxML version 8: a tool for
/// phylogenetic analysis and post-analysis of large phylogenies".
/// Bioinformatics 30(9):1312-1313.
pub struct RaxmlCommandline {
  commandline : AbstractCommandline
}

///|
/// Create a RAxML wrapper; `executable` defaults to "raxmlHPC".
pub fn RaxmlCommandline::new(executable? : String) -> RaxmlCommandline {
  RaxmlCommandline::{
    commandline: AbstractCommandline::new(executable.unwrap_or("raxmlHPC")),
  }
}

///|
/// Append a flag with an optional value to the underlying command line.
fn RaxmlCommandline::flag(
  self : RaxmlCommandline,
  name : String,
  value~ : String?,
) -> Unit {
  self.commandline.arguments.push(name)
  if value is Some(v) {
    self.commandline.arguments.push(v)
  }
}

///|
/// Set the sequence alignment file (-s).
pub fn RaxmlCommandline::set_alignment(
  self : RaxmlCommandline,
  alignment : String,
) -> Unit {
  self.flag("-s", value=Some(alignment))
}

///|
/// Set the run name (-n).
pub fn RaxmlCommandline::set_name(self : RaxmlCommandline, name : String) -> Unit {
  self.flag("-n", value=Some(name))
}

///|
/// Set the substitution model (-m).
pub fn RaxmlCommandline::set_model(
  self : RaxmlCommandline,
  model : String,
) -> Unit {
  self.flag("-m", value=Some(model))
}

///|
/// Set the algorithm type, e.g. "a", "b", "d" (-f).
pub fn RaxmlCommandline::set_algorithm(
  self : RaxmlCommandline,
  algorithm : String,
) -> Unit {
  self.flag("-f", value=Some(algorithm))
}

///|
/// Set the column weight filename (-a).
pub fn RaxmlCommandline::set_weight_filename(
  self : RaxmlCommandline,
  filename : String,
) -> Unit {
  self.flag("-a", value=Some(filename))
}

///|
/// Set the random seed for bootstrapping (-b).
pub fn RaxmlCommandline::set_bootstrap_seed(
  self : RaxmlCommandline,
  seed : Int,
) -> Unit {
  self.flag("-b", value=Some(seed.to_string()))
}

///|
/// Set the number of categories for the per-site rate model (-c).
pub fn RaxmlCommandline::set_num_categories(
  self : RaxmlCommandline,
  count : Int,
) -> Unit {
  self.flag("-c", value=Some(count.to_string()))
}

///|
/// Start from a random tree instead of a parsimony tree (-d).
pub fn RaxmlCommandline::set_random_starting_tree(self : RaxmlCommandline) -> Unit {
  self.flag("-d", value=None)
}

///|
/// Set the likelihood epsilon (-e).
pub fn RaxmlCommandline::set_epsilon(
  self : RaxmlCommandline,
  epsilon : Double,
) -> Unit {
  self.flag("-e", value=Some(epsilon.to_string()))
}

///|
/// Set the file with sites to exclude (-E).
pub fn RaxmlCommandline::set_exclude_filename(
  self : RaxmlCommandline,
  filename : String,
) -> Unit {
  self.flag("-E", value=Some(filename))
}

///|
/// Set the multifurcating constraint file (-g).
pub fn RaxmlCommandline::set_grouping_constraint(
  self : RaxmlCommandline,
  filename : String,
) -> Unit {
  self.flag("-g", value=Some(filename))
}

///|
/// Set the number of rearrangement steps (-i).
pub fn RaxmlCommandline::set_rearrangements(
  self : RaxmlCommandline,
  count : Int,
) -> Unit {
  self.flag("-i", value=Some(count.to_string()))
}

///|
/// Enable checkpointing (-j).
pub fn RaxmlCommandline::set_checkpoints(self : RaxmlCommandline) -> Unit {
  self.flag("-j", value=None)
}

///|
/// Print branch lengths with bootstrap support values (-k).
pub fn RaxmlCommandline::set_bootstrap_branch_lengths(
  self : RaxmlCommandline,
) -> Unit {
  self.flag("-k", value=None)
}

///|
/// Set the clustering threshold (-l).
pub fn RaxmlCommandline::set_cluster_threshold(
  self : RaxmlCommandline,
  threshold : Double,
) -> Unit {
  self.flag("-l", value=Some(threshold.to_string()))
}

///|
/// Set the fast clustering threshold (-L).
pub fn RaxmlCommandline::set_cluster_threshold_fast(
  self : RaxmlCommandline,
  threshold : Double,
) -> Unit {
  self.flag("-L", value=Some(threshold.to_string()))
}

///|
/// Enable per-partition branch length optimization (-M).
pub fn RaxmlCommandline::set_partition_branch_lengths(
  self : RaxmlCommandline,
) -> Unit {
  self.flag("-M", value=None)
}

///|
/// Set the outgroup taxon name (-o).
pub fn RaxmlCommandline::set_outgroup(
  self : RaxmlCommandline,
  outgroup : String,
) -> Unit {
  self.flag("-o", value=Some(outgroup))
}

///|
/// Set the partitioned-model filename (-q).
pub fn RaxmlCommandline::set_partition_filename(
  self : RaxmlCommandline,
  filename : String,
) -> Unit {
  self.flag("-q", value=Some(filename))
}

///|
/// Set the parsimony random seed (-p).
pub fn RaxmlCommandline::set_parsimony_seed(
  self : RaxmlCommandline,
  seed : Int,
) -> Unit {
  self.flag("-p", value=Some(seed.to_string()))
}

///|
/// Set the protein substitution model (-P).
pub fn RaxmlCommandline::set_protein_model(
  self : RaxmlCommandline,
  model : String,
) -> Unit {
  self.flag("-P", value=Some(model))
}

///|
/// Set the binary constraint tree file (-r).
pub fn RaxmlCommandline::set_binary_constraint(
  self : RaxmlCommandline,
  filename : String,
) -> Unit {
  self.flag("-r", value=Some(filename))
}

///|
/// Set the starting tree file (-t).
pub fn RaxmlCommandline::set_starting_tree(
  self : RaxmlCommandline,
  treefile : String,
) -> Unit {
  self.flag("-t", value=Some(treefile))
}

///|
/// Set the number of threads (-T).
pub fn RaxmlCommandline::set_threads(self : RaxmlCommandline, threads : Int) -> Unit {
  self.flag("-T", value=Some(threads.to_string()))
}

///|
/// Set the number of bootstrap searches (-u).
pub fn RaxmlCommandline::set_num_bootstrap_searches(
  self : RaxmlCommandline,
  count : Int,
) -> Unit {
  self.flag("-u", value=Some(count.to_string()))
}

///|
/// Display version information (-v).
pub fn RaxmlCommandline::set_version(self : RaxmlCommandline) -> Unit {
  self.flag("-v", value=None)
}

///|
/// Set the working directory for output files (-w).
pub fn RaxmlCommandline::set_working_dir(
  self : RaxmlCommandline,
  directory : String,
) -> Unit {
  self.flag("-w", value=Some(directory))
}

///|
/// Set the random seed for rapid bootstrapping (-x).
pub fn RaxmlCommandline::set_rapid_bootstrap_seed(
  self : RaxmlCommandline,
  seed : Int,
) -> Unit {
  self.flag("-x", value=Some(seed.to_string()))
}

///|
/// Compute a parsimony starting tree only (-y).
pub fn RaxmlCommandline::set_parsimony(self : RaxmlCommandline) -> Unit {
  self.flag("-y", value=None)
}

///|
/// Set the bipartition filename to draw support values (-z).
pub fn RaxmlCommandline::set_bipartition_filename(
  self : RaxmlCommandline,
  filename : String,
) -> Unit {
  self.flag("-z", value=Some(filename))
}

///|
/// Set the number of bootstrap replicates (-N, also written -#).
pub fn RaxmlCommandline::set_num_replicates(
  self : RaxmlCommandline,
  count : Int,
) -> Unit {
  self.flag("-N", value=Some(count.to_string()))
}

///|
/// Serialize the command line.
pub fn RaxmlCommandline::build(self : RaxmlCommandline) -> String {
  self.commandline.build_command()
}
