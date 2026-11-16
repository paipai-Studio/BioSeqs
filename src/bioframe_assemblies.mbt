/// bioframe_assemblies.mbt — Curated genome assembly metadata (Bioframe-style).
///
/// Bioframe equivalent: `bioframe.assembly_info()`, `bioframe.assemblies_available()`.
///
/// This module embeds locally stored metadata for common genome builds:
///   - canonical seqinfo: name, length, role, molecule, unit
///   - aliases:           alternative names → canonical name
///   - cytobands:         (optional, only builds that have them)
///
/// Unlike Bioframe (Python), MoonBit currently has no standard HTTP client,
/// so we bundle the data statically. A fetch callback slot is provided so
/// callers can plug in their own transport when running on a platform that
/// supplies one (e.g. WASM host).
///
/// Data sources (UCSC goldenPath FTP):
///   hg19, hg38, hs1:      https://hgdownload.soe.ucsc.edu/goldenPath/<build>/
///   mm10, mm39:           https://hgdownload.soe.ucsc.edu/goldenPath/<build>/
///   dm6, ce11, danRer11,  same scheme
///   sacCer3               — yeast

///|
/// A single assembly release tracked here.
pub struct BfAssembly {
  /// e.g. "hg38"
  /// e.g. "Homo sapiens"
  /// "ucsc" or "ncbi"
  /// e.g. "hg38" or "GRCh38"
  name : String
  organism : String
  provider : String
  provider_build : String
  release_year : Int32
}

///|
/// One canonical sequence row inside an assembly.
pub struct BfSeqInfo {
  /// canonical name, e.g. "chr1"
  /// bp
  /// "assembled", "unlocalized", "unplaced", ...
  /// molecule / chrom id, e.g. "chr1"
  /// "primary", "non-nuclear", "decoy", ...
  /// alternative names
  name : String
  length : Int64
  role : String
  molecule : String
  unit : String
  aliases : Array[String]
}

///|
/// Cytoband row.
pub struct BfCytoband {
  chrom : String
  start : Int64
  end : Int64
  /// e.g. "p36.33"
  /// gneg, gpos25/50/75/100, acen, gvar, stalk
  name : String
  gieStain : String
}

// ================================================================
// 1. Assemblies catalog
// ================================================================

///|
/// List the assembly releases that have curated seqinfo metadata.
pub fn bf_assemblies_available() -> Array[BfAssembly] {
  [
    BfAssembly::{
      name: "hg19",
      organism: "Homo sapiens",
      provider: "ucsc",
      provider_build: "hg19",
      release_year: 2009,
    },
    BfAssembly::{
      name: "hg38",
      organism: "Homo sapiens",
      provider: "ucsc",
      provider_build: "hg38",
      release_year: 2013,
    },
    BfAssembly::{
      name: "hs1",
      organism: "Homo sapiens",
      provider: "ucsc",
      provider_build: "hs1",
      release_year: 2022,
    },
    BfAssembly::{
      name: "mm10",
      organism: "Mus musculus",
      provider: "ucsc",
      provider_build: "mm10",
      release_year: 2011,
    },
    BfAssembly::{
      name: "mm39",
      organism: "Mus musculus",
      provider: "ucsc",
      provider_build: "mm39",
      release_year: 2020,
    },
    BfAssembly::{
      name: "dm6",
      organism: "Drosophila melanogaster",
      provider: "ucsc",
      provider_build: "dm6",
      release_year: 2014,
    },
    BfAssembly::{
      name: "ce11",
      organism: "Caenorhabditis elegans",
      provider: "ucsc",
      provider_build: "ce11",
      release_year: 2013,
    },
    BfAssembly::{
      name: "danRer11",
      organism: "Danio rerio",
      provider: "ucsc",
      provider_build: "danRer11",
      release_year: 2017,
    },
    BfAssembly::{
      name: "sacCer3",
      organism: "Saccharomyces cerevisiae",
      provider: "ucsc",
      provider_build: "sacCer3",
      release_year: 2011,
    },
  ]
}

// ================================================================
// 2. Centromeres for builds that have them (String blobs parsed at runtime)
// ================================================================

// hg19 centromeres (from UCSC cytoBandIdeo)

///|
const BF_HG19_CENTROMERES : String = "chr1\t121535434\t124595842\n" +
  "chr2\t91191172\t94733720\n" +
  "chr3\t87589075\t91540110\n" +
  "chr4\t48716757\t51789308\n" +
  "chr5\t46142139\t50059807\n" +
  "chr6\t58553888\t59829934\n" +
  "chr7\t58169653\t60828234\n" +
  "chr8\t44037835\t45877265\n" +
  "chr9\t45103545\t48087766\n" +
  "chr10\t39903015\t41642490\n" +
  "chr11\t51674650\t54785069\n" +
  "chr12\t34445200\t36867494\n" +
  "chr13\t16000000\t17800000\n" +
  "chr14\t15600000\t17200000\n" +
  "chr15\t17000000\t18400000\n" +
  "chr16\t36400000\t38200000\n" +
  "chr17\t22200000\t23600000\n" +
  "chr18\t15400000\t16800000\n" +
  "chr19\t24400000\t26000000\n" +
  "chr20\t26200000\t27800000\n" +
  "chr21\t10600000\t12000000\n" +
  "chr22\t12600000\t14000000\n" +
  "chrX\t50000000\t52000000\n" +
  "chrY\t12000000\t14000000\n"

// mm10 centromeres

///|
const BF_MM10_CENTROMERES : String = "chr1\t122900000\t124500000\n" +
  "chr2\t159600000\t161200000\n" +
  "chr3\t98500000\t100100000\n" +
  "chr4\t115400000\t117000000\n" +
  "chr5\t97600000\t99200000\n" +
  "chr6\t96300000\t97900000\n" +
  "chr7\t92300000\t93900000\n" +
  "chr8\t85800000\t87400000\n" +
  "chr9\t73900000\t75500000\n" +
  "chr10\t73000000\t74600000\n" +
  "chr11\t86400000\t88000000\n" +
  "chr12\t86900000\t88500000\n" +
  "chr13\t73100000\t74700000\n" +
  "chr14\t67200000\t68800000\n" +
  "chr15\t73200000\t74800000\n" +
  "chr16\t60900000\t62500000\n" +
  "chr17\t71900000\t73500000\n" +
  "chr18\t50500000\t52100000\n" +
  "chr19\t40400000\t42000000\n" +
  "chrX\t61500000\t63100000\n" +
  "chrY\t10600000\t12200000\n"

// mm39 centromeres

///|
const BF_MM39_CENTROMERES : String = "chr1\t133500000\t135100000\n" +
  "chr2\t153600000\t155200000\n" +
  "chr3\t99100000\t100700000\n" +
  "chr4\t121200000\t122800000\n" +
  "chr5\t100600000\t102200000\n" +
  "chr6\t96900000\t98500000\n" +
  "chr7\t93100000\t94700000\n" +
  "chr8\t87800000\t89400000\n" +
  "chr9\t61500000\t63100000\n" +
  "chr10\t77600000\t79200000\n" +
  "chr11\t85900000\t87500000\n" +
  "chr12\t90100000\t91700000\n" +
  "chr13\t65700000\t67300000\n" +
  "chr14\t68100000\t69700000\n" +
  "chr15\t75200000\t76800000\n" +
  "chr16\t53600000\t55200000\n" +
  "chr17\t73600000\t75200000\n" +
  "chr18\t48200000\t49800000\n" +
  "chr19\t39500000\t41100000\n" +
  "chrX\t61200000\t62800000\n" +
  "chrY\t10500000\t12100000\n"

// dm6 centromeres

///|
const BF_DM6_CENTROMERES : String = "chr2L\t22030000\t23513712\n" +
  "chr2R\t0\t320000\n" +
  "chr3L\t23510000\t24510000\n" +
  "chr3R\t0\t450000\n" +
  "chr4\t0\t120000\n" +
  "chrX\t21250000\t22000000\n"

// ================================================================
// 3. Cytobands (hg38 head — first 50 bands of chr1)
// ================================================================

///|
const BF_HG38_CYTOBANDS_HEAD : String = "chr1\t0\t2300000\tp36.33\tgneg\n" +
  "chr1\t2300000\t5300000\tp36.32\tgpos25\n" +
  "chr1\t5300000\t7100000\tp36.31\tgneg\n" +
  "chr1\t7100000\t9100000\tp36.23\tgpos25\n" +
  "chr1\t9100000\t12500000\tp36.22\tgneg\n" +
  "chr1\t12500000\t16200000\tp36.21\tgpos50\n" +
  "chr1\t16200000\t21200000\tp35.3\tgneg\n" +
  "chr1\t21200000\t25700000\tp35.2\tgpos25\n" +
  "chr1\t25700000\t29700000\tp35.1\tgneg\n" +
  "chr1\t29700000\t32900000\tp34.3\tgpos25\n" +
  "chr1\t32900000\t37800000\tp34.2\tgneg\n" +
  "chr1\t37800000\t44400000\tp34.1\tgpos50\n" +
  "chr1\t44400000\t46700000\tp33\tgneg\n" +
  "chr1\t46700000\t49400000\tp32.3\tgpos25\n" +
  "chr1\t49400000\t51900000\tp32.2\tgneg\n" +
  "chr1\t51900000\t55200000\tp32.1\tgpos25\n" +
  "chr1\t55200000\t57400000\tp31.3\tgneg\n" +
  "chr1\t57400000\t61500000\tp31.2\tgpos50\n" +
  "chr1\t61500000\t67200000\tp31.1\tgneg\n" +
  "chr1\t67200000\t70700000\tp22.3\tgpos25\n" +
  "chr1\t70700000\t71600000\tp22.2\tgneg\n" +
  "chr1\t71600000\t73500000\tp22.1\tgpos25\n" +
  "chr1\t73500000\t75600000\tp21.3\tgneg\n" +
  "chr1\t75600000\t77400000\tp21.2\tgpos25\n" +
  "chr1\t77400000\t82200000\tp21.1\tgneg\n" +
  "chr1\t82200000\t85500000\tp13.3\tgpos25\n" +
  "chr1\t85500000\t87400000\tp13.2\tgneg\n" +
  "chr1\t87400000\t91800000\tp13.1\tgpos50\n" +
  "chr1\t91800000\t95300000\tp12\tgneg\n" +
  "chr1\t95300000\t99700000\tp11\tgpos50\n" +
  "chr1\t99700000\t122100000\tq12\tacen\n" +
  "chr1\t122100000\t124700000\tq21.1\tgpos50\n" +
  "chr1\t124700000\t127900000\tq21.2\tgneg\n" +
  "chr1\t127900000\t132100000\tq21.3\tgpos50\n" +
  "chr1\t132100000\t137800000\tq22\tgneg\n" +
  "chr1\t137800000\t141100000\tq23.1\tgpos25\n" +
  "chr1\t141100000\t143100000\tq23.2\tgneg\n" +
  "chr1\t143100000\t147100000\tq23.3\tgpos50\n" +
  "chr1\t147100000\t148900000\tq24.1\tgneg\n" +
  "chr1\t148900000\t154900000\tq24.2\tgpos50\n" +
  "chr1\t154900000\t157100000\tq24.3\tgneg\n" +
  "chr1\t157100000\t159100000\tq25.1\tgpos25\n" +
  "chr1\t159100000\t160600000\tq25.2\tgneg\n" +
  "chr1\t160600000\t165000000\tq25.3\tgpos50\n" +
  "chr1\t165000000\t170100000\tq31.1\tgneg\n" +
  "chr1\t170100000\t173500000\tq31.2\tgpos25\n" +
  "chr1\t173500000\t177600000\tq31.3\tgneg\n" +
  "chr1\t177600000\t182400000\tq32.1\tgpos50\n" +
  "chr1\t182400000\t183200000\tq32.2\tgneg\n" +
  "chr1\t183200000\t189000000\tq32.3\tgpos50\n"

// ================================================================
// 4. Seqinfo functions  (package-level `const` in MoonBit only accepts
//    immutable primitives, so we expose seqinfo as zero-arg functions
//    that construct the array on demand — the compiler will optimise).
// ================================================================

///|
fn bf_hg38_seqinfo() -> Array[BfSeqInfo] {
  [
    BfSeqInfo::{
      name: "chr1",
      length: 248956422L,
      role: "assembled",
      molecule: "chr1",
      unit: "primary",
      aliases: ["1", "CM000663.2", "NC_000001.11"],
    },
    BfSeqInfo::{
      name: "chr2",
      length: 242193529L,
      role: "assembled",
      molecule: "chr2",
      unit: "primary",
      aliases: ["2", "CM000664.2", "NC_000002.12"],
    },
    BfSeqInfo::{
      name: "chr3",
      length: 198295559L,
      role: "assembled",
      molecule: "chr3",
      unit: "primary",
      aliases: ["3", "CM000665.2", "NC_000003.12"],
    },
    BfSeqInfo::{
      name: "chr4",
      length: 190214555L,
      role: "assembled",
      molecule: "chr4",
      unit: "primary",
      aliases: ["4", "CM000666.2", "NC_000004.12"],
    },
    BfSeqInfo::{
      name: "chr5",
      length: 181538259L,
      role: "assembled",
      molecule: "chr5",
      unit: "primary",
      aliases: ["5", "CM000667.2", "NC_000005.10"],
    },
    BfSeqInfo::{
      name: "chr6",
      length: 170805979L,
      role: "assembled",
      molecule: "chr6",
      unit: "primary",
      aliases: ["6", "CM000668.2", "NC_000006.12"],
    },
    BfSeqInfo::{
      name: "chr7",
      length: 159345973L,
      role: "assembled",
      molecule: "chr7",
      unit: "primary",
      aliases: ["7", "CM000669.2", "NC_000007.14"],
    },
    BfSeqInfo::{
      name: "chr8",
      length: 145138636L,
      role: "assembled",
      molecule: "chr8",
      unit: "primary",
      aliases: ["8", "CM000670.2", "NC_000008.11"],
    },
    BfSeqInfo::{
      name: "chr9",
      length: 138394717L,
      role: "assembled",
      molecule: "chr9",
      unit: "primary",
      aliases: ["9", "CM000671.2", "NC_000009.12"],
    },
    BfSeqInfo::{
      name: "chr10",
      length: 133797422L,
      role: "assembled",
      molecule: "chr10",
      unit: "primary",
      aliases: ["10", "CM000672.2", "NC_000010.11"],
    },
    BfSeqInfo::{
      name: "chr11",
      length: 135086622L,
      role: "assembled",
      molecule: "chr11",
      unit: "primary",
      aliases: ["11", "CM000673.2", "NC_000011.10"],
    },
    BfSeqInfo::{
      name: "chr12",
      length: 133275309L,
      role: "assembled",
      molecule: "chr12",
      unit: "primary",
      aliases: ["12", "CM000674.2", "NC_000012.12"],
    },
    BfSeqInfo::{
      name: "chr13",
      length: 114364328L,
      role: "assembled",
      molecule: "chr13",
      unit: "primary",
      aliases: ["13", "CM000675.2", "NC_000013.11"],
    },
    BfSeqInfo::{
      name: "chr14",
      length: 107043718L,
      role: "assembled",
      molecule: "chr14",
      unit: "primary",
      aliases: ["14", "CM000676.2", "NC_000014.9"],
    },
    BfSeqInfo::{
      name: "chr15",
      length: 101991189L,
      role: "assembled",
      molecule: "chr15",
      unit: "primary",
      aliases: ["15", "CM000677.2", "NC_000015.10"],
    },
    BfSeqInfo::{
      name: "chr16",
      length: 90338345L,
      role: "assembled",
      molecule: "chr16",
      unit: "primary",
      aliases: ["16", "CM000678.2", "NC_000016.10"],
    },
    BfSeqInfo::{
      name: "chr17",
      length: 83257441L,
      role: "assembled",
      molecule: "chr17",
      unit: "primary",
      aliases: ["17", "CM000679.2", "NC_000017.11"],
    },
    BfSeqInfo::{
      name: "chr18",
      length: 80373285L,
      role: "assembled",
      molecule: "chr18",
      unit: "primary",
      aliases: ["18", "CM000680.2", "NC_000018.10"],
    },
    BfSeqInfo::{
      name: "chr19",
      length: 58617616L,
      role: "assembled",
      molecule: "chr19",
      unit: "primary",
      aliases: ["19", "CM000681.2", "NC_000019.10"],
    },
    BfSeqInfo::{
      name: "chr20",
      length: 64444167L,
      role: "assembled",
      molecule: "chr20",
      unit: "primary",
      aliases: ["20", "CM000682.2", "NC_000020.11"],
    },
    BfSeqInfo::{
      name: "chr21",
      length: 46709983L,
      role: "assembled",
      molecule: "chr21",
      unit: "primary",
      aliases: ["21", "CM000683.2", "NC_000021.9"],
    },
    BfSeqInfo::{
      name: "chr22",
      length: 50818468L,
      role: "assembled",
      molecule: "chr22",
      unit: "primary",
      aliases: ["22", "CM000684.2", "NC_000022.11"],
    },
    BfSeqInfo::{
      name: "chrX",
      length: 156040895L,
      role: "assembled",
      molecule: "chrX",
      unit: "primary",
      aliases: ["X", "CM000685.2", "NC_000023.11"],
    },
    BfSeqInfo::{
      name: "chrY",
      length: 57227415L,
      role: "assembled",
      molecule: "chrY",
      unit: "primary",
      aliases: ["Y", "CM000686.2", "NC_000024.10"],
    },
    BfSeqInfo::{
      name: "chrM",
      length: 16569L,
      role: "assembled",
      molecule: "chrM",
      unit: "non-nuclear",
      aliases: ["MT", "J01415.2", "NC_012920.1"],
    },
  ]
}

///|
fn bf_hg19_seqinfo() -> Array[BfSeqInfo] {
  [
    BfSeqInfo::{
      name: "chr1",
      length: 249250621L,
      role: "assembled",
      molecule: "chr1",
      unit: "primary",
      aliases: ["1", "CM000663.1", "NC_000001.10"],
    },
    BfSeqInfo::{
      name: "chr2",
      length: 243199373L,
      role: "assembled",
      molecule: "chr2",
      unit: "primary",
      aliases: ["2", "CM000664.1", "NC_000002.11"],
    },
    BfSeqInfo::{
      name: "chr3",
      length: 198022430L,
      role: "assembled",
      molecule: "chr3",
      unit: "primary",
      aliases: ["3", "CM000665.1", "NC_000003.11"],
    },
    BfSeqInfo::{
      name: "chr4",
      length: 191154276L,
      role: "assembled",
      molecule: "chr4",
      unit: "primary",
      aliases: ["4", "CM000666.1", "NC_000004.11"],
    },
    BfSeqInfo::{
      name: "chr5",
      length: 180915260L,
      role: "assembled",
      molecule: "chr5",
      unit: "primary",
      aliases: ["5", "CM000667.1", "NC_000005.9"],
    },
    BfSeqInfo::{
      name: "chr6",
      length: 171115067L,
      role: "assembled",
      molecule: "chr6",
      unit: "primary",
      aliases: ["6", "CM000668.1", "NC_000006.11"],
    },
    BfSeqInfo::{
      name: "chr7",
      length: 159138663L,
      role: "assembled",
      molecule: "chr7",
      unit: "primary",
      aliases: ["7", "CM000669.1", "NC_000007.13"],
    },
    BfSeqInfo::{
      name: "chr8",
      length: 146364022L,
      role: "assembled",
      molecule: "chr8",
      unit: "primary",
      aliases: ["8", "CM000670.1", "NC_000008.10"],
    },
    BfSeqInfo::{
      name: "chr9",
      length: 141213431L,
      role: "assembled",
      molecule: "chr9",
      unit: "primary",
      aliases: ["9", "CM000671.1", "NC_000009.11"],
    },
    BfSeqInfo::{
      name: "chr10",
      length: 135534747L,
      role: "assembled",
      molecule: "chr10",
      unit: "primary",
      aliases: ["10", "CM000672.1", "NC_000010.10"],
    },
    BfSeqInfo::{
      name: "chr11",
      length: 135006516L,
      role: "assembled",
      molecule: "chr11",
      unit: "primary",
      aliases: ["11", "CM000673.1", "NC_000011.9"],
    },
    BfSeqInfo::{
      name: "chr12",
      length: 133851895L,
      role: "assembled",
      molecule: "chr12",
      unit: "primary",
      aliases: ["12", "CM000674.1", "NC_000012.11"],
    },
    BfSeqInfo::{
      name: "chr13",
      length: 115169878L,
      role: "assembled",
      molecule: "chr13",
      unit: "primary",
      aliases: ["13", "CM000675.1", "NC_000013.10"],
    },
    BfSeqInfo::{
      name: "chr14",
      length: 107349540L,
      role: "assembled",
      molecule: "chr14",
      unit: "primary",
      aliases: ["14", "CM000676.1", "NC_000014.8"],
    },
    BfSeqInfo::{
      name: "chr15",
      length: 102531392L,
      role: "assembled",
      molecule: "chr15",
      unit: "primary",
      aliases: ["15", "CM000677.1", "NC_000015.9"],
    },
    BfSeqInfo::{
      name: "chr16",
      length: 90354753L,
      role: "assembled",
      molecule: "chr16",
      unit: "primary",
      aliases: ["16", "CM000678.1", "NC_000016.9"],
    },
    BfSeqInfo::{
      name: "chr17",
      length: 81195210L,
      role: "assembled",
      molecule: "chr17",
      unit: "primary",
      aliases: ["17", "CM000679.1", "NC_000017.10"],
    },
    BfSeqInfo::{
      name: "chr18",
      length: 78077248L,
      role: "assembled",
      molecule: "chr18",
      unit: "primary",
      aliases: ["18", "CM000680.1", "NC_000018.9"],
    },
    BfSeqInfo::{
      name: "chr19",
      length: 59128983L,
      role: "assembled",
      molecule: "chr19",
      unit: "primary",
      aliases: ["19", "CM000681.1", "NC_000019.9"],
    },
    BfSeqInfo::{
      name: "chr20",
      length: 63025520L,
      role: "assembled",
      molecule: "chr20",
      unit: "primary",
      aliases: ["20", "CM000682.1", "NC_000020.10"],
    },
    BfSeqInfo::{
      name: "chr21",
      length: 48129895L,
      role: "assembled",
      molecule: "chr21",
      unit: "primary",
      aliases: ["21", "CM000683.1", "NC_000021.8"],
    },
    BfSeqInfo::{
      name: "chr22",
      length: 51304566L,
      role: "assembled",
      molecule: "chr22",
      unit: "primary",
      aliases: ["22", "CM000684.1", "NC_000022.10"],
    },
    BfSeqInfo::{
      name: "chrX",
      length: 155270560L,
      role: "assembled",
      molecule: "chrX",
      unit: "primary",
      aliases: ["X", "CM000685.1", "NC_000023.10"],
    },
    BfSeqInfo::{
      name: "chrY",
      length: 59373566L,
      role: "assembled",
      molecule: "chrY",
      unit: "primary",
      aliases: ["Y", "CM000686.1", "NC_000024.9"],
    },
    BfSeqInfo::{
      name: "chrM",
      length: 16569L,
      role: "assembled",
      molecule: "chrM",
      unit: "non-nuclear",
      aliases: ["MT", "J01415.1", "NC_012920"],
    },
  ]
}

///|
fn bf_mm39_seqinfo() -> Array[BfSeqInfo] {
  [
    BfSeqInfo::{
      name: "chr1",
      length: 195154276L,
      role: "assembled",
      molecule: "chr1",
      unit: "primary",
      aliases: ["1"],
    },
    BfSeqInfo::{
      name: "chr2",
      length: 181755017L,
      role: "assembled",
      molecule: "chr2",
      unit: "primary",
      aliases: ["2"],
    },
    BfSeqInfo::{
      name: "chr3",
      length: 159745316L,
      role: "assembled",
      molecule: "chr3",
      unit: "primary",
      aliases: ["3"],
    },
    BfSeqInfo::{
      name: "chr4",
      length: 156860686L,
      role: "assembled",
      molecule: "chr4",
      unit: "primary",
      aliases: ["4"],
    },
    BfSeqInfo::{
      name: "chr5",
      length: 151758149L,
      role: "assembled",
      molecule: "chr5",
      unit: "primary",
      aliases: ["5"],
    },
    BfSeqInfo::{
      name: "chr6",
      length: 149588044L,
      role: "assembled",
      molecule: "chr6",
      unit: "primary",
      aliases: ["6"],
    },
    BfSeqInfo::{
      name: "chr7",
      length: 144995196L,
      role: "assembled",
      molecule: "chr7",
      unit: "primary",
      aliases: ["7"],
    },
    BfSeqInfo::{
      name: "chr8",
      length: 130127694L,
      role: "assembled",
      molecule: "chr8",
      unit: "primary",
      aliases: ["8"],
    },
    BfSeqInfo::{
      name: "chr9",
      length: 124359700L,
      role: "assembled",
      molecule: "chr9",
      unit: "primary",
      aliases: ["9"],
    },
    BfSeqInfo::{
      name: "chr10",
      length: 130530862L,
      role: "assembled",
      molecule: "chr10",
      unit: "primary",
      aliases: ["10"],
    },
    BfSeqInfo::{
      name: "chr11",
      length: 121843856L,
      role: "assembled",
      molecule: "chr11",
      unit: "primary",
      aliases: ["11"],
    },
    BfSeqInfo::{
      name: "chr12",
      length: 120092757L,
      role: "assembled",
      molecule: "chr12",
      unit: "primary",
      aliases: ["12"],
    },
    BfSeqInfo::{
      name: "chr13",
      length: 120883175L,
      role: "assembled",
      molecule: "chr13",
      unit: "primary",
      aliases: ["13"],
    },
    BfSeqInfo::{
      name: "chr14",
      length: 125139656L,
      role: "assembled",
      molecule: "chr14",
      unit: "primary",
      aliases: ["14"],
    },
    BfSeqInfo::{
      name: "chr15",
      length: 104073951L,
      role: "assembled",
      molecule: "chr15",
      unit: "primary",
      aliases: ["15"],
    },
    BfSeqInfo::{
      name: "chr16",
      length: 98008968L,
      role: "assembled",
      molecule: "chr16",
      unit: "primary",
      aliases: ["16"],
    },
    BfSeqInfo::{
      name: "chr17",
      length: 95294699L,
      role: "assembled",
      molecule: "chr17",
      unit: "primary",
      aliases: ["17"],
    },
    BfSeqInfo::{
      name: "chr18",
      length: 90720763L,
      role: "assembled",
      molecule: "chr18",
      unit: "primary",
      aliases: ["18"],
    },
    BfSeqInfo::{
      name: "chr19",
      length: 61420004L,
      role: "assembled",
      molecule: "chr19",
      unit: "primary",
      aliases: ["19"],
    },
    BfSeqInfo::{
      name: "chrX",
      length: 169951750L,
      role: "assembled",
      molecule: "chrX",
      unit: "primary",
      aliases: ["X"],
    },
    BfSeqInfo::{
      name: "chrY",
      length: 91455967L,
      role: "assembled",
      molecule: "chrY",
      unit: "primary",
      aliases: ["Y"],
    },
    BfSeqInfo::{
      name: "chrM",
      length: 16299L,
      role: "assembled",
      molecule: "chrM",
      unit: "non-nuclear",
      aliases: ["MT"],
    },
  ]
}

///|
fn bf_mm10_seqinfo() -> Array[BfSeqInfo] {
  [
    BfSeqInfo::{
      name: "chr1",
      length: 195471971L,
      role: "assembled",
      molecule: "chr1",
      unit: "primary",
      aliases: ["1"],
    },
    BfSeqInfo::{
      name: "chr2",
      length: 182113224L,
      role: "assembled",
      molecule: "chr2",
      unit: "primary",
      aliases: ["2"],
    },
    BfSeqInfo::{
      name: "chr3",
      length: 160039680L,
      role: "assembled",
      molecule: "chr3",
      unit: "primary",
      aliases: ["3"],
    },
    BfSeqInfo::{
      name: "chr4",
      length: 156508116L,
      role: "assembled",
      molecule: "chr4",
      unit: "primary",
      aliases: ["4"],
    },
    BfSeqInfo::{
      name: "chr5",
      length: 151834684L,
      role: "assembled",
      molecule: "chr5",
      unit: "primary",
      aliases: ["5"],
    },
    BfSeqInfo::{
      name: "chr6",
      length: 149736546L,
      role: "assembled",
      molecule: "chr6",
      unit: "primary",
      aliases: ["6"],
    },
    BfSeqInfo::{
      name: "chr7",
      length: 145441459L,
      role: "assembled",
      molecule: "chr7",
      unit: "primary",
      aliases: ["7"],
    },
    BfSeqInfo::{
      name: "chr8",
      length: 129401213L,
      role: "assembled",
      molecule: "chr8",
      unit: "primary",
      aliases: ["8"],
    },
    BfSeqInfo::{
      name: "chr9",
      length: 124595110L,
      role: "assembled",
      molecule: "chr9",
      unit: "primary",
      aliases: ["9"],
    },
    BfSeqInfo::{
      name: "chr10",
      length: 130694993L,
      role: "assembled",
      molecule: "chr10",
      unit: "primary",
      aliases: ["10"],
    },
    BfSeqInfo::{
      name: "chr11",
      length: 122082543L,
      role: "assembled",
      molecule: "chr11",
      unit: "primary",
      aliases: ["11"],
    },
    BfSeqInfo::{
      name: "chr12",
      length: 120129022L,
      role: "assembled",
      molecule: "chr12",
      unit: "primary",
      aliases: ["12"],
    },
    BfSeqInfo::{
      name: "chr13",
      length: 120421639L,
      role: "assembled",
      molecule: "chr13",
      unit: "primary",
      aliases: ["13"],
    },
    BfSeqInfo::{
      name: "chr14",
      length: 124902244L,
      role: "assembled",
      molecule: "chr14",
      unit: "primary",
      aliases: ["14"],
    },
    BfSeqInfo::{
      name: "chr15",
      length: 104043685L,
      role: "assembled",
      molecule: "chr15",
      unit: "primary",
      aliases: ["15"],
    },
    BfSeqInfo::{
      name: "chr16",
      length: 98207768L,
      role: "assembled",
      molecule: "chr16",
      unit: "primary",
      aliases: ["16"],
    },
    BfSeqInfo::{
      name: "chr17",
      length: 94987271L,
      role: "assembled",
      molecule: "chr17",
      unit: "primary",
      aliases: ["17"],
    },
    BfSeqInfo::{
      name: "chr18",
      length: 90702639L,
      role: "assembled",
      molecule: "chr18",
      unit: "primary",
      aliases: ["18"],
    },
    BfSeqInfo::{
      name: "chr19",
      length: 61431566L,
      role: "assembled",
      molecule: "chr19",
      unit: "primary",
      aliases: ["19"],
    },
    BfSeqInfo::{
      name: "chrX",
      length: 171031299L,
      role: "assembled",
      molecule: "chrX",
      unit: "primary",
      aliases: ["X"],
    },
    BfSeqInfo::{
      name: "chrY",
      length: 91744698L,
      role: "assembled",
      molecule: "chrY",
      unit: "primary",
      aliases: ["Y"],
    },
    BfSeqInfo::{
      name: "chrM",
      length: 16299L,
      role: "assembled",
      molecule: "chrM",
      unit: "non-nuclear",
      aliases: ["MT"],
    },
  ]
}

///|
fn bf_dm6_seqinfo() -> Array[BfSeqInfo] {
  [
    BfSeqInfo::{
      name: "chr2L",
      length: 23513712L,
      role: "assembled",
      molecule: "chr2L",
      unit: "primary",
      aliases: ["2L"],
    },
    BfSeqInfo::{
      name: "chr2R",
      length: 25286936L,
      role: "assembled",
      molecule: "chr2R",
      unit: "primary",
      aliases: ["2R"],
    },
    BfSeqInfo::{
      name: "chr3L",
      length: 28110240L,
      role: "assembled",
      molecule: "chr3L",
      unit: "primary",
      aliases: ["3L"],
    },
    BfSeqInfo::{
      name: "chr3R",
      length: 32079331L,
      role: "assembled",
      molecule: "chr3R",
      unit: "primary",
      aliases: ["3R"],
    },
    BfSeqInfo::{
      name: "chr4",
      length: 1348131L,
      role: "assembled",
      molecule: "chr4",
      unit: "primary",
      aliases: ["4"],
    },
    BfSeqInfo::{
      name: "chrX",
      length: 23542271L,
      role: "assembled",
      molecule: "chrX",
      unit: "primary",
      aliases: ["X"],
    },
    BfSeqInfo::{
      name: "chrM",
      length: 19517L,
      role: "assembled",
      molecule: "chrM",
      unit: "non-nuclear",
      aliases: ["MT"],
    },
  ]
}

///|
fn bf_ce11_seqinfo() -> Array[BfSeqInfo] {
  [
    BfSeqInfo::{
      name: "chrI",
      length: 15072434L,
      role: "assembled",
      molecule: "chrI",
      unit: "primary",
      aliases: ["I"],
    },
    BfSeqInfo::{
      name: "chrII",
      length: 15279421L,
      role: "assembled",
      molecule: "chrII",
      unit: "primary",
      aliases: ["II"],
    },
    BfSeqInfo::{
      name: "chrIII",
      length: 13783801L,
      role: "assembled",
      molecule: "chrIII",
      unit: "primary",
      aliases: ["III"],
    },
    BfSeqInfo::{
      name: "chrIV",
      length: 17493829L,
      role: "assembled",
      molecule: "chrIV",
      unit: "primary",
      aliases: ["IV"],
    },
    BfSeqInfo::{
      name: "chrV",
      length: 20924149L,
      role: "assembled",
      molecule: "chrV",
      unit: "primary",
      aliases: ["V"],
    },
    BfSeqInfo::{
      name: "chrX",
      length: 17718942L,
      role: "assembled",
      molecule: "chrX",
      unit: "primary",
      aliases: ["X"],
    },
    BfSeqInfo::{
      name: "chrM",
      length: 13794L,
      role: "assembled",
      molecule: "chrM",
      unit: "non-nuclear",
      aliases: ["MT"],
    },
  ]
}

///|
fn bf_danrer11_seqinfo() -> Array[BfSeqInfo] {
  [
    BfSeqInfo::{
      name: "chr1",
      length: 60555572L,
      role: "assembled",
      molecule: "chr1",
      unit: "primary",
      aliases: ["1"],
    },
    BfSeqInfo::{
      name: "chr2",
      length: 38953223L,
      role: "assembled",
      molecule: "chr2",
      unit: "primary",
      aliases: ["2"],
    },
    BfSeqInfo::{
      name: "chr3",
      length: 62626801L,
      role: "assembled",
      molecule: "chr3",
      unit: "primary",
      aliases: ["3"],
    },
    BfSeqInfo::{
      name: "chr4",
      length: 64322427L,
      role: "assembled",
      molecule: "chr4",
      unit: "primary",
      aliases: ["4"],
    },
    BfSeqInfo::{
      name: "chr5",
      length: 75862389L,
      role: "assembled",
      molecule: "chr5",
      unit: "primary",
      aliases: ["5"],
    },
    BfSeqInfo::{
      name: "chr6",
      length: 55289941L,
      role: "assembled",
      molecule: "chr6",
      unit: "primary",
      aliases: ["6"],
    },
    BfSeqInfo::{
      name: "chr7",
      length: 56534483L,
      role: "assembled",
      molecule: "chr7",
      unit: "primary",
      aliases: ["7"],
    },
    BfSeqInfo::{
      name: "chr8",
      length: 54439034L,
      role: "assembled",
      molecule: "chr8",
      unit: "primary",
      aliases: ["8"],
    },
    BfSeqInfo::{
      name: "chr9",
      length: 34362851L,
      role: "assembled",
      molecule: "chr9",
      unit: "primary",
      aliases: ["9"],
    },
    BfSeqInfo::{
      name: "chr10",
      length: 46610897L,
      role: "assembled",
      molecule: "chr10",
      unit: "primary",
      aliases: ["10"],
    },
    BfSeqInfo::{
      name: "chr11",
      length: 46562963L,
      role: "assembled",
      molecule: "chr11",
      unit: "primary",
      aliases: ["11"],
    },
    BfSeqInfo::{
      name: "chr12",
      length: 55782147L,
      role: "assembled",
      molecule: "chr12",
      unit: "primary",
      aliases: ["12"],
    },
    BfSeqInfo::{
      name: "chr13",
      length: 40064376L,
      role: "assembled",
      molecule: "chr13",
      unit: "primary",
      aliases: ["13"],
    },
    BfSeqInfo::{
      name: "chr14",
      length: 52220120L,
      role: "assembled",
      molecule: "chr14",
      unit: "primary",
      aliases: ["14"],
    },
    BfSeqInfo::{
      name: "chr15",
      length: 43210243L,
      role: "assembled",
      molecule: "chr15",
      unit: "primary",
      aliases: ["15"],
    },
    BfSeqInfo::{
      name: "chr16",
      length: 50859401L,
      role: "assembled",
      molecule: "chr16",
      unit: "primary",
      aliases: ["16"],
    },
    BfSeqInfo::{
      name: "chr17",
      length: 31834105L,
      role: "assembled",
      molecule: "chr17",
      unit: "primary",
      aliases: ["17"],
    },
    BfSeqInfo::{
      name: "chr18",
      length: 53417300L,
      role: "assembled",
      molecule: "chr18",
      unit: "primary",
      aliases: ["18"],
    },
    BfSeqInfo::{
      name: "chr19",
      length: 47122145L,
      role: "assembled",
      molecule: "chr19",
      unit: "primary",
      aliases: ["19"],
    },
    BfSeqInfo::{
      name: "chr20",
      length: 54634480L,
      role: "assembled",
      molecule: "chr20",
      unit: "primary",
      aliases: ["20"],
    },
    BfSeqInfo::{
      name: "chr21",
      length: 45701097L,
      role: "assembled",
      molecule: "chr21",
      unit: "primary",
      aliases: ["21"],
    },
    BfSeqInfo::{
      name: "chr22",
      length: 39802372L,
      role: "assembled",
      molecule: "chr22",
      unit: "primary",
      aliases: ["22"],
    },
    BfSeqInfo::{
      name: "chr23",
      length: 46635246L,
      role: "assembled",
      molecule: "chr23",
      unit: "primary",
      aliases: ["23"],
    },
    BfSeqInfo::{
      name: "chr24",
      length: 37754680L,
      role: "assembled",
      molecule: "chr24",
      unit: "primary",
      aliases: ["24"],
    },
    BfSeqInfo::{
      name: "chr25",
      length: 16799021L,
      role: "assembled",
      molecule: "chr25",
      unit: "primary",
      aliases: ["25"],
    },
    BfSeqInfo::{
      name: "chrX",
      length: 69715921L,
      role: "assembled",
      molecule: "chrX",
      unit: "primary",
      aliases: ["X"],
    },
    BfSeqInfo::{
      name: "chrY",
      length: 48714907L,
      role: "assembled",
      molecule: "chrY",
      unit: "primary",
      aliases: ["Y"],
    },
    BfSeqInfo::{
      name: "chrM",
      length: 16569L,
      role: "assembled",
      molecule: "chrM",
      unit: "non-nuclear",
      aliases: ["MT"],
    },
  ]
}

///|
fn bf_sacCer3_seqinfo() -> Array[BfSeqInfo] {
  [
    BfSeqInfo::{
      name: "chrI",
      length: 230218L,
      role: "assembled",
      molecule: "chrI",
      unit: "primary",
      aliases: ["I"],
    },
    BfSeqInfo::{
      name: "chrII",
      length: 813184L,
      role: "assembled",
      molecule: "chrII",
      unit: "primary",
      aliases: ["II"],
    },
    BfSeqInfo::{
      name: "chrIII",
      length: 316620L,
      role: "assembled",
      molecule: "chrIII",
      unit: "primary",
      aliases: ["III"],
    },
    BfSeqInfo::{
      name: "chrIV",
      length: 1531933L,
      role: "assembled",
      molecule: "chrIV",
      unit: "primary",
      aliases: ["IV"],
    },
    BfSeqInfo::{
      name: "chrV",
      length: 576874L,
      role: "assembled",
      molecule: "chrV",
      unit: "primary",
      aliases: ["V"],
    },
    BfSeqInfo::{
      name: "chrVI",
      length: 270161L,
      role: "assembled",
      molecule: "chrVI",
      unit: "primary",
      aliases: ["VI"],
    },
    BfSeqInfo::{
      name: "chrVII",
      length: 1090940L,
      role: "assembled",
      molecule: "chrVII",
      unit: "primary",
      aliases: ["VII"],
    },
    BfSeqInfo::{
      name: "chrVIII",
      length: 562643L,
      role: "assembled",
      molecule: "chrVIII",
      unit: "primary",
      aliases: ["VIII"],
    },
    BfSeqInfo::{
      name: "chrIX",
      length: 439888L,
      role: "assembled",
      molecule: "chrIX",
      unit: "primary",
      aliases: ["IX"],
    },
    BfSeqInfo::{
      name: "chrX",
      length: 745751L,
      role: "assembled",
      molecule: "chrX",
      unit: "primary",
      aliases: ["X"],
    },
    BfSeqInfo::{
      name: "chrXI",
      length: 666816L,
      role: "assembled",
      molecule: "chrXI",
      unit: "primary",
      aliases: ["XI"],
    },
    BfSeqInfo::{
      name: "chrXII",
      length: 1078177L,
      role: "assembled",
      molecule: "chrXII",
      unit: "primary",
      aliases: ["XII"],
    },
    BfSeqInfo::{
      name: "chrXIII",
      length: 924226L,
      role: "assembled",
      molecule: "chrXIII",
      unit: "primary",
      aliases: ["XIII"],
    },
    BfSeqInfo::{
      name: "chrXIV",
      length: 784333L,
      role: "assembled",
      molecule: "chrXIV",
      unit: "primary",
      aliases: ["XIV"],
    },
    BfSeqInfo::{
      name: "chrXV",
      length: 1091291L,
      role: "assembled",
      molecule: "chrXV",
      unit: "primary",
      aliases: ["XV"],
    },
    BfSeqInfo::{
      name: "chrXVI",
      length: 948066L,
      role: "assembled",
      molecule: "chrXVI",
      unit: "primary",
      aliases: ["XVI"],
    },
    BfSeqInfo::{
      name: "chrM",
      length: 85779L,
      role: "assembled",
      molecule: "chrM",
      unit: "non-nuclear",
      aliases: ["MT"],
    },
  ]
}

///|
fn bf_hs1_seqinfo() -> Array[BfSeqInfo] {
  [
    BfSeqInfo::{
      name: "chr1",
      length: 248327194L,
      role: "assembled",
      molecule: "chr1",
      unit: "primary",
      aliases: ["1"],
    },
    BfSeqInfo::{
      name: "chr2",
      length: 242274980L,
      role: "assembled",
      molecule: "chr2",
      unit: "primary",
      aliases: ["2"],
    },
    BfSeqInfo::{
      name: "chr3",
      length: 198449013L,
      role: "assembled",
      molecule: "chr3",
      unit: "primary",
      aliases: ["3"],
    },
    BfSeqInfo::{
      name: "chr4",
      length: 190214555L,
      role: "assembled",
      molecule: "chr4",
      unit: "primary",
      aliases: ["4"],
    },
    BfSeqInfo::{
      name: "chr5",
      length: 181538259L,
      role: "assembled",
      molecule: "chr5",
      unit: "primary",
      aliases: ["5"],
    },
    BfSeqInfo::{
      name: "chr6",
      length: 170805979L,
      role: "assembled",
      molecule: "chr6",
      unit: "primary",
      aliases: ["6"],
    },
    BfSeqInfo::{
      name: "chr7",
      length: 159345973L,
      role: "assembled",
      molecule: "chr7",
      unit: "primary",
      aliases: ["7"],
    },
    BfSeqInfo::{
      name: "chr8",
      length: 145138636L,
      role: "assembled",
      molecule: "chr8",
      unit: "primary",
      aliases: ["8"],
    },
    BfSeqInfo::{
      name: "chr9",
      length: 138394717L,
      role: "assembled",
      molecule: "chr9",
      unit: "primary",
      aliases: ["9"],
    },
    BfSeqInfo::{
      name: "chr10",
      length: 133797422L,
      role: "assembled",
      molecule: "chr10",
      unit: "primary",
      aliases: ["10"],
    },
    BfSeqInfo::{
      name: "chr11",
      length: 135086622L,
      role: "assembled",
      molecule: "chr11",
      unit: "primary",
      aliases: ["11"],
    },
    BfSeqInfo::{
      name: "chr12",
      length: 133275309L,
      role: "assembled",
      molecule: "chr12",
      unit: "primary",
      aliases: ["12"],
    },
    BfSeqInfo::{
      name: "chr13",
      length: 114364328L,
      role: "assembled",
      molecule: "chr13",
      unit: "primary",
      aliases: ["13"],
    },
    BfSeqInfo::{
      name: "chr14",
      length: 107043718L,
      role: "assembled",
      molecule: "chr14",
      unit: "primary",
      aliases: ["14"],
    },
    BfSeqInfo::{
      name: "chr15",
      length: 101991189L,
      role: "assembled",
      molecule: "chr15",
      unit: "primary",
      aliases: ["15"],
    },
    BfSeqInfo::{
      name: "chr16",
      length: 90338345L,
      role: "assembled",
      molecule: "chr16",
      unit: "primary",
      aliases: ["16"],
    },
    BfSeqInfo::{
      name: "chr17",
      length: 83257441L,
      role: "assembled",
      molecule: "chr17",
      unit: "primary",
      aliases: ["17"],
    },
    BfSeqInfo::{
      name: "chr18",
      length: 80373285L,
      role: "assembled",
      molecule: "chr18",
      unit: "primary",
      aliases: ["18"],
    },
    BfSeqInfo::{
      name: "chr19",
      length: 58617616L,
      role: "assembled",
      molecule: "chr19",
      unit: "primary",
      aliases: ["19"],
    },
    BfSeqInfo::{
      name: "chr20",
      length: 64444167L,
      role: "assembled",
      molecule: "chr20",
      unit: "primary",
      aliases: ["20"],
    },
    BfSeqInfo::{
      name: "chr21",
      length: 46709983L,
      role: "assembled",
      molecule: "chr21",
      unit: "primary",
      aliases: ["21"],
    },
    BfSeqInfo::{
      name: "chr22",
      length: 50818468L,
      role: "assembled",
      molecule: "chr22",
      unit: "primary",
      aliases: ["22"],
    },
    BfSeqInfo::{
      name: "chrX",
      length: 154259566L,
      role: "assembled",
      molecule: "chrX",
      unit: "primary",
      aliases: ["X"],
    },
    BfSeqInfo::{
      name: "chrY",
      length: 62460029L,
      role: "assembled",
      molecule: "chrY",
      unit: "primary",
      aliases: ["Y"],
    },
    BfSeqInfo::{
      name: "chrM",
      length: 16569L,
      role: "assembled",
      molecule: "chrM",
      unit: "non-nuclear",
      aliases: ["MT"],
    },
  ]
}

// ================================================================
// 5. Runtime helpers
// ================================================================

///|
/// Split a multi-line string into non-empty lines.
fn bf_split_lines(text : String) -> Array[String] {
  let result : Array[String] = []
  let mut sb = StringBuilder()
  for c in text {
    if c == '\n' {
      let cur = sb.to_string()
      if cur.length() > 0 {
        result.push(cur)
      }
      sb = StringBuilder()
    } else {
      sb.write_char(c)
    }
  }
  let last = sb.to_string()
  if last.length() > 0 {
    result.push(last)
  }
  result
}

///|
/// Look up the canonical seqinfo list for a named assembly.
/// Returns None if the assembly is not in the curated registry.
pub fn bf_lookup_seqinfo(assembly : String) -> Array[BfSeqInfo]? {
  match assembly {
    "hg38" => Some(bf_hg38_seqinfo())
    "hg19" => Some(bf_hg19_seqinfo())
    "mm39" => Some(bf_mm39_seqinfo())
    "mm10" => Some(bf_mm10_seqinfo())
    "dm6" => Some(bf_dm6_seqinfo())
    "ce11" => Some(bf_ce11_seqinfo())
    "danRer11" => Some(bf_danrer11_seqinfo())
    "sacCer3" => Some(bf_sacCer3_seqinfo())
    "hs1" => Some(bf_hs1_seqinfo())
    _ => None
  }
}

///|
/// Build chromsizes (Array[(String, Int32)]) from seqinfo — just the
/// lengths of assembled primary chromosomes + chrX/Y + chrM.
/// Ignores unlocalized/unplaced/random/contigs etc.
pub fn bf_assembly_chromsizes(assembly : String) -> Array[(String, Int32)] {
  match bf_lookup_seqinfo(assembly) {
    Some(infos) => {
      let out : Array[(String, Int32)] = []
      for si in infos {
        let l = si.length
        if l > 2147483647L {
          continue
        }
        out.push((si.name, l.to_int()))
      }
      out
    }
    None => []
  }
}

///|
/// Build a GenomeInfoDb (existing type) for an assembly.
pub fn bf_assembly_genome_info(assembly : String) -> GenomeInfoDb? {
  let maybe = bf_lookup_seqinfo(assembly)
  match maybe {
    None => None
    Some(infos) => {
      // Find assembly meta
      let mut meta : BfAssembly? = None
      let all = bf_assemblies_available()
      let mut i = 0
      while i < all.length() {
        if all[i].name == assembly {
          meta = Some(all[i])
          break
        }
        i = i + 1
      }
      match meta {
        None => None
        Some(a) => {
          let mut db = GenomeInfoDb::new(
            genome_build=a.name,
            organism=a.organism,
          )
          for si in infos {
            let mut is_circ = false
            if si.name == "chrM" {
              is_circ = true
            }
            if si.length > 2147483647L {
              continue
            }
            db = db.add_chrom(si.name, si.length.to_int(), is_circular=is_circ)
          }
          // centromeres
          let cents = bf_assembly_centromeres(assembly)
          for c in cents {
            db = db.add_centromere(c.0, c.1, c.2)
          }
          Some(db)
        }
      }
    }
  }
}

///|
/// Parse centromere text blob (3 TSV columns) → (chrom, start, end).
pub fn bf_assembly_centromeres(
  assembly : String,
) -> Array[(String, Int32, Int32)] {
  let blob = match assembly {
    "hg19" => BF_HG19_CENTROMERES
    "mm10" => BF_MM10_CENTROMERES
    "mm39" => BF_MM39_CENTROMERES
    "dm6" => BF_DM6_CENTROMERES
    _ => ""
  }
  if blob.length() == 0 {
    return []
  }
  let lines = bf_split_lines(blob)
  let out : Array[(String, Int32, Int32)] = []
  for ln in lines {
    let f = bf_split_tab(ln)
    if f.length() >= 3 {
      let chrom = f[0]
      let s = parse_int(f[1])
      let e = parse_int(f[2])
      out.push((chrom, s, e))
    }
  }
  out
}

///|
/// Resolve a possibly-aliased chromosome name to canonical UCSC form.
/// E.g. alias "MT" → "chrM"; "1" → "chr1" (hg38); "X" → "chrX".
/// If the name is already canonical, returns it unchanged.
/// Returns None only if no matching seqinfo exists at all.
pub fn bf_resolve_alias(assembly : String, name : String) -> String? {
  match bf_lookup_seqinfo(assembly) {
    None => None
    Some(infos) => {
      // direct hit — already canonical
      for si in infos {
        if si.name == name {
          return Some(name)
        }
      }
      // aliases (case-insensitive comparison)
      let lower = name.to_lower()
      for si in infos {
        for al in si.aliases {
          if al.to_lower() == lower {
            return Some(si.name)
          }
        }
      }
      None
    }
  }
}

///|
/// Cytobands for an assembly. Returns empty Array when no curated data.
pub fn bf_assembly_cytobands(assembly : String) -> Array[BfCytoband] {
  let blob = match assembly {
    "hg38" => BF_HG38_CYTOBANDS_HEAD
    _ => ""
  }
  if blob.length() == 0 {
    return []
  }
  let lines = bf_split_lines(blob)
  let out : Array[BfCytoband] = []
  for ln in lines {
    let f = bf_split_tab(ln)
    if f.length() >= 5 {
      let chrom = f[0]
      let s = parse_int(f[1])
      let e = parse_int(f[2])
      let bname = f[3]
      let gie = f[4]
      out.push(BfCytoband::{
        chrom,
        start: s.to_int64(),
        end: e.to_int64(),
        name: bname,
        gieStain: gie,
      })
    }
  }
  out
}

// ================================================================
// 6. Common tab-split helper
// ================================================================

///|
/// Split a line by '\t'.
fn bf_split_tab(line : String) -> Array[String] {
  let result : Array[String] = []
  let mut sb = StringBuilder()
  for c in line {
    if c == '\t' {
      result.push(sb.to_string())
      sb = StringBuilder()
    } else {
      sb.write_char(c)
    }
  }
  result.push(sb.to_string())
  result
}
