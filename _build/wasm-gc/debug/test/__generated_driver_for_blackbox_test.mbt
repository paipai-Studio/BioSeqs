// #region Define traits with default implementations

///|
/// is_being_cancelled may be unused if there's no async test
#warnings("-unused_value")
trait MoonBit_Async_Test_Driver {
  fn run_async_tests(@moonbitlang/core/prelude.Array[async () -> Unit]) -> Unit = _
  fn is_being_cancelled() -> Bool = _
}

///|
impl MoonBit_Async_Test_Driver with fn run_async_tests(_) {
  ()
}

///|
impl MoonBit_Async_Test_Driver with fn is_being_cancelled() {
  false
}

///|
trait MoonBit_Test_Driver {
  /// returns if the test is found and executed
  fn run_test(
    @moonbitlang/core/prelude.Array[async () -> Unit],
    String,
    Int,
    () -> Unit,
    (String, String, Bool) -> Unit,
    (Error) -> String,
  ) -> Bool raise MoonBitTestDriverInternalSkipTest = _
}

///|
impl MoonBit_Test_Driver with fn run_test(_, _, _, _, _, _) {
  false
}

// #endregion
// #region Define abstract types for namespace with implementation declarations

///|
priv type MoonBit_Test_Driver_Internal_No_Args

///|
impl MoonBit_Test_Driver for MoonBit_Test_Driver_Internal_No_Args

///|
priv type MoonBit_Test_Driver_Internal_With_Args

///|
impl MoonBit_Test_Driver for MoonBit_Test_Driver_Internal_With_Args

///|
priv type MoonBit_Test_Driver_Internal_Async_No_Args

///|
impl MoonBit_Test_Driver for MoonBit_Test_Driver_Internal_Async_No_Args

///|
priv type MoonBit_Test_Driver_Internal_Async_With_Args

///|
impl MoonBit_Test_Driver for MoonBit_Test_Driver_Internal_Async_With_Args

///|
priv type MoonBit_Test_Driver_Internal_With_Bench_Args

///|
impl MoonBit_Test_Driver for MoonBit_Test_Driver_Internal_With_Bench_Args

///|
priv type MoonBit_Async_Test_Driver_Impl

///|
impl MoonBit_Async_Test_Driver for MoonBit_Async_Test_Driver_Impl

// #endregion

// #region Utility Types

///|
/// Maybe unused if there's no tests
#warnings("-unused_type_declaration")
priv struct MoonBitTestDriverInternalTestMap[F](
  @moonbitlang/core/prelude.Map[
    String,
    @moonbitlang/core/prelude.Map[Int, (F, ReadOnlyArray[String])],
  ]
)

///|
priv suberror MoonBitTestDriverInternalSkipTest {
  MoonBitTestDriverInternalSkipTest(name~ : String)
}

// #endregion
///|
fn moonbit_test_driver_internal_error_to_string(x : Error) -> String = "%error.to_string"

///|
#cfg(any(target="native", target="llvm"))
fn string_get(s : String, index : Int) -> Int = "%string_get"

///|
#cfg(not(target="js"))
fn moonbit_unsafe_char_from_int(x : Int) -> Char = "%identity"

///|
#cfg(any(target="native", target="llvm"))
const MOONBIT_TEST_DRIVER_INTERNAL_IS_NATIVE = true

///|
#cfg(not(any(target="native", target="llvm")))
const MOONBIT_TEST_DRIVER_INTERNAL_IS_NATIVE = false

// =========================================================
// ================== BEGIN WASM SPECIFIC ==================
// =========================================================

///|
#cfg(any(target="wasm", target="wasm-gc"))
fn moonbit_test_driver_internal_get_file_name(
  file_name : MoonbitTestDriverInternalExternString,
) -> String {
  moonbit_test_driver_internal_string_from_extern(file_name)
}

///|
#cfg(any(target="wasm", target="wasm-gc"))
#external
type MoonbitTestDriverInternalStringReadHandle

///|
#cfg(any(target="wasm", target="wasm-gc"))
#external
type MoonbitTestDriverInternalExternString

///|
#cfg(any(target="wasm", target="wasm-gc"))
fn moonbit_test_driver_internal_begin_read_string(
  s : MoonbitTestDriverInternalExternString,
) -> MoonbitTestDriverInternalStringReadHandle = "__moonbit_fs_unstable" "begin_read_string"

///|
#cfg(any(target="wasm", target="wasm-gc"))
fn moonbit_test_driver_internal_string_read_char(
  handle : MoonbitTestDriverInternalStringReadHandle,
) -> Int = "__moonbit_fs_unstable" "string_read_char"

///|
#cfg(any(target="wasm", target="wasm-gc"))
fn moonbit_test_driver_internal_finish_read_string(
  handle : MoonbitTestDriverInternalStringReadHandle,
) = "__moonbit_fs_unstable" "finish_read_string"

///|
#cfg(any(target="wasm", target="wasm-gc"))
fn fixedarray_to_bytes(arr : FixedArray[Byte]) -> Bytes = "%identity"

///|
#cfg(any(target="wasm", target="wasm-gc"))
fn moonbit_test_driver_internal_string_from_extern(
  e : MoonbitTestDriverInternalExternString,
) -> String {
  let buf = FixedArray::make(512, Byte::default())
  let mut len = 0
  let handle = moonbit_test_driver_internal_begin_read_string(e)
  while true {
    let ch = moonbit_test_driver_internal_string_read_char(handle)
    if ch == -1 {
      break
    } else {
      let written_bytes_len = buf.set_utf16le_char(
        len,
        moonbit_unsafe_char_from_int(ch),
      )
      len += written_bytes_len
    }
  }
  moonbit_test_driver_internal_finish_read_string(handle)
  fixedarray_to_bytes(buf).to_unchecked_string(offset=0, length=len).to_string()
}

// =========================================================
// =================== END WASM SPECIFIC ===================
// =========================================================

// =========================================================
// ============== BEGIN NATIVE/LLVM SPECIFIC ===============
// =========================================================

///|
#cfg(any(target="native", target="llvm"))
fn moonbit_test_driver_internal_native_parse_args() -> @moonbitlang/core/prelude.Array[
  (String, Int),
] {
  fn moonbit_test_driver_internal_parse_int_(s : String) -> Int {
    let mut res = 0
    let len = s.length()
    for i in 0..<len {
      res = res * 10 + (string_get(s, i) - '0')
    }
    res
  }

  fn moonbit_test_driver_internal_utf8_bytes_to_mbt_string(
    bytes : Bytes,
  ) -> String {
    let res = @moonbitlang/core/prelude.StringBuilder::new()
    let len = bytes.length()
    let mut i = 0
    while i < len {
      let mut c = bytes[i].to_int()
      if c < 0x80 {
        res.write_char(moonbit_unsafe_char_from_int(c))
        i += 1
      } else if c < 0xE0 {
        if i + 1 >= len {
          break
        }
        c = ((c & 0x1F) << 6) | (bytes[i + 1].to_int() & 0x3F)
        res.write_char(moonbit_unsafe_char_from_int(c))
        i += 2
      } else if c < 0xF0 {
        if i + 2 >= len {
          break
        }
        c = ((c & 0x0F) << 12) |
          ((bytes[i + 1].to_int() & 0x3F) << 6) |
          (bytes[i + 2].to_int() & 0x3F)
        res.write_char(moonbit_unsafe_char_from_int(c))
        i += 3
      } else {
        if i + 3 >= len {
          break
        }
        c = ((c & 0x07) << 18) |
          ((bytes[i + 1].to_int() & 0x3F) << 12) |
          ((bytes[i + 2].to_int() & 0x3F) << 6) |
          (bytes[i + 3].to_int() & 0x3F)
        c -= 0x10000
        res.write_char(moonbit_unsafe_char_from_int((c >> 10) + 0xD800))
        res.write_char(moonbit_unsafe_char_from_int((c & 0x3FF) + 0xDC00))
        i += 4
      }
    }
    res.to_string()
  }

  fn moonbit_test_driver_internal_get_cli_args_internal() -> @moonbitlang/core/prelude.Array[
    String,
  ] {
    let tmp = moonbit_test_driver_internal_get_cli_args_ffi()
    let res = @moonbitlang/core/prelude.Array::new(capacity=tmp.length())
    for i in 0..<tmp.length() {
      res.push(moonbit_test_driver_internal_utf8_bytes_to_mbt_string(tmp[i]))
    }
    res
  }

  fn moonbit_test_driver_internal_split_mbt_string(
    s : String,
    sep : Char,
  ) -> @moonbitlang/core/prelude.Array[String] {
    let res = []
    let mut i = 0
    let mut start = 0
    while i < s.length() {
      if moonbit_unsafe_char_from_int(string_get(s, i)) == sep {
        res.push(s.unsafe_substring(start~, end=i))
        start = i + 1
      }
      i += 1
    }
    if start < s.length() {
      res.push(s.unsafe_substring(start~, end=s.length()))
    }
    res
  }

  let file_and_index = []
  let cli_args = moonbit_test_driver_internal_get_cli_args_internal()
  let test_args = moonbit_test_driver_internal_split_mbt_string(
    cli_args[1],
    '/',
  )
  for arg in test_args {
    let file_and_range = moonbit_test_driver_internal_split_mbt_string(arg, ':')
    let file = file_and_range[0]
    let range = file_and_range[1]
    let start_and_end = moonbit_test_driver_internal_split_mbt_string(
      range, '-',
    )
    let start = moonbit_test_driver_internal_parse_int_(start_and_end[0])
    let end = moonbit_test_driver_internal_parse_int_(start_and_end[1])
    for i in start..<end {
      file_and_index.push((file, i))
    }
  }
  file_and_index
}

///|
#cfg(any(target="native", target="llvm"))
fn moonbit_test_driver_internal_get_cli_args_ffi() -> FixedArray[Bytes] = "$moonbit.get_cli_args"

// =========================================================
// =============== END NATIVE/LLVM SPECIFIC ================
// =========================================================

// =========================================================
// -================== BEGIN JS SPECIFIC ===================
// =========================================================

///|
#cfg(target="js")
#external
type MoonBitTestDriverInternalTestEntry

///|
#cfg(target="js")
extern "js" fn MoonBitTestDriverInternalTestEntry::filename(
  self : Self,
) -> String =
  #| (entry) => entry[0]

///|
#cfg(target="js")
extern "js" fn MoonBitTestDriverInternalTestEntry::index(self : Self) -> Int =
  #| (entry) => entry[1]

///|
/// JavaScript error captured by the generated test driver.
pub(all) suberror MoonBitTestDriverInternalJsError {
  MoonBitTestDriverInternalJsError(String)
}

///|
/// Maybe unused if there's no tests
#cfg(target="js")
#warnings("-unused_value")
fn moonbit_test_driver_internal_catch_error(
  f : () -> Unit raise,
  on_ok~ : () -> Unit,
  on_err~ : (Error) -> Unit,
) -> Unit {
  moonbit_test_driver_internal_js_catch(
    on_err=err_str => on_err(MoonBitTestDriverInternalJsError(err_str)),
    () => {
      try f() catch {
        err => on_err(err)
      } noraise {
        _ => on_ok()
      }
    },
  )
}

///|
#cfg(target="js")
extern "js" fn moonbit_test_driver_internal_js_catch(
  f : () -> Unit,
  on_err~ : (String) -> Unit,
) =
  #| (f, on_err) => {
  #|   try {
  #|     f()
  #|   } catch (err) {
  #|     const msg = err.stack.toString()
  #|     on_err(msg)
  #|   }
  #| }

///|
/// Maybe unused if there's no tests
#cfg(not(target="js"))
#warnings("-unused_value")
fn moonbit_test_driver_internal_catch_error(
  f : () -> Unit raise,
  on_ok~ : () -> Unit,
  on_err~ : (Error) -> Unit,
) -> Unit {
  try f() catch {
    err => on_err(err)
  } noraise {
    _ => on_ok()
  }
}

// =========================================================
// -=================== END JS SPECIFIC ====================
// =========================================================
// This file is introduced when tests without args is needed

///|
let moonbit_test_driver_internal_no_args_tests : MoonBitTestDriverInternalTestMap[
  () -> Unit raise,
] = {
  "skbio_tree_test.mbt": {
    0: (__test_736b62696f5f747265655f746573742e6d6274_0, ["tree_node_creation"]),
    1: (__test_736b62696f5f747265655f746573742e6d6274_1, ["tree_node_with_children"]),
    2: (__test_736b62696f5f747265655f746573742e6d6274_2, ["newick_parse_simple"]),
    3: (__test_736b62696f5f747265655f746573742e6d6274_3, ["newick_parse_with_lengths"]),
    4: (__test_736b62696f5f747265655f746573742e6d6274_4, ["newick_parse_nested"]),
    5: (__test_736b62696f5f747265655f746573742e6d6274_5, ["newick_to_newick"]),
    6: (__test_736b62696f5f747265655f746573742e6d6274_6, ["preorder_traversal"]),
    7: (__test_736b62696f5f747265655f746573742e6d6274_7, ["postorder_traversal"]),
    8: (__test_736b62696f5f747265655f746573742e6d6274_8, ["levelorder_traversal"]),
    9: (__test_736b62696f5f747265655f746573742e6d6274_9, ["tips"]),
    10: (__test_736b62696f5f747265655f746573742e6d6274_10, ["count_tips"]),
    11: (__test_736b62696f5f747265655f746573742e6d6274_11, ["count_nodes"]),
    12: (__test_736b62696f5f747265655f746573742e6d6274_12, ["find_node"]),
    13: (__test_736b62696f5f747265655f746573742e6d6274_13, ["find_not_found"]),
    14: (__test_736b62696f5f747265655f746573742e6d6274_14, ["total_length"]),
    15: (__test_736b62696f5f747265655f746573742e6d6274_15, ["lca"]),
    16: (__test_736b62696f5f747265655f746573742e6d6274_16, ["distance"]),
  },
  "skbio_tree_perf_test.mbt": {
    0: (__test_736b62696f5f747265655f706572665f746573742e6d6274_0, ["perf_lca_balanced_1024"]),
    1: (__test_736b62696f5f747265655f706572665f746573742e6d6274_1, ["perf_lca_balanced_4096_deep"]),
    2: (__test_736b62696f5f747265655f706572665f746573742e6d6274_2, ["perf_distance_balanced_1024"]),
    3: (__test_736b62696f5f747265655f706572665f746573742e6d6274_3, ["perf_distance_balanced_4096"]),
    4: (__test_736b62696f5f747265655f706572665f746573742e6d6274_4, ["perf_levelorder_4096"]),
    5: (__test_736b62696f5f747265655f706572665f746573742e6d6274_5, ["perf_find_4096"]),
    6: (__test_736b62696f5f747265655f706572665f746573742e6d6274_6, ["perf_parse_4096"]),
    7: (__test_736b62696f5f747265655f706572665f746573742e6d6274_7, ["perf_to_newick_4096"]),
    8: (__test_736b62696f5f747265655f706572665f746573742e6d6274_8, ["perf_lca_caterpillar_500"]),
    9: (__test_736b62696f5f747265655f706572665f746573742e6d6274_9, ["perf_distance_caterpillar_500"]),
  },
  "skbio_tree_bench.mbt": {
  },
  "skbio_tree.mbt": {
  },
  "skbio_diversity_test.mbt": {
    0: (__test_736b62696f5f6469766572736974795f746573742e6d6274_0, ["observed_count_basic"]),
    1: (__test_736b62696f5f6469766572736974795f746573742e6d6274_1, ["observed_count_all_zero"]),
    2: (__test_736b62696f5f6469766572736974795f746573742e6d6274_2, ["shannon_basic"]),
    3: (__test_736b62696f5f6469766572736974795f746573742e6d6274_3, ["shannon_single_species"]),
    4: (__test_736b62696f5f6469766572736974795f746573742e6d6274_4, ["shannon_empty"]),
    5: (__test_736b62696f5f6469766572736974795f746573742e6d6274_5, ["shannon_base2"]),
    6: (__test_736b62696f5f6469766572736974795f746573742e6d6274_6, ["simpson_basic"]),
    7: (__test_736b62696f5f6469766572736974795f746573742e6d6274_7, ["simpson_single_species"]),
    8: (__test_736b62696f5f6469766572736974795f746573742e6d6274_8, ["inv_simpson_basic"]),
    9: (__test_736b62696f5f6469766572736974795f746573742e6d6274_9, ["pielou_e_basic"]),
    10: (__test_736b62696f5f6469766572736974795f746573742e6d6274_10, ["pielou_e_single"]),
    11: (__test_736b62696f5f6469766572736974795f746573742e6d6274_11, ["simpson_e_basic"]),
    12: (__test_736b62696f5f6469766572736974795f746573742e6d6274_12, ["heip_e_basic"]),
    13: (__test_736b62696f5f6469766572736974795f746573742e6d6274_13, ["chao1_simple"]),
    14: (__test_736b62696f5f6469766572736974795f746573742e6d6274_14, ["chao1_no_singletons"]),
    15: (__test_736b62696f5f6469766572736974795f746573742e6d6274_15, ["chao1_no_doubletons"]),
    16: (__test_736b62696f5f6469766572736974795f746573742e6d6274_16, ["ace_basic"]),
    17: (__test_736b62696f5f6469766572736974795f746573742e6d6274_17, ["berger_parker_basic"]),
    18: (__test_736b62696f5f6469766572736974795f746573742e6d6274_18, ["margalef_basic"]),
    19: (__test_736b62696f5f6469766572736974795f746573742e6d6274_19, ["menhinick_basic"]),
    20: (__test_736b62696f5f6469766572736974795f746573742e6d6274_20, ["fisher_alpha_basic"]),
    21: (__test_736b62696f5f6469766572736974795f746573742e6d6274_21, ["braycurtis_identical"]),
    22: (__test_736b62696f5f6469766572736974795f746573742e6d6274_22, ["braycurtis_maximal"]),
    23: (__test_736b62696f5f6469766572736974795f746573742e6d6274_23, ["braycurtis_basic"]),
    24: (__test_736b62696f5f6469766572736974795f746573742e6d6274_24, ["jaccard_identical"]),
    25: (__test_736b62696f5f6469766572736974795f746573742e6d6274_25, ["jaccard_maximal"]),
    26: (__test_736b62696f5f6469766572736974795f746573742e6d6274_26, ["euclidean_basic"]),
    27: (__test_736b62696f5f6469766572736974795f746573742e6d6274_27, ["manhattan_basic"]),
    28: (__test_736b62696f5f6469766572736974795f746573742e6d6274_28, ["chebyshev_basic"]),
    29: (__test_736b62696f5f6469766572736974795f746573742e6d6274_29, ["canberra_basic"]),
    30: (__test_736b62696f5f6469766572736974795f746573742e6d6274_30, ["cosine_identical"]),
    31: (__test_736b62696f5f6469766572736974795f746573742e6d6274_31, ["cosine_orthogonal"]),
    32: (__test_736b62696f5f6469766572736974795f746573742e6d6274_32, ["faith_pd_simple"]),
    33: (__test_736b62696f5f6469766572736974795f746573742e6d6274_33, ["faith_pd_subset"]),
    34: (__test_736b62696f5f6469766572736974795f746573742e6d6274_34, ["faith_pd_single_taxon"]),
    35: (__test_736b62696f5f6469766572736974795f746573742e6d6274_35, ["gini_simpson_basic"]),
    36: (__test_736b62696f5f6469766572736974795f746573742e6d6274_36, ["dominance_basic"]),
    37: (__test_736b62696f5f6469766572736974795f746573742e6d6274_37, ["compare_with_skbio_alpha"]),
    38: (__test_736b62696f5f6469766572736974795f746573742e6d6274_38, ["compare_with_skbio_beta"]),
    39: (__test_736b62696f5f6469766572736974795f746573742e6d6274_39, ["compare_with_skbio_faith_pd"]),
  },
  "skbio_diversity.mbt": {
  },
  "skbio_alignment_test.mbt": {
    0: (__test_736b62696f5f616c69676e6d656e745f746573742e6d6274_0, ["dna_creation"]),
    1: (__test_736b62696f5f616c69676e6d656e745f746573742e6d6274_1, ["rna_creation"]),
    2: (__test_736b62696f5f616c69676e6d656e745f746573742e6d6274_2, ["protein_creation"]),
    3: (__test_736b62696f5f616c69676e6d656e745f746573742e6d6274_3, ["tabular_msa_creation"]),
    4: (__test_736b62696f5f616c69676e6d656e745f746573742e6d6274_4, ["tabular_msa_column"]),
    5: (__test_736b62696f5f616c69676e6d656e745f746573742e6d6274_5, ["tabular_msa_consensus"]),
    6: (__test_736b62696f5f616c69676e6d656e745f746573742e6d6274_6, ["tabular_msa_conservation"]),
    7: (__test_736b62696f5f616c69676e6d656e745f746573742e6d6274_7, ["global_pairwise_align_nucleotide"]),
    8: (__test_736b62696f5f616c69676e6d656e745f746573742e6d6274_8, ["local_pairwise_align_nucleotide"]),
    9: (__test_736b62696f5f616c69676e6d656e745f746573742e6d6274_9, ["global_pairwise_align_protein"]),
    10: (__test_736b62696f5f616c69676e6d656e745f746573742e6d6274_10, ["local_pairwise_align_protein"]),
    11: (__test_736b62696f5f616c69676e6d656e745f746573742e6d6274_11, ["align_score"]),
    12: (__test_736b62696f5f616c69676e6d656e745f746573742e6d6274_12, ["global_alignment_score"]),
  },
  "skbio_alignment.mbt": {
  },
  "sequtils_test.mbt": {
    0: (__test_7365717574696c735f746573742e6d6274_0, ["crc32_basic"]),
    1: (__test_7365717574696c735f746573742e6d6274_1, ["crc32_simple"]),
    2: (__test_7365717574696c735f746573742e6d6274_2, ["gcg_basic"]),
    3: (__test_7365717574696c735f746573742e6d6274_3, ["seguid_basic"]),
    4: (__test_7365717574696c735f746573742e6d6274_4, ["sha1_basic"]),
    5: (__test_7365717574696c735f746573742e6d6274_5, ["sha1_h0_after_round"]),
    6: (__test_7365717574696c735f746573742e6d6274_6, ["sha1_raw_bytes"]),
    7: (__test_7365717574696c735f746573742e6d6274_7, ["sha1_debug"]),
    8: (__test_7365717574696c735f746573742e6d6274_8, ["base64_test"]),
    9: (__test_7365717574696c735f746573742e6d6274_9, ["string_to_bytes_test"]),
    10: (__test_7365717574696c735f746573742e6d6274_10, ["gc_simple"]),
    11: (__test_7365717574696c735f746573742e6d6274_11, ["gc_full"]),
    12: (__test_7365717574696c735f746573742e6d6274_12, ["protein_molecular_weight"]),
    13: (__test_7365717574696c735f746573742e6d6274_13, ["protein_aromaticity"]),
    14: (__test_7365717574696c735f746573742e6d6274_14, ["protein_gravy"]),
    15: (__test_7365717574696c735f746573742e6d6274_15, ["tm_wallace_simple"]),
    16: (__test_7365717574696c735f746573742e6d6274_16, ["tm_wallace_gc"]),
    17: (__test_7365717574696c735f746573742e6d6274_17, ["tm_gc_simple"]),
    18: (__test_7365717574696c735f746573742e6d6274_18, ["tm_gc_long"]),
  },
  "sequtils.mbt": {
  },
  "seqio.mbt": {
  },
  "seqfeature_test.mbt": {
    0: (__test_736571666561747572655f746573742e6d6274_0, ["feature_location_basic"]),
    1: (__test_736571666561747572655f746573742e6d6274_1, ["feature_location_strand"]),
    2: (__test_736571666561747572655f746573742e6d6274_2, ["seqfeature_basic"]),
    3: (__test_736571666561747572655f746573742e6d6274_3, ["seqfeature_extract"]),
    4: (__test_736571666561747572655f746573742e6d6274_4, ["seqfeature_extract_minus_strand"]),
    5: (__test_736571666561747572655f746573742e6d6274_5, ["seqfeature_qualifiers"]),
    6: (__test_736571666561747572655f746573742e6d6274_6, ["compound_location"]),
    7: (__test_736571666561747572655f746573742e6d6274_7, ["position_types"]),
  },
  "seqfeature.mbt": {
  },
  "seq_record.mbt": {
  },
  "seq.mbt": {
  },
  "phylo_test.mbt": {
    0: (__test_7068796c6f5f746573742e6d6274_0, ["clade_is_terminal"]),
    1: (__test_7068796c6f5f746573742e6d6274_1, ["clade_count_terminals"]),
    2: (__test_7068796c6f5f746573742e6d6274_2, ["clade_get_terminals"]),
    3: (__test_7068796c6f5f746573742e6d6274_3, ["clade_get_nonterminals"]),
    4: (__test_7068796c6f5f746573742e6d6274_4, ["clade_find_clades"]),
    5: (__test_7068796c6f5f746573742e6d6274_5, ["clade_find_any"]),
    6: (__test_7068796c6f5f746573742e6d6274_6, ["clade_get_path"]),
    7: (__test_7068796c6f5f746573742e6d6274_7, ["clade_depths"]),
    8: (__test_7068796c6f5f746573742e6d6274_8, ["clade_is_bifurcating"]),
    9: (__test_7068796c6f5f746573742e6d6274_9, ["tree_basic"]),
    10: (__test_7068796c6f5f746573742e6d6274_10, ["tree_get_terminals"]),
    11: (__test_7068796c6f5f746573742e6d6274_11, ["tree_find_clades"]),
    12: (__test_7068796c6f5f746573742e6d6274_12, ["tree_find_any"]),
    13: (__test_7068796c6f5f746573742e6d6274_13, ["tree_distance"]),
    14: (__test_7068796c6f5f746573742e6d6274_14, ["tree_common_ancestor"]),
    15: (__test_7068796c6f5f746573742e6d6274_15, ["newick_parse_simple"]),
    16: (__test_7068796c6f5f746573742e6d6274_16, ["newick_parse_nested"]),
    17: (__test_7068796c6f5f746573742e6d6274_17, ["newick_parse_with_branch_lengths"]),
    18: (__test_7068796c6f5f746573742e6d6274_18, ["newick_parse_no_semicolon"]),
    19: (__test_7068796c6f5f746573742e6d6274_19, ["newick_roundtrip"]),
    20: (__test_7068796c6f5f746573742e6d6274_20, ["tree_draw_ascii"]),
    21: (__test_7068796c6f5f746573742e6d6274_21, ["tree_repr"]),
  },
  "phylo.mbt": {
  },
  "phylip_io.mbt": {
  },
  "pdb_test.mbt": {
    0: (__test_7064625f746573742e6d6274_0, ["vector3_distance"]),
    1: (__test_7064625f746573742e6d6274_1, ["atom_creation"]),
    2: (__test_7064625f746573742e6d6274_2, ["atom_distance"]),
    3: (__test_7064625f746573742e6d6274_3, ["residue_creation"]),
    4: (__test_7064625f746573742e6d6274_4, ["residue_find_atom"]),
    5: (__test_7064625f746573742e6d6274_5, ["residue_get_ca"]),
    6: (__test_7064625f746573742e6d6274_6, ["residue_center_of_mass"]),
    7: (__test_7064625f746573742e6d6274_7, ["chain_creation"]),
    8: (__test_7064625f746573742e6d6274_8, ["chain_find_residue"]),
    9: (__test_7064625f746573742e6d6274_9, ["model_creation"]),
    10: (__test_7064625f746573742e6d6274_10, ["model_find_chain"]),
    11: (__test_7064625f746573742e6d6274_11, ["structure_creation"]),
    12: (__test_7064625f746573742e6d6274_12, ["pdb_parse_simple"]),
    13: (__test_7064625f746573742e6d6274_13, ["pdb_parse_multiple_residues"]),
    14: (__test_7064625f746573742e6d6274_14, ["pdb_parse_multiple_chains"]),
    15: (__test_7064625f746573742e6d6274_15, ["pdb_parse_with_model"]),
    16: (__test_7064625f746573742e6d6274_16, ["pdb_write"]),
    17: (__test_7064625f746573742e6d6274_17, ["pdb_roundtrip"]),
    18: (__test_7064625f746573742e6d6274_18, ["structure_distance"]),
  },
  "pdb_io.mbt": {
  },
  "pdb.mbt": {
  },
  "newick_io.mbt": {
  },
  "genbank_io.mbt": {
  },
  "fastq_io.mbt": {
  },
  "fasta_io.mbt": {
  },
  "complement.mbt": {
  },
  "codon_table.mbt": {
  },
  "clustal_io.mbt": {
  },
  "bio_seq_test.mbt": {
    0: (__test_62696f5f7365715f746573742e6d6274_0, ["construction and to_string"]),
    1: (__test_62696f5f7365715f746573742e6d6274_1, ["complement and complement_rna"]),
    2: (__test_62696f5f7365715f746573742e6d6274_2, ["reverse_complement"]),
    3: (__test_62696f5f7365715f746573742e6d6274_3, ["transcribe and back_transcribe"]),
    4: (__test_62696f5f7365715f746573742e6d6274_4, ["translate default"]),
    5: (__test_62696f5f7365715f746573742e6d6274_5, ["translate to_stop"]),
    6: (__test_62696f5f7365715f746573742e6d6274_6, ["translate cds valid"]),
    7: (__test_62696f5f7365715f746573742e6d6274_7, ["translate cds invalid start"]),
    8: (__test_62696f5f7365715f746573742e6d6274_8, ["translate cds internal stop"]),
    9: (__test_62696f5f7365715f746573742e6d6274_9, ["translate partial codon"]),
    10: (__test_62696f5f7365715f746573742e6d6274_10, ["translate ambiguous codon"]),
    11: (__test_62696f5f7365715f746573742e6d6274_11, ["translate invalid codon raises"]),
    12: (__test_62696f5f7365715f746573742e6d6274_12, ["translate stop_symbol"]),
    13: (__test_62696f5f7365715f746573742e6d6274_13, ["upper and lower"]),
    14: (__test_62696f5f7365715f746573742e6d6274_14, ["count"]),
    15: (__test_62696f5f7365715f746573742e6d6274_15, ["find and rfind"]),
    16: (__test_62696f5f7365715f746573742e6d6274_16, ["startswith and endswith"]),
    17: (__test_62696f5f7365715f746573742e6d6274_17, ["contains"]),
    18: (__test_62696f5f7365715f746573742e6d6274_18, ["add and mul"]),
    19: (__test_62696f5f7365715f746573742e6d6274_19, ["slice"]),
    20: (__test_62696f5f7365715f746573742e6d6274_20, ["get_char"]),
    21: (__test_62696f5f7365715f746573742e6d6274_21, ["strip lstrip rstrip"]),
    22: (__test_62696f5f7365715f746573742e6d6274_22, ["replace"]),
    23: (__test_62696f5f7365715f746573742e6d6274_23, ["split"]),
    24: (__test_62696f5f7365715f746573742e6d6274_24, ["join"]),
    25: (__test_62696f5f7365715f746573742e6d6274_25, ["equality"]),
    26: (__test_62696f5f7365715f746573742e6d6274_26, ["reverse"]),
    27: (__test_62696f5f7365715f746573742e6d6274_27, ["complement_rna round trip"]),
  },
  "bio_seq.mbt": {
  },
  "alignio.mbt": {
  },
  "align.mbt": {
  },
  "README.mbt.md": {
  },
}


///|
impl MoonBit_Test_Driver for MoonBit_Test_Driver_Internal_No_Args with fn run_test(
  _,
  filename,
  index,
  handle_start,
  handle_result,
  error_to_string,
) {
  if moonbit_test_driver_internal_no_args_tests.0.get(filename)
    is Some(index_map) &&
    index_map.get(index) is Some((f, attrs)) {
    let name = if attrs is [name, ..] { name } else { "" }
    if MOONBIT_TEST_DRIVER_INTERNAL_IS_NATIVE {
      for attr in attrs {
        if attr is [.. "panic", ..] {
          raise MoonBitTestDriverInternalSkipTest(name~)
        }
      }
    }
    handle_start()
    moonbit_test_driver_internal_catch_error(
      f,
      on_ok=() => handle_result(name, "", false),
      on_err=err => handle_result(name, error_to_string(err), false),
    )
    true
  } else {
    false
  }
}
///|
#cfg(any(target="native", target="llvm"))
#test_entry
fn main {
  let async_tests = []
  for arg in moonbit_test_driver_internal_native_parse_args() {
    moonbit_test_driver_internal_do_execute(async_tests, arg.0, arg.1)
  }
  MoonBit_Async_Test_Driver_Impl::run_async_tests(async_tests)
}

///|
#cfg(target="js")
#test_entry
fn main {
  let test_params = moonbit_test_driver_js_parse_args()
  moonbit_test_driver_internal_execute(test_params)
  moonbit_test_driver_finish()
}

///|
#cfg(any(target="wasm", target="wasm-gc"))
#test_entry
fn main {
  ()
}

///|
#cfg(target="js")
extern "js" fn moonbit_test_driver_js_parse_args() -> Array[
  MoonBitTestDriverInternalTestEntry,
] =
  #|() => {
  #|  try {
  #|    let testParams = []
  #|    // check if there is command line argument
  #|    if (process.argv.length > 2) {
  #|      const testArgs = JSON.parse(process.argv[2])
  #|      testParams = testArgs.file_and_index.flatMap(([file, ranges]) => {
  #|        let arr = []
  #|        for (let range of ranges) {
  #|          for (let i = range.start; i < range.end; i++) {
  #|            arr.push([file, i])
  #|          }
  #|        }
  #|        return arr
  #|      })
  #|      return testParams
  #|    }
  #|  } catch (error) {
  #|    console.error('failed to parse args: ', error.message)
  #|    process.exit(1)
  #|  }
  #|}

///|
/// Executes selected JavaScript tests from parsed test entries.
#cfg(target="js")
pub fn moonbit_test_driver_internal_execute(
  tests : Array[MoonBitTestDriverInternalTestEntry],
) -> Unit {
  let async_tests = []
  for i in 0..<tests.length() {
    let entry = tests[i]
    moonbit_test_driver_internal_do_execute(
      async_tests,
      entry.filename(),
      entry.index(),
    )
  }
  MoonBit_Async_Test_Driver_Impl::run_async_tests(async_tests)
}

///|
#cfg(target="wasm")
let moonbit_test_driver_internal_wasm_async_tests : Array[async () -> Unit] = []

///|
/// Executes one selected WebAssembly test and stores async work for finalization.
#cfg(target="wasm")
pub fn moonbit_test_driver_internal_execute(
  filename : MoonbitTestDriverInternalExternString,
  index : Int,
) -> Unit {
  let filename = moonbit_test_driver_internal_get_file_name(filename)
  moonbit_test_driver_internal_do_execute(
    moonbit_test_driver_internal_wasm_async_tests, filename, index,
  )
}

///|
/// Executes one selected WebAssembly GC test and drains async work immediately.
#cfg(target="wasm-gc")
pub fn moonbit_test_driver_internal_execute(
  filename : MoonbitTestDriverInternalExternString,
  index : Int,
) -> Unit {
  let filename = moonbit_test_driver_internal_get_file_name(filename)
  let async_tests = []
  moonbit_test_driver_internal_do_execute(async_tests, filename, index)
  MoonBit_Async_Test_Driver_Impl::run_async_tests(async_tests)
}

///|
/// Finishes pending WebAssembly async tests and coverage collection.
#cfg(target="wasm")
pub fn moonbit_test_driver_finish() -> Unit {
  MoonBit_Async_Test_Driver_Impl::run_async_tests(
    moonbit_test_driver_internal_wasm_async_tests,
  )
  
}

///|
/// Finishes coverage collection for backends without pending async work.
#cfg(any(target="js", target="wasm-gc"))
pub fn moonbit_test_driver_finish() -> Unit {
  
}

///|
fn moonbit_test_driver_internal_do_execute(
  async_tests : Array[async () -> Unit],
  filename : String,
  index : Int,
) -> Unit {
  // These two are internalized to avoid unused warnings

  fn handle_start() {
    let file_name = filename.escape()
    @moonbitlang/core/prelude.println("----- BEGIN MOON TEST RESULT -----")
    @moonbitlang/core/prelude.println(
      "{\"type\":\"start\",\"file\":\{file_name},\"index\":\{index}}",
    )
    @moonbitlang/core/prelude.println("----- END MOON TEST RESULT -----")
  }

  fn handle_result(_testname : String, message : String, skipped : Bool) {
    // In WASM, this function is called by the external test driver, but here
    // in native mode there's nothing that calls it, so we call it manually.
    if !skipped && MOONBIT_TEST_DRIVER_INTERNAL_IS_NATIVE {
      
    }
    let file_name = filename.escape()
    let message = message.escape()
    @moonbitlang/core/prelude.println("----- BEGIN MOON TEST RESULT -----")
    @moonbitlang/core/prelude.println(
      "{\"type\":\"result\",\"file\":\{file_name},\"index\":\{index},\"message\":\{message}}",
    )
    @moonbitlang/core/prelude.println("----- END MOON TEST RESULT -----")
  }

  fn error_to_string(err : Error) -> String {
    match err {
      @moonbitlang/core/prelude.Failure::Failure(e)
      | @moonbitlang/core/prelude.InspectError::InspectError(e)
      | @moonbitlang/core/prelude.SnapshotError::SnapshotError(e)
      | MoonBitTestDriverInternalJsError(e) => e
      e => moonbit_test_driver_internal_error_to_string(e)
    }
  }

  try {
    if !(MoonBit_Test_Driver_Internal_No_Args::run_test(
        async_tests, filename, index, handle_start, handle_result, error_to_string,
      ) ||
      MoonBit_Test_Driver_Internal_With_Args::run_test(
        async_tests, filename, index, handle_start, handle_result, error_to_string,
      ) ||
      MoonBit_Test_Driver_Internal_Async_No_Args::run_test(
        async_tests, filename, index, handle_start, handle_result, error_to_string,
      ) ||
      MoonBit_Test_Driver_Internal_Async_With_Args::run_test(
        async_tests, filename, index, handle_start, handle_result, error_to_string,
      ) ||
      MoonBit_Test_Driver_Internal_With_Bench_Args::run_test(
        async_tests, filename, index, handle_start, handle_result, error_to_string,
      )) {
      raise MoonBitTestDriverInternalSkipTest(name="")
    }
  } catch {
    err => {
      guard err is MoonBitTestDriverInternalSkipTest(name~)
      handle_result(name, "skipped test", true)
    }
  }
}
