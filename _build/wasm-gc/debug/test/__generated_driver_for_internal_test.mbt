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
  "skbio_tree_bench.mbt": {
    0: (__test_736b62696f5f747265655f62656e63682e6d6274_0, ["bench_setup_check"]),
    1: (__test_736b62696f5f747265655f62656e63682e6d6274_1, ["bench_newick_parse_small"]),
    2: (__test_736b62696f5f747265655f62656e63682e6d6274_2, ["bench_newick_parse_large"]),
    3: (__test_736b62696f5f747265655f62656e63682e6d6274_3, ["bench_to_newick_small"]),
    4: (__test_736b62696f5f747265655f62656e63682e6d6274_4, ["bench_to_newick_large"]),
    5: (__test_736b62696f5f747265655f62656e63682e6d6274_5, ["bench_preorder"]),
    6: (__test_736b62696f5f747265655f62656e63682e6d6274_6, ["bench_postorder"]),
    7: (__test_736b62696f5f747265655f62656e63682e6d6274_7, ["bench_levelorder"]),
    8: (__test_736b62696f5f747265655f62656e63682e6d6274_8, ["bench_find"]),
    9: (__test_736b62696f5f747265655f62656e63682e6d6274_9, ["bench_find_last"]),
    10: (__test_736b62696f5f747265655f62656e63682e6d6274_10, ["bench_lca"]),
    11: (__test_736b62696f5f747265655f62656e63682e6d6274_11, ["bench_lca_deep"]),
    12: (__test_736b62696f5f747265655f62656e63682e6d6274_12, ["bench_distance"]),
    13: (__test_736b62696f5f747265655f62656e63682e6d6274_13, ["bench_total_length"]),
    14: (__test_736b62696f5f747265655f62656e63682e6d6274_14, ["bench_count_tips"]),
    15: (__test_736b62696f5f747265655f62656e63682e6d6274_15, ["bench_caterpillar_find"]),
    16: (__test_736b62696f5f747265655f62656e63682e6d6274_16, ["bench_caterpillar_lca"]),
    17: (__test_736b62696f5f747265655f62656e63682e6d6274_17, ["bench_caterpillar_distance"]),
  },
  "skbio_tree.mbt": {
  },
  "skbio_stats.mbt": {
  },
  "skbio_faidx.mbt": {
  },
  "skbio_diversity.mbt": {
  },
  "skbio_common.mbt": {
  },
  "skbio_alignment.mbt": {
  },
  "sequtils.mbt": {
  },
  "seqio.mbt": {
  },
  "seqfeature.mbt": {
  },
  "seq_record.mbt": {
  },
  "seq.mbt": {
  },
  "phylo.mbt": {
  },
  "phylip_io.mbt": {
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
  "bio_seq.mbt": {
  },
  "alignio.mbt": {
  },
  "align.mbt": {
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
