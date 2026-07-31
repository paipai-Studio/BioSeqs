///|
/// Tests for file module (SmartFile).

///|
test "compression_format_to_string" {
  assert_eq(@src.compression_plain().to_string(), "plain")
  assert_eq(@src.compression_gzip().to_string(), "gzip")
  assert_eq(@src.compression_bzip2().to_string(), "bzip2")
}

///|
test "detect_compression_plain" {
  let f = @src.detect_compression("data.txt")
  assert_eq(f.to_string(), "plain")
}

///|
test "detect_compression_gzip" {
  let f = @src.detect_compression("data.txt.gz")
  assert_eq(f.to_string(), "gzip")
}

///|
test "detect_compression_bzip2" {
  let f = @src.detect_compression("data.txt.bz2")
  assert_eq(f.to_string(), "bzip2")
}

///|
test "detect_compression_uppercase" {
  let f = @src.detect_compression("data.GZ")
  assert_eq(f.to_string(), "gzip")
}

///|
test "smart_file_creation" {
  let sf = @src.SmartFile::new("test.txt", mode="r")
  assert_eq(sf.get_path(), "test.txt")
  assert_eq(sf.get_mode(), "r")
  assert_true(sf.is_read_mode())
  assert_false(sf.is_write_mode())
}

///|
test "smart_file_creation_write" {
  let sf = @src.SmartFile::new("out.txt", mode="w")
  assert_eq(sf.get_mode(), "w")
  assert_true(sf.is_write_mode())
  assert_false(sf.is_read_mode())
}

///|
test "smart_file_with_format" {
  let sf = @src.SmartFile::with_format("data.gz", @src.compression_gzip(), mode="r")
  assert_eq(sf.get_format().to_string(), "gzip")
}

///|
test "smart_file_open_close" {
  let sf = @src.SmartFile::new("test.txt", mode="r")
  assert_false(sf.is_open())
  let opened = sf.open()
  assert_true(opened.is_open())
  let closed = opened.close()
  assert_false(closed.is_open())
}

///|
test "smart_file_set_lines" {
  let sf = @src.SmartFile::new("test.txt", mode="r")
  let lines = ["hello", "world"]
  let sf2 = sf.set_lines(lines)
  assert_eq(sf2.lines.length(), 2)
  assert_eq(sf2.lines[0], "hello")
  assert_eq(sf2.lines[1], "world")
}

///|
test "smart_file_read_line" {
  let sf = @src.SmartFile::new("test.txt", mode="r")
  let lines = ["line1", "line2", "line3"]
  let sf2 = sf.set_lines(lines)
  match sf2.read_line() {
    Some(l) => assert_eq(l, "line1")
    None => assert_true(false)
  }
}

///|
test "smart_file_read_line_empty" {
  let sf = @src.SmartFile::new("test.txt", mode="r")
  match sf.read_line() {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "smart_file_read_all" {
  let sf = @src.SmartFile::new("test.txt", mode="r")
  let lines = ["a", "b", "c"]
  let sf2 = sf.set_lines(lines)
  let all = sf2.read_all()
  assert_eq(all, "a\nb\nc")
}

///|
test "smart_file_read_lines" {
  let sf = @src.SmartFile::new("test.txt", mode="r")
  let lines = ["x", "y"]
  let sf2 = sf.set_lines(lines)
  let read = sf2.read_lines()
  assert_eq(read.length(), 2)
  assert_eq(read[0], "x")
  assert_eq(read[1], "y")
}

///|
test "smart_file_write_line" {
  let sf = @src.SmartFile::new("out.txt", mode="w")
  let sf2 = sf.write_line("hello")
  let sf3 = sf2.write_line("world")
  assert_eq(sf3.get_written_content(), "hello\nworld")
}

///|
test "smart_file_write" {
  let sf = @src.SmartFile::new("out.txt", mode="w")
  let sf2 = sf.write("single content")
  assert_eq(sf2.get_written_content(), "single content")
}

///|
test "smart_file_from_string" {
  let sf = @src.smart_file_from_string("test.txt", "a\nb\nc")
  assert_eq(sf.lines.length(), 3)
  assert_eq(sf.lines[0], "a")
  assert_eq(sf.lines[1], "b")
  assert_eq(sf.lines[2], "c")
}

///|
test "smart_file_from_lines" {
  let sf = @src.smart_file_from_lines("test.txt", ["hello", "world"])
  assert_eq(sf.lines.length(), 2)
}

///|
test "smart_file_writer" {
  let sf = @src.smart_file_writer("out.txt", format=@src.compression_gzip())
  assert_eq(sf.get_mode(), "w")
  assert_eq(sf.get_format().to_string(), "gzip")
}

///|
test "has_compression_extension_gz" {
  assert_true(@src.has_compression_extension("data.gz"))
}

///|
test "has_compression_extension_bz2" {
  assert_true(@src.has_compression_extension("data.bz2"))
}

///|
test "has_compression_extension_plain" {
  assert_false(@src.has_compression_extension("data.txt"))
}

///|
test "strip_compression_extension_gz" {
  let result = @src.strip_compression_extension("data.txt.gz")
  assert_eq(result, "data.txt")
}

///|
test "strip_compression_extension_bz2" {
  let result = @src.strip_compression_extension("data.txt.bz2")
  assert_eq(result, "data.txt")
}

///|
test "strip_compression_extension_plain" {
  let result = @src.strip_compression_extension("data.txt")
  assert_eq(result, "data.txt")
}

///|
test "smart_file_summary" {
  let sf = @src.SmartFile::new("test.txt", mode="r")
  let s = sf.summary()
  assert_true(s.length() > 0)
}

///|
test "file_sample_data" {
  let data = @src.file_sample_data()
  assert_eq(data.length(), 4)
}

///|
test "file_sample_lines" {
  let data = @src.file_sample_lines()
  assert_eq(data.length(), 3)
}
