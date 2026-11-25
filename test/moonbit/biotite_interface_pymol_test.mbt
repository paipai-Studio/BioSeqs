///|
/// Tests for biotite_interface_pymol.mbt.
test "biotite_pymol_script_commands" {
  let script = @src.biotite_pymol_script()
  let loaded = @src.biotite_pymol_load(script, "a.pdb")
  let selected = @src.biotite_pymol_select(loaded, "chainA", "chain A")
  let colored = @src.biotite_pymol_color(selected, "red", "chainA")
  let text = @src.biotite_pymol_to_script(colored)
  assert_true(text.contains("cmd.load(\"a.pdb\")"))
  assert_true(text.contains("cmd.select(\"chainA\", \"chain A\")"))
  assert_true(text.contains("cmd.color(\"red\", \"chainA\")"))
}

///|
test "biotite_pymol_show_hide_save" {
  let script = @src.biotite_pymol_show(
    @src.biotite_pymol_hide(@src.biotite_pymol_script(), "lines", "all"),
    "cartoon",
    "all",
  )
  let saved = @src.biotite_pymol_save(script, "out.pse")
  let text = @src.biotite_pymol_to_script(saved)
  assert_true(text.contains("cmd.hide(\"lines\", \"all\")"))
  assert_true(text.contains("cmd.show(\"cartoon\", \"all\")"))
  assert_true(text.contains("cmd.save(\"out.pse\")"))
}

///|
test "biotite_pymol_empty_script" {
  assert_eq(@src.biotite_pymol_to_script(@src.biotite_pymol_script()), "")
}
