///|
/// biotite_interface_pymol.mbt — PyMOL command script construction.
/// Only generates commands; it never starts PyMOL.

///|
pub(all) struct BiotitePymolScript {
  commands : Array[String]
}

///|
pub fn biotite_pymol_script() -> BiotitePymolScript {
  { commands: [], }
}

///|
fn biotite_pymol_append(
  script : BiotitePymolScript,
  command : String,
) -> BiotitePymolScript {
  let commands : Array[String] = Array(capacity=script.commands.length() + 1)
  for c in script.commands {
    commands.push(c)
  }
  commands.push(command)
  { commands, }
}

///|
pub fn biotite_pymol_load(
  script : BiotitePymolScript,
  file : String,
) -> BiotitePymolScript {
  biotite_pymol_append(script, "cmd.load(\"" + file + "\")")
}

///|
pub fn biotite_pymol_select(
  script : BiotitePymolScript,
  name : String,
  expression : String,
) -> BiotitePymolScript {
  biotite_pymol_append(
    script,
    "cmd.select(\"" + name + "\", \"" + expression + "\")",
  )
}

///|
pub fn biotite_pymol_color(
  script : BiotitePymolScript,
  color : String,
  selection : String,
) -> BiotitePymolScript {
  biotite_pymol_append(
    script,
    "cmd.color(\"" + color + "\", \"" + selection + "\")",
  )
}

///|
pub fn biotite_pymol_show(
  script : BiotitePymolScript,
  representation : String,
  selection : String,
) -> BiotitePymolScript {
  biotite_pymol_append(
    script,
    "cmd.show(\"" + representation + "\", \"" + selection + "\")",
  )
}

///|
pub fn biotite_pymol_hide(
  script : BiotitePymolScript,
  representation : String,
  selection : String,
) -> BiotitePymolScript {
  biotite_pymol_append(
    script,
    "cmd.hide(\"" + representation + "\", \"" + selection + "\")",
  )
}

///|
pub fn biotite_pymol_save(
  script : BiotitePymolScript,
  file : String,
) -> BiotitePymolScript {
  biotite_pymol_append(script, "cmd.save(\"" + file + "\")")
}

///|
/// Render the accumulated commands as a newline-separated script.
pub fn biotite_pymol_to_script(script : BiotitePymolScript) -> String {
  let sb = StringBuilder()
  for command in script.commands {
    let _ = sb.write(command + "\n")
  }
  sb.to_string()
}
