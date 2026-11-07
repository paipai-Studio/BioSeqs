/// omm_xml_ff.mbt — OpenMM XML force field parser + residue template matching
/// Item 3: OpenMM native XML force field format

///|
pub struct OmmXmlForceField {
  name : String
  // Atom types with parameters
  atom_types : Array[OmmXmlAtomType]
  // Residue templates
  residues : Array[OmmXmlResidue]
  // Bond force: (type1, type2, k, r0)
  bond_forces : Array[(String, String, Double, Double)]
  // Angle force: (t1, t2, t3, k, theta0)
  angle_forces : Array[(String, String, String, Double, Double)]
  // Periodic torsion: (t1, t2, t3, t4, k, phi0, n)
  torsion_forces : Array[
    (String, String, String, String, Double, Double, Int32),
  ]
  // Nonbonded: (type, charge, sigma, epsilon)
  nonbonded : Array[(String, Double, Double, Double)]
}

///|
pub struct OmmXmlAtomType {
  name : String
  element : String
  mass : Double
}

///|
pub struct OmmXmlResidue {
  name : String
  atoms : Array[(String, String)] // (atom_name, atom_type)
  bonds : Array[(String, String)] // (atom_name1, atom_name2)
}

///|
/// Parse OpenMM XML force field (simplified — extracts key sections).
pub fn omm_parse_xml_ff(content : String) -> OmmXmlForceField {
  let atom_types : Array[OmmXmlAtomType] = Array(capacity=0)
  let residues : Array[OmmXmlResidue] = Array(capacity=0)
  let bond_forces : Array[(String, String, Double, Double)] = Array(capacity=0)
  let angle_forces : Array[(String, String, String, Double, Double)] = Array(
    capacity=0,
  )
  let torsion_forces : Array[
    (String, String, String, String, Double, Double, Int32),
  ] = Array(capacity=0)
  let nonbonded : Array[(String, Double, Double, Double)] = Array(capacity=0)
  let lines = content.split("\n").to_array()
  let mut name = "unknown"
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned().trim().to_owned()
    // Extract force field name
    if line.contains("<ForceField") && line.contains("name=") {
      let start = line.find("name=\"")
      match start {
        Some(s) => {
          let rest = line[s + 6:].to_owned()
          let end = rest.find("\"")
          match end {
            Some(e) => name = rest[:e].to_owned()
            None => ()
          }
        }
        None => ()
      }
    }
    // Atom types
    if line.contains("<Type") && line.contains("name=") {
      let an = omm_xml_attr(line, "name")
      let el = omm_xml_attr(line, "element")
      let m = parse_double(omm_xml_attr(line, "mass")).unwrap_or(12.0)
      atom_types.push(OmmXmlAtomType::{ name: an, element: el, mass: m, })
    }
    // Residue
    if line.contains("<Residue") && line.contains("name=") {
      let rn = omm_xml_attr(line, "name")
      let ratoms : Array[(String, String)] = Array(capacity=0)
      let rbonds : Array[(String, String)] = Array(capacity=0)
      // Read residue content until </Residue>
      let mut j = i + 1
      while j < lines.length() {
        let rl = lines[j].to_owned().trim().to_owned()
        if rl.contains("</Residue>") {
          break
        }
        if rl.contains("<Atom") {
          let an = omm_xml_attr(rl, "name")
          let at = omm_xml_attr(rl, "type")
          ratoms.push((an, at))
        }
        if rl.contains("<Bond") {
          let b1 = omm_xml_attr(rl, "from")
          let b2 = omm_xml_attr(rl, "to")
          rbonds.push((b1, b2))
        }
        j = j + 1
      }
      residues.push(OmmXmlResidue::{ name: rn, atoms: ratoms, bonds: rbonds, })
    }
    // Bond force
    if line.contains("<Bond") && line.contains("k=") {
      let t1 = omm_xml_attr(line, "type1")
      let t2 = omm_xml_attr(line, "type2")
      let k = parse_double(omm_xml_attr(line, "k")).unwrap_or(0.0)
      let r0 = parse_double(omm_xml_attr(line, "r0")).unwrap_or(0.0)
      if t1.length() > 0 {
        bond_forces.push((t1, t2, k, r0))
      }
    }
    // Nonbonded
    if line.contains("<Atom") && line.contains("sigma=") {
      let t = omm_xml_attr(line, "type")
      let q = parse_double(omm_xml_attr(line, "charge")).unwrap_or(0.0)
      let s = parse_double(omm_xml_attr(line, "sigma")).unwrap_or(0.0)
      let e = parse_double(omm_xml_attr(line, "epsilon")).unwrap_or(0.0)
      nonbonded.push((t, q, s, e))
    }
  }
  OmmXmlForceField::{
    name,
    atom_types,
    residues,
    bond_forces,
    angle_forces,
    torsion_forces,
    nonbonded,
  }
}

///|
/// Extract XML attribute value by name.
fn omm_xml_attr(line : String, attr : String) -> String {
  let key = attr + "=\""
  let start = line.find(key)
  match start {
    Some(s) => {
      let rest = line[s + key.length():].to_owned()
      let end = rest.find("\"")
      match end {
        Some(e) => rest[:e].to_owned()
        None => ""
      }
    }
    None => ""
  }
}

///|
/// Match a residue to a force field template and assign atom types.
/// Returns array of (atom_name, atom_type) assignments.
pub fn omm_match_residue_template(
  ff : OmmXmlForceField,
  res_name : String,
  atom_names : Array[String],
) -> Array[(String, String)] {
  // Find matching residue template
  let mut matched : OmmXmlResidue? = None
  for r in ff.residues {
    if r.name.to_lower() == res_name.to_lower() {
      matched = Some(r)
      break
    }
  }
  match matched {
    Some(r) => {
      // Match atoms by name
      let result : Array[(String, String)] = Array(capacity=0)
      for an in atom_names {
        let mut atype = "X"
        for (tn, tt) in r.atoms {
          if tn == an {
            atype = tt
            break
          }
        }
        result.push((an, atype))
      }
      result
    }
    None => {
      // No template found — return unknown types
      let result : Array[(String, String)] = Array(capacity=0)
      for an in atom_names {
        result.push((an, "X"))
      }
      result
    }
  }
}

///|
/// Look up nonbonded parameters for an atom type.
pub fn omm_lookup_nonbonded(
  ff : OmmXmlForceField,
  atom_type : String,
) -> (Double, Double, Double) {
  for (t, q, s, e) in ff.nonbonded {
    if t == atom_type {
      return (q, s, e)
    }
  }
  (0.0, 0.0, 0.0)
}

///|
pub fn omm_xml_ff_dummy() -> Unit {
  ignore(OmmXmlForceField::{
    name: "",
    atom_types: Array(capacity=0),
    residues: Array(capacity=0),
    bond_forces: Array(capacity=0),
    angle_forces: Array(capacity=0),
    torsion_forces: Array(capacity=0),
    nonbonded: Array(capacity=0),
  })
  ignore(OmmXmlAtomType::{ name: "", element: "", mass: 0.0, })
  ignore(OmmXmlResidue::{
    name: "",
    atoms: Array(capacity=0),
    bonds: Array(capacity=0),
  })
}
