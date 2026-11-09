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
    // Angle force: <Angle type1= type2= type3= k= theta0=/> (theta0 in deg)
    if line.contains("<Angle") && line.contains("k=") {
      let t1 = omm_xml_attr(line, "type1")
      let t2 = omm_xml_attr(line, "type2")
      let t3 = omm_xml_attr(line, "type3")
      let k = parse_double(omm_xml_attr(line, "k")).unwrap_or(0.0)
      let th = parse_double(omm_xml_attr(line, "theta0")).unwrap_or(0.0)
      if t1.length() > 0 && t3.length() > 0 {
        angle_forces.push((t1, t2, t3, k, th))
      }
    }
    // Periodic torsions: <Proper/> / <Improper/> lines carry up to 5
    // Fourier terms (k1..k5, phi1..phi5, periodicity1..5); every term with
    // a nonzero force constant becomes one record.
    if line.contains("<Proper") || line.contains("<Improper") {
      let t1 = omm_xml_attr(line, "type1")
      let t4 = omm_xml_attr(line, "type4")
      if t1.length() > 0 && t4.length() > 0 {
        let t2 = omm_xml_attr(line, "type2")
        let t3 = omm_xml_attr(line, "type3")
        for m = 1; m <= 5; m = m + 1 {
          let k = parse_double(omm_xml_attr(line, "k" + m.to_string())).unwrap_or(
            0.0,
          )
          if k != 0.0 {
            let phi = parse_double(omm_xml_attr(line, "phi" + m.to_string())).unwrap_or(
              0.0,
            )
            let per = parse_double(
              omm_xml_attr(line, "periodicity" + m.to_string()),
            ).unwrap_or(3.0)
            torsion_forces.push((t1, t2, t3, t4, k, phi, per.to_int()))
          }
        }
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
/// Look up bond parameters by atom type pair (symmetric). XML units are
/// preserved (k kJ/mol/nm^2, r0 nm).
pub fn omm_xml_ff_lookup_bond(
  ff : OmmXmlForceField,
  t1 : String,
  t2 : String,
) -> (Double, Double)? {
  for (a, b, k, r0) in ff.bond_forces {
    if (a == t1 && b == t2) || (a == t2 && b == t1) {
      return Some((k, r0))
    }
  }
  None
}

///|
/// Look up angle parameters by atom type triple (ends are symmetric).
/// XML units are preserved (k kJ/mol/rad^2, theta0 degrees).
pub fn omm_xml_ff_lookup_angle(
  ff : OmmXmlForceField,
  t1 : String,
  t2 : String,
  t3 : String,
) -> (Double, Double)? {
  for (a, b, c, k, th) in ff.angle_forces {
    if (a == t1 && b == t2 && c == t3) || (a == t3 && b == t2 && c == t1) {
      return Some((k, th))
    }
  }
  None
}

///|
/// Look up all periodic-torsion Fourier terms for an atom type quadruple
/// (chain is symmetric under full reversal). XML units are preserved
/// (k kJ/mol, phi degrees).
pub fn omm_xml_ff_lookup_torsions(
  ff : OmmXmlForceField,
  t1 : String,
  t2 : String,
  t3 : String,
  t4 : String,
) -> Array[(Double, Double, Int32)] {
  let out : Array[(Double, Double, Int32)] = Array(capacity=0)
  for (a, b, c, d, k, phi, n) in ff.torsion_forces {
    let fwd = a == t1 && b == t2 && c == t3 && d == t4
    let rev = a == t4 && b == t3 && c == t2 && d == t1
    if fwd || rev {
      out.push((k, phi, n))
    }
  }
  out
}

///|
/// First atom type whose element matches `element` (template fallback when
/// an atom cannot be typed from a residue template).
pub fn omm_xml_ff_type_for_element(
  ff : OmmXmlForceField,
  element : String,
) -> String {
  for t in ff.atom_types {
    if t.element == element {
      return t.name
    }
  }
  ""
}

///|
/// Element -> standard mass (amu) used when no atom-type mass exists.
fn omm_xml_ff_element_mass(element : String) -> Double {
  if element == "H" {
    1.008
  } else if element == "C" {
    12.011
  } else if element == "N" {
    14.007
  } else if element == "O" {
    15.999
  } else if element == "P" {
    30.974
  } else if element == "S" {
    32.06
  } else if element == "Na" {
    22.990
  } else if element == "Cl" {
    35.453
  } else {
    12.011
  }
}

///|
/// One atom of the molecule to parameterize: name + residue name + element.
pub(all) struct OmmXmlAtomSpec {
  name : String
  resname : String
  element : String
}

///|
/// Build a runnable OmmSystem from a parsed XML force field plus a molecule
/// given as `OmmXmlAtomSpec` list (atom name + residue name + element).
///
/// Workflow per atom: residue-template type match, then element fallback.
/// Bonds come from residue templates (within one residue); angles and
/// torsions are generated from that bond graph. Parameters use OpenMM XML
/// units and are converted to the OmmSystem unit system:
/// k kJ/mol/nm^2 -> kcal/mol/A^2 (bond, /418.4), r0 nm -> A (*10),
/// angle k kJ/mol/rad^2 -> kcal/mol/rad^2 (/4.184), theta0 deg -> rad,
/// torsion k kJ/mol -> kcal/mol (/4.184), phase deg -> rad,
/// sigma nm -> A (*10), epsilon kJ/mol -> kcal/mol (/4.184).
///
/// Fallbacks when no XML parameter exists: bond (300 kcal/mol/A^2, 1.5 A),
/// angle (50 kcal/mol/rad^2, 109.5 deg), torsion (term skipped).
pub fn omm_xml_ff_build_system(
  ff : OmmXmlForceField,
  atoms : Array[OmmXmlAtomSpec],
) -> OmmSystem {
  let n = atoms.length()
  let types : Array[String] = Array::make(n, "")
  let system = OmmSystem::new()
  // Assign atom types + particles.
  for i = 0; i < n; i = i + 1 {
    let assigned = omm_match_residue_template(ff, atoms[i].resname, [
      atoms[i].name,
    ])
    let mut ty = if assigned.length() > 0 { assigned[0].1 } else { "X" }
    if ty == "X" {
      let el_type = omm_xml_ff_type_for_element(ff, atoms[i].element)
      if el_type.length() > 0 {
        ty = el_type
      }
    }
    types[i] = ty
    // Mass: atom type table, else element mass.
    let mut mass = 0.0
    for t in ff.atom_types {
      if t.name == ty {
        mass = t.mass
        break
      }
    }
    if mass == 0.0 {
      mass = omm_xml_ff_element_mass(atoms[i].element)
    }
    // Nonbonded (XML nm/kJ -> A/kcal).
    let (q, s_nm, e_kj) = omm_lookup_nonbonded(ff, ty)
    ignore(
      system.add_particle(
        particle=OmmParticle::new(
          mass~,
          charge=q,
          sigma=s_nm * 10.0,
          epsilon=e_kj / 4.184,
        ),
      ),
    )
  }
  // Residue templates: per-residue name-indexed bond graph.
  fn template_bonds(resname : String) -> Array[(String, String)] {
    for r in ff.residues {
      if r.name.to_lower() == resname.to_lower() {
        let out : Array[(String, String)] = Array(capacity=0)
        for (f, t) in r.bonds {
          let fi = parse_double(f)
          let ti = parse_double(t)
          match (fi, ti) {
            (Some(a), Some(b)) => {
              let ai = a.to_int()
              let bi = b.to_int()
              if ai >= 0 &&
                ai < r.atoms.length() &&
                bi >= 0 &&
                bi < r.atoms.length() {
                out.push((r.atoms[ai].0, r.atoms[bi].0))
              }
            }
            _ => out.push((f, t)) // already name-based
          }
        }
        return out
      }
    }
    []
  }
  // Map (residue, atom name) -> system index.
  fn find_idx(start : Int32, resname : String, name : String) -> Int32 {
    let mut i = start
    while i < n {
      if atoms[i].resname == resname && atoms[i].name == name {
        return i
      }
      i = i + 1
    }
    -1
  }
  // Bonds + per-atom adjacency for angle/torsion generation.
  let adjacency : Array[Array[Int32]] = Array::make(n, [])
  for i = 0; i < n; i = i + 1 {
    adjacency[i] = []
  }
  let bond_pairs : Array[(Int32, Int32)] = Array(capacity=0)
  let mut scan = 0
  while scan < n {
    let resname = atoms[scan].resname
    let mut end = scan
    while end < n && atoms[end].resname == resname {
      end = end + 1
    }
    let tb = template_bonds(resname)
    for (n1, n2) in tb {
      let i1 = find_idx(scan, resname, n1)
      let i2 = find_idx(scan, resname, n2)
      if i1 >= 0 && i2 >= 0 {
        bond_pairs.push((i1, i2))
        adjacency[i1].push(i2)
        adjacency[i2].push(i1)
        // Bond parameters (with element-based type fallback already applied).
        let (k_xml, r0_xml) = match
          omm_xml_ff_lookup_bond(ff, types[i1], types[i2]) {
          Some(v) => v
          None => (1255.2, 0.15) // 300 kcal/mol/A^2, 1.5 A in XML units
        }
        system.add_bond(i1, i2, k_xml / 418.4, r0_xml * 10.0)
      }
    }
    scan = end
  }
  // Angles: every path i-j-k through the bond graph.
  for j = 0; j < n; j = j + 1 {
    let nb = adjacency[j]
    for a = 0; a < nb.length(); a = a + 1 {
      for b = a + 1; b < nb.length(); b = b + 1 {
        let i1 = nb[a]
        let i2 = nb[b]
        let (k_xml, th_xml) = match
          omm_xml_ff_lookup_angle(ff, types[i1], types[j], types[i2]) {
          Some(v) => v
          None => (209.2, 109.5) // 50 kcal/mol/rad^2, 109.5 deg
        }
        system.add_angle(
          i1,
          j,
          i2,
          k_xml / 4.184,
          th_xml * 3.141592653589793 / 180.0,
        )
      }
    }
  }
  // Torsions: every i-j-k-l with j-k a bond (skip terms without parameters).
  for j = 0; j < n; j = j + 1 {
    for kx = 0; kx < adjacency[j].length(); kx = kx + 1 {
      let kk = adjacency[j][kx]
      if kk <= j {
        continue // each center bond once
      }
      for a = 0; a < adjacency[j].length(); a = a + 1 {
        let i1 = adjacency[j][a]
        if i1 == kk {
          continue
        }
        for b = 0; b < adjacency[kk].length(); b = b + 1 {
          let i2 = adjacency[kk][b]
          if i2 == j {
            continue
          }
          let terms = omm_xml_ff_lookup_torsions(
            ff,
            types[i1],
            types[j],
            types[kk],
            types[i2],
          )
          for (k_xml, phi_xml, per) in terms {
            system.add_dihedral(
              i1,
              j,
              kk,
              i2,
              k_xml / 4.184,
              phi_xml * 3.141592653589793 / 180.0,
              per,
            )
          }
        }
      }
    }
  }
  system
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
