///|
/// Bio.PDB.PDBMLParser demo.
///
/// Parses a small PDBML (PDB XML) document and prints the header metadata plus
/// the Model -> Chain -> Residue -> Atom hierarchy.

///|
fn main raise {
  let xml = "<PDBx:datablock xmlns:PDBx='http://pdbml.pdb.org/schema/pdbx-v50.xsd'>" +
  "<PDBx:entryCategory><PDBx:entry id='1ABC'/></PDBx:entryCategory>" +
  "<PDBx:structCategory><PDBx:struct><PDBx:title>TEST STRUCTURE</PDBx:title></PDBx:struct></PDBx:structCategory>" +
  "<PDBx:struct_keywordsCategory><PDBx:struct_keywords><PDBx:text>demo keyword</PDBx:text></PDBx:struct_keywords></PDBx:struct_keywordsCategory>" +
  "<PDBx:pdbx_database_statusCategory><PDBx:pdbx_database_status><PDBx:recvd_initial_deposition_date>2020-01-02</PDBx:recvd_initial_deposition_date></PDBx:pdbx_database_status></PDBx:pdbx_database_statusCategory>" +
  "<PDBx:exptlCategory><PDBx:exptl method='X-RAY DIFFRACTION'/></PDBx:exptlCategory>" +
  "<PDBx:refineCategory><PDBx:refine><PDBx:ls_d_res_high>1.50</PDBx:ls_d_res_high></PDBx:refine></PDBx:refineCategory>" +
  "<PDBx:atom_siteCategory>" +
  "<PDBx:atom_site id='1'><PDBx:group_PDB>ATOM</PDBx:group_PDB><PDBx:label_comp_id>ALA</PDBx:label_comp_id><PDBx:auth_seq_id>1</PDBx:auth_seq_id><PDBx:pdbx_PDB_ins_code>?</PDBx:pdbx_PDB_ins_code><PDBx:auth_asym_id>A</PDBx:auth_asym_id><PDBx:pdbx_PDB_model_num>1</PDBx:pdbx_PDB_model_num><PDBx:label_atom_id>CA</PDBx:label_atom_id><PDBx:type_symbol>C</PDBx:type_symbol><PDBx:label_alt_id>.</PDBx:label_alt_id><PDBx:Cartn_x>1.1</PDBx:Cartn_x><PDBx:Cartn_y>2.2</PDBx:Cartn_y><PDBx:Cartn_z>3.3</PDBx:Cartn_z><PDBx:occupancy>1.00</PDBx:occupancy><PDBx:B_iso_or_equiv>20.5</PDBx:B_iso_or_equiv></PDBx:atom_site>" +
  "<PDBx:atom_site id='2'><PDBx:group_PDB>HETATM</PDBx:group_PDB><PDBx:label_comp_id>HOH</PDBx:label_comp_id><PDBx:auth_seq_id>2</PDBx:auth_seq_id><PDBx:pdbx_PDB_ins_code>?</PDBx:pdbx_PDB_ins_code><PDBx:auth_asym_id>B</PDBx:auth_asym_id><PDBx:pdbx_PDB_model_num>1</PDBx:pdbx_PDB_model_num><PDBx:label_atom_id>O</PDBx:label_atom_id><PDBx:type_symbol>O</PDBx:type_symbol><PDBx:label_alt_id>.</PDBx:label_alt_id><PDBx:Cartn_x>4.0</PDBx:Cartn_x><PDBx:Cartn_y>5.0</PDBx:Cartn_y><PDBx:Cartn_z>6.0</PDBx:Cartn_z><PDBx:occupancy>0.50</PDBx:occupancy><PDBx:B_iso_or_equiv>30.0</PDBx:B_iso_or_equiv></PDBx:atom_site>" +
  "</PDBx:atom_siteCategory></PDBx:datablock>"

  let s = @src.pdbml_get_structure(xml)
  println("idcode       : \{s.header.idcode}")
  println("name         : \{s.header.name}")
  println("head         : \{s.header.head}")
  println("deposition   : \{s.header.deposition_date}")
  println("method       : \{s.header.structure_method}")
  println("resolution   : \{s.header.resolution}")

  for m in s.models {
    println("model serial : \{m.serial_number}")
    for c in m.chains {
      println("  chain      : \{c.chain_id}")
      for r in c.residues {
        println("    residue  : \{r.component_id} hetero='\{r.hetero_field}' seq=\{r.sequence_id} ins='\{r.insertion_code}'")
        for a in r.atoms {
          println("      atom   : \{a.name} serial=\{a.serial_number} element=\{a.element} coord=(\{a.coord[0]},\{a.coord[1]},\{a.coord[2]})")
        }
      }
    }
  }
}