// Melting temperature via nearest-neighbor thermodynamics - faithful port of
// Biopython's Bio.SeqUtils.MeltingTemp (Tm_Wallace / Tm_GC / Tm_NN,
// salt_correction methods 1-7, chem_correction for DMSO/formamide).
//
// Tm_NN implements the full Biopython algorithm: eight published NN
// thermodynamic table sets (DNA Breslauer'86/Sugimoto'96/Allawi-SantaLucia'97/
// SantaLucia-Hicks'04; RNA Freier'86/Xia'98/Chen'12; RNA-DNA Sugimoto'95),
// initiation terms, terminal/internal mismatch tables, dangling-end tables,
// self-complement symmetry correction, dnac1/dnac2 strand concentrations,
// salt corrections (Schildkraut-Lifson, Wetmur, SantaLucia entropy, Owczarzy
// 2004/2008 with Mg2+ and dNTP binding) and DMSO/formamide corrections.

///|
/// One thermodynamic table: keys are initiation names ("init", "sym",
/// "init_A/T", ...) or duplex NN keys ("AA/TT", ...); values are (dH in
/// kcal/mol, dS in cal/(mol*K)).
pub type ThermoTable = Map[String, (Double, Double)]

///|
fn thermo(pairs : Array[(String, Double, Double)]) -> ThermoTable {
  let m : ThermoTable = Map([])
  for p in pairs {
    m.set(p.0, (p.1, p.2))
  }
  m
}

///|
/// Return one of the built-in thermodynamic tables by Biopython name:
/// "DNA_NN1".."DNA_NN4", "RNA_NN1".."RNA_NN3", "R_DNA_NN1", "DNA_IMM1",
/// "DNA_TMM1", "DNA_DE1", "RNA_DE1".
pub fn tm_thermo_table(name : String) -> ThermoTable raise SeqError {
  match name {
    "DNA_NN1" =>
      thermo([
        ("init", 0.0, 0.0),
        ("init_A/T", 0.0, 0.0),
        ("init_G/C", 0.0, 0.0),
        ("init_oneG/C", 0.0, -16.8),
        ("init_allA/T", 0.0, -20.1),
        ("init_5T/A", 0.0, 0.0),
        ("sym", 0.0, -1.3),
        ("AA/TT", -9.1, -24.0),
        ("AT/TA", -8.6, -23.9),
        ("TA/AT", -6.0, -16.9),
        ("CA/GT", -5.8, -12.9),
        ("GT/CA", -6.5, -17.3),
        ("CT/GA", -7.8, -20.8),
        ("GA/CT", -5.6, -13.5),
        ("CG/GC", -11.9, -27.8),
        ("GC/CG", -11.1, -26.7),
        ("GG/CC", -11.0, -26.6),
      ])
    "DNA_NN2" =>
      thermo([
        ("init", 0.6, -9.0),
        ("init_A/T", 0.0, 0.0),
        ("init_G/C", 0.0, 0.0),
        ("init_oneG/C", 0.0, 0.0),
        ("init_allA/T", 0.0, 0.0),
        ("init_5T/A", 0.0, 0.0),
        ("sym", 0.0, -1.4),
        ("AA/TT", -8.0, -21.9),
        ("AT/TA", -5.6, -15.2),
        ("TA/AT", -6.6, -18.4),
        ("CA/GT", -8.2, -21.0),
        ("GT/CA", -9.4, -25.5),
        ("CT/GA", -6.6, -16.4),
        ("GA/CT", -8.8, -23.5),
        ("CG/GC", -11.8, -29.0),
        ("GC/CG", -10.5, -26.4),
        ("GG/CC", -10.9, -28.4),
      ])
    "DNA_NN3" =>
      thermo([
        ("init", 0.0, 0.0),
        ("init_A/T", 2.3, 4.1),
        ("init_G/C", 0.1, -2.8),
        ("init_oneG/C", 0.0, 0.0),
        ("init_allA/T", 0.0, 0.0),
        ("init_5T/A", 0.0, 0.0),
        ("sym", 0.0, -1.4),
        ("AA/TT", -7.9, -22.2),
        ("AT/TA", -7.2, -20.4),
        ("TA/AT", -7.2, -21.3),
        ("CA/GT", -8.5, -22.7),
        ("GT/CA", -8.4, -22.4),
        ("CT/GA", -7.8, -21.0),
        ("GA/CT", -8.2, -22.2),
        ("CG/GC", -10.6, -27.2),
        ("GC/CG", -9.8, -24.4),
        ("GG/CC", -8.0, -19.9),
      ])
    "DNA_NN4" =>
      thermo([
        ("init", 0.2, -5.7),
        ("init_A/T", 2.2, 6.9),
        ("init_G/C", 0.0, 0.0),
        ("init_oneG/C", 0.0, 0.0),
        ("init_allA/T", 0.0, 0.0),
        ("init_5T/A", 0.0, 0.0),
        ("sym", 0.0, -1.4),
        ("AA/TT", -7.6, -21.3),
        ("AT/TA", -7.2, -20.4),
        ("TA/AT", -7.2, -21.3),
        ("CA/GT", -8.5, -22.7),
        ("GT/CA", -8.4, -22.4),
        ("CT/GA", -7.8, -21.0),
        ("GA/CT", -8.2, -22.2),
        ("CG/GC", -10.6, -27.2),
        ("GC/CG", -9.8, -24.4),
        ("GG/CC", -8.0, -19.9),
      ])
    "RNA_NN1" =>
      thermo([
        ("init", 0.0, -10.8),
        ("init_A/T", 0.0, 0.0),
        ("init_G/C", 0.0, 0.0),
        ("init_oneG/C", 0.0, 0.0),
        ("init_allA/T", 0.0, 0.0),
        ("init_5T/A", 0.0, 0.0),
        ("sym", 0.0, -1.4),
        ("AA/TT", -6.6, -18.4),
        ("AT/TA", -5.7, -15.5),
        ("TA/AT", -8.1, -22.6),
        ("CA/GT", -10.5, -27.8),
        ("GT/CA", -10.2, -26.2),
        ("CT/GA", -7.6, -19.2),
        ("GA/CT", -13.3, -35.5),
        ("CG/GC", -8.0, -19.4),
        ("GC/CG", -14.2, -34.9),
        ("GG/CC", -12.2, -29.7),
      ])
    "RNA_NN2" =>
      thermo([
        ("init", 3.61, -1.5),
        ("init_A/T", 3.72, 10.5),
        ("init_G/C", 0.0, 0.0),
        ("init_oneG/C", 0.0, 0.0),
        ("init_allA/T", 0.0, 0.0),
        ("init_5T/A", 0.0, 0.0),
        ("sym", 0.0, -1.4),
        ("AA/TT", -6.82, -19.0),
        ("AT/TA", -9.38, -26.7),
        ("TA/AT", -7.69, -20.5),
        ("CA/GT", -10.44, -26.9),
        ("GT/CA", -11.40, -29.5),
        ("CT/GA", -10.48, -27.1),
        ("GA/CT", -12.44, -32.5),
        ("CG/GC", -10.64, -26.7),
        ("GC/CG", -14.88, -36.9),
        ("GG/CC", -13.39, -32.7),
      ])
    "RNA_NN3" =>
      thermo([
        ("init", 6.40, 6.99),
        ("init_A/T", 3.85, 11.04),
        ("init_G/C", 0.0, 0.0),
        ("init_oneG/C", 0.0, 0.0),
        ("init_allA/T", 0.0, 0.0),
        ("init_5T/A", 0.0, 0.0),
        ("sym", 0.0, -1.4),
        ("AA/TT", -7.09, -19.8),
        ("AT/TA", -9.11, -25.8),
        ("TA/AT", -8.50, -22.9),
        ("CA/GT", -11.03, -28.8),
        ("GT/CA", -11.98, -31.3),
        ("CT/GA", -10.90, -28.5),
        ("GA/CT", -13.21, -34.9),
        ("CG/GC", -10.88, -27.4),
        ("GC/CG", -16.04, -40.6),
        ("GG/CC", -14.18, -35.0),
        ("GT/TG", -13.83, -46.9),
        ("GG/TT", -17.82, -56.7),
        ("AG/TT", -3.96, -11.6),
        ("TG/AT", -0.96, -1.8),
        ("TT/AG", -10.38, -31.8),
        ("TG/GT", -12.64, -38.9),
        ("AT/TG", -7.39, -21.0),
        ("CG/GT", -5.56, -13.9),
        ("CT/GG", -9.44, -24.7),
        ("GG/CT", -7.03, -16.8),
        ("GT/CG", -11.09, -28.8),
      ])
    "R_DNA_NN1" =>
      thermo([
        ("init", 1.9, -3.9),
        ("init_A/T", 0.0, 0.0),
        ("init_G/C", 0.0, 0.0),
        ("init_oneG/C", 0.0, 0.0),
        ("init_allA/T", 0.0, 0.0),
        ("init_5T/A", 0.0, 0.0),
        ("sym", 0.0, 0.0),
        ("TT/AA", -11.5, -36.4),
        ("GT/CA", -7.8, -21.6),
        ("CT/GA", -7.0, -19.7),
        ("AT/TA", -8.3, -23.9),
        ("TG/AC", -10.4, -28.4),
        ("GG/CC", -12.8, -31.9),
        ("CG/GC", -16.3, -47.1),
        ("AG/TC", -9.1, -23.5),
        ("TC/AG", -8.6, -22.9),
        ("GC/CG", -8.0, -17.1),
        ("CC/GG", -9.3, -23.2),
        ("AC/TG", -5.9, -12.3),
        ("TA/AT", -7.8, -23.2),
        ("GA/CT", -5.5, -13.5),
        ("CA/GT", -9.0, -26.1),
        ("AA/TT", -7.8, -21.9),
      ])
    "DNA_IMM1" =>
      thermo([
        ("AG/TT", 1.0, 0.9),
        ("AT/TG", -2.5, -8.3),
        ("CG/GT", -4.1, -11.7),
        ("CT/GG", -2.8, -8.0),
        ("GG/CT", 3.3, 10.4),
        ("GG/TT", 5.8, 16.3),
        ("GT/CG", -4.4, -12.3),
        ("GT/TG", 4.1, 9.5),
        ("TG/AT", -0.1, -1.7),
        ("TG/GT", -1.4, -6.2),
        ("TT/AG", -1.3, -5.3),
        ("AA/TG", -0.6, -2.3),
        ("AG/TA", -0.7, -2.3),
        ("CA/GG", -0.7, -2.3),
        ("CG/GA", -4.0, -13.2),
        ("GA/CG", -0.6, -1.0),
        ("GG/CA", 0.5, 3.2),
        ("TA/AG", 0.7, 0.7),
        ("TG/AA", 3.0, 7.4),
        ("AC/TT", 0.7, 0.2),
        ("AT/TC", -1.2, -6.2),
        ("CC/GT", -0.8, -4.5),
        ("CT/GC", -1.5, -6.1),
        ("GC/CT", 2.3, 5.4),
        ("GT/CC", 5.2, 13.5),
        ("TC/AT", 1.2, 0.7),
        ("TT/AC", 1.0, 0.7),
        ("AA/TC", 2.3, 4.6),
        ("AC/TA", 5.3, 14.6),
        ("CA/GC", 1.9, 3.7),
        ("CC/GA", 0.6, -0.6),
        ("GA/CC", 5.2, 14.2),
        ("GC/CA", -0.7, -3.8),
        ("TA/AC", 3.4, 8.0),
        ("TC/AA", 7.6, 20.2),
        ("AA/TA", 1.2, 1.7),
        ("CA/GA", -0.9, -4.2),
        ("GA/CA", -2.9, -9.8),
        ("TA/AA", 4.7, 12.9),
        ("AC/TC", 0.0, -4.4),
        ("CC/GC", -1.5, -7.2),
        ("GC/CC", 3.6, 8.9),
        ("TC/AC", 6.1, 16.4),
        ("AG/TG", -3.1, -9.5),
        ("CG/GG", -4.9, -15.3),
        ("GG/CG", -6.0, -15.8),
        ("TG/AG", 1.6, 3.6),
        ("AT/TT", -2.7, -10.8),
        ("CT/GT", -5.0, -15.8),
        ("GT/CT", -2.2, -8.4),
        ("TT/AT", 0.2, -1.5),
        ("AI/TC", -8.9, -25.5),
        ("TI/AC", -5.9, -17.4),
        ("AC/TI", -8.8, -25.4),
        ("TC/AI", -4.9, -13.9),
        ("CI/GC", -5.4, -13.7),
        ("GI/CC", -6.8, -19.1),
        ("CC/GI", -8.3, -23.8),
        ("GC/CI", -5.0, -12.6),
        ("AI/TA", -8.3, -25.0),
        ("TI/AA", -3.4, -11.2),
        ("AA/TI", -0.7, -2.6),
        ("TA/AI", -1.3, -4.6),
        ("CI/GA", 2.6, 8.9),
        ("GI/CA", -7.8, -21.1),
        ("CA/GI", -7.0, -20.0),
        ("GA/CI", -7.6, -20.2),
        ("AI/TT", 0.49, -0.7),
        ("TI/AT", -6.5, -22.0),
        ("AT/TI", -5.6, -18.7),
        ("TT/AI", -0.8, -4.3),
        ("CI/GT", -1.0, -2.4),
        ("GI/CT", -3.5, -10.6),
        ("CT/GI", 0.1, -1.0),
        ("GT/CI", -4.3, -12.1),
        ("AI/TG", -4.9, -15.8),
        ("TI/AG", -1.9, -8.5),
        ("AG/TI", 0.1, -1.8),
        ("TG/AI", 1.0, 1.0),
        ("CI/GG", 7.1, 21.3),
        ("GI/CG", -1.1, -3.2),
        ("CG/GI", 5.8, 16.9),
        ("GG/CI", -7.6, -22.0),
        ("AI/TI", -3.3, -11.9),
        ("TI/AI", 0.1, -2.3),
        ("CI/GI", 1.3, 3.0),
        ("GI/CI", -0.5, -1.3),
      ])
    "DNA_TMM1" =>
      thermo([
        ("AA/TA", -3.1, -7.8),
        ("TA/AA", -2.5, -6.3),
        ("CA/GA", -4.3, -10.7),
        ("GA/CA", -8.0, -22.5),
        ("AC/TC", -0.1, 0.5),
        ("TC/AC", -0.7, -1.3),
        ("CC/GC", -2.1, -5.1),
        ("GC/CC", -3.9, -10.6),
        ("AG/TG", -1.1, -2.1),
        ("TG/AG", -1.1, -2.7),
        ("CG/GG", -3.8, -9.5),
        ("GG/CG", -0.7, -19.2),
        ("AT/TT", -2.4, -6.5),
        ("TT/AT", -3.2, -8.9),
        ("CT/GT", -6.1, -16.9),
        ("GT/CT", -7.4, -21.2),
        ("AA/TC", -1.6, -4.0),
        ("AC/TA", -1.8, -3.8),
        ("CA/GC", -2.6, -5.9),
        ("CC/GA", -2.7, -6.0),
        ("GA/CC", -5.0, -13.8),
        ("GC/CA", -3.2, -7.1),
        ("TA/AC", -2.3, -5.9),
        ("TC/AA", -2.7, -7.0),
        ("AC/TT", -0.9, -1.7),
        ("AT/TC", -2.3, -6.3),
        ("CC/GT", -3.2, -8.0),
        ("CT/GC", -3.9, -10.6),
        ("GC/CT", -4.9, -13.5),
        ("GT/CC", -3.0, -7.8),
        ("TC/AT", -2.5, -6.3),
        ("TT/AC", -0.7, -1.2),
        ("AA/TG", -1.9, -4.4),
        ("AG/TA", -2.5, -5.9),
        ("CA/GG", -3.9, -9.6),
        ("CG/GA", -6.0, -15.5),
        ("GA/CG", -4.3, -11.1),
        ("GG/CA", -4.6, -11.4),
        ("TA/AG", -2.0, -4.7),
        ("TG/AA", -2.4, -5.8),
        ("AG/TT", -3.2, -8.7),
        ("AT/TG", -3.5, -9.4),
        ("CG/GT", -3.8, -9.0),
        ("CT/GG", -6.6, -18.7),
        ("GG/CT", -5.7, -15.9),
        ("GT/CG", -5.9, -16.1),
        ("TG/AT", -3.9, -10.5),
        ("TT/AG", -3.6, -9.8),
      ])
    "DNA_DE1" =>
      thermo([
        ("AA/.T", 0.2, 2.3),
        ("AC/.G", -6.3, -17.1),
        ("AG/.C", -3.7, -10.0),
        ("AT/.A", -2.9, -7.6),
        ("CA/.T", 0.6, 3.3),
        ("CC/.G", -4.4, -12.6),
        ("CG/.C", -4.0, -11.9),
        ("CT/.A", -4.1, -13.0),
        ("GA/.T", -1.1, -1.6),
        ("GC/.G", -5.1, -14.0),
        ("GG/.C", -3.9, -10.9),
        ("GT/.A", -4.2, -15.0),
        ("TA/.T", -6.9, -20.0),
        ("TC/.G", -4.0, -10.9),
        ("TG/.C", -4.9, -13.8),
        ("TT/.A", -0.2, -0.5),
        (".A/AT", -0.7, -0.8),
        (".C/AG", -2.1, -3.9),
        (".G/AC", -5.9, -16.5),
        (".T/AA", -0.5, -1.1),
        (".A/CT", 4.4, 14.9),
        (".C/CG", -0.2, -0.1),
        (".G/CC", -2.6, -7.4),
        (".T/CA", 4.7, 14.2),
        (".A/GT", -1.6, -3.6),
        (".C/GG", -3.9, -11.2),
        (".G/GC", -3.2, -10.4),
        (".T/GA", -4.1, -13.1),
        (".A/TT", 2.9, 10.4),
        (".C/TG", -4.4, -13.1),
        (".G/TC", -5.2, -15.0),
        (".T/TA", -3.8, -12.6),
      ])
    "RNA_DE1" =>
      thermo([
        (".T/AA", -4.9, -13.2),
        (".T/CA", -0.9, -1.3),
        (".T/GA", -5.5, -15.1),
        (".T/TA", -2.3, -5.5),
        (".G/AC", -9.0, -23.5),
        (".G/CC", -4.1, -10.6),
        (".G/GC", -8.6, -22.2),
        (".G/TC", -7.5, -20.31),
        (".C/AG", -7.4, -20.3),
        (".C/CG", -2.8, -7.7),
        (".C/GG", -6.4, -16.4),
        (".C/TG", -3.6, -9.7),
        (".T/AG", -4.9, -13.2),
        (".T/CG", -0.9, -1.3),
        (".T/GG", -5.5, -15.1),
        (".T/TG", -2.3, -5.5),
        (".A/AT", -5.7, -16.1),
        (".A/CT", -0.7, -1.9),
        (".A/GT", -5.8, -16.4),
        (".A/TT", -2.2, -6.8),
        (".G/AT", -5.7, -16.1),
        (".G/CT", -0.7, -1.9),
        (".G/GT", -5.8, -16.4),
        (".G/TT", -2.2, -6.8),
        ("AT/.A", -0.5, -0.6),
        ("CT/.A", 6.9, 22.6),
        ("GT/.A", 0.6, 2.6),
        ("TT/.A", 0.6, 2.6),
        ("AG/.C", -1.6, -4.5),
        ("CG/.C", 0.7, 3.2),
        ("GG/.C", -4.6, -14.8),
        ("TG/.C", -0.4, -1.3),
        ("AC/.G", -2.4, -6.1),
        ("CC/.G", 3.3, 11.6),
        ("GC/.G", 0.8, 3.2),
        ("TC/.G", -1.4, -4.2),
        ("AT/.G", -0.5, -0.6),
        ("CT/.G", 6.9, 22.6),
        ("GT/.G", 0.6, 2.6),
        ("TT/.G", 0.6, 2.6),
        ("AA/.T", 1.6, 6.1),
        ("CA/.T", 2.2, 8.1),
        ("GA/.T", 0.7, 3.5),
        ("TA/.T", 3.1, 10.6),
        ("AG/.T", 1.6, 6.1),
        ("CG/.T", 2.2, 8.1),
        ("GG/.T", 0.7, 3.5),
        ("TG/.T", 3.1, 10.6),
      ])
    _ => raise SeqError("unknown thermodynamic table: \{name}")
  }
}

///|
/// Create a new thermodynamic table, either from a built-in table (`base`
/// name, "" for a blank zero table) with overridden entries. Port of
/// Biopython `make_table`.
pub fn tm_make_table(
  base? : String = "",
  values? : ThermoTable = Map([]),
) -> ThermoTable raise SeqError {
  let table : ThermoTable = if base == "" {
    thermo([
      ("init", 0.0, 0.0),
      ("init_A/T", 0.0, 0.0),
      ("init_G/C", 0.0, 0.0),
      ("init_oneG/C", 0.0, 0.0),
      ("init_allA/T", 0.0, 0.0),
      ("init_5T/A", 0.0, 0.0),
      ("sym", 0.0, 0.0),
      ("AA/TT", 0.0, 0.0),
      ("AT/TA", 0.0, 0.0),
      ("TA/AT", 0.0, 0.0),
      ("CA/GT", 0.0, 0.0),
      ("GT/CA", 0.0, 0.0),
      ("CT/GA", 0.0, 0.0),
      ("GA/CT", 0.0, 0.0),
      ("CG/GC", 0.0, 0.0),
      ("GC/CG", 0.0, 0.0),
      ("GG/CC", 0.0, 0.0),
    ])
  } else {
    tm_thermo_table(base)
  }
  for k, v in values {
    table.set(k, v)
  }
  table
}

///|
fn mt_dots(n : Int) -> String {
  let sb = StringBuilder()
  let mut i = 0
  while i < n {
    sb.write_string(".")
    i = i + 1
  }
  sb.to_string()
}

///|
/// Complement used when no c_seq is given (same semantics as Biopython,
/// which complements the raw input before the Tm_NN checks).
fn tm_complement(seq : String) -> String {
  let out = StringBuilder()
  let n = seq.length()
  let mut i = 0
  while i < n {
    let c = substring(seq, i, 1)
    out.write_string(
      match c {
        "A" => "T"
        "T" => "A"
        "U" => "A"
        "G" => "C"
        "C" => "G"
        _ => c
      },
    )
    i = i + 1
  }
  out.to_string()
}

///|
/// Normalize a sequence for a Tm method: drop whitespace, uppercase,
/// back-transcribe RNA (U->T), then keep only bases valid for the method.
/// Port of Biopython MeltingTemp._check.
fn tm_check(seq : String, method_name : String) -> String {
  // remove all whitespace then uppercase
  let s = seq.to_upper()
  let kept = StringBuilder()
  let baseset : Array[String] = if method_name == "Tm_GC" {
    [
      "A", "B", "C", "D", "G", "H", "I", "K", "M", "N", "R", "S", "T", "V", "W",
      "X", "Y",
    ]
  } else if method_name == "Tm_NN" {
    ["A", "C", "G", "T", "I"]
  } else {
    []
  }
  let n = s.length()
  let mut i = 0
  while i < n {
    let c = substring(s, i, 1)
    if c == " " || c == "\t" || c == "\n" || c == "\r" {
      i = i + 1
      continue
    }
    // back-transcribe
    let base = if c == "U" { "T" } else { c }
    if method_name == "Tm_Wallace" {
      kept.write_string(base)
    } else {
      for b in baseset {
        if b == base {
          kept.write_string(base)
          break
        }
      }
    }
    i = i + 1
  }
  kept.to_string()
}

///|
fn mt_reverse(s : String) -> String {
  let out = StringBuilder()
  let mut i = s.length() - 1
  while i >= 0 {
    out.write_string(substring(s, i, 1))
    i = i - 1
  }
  out.to_string()
}

///|
/// Salt correction term for Tm calculations. Port of Biopython
/// `salt_correction`; methods 1-4 return an additive Tm term, method 5 an
/// additive entropy (cal/(mol*K)) term, methods 6-7 a term to add to
/// 1/Tm(K). Na, K, Tris, Mg, dNTPs are millimolar concentrations.
pub fn tm_salt_correction(
  na? : Double = 0.0,
  k? : Double = 0.0,
  tris? : Double = 0.0,
  mg? : Double = 0.0,
  dntps? : Double = 0.0,
  salt_method? : Int = 1,
  seq? : String = "",
) -> Double raise SeqError {
  if salt_method == 5 || salt_method == 6 || salt_method == 7 {
    if seq == "" {
      raise SeqError(
        "sequence is missing (is needed to calculate GC content or sequence length).",
      )
    }
  }
  if salt_method > 7 {
    raise SeqError("Allowed values for parameter 'salt_method' are 1-7.")
  }
  if salt_method == 0 {
    return 0.0
  }
  let mut mon_val = na + k + tris / 2.0 // all millimolar
  let mg_lower = mg * 1.0e-3 // molar
  if k + mg + tris + dntps > 0.0 && salt_method != 7 && dntps < mg {
    // dNTPs bind Mg2+ strongly; sodium-equivalent per von Ahsen (2001)
    mon_val = mon_val + 120.0 * (mg - dntps).sqrt()
  }
  let mon = mon_val * 1.0e-3
  if salt_method >= 1 && salt_method <= 6 && mon == 0.0 {
    raise SeqError(
      "Total ion concentration of zero is not allowed in this method.",
    )
  }
  let ln10 = @math.ln(10.0)
  let log10 = fn(x : Double) { @math.ln(x) / ln10 }
  if salt_method == 1 {
    16.6 * log10(mon)
  } else if salt_method == 2 {
    16.6 * log10(mon / (1.0 + 0.7 * mon))
  } else if salt_method == 3 {
    12.5 * log10(mon)
  } else if salt_method == 4 {
    11.7 * log10(mon)
  } else if salt_method == 5 {
    0.368 * (seq.length() - 1).to_double() * @math.ln(mon)
  } else if salt_method == 6 {
    let gcf = gc_fraction(seq, ambiguous="ignore")
    let ln_m = @math.ln(mon)
    (4.29 * gcf - 3.95) * 1.0e-5 * ln_m + 9.40e-6 * ln_m * ln_m
  } else {
    // salt_method 7: Owczarzy et al. (2008)
    let mut a = 3.92
    let b = -0.911
    let c = 6.26
    let mut d = 1.42
    let e = -48.2
    let f = 52.5
    let mut g = 8.31
    let mut mg_free = mg_lower
    if dntps > 0.0 {
      let dntps_m = dntps * 1.0e-3
      let ka = 3.0e4 // dissociation constant for Mg:dNTP
      let q = ka * dntps_m - ka * mg_lower + 1.0
      mg_free = (-q + (q * q + 4.0 * ka * mg_lower).sqrt()) / (2.0 * ka)
    }
    if mon_val > 0.0 {
      let ratio = mg_free.sqrt() / mon
      let ln_m = @math.ln(mon)
      if ratio < 0.22 {
        let gcf = gc_fraction(seq, ambiguous="ignore")
        return (4.29 * gcf - 3.95) * 1.0e-5 * ln_m + 9.40e-6 * ln_m * ln_m
      } else if ratio < 6.0 {
        a = 3.92 * (0.843 - 0.352 * mon.sqrt() * ln_m)
        d = 1.42 * (1.279 - 4.03e-3 * ln_m - 8.03e-3 * ln_m * ln_m)
        g = 8.31 * (0.486 - 0.258 * ln_m + 5.25e-3 * ln_m * ln_m * ln_m)
      }
    } else if mg_free <= 0.0 {
      raise SeqError(
        "Total ion concentration of zero is not allowed in this method.",
      )
    }
    let gcf = gc_fraction(seq, ambiguous="ignore")
    let ln_mg = @math.ln(mg_free)
    (
      a +
      b * ln_mg +
      gcf * (c + d * ln_mg) +
      1.0 /
      (2.0 * (seq.length() - 1).to_double()) *
      (e + f * ln_mg + g * ln_mg * ln_mg)
    ) *
    1.0e-5
  }
}

///|
/// Chemical additive correction (DMSO / formamide) for a Tm in degrees C.
/// Port of Biopython `chem_correction`. DMSO is given in percent, formamide
/// in percent (fmdmethod=1) or molar (fmdmethod=2, GC in percent).
pub fn tm_chem_correction(
  melting_temp : Double,
  dmso? : Double = 0.0,
  fmd? : Double = 0.0,
  dmso_factor? : Double = 0.75,
  fmd_factor? : Double = 0.65,
  fmd_method? : Int = 1,
  gc? : Double = -1.0,
) -> Double raise SeqError {
  let mut tm = melting_temp
  if dmso != 0.0 {
    tm = tm - dmso_factor * dmso
  }
  if fmd != 0.0 {
    if fmd_method == 1 {
      tm = tm - fmd_factor * fmd
    } else if fmd_method == 2 {
      if gc < 0.0 {
        raise SeqError("'GC' is missing or negative")
      }
      tm = tm + (0.453 * (gc / 100.0) - 2.88) * fmd
    } else {
      raise SeqError("'fmdmethod' must be 1 or 2")
    }
  }
  tm
}

///|
/// Options for `tm_wallace_full` / `tm_gc_full` / `tm_nn_full`, mirroring the
/// keyword arguments of the Biopython Tm methods.
struct TmOptions {
  check : Bool
  strict : Bool
  c_seq : String
  shift : Int
  nn_table : String
  tmm_table : String
  imm_table : String
  de_table : String
  nn_custom : ThermoTable
  tmm_custom : ThermoTable
  imm_custom : ThermoTable
  de_custom : ThermoTable
  dnac1 : Double
  dnac2 : Double
  selfcomp : Bool
  na : Double
  k : Double
  tris : Double
  mg : Double
  dntps : Double
  saltcorr : Int
}

///|
/// Default options (Biopython defaults: Na=50 mM, saltcorr=5 for Tm_NN).
pub fn TmOptions::default() -> TmOptions {
  {
    check: true,
    strict: true,
    c_seq: "",
    shift: 0,
    nn_table: "DNA_NN3",
    tmm_table: "DNA_TMM1",
    imm_table: "DNA_IMM1",
    de_table: "DNA_DE1",
    nn_custom: Map([]),
    tmm_custom: Map([]),
    imm_custom: Map([]),
    de_custom: Map([]),
    dnac1: 25.0,
    dnac2: 25.0,
    selfcomp: false,
    na: 50.0,
    k: 0.0,
    tris: 0.0,
    mg: 0.0,
    dntps: 0.0,
    saltcorr: 5,
  }
}

///|
/// Construct options with explicit overrides; unset parameters keep the
/// Biopython defaults (see `TmOptions::default`).
pub fn TmOptions::of(
  check? : Bool = true,
  strict? : Bool = true,
  c_seq? : String = "",
  shift? : Int = 0,
  nn_table? : String = "DNA_NN3",
  tmm_table? : String = "DNA_TMM1",
  imm_table? : String = "DNA_IMM1",
  de_table? : String = "DNA_DE1",
  nn_custom? : ThermoTable = Map([]),
  tmm_custom? : ThermoTable = Map([]),
  imm_custom? : ThermoTable = Map([]),
  de_custom? : ThermoTable = Map([]),
  dnac1? : Double = 25.0,
  dnac2? : Double = 25.0,
  selfcomp? : Bool = false,
  na? : Double = 50.0,
  k? : Double = 0.0,
  tris? : Double = 0.0,
  mg? : Double = 0.0,
  dntps? : Double = 0.0,
  saltcorr? : Int = 5,
) -> TmOptions {
  {
    check,
    strict,
    c_seq,
    shift,
    nn_table,
    tmm_table,
    imm_table,
    de_table,
    nn_custom,
    tmm_custom,
    imm_custom,
    de_custom,
    dnac1,
    dnac2,
    selfcomp,
    na,
    k,
    tris,
    mg,
    dntps,
    saltcorr,
  }
}

///|
/// Resolve a table: custom override wins, otherwise load the built-in name.
fn tm_resolve(
  custom : ThermoTable,
  name : String,
) -> ThermoTable raise SeqError {
  if custom.length() > 0 {
    custom
  } else {
    tm_thermo_table(name)
  }
}

///|
/// Look up a duplex key in a table, trying the key and its reverse; returns
/// Some((dH, dS)) or None.
fn tm_table_lookup(table : ThermoTable, key : String) -> (Double, Double)? {
  match table.get(key) {
    Some(v) => Some(v)
    None => table.get(mt_reverse(key))
  }
}

///|
/// Raise or silently ignore a missing-neighbor condition (strict mode).
fn tm_key_error(neighbors : String, strict : Bool) -> Unit raise SeqError {
  if strict {
    raise SeqError(
      "no thermodynamic data for neighbors '\{neighbors}' available",
    )
  }
}

///|
/// Tm via the Wallace rule of thumb: 4*GC + 2*AT. Ambiguous bases contribute
/// expected values unless `strict` is set. Port of Biopython `Tm_Wallace`.
pub fn tm_wallace_full(
  seq : String,
  options? : TmOptions = TmOptions::default(),
) -> Double raise SeqError {
  let s = if options.check {
    tm_check(seq, "Tm_Wallace")
  } else {
    seq.to_upper()
  }
  let count = fn(b : String) -> Int {
    let mut n = 0
    let mut i = 0
    while i < s.length() {
      if substring(s, i, 1) == b {
        n = n + 1
      }
      i = i + 1
    }
    n
  }
  let count_any = fn(bs : Array[String]) -> Int {
    let mut n = 0
    for b in bs {
      n = n + count(b)
    }
    n
  }
  let tm = 2.0 * count_any(["A", "T", "W"]).to_double() +
    4.0 * count_any(["C", "G", "S"]).to_double()
  let tmp = 3.0 * count_any(["K", "M", "N", "R", "Y"]).to_double() +
    10.0 / 3.0 * count_any(["B", "V"]).to_double() +
    8.0 / 3.0 * count_any(["D", "H"]).to_double()
  if options.strict && tmp > 0.0 {
    raise SeqError(
      "ambiguous bases B, D, H, K, M, N, R, V, Y not allowed when strict=True",
    )
  }
  tm + tmp
}

///|
/// Empirical GC-content Tm formulas: Tm = A + B(%GC) - C/N (+ salt
/// correction) - D(%mismatch). Port of Biopython `Tm_GC`; `valueset` picks
/// one of the eight published constant sets (default 7), and `userset`
/// supplies a custom (A, B, C, D) tuple.
pub fn tm_gc_full(
  seq : String,
  valueset? : Int = 7,
  userset? : (Double, Double, Double, Double)? = None,
  options? : TmOptions = TmOptions::default(),
) -> Double raise SeqError {
  let opts = options
  // Tm_GC uses the valueset-assigned salt correction; saltcorr 5 (the Tm_NN
  // entropy method) is invalid here. TmOptions defaults to 5 for Tm_NN, so
  // treat it as "unset" (0) in Tm_GC; valuesets override below.
  let mut saltcorr = if opts.saltcorr == 5 { 0 } else { opts.saltcorr }
  let s = if opts.check { tm_check(seq, "Tm_GC") } else { seq.to_upper() }
  if opts.strict {
    let ambiguous = ["K", "M", "N", "R", "Y", "B", "V", "D", "H"]
    for b in ambiguous {
      let mut i = 0
      while i < s.length() {
        if substring(s, i, 1) == b {
          raise SeqError(
            "ambiguous bases B, D, H, K, M, N, R, V, Y not allowed when 'strict=True'",
          )
        }
        i = i + 1
      }
    }
  }
  let mut percent_gc = gc_fraction(s, ambiguous="weighted") * 100.0
  let mut n_mismatch = 0
  let mut i = 0
  while i < s.length() {
    if substring(s, i, 1) == "X" {
      n_mismatch = n_mismatch + 1
    }
    i = i + 1
  }
  // gc_fraction counts X as 0.5
  percent_gc = percent_gc -
    n_mismatch.to_double() * 50.0 / s.length().to_double()
  let mut a = 0.0
  let mut b = 0.0
  let mut c = 0.0
  let mut d = 0.0
  match userset {
    Some(u) => {
      a = u.0
      b = u.1
      c = u.2
      d = u.3
    }
    None =>
      match valueset {
        1 => {
          a = 69.3
          b = 0.41
          c = 650.0
          d = 1.0
          saltcorr = 0
        }
        2 => {
          a = 81.5
          b = 0.41
          c = 675.0
          d = 1.0
          saltcorr = 0
        }
        3 => {
          a = 81.5
          b = 0.41
          c = 675.0
          d = 1.0
          saltcorr = 1
        }
        4 => {
          a = 81.5
          b = 0.41
          c = 500.0
          d = 1.0
          saltcorr = 2
        }
        5 => {
          a = 78.0
          b = 0.7
          c = 500.0
          d = 1.0
          saltcorr = 2
        }
        6 => {
          a = 67.0
          b = 0.8
          c = 500.0
          d = 1.0
          saltcorr = 2
        }
        7 => {
          a = 81.5
          b = 0.41
          c = 600.0
          d = 1.0
          saltcorr = 1
        }
        8 => {
          a = 77.1
          b = 0.41
          c = 528.0
          d = 1.0
          saltcorr = 4
        }
        _ => raise SeqError("allowed values for parameter 'valueset' are 0-8.")
      }
  }
  if valueset > 8 {
    raise SeqError("allowed values for parameter 'valueset' are 0-8.")
  }
  let mut tm = a + b * percent_gc - c / s.length().to_double()
  if saltcorr != 0 {
    tm = tm +
      tm_salt_correction(
        na=opts.na,
        k=opts.k,
        tris=opts.tris,
        mg=opts.mg,
        dntps=opts.dntps,
        salt_method=saltcorr,
        seq=s,
      )
  }
  tm = tm - d * (n_mismatch.to_double() * 100.0 / s.length().to_double())
  tm
}

///|
/// Tm via nearest-neighbor thermodynamics. Port of Biopython `Tm_NN`.
/// Supports the eight built-in NN tables (or custom tables via the
/// `*_custom` options), dangling ends (`shift` + `c_seq`), terminal and
/// internal mismatches (including inosine), self-complementary sequences,
/// strand concentrations (dnac1/dnac2) and salt corrections 1-7.
pub fn tm_nn_full(
  seq : String,
  options? : TmOptions = TmOptions::default(),
) -> Double raise SeqError {
  let opts = options
  let nn_table = tm_resolve(opts.nn_custom, opts.nn_table)
  let tmm_table = tm_resolve(opts.tmm_custom, opts.tmm_table)
  let imm_table = tm_resolve(opts.imm_custom, opts.imm_table)
  let de_table = tm_resolve(opts.de_custom, opts.de_table)

  let raw = seq.to_upper()
  let mut cseq = if opts.c_seq == "" {
    tm_complement(raw)
  } else {
    opts.c_seq.to_upper()
  }
  let mut tmp_seq = if opts.check { tm_check(raw, "Tm_NN") } else { raw }
  let mut tmp_cseq = if opts.check { tm_check(cseq, "Tm_NN") } else { cseq }
  cseq = tmp_cseq
  // checked original sequences: initiation terms and salt correction use these
  // (Biopython uses `seq`/`c_seq` here, not the dangling-end-trimmed copies)
  let seq_nn = tmp_seq
  let mut delta_h = 0.0
  let mut delta_s = 0.0

  // Dangling ends: alignment with shift and length differences
  if opts.shift != 0 || tmp_seq.length() != tmp_cseq.length() {
    if opts.shift > 0 {
      tmp_seq = mt_dots(opts.shift) + tmp_seq
    }
    if opts.shift < 0 {
      tmp_cseq = mt_dots(-opts.shift) + tmp_cseq
    }
    if tmp_cseq.length() > tmp_seq.length() {
      tmp_seq = tmp_seq + mt_dots(tmp_cseq.length() - tmp_seq.length())
    }
    if tmp_cseq.length() < tmp_seq.length() {
      tmp_cseq = tmp_cseq + mt_dots(tmp_seq.length() - tmp_cseq.length())
    }
    // remove over-dangling ends
    let mut go = true
    while go {
      if (tmp_seq.length() >= 2 && substring(tmp_seq, 0, 2) == "..") ||
        (tmp_cseq.length() >= 2 && substring(tmp_cseq, 0, 2) == "..") {
        tmp_seq = substring(tmp_seq, 1, tmp_seq.length() - 1)
        tmp_cseq = substring(tmp_cseq, 1, tmp_cseq.length() - 1)
      } else {
        go = false
      }
    }
    let mut go2 = true
    while go2 {
      let n = tmp_seq.length()
      if n >= 2 &&
        (
          substring(tmp_seq, n - 2, 2) == ".." ||
          substring(tmp_cseq, n - 2, 2) == ".."
        ) {
        tmp_seq = substring(tmp_seq, 0, n - 1)
        tmp_cseq = substring(tmp_cseq, 0, n - 1)
      } else {
        go2 = false
      }
    }
    // the dangling ends themselves
    if substring(tmp_seq, 0, 1) == "." || substring(tmp_cseq, 0, 1) == "." {
      let left_de = substring(tmp_seq, 0, 2) + "/" + substring(tmp_cseq, 0, 2)
      match de_table.get(left_de) {
        Some(v) => {
          delta_h = delta_h + v.0
          delta_s = delta_s + v.1
        }
        None => tm_key_error(left_de, opts.strict)
      }
      tmp_seq = substring(tmp_seq, 1, tmp_seq.length() - 1)
      tmp_cseq = substring(tmp_cseq, 1, tmp_cseq.length() - 1)
    }
    let n2 = tmp_seq.length()
    if substring(tmp_seq, n2 - 1, 1) == "." ||
      substring(tmp_cseq, n2 - 1, 1) == "." {
      let right_de = mt_reverse(substring(tmp_cseq, n2 - 2, 2)) +
        "/" +
        mt_reverse(substring(tmp_seq, n2 - 2, 2))
      match de_table.get(right_de) {
        Some(v) => {
          delta_h = delta_h + v.0
          delta_s = delta_s + v.1
        }
        None => tm_key_error(right_de, opts.strict)
      }
      tmp_seq = substring(tmp_seq, 0, n2 - 1)
      tmp_cseq = substring(tmp_cseq, 0, n2 - 1)
    }
  }

  // Terminal mismatches
  let left_tmm = mt_reverse(substring(tmp_cseq, 0, 2)) +
    "/" +
    mt_reverse(substring(tmp_seq, 0, 2))
  match tmm_table.get(left_tmm) {
    Some(v) => {
      delta_h = delta_h + v.0
      delta_s = delta_s + v.1
      tmp_seq = substring(tmp_seq, 1, tmp_seq.length() - 1)
      tmp_cseq = substring(tmp_cseq, 1, tmp_cseq.length() - 1)
    }
    None => ()
  }
  let nt = tmp_seq.length()
  let right_tmm = substring(tmp_seq, nt - 2, 2) +
    "/" +
    substring(tmp_cseq, nt - 2, 2)
  match tmm_table.get(right_tmm) {
    Some(v) => {
      delta_h = delta_h + v.0
      delta_s = delta_s + v.1
      tmp_seq = substring(tmp_seq, 0, nt - 1)
      tmp_cseq = substring(tmp_cseq, 0, nt - 1)
    }
    None => ()
  }

  // Initiation terms
  let init_add = fn(key : String) -> Unit raise SeqError {
    match nn_table.get(key) {
      Some(v) => {
        delta_h = delta_h + v.0
        delta_s = delta_s + v.1
      }
      None => tm_key_error(key, opts.strict)
    }
  }
  init_add("init")
  let gcf = gc_fraction(seq_nn, ambiguous="ignore")
  if gcf == 0.0 {
    init_add("init_allA/T")
  } else {
    init_add("init_oneG/C")
  }
  // Biopython adds the 5T/A penalty once for a 5' T and again for a 3' A
  if substring(seq_nn, 0, 1) == "T" {
    init_add("init_5T/A")
  }
  if substring(seq_nn, seq_nn.length() - 1, 1) == "A" {
    init_add("init_5T/A")
  }
  let ends = substring(seq_nn, 0, 1) + substring(seq_nn, seq_nn.length() - 1, 1)
  let mut at = 0
  let mut gc = 0
  for e in [substring(ends, 0, 1), substring(ends, 1, 1)] {
    if e == "A" || e == "T" {
      at = at + 1
    } else if e == "G" || e == "C" {
      gc = gc + 1
    }
  }
  let atv = nn_table.get("init_A/T")
  match atv {
    Some(v) => {
      delta_h = delta_h + v.0 * at.to_double()
      delta_s = delta_s + v.1 * at.to_double()
    }
    None => ()
  }
  let gcv = nn_table.get("init_G/C")
  match gcv {
    Some(v) => {
      delta_h = delta_h + v.0 * gc.to_double()
      delta_s = delta_s + v.1 * gc.to_double()
    }
    None => ()
  }

  // Zipping: internal NNs and internal mismatches
  let mut bn = 0
  while bn < tmp_seq.length() - 1 {
    let neighbors = substring(tmp_seq, bn, 2) + "/" + substring(tmp_cseq, bn, 2)
    let picked = match imm_table.get(neighbors) {
      Some(v) => Some(v)
      None =>
        match imm_table.get(mt_reverse(neighbors)) {
          Some(v) => Some(v)
          None =>
            match nn_table.get(neighbors) {
              Some(v) => Some(v)
              None =>
                match nn_table.get(mt_reverse(neighbors)) {
                  Some(v) => Some(v)
                  None => None
                }
            }
        }
    }
    match picked {
      Some(v) => {
        delta_h = delta_h + v.0
        delta_s = delta_s + v.1
      }
      None => tm_key_error(neighbors, opts.strict)
    }
    bn = bn + 1
  }

  // Strand concentration and symmetry
  let mut kk = (opts.dnac1 - opts.dnac2 / 2.0) * 1.0e-9
  if opts.selfcomp {
    kk = opts.dnac1 * 1.0e-9
    init_add("sym")
  }
  let r_gas = 1.987 // cal/(mol*K)
  let corr = if opts.saltcorr != 0 {
    tm_salt_correction(
      na=opts.na,
      k=opts.k,
      tris=opts.tris,
      mg=opts.mg,
      dntps=opts.dntps,
      salt_method=opts.saltcorr,
      seq=seq_nn,
    )
  } else {
    0.0
  }
  if opts.saltcorr == 5 {
    delta_s = delta_s + corr
  }
  let mut tm = 1000.0 * delta_h / (delta_s + r_gas * @math.ln(kk)) - 273.15
  if opts.saltcorr >= 1 && opts.saltcorr <= 4 {
    tm = tm + corr
  }
  if opts.saltcorr == 6 || opts.saltcorr == 7 {
    tm = 1.0 / (1.0 / (tm + 273.15) + corr) - 273.15
  }
  tm
}
