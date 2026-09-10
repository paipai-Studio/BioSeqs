///|
/// ProtParam amino acid scale tables ported from Biopython 1.88
/// Bio.SeqUtils.ProtParamData (28 hydrophobicity scales + Vihinen Flex).

///|
/// Normalized flexibility parameters (Vihinen et al. 1994).
let flex_values : Map[Char, Double] = Map(
  [
    ('A', 0.984), ('R', 1.008), ('N', 1.048), ('D', 1.068), ('C', 0.906),
    ('Q', 1.037), ('E', 1.094), ('G', 1.031), ('H', 0.95), ('I', 0.927),
    ('L', 0.935), ('K', 1.102), ('M', 0.952), ('F', 0.915), ('P', 1.049),
    ('S', 1.046), ('T', 0.997), ('W', 0.904), ('Y', 0.929), ('V', 0.931),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "KyteDoolitle" (ProtParamData).
let gravy_scale_0 : Map[Char, Double] = Map(
  [
    ('A', 1.8), ('R', -4.5), ('N', -3.5), ('D', -3.5), ('C', 2.5),
    ('Q', -3.5), ('E', -3.5), ('G', -0.4), ('H', -3.2), ('I', 4.5),
    ('L', 3.8), ('K', -3.9), ('M', 1.9), ('F', 2.8), ('P', -1.6),
    ('S', -0.8), ('T', -0.7), ('W', -0.9), ('Y', -1.3), ('V', 4.2),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Aboderin" (ProtParamData).
let gravy_scale_1 : Map[Char, Double] = Map(
  [
    ('A', 5.1), ('R', 2.0), ('N', 0.6), ('D', 0.7), ('C', 0.0),
    ('Q', 1.4), ('E', 1.8), ('G', 4.1), ('H', 1.6), ('I', 9.3),
    ('L', 10.0), ('K', 1.3), ('M', 8.7), ('F', 9.6), ('P', 4.9),
    ('S', 3.1), ('T', 3.5), ('W', 9.2), ('Y', 8.0), ('V', 8.5),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "AbrahamLeo" (ProtParamData).
let gravy_scale_2 : Map[Char, Double] = Map(
  [
    ('A', 0.44), ('R', -2.42), ('N', -1.32), ('D', -0.31), ('C', 0.58),
    ('Q', -0.71), ('E', -0.34), ('G', 0.0), ('H', -0.01), ('I', 2.46),
    ('L', 2.46), ('K', -2.45), ('M', 1.1), ('F', 2.54), ('P', 1.29),
    ('S', -0.84), ('T', -0.41), ('W', 2.56), ('Y', 1.63), ('V', 1.73),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Argos" (ProtParamData).
let gravy_scale_3 : Map[Char, Double] = Map(
  [
    ('A', 0.61), ('R', 0.6), ('N', 0.06), ('D', 0.46), ('C', 1.07),
    ('Q', 0.0), ('E', 0.47), ('G', 0.07), ('H', 0.61), ('I', 2.22),
    ('L', 1.53), ('K', 1.15), ('M', 1.18), ('F', 2.02), ('P', 1.95),
    ('S', 0.05), ('T', 0.05), ('W', 2.65), ('Y', 1.88), ('V', 1.32),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "BlackMould" (ProtParamData).
let gravy_scale_4 : Map[Char, Double] = Map(
  [
    ('A', 0.616), ('R', 0.0), ('N', 0.236), ('D', 0.028), ('C', 0.68),
    ('Q', 0.251), ('E', 0.043), ('G', 0.501), ('H', 0.165), ('I', 0.943),
    ('L', 0.943), ('K', 0.283), ('M', 0.738), ('F', 1.0), ('P', 0.711),
    ('S', 0.359), ('T', 0.45), ('W', 0.878), ('Y', 0.88), ('V', 0.825),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "BullBreese" (ProtParamData).
let gravy_scale_5 : Map[Char, Double] = Map(
  [
    ('A', 0.61), ('R', 0.69), ('N', 0.89), ('D', 0.61), ('C', 0.36),
    ('Q', 0.97), ('E', 0.51), ('G', 0.81), ('H', 0.69), ('I', -1.45),
    ('L', -1.65), ('K', 0.46), ('M', -0.66), ('F', -1.52), ('P', -0.17),
    ('S', 0.42), ('T', 0.29), ('W', -1.2), ('Y', -1.43), ('V', -0.75),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Casari" (ProtParamData).
let gravy_scale_6 : Map[Char, Double] = Map(
  [
    ('A', 0.2), ('R', -0.7), ('N', -0.5), ('D', -1.4), ('C', 1.9),
    ('Q', -1.1), ('E', -1.3), ('G', -0.1), ('H', 0.4), ('I', 1.4),
    ('L', 0.5), ('K', -1.6), ('M', 0.5), ('F', 1.0), ('P', -1.0),
    ('S', -0.7), ('T', -0.4), ('W', 1.6), ('Y', 0.5), ('V', 0.7),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Cid" (ProtParamData).
let gravy_scale_7 : Map[Char, Double] = Map(
  [
    ('A', 0.02), ('R', -0.42), ('N', -0.77), ('D', -1.04), ('C', 0.77),
    ('Q', -1.1), ('E', -1.14), ('G', -0.8), ('H', 0.26), ('I', 1.81),
    ('L', 1.14), ('K', -0.41), ('M', 1.0), ('F', 1.35), ('P', -0.09),
    ('S', -0.97), ('T', -0.77), ('W', 1.71), ('Y', 1.11), ('V', 1.13),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Cowan3.4" (ProtParamData).
let gravy_scale_8 : Map[Char, Double] = Map(
  [
    ('A', 0.42), ('R', -1.56), ('N', -1.03), ('D', -0.51), ('C', 0.84),
    ('Q', -0.96), ('E', -0.37), ('G', 0.0), ('H', -2.28), ('I', 1.81),
    ('L', 1.8), ('K', -2.03), ('M', 1.18), ('F', 1.74), ('P', 0.86),
    ('S', -0.64), ('T', -0.26), ('W', 1.46), ('Y', 0.51), ('V', 1.34),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Cowan7.5" (ProtParamData).
let gravy_scale_9 : Map[Char, Double] = Map(
  [
    ('A', 0.35), ('R', -1.5), ('N', -0.99), ('D', -2.15), ('C', 0.76),
    ('Q', -0.93), ('E', -1.95), ('G', 0.0), ('H', -0.65), ('I', 1.83),
    ('L', 1.8), ('K', -1.54), ('M', 1.1), ('F', 1.69), ('P', 0.84),
    ('S', -0.63), ('T', -0.27), ('W', 1.35), ('Y', 0.39), ('V', 1.32),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Eisenberg" (ProtParamData).
let gravy_scale_10 : Map[Char, Double] = Map(
  [
    ('A', 0.62), ('R', -2.53), ('N', -0.78), ('D', -0.9), ('C', 0.29),
    ('Q', -0.85), ('E', -0.74), ('G', 0.48), ('H', -0.4), ('I', 1.38),
    ('L', 1.06), ('K', -1.5), ('M', 0.64), ('F', 1.19), ('P', 0.12),
    ('S', -0.18), ('T', -0.05), ('W', 0.81), ('Y', 0.26), ('V', 1.08),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Engelman" (ProtParamData).
let gravy_scale_11 : Map[Char, Double] = Map(
  [
    ('A', -1.6), ('R', 12.3), ('N', 4.8), ('D', 9.2), ('C', -2.0),
    ('Q', 4.1), ('E', 8.2), ('G', -1.0), ('H', 3.0), ('I', -3.1),
    ('L', -2.8), ('K', 8.8), ('M', -3.4), ('F', -3.7), ('P', 0.2),
    ('S', -0.6), ('T', -1.2), ('W', -1.9), ('Y', 0.7), ('V', -2.6),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Fasman" (ProtParamData).
let gravy_scale_12 : Map[Char, Double] = Map(
  [
    ('A', -0.21), ('R', 2.11), ('N', 0.96), ('D', 1.36), ('C', -6.04),
    ('Q', 1.52), ('E', 2.3), ('G', 0.0), ('H', -1.23), ('I', -4.81),
    ('L', -4.68), ('K', 3.88), ('M', -3.66), ('F', -4.65), ('P', 0.75),
    ('S', 1.74), ('T', 0.78), ('W', -3.32), ('Y', -1.01), ('V', -3.5),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Fauchere" (ProtParamData).
let gravy_scale_13 : Map[Char, Double] = Map(
  [
    ('A', 0.31), ('R', -1.01), ('N', -0.6), ('D', -0.77), ('C', 1.54),
    ('Q', -0.22), ('E', -0.64), ('G', 0.0), ('H', 0.13), ('I', 1.8),
    ('L', 1.7), ('K', -0.99), ('M', 1.23), ('F', 1.79), ('P', 0.72),
    ('S', -0.04), ('T', 0.26), ('W', 2.25), ('Y', 0.96), ('V', 1.22),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "GoldSack" (ProtParamData).
let gravy_scale_14 : Map[Char, Double] = Map(
  [
    ('A', 0.75), ('R', 0.75), ('N', 0.69), ('D', 0.0), ('C', 1.0),
    ('Q', 0.59), ('E', 0.0), ('G', 0.0), ('H', 0.0), ('I', 2.95),
    ('L', 2.4), ('K', 1.5), ('M', 1.3), ('F', 2.65), ('P', 2.6),
    ('S', 0.0), ('T', 0.45), ('W', 3.0), ('Y', 2.85), ('V', 1.7),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Guy" (ProtParamData).
let gravy_scale_15 : Map[Char, Double] = Map(
  [
    ('A', 0.1), ('R', 1.91), ('N', 0.48), ('D', 0.78), ('C', -1.42),
    ('Q', 0.95), ('E', 0.83), ('G', 0.33), ('H', -0.5), ('I', -1.13),
    ('L', -1.18), ('K', 1.4), ('M', -1.59), ('F', -2.12), ('P', 0.73),
    ('S', 0.52), ('T', 0.07), ('W', -0.51), ('Y', -0.21), ('V', -1.27),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Jones" (ProtParamData).
let gravy_scale_16 : Map[Char, Double] = Map(
  [
    ('A', 0.87), ('R', 0.85), ('N', 0.09), ('D', 0.66), ('C', 1.52),
    ('Q', 0.0), ('E', 0.67), ('G', 0.1), ('H', 0.87), ('I', 3.15),
    ('L', 2.17), ('K', 1.64), ('M', 1.67), ('F', 2.87), ('P', 2.77),
    ('S', 0.07), ('T', 0.07), ('W', 3.77), ('Y', 2.67), ('V', 1.87),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Juretic" (ProtParamData).
let gravy_scale_17 : Map[Char, Double] = Map(
  [
    ('A', 1.1), ('R', -5.1), ('N', -3.5), ('D', -3.6), ('C', 2.5),
    ('Q', -3.68), ('E', -3.2), ('G', -0.64), ('H', -3.2), ('I', 4.5),
    ('L', 3.8), ('K', -4.11), ('M', 1.9), ('F', 2.8), ('P', -1.9),
    ('S', -0.5), ('T', -0.7), ('W', -0.46), ('Y', -1.3), ('V', 4.2),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Kidera" (ProtParamData).
let gravy_scale_18 : Map[Char, Double] = Map(
  [
    ('A', -0.27), ('R', 1.87), ('N', 0.81), ('D', 0.81), ('C', -1.05),
    ('Q', 1.1), ('E', 1.17), ('G', -0.16), ('H', 0.28), ('I', -0.77),
    ('L', -1.1), ('K', 1.7), ('M', -0.73), ('F', -1.43), ('P', -0.75),
    ('S', 0.42), ('T', 0.63), ('W', -1.57), ('Y', -0.56), ('V', -0.4),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Miyazawa" (ProtParamData).
let gravy_scale_19 : Map[Char, Double] = Map(
  [
    ('A', 5.33), ('R', 4.18), ('N', 3.71), ('D', 3.59), ('C', 7.93),
    ('Q', 3.87), ('E', 3.65), ('G', 4.48), ('H', 5.1), ('I', 8.83),
    ('L', 8.47), ('K', 2.95), ('M', 8.95), ('F', 9.03), ('P', 3.87),
    ('S', 4.09), ('T', 4.49), ('W', 7.66), ('Y', 5.89), ('V', 7.63),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Parker" (ProtParamData).
let gravy_scale_20 : Map[Char, Double] = Map(
  [
    ('A', 2.1), ('R', 4.2), ('N', 7.0), ('D', 10.0), ('C', 1.4),
    ('Q', 6.0), ('E', 7.8), ('G', 5.7), ('H', 2.1), ('I', -8.0),
    ('L', -9.2), ('K', 5.7), ('M', -4.2), ('F', -9.2), ('P', 2.1),
    ('S', 6.5), ('T', 5.2), ('W', -10.0), ('Y', -1.9), ('V', -3.7),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Ponnuswamy" (ProtParamData).
let gravy_scale_21 : Map[Char, Double] = Map(
  [
    ('A', 0.85), ('R', 0.2), ('N', -0.48), ('D', -1.1), ('C', 2.1),
    ('Q', -0.42), ('E', -0.79), ('G', 0.0), ('H', 0.22), ('I', 3.14),
    ('L', 1.99), ('K', -1.19), ('M', 1.42), ('F', 1.69), ('P', -1.14),
    ('S', -0.52), ('T', -0.08), ('W', 1.76), ('Y', 1.37), ('V', 2.53),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Rose" (ProtParamData).
let gravy_scale_22 : Map[Char, Double] = Map(
  [
    ('A', 0.74), ('R', 0.64), ('N', 0.63), ('D', 0.62), ('C', 0.91),
    ('Q', 0.62), ('E', 0.62), ('G', 0.72), ('H', 0.78), ('I', 0.88),
    ('L', 0.85), ('K', 0.52), ('M', 0.85), ('F', 0.88), ('P', 0.64),
    ('S', 0.66), ('T', 0.7), ('W', 0.85), ('Y', 0.76), ('V', 0.86),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Roseman" (ProtParamData).
let gravy_scale_23 : Map[Char, Double] = Map(
  [
    ('A', 0.39), ('R', -3.95), ('N', -1.91), ('D', -3.81), ('C', 0.25),
    ('Q', -1.3), ('E', -2.91), ('G', 0.0), ('H', -0.64), ('I', 1.82),
    ('L', 1.82), ('K', -2.77), ('M', 0.96), ('F', 2.27), ('P', 0.99),
    ('S', -1.24), ('T', -1.0), ('W', 2.13), ('Y', 1.47), ('V', 1.3),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Sweet" (ProtParamData).
let gravy_scale_24 : Map[Char, Double] = Map(
  [
    ('A', -0.4), ('R', -0.59), ('N', -0.92), ('D', -1.31), ('C', 0.17),
    ('Q', -0.91), ('E', -1.22), ('G', -0.67), ('H', -0.64), ('I', 1.25),
    ('L', 1.22), ('K', -0.67), ('M', 1.02), ('F', 1.92), ('P', -0.49),
    ('S', -0.55), ('T', -0.28), ('W', 0.5), ('Y', 1.67), ('V', 0.91),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Tanford" (ProtParamData).
let gravy_scale_25 : Map[Char, Double] = Map(
  [
    ('A', 0.62), ('R', -2.53), ('N', -0.78), ('D', -0.09), ('C', 0.29),
    ('Q', -0.85), ('E', -0.74), ('G', 0.48), ('H', -0.4), ('I', 1.38),
    ('L', 1.53), ('K', -1.5), ('M', 0.64), ('F', 1.19), ('P', 0.12),
    ('S', -0.18), ('T', -0.05), ('W', 0.81), ('Y', 0.26), ('V', 1.8),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Wilson" (ProtParamData).
let gravy_scale_26 : Map[Char, Double] = Map(
  [
    ('A', -0.3), ('R', -1.1), ('N', -0.2), ('D', -1.4), ('C', 6.3),
    ('Q', -0.2), ('E', 0.0), ('G', 1.2), ('H', -1.3), ('I', 4.3),
    ('L', 6.6), ('K', -3.6), ('M', 2.5), ('F', 7.5), ('P', 2.2),
    ('S', -0.6), ('T', -2.2), ('W', 7.9), ('Y', 7.1), ('V', 5.9),
  ],
  capacity=20,
)

///|
/// Hydrophobicity scale "Zimmerman" (ProtParamData).
let gravy_scale_27 : Map[Char, Double] = Map(
  [
    ('A', 0.83), ('R', 0.83), ('N', 0.09), ('D', 0.64), ('C', 1.48),
    ('Q', 0.0), ('E', 0.65), ('G', 0.1), ('H', 1.1), ('I', 3.07),
    ('L', 2.52), ('K', 1.6), ('M', 1.4), ('F', 2.75), ('P', 2.7),
    ('S', 0.14), ('T', 0.54), ('W', 0.31), ('Y', 2.97), ('V', 1.79),
  ],
  capacity=20,
)

///|
/// Return the registered hydrophobicity scale table with the given name.
/// Raises SeqError for unknown scale names (mirrors Biopython ValueError).
pub fn gravy_scale_table(name : String) -> Map[Char, Double] raise SeqError {
  match name {
    "KyteDoolitle" => gravy_scale_0
    "Aboderin" => gravy_scale_1
    "AbrahamLeo" => gravy_scale_2
    "Argos" => gravy_scale_3
    "BlackMould" => gravy_scale_4
    "BullBreese" => gravy_scale_5
    "Casari" => gravy_scale_6
    "Cid" => gravy_scale_7
    "Cowan3.4" => gravy_scale_8
    "Cowan7.5" => gravy_scale_9
    "Eisenberg" => gravy_scale_10
    "Engelman" => gravy_scale_11
    "Fasman" => gravy_scale_12
    "Fauchere" => gravy_scale_13
    "GoldSack" => gravy_scale_14
    "Guy" => gravy_scale_15
    "Jones" => gravy_scale_16
    "Juretic" => gravy_scale_17
    "Kidera" => gravy_scale_18
    "Miyazawa" => gravy_scale_19
    "Parker" => gravy_scale_20
    "Ponnuswamy" => gravy_scale_21
    "Rose" => gravy_scale_22
    "Roseman" => gravy_scale_23
    "Sweet" => gravy_scale_24
    "Tanford" => gravy_scale_25
    "Wilson" => gravy_scale_26
    "Zimmerman" => gravy_scale_27
    _ => raise SeqError("scale: " + name + " not known")
  }
}

///|
/// Hopp & Wood (1981) hydrophilicity scale (ProtParamData.hw).
let hw_values : Map[Char, Double] = Map(
  [
    ('A', -0.5),
    ('R', 3.0),
    ('N', 0.2),
    ('D', 3.0),
    ('C', -1.0),
    ('Q', 0.2),
    ('E', 3.0),
    ('G', 0.0),
    ('H', -0.5),
    ('I', -1.8),
    ('L', -1.8),
    ('K', 3.0),
    ('M', -1.3),
    ('F', -2.5),
    ('P', 0.0),
    ('S', 0.3),
    ('T', -0.4),
    ('W', -3.4),
    ('Y', -2.3),
    ('V', -1.5),
  ],
  capacity=20,
)

///|
/// Return the Hopp & Wood hydrophilicity table, ready for protein_scale.
pub fn hydrophilicity_table() -> Map[Char, Double] {
  hw_values
}
