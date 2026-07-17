///|
fn main {
  println("=== dplyr Data Manipulation Demo ===")
  println("")

  let cols = [
    @bio.column("gene", @bio.data_type_string(), [
      @bio.data_value_string("GeneA"),
      @bio.data_value_string("GeneB"),
      @bio.data_value_string("GeneC"),
      @bio.data_value_string("GeneD"),
      @bio.data_value_string("GeneE"),
      @bio.data_value_string("GeneF"),
    ]),
    @bio.column("condition", @bio.data_type_string(), [
      @bio.data_value_string("Control"),
      @bio.data_value_string("Control"),
      @bio.data_value_string("Treatment"),
      @bio.data_value_string("Treatment"),
      @bio.data_value_string("Control"),
      @bio.data_value_string("Treatment"),
    ]),
    @bio.column("expression", @bio.data_type_int(), [
      @bio.data_value_int(100),
      @bio.data_value_int(120),
      @bio.data_value_int(200),
      @bio.data_value_int(220),
      @bio.data_value_int(50),
      @bio.data_value_int(150),
    ]),
    @bio.column("pvalue", @bio.data_type_double(), [
      @bio.data_value_double(0.05),
      @bio.data_value_double(0.01),
      @bio.data_value_double(0.001),
      @bio.data_value_double(0.0001),
      @bio.data_value_double(0.1),
      @bio.data_value_double(0.02),
    ]),
  ]

  let df = @bio.dataframe(cols)

  println("1. Original DataFrame:")
  println("   Rows: \{@bio.dataframe_nrow(df)}")
  println("   Cols: \{@bio.dataframe_ncol(df)}")
  println("   Column names: \{@bio.dataframe_colnames(df)}")
  println("")

  println("2. Select columns (gene, expression):")
  let selected = @bio.dataframe_select(df, ["gene", "expression"])
  println("   Column names: \{@bio.dataframe_colnames(selected)}")
  println("")

  println("3. Filter rows (pvalue < 0.05):")
  let filtered = @bio.dataframe_filter(df, fn(
    df : @bio.DataFrame,
    i : Int,
  ) -> Bool {
    match df.columns[3].values[i] {
      @bio.DataValue::DoubleValue(p) => p < 0.05
      _ => false
    }
  })
  println("   Rows after filter: \{@bio.dataframe_nrow(filtered)}")
  println("")

  println("4. Mutate (add doubled_expression column):")
  let mutated = @bio.dataframe_mutate(df, [
    (
      "doubled_expression",
      @bio.data_type_int(),
      fn(df : @bio.DataFrame, i : Int) -> @bio.DataValue {
        match df.columns[2].values[i] {
          @bio.DataValue::IntValue(expr) => @bio.data_value_int(expr * 2)
          _ => @bio.data_value_int(0)
        }
      },
    ),
  ])
  println("   Column names: \{@bio.dataframe_colnames(mutated)}")
  println("")

  println("5. Arrange by expression (descending):")
  let arranged = @bio.dataframe_arrange(df, "expression", true)
  println("   First row expression (highest): ")
  match arranged.columns[2].values[0] {
    @bio.DataValue::IntValue(expr) => println("     \{expr}")
    _ => println("     N/A")
  }
  println("")

  println("6. Group by condition:")
  let groups = @bio.dataframe_group_by(df, "condition")
  println("   Number of groups: \{groups.length()}")
  let mut i = 0
  while i < groups.length() {
    let (key, group_df) = groups[i]
    println("   Group \{key}: \{@bio.dataframe_nrow(group_df)} rows")
    i = i + 1
  }
  println("")

  println("7. Summarize (mean expression by condition):")
  let mean_func = fn(values : Array[@bio.DataValue]) -> @bio.DataValue {
    let mut sum = 0.0
    let mut count = 0
    let mut j = 0
    while j < values.length() {
      match values[j] {
        @bio.DataValue::IntValue(v) => {
          sum = sum + v.to_double()
          count = count + 1
        }
        @bio.DataValue::DoubleValue(v) => {
          sum = sum + v
          count = count + 1
        }
        _ => ()
      }
      j = j + 1
    }
    @bio.data_value_double(
      if count > 0 {
        sum / count.to_double()
      } else {
        0.0
      },
    )
  }
  let summarized = @bio.dataframe_summarize(df, "condition", [
    ("mean_expression", mean_func),
  ])
  println("   Column names: \{@bio.dataframe_colnames(summarized)}")
  println("   Rows: \{@bio.dataframe_nrow(summarized)}")
  println("")

  println("8. Create second DataFrame for join:")
  let cols2 = [
    @bio.column("gene", @bio.data_type_string(), [
      @bio.data_value_string("GeneA"),
      @bio.data_value_string("GeneB"),
      @bio.data_value_string("GeneG"),
    ]),
    @bio.column("function", @bio.data_type_string(), [
      @bio.data_value_string("Kinase"),
      @bio.data_value_string("Transcription Factor"),
      @bio.data_value_string("Unknown"),
    ]),
  ]
  let df2 = @bio.dataframe(cols2)

  println("9. Inner join (by gene):")
  let inner_joined = @bio.dataframe_inner_join(df, df2, "gene")
  println("   Rows after inner join: \{@bio.dataframe_nrow(inner_joined)}")
  println("")

  println("10. Left join (by gene):")
  let left_joined = @bio.dataframe_left_join(df, df2, "gene")
  println("   Rows after left join: \{@bio.dataframe_nrow(left_joined)}")
  println("")

  println("=== Analysis Complete ===")
}
