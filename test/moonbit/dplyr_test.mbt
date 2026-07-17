///|
test "dplyr_create_dataframe" {
  let df_cols = [
    @src.column("name", @src.data_type_string(), [
      @src.data_value_string("Alice"),
      @src.data_value_string("Bob"),
      @src.data_value_string("Charlie"),
    ]),
    @src.column("age", @src.data_type_int(), [
      @src.data_value_int(25),
      @src.data_value_int(30),
      @src.data_value_int(35),
    ]),
  ]
  let df = @src.dataframe(df_cols)
  assert_eq(@src.dataframe_nrow(df), 3)
  assert_eq(@src.dataframe_ncol(df), 2)
}

///|
test "dplyr_select" {
  let sel_cols = [
    @src.column("name", @src.data_type_string(), [
      @src.data_value_string("Alice"),
      @src.data_value_string("Bob"),
      @src.data_value_string("Charlie"),
    ]),
    @src.column("age", @src.data_type_int(), [
      @src.data_value_int(25),
      @src.data_value_int(30),
      @src.data_value_int(35),
    ]),
    @src.column("score", @src.data_type_double(), [
      @src.data_value_double(85.5),
      @src.data_value_double(90.0),
      @src.data_value_double(78.0),
    ]),
  ]
  let sel_df = @src.dataframe(sel_cols)
  let sel_result = @src.dataframe_select(sel_df, ["name", "score"])
  assert_eq(@src.dataframe_ncol(sel_result), 2)
}

///|
test "dplyr_filter" {
  let fil_cols = [
    @src.column("name", @src.data_type_string(), [
      @src.data_value_string("Alice"),
      @src.data_value_string("Bob"),
      @src.data_value_string("Charlie"),
    ]),
    @src.column("age", @src.data_type_int(), [
      @src.data_value_int(25),
      @src.data_value_int(30),
      @src.data_value_int(35),
    ]),
  ]
  let fil_df = @src.dataframe(fil_cols)
  let fil_result = @src.dataframe_filter(fil_df, fn(
    d : @src.DataFrame,
    i : Int,
  ) -> Bool {
    match d.columns[1].values[i] {
      @src.DataValue::IntValue(a) => a > 28
      _ => false
    }
  })
  assert_eq(@src.dataframe_nrow(fil_result), 2)
}

///|
test "dplyr_mutate" {
  let mut_cols = [
    @src.column("age", @src.data_type_int(), [
      @src.data_value_int(25),
      @src.data_value_int(30),
      @src.data_value_int(35),
    ]),
  ]
  let mut_df = @src.dataframe(mut_cols)
  let mut_result = @src.dataframe_mutate(mut_df, [
    (
      "age_plus_10",
      @src.data_type_int(),
      fn(d : @src.DataFrame, i : Int) -> @src.DataValue {
        match d.columns[0].values[i] {
          @src.DataValue::IntValue(a) => @src.data_value_int(a + 10)
          _ => @src.data_value_int(0)
        }
      },
    ),
  ])
  assert_eq(@src.dataframe_ncol(mut_result), 2)
}

///|
test "dplyr_arrange" {
  let arr_cols = [
    @src.column("name", @src.data_type_string(), [
      @src.data_value_string("Bob"),
      @src.data_value_string("Alice"),
      @src.data_value_string("Charlie"),
    ]),
    @src.column("age", @src.data_type_int(), [
      @src.data_value_int(30),
      @src.data_value_int(25),
      @src.data_value_int(35),
    ]),
  ]
  let arr_df = @src.dataframe(arr_cols)
  let arr_result = @src.dataframe_arrange(arr_df, "age", false)
  match arr_result.columns[1].values[0] {
    @src.DataValue::IntValue(a) => assert_eq(a, 25)
    _ => assert_true(false)
  }
}

///|
test "dplyr_group_by" {
  let grp_cols = [
    @src.column("group", @src.data_type_string(), [
      @src.data_value_string("A"),
      @src.data_value_string("A"),
      @src.data_value_string("B"),
    ]),
    @src.column("value", @src.data_type_int(), [
      @src.data_value_int(10),
      @src.data_value_int(20),
      @src.data_value_int(30),
    ]),
  ]
  let grp_df = @src.dataframe(grp_cols)
  let grp_result = @src.dataframe_group_by(grp_df, "group")
  assert_true(grp_result.length() >= 1)
}

///|
test "dplyr_summarize" {
  let smr_cols = [
    @src.column("group", @src.data_type_string(), [
      @src.data_value_string("A"),
      @src.data_value_string("A"),
      @src.data_value_string("B"),
      @src.data_value_string("B"),
    ]),
    @src.column("value", @src.data_type_int(), [
      @src.data_value_int(10),
      @src.data_value_int(20),
      @src.data_value_int(30),
      @src.data_value_int(40),
    ]),
  ]
  let smr_df = @src.dataframe(smr_cols)
  let smr_mean = fn(vals : Array[@src.DataValue]) -> @src.DataValue {
    let mut s = 0.0
    let mut c = 0
    let mut i = 0
    while i < vals.length() {
      match vals[i] {
        @src.DataValue::IntValue(v) => {
          s = s + v.to_double()
          c = c + 1
        }
        @src.DataValue::DoubleValue(v) => {
          s = s + v
          c = c + 1
        }
        _ => ()
      }
      i = i + 1
    }
    @src.data_value_double(if c > 0 { s / c.to_double() } else { 0.0 })
  }
  let smr_result = @src.dataframe_summarize(smr_df, "group", [
    ("mean_value", smr_mean),
  ])
  assert_eq(@src.dataframe_nrow(smr_result), 2)
}

///|
test "dplyr_inner_join" {
  let inj_cols1 = [
    @src.column("id", @src.data_type_int(), [
      @src.data_value_int(1),
      @src.data_value_int(2),
      @src.data_value_int(3),
    ]),
    @src.column("name", @src.data_type_string(), [
      @src.data_value_string("Alice"),
      @src.data_value_string("Bob"),
      @src.data_value_string("Charlie"),
    ]),
  ]
  let inj_df1 = @src.dataframe(inj_cols1)
  let inj_cols2 = [
    @src.column("id", @src.data_type_int(), [
      @src.data_value_int(1),
      @src.data_value_int(2),
      @src.data_value_int(4),
    ]),
    @src.column("score", @src.data_type_double(), [
      @src.data_value_double(85.5),
      @src.data_value_double(90.0),
      @src.data_value_double(78.0),
    ]),
  ]
  let inj_df2 = @src.dataframe(inj_cols2)
  let inj_result = @src.dataframe_inner_join(inj_df1, inj_df2, "id")
  assert_eq(@src.dataframe_nrow(inj_result), 2)
}

///|
test "dplyr_left_join" {
  let lfj_cols1 = [
    @src.column("id", @src.data_type_int(), [
      @src.data_value_int(1),
      @src.data_value_int(2),
      @src.data_value_int(3),
    ]),
    @src.column("name", @src.data_type_string(), [
      @src.data_value_string("Alice"),
      @src.data_value_string("Bob"),
      @src.data_value_string("Charlie"),
    ]),
  ]
  let lfj_df1 = @src.dataframe(lfj_cols1)
  let lfj_cols2 = [
    @src.column("id", @src.data_type_int(), [
      @src.data_value_int(1),
      @src.data_value_int(4),
    ]),
    @src.column("score", @src.data_type_double(), [
      @src.data_value_double(85.5),
      @src.data_value_double(78.0),
    ]),
  ]
  let lfj_df2 = @src.dataframe(lfj_cols2)
  let lfj_result = @src.dataframe_left_join(lfj_df1, lfj_df2, "id")
  assert_eq(@src.dataframe_nrow(lfj_result), 3)
}
