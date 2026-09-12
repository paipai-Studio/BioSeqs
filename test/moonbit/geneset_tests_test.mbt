///|
/// Tests for the limma gene-set test port (rankSumTestWithCorrelation,
/// cameraPR.default, camera.default). All expected values were frozen from
/// limma 3.68.5 under R 4.3.3 using /tmp/gen_geneset_expected.R (the R
/// sources were sourced directly; no limma package install required).

///|
fn gstt_approx(
  left : Double,
  right : Double,
  tolerance? : Double = 1.0e-9,
) -> Bool {
  (left - right).abs() <= tolerance
}

///|
fn gstt_raises(operation : () -> Unit raise @src.GeneSetTestError) -> Bool {
  try {
    operation()
    false
  } catch {
    _ => true
  }
}

///|
/// 40 pre-ranked statistics, set.seed(77); round(rnorm(40), 4) in R.
fn gstt_statistic_s() -> Array[Double] {
  [
    -0.54959999999999998, 1.0911, 0.63980000000000004, 1.0426, 0.16969999999999999,
    1.1377999999999999, -0.97060000000000002, -0.1318, 0.1462, 1.4413, -2.9413999999999998,
    -0.24279999999999999, -0.1406, -0.0327, 0.27979999999999999, 0.59009999999999996,
    1.0243, 2.1073, 0.15459999999999999, 0.91310000000000002, -0.25419999999999998,
    1.5194000000000001, 1.7813000000000001, -0.87880000000000003, -1.5291999999999999,
    0.13589999999999999, -0.70889999999999997, -1.4096, 1.8309, 1.2902, -2.3624000000000001,
    -0.55100000000000005, -0.30470000000000003, -0.75, 0.14369999999999999, -0.5494,
    0.15989999999999999, -0.0877, 0.081100000000000005, 0.89910000000000001,
  ]
}

///|
/// R index c(1,3,5,7,9,11), converted to zero-based.
fn gstt_index_six() -> Array[Int] {
  [0, 2, 4, 6, 8, 10]
}

///|
test "rankSumTestWithCorrelation matches limma (normal reference)" {
  let s = gstt_statistic_s()
  let idx = gstt_index_six()
  let (down0, up0) = @src.rank_sum_test_with_correlation(idx, s)
  assert_true(gstt_approx(down0, 0.11640579004813852))
  assert_true(gstt_approx(up0, 0.89084329546327112))
  let (downc, upc) = @src.rank_sum_test_with_correlation(
    idx,
    s,
    correlation=0.25,
  )
  assert_true(gstt_approx(downc, 0.19933876523395497))
  assert_true(gstt_approx(upc, 0.80806248747259279))
}

///|
test "rankSumTestWithCorrelation finite df uses Student t" {
  let s = gstt_statistic_s()
  let idx = gstt_index_six()
  let (down, up) = @src.rank_sum_test_with_correlation(
    idx,
    s,
    correlation=0.25,
    df=Some(38.0),
  )
  assert_true(gstt_approx(down, 0.20197976385656841))
  assert_true(gstt_approx(up, 0.80532847404424568))
}

///|
test "rankSumTestWithCorrelation average-rank tie correction" {
  // R: statistic <- rep(1:8, 5); index <- 1:10
  let stats : Array[Double] = []
  let mut rep = 0
  while rep < 5 {
    let mut v = 1
    while v <= 8 {
      stats.push(v.to_double())
      v += 1
    }
    rep += 1
  }
  let idx : Array[Int] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
  let (down0, up0) = @src.rank_sum_test_with_correlation(idx, stats)
  assert_true(gstt_approx(down0, 0.17659512415857861))
  assert_true(gstt_approx(up0, 0.83144492834555495))
  let (downc, upc) = @src.rank_sum_test_with_correlation(
    idx,
    stats,
    correlation=0.3,
    df=Some(20.0),
  )
  assert_true(gstt_approx(downc, 0.29607158366417191))
  assert_true(gstt_approx(upc, 0.71012646798155066))
}

///|
/// Gene sets (R 1-based positions converted to zero-based):
/// UP c(2,5,9,14,20,27,31,38), DOWN c(3,11,17,23,33),
/// NULLSET c(1,4,8,15,22), SINGLE 20.
fn gstt_sets() -> Array[@src.CameraGeneSet] raise @src.GeneSetTestError {
  [
    @src.CameraGeneSet::create("UP", [1, 4, 8, 13, 19, 26, 30, 37]),
    @src.CameraGeneSet::create("DOWN", [2, 10, 16, 22, 32]),
    @src.CameraGeneSet::create("NULLSET", [0, 3, 7, 14, 21]),
    @src.CameraGeneSet::create("SINGLE", [19]),
  ]
}

///|
fn gstt_check_row(
  row : @src.CameraResult,
  name~ : String,
  n_genes~ : Int,
  direction~ : String,
  p_value~ : Double,
  fdr~ : Double,
  tolerance? : Double = 1.0e-9,
) -> Unit raise Error {
  assert_eq(row.name(), name)
  assert_eq(row.n_genes(), n_genes)
  assert_eq(row.direction(), direction)
  assert_true(gstt_approx(row.p_value(), p_value, tolerance~))
  assert_true(gstt_approx(row.fdr().unwrap_or(-1.0), fdr, tolerance~))
}

///|
test "cameraPR default parametric test matches limma" {
  let results = @src.camera_pr(gstt_statistic_s(), gstt_sets(), sort=false)
  assert_eq(results.length(), 4)
  assert_true(results[0].correlation().is_none())
  gstt_check_row(
    results[0],
    name="UP",
    n_genes=8,
    direction="Down",
    p_value=0.55333238653511674,
    fdr=0.73777651538015565,
  )
  gstt_check_row(
    results[1],
    name="DOWN",
    n_genes=5,
    direction="Down",
    p_value=0.89113368413850469,
    fdr=0.89113368413850469,
  )
  gstt_check_row(
    results[2],
    name="NULLSET",
    n_genes=5,
    direction="Up",
    p_value=0.48724286541348771,
    fdr=0.73777651538015565,
  )
  gstt_check_row(
    results[3],
    name="SINGLE",
    n_genes=1,
    direction="Up",
    p_value=0.4596798481474853,
    fdr=0.73777651538015565,
  )
}

///|
test "cameraPR ranks match limma" {
  let results = @src.camera_pr(
    gstt_statistic_s(),
    gstt_sets(),
    use_ranks=true,
    sort=false,
  )
  gstt_check_row(
    results[0],
    name="UP",
    n_genes=8,
    direction="Down",
    p_value=0.77935700783756134,
    fdr=0.77935700783756134,
  )
  gstt_check_row(
    results[1],
    name="DOWN",
    n_genes=5,
    direction="Up",
    p_value=0.77817726349346572,
    fdr=0.77935700783756134,
  )
  gstt_check_row(
    results[2],
    name="NULLSET",
    n_genes=5,
    direction="Up",
    p_value=0.46884536925434778,
    fdr=0.77935700783756134,
  )
  gstt_check_row(
    results[3],
    name="SINGLE",
    n_genes=1,
    direction="Up",
    p_value=0.43558717431803001,
    fdr=0.77935700783756134,
  )
}

///|
test "cameraPR non-directional uses absolute statistics and ranks" {
  let results = @src.camera_pr(
    gstt_statistic_s(),
    gstt_sets(),
    use_ranks=true,
    directional=false,
    sort=false,
  )
  assert_eq(results[0].direction(), "")
  gstt_check_row(
    results[0],
    name="UP",
    n_genes=8,
    direction="",
    p_value=0.81763539867562729,
    fdr=0.81763539867562729,
  )
  gstt_check_row(
    results[1],
    name="DOWN",
    n_genes=5,
    direction="",
    p_value=0.07949511151233471,
    fdr=0.31798044604933884,
  )
  gstt_check_row(
    results[2],
    name="NULLSET",
    n_genes=5,
    direction="",
    p_value=0.59539797913521586,
    fdr=0.79386397218028781,
  )
  gstt_check_row(
    results[3],
    name="SINGLE",
    n_genes=1,
    direction="",
    p_value=0.36447655469089063,
    fdr=0.72895310938178126,
  )
}

///|
test "cameraPR per-set correlations add Correlation column" {
  let results = @src.camera_pr(
    gstt_statistic_s(),
    gstt_sets(),
    per_set_inter_gene_cor=Some([0.0, 0.02, 0.05, 0.01]),
    sort=false,
  )
  assert_true(gstt_approx(results[0].correlation().unwrap(), 0.0))
  assert_true(gstt_approx(results[1].correlation().unwrap(), 0.02))
  assert_true(gstt_approx(results[2].correlation().unwrap(), 0.05))
  assert_true(gstt_approx(results[3].correlation().unwrap(), 0.01))
  gstt_check_row(
    results[0],
    name="UP",
    n_genes=8,
    direction="Down",
    p_value=0.54248057689052853,
    fdr=0.72330743585403801,
  )
  gstt_check_row(
    results[1],
    name="DOWN",
    n_genes=5,
    direction="Down",
    p_value=0.89291766898604674,
    fdr=0.89291766898604674,
  )
  gstt_check_row(
    results[2],
    name="NULLSET",
    n_genes=5,
    direction="Up",
    p_value=0.51423951305031756,
    fdr=0.72330743585403801,
  )
  gstt_check_row(
    results[3],
    name="SINGLE",
    n_genes=1,
    direction="Up",
    p_value=0.4596798481474853,
    fdr=0.72330743585403801,
  )
}

///|
test "cameraPR single set has no FDR or Correlation" {
  let results = @src.camera_pr(gstt_statistic_s(), [
    @src.CameraGeneSet::create("one", gstt_index_six()),
  ])
  assert_eq(results.length(), 1)
  assert_true(results[0].correlation().is_none())
  assert_true(results[0].fdr().is_none())
  gstt_check_row(
    results[0],
    name="one",
    n_genes=6,
    direction="Down",
    p_value=0.10007499413141692,
    fdr=-1.0,
  )
}

///|
test "cameraPR sorts by p-value and validates inputs" {
  let results = @src.camera_pr(gstt_statistic_s(), gstt_sets())
  assert_eq(results[0].name(), "SINGLE")
  assert_eq(results[1].name(), "NULLSET")
  assert_eq(results[2].name(), "UP")
  assert_eq(results[3].name(), "DOWN")
  assert_true(
    gstt_raises(fn() raise @src.GeneSetTestError {
      @src.camera_pr([1.0, 0.0], gstt_sets()) |> ignore
    }),
  )
  assert_true(
    gstt_raises(fn() raise @src.GeneSetTestError {
      @src.camera_pr(gstt_statistic_s(), gstt_sets(), fixed_inter_gene_cor=1.0)
      |> ignore
    }),
  )
  assert_true(
    gstt_raises(fn() raise @src.GeneSetTestError {
      @src.camera_pr(
        gstt_statistic_s(),
        gstt_sets(),
        per_set_inter_gene_cor=Some([0.01, 0.01]),
      )
      |> ignore
    }),
  )
  assert_true(
    gstt_raises(fn() raise @src.GeneSetTestError {
      let bad = @src.CameraGeneSet::create("bad", [0, 1, 40])
      @src.camera_pr(gstt_statistic_s(), [bad]) |> ignore
    }),
  )
}

///|
/// 40 x 8 expression matrix from the R generator (set.seed(2026), group
/// effects added on columns 5-8 for the UP/DOWN sets).
fn gstt_y_data() -> Array[Double] {
  [
    [
      0.52058907291852285, -1.1572538576234426, -1.1555774607705176, 1.9062041885624212,
      1.4054062988468312, -0.43202637747734046, 0.26324498721500261, 0.6197125993564252,
    ],
    [
      -1.07969076235228, 0.77799877400621342, 1.5089163463926982, 2.1317977184788686,
      1.134730147637365, 2.4585911476867994, 2.1104211894466625, 3.3596632579562025,
    ],
    [
      0.13923811501927341, -1.2042217859730673, -1.1554961602573557, 0.23145636301535205,
      -2.8739539566200927, -2.2293519860593207, -0.80423637871259657, -3.6210647627137966,
    ],
    [
      -0.08474878494857653, 0.30671410214768002, -1.5568118102066841, 0.896945737827196,
      0.77431820082489955, -1.8569622098408307, 0.73330951363310293, 0.079559499355537452,
    ],
    [
      -0.66663961528459614, -0.83365854064205591, -0.055247925223556642, -1.7387927123515201,
      1.5981097584827479, 1.7103258270188468, 1.9954883084583575, 2.5125931373067534,
    ],
    [
      -2.5160890320094564, 1.4181411635037484, 0.84919316270271517, 0.46884794013553166,
      -0.21513783053589386, 1.7671412378691886, -0.021215499595033872, -0.21910193815088014,
    ],
    [
      -0.73514679745667699, 0.71178611785071821, -0.011096722503659213, -0.54414767319773871,
      0.60525613051401106, -0.44268289730081639, 0.13850173298157384, -0.82593869467054037,
    ],
    [
      -1.020122263135093, -0.40249772286237345, -0.76031383178071854, -0.1654141538084471,
      -0.61456768364747549, -0.58886440200156642, -1.0701270699532761, -1.0135782441118413,
    ],
    [
      0.11355444129730689, 0.79936221871450863, 1.1757992516762399, 0.55216656257281327,
      1.2755770685755883, 1.8794117416470988, 3.2360620669768574, 3.282066246028311,
    ],
    [
      -0.47379098184009499, 0.42648443555694809, 2.4144497090347583, 1.0333302116034961,
      -0.25694263718177718, 1.6530612818618051, -0.53981009473541886, -0.21214449917961803,
    ],
    [
      -0.40821470433792817, -1.1698986699356224, -1.9585127682261312, -0.046178804783545607,
      -2.3920727363495593, -2.8713545686297435, -2.2962289146560675, -0.81512276337359468,
    ],
    [
      -0.7304332785936144, -0.20620016029301858, 1.4735898934298002, 2.6387073125957894,
      -0.68670196174465625, 0.78066827085778212, 0.68821880334963759, -0.65061243677228664,
    ],
    [
      -0.22143659940617352, -0.93061030924361188, -0.4748839614987303, 0.58900574922501769,
      0.87584034013081791, -0.61387786931866162, -0.39900730372786841, 1.2873672700585292,
    ],
    [
      -0.22581652442895145, 0.44939756579173729, 0.98117098378587364, -0.2023773342292117,
      1.5189420139575078, 1.6724086886939511, 1.8965151471505128, 2.5389332156633175,
    ],
    [
      -2.5468814461282956, -0.64480645456231322, -1.8243552195365609, 0.44156013800238225,
      0.068357590868096174, -0.35521601590155316, 1.529534040595421, 0.025540553803906248,
    ],
    [
      1.3470014992970318, -0.23142268981095807, -0.26061081735226788, -0.10025770208283445,
      0.13702489281584826, 0.78093407436683959, -0.63463645844883376, 1.1639475327246045,
    ],
    [
      0.61640814584996906, -1.2363680114453701, -0.95884531257366012, -1.0934039919911085,
      -3.8587685188705816, -1.3913298282145736, -2.2268545890977283, -2.3676257321171712,
    ],
    [
      0.21756433830788766, -0.96095529822301518, -0.49029571006652489, 0.50324733858265003,
      0.34627657136052015, 1.0761790174734789, 1.5010840927328022, -0.36336956950137561,
    ],
    [
      -0.80471883040026315, 0.1339562815695613, -1.0772967136006564, 0.9492400911176353,
      1.5016338199469252, 1.0655576523132153, -0.34246348237863444, -0.0067533647454471228,
    ],
    [
      0.68974677762581926, -0.99905272296928349, 0.36942591054444629, 0.38280205602231232,
      1.9862618796705998, 4.161646412829259, 1.4462724439885175, 1.8995275011353587,
    ],
    [
      -0.32867201409023966, -0.14147069081321334, 1.2162664811393504, 0.37157821089805299,
      -0.77524519751102983, 0.070456482052291131, 0.25762687069895868, -0.14859207836111962,
    ],
    [
      -0.16468157692584029, 0.16732989964881301, -0.49384323561410304, 0.15721640165142459,
      1.3111581356761979, -2.5352252195931748, 0.61312175865245711, -0.18295275543693362,
    ],
    [
      -1.3920288797713043, -0.19878861267168074, -0.22778472097718763, 0.84793988997537484,
      -1.739475443445158, -2.5413341447233915, -2.6522873667775744, -2.8581478376537497,
    ],
    [
      1.4658247612013915, -0.29120718809209706, -1.1116488870472734, -1.2833926565051761,
      1.0140170076257526, -0.77593957354114551, -0.13828102908689452, 0.958159182178004,
    ],
    [
      0.048206808225443779, -1.7343970342492883, 0.99583308940962401, 1.1582786080896008,
      0.21559420847465305, -0.29595541823976274, -0.12363800186919924, -0.83907694183230297,
    ],
    [
      1.9080838346419888, -0.27272833173260363, 0.56182889094608723, -0.90910621316598828,
      0.6134046260919882, 1.1513729054233763, -0.9647052163833707, -0.75972427931168318,
    ],
    [
      1.7309445272373212, -1.7999294816968234, -1.5642785565158785, 0.33499070516858143,
      3.4459884810620274, 1.4311121324684843, 2.0190431840153162, 1.133485180463786,
    ],
    [
      0.058145837259296215, 1.1527409768030377, 1.4810770191053098, 0.77585408191640737,
      0.65808590074518103, -0.80578634002503213, 1.6916923531990611, 1.6590256686813898,
    ],
    [
      0.64532816177517105, -1.003319485592004, -1.1692051501122338, 0.23714015227600846,
      0.37523449649149093, 0.492593861841638, 0.88580907147416865, 0.55694227721051914,
    ],
    [
      1.7256289865063026, 0.14821004429253856, -1.0520135374663548, -1.5598701345871997,
      -0.67468827593817926, 0.70075833733224313, -0.72993014548252566, 0.60138612966221472,
    ],
    [
      -0.52896691745624591, 0.51949668074917621, -1.471845524829025, 0.024326751188882084,
      2.5804805732962741, 1.7568747084499159, 2.3518846830980129, 0.030213178895757142,
    ],
    [
      0.16639202561255087, 0.0054362944757612797, -1.1934369207641855, 0.33463518915828694,
      -0.41422933222661223, 0.41248117692247616, -0.24733089326019458, -0.38010774897368182,
    ],
    [
      -0.25472357476394636, 1.3470208346539403, -0.34797942469193677, 0.98941099914713859,
      0.18050398539241863, -1.501834690641201, -2.1635918317233216, -0.84222645593767176,
    ],
    [
      0.3327823594379094, -0.84703341729599568, 0.53282254314941857, 0.50037511354159758,
      0.34285764365025712, 0.36767911462303332, 0.87236056020824637, 1.3896346541629399,
    ],
    [
      0.1824328128276346, 0.44339801731577189, -0.71098857585084252, 0.57690786545545159,
      -0.79802466190715593, -0.57060286918970027, 0.60593821388331714, -2.0115084247761899,
    ],
    [
      1.1645937495180019, -0.9771494683233134, -0.19898913763240356, 1.1426909401191054,
      1.6410361662536537, 0.73881651771436818, 0.92724081261281954, -2.336308503763895,
    ],
    [
      0.5934928507278675, 2.124491135513606, -0.2813733235184645, -0.70744781442297744,
      -0.4600051116512518, -2.2682913256649466, 0.4580444366512087, -1.3684597701113821,
    ],
    [
      -0.89177972921736914, 0.68769896056154145, 1.1709262752049787, -0.61855964098112792,
      0.41352901139382192, 2.6140410400524661, 1.1988859757851043, 0.62690573964638752,
    ],
    [
      0.57725331035740268, -0.34336822018039642, 2.28873318592, -0.0083141961424556977,
      -1.0449156740717873, 0.40363580599992022, -0.49712315299553933, 0.73738855008549076,
    ],
    [
      -0.8249055810580429, 0.78516947267394588, -1.0688566144990765, 0.33891973395268893,
      0.041653665885759195, 2.0515793266993096, -1.2129678148543903, 0.45470934175017702,
    ],
  ].flatten()
}

///|
fn gstt_design() -> @src.Matrix {
  @src.Matrix::new(8, 2, [
    1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
  ])
}

///|
test "camera.default fixed correlation parametric test matches limma" {
  let y = @src.Matrix::new(40, 8, gstt_y_data())
  let results = @src.camera_default(y, gstt_sets(), gstt_design(), sort=false)
  assert_eq(results.length(), 4)
  assert_true(results[0].correlation().is_none())
  gstt_check_row(
    results[0],
    name="UP",
    n_genes=8,
    direction="Up",
    p_value=2.6775554837069889e-07,
    fdr=1.0710221934827956e-06,
    tolerance=1.0e-11,
  )
  gstt_check_row(
    results[1],
    name="DOWN",
    n_genes=5,
    direction="Down",
    p_value=2.3149688898764702e-05,
    fdr=4.6299377797529404e-05,
    tolerance=1.0e-10,
  )
  gstt_check_row(
    results[2],
    name="NULLSET",
    n_genes=5,
    direction="Up",
    p_value=0.83647204836549549,
    fdr=0.83647204836549549,
  )
  gstt_check_row(
    results[3],
    name="SINGLE",
    n_genes=1,
    direction="Up",
    p_value=0.091703973041782433,
    fdr=0.12227196405570991,
  )
}

///|
test "camera.default ranks match limma" {
  let y = @src.Matrix::new(40, 8, gstt_y_data())
  let results = @src.camera_default(
    y,
    gstt_sets(),
    gstt_design(),
    use_ranks=true,
    sort=false,
  )
  gstt_check_row(
    results[0],
    name="UP",
    n_genes=8,
    direction="Up",
    p_value=3.5274072974447538e-05,
    fdr=0.00014109629189779015,
    tolerance=1.0e-10,
  )
  gstt_check_row(
    results[1],
    name="DOWN",
    n_genes=5,
    direction="Down",
    p_value=0.00053852017834995819,
    fdr=0.0010770403566999164,
    tolerance=1.0e-9,
  )
  gstt_check_row(
    results[2],
    name="NULLSET",
    n_genes=5,
    direction="Up",
    p_value=0.87211808067269936,
    fdr=0.87211808067269936,
  )
  gstt_check_row(
    results[3],
    name="SINGLE",
    n_genes=1,
    direction="Up",
    p_value=0.11891756267142568,
    fdr=0.15855675022856758,
  )
}

///|
test "camera.default estimated correlation matches limma" {
  let y = @src.Matrix::new(40, 8, gstt_y_data())
  let results = @src.camera_default(
    y,
    gstt_sets(),
    gstt_design(),
    inter_gene_cor=None,
    sort=false,
  )
  assert_true(
    gstt_approx(
      results[0].correlation().unwrap(),
      -0.054866578412780087,
      tolerance=1.0e-8,
    ),
  )
  assert_true(
    gstt_approx(
      results[1].correlation().unwrap(),
      -0.07911784914579284,
      tolerance=1.0e-8,
    ),
  )
  assert_true(
    gstt_approx(
      results[2].correlation().unwrap(),
      0.4099441607393175,
      tolerance=1.0e-8,
    ),
  )
  assert_true(results[3].correlation().is_none())
  gstt_check_row(
    results[0],
    name="UP",
    n_genes=8,
    direction="Up",
    p_value=0.00067928066255676844,
    fdr=0.0027171226502270738,
    tolerance=1.0e-8,
  )
  gstt_check_row(
    results[1],
    name="DOWN",
    n_genes=5,
    direction="Down",
    p_value=0.0026968142877907718,
    fdr=0.0053936285755815435,
    tolerance=1.0e-8,
  )
  gstt_check_row(
    results[2],
    name="NULLSET",
    n_genes=5,
    direction="Up",
    p_value=0.89664681919393296,
    fdr=0.89664681919393296,
    tolerance=1.0e-8,
  )
  gstt_check_row(
    results[3],
    name="SINGLE",
    n_genes=1,
    direction="Up",
    p_value=0.13431177304035277,
    fdr=0.17908236405380368,
    tolerance=1.0e-8,
  )
}

///|
test "camera.default estimated correlation with ranks matches limma" {
  let y = @src.Matrix::new(40, 8, gstt_y_data())
  let results = @src.camera_default(
    y,
    gstt_sets(),
    gstt_design(),
    inter_gene_cor=None,
    use_ranks=true,
    sort=false,
  )
  gstt_check_row(
    results[0],
    name="UP",
    n_genes=8,
    direction="Up",
    p_value=0.0054189174500864938,
    fdr=0.021675669800345975,
    tolerance=1.0e-8,
  )
  gstt_check_row(
    results[1],
    name="DOWN",
    n_genes=5,
    direction="Down",
    p_value=0.012563483807419572,
    fdr=0.025126967614839145,
    tolerance=1.0e-8,
  )
  gstt_check_row(
    results[2],
    name="NULLSET",
    n_genes=5,
    direction="Up",
    p_value=0.91860110442411647,
    fdr=0.91860110442411647,
    tolerance=1.0e-8,
  )
  gstt_check_row(
    results[3],
    name="SINGLE",
    n_genes=1,
    direction="Up",
    p_value=0.16993330033451662,
    fdr=0.22657773377935547,
    tolerance=1.0e-8,
  )
}

///|
test "camera.default non-directional ranks match limma" {
  let y = @src.Matrix::new(40, 8, gstt_y_data())
  let results = @src.camera_default(
    y,
    gstt_sets(),
    gstt_design(),
    use_ranks=true,
    directional=false,
    sort=false,
  )
  assert_eq(results[0].direction(), "")
  gstt_check_row(
    results[0],
    name="UP",
    n_genes=8,
    direction="",
    p_value=0.00019767924812788958,
    fdr=0.0007907169925115583,
    tolerance=1.0e-9,
  )
  gstt_check_row(
    results[1],
    name="DOWN",
    n_genes=5,
    direction="",
    p_value=0.0062976026456839571,
    fdr=0.012595205291367914,
    tolerance=1.0e-8,
  )
  gstt_check_row(
    results[2],
    name="NULLSET",
    n_genes=5,
    direction="",
    p_value=0.94626685929890586,
    fdr=0.94626685929890586,
  )
  gstt_check_row(
    results[3],
    name="SINGLE",
    n_genes=1,
    direction="",
    p_value=0.07041550415687195,
    fdr=0.09388733887338873,
    tolerance=1.0e-8,
  )
}

///|
test "camera.default validates design and sets" {
  let y = @src.Matrix::new(40, 8, gstt_y_data())
  let bad_design = @src.Matrix::new(7, 2, Array::make(14, 1.0))
  assert_true(
    gstt_raises(fn() raise @src.GeneSetTestError {
      @src.camera_default(y, gstt_sets(), bad_design) |> ignore
    }),
  )
  assert_true(
    gstt_raises(fn() raise @src.GeneSetTestError {
      @src.camera_default(y, gstt_sets(), gstt_design(), contrast=Some(5))
      |> ignore
    }),
  )
  assert_true(
    gstt_raises(fn() raise @src.GeneSetTestError {
      @src.camera_default(
        y,
        gstt_sets(),
        gstt_design(),
        inter_gene_cor=Some(1.0),
      )
      |> ignore
    }),
  )
}
