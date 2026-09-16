///|
/// WGCNA frozen-value tests
/// All values from Python (np.random.seed(42) expr, seed(123) traits).

fn _wgcna_approx_eq(a : Double, b : Double) -> Bool {
  let d = a - b
  let ad = if d < 0.0 { -d } else { d }
  ad <= 0.0001
}

fn _wgcna_close_eq(a : Double, b : Double) -> Bool {
  let d = a - b
  let ad = if d < 0.0 { -d } else { d }
  ad <= 0.001
}

let _EXPR : Array[Array[Double]] = [
  [0.1543384306, 1.3635813302, -1.1547448005, 1.5175847015, -0.9124631149, -1.3063777269, 0.6290608678, -1.1691461191, -1.0652111510, -1.5255368171, 2.7153788530, -0.1799695488, -0.9286071262, -0.6356702562, -1.0309509908],
  [-0.3751624105, 0.9262979206, -1.3214654032, 0.3020577355, -1.5015382782, -0.4564848639, -0.2734336902, -1.5190993813, -1.5428480511, -2.0690082011, 3.0025489942, -0.5226014608, -0.8617392528, 0.0206505211, -0.3143389494],
  [-0.3948719520, 2.3588525730, -1.4853376819, 0.2272143088, -0.3841139715, -1.7997310734, -0.0561137424, -2.5326985456, -1.4945669825, -1.6983862209, 3.3163209894, 0.1385794686, -1.1074124032, 0.1704180183, -0.9077530697],
  [-0.4539407501, 1.2023940953, -0.9500279564, 0.9278789181, -1.6769065052, -1.0272672637, -0.3530866801, -1.8913244838, -0.5246358136, -1.2813170776, 3.4127277589, -0.3667134337, -1.2041944499, 0.4866015818, 0.3192804890],
  [-0.3336057648, 1.3398839924, -2.0317565565, 0.1579664612, -0.3891235164, -0.5111892341, -0.1965506007, -1.0510970347, -0.6496559455, -2.1193767162, 3.1277855021, 0.8219136112, -1.0675012815, 1.1032916940, -1.4783646266],
  [0.1457877327, 0.7805235742, 0.7595055768, -0.5907209889, -1.6847663461, -0.4398979805, -2.5837491119, -1.0792822523, 1.3811140118, 1.5697144548, -0.5497217693, 1.1495694145, 0.5378757403, 3.4297768683, 1.6001236103],
  [-0.2166247448, 1.2213225353, 0.5579827051, -0.8004324504, -0.8870359653, -1.0618195106, -2.6142452592, -1.6877016385, 1.6428058493, 1.8566676895, -1.0065286185, 0.4815456943, 0.2021429272, 3.2935183355, 1.2628470379],
  [-0.0631380910, 1.6800929907, 0.9962981585, -0.5078261818, -0.7282048467, -1.2894476442, -2.7755623354, -1.7881141697, 2.8718701772, 1.8777807738, -0.1480695763, 0.6745124708, -0.2108388334, 4.2660683774, 1.7194564101],
  [0.1303524541, 0.2823063127, 1.6104064075, -1.3375269085, -0.3975533419, 0.7651657763, -3.2575735602, -2.1013781394, 1.6900748035, 1.7222234291, -1.0741749630, 0.7261498431, -0.1576516714, 3.9314531855, 0.8837777766],
  [0.5098036831, 0.3453733939, 0.7479784940, -0.2298427684, -1.3064140470, -0.2163320693, -2.1087340205, -2.6219708919, 1.7325660502, 2.1039026533, 0.0920681884, 0.0733930002, -0.2867281211, 3.9556277530, 1.4919822303],
  [0.2915392941, 1.4154266775, 1.7483678555, -0.1483758780, 1.1297922055, 1.4917541135, 1.7734320051, 1.6429288652, -2.7889358742, -0.8178548797, -0.4769764014, -1.5303850787, 0.4787189896, 0.7093283026, -1.4983732358],
  [0.3726833324, 1.6532326528, 3.0367767077, -0.3871967846, 0.6063878867, 1.4041726077, 0.4326396070, 1.3674615500, -2.0227081382, -1.0077862844, 0.4239560582, -1.9174264346, 0.6261882386, 0.9873410451, -0.6199767172],
  [0.4791265428, 0.8136237946, 1.5529339673, -0.0232665190, 0.8715245762, 2.2059300695, 1.0771635618, 1.3695979484, -2.6166809846, -1.9035552964, -0.2128930423, -1.4957279834, 0.0064680721, 0.4968000724, -1.8934708375],
  [0.3589515588, 0.8002738547, 2.1652427693, -0.2353983674, 0.4117708200, 2.0278235027, 1.1209370127, 1.9475380263, -1.6663830996, -1.8349663680, -0.4585480862, -1.6664097470, 0.1563141755, 1.3771933049, -0.0536955551],
  [0.4517381243, 1.8099853929, 2.5653810981, 0.0611928990, 0.8256213465, 2.2284144327, 0.4541321422, 1.2876031013, -2.4359658995, -1.1051946143, 1.1676937170, -2.8575599769, 0.2425512952, 0.3133115261, -2.2160272333],
  [0.0055652548, 1.0693861639, -0.3485337109, 0.5179663988, 1.7225986054, -1.2913608634, -1.8011873457, 1.4213534938, -0.8887694105, 1.7226146138, 0.8040848949, -0.4201773691, 0.9572179271, -0.3794148603, 1.2571667838],
  [-0.9351704129, 0.9798779336, 0.4428323176, 1.3084958505, 0.7826515274, -1.0934281655, -2.1468892961, 1.0719029575, 0.3197578834, 0.8531334246, -0.1433165934, 1.0513248978, 1.9250688682, 0.4677110023, 0.0712651334],
  [-0.7810270801, 1.6707017289, -0.1634960547, 1.0975279675, 1.7701167577, -1.3896427833, -1.9391793185, -0.2220660962, -1.0751630574, 0.5243584934, -0.1367662072, 1.4085995762, 0.1489200807, -0.2685438713, 0.8963204051],
  [0.1817266008, 0.3193150787, 0.7719205542, 0.8807347839, 0.8920454054, -0.6951258104, -1.8098867927, 1.0984591352, -0.5280681942, 0.4579857707, 0.5438840564, 0.9234592615, 1.6569991778, -0.6674293774, 1.8974668037],
  [-1.5149539434, 0.9613536068, 0.4844972813, 1.0161141873, 1.0714499711, -1.0302386727, -2.1559171078, 1.1038851953, -0.1381681882, 0.8291503121, 0.1406705861, 1.0421938620, 1.0176405301, 0.3579094315, 1.1457645374],
  [-0.3333604577, -0.2951619529, 0.3019623885, 0.0329160341, -0.1834821962, 0.2852810594, 0.7891076042, -0.5890522679, 0.3033483822, 0.1243198108, -0.6096621860, 0.7123964837, 0.2139210468, 0.5779344186, 0.4237616566],
  [0.0916389686, 0.3259050528, -0.2268177924, -0.1101859224, -0.2381029264, 0.2751153502, 0.4478536690, -0.7023769151, 1.0759933295, -0.2775925537, -1.1079158907, 0.7420634945, 0.1970442193, 0.4832395091, 0.0851950077],
  [0.0750136613, -0.4636986186, -0.0337821350, -0.6108499544, 0.3145284675, 0.1530886770, -0.2624734421, -0.4539594943, 0.2362654188, -0.0564461394, -0.9119317821, 0.2848516635, -0.0763038421, -0.0822919871, -0.4644968998],
  [0.1971620164, -0.7391136036, -0.7754163012, -0.6314912092, -0.2797549750, 0.3820711505, 0.8879532648, 0.1355632382, -0.0501695734, 0.2159080331, -1.0020862667, 0.1537514897, -0.3431164471, 0.3325388807, -0.6425932187],
  [0.3408103049, 0.7512980237, -0.1260644883, -0.0714132376, 0.1720405967, 0.0260071318, 0.2623213972, -0.2869703730, 0.0786377409, -0.1610887549, -0.4885664972, 0.4120072033, 0.5267846762, 0.6508150136, 0.8476134818],
  [0.3382790588, 1.0871772935, 0.0724370242, 0.9910124659, -0.0424224926, -0.6752979046, 0.3649117022, -1.0355159427, 0.2215289691, -0.1188318877, 0.0562374247, 0.5150433340, 1.0919076373, 0.0726947954, -0.0576950800],
  [0.2727455046, 0.8969765609, -0.6793505822, 0.8118403822, 0.9514467104, -0.4900248096, -0.1919592396, 0.7033681065, 0.4271365573, 0.8796425991, -0.8161732979, 0.0444778346, 0.1564440676, -0.3790268268, 0.0057240121],
  [1.0333778064, 0.1172067605, -0.0904237212, -0.0437411848, 0.6189460671, 0.1003704555, 0.1022869790, -0.7406251659, 1.1232448909, 0.4268914491, -0.3932027368, 1.1197413339, 0.2116595348, 0.1871314683, 0.2645159885],
  [1.7523268027, 1.5286873964, -0.1437160529, 0.3818964748, 0.6844146248, 0.4288787953, 0.1821462946, 0.3694576495, 1.0136187233, -0.6186353001, -0.6105887249, -0.6754507428, 0.0191188005, -0.0437459958, 0.9819352739],
  [0.7590002305, 1.4653247479, -0.7092847078, -0.9555802204, 0.3339528006, -0.0634042590, 0.6482606509, -1.0072891306, 0.4398464601, -0.3915003071, 0.3158768059, 0.5274644692, -0.3161176755, -0.6594505824, -0.2988500131],
]
let _GENES : Array[String] = [
  "Gene0","Gene1","Gene2","Gene3","Gene4","Gene5","Gene6","Gene7","Gene8","Gene9",
  "Gene10","Gene11","Gene12","Gene13","Gene14","Gene15","Gene16","Gene17","Gene18","Gene19",
  "Gene20","Gene21","Gene22","Gene23","Gene24","Gene25","Gene26","Gene27","Gene28","Gene29",
]
let _SAMPLES : Array[String] = [
  "S0","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14",
]
let _POWERS : Array[Double] = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
let _TRAITS : Array[Array[Double]] = [
  [-0.5428153017, 0.4986727233, 0.1414892490, -0.7531473570, -0.2893001260, 0.8257182685, -1.2133396217, -0.2144563144, 0.6329681294, -0.4333702011, -0.3394430758, -0.0473544845, 0.7456948131, -0.3194509983, -0.2219909798],
  [-0.2171756378, 1.1029650414, 1.0933930445, 0.5020269489, 0.1930931996, 0.3686842879, 0.7453660141, -0.4679169342, 0.5879145224, -0.6269403339, -0.3188757512, 0.4535525979, -0.7143403501, -0.0700343601, -0.4308774479],
  [-0.1278096853, -1.3992945527, -0.8857665523, -0.3499386173, 0.4637312159, -0.0868178414, 0.0014229579, 0.3441113556, -0.4397681715, 0.1418136619, -0.4026832590, -0.8638347471, -0.1954498969, 0.2869029312, 0.1692945255],
]
let _TRAIT_NAMES : Array[String] = ["Weight", "Height", "Age"]

// ============ Adjacency ============
test "wgcna_adjacency_power2_frozen" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let wd = @src.adjacency(wd, 2.0)
  assert_true(_wgcna_approx_eq(wd.adjacency_matrix[0][1], 0.7824233087))
  assert_true(_wgcna_approx_eq(wd.adjacency_matrix[5][6], 0.9246171560))
  assert_true(_wgcna_approx_eq(wd.adjacency_matrix[0][0], 1.0))
  assert_true(_wgcna_approx_eq(wd.adjacency_matrix[29][29], 1.0))
}

test "wgcna_adjacency_symmetric" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let wd = @src.adjacency(wd, 6.0)
  let mut i = 0
  let mut ok = true
  while i < 30 {
    let mut j = 0
    while j < 30 {
      if !_wgcna_approx_eq(wd.adjacency_matrix[i][j], wd.adjacency_matrix[j][i]) {
        ok = false
        break
      }
      j = j + 1
    }
    if !ok { break }
    i = i + 1
  }
  assert_true(ok)
}

// ============ TOM ============
test "wgcna_tom_power2_frozen" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let wd = @src.adjacency(wd, 2.0)
  let wd = @src.tom_similarity(wd)
  assert_true(_wgcna_approx_eq(wd.tom_matrix[0][1], 0.3304018720))
  assert_true(_wgcna_approx_eq(wd.tom_matrix[5][6], 0.3500655610))
  assert_true(_wgcna_approx_eq(wd.tom_matrix[0][0], 1.0))
}

test "wgcna_tom_bounds" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let wd = @src.adjacency(wd, 4.0)
  let wd = @src.tom_similarity(wd)
  let mut ok = true
  let mut i = 0
  while i < 30 && ok {
    let mut j = 0
    while j < 30 && ok {
      if wd.tom_matrix[i][j] < -0.001 || wd.tom_matrix[i][j] > 1.001 {
        ok = false
      }
      j = j + 1
    }
    i = i + 1
  }
  assert_true(ok)
}

// ============ Pick Soft Threshold ============
test "wgcna_pick_soft_threshold_basic" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let best = @src.pick_soft_threshold(wd, _POWERS, 0.8)
  assert_true(best >= 2.0)
  assert_true(best <= 20.0)
}

test "wgcna_pick_soft_threshold_low_r2" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let low = [1.0, 2.0, 4.0, 6.0, 8.0]
  let best = @src.pick_soft_threshold(wd, low, 0.5)
  assert_true(best >= 1.0)
}

test "wgcna_pick_soft_threshold_small_dataset" {
  let expr : Array[Array[Double]] = [[1.0,2.0,3.0],[2.0,3.0,4.0]]
  let wd = @src.create_wgcna_data(expr, ["A","B"], ["s1","s2","s3"])
  let best = @src.pick_soft_threshold(wd, [1.0, 4.0, 6.0], 0.8)
  assert_true(best == 1.0)
}

// ============ Blockwise Modules ============
test "wgcna_blockwise_modules_count" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let wd = @src.blockwise_modules(wd, 2.0, 10)
  assert_true(wd.modules.length() == 3)
  assert_true(wd.module_colors.length() == 3)
  match wd.modules.get("Module_0") {
    Option::None => assert_true(false)
    Some(g) => assert_true(g.length() == 10)
  }
}

test "wgcna_blockwise_modules_colors" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let wd = @src.blockwise_modules(wd, 6.0, 10)
  match wd.module_colors.get("Module_0") {
    Option::None => assert_true(false)
    Some(c) => assert_true(c == "red")
  }
}

// ============ Module Eigengenes ============
test "wgcna_module_eigengenes_count" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let wd = @src.blockwise_modules(wd, 2.0, 10)
  let mes = @src.module_eigengenes(wd)
  assert_true(mes.length() == 3)
}

test "wgcna_module_eigengenes_normalized" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let wd = @src.blockwise_modules(wd, 2.0, 10)
  let mes = @src.module_eigengenes(wd)
  let mut i = 0
  while i < mes.length() {
    let me = mes[i].eigengene
    let mut nrm = 0.0
    let mut j = 0
    while j < me.length() { nrm = nrm + me[j]*me[j]; j = j + 1 }
    assert_true(_wgcna_close_eq(nrm.sqrt(), 1.0))
    i = i + 1
  }
}

test "wgcna_module_eigengenes_frozen" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let wd = @src.blockwise_modules(wd, 2.0, 10)
  let mes = @src.module_eigengenes(wd)
  assert_true(mes[0].module_name == "Module_0")
  assert_true(_wgcna_close_eq(mes[0].eigengene[0], -0.023894))
  assert_true(_wgcna_close_eq(mes[0].eigengene[1], 0.306334))
  assert_true(mes[1].module_name == "Module_1")
  assert_true(_wgcna_close_eq(mes[1].eigengene[0], -0.036700))
}

// ============ Module-Trait Corr ============
test "wgcna_module_trait_corr_count" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let wd = @src.blockwise_modules(wd, 2.0, 10)
  let corr = @src.module_trait_correlation(wd, _TRAITS, _TRAIT_NAMES)
  assert_true(corr.length() == 3)
}

test "wgcna_module_trait_corr_range" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let wd = @src.blockwise_modules(wd, 2.0, 10)
  let corr = @src.module_trait_correlation(wd, _TRAITS, _TRAIT_NAMES)
  let mod_names = corr.keys().to_array()
  let mut i = 0
  while i < mod_names.length() {
    match corr.get(mod_names[i]) {
      Option::None => assert_true(false)
      Some(tmap) => {
        let tnames = tmap.keys().to_array()
        let mut j = 0
        while j < tnames.length() {
          match tmap.get(tnames[j]) {
            Option::None => assert_true(false)
            Some(c) => assert_true(c >= -1.001 && c <= 1.001)
          }
          j = j + 1
        }
      }
    }
    i = i + 1
  }
}

test "wgcna_module_trait_corr_frozen" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let wd = @src.blockwise_modules(wd, 2.0, 10)
  let corr = @src.module_trait_correlation(wd, _TRAITS, _TRAIT_NAMES)
  match corr.get("Module_0") {
    Option::None => assert_true(false)
    Some(t) => {
      match t.get("Weight") {
        Option::None => assert_true(false)
        Some(c) => assert_true(_wgcna_close_eq(c, 0.122107))
      }
    }
  }
}

// ============ kME ============
test "wgcna_module_membership_count" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let wd = @src.adjacency(wd, 2.0)
  let wd = @src.blockwise_modules(wd, 2.0, 10)
  let members = @src.get_module_membership(wd)
  assert_true(members.length() == 30)
}

test "wgcna_module_membership_kme_range" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let wd = @src.adjacency(wd, 2.0)
  let wd = @src.blockwise_modules(wd, 2.0, 10)
  let members = @src.get_module_membership(wd)
  let mut i = 0
  while i < members.length() {
    assert_true(members[i].kme >= -0.001 && members[i].kme <= 1.001)
    assert_true(members[i].connectivity >= -0.001)
    i = i + 1
  }
}

test "wgcna_module_membership_frozen" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let wd = @src.adjacency(wd, 2.0)
  let wd = @src.blockwise_modules(wd, 2.0, 10)
  let members = @src.get_module_membership(wd)
  assert_true(members[0].gene_name == "Gene0")
  assert_true(members[0].module_name == "Module_0")
  assert_true(_wgcna_close_eq(members[0].kme, 0.370947))
  assert_true(_wgcna_close_eq(members[0].connectivity, 4.779160))
}

// ============ Summary ============
test "wgcna_summary" {
  let wd = @src.create_wgcna_data(_EXPR, _GENES, _SAMPLES)
  let wd = @src.blockwise_modules(wd, 6.0, 10)
  let s = @src.summarize_wgcna(wd)
  assert_true(s.contains("30"))
  assert_true(s.contains("15"))
  assert_true(s.contains("3"))
}
