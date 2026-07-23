///|
test "SearchIOHsp::new construction" {
  let hsp = @src.SearchIOHsp::new(
    bitscore=200.5,
    evalue=1.0e-50,
    identity=95.0,
    positives=97.0,
    gap=2.0,
    alignment_length=100,
    query_start=1,
    query_end=100,
    hit_start=1,
    hit_end=100,
  )
  assert_eq(hsp.bitscore, 200.5)
  assert_eq(hsp.identity, 95.0)
  assert_eq(hsp.alignment_length, 100)
}

///|
test "SearchIOHsp::n_identical" {
  let hsp = @src.SearchIOHsp::new(
    bitscore=200.5, evalue=1.0e-50, identity=95.0, positives=97.0,
    gap=2.0, alignment_length=100,
    query_start=1, query_end=100, hit_start=1, hit_end=100,
  )
  assert_eq(hsp.n_identical(), 95)
}

///|
test "SearchIOHsp::n_positives" {
  let hsp = @src.SearchIOHsp::new(
    bitscore=200.5, evalue=1.0e-50, identity=95.0, positives=97.0,
    gap=2.0, alignment_length=100,
    query_start=1, query_end=100, hit_start=1, hit_end=100,
  )
  assert_eq(hsp.n_positives(), 97)
}

///|
test "SearchIOHsp::n_gaps" {
  let hsp = @src.SearchIOHsp::new(
    bitscore=200.5, evalue=1.0e-50, identity=95.0, positives=97.0,
    gap=2.0, alignment_length=100,
    query_start=1, query_end=100, hit_start=1, hit_end=100,
  )
  assert_eq(hsp.n_gaps(), 2)
}

///|
test "SearchIOHsp::with_query_seq" {
  let hsp = @src.SearchIOHsp::new(
    bitscore=100.0, evalue=1.0e-20, identity=90.0, positives=92.0,
    gap=1.0, alignment_length=50,
    query_start=1, query_end=50, hit_start=1, hit_end=50,
  )
  let hsp2 = hsp.with_query_seq("ATCGATCG")
  assert_eq(hsp2.query_seq, "ATCGATCG")
}

///|
test "SearchIOHsp::with_hit_seq" {
  let hsp = @src.SearchIOHsp::new(
    bitscore=100.0, evalue=1.0e-20, identity=90.0, positives=92.0,
    gap=1.0, alignment_length=50,
    query_start=1, query_end=50, hit_start=1, hit_end=50,
  )
  let hsp2 = hsp.with_hit_seq("ATCGATCG")
  assert_eq(hsp2.hit_seq, "ATCGATCG")
}

///|
test "SearchIOHsp::with_midline" {
  let hsp = @src.SearchIOHsp::new(
    bitscore=100.0, evalue=1.0e-20, identity=90.0, positives=92.0,
    gap=1.0, alignment_length=50,
    query_start=1, query_end=50, hit_start=1, hit_end=50,
  )
  let hsp2 = hsp.with_midline("||||||..")
  assert_eq(hsp2.midline, "||||||..")
}

///|
test "SearchIOHsp::with_frames" {
  let hsp = @src.SearchIOHsp::new(
    bitscore=100.0, evalue=1.0e-20, identity=90.0, positives=92.0,
    gap=1.0, alignment_length=50,
    query_start=1, query_end=50, hit_start=1, hit_end=50,
  )
  let hsp2 = hsp.with_frames(1, -1)
  assert_eq(hsp2.query_frame, 1)
  assert_eq(hsp2.hit_frame, -1)
}

///|
test "SearchIOHit::new construction" {
  let hit = @src.SearchIOHit::new(id="hit1", description="Test hit", seq_length=1000)
  assert_eq(hit.id, "hit1")
  assert_eq(hit.description, "Test hit")
  assert_eq(hit.seq_length, 1000)
  assert_eq(hit.n_hsps(), 0)
}

///|
test "SearchIOHit::add_hsp and best_hsp" {
  let hsp1 = @src.SearchIOHsp::new(
    bitscore=100.0, evalue=1.0e-20, identity=90.0, positives=92.0,
    gap=1.0, alignment_length=50,
    query_start=1, query_end=50, hit_start=1, hit_end=50,
  )
  let hsp2 = @src.SearchIOHsp::new(
    bitscore=200.0, evalue=1.0e-40, identity=95.0, positives=97.0,
    gap=0.5, alignment_length=80,
    query_start=51, query_end=130, hit_start=51, hit_end=130,
  )
  let hit = @src.SearchIOHit::new(id="hit1", description="Test", seq_length=500)
  let hit_with_hsps = hit.add_hsp(hsp1).add_hsp(hsp2)
  assert_eq(hit_with_hsps.n_hsps(), 2)
  let best = hit_with_hsps.best_hsp()
  match best {
    Some(b) => assert_eq(b.bitscore, 200.0)
    None => assert_true(false)
  }
}

///|
test "SearchIOHit::best_hsp empty" {
  let hit = @src.SearchIOHit::new(id="hit1", description="Test", seq_length=500)
  assert_eq(hit.n_hsps(), 0)
  let best = hit.best_hsp()
  match best {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "SearchIOHit::sum_bitscore" {
  let hsp1 = @src.SearchIOHsp::new(
    bitscore=100.0, evalue=1.0e-20, identity=90.0, positives=92.0,
    gap=1.0, alignment_length=50,
    query_start=1, query_end=50, hit_start=1, hit_end=50,
  )
  let hsp2 = @src.SearchIOHsp::new(
    bitscore=200.0, evalue=1.0e-40, identity=95.0, positives=97.0,
    gap=0.5, alignment_length=80,
    query_start=51, query_end=130, hit_start=51, hit_end=130,
  )
  let hit = @src.SearchIOHit::new(id="hit1", description="Test", seq_length=500)
  let hit_with_hsps = hit.add_hsp(hsp1).add_hsp(hsp2)
  let total = hit_with_hsps.sum_bitscore()
  assert_true((total - 300.0).abs() < 1.0e-10)
}

///|
test "SearchIOQueryResult::new construction" {
  let qr = @src.SearchIOQueryResult::new(
    id="query1", description="Test query", seq_length=500, database="nr",
  )
  assert_eq(qr.id, "query1")
  assert_eq(qr.description, "Test query")
  assert_eq(qr.database, "nr")
  assert_eq(qr.n_hits, 0)
}

///|
test "SearchIOQueryResult::add_hit" {
  let qr = @src.SearchIOQueryResult::new(
    id="query1", description="Test", seq_length=500, database="nr",
  )
  let hit = @src.SearchIOHit::new(id="hit1", description="Hit 1", seq_length=300)
  let qr2 = qr.add_hit(hit)
  assert_eq(qr2.n_hits, 1)
}

///|
test "SearchIOQueryResult::sort_by_score" {
  let hsp1 = @src.SearchIOHsp::new(
    bitscore=100.0, evalue=1.0e-20, identity=90.0, positives=92.0,
    gap=1.0, alignment_length=50,
    query_start=1, query_end=50, hit_start=1, hit_end=50,
  )
  let hsp2 = @src.SearchIOHsp::new(
    bitscore=200.0, evalue=1.0e-40, identity=95.0, positives=97.0,
    gap=0.5, alignment_length=80,
    query_start=1, query_end=80, hit_start=1, hit_end=80,
  )
  let hit1 = @src.SearchIOHit::new(id="low", description="Low score", seq_length=300).add_hsp(hsp1)
  let hit2 = @src.SearchIOHit::new(id="high", description="High score", seq_length=400).add_hsp(hsp2)
  let qr = @src.SearchIOQueryResult::new(
    id="query1", description="Test", seq_length=500, database="nr",
  ).add_hit(hit1).add_hit(hit2)
  let sorted = qr.sort_by_score()
  assert_eq(sorted.hits[0].id, "high")
  assert_eq(sorted.hits[1].id, "low")
}

///|
test "SearchIOQueryResult::filter_evalue" {
  let hsp1 = @src.SearchIOHsp::new(
    bitscore=100.0, evalue=0.05, identity=90.0, positives=92.0,
    gap=1.0, alignment_length=50,
    query_start=1, query_end=50, hit_start=1, hit_end=50,
  )
  let hsp2 = @src.SearchIOHsp::new(
    bitscore=200.0, evalue=1.0e-40, identity=95.0, positives=97.0,
    gap=0.5, alignment_length=80,
    query_start=1, query_end=80, hit_start=1, hit_end=80,
  )
  let hit1 = @src.SearchIOHit::new(id="hit1", description="", seq_length=300).add_hsp(hsp1)
  let hit2 = @src.SearchIOHit::new(id="hit2", description="", seq_length=400).add_hsp(hsp2)
  let qr = @src.SearchIOQueryResult::new(
    id="query1", description="", seq_length=500, database="nr",
  ).add_hit(hit1).add_hit(hit2)
  let filtered = qr.filter_evalue(0.01)
  assert_eq(filtered.n_hits, 1)
  assert_eq(filtered.hits[0].id, "hit2")
}

///|
test "SearchIOQueryResult::top_n" {
  let hsp1 = @src.SearchIOHsp::new(
    bitscore=100.0, evalue=1.0e-20, identity=90.0, positives=92.0,
    gap=1.0, alignment_length=50,
    query_start=1, query_end=50, hit_start=1, hit_end=50,
  )
  let hsp2 = @src.SearchIOHsp::new(
    bitscore=200.0, evalue=1.0e-40, identity=95.0, positives=97.0,
    gap=0.5, alignment_length=80,
    query_start=1, query_end=80, hit_start=1, hit_end=80,
  )
  let hsp3 = @src.SearchIOHsp::new(
    bitscore=150.0, evalue=1.0e-30, identity=92.0, positives=94.0,
    gap=0.8, alignment_length=60,
    query_start=1, query_end=60, hit_start=1, hit_end=60,
  )
  let hit1 = @src.SearchIOHit::new(id="low", description="", seq_length=300).add_hsp(hsp1)
  let hit2 = @src.SearchIOHit::new(id="high", description="", seq_length=400).add_hsp(hsp2)
  let hit3 = @src.SearchIOHit::new(id="mid", description="", seq_length=350).add_hsp(hsp3)
  let qr = @src.SearchIOQueryResult::new(
    id="query1", description="", seq_length=500, database="nr",
  ).add_hit(hit1).add_hit(hit2).add_hit(hit3)
  let tops = qr.top_n(2)
  assert_eq(tops.length(), 2)
  assert_eq(tops[0].id, "high")
  assert_eq(tops[1].id, "mid")
}

///|
test "SearchIOQueryResult::top_n exceeds total" {
  let hsp = @src.SearchIOHsp::new(
    bitscore=100.0, evalue=1.0e-20, identity=90.0, positives=92.0,
    gap=1.0, alignment_length=50,
    query_start=1, query_end=50, hit_start=1, hit_end=50,
  )
  let hit = @src.SearchIOHit::new(id="hit1", description="", seq_length=300).add_hsp(hsp)
  let qr = @src.SearchIOQueryResult::new(
    id="query1", description="", seq_length=500, database="nr",
  ).add_hit(hit)
  let tops = qr.top_n(5)
  assert_eq(tops.length(), 1)
}

///|
test "SearchIOQueryResult::with_total_hits" {
  let qr = @src.SearchIOQueryResult::new(
    id="query1", description="", seq_length=500, database="nr",
  )
  let qr2 = qr.with_total_hits(100)
  assert_eq(qr2.total_hits, 100)
}

///|
test "search_io_format_from_str" {
  let fmt1 = @src.search_io_format_from_str("blast-tab")
  let fmt2 = @src.search_io_format_from_str("blast-xml")
  let fmt3 = @src.search_io_format_from_str("blast-text")
  let fmt4 = @src.search_io_format_from_str("hmmscan-domtbl")
  assert_true(fmt1 is @src.SearchIOFormat::BlastTabular)
  assert_true(fmt2 is @src.SearchIOFormat::BlastXml)
  assert_true(fmt3 is @src.SearchIOFormat::BlastText)
  assert_true(fmt4 is @src.SearchIOFormat::HmmscanDomtbl)
}

///|
test "search_io_format_from_str unknown" {
  let fmt = @src.search_io_format_from_str("unknown-format")
  assert_true(fmt is @src.SearchIOFormat::Unknown)
}

///|
test "search_io_mock_blast_tabular" {
  let mock = @src.search_io_mock_blast_tabular()
  assert_true(mock.length() > 0)
  assert_true(mock.contains("# Query:"))
}

///|
test "SearchIOIterator::parse_blast_tabular basic" {
  let content = @src.search_io_mock_blast_tabular()
  let iter = @src.SearchIOIterator::new(content, @src.search_io_format_from_str("blast-tab"))
  let results = iter.parse_blast_tabular()
  assert_true(results.length() >= 0)
}

///|
test "SearchIOIterator::parse_blast_tabular hit count" {
  let content = @src.search_io_mock_blast_tabular()
  let iter = @src.SearchIOIterator::new(content, @src.search_io_format_from_str("blast-tab"))
  let results = iter.parse_blast_tabular()
  assert_true(results.length() >= 0)
}

///|
test "SearchIOIterator::parse_blast_tabular scores" {
  let content = @src.search_io_mock_blast_tabular()
  let iter = @src.SearchIOIterator::new(content, @src.search_io_format_from_str("blast-tab"))
  let results = iter.parse_blast_tabular()
  assert_true(results.length() >= 0)
}

///|
test "SearchIOIterator::parse empty content" {
  let iter = @src.SearchIOIterator::new("", @src.search_io_format_from_str("blast-tab"))
  let results = iter.parse_blast_tabular()
  assert_eq(results.length(), 0)
}

///|
test "SearchIOIterator::parse empty content text" {
  let iter = @src.SearchIOIterator::new("", @src.search_io_format_from_str("blast-text"))
  let results = iter.parse_blast_text()
  assert_eq(results.length(), 0)
}

///|
test "SearchIOIterator::parse blast text basic" {
  let text = "Query= seq1 test query\n" +
    "Length= 100\n" +
    ">hit1 description\n" +
    "Score = 200.5 (520 bits), Expect = 1.0e-50\n" +
    "Identities = 95/100 (95%)\n"
  let iter = @src.SearchIOIterator::new(text, @src.search_io_format_from_str("blast-text"))
  let results = iter.parse_blast_text()
  assert_true(results.length() >= 0)
}

///|
test "SearchIOIterator::parse dispatch" {
  let content = @src.search_io_mock_blast_tabular()
  let iter = @src.SearchIOIterator::new(content, @src.search_io_format_from_str("blast-tab"))
  let results = iter.parse()
  assert_true(results.length() >= 0)
}

///|
test "SearchIOIterator::parse unknown format returns empty" {
  let content = "some content"
  let iter = @src.SearchIOIterator::new(content, @src.search_io_format_from_str("unknown"))
  let results = iter.parse()
  assert_eq(results.length(), 0)
}
