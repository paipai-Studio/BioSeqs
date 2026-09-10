///|
/// Bio.SearchIO.BlatIO demo: parse BLAT `blat-psl` / PSLX output, inspect
/// grouped QueryResult / Hit / HSP (block) objects, and write a round trip.
/// Run with: moon run examples/blat_psl_demo

///|
fn strand_label(value : Int) -> String {
  if value > 0 {
    "+"
  } else if value < 0 {
    "-"
  } else {
    "."
  }
}

///|
fn parse(content : String, pslx : Bool) -> Array[@src.BlatQueryResult] {
  @src.blat_psl_parse(content, pslx~) catch {
    SeqError(message) => abort("PSL parsing failed: " + message)
  }
}

///|
fn write(
  content : Array[@src.BlatQueryResult],
  header : Bool,
  pslx : Bool,
) -> String {
  @src.blat_psl_write(content, header~, pslx~) catch {
    SeqError(message) => abort("PSL writing failed: " + message)
  }
}

///|
fn print_qresults(qresults : Array[@src.BlatQueryResult]) -> Unit {
  for qresult in qresults {
    println(
      "query " +
      qresult.id +
      " (len=" +
      qresult.seq_len.to_string() +
      ", program=" +
      qresult.program +
      "): " +
      qresult.hits.length().to_string() +
      " hit(s)",
    )
    for hit in qresult.hits {
      println(
        "  hit " +
        hit.id +
        " (len=" +
        hit.seq_len.to_string() +
        "): " +
        hit.hsps.length().to_string() +
        " HSP(s)",
      )
      for hsp in hit.hsps {
        let kind = if hsp.query_is_protein { "protein" } else { "nucleotide" }
        println(
          "    HSP [" +
          kind +
          "] score=" +
          hsp.score.to_string() +
          " identity=" +
          hsp.ident_pct.to_string() +
          "% identical=" +
          hsp.ident_num.to_string() +
          " gapopens=" +
          hsp.gapopen_num.to_string() +
          " (" +
          hsp.blocks.length().to_string() +
          " blocks)",
        )
        for block in hsp.blocks {
          println(
            "      query " +
            strand_label(block.query_strand) +
            " [" +
            block.query_start.to_string() +
            ", " +
            block.query_end.to_string() +
            ")  hit " +
            strand_label(block.hit_strand) +
            " [" +
            block.hit_start.to_string() +
            ", " +
            block.hit_end.to_string() +
            ")",
          )
        }
      }
    }
  }
}

///|
fn main {
  println("=== Biopython Bio.SearchIO.BlatIO (blat-psl) Demo ===")

  // A psLayout version 3 header followed by one plus-strand DNA row, one
  // minus-strand DNA row (tStarts listed in descending file order) and one
  // protein query row (hit blocks span 3 nucleotides per residue).
  let psl =
    #|psLayout version 3
    #|
    #|match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
    #|     	match	match	  	count	bases	count	bases	     	name     	size	start	end	name     	size	start	end	count
    #|---------------------------------------------------------------------------------------------------------------------------------------------------------------
    #|40	2	1	0	1	5	1	8	+	q1	100	10	45	chr1	1000	100	138	2	20,10,	10,35,	100,128,
    #|27	2	1	0	2	25	2	20	--	q1	100	10	65	chr2	500	100	150	3	10,10,10,	80,60,35,	390,370,350,
    #|28	1	1	0	1	0	1	0	++	q3	40	5	35	chrP	300	30	150	2	10,20,	5,15,	30,90,
  let qresults = parse(psl, false)

  println("\n1. Parsed structure (coordinates are 0-based, half-open)")
  print_qresults(qresults)

  let stats = @src.blat_psl_write_stats(qresults)
  println(
    "\n2. Write stats: " +
    stats.qresults.to_string() +
    " queries, " +
    stats.hits.to_string() +
    " hits, " +
    stats.hsps.to_string() +
    " HSPs, " +
    stats.fragments.to_string() +
    " fragments",
  )

  println("\n3. Round trip through the writer (header re-emitted):")
  let written = write(qresults, true, false)
  println(written)
  let reparsed = parse(written, false)
  println(
    "re-parsed " +
    reparsed.length().to_string() +
    " query result(s) from writer output",
  )

  // PSLX adds two comma-separated sequence columns per row.
  let pslx =
    #|18	1	1	0	1	4	0	0	+	qx	60	5	29	chrA	900	50	74	2	10,10,	5,19,	50,64,	AAAAATTTTT,CCCCCGGGGG,	AAAAATTT-T,CCCCCGGGGG,
    #|16	2	0	0	0	0	1	6	-+	qx	60	38	58	chrB	700	100	120	2	8,8,	2,14,	100,112,	ACGTACGT,TGCATGCA,	AC-TACGT,TGC-TGCA,
  let with_seqs = parse(pslx, true)
  println("4. PSLX block sequences:")
  let first_block = with_seqs[0].hits[0].hsps[0].blocks[0]
  println("    query=" + first_block.query_seq + " hit=" + first_block.hit_seq)
  let written_pslx = write(with_seqs, false, true)
  println("  PSLX round trip:\n" + written_pslx)

  // The generic SearchIO model uses 1-based inclusive coordinates and maps
  // the UCSC score to bitscore and identity percent to identity.
  let generic = qresults[0].to_search_io()
  let hsp = generic.hits[0].hsps[0]
  println("5. Generic SearchIO view of the first HSP:")
  println(
    "    " +
    generic.id +
    " vs " +
    generic.hits[0].id +
    " query [" +
    hsp.query_start.to_string() +
    ", " +
    hsp.query_end.to_string() +
    "] hit [" +
    hsp.hit_start.to_string() +
    ", " +
    hsp.hit_end.to_string() +
    "] bitscore=" +
    hsp.bitscore.to_string() +
    " identity=" +
    hsp.identity.to_string() +
    "%",
  )
}
