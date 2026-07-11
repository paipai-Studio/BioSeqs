#!/usr/bin/env python3
import json
from Bio.SeqUtils import CheckSum, MeltingTemp
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import gc_fraction

results = {}

# CheckSum tests
seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
results["crc32"] = CheckSum.crc32(seq)
results["gcg"] = CheckSum.gcg(seq)
results["seguid"] = CheckSum.seguid(seq)

# GC content tests
results["gc_simple"] = gc_fraction("AGCT") * 100
results["gc_full"] = gc_fraction(seq) * 100

# ProteinAnalysis tests
protein_seq = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTILTRLRLY"
protein = ProteinAnalysis(protein_seq)
results["protein_molecular_weight"] = protein.molecular_weight()
results["protein_aromaticity"] = protein.aromaticity()
results["protein_gravy"] = protein.gravy()

# Melting temperature tests
results["tm_wallace_simple"] = MeltingTemp.Tm_Wallace("AAAAATTTTT")
results["tm_wallace_gc"] = MeltingTemp.Tm_Wallace("GGCCGGCC")
results["tm_gc_simple"] = MeltingTemp.Tm_GC("AGCTAGCT")
results["tm_gc_long"] = MeltingTemp.Tm_GC(seq)

print(json.dumps(results, indent=2))